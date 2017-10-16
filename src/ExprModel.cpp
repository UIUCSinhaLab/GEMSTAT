
#include <stdexcept>      // std::invalid_argument

#include "ExprPredictor.h"
#include "ExprModel.h"

ExprModel::ExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs, FactorIntFunc* _intFunc, int _maxContact, IntMatrix& _coopMat, vector< bool >& _actIndicators, vector< bool>& _repIndicators, IntMatrix& _repressionMat, double _repressionDistThr ) : modelOption( _modelOption), one_qbtm_per_crm( _one_qbtm_per_crm), motifs( _motifs), intFunc( _intFunc), maxContact( _maxContact), coopMat( _coopMat ), actIndicators( _actIndicators), repIndicators( _repIndicators), repressionMat( _repressionMat), repressionDistThr( _repressionDistThr)
{
  shared_scaling = false;
  //A QUENCHING model shall have repIndicators all false.
  if(_modelOption == QUENCHING){
      for(int i = 0;i<repIndicators.size();i++){
          if(repIndicators[i])
            throw std::invalid_argument("A quenching model should not have repIndicators that are all false.");
      }
  }
}


ModelType getModelOption( const string& modelOptionStr )
{
    if ( toupperStr( modelOptionStr ) == "LOGISTIC" ) return LOGISTIC;
    if ( toupperStr( modelOptionStr ) == "DIRECT" ) return DIRECT;
    if ( toupperStr( modelOptionStr ) == "QUENCHING" ) return QUENCHING;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_UNLIMITED" ) return CHRMOD_UNLIMITED;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_LIMITED" ) return CHRMOD_LIMITED;
    if ( toupperStr( modelOptionStr ) == "RATES" ) return RATES;
    if ( toupperStr( modelOptionStr ) == "MARKOV" ) return MARKOV;

    cerr << "modelOptionStr is not a valid model option" << endl;
    exit(1);
}


string getModelOptionStr( ModelType modelOption )
{
    if ( modelOption == LOGISTIC ) return "Logisitic";
    if ( modelOption == DIRECT ) return "Direct";
    if ( modelOption == QUENCHING ) return "Quenching";
    if ( modelOption == CHRMOD_UNLIMITED ) return "ChrMod_Unlimited";
    if ( modelOption == CHRMOD_LIMITED ) return "ChrMod_Limited";
    if ( modelOption == RATES ) return "Rates";
    if ( modelOption == MARKOV ) return "Markov";

    return "Invalid";
}

int ExprModel::getNumCoop() const {
  int num_of_coop_pairs = 0;
  //Calculate the number of cooperative pairs.
  for(int i = 0;i < getNFactors();i++)
    for(int j=i;j< getNFactors();j++)
      num_of_coop_pairs += coopMat.getElement(i,j);

  return num_of_coop_pairs;
}

ExprFunc* ExprModel::createNewExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const
{
  ExprPar parToPass;
  ExprFunc* return_exprfunc = NULL;
  switch(this->modelOption) {
    case LOGISTIC :
      parToPass = par.my_factory->changeSpace(par, ENERGY_SPACE);
      return_exprfunc = new Logistic_ExprFunc(sites_,seq_length,seq_num,
                          this->motifs,
                          this->intFunc,
                          this->actIndicators,
                          this->maxContact,
                          this->repIndicators,
                          this->repressionMat,
                          this->repressionDistThr,
                          parToPass );
      break;
    case DIRECT :
      parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
      return_exprfunc = new Direct_ExprFunc(sites_,seq_length,seq_num,
                        this->motifs,
                        this->intFunc,
                        this->actIndicators,
                        this->maxContact,
                        this->repIndicators,
                        this->repressionMat,
                        this->repressionDistThr,
                        parToPass );
      break;
    case QUENCHING :
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new Quenching_ExprFunc(sites_,seq_length,seq_num,
                          this->motifs,
                          this->intFunc,
                          this->actIndicators,
                          this->maxContact,
                          this->repIndicators,
                          this->repressionMat,
                          this->repressionDistThr,
                          parToPass );
        break;
    case CHRMOD_LIMITED :
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new ChrModLimited_ExprFunc(sites_,seq_length,seq_num,
                          this->motifs,
                          this->intFunc,
                          this->actIndicators,
                          this->maxContact,
                          this->repIndicators,
                          this->repressionMat,
                          this->repressionDistThr,
                          parToPass );
        break;
    case CHRMOD_UNLIMITED :
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new ChrModUnlimited_ExprFunc(sites_,seq_length,seq_num,
                          this->motifs,
                          this->intFunc,
                          this->actIndicators,
                          this->maxContact,
                          this->repIndicators,
                          this->repressionMat,
                          this->repressionDistThr,
                          parToPass );
        break;
    case MARKOV:
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new Markov_ExprFunc(sites_,seq_length,seq_num,
                          this->motifs,
                          this->intFunc,
                          this->actIndicators,
                          this->maxContact,
                          this->repIndicators,
                          this->repressionMat,
                          this->repressionDistThr,
                          parToPass );
        break;
    default :
        cerr << "Somehow, an invalid model argument was passed. " << endl;
        assert(false);//Should never reach here.
  }

  return return_exprfunc;
}
