
#include <stdexcept>      // std::invalid_argument

#include "ExprPredictor.h"
#include "ExprModel.h"

ExprModel::ExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs,int _maxContact, vector< bool >& _actIndicators, vector< bool>& _repIndicators, IntMatrix& _repressionMat, double _repressionDistThr ) : modelOption( _modelOption), one_qbtm_per_crm( _one_qbtm_per_crm), motifs( _motifs), maxContact( _maxContact), actIndicators( _actIndicators), repIndicators( _repIndicators), repressionMat( _repressionMat), repressionDistThr( _repressionDistThr)
{
    int n_tfs = _motifs.size();//might change later? (multi-motif)

    coop_setup = new CoopInfo(n_tfs);

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

ExprFunc* ExprModel::createNewExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const
{
  ExprPar parToPass;
  ExprFunc* return_exprfunc = NULL;
  switch(this->modelOption) {
    case LOGISTIC :
      parToPass = par.my_factory->changeSpace(par, ENERGY_SPACE);
      return_exprfunc = new Logistic_ExprFunc(this, parToPass, sites_,seq_length,seq_num);
      break;
    case DIRECT :
      parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
      return_exprfunc = new Direct_ExprFunc(this, parToPass, sites_,seq_length,seq_num);
      break;
    case QUENCHING :
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new Quenching_ExprFunc(this, parToPass, sites_,seq_length,seq_num);
        break;
    case CHRMOD_LIMITED :
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new ChrModLimited_ExprFunc(this, parToPass, sites_,seq_length,seq_num);
        break;
    case CHRMOD_UNLIMITED :
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new ChrModUnlimited_ExprFunc(this, parToPass, sites_,seq_length,seq_num);
        break;
    case MARKOV:
        parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
        return_exprfunc = new Markov_ExprFunc(this, parToPass, sites_,seq_length,seq_num);
        break;
    default :
        cerr << "Somehow, an invalid model argument was passed. " << endl;
        assert(false);//Should never reach here.
  }

  return return_exprfunc;
}


CoopInfo::CoopInfo(int n_motifs) : coop_matrix( n_motifs, n_motifs, false)
{
    int_funcs.clear();
    int_funcs.push_back(new Null_FactorIntFunc()); //Interaction between sites not otherwise listed.
    int_funcs.push_back(new FactorIntFuncBinary( 20 ));


}

void CoopInfo::set_default_interaction( FactorIntFunc* new_default){
    assert(NULL != new_default);
    FactorIntFunc* tmp = int_funcs[1];
    delete tmp;
    int_funcs[1] = new_default;

}

int CoopInfo::get_longest_coop_thr() const {
    int longest = 0;
    for(int i = 1;i<int_funcs.size();i++){
        int m_dist = int_funcs[i]->getMaxDist();
        if(m_dist > longest){ longest = m_dist; }
    }
    return longest;
}

void CoopInfo::read_coop_file(string filename, map<string, int> factorIdxMap){

        ifstream fin( filename.c_str() );

        std:string line;
        std::istringstream line_ss;
        vector<string> tokens;
        #define LOCAL_TOKENIZE(M_TOK_VECT,M_LINE_STR,M_SS) M_TOK_VECT.clear();\
                                M_SS.clear();\
                                M_SS.str(M_LINE_STR);\
                                copy(istream_iterator<string>(M_SS),\
                                istream_iterator<string>(),\
                                back_inserter(M_TOK_VECT))

        while(std::getline(fin,line)){
                LOCAL_TOKENIZE(tokens,line,line_ss);
                int tf_i = factorIdxMap[tokens[0]];
                int tf_j = factorIdxMap[tokens[1]];

                coop_matrix.setElement(tf_i,tf_j,true);
                coop_matrix.setElement(tf_j,tf_i,true);

                //cerr << "DEBUG read" << tokens[0] << tf_i << " " << tokens[1] << tf_i << endl;
        }

        //cerr << coop_matrix << endl;

        //exit(1);
        fin.close();
}
