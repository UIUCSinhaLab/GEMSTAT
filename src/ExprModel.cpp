#include "ExprPredictor.h"

#include "ExprModel.h"

ExprModel::ExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs, FactorIntFunc* _intFunc, int _maxContact, IntMatrix& _coopMat, vector< bool >& _actIndicators, vector< bool>& _repIndicators, IntMatrix& _repressionMat, double _repressionDistThr ) : modelOption( _modelOption), one_qbtm_per_crm( _one_qbtm_per_crm), motifs( _motifs), intFunc( _intFunc), maxContact( _maxContact), coopMat( _coopMat ), actIndicators( _actIndicators), repIndicators( _repIndicators), repressionMat( _repressionMat), repressionDistThr( _repressionDistThr)
{
  shared_scaling = false;
}


ModelType getModelOption( const string& modelOptionStr )
{
    if ( toupperStr( modelOptionStr ) == "LOGISTIC" ) return LOGISTIC;
    if ( toupperStr( modelOptionStr ) == "DIRECT" ) return DIRECT;
    if ( toupperStr( modelOptionStr ) == "QUENCHING" ) return QUENCHING;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_UNLIMITED" ) return CHRMOD_UNLIMITED;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_LIMITED" ) return CHRMOD_LIMITED;

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
