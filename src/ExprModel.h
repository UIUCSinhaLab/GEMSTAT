#ifndef EXPR_MODEL_H
#define EXPR_MODEL_H

#include "FactorIntFunc.h"
//#include "ExprPredictor.h"
#include "SeqAnnotator.h"

class ExprFunc; // forward reference to avoid circular issues. #include "ExprFunc.h"
class ExprPar;

typedef enum ModelType
{
    LOGISTIC,                                     // logistic regression
    DIRECT,                                       // direct interaction between TF and BTM, repressor works through BTM
    QUENCHING,                                    // repressor stops activator from interacting with BTM
    CHRMOD_UNLIMITED,                             // repressor works by chromatin modification (making it unaccessible), unlimited activation
    CHRMOD_LIMITED,                                // repressor works by chromatin modification (making it unaccessible), limited activation
    RATES,
    MARKOV
} ModelType;

#include "ExprPar.h"

ModelType getModelOption( const string& modelOptionStr );
string getModelOptionStr( ModelType modelOption );


class ExprModel {

public: //TODO: Implement good accessors / mutators instead.
	ExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs, FactorIntFunc* _intFunc, int _maxContact, IntMatrix& _coopMat, vector< bool >& _actIndicators, vector< bool>& _repIndicators, IntMatrix& _repressionMat, double _repressionDistThr );


	ModelType modelOption;			//model option TODO: remove this later and use inheritance for that.
	bool one_qbtm_per_crm;		//Each crm has its own different binding weight for the polymerase
  bool shared_scaling;

	vector< Motif >& motifs;		//The motifs

	FactorIntFunc* intFunc;			// function to compute distance-dependent TF-TF interactions

	int maxContact;

	IntMatrix& coopMat;		// Cooperativities

	//The roles of TFs. Consider creating a "TF Role" object.
	//Then there would be one vector of roles, and that TFRole object could be sub-classed for different models.
	vector< bool >& actIndicators;		// 1 if the TF is in the activator set
	vector< bool >& repIndicators;		// 1 if the TF is in the repressor set

	//repression related stuff
	IntMatrix& repressionMat;		// repression matrix: R(f,f') = 1 if f can repress f' (Matrix representation of a directed graph.)
	double repressionDistThr;		// distance threshold for repression: d_R

  int getNFactors() const {return motifs.size();}
  int getNumCoop() const;

  ExprFunc* createNewExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const;
};



#endif
