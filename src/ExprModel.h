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


class CoopInfo {
    public:
        CoopInfo(int n_motifs);

        //READ Cooperativities
        void read_coop_file(string filename, map<string, int> factorIdxMap);

        //STORE COOPS
        IntMatrix coop_matrix;
        vector< FactorIntFunc* > int_funcs;

        FactorIntFunc* coop_func_for(int i, int j){ return int_funcs[coop_matrix(i,j)]; }


        void set_default_interaction( FactorIntFunc* new_default);
        int get_num_coops() const {return num_coops;}

        int get_longest_coop_thr() const;

        const IntMatrix& get_coop_mat_immutable() const { return coop_matrix;}
        bool has_coop(int i, int j) { return coop_matrix(i,j) > 0;}
private:
        int num_coops;//cached
};

class ExprModel {

public: //TODO: Implement good accessors / mutators instead.
	ExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs, vector< string>& _motif_names, int _maxContact, vector< bool >& _actIndicators, vector< bool>& _repIndicators, IntMatrix& _repressionMat, double _repressionDistThr );


	ModelType modelOption;			//model option TODO: remove this later and use inheritance for that.
	bool one_qbtm_per_crm;		//Each crm has its own different binding weight for the polymerase
    bool shared_scaling;

	vector< Motif >& motifs;		//The motifs
    vector< string > motifnames;

    CoopInfo* coop_setup;

	int maxContact;        //maximum number of TFs that can interact with the BTM at any one time.

	//The roles of TFs. Consider creating a "TF Role" object.
	//Then there would be one vector of roles, and that TFRole object could be sub-classed for different models.
	vector< bool >& actIndicators;		// 1 if the TF is in the activator set
	vector< bool >& repIndicators;		// 1 if the TF is in the repressor set

	//repression related stuff
	IntMatrix& repressionMat;		// repression matrix: R(f,f') = 1 if f can repress f' (Matrix representation of a directed graph.)
	double repressionDistThr;		// distance threshold for repression: d_R

  int getNFactors() const {return motifs.size();}
  int getNumCoop() const { return (coop_setup == NULL ? 0 : coop_setup->get_num_coops());}

  const IntMatrix& get_coop_mat_immutable() const { return coop_setup->get_coop_mat_immutable();}
  int get_longest_coop_thr() const {return coop_setup->get_longest_coop_thr();}

  void set_coop_setup(CoopInfo* in_pointer){if(coop_setup != NULL){delete coop_setup;} coop_setup = in_pointer;}

  ExprFunc* createNewExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const;
};





#endif
