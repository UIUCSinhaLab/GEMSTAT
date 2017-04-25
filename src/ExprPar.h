#ifndef EXPRPAR_H
#define EXPRPAR_H

#include "ExprModel.h"
#include "PredictorTrainer.h"

class ParFactory;
class ExprModel;

enum ThermodynamicParameterSpace {
  PROB_SPACE, //GEMSTAT dynamic programming algorithms are written for this to be the native parameter space.
              //It is also easier for humans to understand. The exp() of ENERGY_SPACE.
  ENERGY_SPACE,//Most thermodynamic textbooks and concepts happen in this parameter space, (-inf, +inf)
                //In particular, regularization, machine learning, etc. generally happen in this space. The log() of the PROB_SPACE.
  CONSTRAINED_SPACE//ENERGY_SPACE parameters are transformed via a Sinusoidal, sigmoidal, or other function so that an unconstrained
                //optimizer can be used to perform constrained optimization. (For user understandability, constraints are specified in PROB_SPACE, _in input, output, and config files_.)
};

string parameterSpaceStr(ThermodynamicParameterSpace in);

/* ExprPar class: the parameters of the expression model */
class ExprPar
{
    friend class ParFactory;
    public:
        // constructors
        ExprPar() : factorIntMat() {}
        ExprPar( int _nFactors, int _nSeqs );     // default values of parameters
        ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >&  _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector< double >& _energyThrFactors );
                                                  // construct from a "flat" vector of free parameters (assuming they are in the correct/uniform scale)
        ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs );
        void copy( const ExprPar& other ) { maxBindingWts = other.maxBindingWts; factorIntMat = other.factorIntMat; txpEffects = other.txpEffects; repEffects = other.repEffects; basalTxps = other.basalTxps; pis = other.pis, betas = other.betas, energyThrFactors = other.energyThrFactors, nSeqs = basalTxps.size(); my_space = other.my_space; my_factory = other.my_factory; }
        ExprPar( const ExprPar& other ) { copy( other ); }

        // assignment
        ExprPar& operator=( const ExprPar& other ) { copy( other ); return *this; }

        // access methods
        int nFactors() const { return maxBindingWts.size(); }

        double getBetaForSeq(int seqID) const; //Returns the appropriate value of beta for this sequence. This allows some sequences to share one beta value, while others share another, or each have their own.

        //get the parameters into a vector.
        void getRawPars(vector< double >& pars) const;//Temporary until these edits are all done.

        // print the parameters
        void print(  ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const;

        //Which parameter space this ExprPar lives in
        ThermodynamicParameterSpace my_space;
        const ParFactory* my_factory;

        // parameters
        vector < double > maxBindingWts;          // binding weight of the strongest site for each TF: K(S_max) [TF_max]
        Matrix factorIntMat;                      // (maximum) interactions between pairs of factors: omega(f,f')
        vector < double > txpEffects;             // transcriptional effects: alpha for Direct and Quenching model, exp(alpha) for Logistic model (so that the same default values can be used). Equal to 1 if a TF is not an activator under the Quenching model
        vector < double > repEffects;             // repression effects: beta under ChrMod models (the equlibrium constant of nucleosome association with chromatin). Equal to 0 if a TF is not a repressor.
        vector < double > basalTxps;              // basal transcription: q_p for Direct and Quenching model, exp(alpha_0) for Logistic model (so that the same default value can be used)
        vector < double > pis;
        //     double expRatio; 		// constant factor of measurement to prediction

        vector < double > betas;
        vector < double > energyThrFactors;
        int nSeqs;
        
        static ModelType modelOption;             // model option
        static SearchType searchOption;           // search option: 0 - unconstrained search; 1 - constrained search
        static int estBindingOption;              // whether to estimate binding parameters
        static bool one_qbtm_per_crm;

        static double default_weight;             // default binding weight
        static double default_interaction;        // default factor interaction
        static double default_effect_Logistic;    // default transcriptional effect under Logistic model
        static double default_effect_Thermo;      // default transcriptional effect under thermo. models
        static double default_repression;         // default repression
        static double default_basal_Logistic;     // default basal transcription under Logistic model
        static double default_basal_Thermo;       // default basal transcriptional under thermo. models
        static double default_pi;
        static double min_pi;
        static double max_pi;
        static double min_weight;                 // min. binding weight
        static double max_weight;                 // max. binding weight
        static double min_interaction;            // min. interaction
        static double max_interaction;            // max. interaction
        static double min_effect_Logistic;        // min. transcriptional effect under Logistic model
        static double max_effect_Logistic;        // max. transcriptional effect under Logistic model
        //     static double min_effect_Direct;   // min. transcriptional effect under Direct model
        static double min_effect_Thermo;          // min. transcriptional effect under thermo. models
        static double max_effect_Thermo;          // max. transcriptional effect under thermo. models
        static double min_repression;             // min. repression
        static double max_repression;             // max. repression
        static double min_basal_Logistic;         // min. basal transcription under Logistic model
        static double max_basal_Logistic;         // max. basal transcription under Logistic model
        static double min_basal_Thermo;           // min. basal transcription under thermo. models
        static double max_basal_Thermo;           // max. basal transcription under thermo. models
        static double min_energyThrFactors;
        static double max_energyThrFactors;
        static double default_energyThrFactors;
        static double delta;                      // a small number for testing the range of parameter values
        static double default_beta;
        static double min_beta;                   // a small number for testing the range of parameter values
        static double max_beta;                   // a small number for testing the range of parameter values
        // 	static double wt_step;		// step of maxExprWt (log10)
        // 	static double int_step;		// step of interaction (log10)
        // 	static double ratio_step;	// step of expRatio
};

class ParFactory
{
    friend class ExprPar;
    public:
      ParFactory( const ExprModel& in_model, int in_nSeqs, const vector<bool>& in_ff);
      ~ParFactory(){};

      //the raison d'etre for this class
      virtual ExprPar create_expr_par() const;
      virtual ExprPar create_expr_par(const vector<double>& pars, const ThermodynamicParameterSpace in_space) const;




      void setFreeFix(const vector<bool>& in_free_fix);

      ExprPar createDefaultMinMax(bool min_or_max) const;

      ExprPar changeSpace(const ExprPar& in_par, const ThermodynamicParameterSpace new_space) const;

      //Code for separating parameters to optimize from those that we don't want to optimize.
      void joinParams(const vector<double>& freepars, const vector<double>& fixpars, vector<double>& output, const vector<bool>& indicator_bool) const;//TODO: Will become unnecessary when we switch to a natrually constrained optimizer.
      void separateParams(const ExprPar& input, vector<double>& free_output, vector<double>& fixed_output, const vector<bool>& indicator_bool) const;

      ExprPar truncateToBounds(const ExprPar& in_par, const vector<bool>& indicator_bool) const;
      bool testWithinBounds(const ExprPar& in_par) const;

      ExprPar randSamplePar( const gsl_rng* rng) const;

      int nFactors() const;

      void setMaximums(const ExprPar& in_maximums){assert(in_maximums.my_space == ENERGY_SPACE); maximums.copy(in_maximums);}
      void setMinimums(const ExprPar& in_minimums){assert(in_minimums.my_space == ENERGY_SPACE); minimums.copy(in_minimums);}
      const ExprPar& getMaximums(){return maximums;}
      const ExprPar& getMinimums(){return minimums;}
      const ExprPar& getDefaults(){return defaults;}

      ExprPar load(const string& file);//Use this.

      const ExprModel& expr_model;
    private:
      const vector<bool>& indicator_bool;//TODO: Later work on making this non-const. Original author const'ed everything.
      int nSeqs;

      ExprPar maximums; //Should be in the ENERGY_SPACE
      ExprPar minimums; //Should be in the ENERGY_SPACE
      ExprPar defaults; //In the ENERGY_SPACE

      void constrained_to_energy_helper(const vector<double>& pars, vector<double>& output, const vector<double>& low, const vector<double>& high) const;
      void energy_to_constrained_helper(const vector<double>& pars, vector<double>& output, const vector<double>& low, const vector<double>& high) const;
      void energy_to_prob_helper(const vector<double>& pars, vector<double>& output) const;
      void prob_to_energy_helper(const vector<double>& pars, vector<double>& output) const;
};

#endif
