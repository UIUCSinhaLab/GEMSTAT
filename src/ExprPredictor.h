#ifndef EXPR_PREDICTOR_H
#define EXPR_PREDICTOR_H

#include "SeqAnnotator.h"
enum ModelType
{
    LOGISTIC,                                     // logistic regression
    DIRECT,                                       // direct interaction between TF and BTM, repressor works through BTM
    QUENCHING,                                    // repressor stops activator from interacting with BTM
    CHRMOD_UNLIMITED,                             // repressor works by chromatin modification (making it unaccessible), unlimited activation
    CHRMOD_LIMITED                                // repressor works by chromatin modification (making it unaccessible), limited activation
};

ModelType getModelOption( const string& modelOptionStr );
string getModelOptionStr( ModelType modelOption );

enum FactorIntType
{
    BINARY,                                       // Binary model of interaction
    GAUSSIAN                                      // Gaussian model of interaction
};

string getIntOptionStr( FactorIntType intOption );

enum ObjType
{
    SSE,                                          // sum of squared error
    CORR,                                         // Pearson correlation
    CROSS_CORR,                                   // cross correlation (maximum in a range of shifts)
    PGP                                           // PGP score
};

ObjType getObjOption( const string& objOptionStr );
string getObjOptionStr( ObjType objOption );

enum SearchType
{
    UNCONSTRAINED,                                // unconstrained search
    CONSTRAINED                                   // constrained search
};

string getSearchOptionStr( SearchType searchOption );
/*****************************************************
 * Factor-Factor Interactions
 ******************************************************/

/* FactorIntFunc class: distance-dependent function of TF-TF interaction  */
class FactorIntFunc
{
    public:
        // compute the factor interaction, given the normal interaction (when they are close enough)
        virtual double compFactorInt( double normalInt, double dist, bool orientation ) const = 0;

        // the maximum distance beyond which there is no interaction
        virtual double getMaxDist() const = 0;
};

/* FactorIntFuncBinary class: binary distance function */
class FactorIntFuncBinary : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncBinary( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyond which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // if distance < thr, the "normal" value; otherwise 1 (no interaction)
        double orientationEffect;                 // the effect of orientation: if at different strands, the effect should be multiplied this value
};

/* FactorIntFuncGaussian class: Gaussian distance function*/
class FactorIntFuncGaussian : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncGaussian( double _distThr, double _sigma ) : distThr( _distThr ), sigma( _sigma )
        {
            assert( distThr > 0 && sigma > 0 );
        }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyone which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // no interaction if distance is greater than thr.
        double sigma;                             // standard deviation of
};

/* FactorIntFuncGeometric class: distance function decays geometrically (but never less than 1) */
class FactorIntFuncGeometric : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncGeometric( double _distThr, double _spacingEffect, double _orientationEffect ) : distThr( _distThr ), spacingEffect( _spacingEffect ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyond which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // if distance < thr, the "normal" value; otherwise decay with distance (by parameter spacingEffect)
        double spacingEffect;                     // the effect of spacing
        double orientationEffect;                 // the effect of orientation: if at different strands, the effect should be multiplied this value
};

/*****************************************************
 * Expression Model and Parameters
 ******************************************************/

/* ExprPar class: the parameters of the expression model */
class ExprPar
{
    public:
        // constructors
        ExprPar() : factorIntMat() {}
        ExprPar( int _nFactors, int _nSeqs );     // default values of parameters
        ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >&  _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector< double >& _energyThrFactors );
                                                  // construct from a "flat" vector of free parameters (assuming they are in the correct/uniform scale)
        ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs );
        void copy( const ExprPar& other ) { maxBindingWts = other.maxBindingWts; factorIntMat = other.factorIntMat; txpEffects = other.txpEffects; repEffects = other.repEffects; basalTxps = other.basalTxps; pis = other.pis, betas = other.betas, energyThrFactors = other.energyThrFactors, nSeqs = basalTxps.size();  }
        ExprPar( const ExprPar& other ) { copy( other ); }

        // assignment
        ExprPar& operator=( const ExprPar& other ) { copy( other ); return *this; }

        // access methods
        int nFactors() const { return maxBindingWts.size(); }

        // get the free parameters (in the correct/uniform scale)
        virtual void getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const;

        // print the parameters
        virtual void print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const;

        // load the parameter values from a file, assuming the parameter has the correct dimensions (and initialized)
        virtual int load( const string& file, const int num_of_coop_pairs );

        // adjust the values of parameters: if the value is close to min or max allowed value, slightly change it s.t. it is away from the boundary
        virtual void adjust( const IntMatrix& coopMat  );

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

/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class ExprFunc
{
    public:
        // constructors
        ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par );

        // access methods
        const vector< Motif >& getMotifs() const
        {
            return motifs;
        }

        // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
        virtual double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num );
        virtual double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num, int TFid );
        virtual double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num, std::ofstream& fout );
        const ExprPar& getPar() const { return par; }

        static ModelType modelOption;             // model option
        static bool one_qbtm_per_crm;
    protected:
        // TF binding motifs
        const vector< Motif >& motifs;

        // control parameters
        const FactorIntFunc* intFunc;             // function to compute distance-dependent TF-TF interactions
        const vector< bool >& actIndicators;      // 1 if the TF is in the activator set
        int maxContact;                           // the maximum contact
        const vector< bool >& repIndicators;      // 1 if the TF is in the repressor set
        const IntMatrix& repressionMat;           // repression matrix: R(f,f') = 1 if f can repress f'
        double repressionDistThr;                 // distance threshold for repression: d_R

        // model parameters
        const ExprPar& par;

        // the sequence whose expression is to be predicted
        SiteVec sites;
        vector< int > boundaries;                 // left boundary of each site beyond which there is no interaction

        // intermediate computational results
        vector< double > bindingWts;

        // compute the partition function when the basal transcriptional machinery (BTM) is not bound
        virtual double compPartFuncOff() const;

        // compute the partition function when the BTM is not bound: ChrMod model
        virtual double compPartFuncOffChrMod() const;

        // compute the partition function when the BTM is bound
        virtual double compPartFuncOn() const;

        // compute the paritition function when the BTM is bound: Direct model
        virtual double compPartFuncOnDirect() const;

        // compute the paritition function when the BTM is bound: Quenching model
        virtual double compPartFuncOnQuenching() const;

        // compute the paritition function when the BTM is bound: ChrMod_Unlimited model
        virtual double compPartFuncOnChrMod_Unlimited() const;

        // compute the paritition function when the BTM is bound: ChrMod_Limited model
        virtual double compPartFuncOnChrMod_Limited() const;

        // compute the TF-TF interaction between two occupied sites
        virtual double compFactorInt( const Site& a, const Site& b ) const;

        double compDen() const;
        double compNum( int TFid ) const;

        // test if one site represses another site
        bool testRepression( const Site& a, const Site& b ) const;
};

/*****************************************************
 * Model Training and Testing
 ******************************************************/

/* ExprPredictor class: the thermodynamic sequence-to-expression predictor */
class ExprPredictor
{
    public:
        // constructors
        ExprPredictor( const vector < Sequence >& _seqs, const vector< SiteVec >& _seqSites, const vector< SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts  );

        // access methods
        int nSeqs() const
        {
            return seqs.size();
        }
        int nFactors() const
        {
            return motifs.size();
        }
        int nConds() const
        {
            return exprData.nCols();
        }
        const IntMatrix& getCoopMat() const
        {
            return coopMat;
        }
        const vector< bool >& getActIndicators() const
        {
            return actIndicators;
        }
        const vector< bool >& getRepIndicators() const
        {
            return repIndicators;
        }
        const IntMatrix& getRepressionMat() const
        {
            return repressionMat;
        }
        const ExprPar& getPar() const { return par_model; }
        double getObj() const { return obj_model; }

        // the objective function to be minimized
        double objFunc( const ExprPar& par ) ;

        // training the model
        int train( const ExprPar& par_init );     // training with the initial values given
                                                  // training with the initial values and allowing random starts
        int train( const ExprPar& par_init, const gsl_rng* rng );
        int train();                              // automatic training: first estimate the initial values, then train

        // predict expression values of a sequence (across the same conditions)
        int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const;

        // test the model, perfOption = 0: RMSE
        // 	double test( const vector< Sequence  >& testSeqs, const Matrix& testExprData, Matrix& predictions ) const;

        //std::ofstream gene_crm_fout;

        static ModelType modelOption;             // model option
        static int estBindingOption;              // whether estimate binding parameters
        static ObjType objOption;                 // option of the objective function

        // the similarity between two expression patterns, using cross-correlation
        static double exprSimCrossCorr( const vector< double >& x, const vector< double >& y );
        static int maxShift;                      // maximum shift when computing cross correlation
        static double shiftPenalty;               // the penalty for shift (when weighting different positions)

        // the parameters for the optimizer
        static int nAlternations;                 // number of alternations (between two optimization methods)
        static int nRandStarts;                   // number of random starts
        static double min_delta_f_SSE;            // the minimum change of the objective function under SSE
        static double min_delta_f_Corr;           // the minimum change of the objective function under correlation
        static double min_delta_f_CrossCorr;      // the minimum change of the objective function under cross correlation
        static double min_delta_f_PGP;            // the minimum change of the objective function under PGP
        static int nSimplexIters;                 // maximum number of iterations for Simplex optimizer
        static int nGradientIters;                // maximum number of iterations for Gradient optimizer
        static bool one_qbtm_per_crm;
        vector < bool > indicator_bool;
        vector <string> motifNames;
        vector < double > fix_pars;
        vector < double > free_pars;
        vector < Sequence > seqs;
    protected:
        // training data
        const vector< SiteVec >& seqSites;        // the extracted sites for all sequences
        const vector< int >& seqLengths;          // lengths of all sequences
        //TODO: R_SEQ Either remove this dead feature or revive it and make it conditional.
	const vector <SiteVec>& r_seqSites;
        const vector< int >& r_seqLengths;        // lengths of all sequences
        const Matrix& exprData;                   // expressions of the corresponding sequences across multiple conditions
        const vector< Motif >& motifs;            // TF binding motifs
        const Matrix& factorExprData;             // [TF] of all factors over multiple conditions
        const vector < int >& axis_start;
        const vector < int >& axis_end;
        const vector < double >& axis_wts;

        // control parameters
        const FactorIntFunc* intFunc;             // function to compute distance-dependent TF-TF interactions
        const IntMatrix& coopMat;                 // cooperativity matrix: C(f,f') = 1 if f and f' bind cooperatively
        const vector< bool >& actIndicators;      // 1 if the TF is in the activator set
        int maxContact;                           // the maximum contact
        const vector< bool >& repIndicators;      // 1 if the TF is in the repressor set
        const IntMatrix& repressionMat;           // repression matrix: R(f,f') = 1 if f can repress f'
        double repressionDistThr;                 // distance threshold for repression: d_R

        // model parameters and the value of the objective function
        ExprPar par_model;
        double obj_model;

        // randomly sample parameter values (only those free parameters), the parameters should be initialized
        int randSamplePar( const gsl_rng* rng, ExprPar& par ) const;

        // check if some parameter combination is valid
        virtual bool testPar( const ExprPar& par ) const;

        // print the parameter values (the ones that are estimated) in a single line
        virtual void printPar( const ExprPar& par ) const;

        // create the expression function
        virtual ExprFunc* createExprFunc( const ExprPar& par ) const;

        // objective functions
        double compRMSE( const ExprPar& par );    // root mean square error between predicted and observed expressions
        double compAvgCorr( const ExprPar& par ); // the average Pearson correlation
                                                  // the average cross correlation -based similarity
        double compAvgCrossCorr( const ExprPar& par );
        double compPGP( const ExprPar& par );     // the average cross correlation -based similarity

        // minimize the objective function, using the current model parameters as initial values
                                                  // simplex
        int simplex_minimize( ExprPar& par_result, double& obj_result );
                                                  // gradient: BFGS or conjugate gradient
        int gradient_minimize( ExprPar& par_result, double& obj_result );
        //  	int SA_minimize( ExprPar& par_result, double& obj_result ) const;	// simulated annealing
};

// the objective function and its gradient of ExprPredictor::simplex_minimize or gradient_minimize
double gsl_obj_f( const gsl_vector* v, void* params );
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad );
void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad );
#endif
