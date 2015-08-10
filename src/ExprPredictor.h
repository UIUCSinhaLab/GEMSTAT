#ifndef EXPR_PREDICTOR_H
#define EXPR_PREDICTOR_H

#include "ExprModel.h"
#include "FactorIntFunc.h"
#include "SeqAnnotator.h"
#include "PredictorTrainer.h"
#include "ExprPar.h"

/*****************************************************
 * Expression Model and Parameters
 ******************************************************/


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
        double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num );
        const ExprPar& getPar() const { return par; }

        static ModelType modelOption;             // model option
        static bool one_qbtm_per_crm;
    private:
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
        double compPartFuncOff() const;

        // compute the partition function when the BTM is not bound: ChrMod model
        double compPartFuncOffChrMod() const;

        // compute the partition function when the BTM is bound
        double compPartFuncOn() const;

        // compute the paritition function when the BTM is bound: Direct model
        double compPartFuncOnDirect() const;

        // compute the paritition function when the BTM is bound: Quenching model
        double compPartFuncOnQuenching() const;

        // compute the paritition function when the BTM is bound: ChrMod_Unlimited model
        double compPartFuncOnChrMod_Unlimited() const;

        // compute the paritition function when the BTM is bound: ChrMod_Limited model
        double compPartFuncOnChrMod_Limited() const;

        // compute the TF-TF interaction between two occupied sites
        double compFactorInt( const Site& a, const Site& b ) const;

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
        ExprPredictor( const vector < Sequence >& _seqs, const vector< SiteVec >& _seqSites, const vector< SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const ExprModel& _expr_model, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts  );
        ~ExprPredictor();
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
            return expr_model.coopMat;
        }
        const vector< bool >& getActIndicators() const
        {
            return expr_model.actIndicators;
        }
        const vector< bool >& getRepIndicators() const
        {
            return expr_model.repIndicators;
        }
        const IntMatrix& getRepressionMat() const
        {
            return expr_model.repressionMat;
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
        vector < bool > indicator_bool;
        vector <string> motifNames;
        vector < double > fix_pars;
        vector < double > free_pars;
        vector < Sequence > seqs;

        //TODO: decide if this needs to be made private
        // Factory for Parameter vectors;
        ParFactory *param_factory;
    private:
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
	      const ExprModel& expr_model;

        // model parameters and the value of the objective function
        ExprPar par_model;
        double obj_model;

        // randomly sample parameter values (only those free parameters), the parameters should be initialized
        int randSamplePar( const gsl_rng* rng, ExprPar& par ) const;

        // check if some parameter combination is valid
        bool testPar( const ExprPar& par ) const;

        // print the parameter values (the ones that are estimated) in a single line
        void printPar( const ExprPar& par ) const;

        // create the expression function
        ExprFunc* createExprFunc( const ExprPar& par ) const;

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
	int optimize_beta( ExprPar& par_result, double& obj_result); // find the current best beta with one-step otimization.
};

// the objective function and its gradient of ExprPredictor::simplex_minimize or gradient_minimize
double gsl_obj_f( const gsl_vector* v, void* params );
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad );
void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad );
#endif
