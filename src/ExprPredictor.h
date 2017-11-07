#ifndef EXPR_PREDICTOR_H
#define EXPR_PREDICTOR_H

#include "ExprModel.h"
#include "FactorIntFunc.h"
#include "SeqAnnotator.h"
#include "PredictorTrainer.h"
#include "ExprPar.h"
#include "ObjFunc.h"
#include "ExprFunc.h"
#include "DataSet.h"

/*****************************************************
 * Model Training and Testing
 ******************************************************/

/* ExprPredictor class: the thermodynamic sequence-to-expression predictor */
class ExprPredictor
{
    public:
        // constructors
        ExprPredictor( const vector < Sequence >& _seqs, const vector< SiteVec >& _seqSites,  const vector< int >& _seqLengths, const DataSet& _training_data, const vector< Motif >& _motifs, const ExprModel& _expr_model, const vector < bool >& _indicator_bool, const vector <string>& _motifNames);
        ~ExprPredictor();
        // access methods
        int nSeqs() const { return seqs.size(); }
        int nFactors() const { return expr_model.motifs.size(); }
        int nConds() const { return training_data.nConds(); }
        const IntMatrix& getCoopMat() const { return expr_model.get_coop_mat_immutable(); }
        const vector< bool >& getActIndicators() const { return expr_model.actIndicators; }
        const vector< bool >& getRepIndicators() const { return expr_model.repIndicators; }
        const IntMatrix& getRepressionMat() const { return expr_model.repressionMat; }
        const vector< SiteVec >& getSeqSites(){ return seqSites; }
        const ExprPar& getPar() const { return par_model; }
        double getObj() const { return obj_model; }

        void set_objective_option( ObjType in_obj_option );

        // the objective function to be minimized
        double objFunc( const ExprPar& par ) ;

        // training the model
        int train( const ExprPar& par_init );     // training with the initial values given
                                                  // training with the initial values and allowing random starts
        int train( const ExprPar& par_init, const gsl_rng* rng );
        int train();                              // automatic training: first estimate the initial values, then train

        // predict expression values of a sequence (across all conditions)
        int predict( const ExprPar& par, const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num) const;
        int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num) const;
        //TODO: Implement this such that the previous calls it, it is not slow, and it is DRY and KISS
        //int predict( const SiteVec& targetSites, int targetSeqLength, vector<int> condition_index_list, vector< double >& targetExprs, int seq_num ) const;

        // test the model, perfOption = 0: RMSE
        // 	double test( const vector< Sequence  >& testSeqs, const Matrix& testExprData, Matrix& predictions ) const;

        //std::ofstream gene_crm_fout;

        ObjType objOption;                 // option of the objective function

        // the similarity between two expression patterns, using cross-correlation
        static double exprSimCrossCorr( const vector< double >& x, const vector< double >& y );
        static int maxShift;                      // maximum shift when computing cross correlation
        static double shiftPenalty;               // the penalty for shift (when weighting different positions)

        // the parameters for the optimizer
        int n_alternations;                 // number of alternations (between two optimization methods)
        int n_random_starts;                   // number of random starts
        static double min_delta_f_SSE;            // the minimum change of the objective function under SSE
        static double min_delta_f_Corr;           // the minimum change of the objective function under correlation
        static double min_delta_f_CrossCorr;      // the minimum change of the objective function under cross correlation
        static double min_delta_f_PGP;            // the minimum change of the objective function under PGP
        int max_simplex_iterations;               // maximum number of iterations for Simplex optimizer (default = 200)
        int max_gradient_iterations;              // maximum number of iterations for Gradient optimizer (default = 50)
        vector < bool > indicator_bool;
        vector <string> motifNames;
        vector < double > fix_pars;
        vector < double > free_pars;
        vector < Sequence > seqs;

        //TODO: decide if this needs to be made private
        // Factory for Parameter vectors;
        ParFactory *param_factory;
        ObjFunc *trainingObjective;
        SearchType search_option;
    private:
        //***** training data ******
        const DataSet& training_data;               //input and output curves
        //
        const vector< SiteVec >& seqSites;        // the extracted sites for all sequences
        const vector< int >& seqLengths;          // lengths of all sequences
        //TODO: R_SEQ Either remove this dead feature or revive it and make it conditional.
        //const vector <SiteVec>& r_seqSites;
        //const vector< int >& r_seqLengths;        // lengths of all sequences


        //***** MODEL *******
        //The model, controls what kinds of interactions, who can interact, what kind of BTM, etc.
    	const ExprModel& expr_model;

        // model parameters and the value of the objective function
        ExprPar par_model;
        double obj_model;

        // randomly sample parameter values (only those free parameters), the parameters should be initialized
        //int randSamplePar( const gsl_rng* rng, ExprPar& par ) const;

        // print the parameter values (the ones that are estimated) in a single line
        void printPar( const ExprPar& par ) const;

        // create the expression function
        ExprFunc* createExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const;

        // objective functions
        double evalObjective( const ExprPar& par );

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
