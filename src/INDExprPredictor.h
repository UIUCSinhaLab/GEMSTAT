#ifndef IND_EXPR_PREDICTOR_H
#define IND_EXPR_PREDICTOR_H

#include "ExprPredictor.h"

/*****************************************************
* Expression Model and Parameters
******************************************************/

/* ExprPar class: the parameters of the expression model */
class INDExprPar : public ExprPar {
public:
        // constructors
        INDExprPar() : ExprPar() {}
	INDExprPar( const ExprPar& other) : ExprPar(other) {}
        INDExprPar( int _nFactors, int _nSeqs );     // default values of parameters
    INDExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >&  _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector< double >& _energyThrFactors, double _cic_att );
    INDExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs );	// construct from a "flat" vector of free parameters (assuming they are in the correct/uniform scale)
    void copy( const INDExprPar& other ) { maxBindingWts = other.maxBindingWts; factorIntMat = other.factorIntMat; txpEffects = other.txpEffects; repEffects = other.repEffects; basalTxps = other.basalTxps; pis = other.pis; betas = other.betas; energyThrFactors = other.energyThrFactors; cic_att = other.cic_att; nSeqs = basalTxps.size();  }
        INDExprPar( const INDExprPar& other ) { copy( other ); }

        // assignment
        INDExprPar& operator=( const INDExprPar& other ) { copy( other ); return *this; }

	 virtual void getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const;

	// print the parameters
        virtual void print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const;

        // load the parameter values from a file, assuming the parameter has the correct dimensions (and initialized)
        virtual int load( const string& file, const int num_of_coop_pairs );

        // adjust the values of parameters: if the value is close to min or max allowed value, slightly change it s.t. it is away from the boundary
        virtual void adjust( const IntMatrix& coopMat  );	


	double cic_att;
    static double min_cic_att;
    static double max_cic_att;
    static double default_cic_att;
};

/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class INDExprFunc : public ExprFunc{
public:
        // constructors
        INDExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const INDExprPar& _par );

        // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
    	double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num, double dperk_conc );
        double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num );
        const INDExprPar& getPar() const { return par; }

	void set_dperk_expr( const Matrix& _dperk_ExprData );

	protected:
		// model parameters
        const INDExprPar& par;
	Matrix& dperk_ExprData;
};

/*****************************************************
* Model Training and Testing
******************************************************/

/* ExprPredictor class: the thermodynamic sequence-to-expression predictor */
class INDExprPredictor : public ExprPredictor {
public:
        // constructors
    INDExprPredictor( const vector < Sequence >& _seqs, const vector< SiteVec >& _seqSites, const vector< SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const Matrix& _dperk_ExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts  );

	//TODO: If we need these because of the C++ object model, alter them, otherwise remove
        /*
        const ExprPar& getPar() const { return par_model; }
        double getObj() const { return obj_model; }

        // the objective function to be minimized
        double objFunc( const ExprPar& par ) ;

        // training the model
        int train( const ExprPar& par_init );     // training with the initial values given
    int train( const ExprPar& par_init, const gsl_rng* rng );   // training with the initial values and allowing random starts
        int train();                              // automatic training: first estimate the initial values, then train

        // predict expression values of a sequence (across the same conditions)
        int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const;

        // test the model, perfOption = 0: RMSE
// 	double test( const vector< Sequence  >& testSeqs, const Matrix& testExprData, Matrix& predictions ) const;    

	std::ofstream gene_crm_fout;

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
	*/
    const Matrix& dperk_ExprData;
private:
        // model parameters and the value of the objective function
        INDExprPar par_model;

	/*
        // randomly sample parameter values (only those free parameters), the parameters should be initialized
        int randSamplePar( const gsl_rng* rng, ExprPar& par ) const;
	*/
        // check if some parameter combination is valid
        virtual bool testPar( const ExprPar& par ) const;

	
        // print the parameter values (the ones that are estimated) in a single line
        virtual void printPar( const ExprPar& par ) const;

        // create the expression function
        virtual ExprFunc* createExprFunc( const ExprPar& par ) const;
	
	/*
        // objective functions
        double compRMSE( const ExprPar& par );    // root mean square error between predicted and observed expressions
        double compAvgCorr( const ExprPar& par ); // the average Pearson correlation
    double compAvgCrossCorr( const ExprPar& par );    // the average cross correlation -based similarity
        double compPGP( const ExprPar& par );     // the average cross correlation -based similarity

        // minimize the objective function, using the current model parameters as initial values
    int simplex_minimize( ExprPar& par_result, double& obj_result );	// simplex	
    int gradient_minimize( ExprPar& par_result, double& obj_result );	// gradient: BFGS or conjugate gradient
//  	int SA_minimize( ExprPar& par_result, double& obj_result ) const;	// simulated annealing 		
    int optimize_beta( ExprPar& par_result, double& obj_result);	// find the current best beta with one-step otimization.
	*/
};

#endif
