#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include <nlopt.hpp>

#include "ExprPredictor.h"
#include "ExprPar.h"
#include "ExprFunc.h"

double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data);

ExprPredictor::ExprPredictor( const vector <Sequence>& _seqs, const vector< SiteVec >& _seqSites, const vector < SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const DataSet& _training_data, const vector< Motif >& _motifs, const ExprModel& _expr_model,
		const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts ) : seqs(_seqs), seqSites( _seqSites ), r_seqSites( _r_seqSites ), seqLengths( _seqLengths ), r_seqLengths( _r_seqLengths ), training_data( _training_data ),
	expr_model( _expr_model),
	indicator_bool ( _indicator_bool ), motifNames ( _motifNames ), axis_start ( _axis_start ), axis_end( _axis_end ), axis_wts( _axis_wts )
{
    //TODO: Move appropriate lines from this block to the ExprModel class.
		cerr << "exprData size: " << training_data.exprData.nRows() << "  " << nSeqs() << endl;
    assert( training_data.exprData.nRows() == nSeqs() );
    assert( training_data.factorExprData.nRows() == nFactors() && training_data.factorExprData.nCols() == nConds() );
    assert( expr_model.coopMat.isSquare() && expr_model.coopMat.isSymmetric() && expr_model.coopMat.nRows() == nFactors() );
    assert( expr_model.actIndicators.size() == nFactors() );
    assert( expr_model.maxContact > 0 );
    assert( expr_model.repIndicators.size() == nFactors() );
    assert( expr_model.repressionMat.isSquare() && expr_model.repressionMat.nRows() == nFactors() );
    assert( expr_model.repressionDistThr >= 0 );

    //gene_crm_fout.open( "gene_crm_fout.txt" );

    // set the model option for ExprPar and ExprFunc
    ExprPar::modelOption = expr_model.modelOption;//TODO: Remove both of these.
    ExprFunc::modelOption = expr_model.modelOption;

    // set the values of the parameter range according to the model option
    if ( expr_model.modelOption != LOGISTIC && expr_model.modelOption != DIRECT )
    {
        //ExprPar::min_effect_Thermo = 0.99;
        //ExprPar::min_interaction = 0.99;
    }

    //expr_model was already initialized. Setup the parameter factory.
    param_factory = new ParFactory(expr_model, nSeqs());

    //TODO: Move this to the front-end or something?
    //Maybe make it have a default SSE score objective, but anything else gets specified in the front-end.
    switch(objOption){
      case CORR:
        trainingObjective = new AvgCorrObjFunc();
        break;
      case PGP:
        trainingObjective = new PGPObjFunc();
        break;
      case CROSS_CORR:
        trainingObjective = new AvgCrossCorrObjFunc(ExprPredictor::maxShift, ExprPredictor::shiftPenalty);
        break;
      case LOGISTIC_REGRESSION:
        trainingObjective = new LogisticRegressionObjFunc();
        break;
			case PEAK_WEIGHTED:
				trainingObjective = new PeakWeightedObjFunc();
				break;
      case SSE:
      default:
        trainingObjective = new RMSEObjFunc();
        break;
    }

    /* DEBUG
    cout << setprecision(10);
    ExprPar foo = param_factory->createDefaultMinMax(true);
    cout << " MAXIMUMS " << endl;
    printPar(foo);
    foo = param_factory->createDefaultMinMax(false);
    cout << " MINIMUMS " << endl;
    printPar(foo);
    */
}

ExprPredictor::~ExprPredictor()
{
  delete param_factory;
  delete trainingObjective;
}

double ExprPredictor::objFunc( const ExprPar& par )
{
    double objective_value = evalObjective( par );

    return objective_value;
}


int ExprPredictor::train( const ExprPar& par_init )
{
    par_model = par_init;

    cout << "*** Diagnostic printing BEFORE adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ){
      par_model = param_factory->truncateToBounds(par_model, indicator_bool);

    }
    obj_model = objFunc( par_model );

    cout << "*** Diagnostic printing AFTER adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( nAlternations == 0 ) return 0;

    // alternate between two different methods
    ExprPar par_result;
    double obj_result;
    for ( int i = 0; i < nAlternations; i++ )
    {
        simplex_minimize( par_result, obj_result );
        par_model = par_result;

        gradient_minimize( par_result, obj_result );
        par_model = par_result;
    }

    #ifdef BETAOPTBROKEN
    optimize_beta( par_model, obj_result );
    #endif

    // commit the parameters and the value of the objective function
    //par_model = par_result;
    obj_model = obj_result;

    return 0;
}


int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
    /*
        //for random starts:
        ExprPar par_rand_start = par_init;
        par_rand_start = param_factor->randSamplePar( rng );
        train( par_rand_start );*/
    // training using the initial values
    train( par_init );

    cout << "Initial training:\tParameters = "; printPar( par_model );
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    // training with random starts
    ExprPar par_best = par_model;
    double obj_best = obj_model;
    for ( int i = 0; i < nRandStarts; i++ )
    {
        ExprPar par_curr = par_init;
        par_curr = param_factory->randSamplePar( rng );
        train( par_curr );
        cout << "Random start " << i + 1 << ":\tParameters = "; printPar( par_model );
        cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
        if ( obj_model < obj_best )
        {
            par_best = par_model;
            obj_best = obj_model;
        }
    }

    // training using the best parameters so far
    if ( nRandStarts ) train( par_best );
    cout << "Final training:\tParameters = "; printPar( par_model );
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    //gene_crm_fout.close();

    return 0;
}


int ExprPredictor::train()
{
    // random number generator
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;     // create rng type
    rng = gsl_rng_alloc( T );
    gsl_rng_set( rng, time( 0 ) );                // set the seed equal to simulTime(0)

    // training using the default initial values with random starts
    ExprPar par_default( nFactors(), nSeqs() );
    train( par_default, rng );

    return 0;
}

//Called fairly rarely, don't worry about optimality.
int ExprPredictor::predict( const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num) const
{
	return this->predict(par_model,targetSites_,targetSeqLength,targetExprs,seq_num);
}

int ExprPredictor::predict( const ExprPar& par, const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const
{
	// predict the expression
    ExprFunc* func = createExprFunc( par , targetSites_, targetSeqLength, seq_num);
	targetExprs.resize(nConds());
    for ( int j = 0; j < nConds(); j++ )
    {
		Condition concs = training_data.getCondition( j );
        double predicted = func->predictExpr( concs );
        targetExprs[j] = ( predicted );
    }

    delete func;
    return 0;
}


ObjType ExprPredictor::objOption = SSE;

int ExprPredictor::maxShift = 5;
double ExprPredictor::shiftPenalty = 0.8;

int ExprPredictor::nAlternations = 4;
int ExprPredictor::nRandStarts = 5;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
double ExprPredictor::min_delta_f_PGP = 1.0E-8;
int ExprPredictor::nSimplexIters = 200;
int ExprPredictor::nGradientIters = 50;



void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 8 );
    //     cout.width( 8 );

    // print binding weights
    cout << "MAXBIND : " << par.maxBindingWts << endl;
    cout << "INTER : " ;
    // print the interaction matrix
    for ( int i = 0; i < nFactors(); i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
           cout << par.factorIntMat( i, j ) << "\t";
        }
    }
    cout << endl;

    // print the transcriptional effects
    cout << "TXP : " << par.txpEffects << endl;

    // print the repression effects
    cout << "REP : " << par.repEffects << endl;

    // print the basal transcriptions
    cout << "BASAL : " << par.basalTxps << endl;

    //print the pi values
    cout << "PIS : " << par.pis << endl;

    //print the beta values
    cout << "BETAS : " << par.betas << endl;
    //assert( par.betas.size() == nSeqs() );

    cout << "THRESH : " << par.energyThrFactors << endl;
    cout << flush;
}


ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const
{

    return expr_model.createNewExprFunc( par, sites_, seq_length, seq_num );
}


int indices_of_crm_in_gene[] =
{
};

double ExprPredictor::evalObjective( const ExprPar& par )
{

    vector< int > seqLengths( seqs.size() );
    for( int i = 0; i < seqs.size(); i++ ){
      seqLengths[i] = seqs[i].size();
    }

    vector< SiteVec > seqSites( seqs.size() ); //
    #ifdef REANNOTATE_EACH_PREDICTION
    SeqAnnotator ann( expr_model.motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ ) {
       	ann.annot( seqs[ i ], seqSites[ i ] );
    }
    #else
    seqSites = this->seqSites;
    #endif

    vector<vector<double> > ground_truths;
    vector<vector<double> > predictions;

    //Create predictions for every sequence and condition
    for ( int i = 0; i < nSeqs(); i++ ) {
        ground_truths.push_back(training_data.exprData.getRow(i));
		vector<double> one_seq_predictions(nConds());

		this->predict(par, seqSites[i], seqLengths[i], one_seq_predictions, i );

		predictions.push_back(one_seq_predictions);
    }

    //Evaluate the objective function on that.
    double ret_val = trainingObjective->eval(ground_truths, predictions, &par);

    return ret_val;

}

int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector < double > pars;

    //ExprPar tmp_par_model = param_factory->changeSpace(par_model, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar tmp_par_model = param_factory->changeSpace(par_model, ENERGY_SPACE);
    param_factory->separateParams(tmp_par_model, free_pars, fix_pars, indicator_bool );

    pars.clear();
    pars = free_pars;

    //SIMPLEX MINIMIZATION with NLOPT
    nlopt::opt optimizer(nlopt::LN_NELDERMEAD, pars.size());
    optimizer.set_min_objective(nlopt_obj_func, this);
    //optimizer.set_maxmaxeval(nSimplexIters);
    optimizer.set_initial_step(1.0);//TODO: enforce simplex staring size.
    optimizer.set_maxeval(nSimplexIters);//TODO: enforce nSimplexIters

    if(ExprPar::searchOption == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }

    nlopt::result result = optimizer.optimize(free_pars, obj_result);
    obj_result = optimizer.last_optimum_value();
    //Done Minimizing

    param_factory->joinParams(free_pars, fix_pars, pars, indicator_bool);
    //tmp_par_model = param_factory->create_expr_par(pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    tmp_par_model = param_factory->create_expr_par(pars, ENERGY_SPACE);
    par_result = param_factory->changeSpace(tmp_par_model, PROB_SPACE);

    printPar( par_result );

    return 0;
}


int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    //cout << "DEBUG: in getFreePars()" << endl;
    //par_model.getFreePars( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators );
    //cout << "DEBUG: out getFreePars()" << endl;
    //ExprPar tmp_par_model = param_factory->changeSpace(par_model, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar tmp_par_model = param_factory->changeSpace(par_model, ENERGY_SPACE);

    param_factory->separateParams(tmp_par_model, free_pars, fix_pars, indicator_bool );

    pars.clear();
    pars = free_pars;

    //GRADIENT MINIMIZATION with NLOPT
    nlopt::opt optimizer(nlopt::LD_LBFGS, pars.size());
    optimizer.set_min_objective(nlopt_obj_func, this);

    //TODO: Move this to a nice lookup table or something.
    //Set the stopping criterion
    double ftol;
    switch(objOption){
      SSE:
        ftol = min_delta_f_SSE;
        break;
      CORR:
        ftol = min_delta_f_Corr;
        break;
      CROSS_CORR:
        ftol = min_delta_f_CrossCorr;
        break;
      PGP:
        ftol = min_delta_f_PGP;
        break;
      default:
        ftol = 1e-5;
        break;
    }
    optimizer.set_ftol_abs(ftol);

    if(ExprPar::searchOption == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }


    //TODO: enforce nGradientIters
    optimizer.set_maxeval(nGradientIters);
    try{
      nlopt::result result = optimizer.optimize(free_pars, obj_result);
      obj_result = optimizer.last_optimum_value();
    }catch(std::runtime_error){
      cerr << "There was an exception in the gradient descent!" << endl;
    }

    //Done Minimizing
    //pars now contains the optimal parameters

    param_factory->joinParams(free_pars, fix_pars, pars, indicator_bool);
    //tmp_par_model = param_factory->create_expr_par(pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    tmp_par_model = param_factory->create_expr_par(pars, ENERGY_SPACE);

    par_result = param_factory->changeSpace(tmp_par_model, PROB_SPACE);
    cout << "DEBUG" << endl;
    return 0;
}


double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data){
        gsl_vector *xv = vector2gsl(x); //TODO: Ugly, remove (Make all objective functions use native STL vectors)
        double objective = gsl_obj_f(xv, f_data);

        if(!grad.empty()){
                gsl_vector *dxv = vector2gsl(grad);

                gsl_obj_df(xv,f_data,dxv);

                for(int i = 0;i< grad.size();i++){
                        grad[i] = dxv->data[i];
                }
                cerr << " obj called, derivative : " << grad << endl;
                gsl_vector_free(dxv);
        }

        cerr << " obj called " << objective << endl;

        gsl_vector_free(xv);
        return objective;
}

double gsl_obj_f( const gsl_vector* v, void* params )
{
    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;

    // parse the variables (parameters to be optimized)
    //     vector< double > expv;
    //     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
    vector <double> temp_free_pars = gsl2vector(v);
    vector < double > all_pars;

    predictor->param_factory->joinParams(temp_free_pars, predictor->fix_pars, all_pars, predictor->indicator_bool);
    //ExprPar par = predictor->param_factory->create_expr_par(all_pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar par = predictor->param_factory->create_expr_par(all_pars, ENERGY_SPACE);
    par = predictor->param_factory->changeSpace(par, PROB_SPACE); //TODO: WTF? This shouldn't be required because it's done in the createExprFunc method. Stack corruption or something?


    // call the ExprPredictor object to evaluate the objective function
    double obj = predictor->objFunc( par );
    return obj;
}

void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    double step = 1.0E-4;
    numeric_deriv( grad, gsl_obj_f, v, params, step );
}


void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params );
    gsl_obj_df( v, params, grad );
}
