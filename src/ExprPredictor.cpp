#include <typeinfo>


#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include <nlopt.hpp>

#include "ExprPredictor.h"
#include "ExprPar.h"
#include "ExprFunc.h"

double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data);

ExprPredictor::ExprPredictor( const vector <Sequence>& _seqs, const vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, TrainingDataset* _training_data, const vector< Motif >& _motifs, const ExprModel& _expr_model,
		const vector < bool >& _indicator_bool, const vector <string>& _motifNames) : TrainingAware(), seqs(_seqs), seqSites( _seqSites ), seqLengths( _seqLengths ), training_data( _training_data ),
	expr_model( _expr_model),
	indicator_bool ( _indicator_bool ), motifNames ( _motifNames ),
	search_option(UNCONSTRAINED)
{
    //TODO: Move appropriate lines from this block to the ExprModel class.
	cerr << "exprData size: " << training_data->n_rows_output() << "  " << nSeqs() << endl;
    assert( training_data->n_rows_output() == nSeqs() );
    //assert( training_data->factorExprData.nRows() == nFactors() && training_data->factorExprData.nCols() == nConds() );
    //assert( expr_model.coopMat.isSquare() && expr_model.coopMat.isSymmetric() && expr_model.coopMat.nRows() == nFactors() );
    assert( expr_model.actIndicators.size() == nFactors() );
    assert( expr_model.maxContact > 0 );
    assert( expr_model.repIndicators.size() == nFactors() );
    assert( expr_model.repressionMat.isSquare() && expr_model.repressionMat.nRows() == nFactors() );
    assert( expr_model.repressionDistThr >= 0 );

	//****** DEFAULT VALUES *********
	objOption = SSE;

	n_alternations = 4;
	n_random_starts = 5;


	max_simplex_iterations = 200;
	max_gradient_iterations = 50;


    //gene_crm_fout.open( "gene_crm_fout.txt" );

    // set the model option for ExprPar and ExprFunc
    //ExprPar::modelOption = expr_model.modelOption;//TODO: Remove both of these.
    //ExprFunc::modelOption = expr_model.modelOption;

    // set the values of the parameter range according to the model option
    if ( expr_model.modelOption != LOGISTIC && expr_model.modelOption != DIRECT )
    {
        //ExprPar::min_effect_Thermo = 0.99;
        //ExprPar::min_interaction = 0.99;
    }

    //expr_model was already initialized. Setup the parameter factory.
    param_factory = new ParFactory(expr_model, nSeqs());

	trainingObjective = NULL;
	set_objective_option(objOption);

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

void ExprPredictor::set_objective_option( ObjType in_obj_option ){
	//TODO: Move this to the front-end or something?
    //Maybe make it have a default SSE score objective, but anything else gets specified in the front-end.
	objOption = in_obj_option;

	if(NULL != trainingObjective){
		delete trainingObjective;
		trainingObjective = NULL;
	}

    switch(in_obj_option){
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

    if ( n_alternations > 0 && this->search_option == CONSTRAINED ){
      par_model = param_factory->truncateToBounds(par_model, indicator_bool);

    }
    obj_model = objFunc( par_model );

    cout << "*** Diagnostic printing AFTER adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( n_alternations == 0 ) return 0;

    // alternate between two different methods
    ExprPar par_result = param_factory->create_expr_par();
    double obj_result;
	this->start_training();

    for ( int i = 1; i <= n_alternations; i++ )
    {
		this->begin_epoch(i);
        simplex_minimize( par_result, obj_result );
        par_model = par_result;

        gradient_minimize( par_result, obj_result );
        par_model = par_result;
    }

    #ifdef BETAOPTBROKEN
    optimize_beta( par_model, obj_result );
    #endif

	this->end_training();

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
    for ( int i = 0; i < n_random_starts; i++ )
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
    if ( n_random_starts ) train( par_best );
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

/**
In training mode, will skip zero-weighted bins.
*/
int ExprPredictor::predict( const ExprPar& par, const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const
{
	// predict the expression


		//Code for skipping during training BEGIN_SKIPPING
		/*TODO: dynamic_cast is slow, maybe it would be better to move this code that decides
		which bins to predict out to some pre-epoch place so it only gets called once.
		For now, we value correctness above efficiency.
		*/
		Matrix *weights = NULL;
		if( NULL != dynamic_cast<const Weighted_ObjFunc_Mixin*>(this->trainingObjective) ){
			weights = ((Weighted_ObjFunc_Mixin*)this->trainingObjective)->get_weights();
		}
		//End of skipping code.	END_SKIPPING


    ExprFunc* func = createExprFunc( par , targetSites_, targetSeqLength, seq_num);
		targetExprs.resize(nConds());
    for ( int j = 0; j < nConds(); j++ )
    {
				//Code for skipping during training BEGIN_SKIPPING
				if( in_training && weights != NULL && weights->getElement(seq_num,j) <= 0.0){
					//cerr << "TEMPORARY DEBUG CODE, skipping unweighted bin (" << seq_num << "," << j << ")." << endl;
					targetExprs[j] = 0.0;
					continue;
				}
				//cerr << "Unskipped bin, making prediction." << endl;
				//End of skipping code. END_SKIPPING


				Condition concs = training_data->getCondition( j , par );
        double predicted = func->predictExpr( concs );
        targetExprs[j] = ( predicted );
    }

    delete func;
    return 0;
}

/**
While in training mode (private in_training variable == true), this method will skip the prediction of zero-weighted positions.
*/
int ExprPredictor::predict_all( const ExprPar& par , vector< vector< double > > &targetExprs ) const
{
	vector< int > seqLengths( seqs.size() );
	vector< SiteVec > seqSites( seqs.size() ); //
	targetExprs.clear();

    for( int i = 0; i < seqs.size(); i++ ){
      seqLengths[i] = seqs[i].size();
    }


    #ifdef REANNOTATE_EACH_PREDICTION
    SeqAnnotator ann( expr_model.motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ ) {
       	ann.annot( seqs[ i ], seqSites[ i ] );
    }
    #else
    seqSites = this->seqSites;
    #endif

    //Create predictions for every sequence and condition
    for ( int i = 0; i < nSeqs(); i++ ) {
			vector<double> one_seq_predictions(nConds());

			this->predict(par, seqSites[i], seqLengths[i], one_seq_predictions, i );

			targetExprs.push_back(one_seq_predictions);
    }

	return 0;
}

int ExprPredictor::maxShift = 5;
double ExprPredictor::shiftPenalty = 0.8;

double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
double ExprPredictor::min_delta_f_PGP = 1.0E-8;



void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 8 );
    //     cout.width( 8 );
	cout << par.my_pars;
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
	vector<vector<double> > ground_truths;
	vector<vector<double> > predictions;

	for(int i = 0;i< nSeqs();i++){//Populate ground truths
		ground_truths.push_back(training_data->get_output_row(i));
	}

	this->predict_all(par, predictions);

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
    optimizer.set_initial_step(1.0);//TODO: enforce simplex starting size.
	if(max_simplex_iterations > -1){ optimizer.set_maxeval(max_simplex_iterations); }

    if(this->search_option == CONSTRAINED){
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

    if(this->search_option == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }


    //TODO: enforce nGradientIters
	if(max_gradient_iterations > -1){ optimizer.set_maxeval(max_gradient_iterations); }

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
    double step = 1.0E-6;
    numeric_deriv( grad, gsl_obj_f, v, params, step );
}


void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params );
    gsl_obj_df( v, params, grad );
}
