/*
 * PredictorTrainer.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: lunt
 */


#include "ExprPar.h"
#include "DataSet.h"


#include "PredictorTrainer.h"

#include "ExprPredictor.h"


ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    if ( toupperStr( objOptionStr ) == "PGP" ) return PGP;
    if ( toupperStr( objOptionStr ) == "LOGISTIC_REGRESSION") return LOGISTIC_REGRESSION;
    if ( toupperStr( objOptionStr ) == "PEAK_WEIGHTED") return PEAK_WEIGHTED;
    if ( toupperStr( objOptionStr ) == "WEIGHTED_SSE") return WEIGHTED_SSE;


    cerr << "objOptionStr is not a valid option of objective function" << endl;
    exit(1);
}


string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    if ( objOption == PGP ) return "PGP";
    if ( objOption == LOGISTIC_REGRESSION ) return "LOGISTIC_REGRESSION";
    if ( objOption == PEAK_WEIGHTED ) return "PEAK_WEIGHTED";
    if ( objOption == WEIGHTED_SSE ) return "WEIGHTED_SSE";

    return "Invalid";
}


string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}

ExprPar TrainingPipeline::train(const ExprPredictor* predictor, const TrainingDataset* training_data, const ExprPar& par_start ){
	ExprPar current = par_start;
	for(int i = 0;i<trainers.size();i++){
		current = trainers[i].get()->train(predictor, training_data, current);
	}
	return current;
}


double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data){
		ExprPredictor* predictor = (ExprPredictor*)f_data;
		predictor->begin_batch();

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
