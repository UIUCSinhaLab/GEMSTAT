#include <cstring>


#include "regularization.hpp"

RegularizedObjFunc::RegularizedObjFunc(ObjFunc* wrapped_obj_func, const ExprPar& centers, const ExprPar& l1, const ExprPar& l2)
{
  my_wrapped_obj_func = wrapped_obj_func;
  ExprPar tmp_energy_space = centers.my_factory->changeSpace(centers, ENERGY_SPACE);
  tmp_energy_space.getRawPars(my_centers );

  //It doesn't matter what space these are in, they are just storage for values.
  l1.getRawPars(lambda1 );
  l2.getRawPars(lambda2 );
  cache_pars = vector<double>(my_centers.size(),0.0);
  //cache_diffs(my_centers.size(),0.0);
  //cache_sq_diffs(my_centers.size(),0.0);
}

double RegularizedObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par){

  double objective_value = my_wrapped_obj_func->eval( ground_truth, prediction, par );



  double l1_running_total = 0.0;
  double l2_running_total = 0.0;

  ExprPar tmp_energy_space = par->my_factory->changeSpace(*par, ENERGY_SPACE);
  tmp_energy_space.getRawPars(cache_pars );

  for(int i = 0;i<cache_pars.size();i++){
    double the_diff = abs(cache_pars[i] - my_centers[i]);
    l1_running_total += lambda1[i]*the_diff;
    l2_running_total += lambda2[i]*pow(the_diff,2.0);
  }

  objective_value += l1_running_total + l2_running_total;

  return objective_value;
}


void regularization_cmndline_init(ExprPredictor *predictor, int argc, char* argv[] ){
    
    
    double l1 = 0.0;
    double l2 = 0.0;
    
    
    for ( int i = 1; i < argc; i++ )
    {
        if ( !strcmp("-l1", argv[ i ]))         { l1 = atof(argv[ ++i ]); }
        else if ( !strcmp("-l2", argv[ i ]))    { l2 = atof(argv[ ++i ]); }
    }
    
    //Setup regularization objective function
    ExprPar tmp_centers, tmp_l1, tmp_l2;
    bool setup_regularization = false;
    if(0.0 != l1 || 0.0 != l2){
        setup_regularization = true;
      cerr << "INFO: Regularization was turned on and will be used. l1 = " << l1 << " l2 = " << l2 << " ."<< endl;

      tmp_centers = predictor->param_factory->create_expr_par();
      tmp_l1 = predictor->param_factory->create_expr_par();
      tmp_l2 = predictor->param_factory->create_expr_par();

      //TODO: add an option to read l1 and l2 values from a file.
      vector< double > tmp_l12_vector;
      tmp_l1.getRawPars(tmp_l12_vector);
      std::fill(tmp_l12_vector.begin(),tmp_l12_vector.end(),l1);
      tmp_l1 = predictor->param_factory->create_expr_par(tmp_l12_vector, ENERGY_SPACE);

      tmp_l2.getRawPars(tmp_l12_vector);
      std::fill(tmp_l12_vector.begin(),tmp_l12_vector.end(),l2);
      tmp_l2 = predictor->param_factory->create_expr_par(tmp_l12_vector, ENERGY_SPACE);

      RegularizedObjFunc *tmp_reg_obj_func = new RegularizedObjFunc(predictor->trainingObjective,
                                              tmp_centers,
                                              tmp_l1,
                                              tmp_l2
                                            );
      predictor->trainingObjective = tmp_reg_obj_func;
    }
    
    //DONE setting up regularization. Do some checking next
    
    /*TODO: Break checking code into its own module.
    * THis code is copied from main. It is absolute garbage.
    */
    
    //**** CHECK that all loaded pars have the same ordering and parameter names.
    
    ExprPar par_init = predictor->param_factory->create_expr_par();//Should only be scoped to this function.
    
    vector<string> path_vector;
    for(gsparams::DictList::iterator itr = par_init.my_pars.begin();itr != par_init.my_pars.end();++itr){
        path_vector.push_back(itr.get_path());
    }
    vector< std::pair<std::string , gsparams::DictList* > > all_loaded_params;

    /*Only if using l1/l2 regularization*/
    if(setup_regularization){
        all_loaded_params.push_back(std::make_pair("reg_centers",&(tmp_centers.my_pars)));
        all_loaded_params.push_back(std::make_pair("l1_weights",&(tmp_l1.my_pars)));
        all_loaded_params.push_back(std::make_pair("l2_weights",&(tmp_l2.my_pars)));
    }

    for(int i = 0;i<all_loaded_params.size();i++){
        gsparams::DictList::iterator itr = all_loaded_params[i].second->begin();
        int j = 0;
        while(itr!=all_loaded_params[i].second->end()){
            if(0 != path_vector[j].compare(itr.get_path())){
                
                throw std::invalid_argument("Error in one of the input parameter files (.par, free_fix, lower/upper bounds, etc.) \nTrying to compare SNOT objects, it seems that one of the inputs is misordered.  starting_parameters" + path_vector[j] + " is not " + all_loaded_params[i].first + itr.get_path() + " ");
            }
            ++itr;
            j++;
        }
    }
    
    
}
