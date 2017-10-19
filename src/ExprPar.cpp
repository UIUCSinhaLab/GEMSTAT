#include "ExprPar.h"
#include "ExprPredictor.h"
#include "conf/ExprParConf.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>


string parameterSpaceStr(ThermodynamicParameterSpace in){
    if(in == CONSTRAINED_SPACE)
      return "CONSTRAINED";
    if(in == ENERGY_SPACE)
      return "ENERGY_SPACE";
    if(in == PROB_SPACE)
      return "PROB_SPACE";
    assert(false);
}


ParFactory::ParFactory( const ExprModel& in_model, int in_nSeqs) : expr_model(in_model), nSeqs(in_nSeqs)
{
  //maximums = ExprPar( expr_model.motifs.size(), nSeqs );
  maximums = createDefaultMinMax(true);
  //minimums = ExprPar( expr_model.motifs.size(), nSeqs );
  minimums = createDefaultMinMax(false);

  defaults = create_expr_par();
}

int ParFactory::nFactors() const
{
  return expr_model.motifs.size();
}

void ParFactory::separateParams(const ExprPar& input, vector<double>& free_output, vector<double>& fixed_output, const vector<bool>& indicator_bool) const
{
  //Hassan's code for separating parameters
  vector<double> pars;

  input.getRawPars( pars );
  int pars_size = pars.size();
  free_output.clear();
  fixed_output.clear();
  for( int index = 0; index < pars_size; index++ )
  {
      if( indicator_bool[ index ] )
      {
          free_output.push_back( pars[ index ]);
      }
      else
      {
          fixed_output.push_back( pars[ index ] );
      }
  }
}

void ParFactory::joinParams(const vector<double>& free_pars, const vector<double>& fix_pars, vector<double>& output, const vector<bool>& indicator_bool) const
{
  output.clear();
  int pars_size = free_pars.size() + fix_pars.size();
  //assert(indicator_bool.size() == pars_size)
  int free_par_counter = 0;
  int fix_par_counter = 0;
  for( int index = 0; index < pars_size; index ++ )
  {
      if( indicator_bool[ index ] )
      {
          output.push_back( free_pars[ free_par_counter ++ ]);
      }
      else
      {
          output.push_back( fix_pars[ fix_par_counter ++ ]);
      }
  }

}

ExprPar ParFactory::truncateToBounds(const ExprPar& in_par, const vector<bool>& indicator_bool) const
{
  //giving an empty indicator bool shall mean that all parameters are subject to truncation
  //otherwise, only those for which indicator_bool is true shall be.
  vector<bool> use_indicator_bool = indicator_bool;
  if(use_indicator_bool.empty()){
    vector<double> vector_tmpfoo;
    in_par.getRawPars(vector_tmpfoo);
    use_indicator_bool = vector<bool>(true,vector_tmpfoo.size());
  }

  ThermodynamicParameterSpace original_space = in_par.my_space;

  vector< double > adjust_free_pars;
  vector< double > adjust_fix_pars;
  vector< double > post_adjust_fix_pars;//Used for penultimate output too

  ExprPar tmp_model_es = this->changeSpace(in_par, ENERGY_SPACE);//Be sure we are in the right parameter space
  this->separateParams(tmp_model_es, adjust_free_pars, adjust_fix_pars, use_indicator_bool );

  ExprPar tmp_constrained_es = this->changeSpace(
                                  this->changeSpace(tmp_model_es,CONSTRAINED_SPACE),
                                  ENERGY_SPACE);
  this->separateParams(tmp_constrained_es, adjust_free_pars, post_adjust_fix_pars, use_indicator_bool);
  this->joinParams(adjust_free_pars, adjust_fix_pars, post_adjust_fix_pars, use_indicator_bool);

  ExprPar final_return = this->create_expr_par(post_adjust_fix_pars, ENERGY_SPACE);
  final_return = this->changeSpace(final_return, original_space);
  return final_return;
}

bool ParFactory::testWithinBounds(const ExprPar& in_par) const
{
  vector< double > mins;
  vector< double > maxes;
  vector< double > totest;


  this->minimums.getRawPars(mins);
  this->maximums.getRawPars(maxes);

  assert(this->minimums.my_space == this->maximums.my_space);
  ExprPar to_test_es = this->changeSpace(in_par,this->minimums.my_space);//should be ENERGY_SPACE
  to_test_es.getRawPars(totest);

  int iter_end = totest.size();
  for(int i=0;i<iter_end;i++){
    if(totest[i] < mins[i] || totest[i] > maxes[i]){return false;}
  }
  return true;
}

ExprPar ParFactory::createDefaultMinMax(bool min_or_max) const
{
  //TODO: the model really should know about the number of factors in the model.
  int _nFactors = expr_model.motifs.size();
  assert( _nFactors > 0 );

  ExprPar tmp_par = create_expr_par();

  tmp_par.maxBindingWts.assign(tmp_par.maxBindingWts.size(), min_or_max ? log(ExprPar::max_weight) : log(ExprPar::min_weight)); // ExprPar::min_weight
  //set the interaction maximums
  tmp_par.factorIntMat.setAll(min_or_max ? log(ExprPar::max_interaction) : log(ExprPar::min_interaction));
  tmp_par.txpEffects.assign(tmp_par.txpEffects.size(), min_or_max ? log(ExprPar::max_effect_Thermo) : log(ExprPar::min_effect_Thermo));//TODO: Doesn't handle Logistic model
  tmp_par.repEffects.assign(tmp_par.repEffects.size(), min_or_max ? log(ExprPar::max_repression) : log(ExprPar::min_repression));
  tmp_par.basalTxps.assign(tmp_par.basalTxps.size(), min_or_max ? log(ExprPar::max_basal_Thermo) : log(ExprPar::min_basal_Thermo));//TODO: Doesn't handle Logistic model
  tmp_par.pis.assign(tmp_par.pis.size(), min_or_max ? log(ExprPar::max_pi) : log(ExprPar::min_pi));
  tmp_par.betas.assign(tmp_par.betas.size(), min_or_max ? log(ExprPar::max_beta) : log(ExprPar::min_beta));
  tmp_par.energyThrFactors.assign(tmp_par.energyThrFactors.size(), min_or_max ? log(ExprPar::max_energyThrFactors) : log(ExprPar::min_energyThrFactors));
  tmp_par.my_space = ENERGY_SPACE;
  return tmp_par;
}

ExprPar ParFactory::create_expr_par() const
{
  ExprPar tmp_par = ExprPar(expr_model.getNFactors(), this->nSeqs);

  int _nFactors = expr_model.motifs.size();
  assert( _nFactors > 0 );

  tmp_par.maxBindingWts.assign( _nFactors, log( ExprPar::default_weight ) );
  tmp_par.factorIntMat.setDimensions( _nFactors, _nFactors );
  tmp_par.factorIntMat.setAll( log( ExprPar::default_interaction ) );

  double defaultEffect = expr_model.modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : log(ExprPar::default_effect_Thermo);
  tmp_par.txpEffects.assign( _nFactors, defaultEffect );
  tmp_par.repEffects.assign( _nFactors, log( ExprPar::default_repression ) );

  int numBTMS = expr_model.one_qbtm_per_crm ? nSeqs : 1;

  double basalTxp_val = expr_model.modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : log( ExprPar::default_basal_Thermo );
  tmp_par.basalTxps.assign( numBTMS, basalTxp_val );
  //for the pausing parameters
  tmp_par.pis.assign( expr_model.shared_scaling ? 1 : nSeqs, log( ExprPar::default_pi ) );

  //for the beta parameters
  tmp_par.betas.assign( expr_model.shared_scaling ? 1 : nSeqs, log( ExprPar::default_beta ) );

  tmp_par.energyThrFactors.assign( _nFactors, log( ExprPar::default_energyThrFactors ) );
  tmp_par.my_space = ENERGY_SPACE;
  tmp_par.my_factory = this;
  return tmp_par;
}

ExprPar ParFactory::create_expr_par(const vector<double>& pars, const ThermodynamicParameterSpace in_space) const
{
      //TODO: Most of this code can be further simplified by using the STL vector's nice functions.

      ExprPar tmp_par = this->create_expr_par();

      int _nFactors = expr_model.getNFactors();

      int counter = 0;

      // set maxBindingWts
      tmp_par.maxBindingWts.clear();
      for ( int i = 0; i < _nFactors; i++ )
      {
          double weight = pars[counter++] ;
          tmp_par.maxBindingWts.push_back( weight );
      }

      // set the interaction matrix
      //Which was created previously.
      double scaled_default = ExprPar::default_interaction;
      if(in_space == ENERGY_SPACE)
        scaled_default = log(ExprPar::default_interaction);
      if(in_space == CONSTRAINED_SPACE)
        scaled_default = infty_transform(log(ExprPar::default_interaction), log( ExprPar::min_interaction ), log( ExprPar::max_interaction ));
      tmp_par.factorIntMat.setAll(scaled_default);


      for ( int i = 0; i < _nFactors; i++ )
      {
          for ( int j = 0; j <= i; j++ )
          {
              if ( expr_model.coop_setup->has_coop( i, j ) )
              {
                  double interaction = pars[counter++] ;
                  tmp_par.factorIntMat( i, j ) = interaction;
                  tmp_par.factorIntMat( j, i ) = interaction;
              }
          }
      }

      //TODO: For the logistic model, actIndicators should always be 1 and actRepressors always 0

      // set the transcriptional effects
      scaled_default = ExprPar::default_effect_Thermo;
      if(in_space == ENERGY_SPACE)
        scaled_default = log(ExprPar::default_effect_Thermo);
      if(in_space == CONSTRAINED_SPACE)
        scaled_default = infty_transform(log(ExprPar::default_effect_Thermo), log( ExprPar::min_effect_Thermo ), log( ExprPar::max_effect_Thermo ));

      for ( int i = 0; i < _nFactors; i++ )
      {
              double effect;

              if ( expr_model.actIndicators[i] )
                effect = pars[counter++];
              else
                effect = scaled_default;//TODO: Which space are we creating in?

              tmp_par.txpEffects[i] = effect ;
      }

      // set the repression effects
      scaled_default = ExprPar::default_repression;
      if(in_space == ENERGY_SPACE)
        scaled_default = log(ExprPar::default_repression);
      if(in_space == CONSTRAINED_SPACE)
        scaled_default = infty_transform(log(ExprPar::default_repression), log( ExprPar::min_repression ), log( ExprPar::max_repression ));

      tmp_par.repEffects.assign(_nFactors, scaled_default);//TODO: Which space are we creating in?
      if ( expr_model.modelOption == CHRMOD_UNLIMITED || expr_model.modelOption == CHRMOD_LIMITED || expr_model.modelOption == DIRECT || expr_model.modelOption == MARKOV )
      {
          for ( int i = 0; i < _nFactors; i++ )
          {
              if ( expr_model.repIndicators[i] )
              {
                  double repression = pars[counter++];
                  tmp_par.repEffects[i] = repression;
              }
          }
      }

      // set the basal transcription

      int num_qbtm = expr_model.one_qbtm_per_crm ? nSeqs : 1;
      tmp_par.basalTxps.clear();
      for( int i = 0; i < num_qbtm; i++ )
      {
            double basal = pars[counter++];
            tmp_par.basalTxps.push_back( basal );
      }

      //for( int i = 0; i < num_qbtm; i++ )//TODO: should be this but that triggers an exception somewhere I don7t want to deal with.
      for( int i = 0; i < (expr_model.shared_scaling ? 1 : nSeqs); i++ )
      {
          //cout << "sending for:" << pars[counter] <<"\t" << log(min_pi) << "\t" <<log(max_pi) << endl;
          double pi_val = pars[ counter++ ] ;
          tmp_par.pis[i] = pi_val ;
      }

      //cout << "Counter before beta: " << counter << endl;
      //cout << "nSeqs: " << nSeqs << "\t _nSeqs: " << _nSeqs << endl;
      //put in the values of the beta params

      //for( int i = 0; i < num_qbtm; i++ )//TODO: should be this but that triggers an exception somewhere I don7t want to deal with.

      for( int i = 0; i < (expr_model.shared_scaling ? 1 : nSeqs); i++ )
      {
          double beta_val = pars[ counter ++];
          tmp_par.betas[i] = beta_val ;
      }

      //cout << "Counter after beta: " << counter << endl;
      for( int i = 0; i < _nFactors; i++ )
      {
          double energyThrFactor_val = pars[ counter++ ];
          tmp_par.energyThrFactors[i] = energyThrFactor_val ;
      }
      assert(counter == pars.size());

      tmp_par.my_space = in_space;
      return tmp_par;
}



ExprPar ParFactory::changeSpace(const ExprPar& in_par, const ThermodynamicParameterSpace new_space) const
{
  if( in_par.my_space == new_space){
    return in_par;
  }

  //Since the spaces are related by CONSTRAINED_SPACE <-> ENERGY_SPACE <-> PROB_SPACE
  //and we already handled a transformation to the same space type above, we know we at least need to go to the energy space.

  vector<double> as_energy_space;
  vector<double> original_pars;


  vector<double> low_vect;//TODO: Cache this when the high and low exprPar are set.
  vector<double> high_vect;
  if(in_par.my_space == CONSTRAINED_SPACE || new_space == CONSTRAINED_SPACE){
        maximums.getRawPars(high_vect );
        minimums.getRawPars(low_vect );
  }


  in_par.getRawPars(original_pars );

  if(in_par.my_space == CONSTRAINED_SPACE){
    constrained_to_energy_helper(original_pars, as_energy_space,low_vect, high_vect);
  }else if(in_par.my_space == PROB_SPACE){//Inputs were in the PROB_SPACE
    prob_to_energy_helper(original_pars, as_energy_space);
  }else{ // Intputs were in ENERGY_SPACE
    as_energy_space = original_pars;
  }

  //Now convert that into whichever the target parameter space was.
  if(new_space == ENERGY_SPACE){
    return create_expr_par(as_energy_space,ENERGY_SPACE);
  }

  if(new_space == PROB_SPACE){
    vector<double> as_prob_space;
    energy_to_prob_helper(as_energy_space,as_prob_space);
    return create_expr_par(as_prob_space,PROB_SPACE);
  }

  if(new_space == CONSTRAINED_SPACE){
    vector<double> as_constrained;
    energy_to_constrained_helper(as_energy_space, as_constrained, low_vect, high_vect);
    return create_expr_par(as_constrained, CONSTRAINED_SPACE);
  }

  assert(false);
}

void ParFactory::constrained_to_energy_helper(const vector<double>& pars, vector<double>& output, const vector<double>& low, const vector<double>& high) const
{
  output.clear();
  assert(pars.size() == low.size() && pars.size() == high.size());

  int N = pars.size();
  for(int i = 0; i<N;i++){
    double value = inverse_infty_transform( pars[i], low[i], high[i] );
    output.push_back(value);
  }


}

void ParFactory::energy_to_constrained_helper(const vector<double>& pars, vector<double>& output, const vector<double>& low, const vector<double>& high) const
{
  output.clear();
  assert(pars.size() == low.size() && pars.size() == high.size());

  int N = pars.size();
  for(int i = 0; i<N;i++){
    double value = infty_transform( pars[i], low[i], high[i] );
    output.push_back(value);
  }

}

void ParFactory::energy_to_prob_helper(const vector<double>& pars, vector<double>& output) const
{
  output.assign(pars.size(), 1.0);

  int N = pars.size();
  for(int i = 0; i<N;i++){
    output[i] = exp(pars[i]);
  }
}

void ParFactory::prob_to_energy_helper(const vector<double>& pars, vector<double>& output) const
{
  output.assign(pars.size(), 0.0);

  int N = pars.size();
  for(int i = 0; i<N;i++){
    output[i] = log(pars[i]);
  }
}

ExprPar ParFactory::randSamplePar( const gsl_rng* rng) const
{
    ExprPar tmp_par = create_expr_par();
    tmp_par.my_space = ENERGY_SPACE;
    // sample binding weights
    //estBindingOption is now always true, fix the bindings to 1.0 if you don't want to estimate them.
    for ( int i = 0; i < nFactors(); i++ )
    {
      tmp_par.maxBindingWts[i] = gsl_ran_flat( rng, minimums.maxBindingWts[i], maximums.maxBindingWts[i] );
    }


    // sample the interaction matrix
    for ( int i = 0; i < nFactors(); i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
            if ( expr_model.coop_setup->has_coop( i, j ) ){
              double rand_interaction = gsl_ran_flat( rng, minimums.factorIntMat(i, j), maximums.factorIntMat(i,j) );
              tmp_par.factorIntMat( i, j ) = rand_interaction;
              tmp_par.factorIntMat( j, i ) = rand_interaction;
           }
        }
    }

    // sample the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
      if( expr_model.actIndicators[i]){//TODO: Consider removing
            double rand_effect = gsl_ran_flat( rng, minimums.txpEffects[i], maximums.txpEffects[i] );
            tmp_par.txpEffects[i] = rand_effect;
      }
    }

    // sample the repression effects
    for ( int i = 0; i < nFactors(); i++ )
    {
      if ( expr_model.repIndicators[i] )//TODO: Consider removing
      {
          double rand_repression = gsl_ran_flat( rng, minimums.repEffects[i], maximums.repEffects[i] );
          tmp_par.repEffects[i] = rand_repression;
      }
    }

    // sample the basal transcription
    int num_qbtm = expr_model.one_qbtm_per_crm ? nSeqs : 1;
    double rand_basal;
    for( int i = 0; i < num_qbtm; i ++ )
    {
        rand_basal = gsl_ran_flat( rng, minimums.basalTxps[i], maximums.basalTxps[i] );
        tmp_par.basalTxps[ i ] = rand_basal;
    }


    return tmp_par;
}


ExprPar ParFactory::load(const string& file){
  ExprPar ret_par;
  // open the file
  ifstream fin( file.c_str() );
  if ( !fin ){ cerr << "Cannot open parameter file " << file << endl; exit( 1 ); }

  std::string header;
  std::getline(fin,header);
  fin.seekg(0,fin.beg);//Every loader can expect to be at the beginning of the file.

  if(0 == header.compare("#GSPAR1.6a")) {//TODO: remove any whitespace off the end of header.
    ret_par = load_1_6a(fin);
  }else{
    ret_par = load_old(fin);
  }
  fin.close();

  return ret_par;
}

ExprPar ParFactory::load_1_6a(istream& fin){
    ExprPar tmp_par = create_expr_par();
    tmp_par = changeSpace(tmp_par, expr_model.modelOption == LOGISTIC ? ENERGY_SPACE : PROB_SPACE );//TODO: get rid of this so that logistic models are stored in the same space with the other models.

    vector< string > motifNames;


    std::string header;
    std::getline(fin,header);
    assert(0 == header.compare("#GSPAR1.6a"));

    cerr << "LOADING 1.6a parfile" << endl;
    //https://stackoverflow.com/questions/7868936/read-file-line-by-line
    //https://stackoverflow.com/questions/236129/split-a-string-in-c
    std:string line;
    std::istringstream line_ss;
    vector<string> tokens;
    #define LOCAL_TOKENIZE(M_TOK_VECT,M_LINE_STR,M_SS) M_TOK_VECT.clear();\
                            M_SS.clear();\
                            M_SS.str(M_LINE_STR);\
                            copy(istream_iterator<string>(M_SS),\
                            istream_iterator<string>(),\
                            back_inserter(M_TOK_VECT))

    //Read the TF lines.
    int tf_i = 0;
    for(tf_i = 0; tf_i < expr_model.getNFactors() && std::getline(fin,line);tf_i++){
        LOCAL_TOKENIZE(tokens,line,line_ss);

        //Read a TF line.
        assert(tokens.size() == 3 || tokens.size() == 4);
        motifNames.push_back(string(tokens[0]));
        tmp_par.maxBindingWts[tf_i] = atof(tokens[1].c_str());
        tmp_par.txpEffects[tf_i] = atof(tokens[2].c_str());
        if(tokens.size() >= 4){
          tmp_par.repEffects[tf_i] = atof(tokens[3].c_str());
        }

    }

    //Expect the basal transcription line
    std::getline(fin,line);
    LOCAL_TOKENIZE(tokens,line,line_ss);

    vector<double> reading_basal;
    reading_basal.clear();
    for(tf_i=2;tf_i<tokens.size();tf_i++){
      reading_basal.push_back(atof(tokens[tf_i].c_str()));
    }

    assert(tokens.size() >= 3);
    assert(0 == tokens[0].compare("basal_transcription") &&
            0 == tokens[1].compare("="));
    //assert(tokens.size() == 2 + tmp_par.basalTxps.size());
    //If there is one provided, but we are using multiple qbtm_per (or one and one)
    if( tmp_par.basalTxps.size() == reading_basal.size()){
      std::copy(reading_basal.begin(),reading_basal.end(),tmp_par.basalTxps.begin());
    }else{//Mismatch size
      double avg = 0.0;
      for(tf_i = 0;tf_i<reading_basal.size();tf_i++){
        avg += reading_basal[tf_i];
      }
      avg /= (reading_basal.size());
      tmp_par.basalTxps.assign(tmp_par.basalTxps.size(),avg);
    }

    //Done reading the basal transcription
    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < expr_model.getNFactors(); i++ )
    {
        factorIdxMap[motifNames[i]] = i;
    }

    //PIS
    std::getline(fin,line);
    LOCAL_TOKENIZE(tokens,line,line_ss);
    vector< double > the_pis;
    the_pis.clear();
    for(tf_i=0;tf_i<tokens.size();tf_i++){
      the_pis.push_back(atof(tokens[tf_i].c_str()));
    }
    if(tmp_par.pis.size() == the_pis.size()){
      std::copy(the_pis.begin(),the_pis.end(),tmp_par.pis.begin());
    }else{//Mismatch size
      double avg = 0.0;
      for(tf_i = 0;tf_i<the_pis.size();tf_i++){
        avg += the_pis[tf_i];
      }
      avg /= (the_pis.size());
      tmp_par.pis.assign(tmp_par.pis.size(),avg);
    }

    //BETA
    std::getline(fin,line);
    LOCAL_TOKENIZE(tokens,line,line_ss);
    vector< double > the_betas;
    the_betas.clear();
    for(tf_i=0;tf_i<tokens.size();tf_i++){
      the_betas.push_back(atof(tokens[tf_i].c_str()));
    }
    if(tmp_par.betas.size() == the_betas.size()){
      std::copy(the_betas.begin(),the_betas.end(),tmp_par.betas.begin());
    }else{//Mismatch size
      double avg = 0.0;
      for(tf_i = 0;tf_i<the_betas.size();tf_i++){
        avg += the_betas[tf_i];
      }
      avg /= (the_betas.size());
      tmp_par.betas.assign(tmp_par.betas.size(),avg);
    }

  //TODO: load Cooperativities
  // read the cooperative interactions
  for( int i = 0; i < expr_model.getNumCoop(); i++ )
  {
      double coopVal;
      string factor1, factor2;

      std::getline(fin,line);
      LOCAL_TOKENIZE(tokens,line,line_ss);

      factor1 = tokens[0];
      factor2 = tokens[1];
      coopVal = atof(tokens[2].c_str());

      if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) throw RET_ERROR;
      int idx1 = factorIdxMap[factor1];
      int idx2 = factorIdxMap[factor2];
      tmp_par.factorIntMat( idx1, idx2 ) = coopVal;
      tmp_par.factorIntMat( idx2, idx1 ) = coopVal;
      //cout << factor1 << "\t" << factor2 << "\t" << idx1 << "\t" << idx2 << endl;
  }

  //TODO: load annotation thresholds
  double factor_thr_val;
  tmp_par.energyThrFactors.clear();
  std::getline(fin,line);
  LOCAL_TOKENIZE(tokens,line,line_ss);
  for(int i = 0; i < tokens.size();i++){
      tmp_par.energyThrFactors.push_back( atof(tokens[i].c_str()) );
  }
    //TODO: write output for this format.


    //tmp_par.print(cerr,motifNames,expr_model.coopMat);
    //cerr << "This feature is not yet fully implemented. Sorry." << endl;
    //exit(1);

    #undef LOCAL_TOKENIZE

    return tmp_par;
}

ExprPar ParFactory::load_old(istream& fin){

  ExprPar tmp_par = create_expr_par();
  tmp_par = changeSpace(tmp_par, expr_model.modelOption == LOGISTIC ? ENERGY_SPACE : PROB_SPACE );//TODO: get rid of this so that logistic models are stored in the same space with the other models.



  // read the factor information
  vector< string > motifNames( expr_model.getNFactors() );
  for ( int i = 0; i < expr_model.getNFactors(); i++ )
  {
      fin >> motifNames[i] >> tmp_par.maxBindingWts[i] >> tmp_par.txpEffects[i];
      if ( expr_model.modelOption == CHRMOD_UNLIMITED || expr_model.modelOption == CHRMOD_LIMITED || expr_model.modelOption == DIRECT || expr_model.modelOption == MARKOV) fin >> tmp_par.repEffects[i];
  }

  // factor name to index mapping
  map< string, int > factorIdxMap;
  for ( int i = 0; i < expr_model.getNFactors(); i++ )
  {
      factorIdxMap[motifNames[i]] = i;
  }

  // read the basal transcription
  string symbol, eqSign, value;
  fin >> symbol >> eqSign >> value;
  if ( symbol != "basal_transcription" || eqSign != "=" ) throw RET_ERROR;
  double basalTxp_val = atof( value.c_str() );
  tmp_par.basalTxps[ 0 ] =  basalTxp_val ;
  if( expr_model.one_qbtm_per_crm )
  {
      for( int _i = 1; _i < nSeqs; _i++ )
      {
          fin >> value;
          double basalTxp_val = atof( value.c_str() );
          tmp_par.basalTxps[ _i ] = basalTxp_val;
      }
  }
  //read the pi values
  for( int i = 0; i < (expr_model.shared_scaling ? 1 : nSeqs); i++ )
  {
      fin >> value;
      tmp_par.pis[ i ] = atof( value.c_str() );
  }

  //read the beta values
  for( int i = 0; i < (expr_model.shared_scaling ? 1 : nSeqs); i++ )
  {
      fin >> value;
      tmp_par.betas[ i ] = atof( value.c_str() );
  }

  // read the cooperative interactions
  string factor1, factor2;
  double coopVal;
  for( int i = 0; i < expr_model.getNumCoop(); i++ )
  {
      fin >> factor1 >> factor2 >> coopVal;
      if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) throw RET_ERROR;
      int idx1 = factorIdxMap[factor1];
      int idx2 = factorIdxMap[factor2];
      tmp_par.factorIntMat( idx1, idx2 ) = coopVal;
      tmp_par.factorIntMat( idx2, idx1 ) = coopVal;
      //cout << factor1 << "\t" << factor2 << "\t" << idx1 << "\t" << idx2 << endl;
  }
  double factor_thr_val;
  tmp_par.energyThrFactors.clear();
  while( fin >> factor_thr_val )
  {
      tmp_par.energyThrFactors.push_back( factor_thr_val );
  }

  return tmp_par;
}


ExprPar::ExprPar( int _nFactors, int _nSeqs ) : factorIntMat()
{
    assert( _nFactors > 0 );

    maxBindingWts.assign( _nFactors, ExprPar::default_weight );

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );

    double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
    txpEffects.assign( _nFactors, defaultEffect );
    repEffects.assign( _nFactors, ExprPar::default_repression );

    nSeqs = _nSeqs;

    int numBTMS = one_qbtm_per_crm ? nSeqs : 1;

    double basalTxp_val = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo;
    basalTxps.assign( numBTMS, basalTxp_val );
    //for the pausing parameters
    pis.assign( nSeqs, ExprPar::default_pi );

    //for the beta parameters
    betas.assign( nSeqs, ExprPar::default_beta );

    energyThrFactors.assign( _nFactors, ExprPar::default_energyThrFactors );

    my_space = modelOption == LOGISTIC ? ENERGY_SPACE : PROB_SPACE;
}


ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector <double>& _energyThrFactors ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxps( _basalTxps ), pis( _pis), betas( _betas ), nSeqs( _nSeqs ), energyThrFactors( _energyThrFactors )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() );
    assert( txpEffects.size() == maxBindingWts.size() && repEffects.size() == maxBindingWts.size() );
    assert( basalTxps.size() == one_qbtm_per_crm ? nSeqs : 1);
}


ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs ) : factorIntMat()
{
    assert(false);
    int _nFactors = actIndicators.size();
    nSeqs = _nSeqs;
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
    //     assert( pars.size() == ( _nFactors * ( _nFactors + 1 ) / 2 + 2 * _nFactors + 2 );
    int counter = 0;

    for ( int i = 0; i < _nFactors; i++ )
    {
        double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) : exp( pars[counter++] );
        maxBindingWts.push_back( weight );
    }

    // set the interaction matrix
    factorIntMat.setDimensions( _nFactors, _nFactors );
    for ( int i = 0; i < _nFactors; i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
            if ( coopMat( i, j ) )
            {
                double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );
                factorIntMat( i, j ) = interaction;
            }
            else factorIntMat( i, j ) = ExprPar::default_interaction;
        }
    }
    for ( int i = 0; i < _nFactors; i++ )
    {
        for ( int j = i + 1; j < _nFactors; j++ )
        {
            factorIntMat( i, j ) = factorIntMat( j, i );
        }
    }

    // set the transcriptional effects
    for ( int i = 0; i < _nFactors; i++ )
    {
        //         double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        if ( modelOption == LOGISTIC )
        {
            double effect = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_effect_Logistic, max_effect_Logistic ) : pars[counter++];
            txpEffects.push_back( effect );
        }
        /* else if ( modelOption == DIRECT ) {
                    double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : exp( pars[counter++] );
                    txpEffects.push_back( effect );
                }*/
        else
        {
            if ( actIndicators[i] )
            {
                double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : exp( pars[counter++] );
                txpEffects.push_back( effect );
            }
            else
            {
                txpEffects.push_back( ExprPar::default_effect_Thermo );
            }
        }
    }

    // set the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT || modelOption == MARKOV )
    {
        for ( int i = 0; i < _nFactors; i++ )
        {
            if ( repIndicators[i] )
            {
                double repression = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_repression ), log( max_repression ) ) ) : exp( pars[counter++] );
                repEffects.push_back( repression );
            }
            else
            {
                repEffects.push_back( ExprPar::default_repression );
            }
        }
    }
    else
    {
        for ( int i = 0; i < _nFactors; i++ ) repEffects.push_back( ExprPar::default_repression );
    }

    // set the basal transcription
    if( one_qbtm_per_crm )
    {
        nSeqs = _nSeqs;
        for( int i = 0; i < nSeqs; i++ )
        {
            if ( modelOption == LOGISTIC )
            {
                double basal = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_basal_Logistic, max_basal_Logistic ) : pars[counter++];
                basalTxps.push_back( basal );
            }
            else
            {
                double basal = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_basal_Thermo ), log( max_basal_Thermo ) ) ) : exp( pars[counter++] );
                basalTxps.push_back( basal );
            }
        }
    }
    else
    {

        if ( modelOption == LOGISTIC )
        {
            double basal = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_basal_Logistic, max_basal_Logistic ) : pars[counter++];
            basalTxps.push_back( basal );
        }
        else
        {
            double basal = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_basal_Thermo ), log( max_basal_Thermo ) ) ) : exp( pars[counter++] );
            basalTxps.push_back( basal );
        }
    }

    for( int i = 0; i < nSeqs; i++ )
    {
        //cout << "sending for:" << pars[counter] <<"\t" << log(min_pi) << "\t" <<log(max_pi) << endl;
        double pi_val = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[ counter++ ], log( min_pi ), log( max_pi ) ) ):exp( pars[ counter++ ] ) ;
        pis.push_back( pi_val );
    }
    //cout << "Counter before beta: " << counter << endl;
    //cout << "nSeqs: " << nSeqs << "\t _nSeqs: " << _nSeqs << endl;
    //put in the values of the beta params
    for( int i = 0; i < nSeqs; i++ )
    {
        double beta_val = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[ counter ++], log(min_beta), log( max_beta) ) ) : exp( pars[ counter++ ] );
        betas.push_back( beta_val );
    }
    //cout << "Counter after beta: " << counter << endl;
    for( int i = 0; i < _nFactors; i++ )
    {
        double energyThrFactor_val = searchOption == CONSTRAINED ? exp ( inverse_infty_transform( pars[ counter++ ], log( min_energyThrFactors ), log( max_energyThrFactors ) ) ) : exp ( pars[ counter++ ] );
        energyThrFactors.push_back( energyThrFactor_val );
    }
}



void ExprPar::getRawPars( vector< double >& pars) const
{
    IntMatrix& coopMat = this->my_factory->expr_model.coop_setup->coop_matrix;//TODO: I bet this results in a bug.
    vector< bool >& actIndicators = this->my_factory->expr_model.actIndicators;
    vector< bool >& repIndicators = this->my_factory->expr_model.repIndicators;

    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();

    // write maxBindingWts
    for ( int i = 0; i < nFactors(); i++ )
    {
      pars.push_back( maxBindingWts[ i ] );
    }

    // write the interaction matrix
    if ( modelOption != LOGISTIC )//TODO: just make the coopMat all false for LOGISTIC models.
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            for ( int j = 0; j <= i; j++ )
            {
                if ( coopMat( i, j ) )
                {
                      pars.push_back( factorIntMat( i, j ) );
                }
            }
        }
    }

    // write the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( modelOption == LOGISTIC )//TODO: just make actIndicators all true for LOGISTIC model.
        {
          pars.push_back( txpEffects[i] );
        }
        /*else if ( modelOption == DIRECT ) {
                    double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
                    pars.push_back( effect );
                }*/
        else
        {
            if ( actIndicators[i] ) //TODO: make actIndicators all true for DIRECT model
            {
                pars.push_back( txpEffects[i] );
            }
        }
    }

    // write the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT || modelOption == MARKOV)
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( repIndicators[i] )
            {
              pars.push_back( repEffects[i] );
            }
        }
    }

    for( int i = 0; i < basalTxps.size(); i++ )
    {
        pars.push_back( basalTxps[ i ] );
    }

    //write the pis
    for( int i = 0; i < pis.size(); i++ )
    {
        pars.push_back( pis[ i ] );
    }

    //write the betas
    for( int i = 0; i < betas.size(); i++ )
    {
        pars.push_back( betas[ i ] );
    }
    for( int i = 0; i < energyThrFactors.size(); i++ )
    {
        pars.push_back( energyThrFactors[ i ] );
    }
}

GEMSTAT_PAR_FLOAT_T ExprPar::getBetaForSeq(int enhancer_ID) const {
    int use_enhancerID = (this->my_factory->expr_model.shared_scaling ? 0 : enhancer_ID);
    return betas[ use_enhancerID ];
}

GEMSTAT_PROMOTER_DATA_T ExprPar::getPromoterData(int enhancer_ID) const {
    GEMSTAT_PROMOTER_DATA_T the_return_value;

    //TODO: we can create a more complicated mapping later.
    int use_enhancerID = (this->my_factory->expr_model.shared_scaling ? 0 : enhancer_ID);
    int use_basal = (this->my_factory->expr_model.one_qbtm_per_crm ? use_enhancerID : 0);
    the_return_value.basal_trans = basalTxps[ use_basal ];
    the_return_value.pi = pis[ use_enhancerID ];
    the_return_value.beta = betas[ use_enhancerID ];

    return the_return_value;
}

void ExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const
{
//    os.setf( ios::fixed );
//    os.precision( 50 );
    os << "#GSPAR1.6a" << endl; //par file version header.
    // print the factor information
    for ( int i = 0; i < nFactors(); i++ )
    {
        os << motifNames[i] << "\t" << maxBindingWts[i] << "\t" << txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT || modelOption == MARKOV ) os << "\t" << repEffects[i];
        os << endl;
    }

    // print the basal transcription
    os << "basal_transcription = " << basalTxps[ 0 ];
    for( int _i = 1; _i < basalTxps.size(); _i++ )
    {
        os << "\t" << basalTxps[ _i ];
    }
    os << endl;

    //print the pi vals
    for( int _i = 0; _i < pis.size(); _i++ )
    {
        if(_i > 0){
          os << "\t";
        }
        os << pis[ _i ];
    }
    os << endl;

    //print the beta values
    for( int i = 0; i < betas.size(); i++ )
    {
      if(i > 0){
        os << "\t";
      }
        os << betas[ i ];
    }
    os << endl;

    // print the cooperative interactions
    for ( int i = 0; i < nFactors(); i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
            if ( coopMat( i, j ) ) os << motifNames[i] << "\t" << motifNames[j] << "\t" << factorIntMat( i, j ) << endl;
        }
    }

    for( int i = 0; i < energyThrFactors.size(); i++ )
    {
        os << energyThrFactors[ i ] << "\t";
    }
    os << endl;

}
