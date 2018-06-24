#include <gsl/gsl_math.h>

#include "ExprPredictor.h"
#include "ExprPar.h"

//#define DEBUG

//It is a precondition that site_a.start <= site_b.start
#define ORDERED_SITE_OVERLAP(site_a, site_b) site_a.end < site_b.start

//It is a precondition that the sites do not overlap, and that site_b comes after site_a
#define SITE_DISTANCE(site_a, site_b) site_b.start - site_a.end


ExprFunc::ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num): expr_model(_model), par(_par), motifs( _model->motifs ), actIndicators( _model->actIndicators ), maxContact( _model->maxContact ), repIndicators( _model->repIndicators ), repressionMat( _model->repressionMat ), repressionDistThr( _model->repressionDistThr ), factorIntMat(_model->motifs.size(),_model->motifs.size(),1.0)
{
    //par = _par;//NOTE: made this const, and that solved a memory leak.

    int nFactors = par.nFactors();
    if(motifs.size() != nFactors){
        //std::cerr << " MOTFIS " << motifs.size() << " : factors par " << nFactors << std::endl;
        throw std::runtime_error(" Motifs and factors did not match");
    }
    assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );

    n_sites = sites_.size();//number of non-psudo sites.
    seq_number = seq_num;
    seq_length = seq_len;

    //setup legacy parameters
    maxBindingWts.clear();
    maxBindingWts.assign(nFactors,1.0);
    //maxBindingWts = vector < GEMSTAT_PAR_FLOAT_T >(nFactors);          // binding weight of the strongest site for each TF: K(S_max) [TF_max]
    txpEffects.clear();
    txpEffects.resize(nFactors,1.0);
    repEffects.clear();
    repEffects.resize(nFactors,1.0);

    std::map<std::string, int> tf_names_to_ids;
    for(int i = 0;i<expr_model->motifnames.size();i++){
        tf_names_to_ids[expr_model->motifnames.at(i)] = i;
    }

    //Do not assume that the tfs dictionary is in our internal order.
    for(int i = 0;i<nFactors;i++){
        const std::string &which_tf = expr_model->motifnames.at(i);
        maxBindingWts[i] = ((gsparams::DictList&)par.my_pars)["tfs"][which_tf]["maxbind"] ;
        txpEffects[i] = ((gsparams::DictList&)par.my_pars)["tfs"][which_tf]["alpha_a"] ;
        repEffects[i] = ((gsparams::DictList&)par.my_pars)["tfs"][which_tf]["alpha_r"] ;
    }

    //factorIntMat = Matrix(nFactors, nFactors);                      // (maximum) interactions between pairs of factors: omega(f,f')
    //Order doesn't matter for factor interactions, positions will be looked up.
    for(int k = 0;k<((gsparams::DictList&)par.my_pars)["inter"].size();k++){
        //need to split the name.
        std::string key = ((gsparams::DictList&)par.my_pars)["inter"].map_key_storage.at(k);
        double value = ((gsparams::DictList&)par.my_pars)["inter"].list_storage.at(k);

        int split_pos = key.find(":");
        int i = tf_names_to_ids[key.substr(0,split_pos)];//TODO: more defensive here. Make sure the value actually existed.
        int j = tf_names_to_ids[key.substr(split_pos+1,key.size())];
        factorIntMat.setElement(i,j,value);
        factorIntMat.setElement(j,i,value);
    }


    this->setupSitesAndBoundaries(sites_,seq_length, seq_num);

}

void ExprFunc::setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num){
  #ifdef DEBUG
  cerr << "running ExprFunc::setupSitesAndBoundaries(...)" << endl;
  #endif

	
  int n = _sites.size();
  sites = SiteVec(_sites);

  int pseudo_start = -1000;
  int pseudo_end = length+1000;
  if(sites.size() > 0){
	pseudo_start = _sites[0].start-1000;
	pseudo_end = _sites[_sites.size()-1].end+1000;
  }

  sites.insert( sites.begin(), Site(pseudo_start,pseudo_start,true,-1,0.0,1.0) );        // start with a pseudo-site at position 0
  sites.push_back( Site(pseudo_end,pseudo_end,true,-1,0.0,1.0) );       //and another pseudo-site at the end

  boundaries.resize(n_sites+2);
  boundaries[0] = 0;//value for starting pseudosite
  double range = max( (double)expr_model->get_longest_coop_thr(), (double)repressionDistThr );
  for ( int i = 1; i <= n; i++ )
  {
      int j;
      for ( j = i - 1; j >= 1; j-- )
      {
          if ( SITE_DISTANCE(sites[j],sites[i]) > range ) break;
      }//If the loop never broke, j will leak with value 0.
      int boundary = j;
      boundaries[i] = ( boundary );
  }
  boundaries[n_sites+1] = n_sites;//last-non-psudoe-site position boundary for ending pseudosite?

}

void ExprFunc::setupBindingWeights(const vector< double >& factorConcs){
  bindingWts.resize(sites.size());
  bindingWts[0] = 1.0;  //for first pseudosite
  bindingWts[bindingWts.size()-1] = 1.0; //might or might not be a pseudosite, gets overwritten if not
  int n = n_sites;
  for ( int i = 1; i <= n; i++ )
  {
      bindingWts[i] = maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio ;
      /*
      double samee = maxBindingWts[sites[i].factorIdx]*factorConcs[sites[i].factorIdx]*sites[i].prior_probability*sites[i].wtRatio;
      if(samee != samee)
      {
          cout << "DEBUG: samee for " << i << "\t" << sites[i].factorIdx << "\t" <<  maxBindingWts[sites[i].factorIdx] <<"\t" << factorConcs[sites[i].factorIdx] <<"\t" << sites[i].prior_probability <<"\t" << sites[i].wtRatio << endl;
          exit(1);
      }
      */
  }
}

double ExprFunc::predictExpr(const Condition& in_condition){
  double to_return = this->predictExpr( in_condition.concs );
  if(to_return < 0.0 || to_return != to_return){
	cerr << "returning a nonsense expression prediction " << to_return << endl;
	exit(1);
  }
  assert(to_return >= 0.0); //Positive expression
  assert(!(to_return != to_return)); //not NaN
  return to_return;
}

double ExprFunc::predictExpr( const vector< double >& factorConcs )
{

    // compute the Boltzman weights of binding for all sites
    setupBindingWeights(factorConcs);

    // Thermodynamic models: Direct, Quenching, ChrMod_Unlimited and ChrMod_Limited
    // compute the partition functions
    gemstat_dp_t Z_off = compPartFuncOff();
    //cout << "Z_off = " << Z_off << endl;
    gemstat_dp_t Z_on = compPartFuncOn();
    //cout << "Z_on = " << Z_on << endl;

    // compute the expression (promoter occupancy)
    gemstat_dp_t efficiency = Z_on / Z_off;
    //cout << "efficiency = " << efficiency << endl;
    //cout << "basalTxp = " << basalTxps[ seq_num ] << endl;

    GEMSTAT_PROMOTER_DATA_T my_promoter = par.getPromoterData( this->seq_number );

    gemstat_dp_t promoterOcc = efficiency * my_promoter.basal_trans / ( 1.0 + efficiency * my_promoter.basal_trans /** ( 1 + my_promoter.pi )*/ );
    #ifdef DEBUG
    if(promoterOcc < 0.0 || promoterOcc != promoterOcc){
	cerr << "Ridiculous in Direct!" << endl;
	cerr << "efficiency " << efficiency << endl;
	cerr << "Z_on" << Z_on << endl;
	cerr << "Z_off" << Z_off << endl;
	cerr << "basal " << my_promoter.basal_trans << endl; //TODO: I think I just found the bug.
	cerr << "=====" << endl;
	}
    #else
    assert(promoterOcc >= 0.0 && promoterOcc == promoterOcc);
    #endif
    return promoterOcc;
}

double Logistic_ExprFunc::predictExpr( const vector< double >& factorConcs ){

  // compute the Boltzman weights of binding for all sites
  setupBindingWeights(factorConcs);

  GEMSTAT_PROMOTER_DATA_T my_promoter = par.getPromoterData( this->seq_number );

  // total occupancy of each factor
  vector< gemstat_dp_t > factorOcc( motifs.size(), 0 );
  for ( int i = 1; i <= n_sites; i++ )
  {
      factorOcc[ sites[i].factorIdx ] += bindingWts[i] / ( 1.0 + bindingWts[i] );
  }
  gemstat_dp_t totalEffect = 0;
  //         cout << "factor\toccupancy\ttxp_effect" << endl;
  for ( int i = 0; i < motifs.size(); i++ )
  {
      gemstat_dp_t effect = txpEffects[i] * factorOcc[i];
      totalEffect += effect;
      //             cout << i << "\t" << factorOcc[i] << "\t" << effect << endl;

      // length correction
      //             totalEffect = totalEffect / (double)length;
  }
  //         return expRatio * logistic( log( my_promoter.basal_trans ) + totalEffect );
  return logistic( my_promoter.basal_trans + totalEffect );
}

Markov_ExprFunc::Markov_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){

    //Additional setup for a Markov_ExprFunc
    //void Markov_ExprFunc::setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num){

    #ifdef DEBUG
    cerr << "running Markov_ExprFunc::setupSitesAndBoundaries(...)" << endl;
    #endif

    double range = max( (double)expr_model->get_longest_coop_thr(), (double)repressionDistThr );
    rev_bounds.resize(n_sites+2);
    rev_bounds[n_sites] = n_sites+1; //last true site points to pseudosite
    rev_bounds[0] = 1; // first pseudosite has reverse boundary pointing to first true site
    for ( int i = n_sites; i > 0; i-- )
    {
        int j;
        for ( j = i+1; j <= n_sites; j++ )
        {
				//site_distance arugments backwards because we are going backwards.
            if ( SITE_DISTANCE(sites[i], sites[j]) > range ) break;
        }//If the loop never broke, j will leak with value n+1
        int boundary = j;
        rev_bounds[i] = ( boundary );
    }
    assert(rev_bounds[n_sites] == n_sites+1);
}

double Markov_ExprFunc::predictExpr( const vector< double >& factorConcs )
{
  int seq_num = this->seq_number;

  int n = n_sites;
  #ifdef DEBUG
  cout << "SITES size : " << sites.size() << " : n_sites : " << n_sites << endl;
  assert(sites.size() == n_sites+2);
  #endif

  /*
  cerr << "BOUNDARIES " << n << endl;
  cerr << boundaries << endl;
  cerr << "rev bounds " << endl;
  cerr << rev_bounds << endl;
  */
  setupBindingWeights(factorConcs);

    #ifdef DEBUG
    cerr << "Done setting bindingWts" << endl;
    #endif
    // initialization
    vector< gemstat_dp_t > Z( n + 2 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 2 );
    Zt[0] = 1.0;

    vector< gemstat_dp_t > backward_Z(n+2,0.0);
    backward_Z[backward_Z.size()-1] = 1.0;
    vector< gemstat_dp_t > backward_Z_sum(n+1,0.0);
    vector< gemstat_dp_t > backward_Zt(n+2,0.0);
    backward_Zt[backward_Zt.size()-1] = 1.0;

    // recurrence forward
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    // recurrence backward
    for ( int i = n; i >= 1; i-- )
    {
        gemstat_dp_t sum = backward_Zt[rev_bounds[i]];
        for ( int j = rev_bounds[i] - 1; j > i; j-- )
        {
			//Arguments to siteOverlap and interaction are backwards because we go backwards.
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * backward_Z[ j ];
        }
        backward_Z_sum[i] = sum;
        backward_Z[ i ] =  sum*bindingWts[i];
        backward_Zt[i] = backward_Z[i] + backward_Zt[i + 1] ;
    }




    #ifdef DEBUG
    vector< gemstat_dp_t > final_Z(n_sites+2,0.0);
    vector< gemstat_dp_t > final_Zt(n_sites+2,0.0);//not used, for debug only.
    bool problem = false;
    #endif

    vector< double > bindprobs(n_sites+2,0.0);

    for(int i = 1;i<=n_sites;i++){
      //Notice the i+1, we are skipping the pseudosite.
      gemstat_dp_t one_final_Z = Z[i] * backward_Z_sum[i];
      #ifdef DEBUG
      final_Z[i] = one_final_Z;
      final_Zt[i] = Zt[i] * backward_Zt[i];
      #endif
      //bindprobs[i] = final_Z[i] / final_Zt[i];
      bindprobs[i] = one_final_Z / backward_Zt[1];
      #ifdef DEBUG
      if( bindprobs[i] <= 0.0 || bindprobs[i] >= 1.0){
        problem = true;
      }
      #else
      //assert(bindprobs[i] >= 0.0);
      //assert(bindprobs[i] <= 1.0);
      if(bindprobs[i] < 0.0){
	bindprobs[i] = 0.0;
      } else if(bindprobs[i] > 1.0){
	bindprobs[i] = 1.0;
	}
      #endif
    }

    #ifdef DEBUG
    if(problem){
    cout << endl;
    cerr << "=====DEBUG=====" << endl;
    cerr << "Forward_Zt " << endl << Zt << endl;
    cerr << "====" << endl;
    cerr << "Backward_Zt " << endl << backward_Zt << endl;
    cerr << "====" << endl;
    cerr << "final_Zt " << endl << final_Zt << endl;
    cerr << "====" << endl;
    cerr << "final_Z " << endl << final_Z << endl;
    cerr << "====" << endl;
    cerr << "=====END=======" << endl;
    }
    #endif
    return this->expr_from_config(bindprobs);
}

double Markov_ExprFunc::expr_from_config(const vector< double >& marginals){
  double sum_total = 0.0;
  int n = n_sites;

  GEMSTAT_PROMOTER_DATA_T my_promoter = par.getPromoterData( seq_number );

  assert(n_sites + 2 == marginals.size());

  for(int i = 1; i <= n_sites; i++){
    //Iterating over the sites, not counting beginning or ending pseudosites

    double log_effect = 0.0;

    if( actIndicators[ sites[ i ].factorIdx ] )
    {
        log_effect = log(txpEffects[ sites[ i ].factorIdx ]);
    }
    if( repIndicators[ sites[ i ].factorIdx ] )
    {
        log_effect = log(repEffects[ sites[ i ].factorIdx ]);
    }

    sum_total += log_effect*marginals[i];
  }

  //TODO: do I need to make it negative?....
  //TODO: check one_qbtm_per_crm only once, higher up.
  double Z_on = exp(sum_total)*my_promoter.basal_trans;
  return Z_on / (1.0 + Z_on);

}

//ModelType ExprFunc::modelOption = QUENCHING;

gemstat_dp_t ExprFunc::compPartFuncOff() const
{
    #ifdef DEBUG
      //assert(modelOption != CHRMOD_UNLIMITED && modelOption != CHRMOD_LIMITED );
    #endif

    int n = n_sites;
    // initialization
    vector< gemstat_dp_t > Z( n + 1 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        if( sum != sum )
        {
            cout << "DEBUG: sum nan" << "\t" << Zt[ boundaries[i] ] <<  endl;
            exit(1);
        }
        //cout << "DEBUG: sum = " << n << endl;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            //cout << "compFactorInt: " << compFactorInt( sites[ j ], sites[ i ] ) << "\t";
            //cout << "Z[j]: " << Z[ j ] << endl;
            gemstat_dp_t old_sum = sum;
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
            if( sum != sum || isinf( sum ))
            {
                cout << "Old sum:\t" << old_sum << endl;
                cout << "Factors:\t" << sites[ i ].factorIdx << "\t" << sites[ j ].factorIdx << endl;
                cout << "compFactorInt:\t" << compFactorInt( sites[ j ], sites[ i ] ) << endl;
                cout << "Z[j]:\t" << Z[ j ] << endl;
                cout << i << "\t" << j << "\t" << factorIntMat( (sites[i]).factorIdx, (sites[j]).factorIdx ) << endl;
                cout << "DEBUG: sum nan/inf\t"<< sum << endl;
                exit(1);
            }
        }

        Z[i] = bindingWts[ i ] * sum;
        if( Z[i]!=Z[i] )
        {
            cout << "DEBUG: Z bindingWts[i]: " << sites[i].factorIdx << "\t" << bindingWts[ sites[i].factorIdx ] <<"\t" << sum << endl;
            exit(1);
        }
        Zt[i] = Z[i] + Zt[i - 1];
        //cout << "debug: Zt[i] = " << Zt[i] << endl;
    }

    // the partition function
    // 	gemstat_dp_t Z_bind = 1;
    // 	for ( int i = 0; i < sites.size(); i++ ) {
    // 		Z_bind += Z[ i ];
    // 	}
    return Zt[n];
}


gemstat_dp_t ChrMod_ExprFunc::compPartFuncOff() const
{
    int n = n_sites;

    // initialization
    vector< gemstat_dp_t > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< gemstat_dp_t > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        gemstat_dp_t sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            gemstat_dp_t dist = SITE_DISTANCE(sites[i], sites[j]);

            // sum for Z0
            sum0 += compFactorInt( sites[j], sites[i] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j];

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] )
            {
                sum1 += compFactorInt( sites[j], sites[i] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];
}


gemstat_dp_t ExprFunc::compPartFuncOn() const
{
    /*
    if ( modelOption == DIRECT ) assert(false);//should never make it here.
    if ( modelOption == QUENCHING ) assert(false);
    if ( modelOption == CHRMOD_UNLIMITED) assert(false);//return compPartFuncOnChrMod_Unlimited();
    if ( modelOption == CHRMOD_LIMITED ) assert(false);//return compPartFuncOnChrMod_Limited();
    */
//TODO: A compiler warning is generated here. Shouldn't there be some defensive coding?
    assert(false);
    return 0.0;
}


gemstat_dp_t Direct_ExprFunc::compPartFuncOn() const
{
    int n = n_sites;

    // initialization
    vector< gemstat_dp_t > Z( n + 1 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
        }
        //Z[i] = bindingWts[ i ] * txpEffects[ sites[i].factorIdx ] * sum;
        if( actIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = bindingWts[ i ] * txpEffects[ sites[ i ].factorIdx ] * sum;
            //cout << "1: " << txpEffects[ sites[ i ].factorIdx ] << endl;
        }
        if( repIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = bindingWts[ i ] * repEffects[ sites[ i ].factorIdx ] * sum;
            //cout << "2: " << repEffects[ sites[ i ].factorIdx ] << endl;
        }
        //cout << "DEBUG 0: " << sum << "\t" << Zt[ i - 1] << endl;
        Zt[i] = Z[i] + Zt[i - 1];
        /*if( actIndicators[ sites[ i ].factorIdx ] )
            cout << "DEBUG 1: " << Zt[i] << "\t" << bindingWts[i]*txpEffects[sites[i].factorIdx]*(Zt[ i - 1] + 1) << endl;
        if( repIndicators[ sites[ i ].factorIdx ] )
            cout << "DEBUG 2: " << Zt[i] << "\t" << bindingWts[i]*repEffects[sites[i].factorIdx]*(Zt[ i - 1] + 1) << endl;*/
    }

    return Zt[n];
}


gemstat_dp_t Quenching_ExprFunc::compPartFuncOn() const
{
    int n = n_sites;
    int N0 = maxContact;
    Matrix Z1(n+1, N0+1);
    Matrix Z0(n+1, N0+1);

    // k = 0
    for ( int i = 0; i <= n; i++ )
    {
        gemstat_dp_t sum1 = 1, sum0 = 0;
        for ( int j = 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            bool R = testRepression( sites[j], sites[i] );
            gemstat_dp_t term = compFactorInt( sites[ j ], sites[ i ] ) * ( Z1.getElement(j,0) + Z0.getElement(j,0) );
            sum1 += ( 1 - R )* term;
            sum0 += R * term;
        }
	Z1.setElement(i,0, bindingWts[i] * sum1);
	Z0.setElement(i,0, bindingWts[i] * sum0);
    }

    // k >= 1
    for ( int k = 1; k <= N0; k++ )
    {
        for ( int i = 0; i <= n; i++ )
        {
            if ( i < k )
            {
                Z1.setElement(i,k,0.0);
                Z0.setElement(i,k,0.0);
                continue;
            }
            gemstat_dp_t sum1 = 0, sum0 = 0;
            for ( int j = 1; j < i; j++ )
            {
                if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
                bool R = testRepression( sites[j], sites[i] );
                gemstat_dp_t effect = actIndicators[sites[j].factorIdx] * ( 1 - testRepression( sites[i], sites[j] ) ) * Z1.getElement(j,k-1) * txpEffects[sites[j].factorIdx];
                gemstat_dp_t term = compFactorInt( sites[ j ], sites[ i ] ) * ( Z1.getElement(j,k) + Z0.getElement(j,k) + effect );
                sum1 += ( 1 - R )* term;
                sum0 += R * term;
            }
            Z1.setElement(i,k,bindingWts[i] * sum1);
            Z0.setElement(i,k,bindingWts[i] * sum0);
        }
    }

    //     for ( int i = 1; i <= n; i++ ) {
    //         for ( int k = 0; k <= N0; k++ ) {
    //             cout << "Z1(" << i << ", " << k << ") = " << Z1[i][k] << "\t";
    //             cout << "Z0(" << i << ", " << k << ") = " << Z0[i][k] << endl;
    //         }
    //         cout << endl;
    //     }

    // the partition function
    gemstat_dp_t Z_on = 1;
    for ( int i = 1; i <= n; i++ )
    {
        for ( int k = 0; k <= N0; k++ )
        {
            gemstat_dp_t term = Z1.getElement(i,k) + Z0.getElement(i,k);
            Z_on += term;
        }
        for ( int k = 0; k <= N0 - 1; k++ )
        {
            Z_on += actIndicators[sites[i].factorIdx] * Z1.getElement(i,k) * txpEffects[sites[i].factorIdx];
        }
    }
    return Z_on;
}


gemstat_dp_t ChrModUnlimited_ExprFunc::compPartFuncOn() const
{
    int n = n_sites;

    // initialization
    vector< gemstat_dp_t > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< gemstat_dp_t > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        gemstat_dp_t sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
			//Remember that for current site i, we are comparing to previous sites j.
            gemstat_dp_t dist = SITE_DISTANCE(sites[j], sites[i]);
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;

            // sum for Z0
            sum0 += compFactorInt( sites[j], sites[i] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j];

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] )
            {
                sum1 += compFactorInt( sites[j], sites[i] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * txpEffects[ sites[i].factorIdx ] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];
}


gemstat_dp_t ChrModLimited_ExprFunc::compPartFuncOn() const
{
    int n = n_sites;

    // initialization
    int N0 = maxContact;
    Matrix Z0( n + 1, N0 + 1 );
    Matrix Z1( n + 1, N0 + 1 );
    Matrix Zt( n + 1, N0 + 1 );
    Z0.setElement(0,0,0.0);
    Z1.setElement(0,0,0.0);
    Zt.setElement(0,0,1.0);
    for ( int k = 1; k <= N0; k++ )
    {
        Z0.setElement(0,k,0.0);
        Z1.setElement(0,k,0.0);
        Zt.setElement(0,k,0.0);
    }

    // recurrence
    for ( int k = 0; k <= N0; k++ )
    {
        for ( int i = 1; i <= n; i++ )
        {
            //             cout << "k = " << k << " i = " << i << endl;
            gemstat_dp_t sum0 = Zt.getElement(boundaries[i],k);
	    gemstat_dp_t sum0A = k > 0 ? Zt.getElement(boundaries[i],k-1) : 0.0;
	    gemstat_dp_t sum1 = sum0;

            for ( int j = boundaries[i] + 1; j < i; j++ )
            {
                double dist = SITE_DISTANCE(sites[j], sites[i]);//Remeber order of sites.
                if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;

                // sum for Z0
                sum0 += compFactorInt( sites[j], sites[i] ) * Z0.getElement(j,k);
                sum0A += k > 0 ? compFactorInt( sites[j], sites[i] ) * Z0.getElement(j,k-1) : 0;
                if ( dist > repressionDistThr )
                {
                    sum0 += Z1.getElement(j,k);
                    sum0A += k > 0 ? Z1.getElement(j,k-1) : 0;
                }

                // sum for Z1
                if ( repIndicators[ sites[i].factorIdx ] )
                {
                    sum1 += compFactorInt( sites[j], sites[i] ) * Z1.getElement(j,k);
                    if ( dist > repressionDistThr ) sum1 += Z0.getElement(j,k);
                }
            }
            Z0.setElement(i,k,bindingWts[i] * sum0);
            if ( actIndicators[sites[i].factorIdx] ) Z0(i,k) += k > 0 ? bindingWts[i] * txpEffects[sites[i].factorIdx] * sum0A : 0;
            if ( repIndicators[ sites[i].factorIdx ] ) Z1(i,k) = bindingWts[i] * repEffects[ sites[i].factorIdx ] * sum1;
            else Z1.setElement(i,k,0.0);
            Zt.setElement(i,k,Z0.getElement(i,k) + Z1.getElement(i,k) + Zt.getElement(i - 1,k));
            //             cout << "i = " << i << " k = " << k << " Z0 = " << Z0[i][k] << " Z1 = " << Z1[i][k] << " Zt = " << Zt[i][k] << endl;
        }
    }

    // the partition function
    //     cout << "Zt[n] = " << Zt[n] << endl;
    return sum( Zt.getRow(n) );//And we end up with a vector anyway. See about fixing this.
}

/**
	PRE: Site a must start at or before the same base pair as Site b.
*/
double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
	#ifdef DEBUG
		assert( a.start <= b.start);
	#endif
    // 	assert( !siteOverlap( a, b, motifs ) );
    double maxInt = factorIntMat( a.factorIdx, b.factorIdx );
    double dist = SITE_DISTANCE(a, b);
    //bool orientation = ( a.strand == b.strand );

    FactorIntFunc* an_int_func = expr_model->coop_setup->coop_func_for(a.factorIdx, b.factorIdx);
    return an_int_func->compFactorInt( maxInt, dist, a.strand, b.strand );

    //TODO: we need to get this information from the expr_model.
    //This is going to be very slow, we should have cached it.
}

/**
	PRE: Site a must start at or before the same base pair as Site b.
*/
bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );
	#ifdef DEBUG
		assert( a.start <= b.start);
	#endif

    double dist = SITE_DISTANCE(a, b);
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}
