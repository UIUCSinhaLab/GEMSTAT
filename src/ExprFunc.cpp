#include <gsl/gsl_math.h>

#include "ExprPredictor.h"
#include "ExprPar.h"


ExprFunc::ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : motifs( _motifs ), intFunc( _intFunc ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr )
{
    par = _par;

    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );

}

void ExprFunc::setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num){
  boundaries.clear();

  if( !one_qbtm_per_crm )
      seq_num = 0;
  // store the sequence
  int n = _sites.size();
  sites = SiteVec(_sites);
  sites.insert( sites.begin(), Site() );        // start with a pseudo-site at position 0
  boundaries.push_back( 0 );
  double range = max( intFunc->getMaxDist(), repressionDistThr );
  for ( int i = 1; i <= n; i++ )
  {
      int j;
      for ( j = i - 1; j >= 1; j-- )
      {
          if ( ( sites[i].start - sites[j].start ) > range ) break;
      }
      int boundary = j;
      boundaries.push_back( boundary );
  }
}

void ExprFunc::setupBindingWeights(const vector< double >& factorConcs){
  bindingWts.clear();
  bindingWts.push_back( 1.0 );
  int n = sites.size();
  for ( int i = 1; i < n; i++ )
  {
      bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio );
      /*
      double samee = par.maxBindingWts[sites[i].factorIdx]*factorConcs[sites[i].factorIdx]*sites[i].prior_probability*sites[i].wtRatio;
      if(samee != samee)
      {
          cout << "DEBUG: samee for " << i << "\t" << sites[i].factorIdx << "\t" <<  par.maxBindingWts[sites[i].factorIdx] <<"\t" << factorConcs[sites[i].factorIdx] <<"\t" << sites[i].prior_probability <<"\t" << sites[i].wtRatio << endl;
          exit(1);
      }
      */
  }
}

double ExprFunc::predictExpr( const SiteVec& _sites, int length, const Condition& in_condition, int seq_num ){
  return this->predictExpr(_sites,length, in_condition.concs, seq_num);
}

double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
    //initialize the sites and boundaries and whatnot.
    setupSitesAndBoundaries(_sites,length,seq_num);

    // compute the Boltzman weights of binding for all sites
    setupBindingWeights(factorConcs);

    // Thermodynamic models: Direct, Quenching, ChrMod_Unlimited and ChrMod_Limited
    // compute the partition functions
    double Z_off = compPartFuncOff();
    //cout << "Z_off = " << Z_off << endl;
    double Z_on = compPartFuncOn();
    //cout << "Z_on = " << Z_on << endl;

    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    //cout << "efficiency = " << efficiency << endl;
    //cout << "basalTxp = " << par.basalTxps[ seq_num ] << endl;
    double promoterOcc = efficiency * par.basalTxps[ seq_num ] / ( 1.0 + efficiency * par.basalTxps[ seq_num ] /** ( 1 + par.pis[ seq_num ] )*/ );
    return promoterOcc;
}

double Logistic_ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num ){
  //initialize the sites and boundaries and whatnot.
  setupSitesAndBoundaries(_sites,length,seq_num);

  // compute the Boltzman weights of binding for all sites
  setupBindingWeights(factorConcs);


  // total occupancy of each factor
  vector< double > factorOcc( motifs.size(), 0 );
  for ( int i = 1; i < sites.size(); i++ )
  {
      factorOcc[ sites[i].factorIdx ] += bindingWts[i] / ( 1.0 + bindingWts[i] );
  }
  double totalEffect = 0;
  //         cout << "factor\toccupancy\ttxp_effect" << endl;
  for ( int i = 0; i < motifs.size(); i++ )
  {
      double effect = par.txpEffects[i] * factorOcc[i];
      totalEffect += effect;
      //             cout << i << "\t" << factorOcc[i] << "\t" << effect << endl;

      // length correction
      //             totalEffect = totalEffect / (double)length;
  }
  //         return par.expRatio * logistic( log( par.basalTxp ) + totalEffect );
  return logistic( par.basalTxps[ seq_num ] + totalEffect );
}

double Markov_ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
  boundaries.clear();

  if( !one_qbtm_per_crm )
      seq_num = 0;
  // store the sequence
  int n = _sites.size();
  sites = SiteVec(_sites);
  sites.insert( sites.begin(), Site() );        // start with a pseudo-site at position 0
  sites.push_back(Site());//and another at the end
  sites[sites.size()-1].start = length + 1000;
  boundaries.push_back( 0 );
  double range = max( intFunc->getMaxDist(), repressionDistThr );
  for ( int i = 1; i <= n; i++ )
  {
      int j;
      for ( j = i - 1; j >= 1; j-- )
      {
          if ( ( sites[i].start - sites[j].start ) > range ) break;
      }
      int boundary = j;
      boundaries.push_back( boundary );
  }
  boundaries.push_back(n);//self-boundary for ending pseudosite?

  vector<int> rev_bounds(boundaries.size(), n+1);
  rev_bounds[rev_bounds.size()-1] = n+1;
  rev_bounds[0] = 1;
  for ( int i = n; i > 0; i-- )
  {
      int j;
      for ( j = i+1; j <= n; j++ )
      {
          if ( ( sites[j].start - sites[i].start) > range ) break;
      }
      int boundary = j;
      rev_bounds[i] = ( boundary );
  }
  rev_bounds.push_back(n+1);
  assert(rev_bounds[n] == n+1);

  /*
  cerr << "BOUNDARIES " << n << endl;
  cerr << boundaries << endl;
  cerr << "rev bounds " << endl;
  cerr << rev_bounds << endl;
  */

  bindingWts.clear();
  bindingWts.push_back( 1.0 );
  for ( int i = 1; i <= n; i++ )
  {
      bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio );
  }
  bindingWts.push_back(1.0);

    // initialization
    vector< double > Z( n + 2 );
    Z[0] = 1.0;
    vector< double > Zt( n + 2 );
    Zt[0] = 1.0;

    vector< double > backward_Z(n+2,0.0);
    backward_Z[backward_Z.size()-1] = 1.0;
    vector< double > backward_Zt(n+2,0.0);
    backward_Zt[backward_Zt.size()-1] = 1.0;

    // recurrence forward
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    // recurrence backward
    for ( int i = n; i >= 1; i-- )
    {
        double sum = backward_Zt[rev_bounds[i]];
        for ( int j = rev_bounds[i] - 1; j > i; j-- )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * backward_Z[ j ];
        }
        backward_Z[ i ] =  sum*bindingWts[i];
        backward_Zt[i] = backward_Z[i] + backward_Zt[i + 1] ;
    }



    vector< double > final_Z(Zt.size()-2,0.0);
    vector< double > final_Zt(Zt.size()-2,0.0);

    vector< double > bindprobs(Zt.size()-2,0.0);

    for(int i = 0;i<n;i++){
      final_Z[i] = Z[i+1] * backward_Z[i+1];
      final_Zt[i] = Zt[i+1] * backward_Zt[i+1];
      bindprobs[i] = final_Z[i] / final_Zt[i];
    }

    /*
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
    /**/
    return this->expr_from_config(_sites, length, seq_num, bindprobs);
}

double Markov_ExprFunc::expr_from_config(const SiteVec& _sites, int length, int seq_num, const vector< double >& marginals){
  double sum_total = 0.0;
  int n = _sites.size();


  assert(_sites.size() == marginals.size());

  for(int i = 0; i < n; i++){

    double log_effect = 0.0;

    if( actIndicators[ _sites[ i ].factorIdx ] )
    {
        log_effect = log(par.txpEffects[ _sites[ i ].factorIdx ]);
    }
    if( repIndicators[ _sites[ i ].factorIdx ] )
    {
        log_effect = log(par.repEffects[ _sites[ i ].factorIdx ]);
    }

    sum_total += log_effect*marginals[i];
  }

  //TODO: do I need to make it negative?....
  //TODO: check one_qbtm_per_crm only once, higher up.
  double Z_on = exp(sum_total)*par.basalTxps[ one_qbtm_per_crm ? 0 : seq_num ];
  return Z_on / (1.0 + Z_on);

}

ModelType ExprFunc::modelOption = QUENCHING;

double ExprFunc::compPartFuncOff() const
{
    #ifdef DEBUG
      assert(modelOption != CHRMOD_UNLIMITED && modelOption != CHRMOD_LIMITED );
    #endif

    int n = sites.size() - 1;
    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        if( sum != sum )
        {
            cout << "DEBUG: sum nan" << "\t" << Zt[ boundaries[i] ] <<  endl;
            exit(1);
        }
        //cout << "DEBUG: sum = " << n << endl;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            //cout << "compFactorInt: " << compFactorInt( sites[ i ], sites[ j ] ) << "\t";
            //cout << "Z[j]: " << Z[ j ] << endl;
            double old_sum = sum;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
            if( sum != sum || isinf( sum ))
            {
                cout << "Old sum:\t" << old_sum << endl;
                cout << "Factors:\t" << sites[ i ].factorIdx << "\t" << sites[ j ].factorIdx << endl;
                cout << "compFactorInt:\t" << compFactorInt( sites[ i ], sites[ j ] ) << endl;
                cout << "Z[j]:\t" << Z[ j ] << endl;
                cout << i << "\t" << j << "\t" << par.factorIntMat( (sites[i]).factorIdx, (sites[j]).factorIdx ) << endl;
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
    // 	double Z_bind = 1;
    // 	for ( int i = 0; i < sites.size(); i++ ) {
    // 		Z_bind += Z[ i ];
    // 	}
    return Zt[n];
}


double ChrMod_ExprFunc::compPartFuncOff() const
{
    int n = sites.size()- 1;

    // initialization
    vector< double > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< double > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        double sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            double dist = sites[i].start - sites[j].start;

            // sum for Z0
            sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j];

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] )
            {
                sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];
}


double ExprFunc::compPartFuncOn() const
{
    if ( modelOption == DIRECT ) assert(false);//should never make it here.
    if ( modelOption == QUENCHING ) assert(false);
    if ( modelOption == CHRMOD_UNLIMITED) assert(false);//return compPartFuncOnChrMod_Unlimited();
    if ( modelOption == CHRMOD_LIMITED ) assert(false);//return compPartFuncOnChrMod_Limited();
//TODO: A compiler warning is generated here. Shouldn't there be some defensive coding?
    assert(false);
}


double Direct_ExprFunc::compPartFuncOn() const
{
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        //Z[i] = bindingWts[ i ] * par.txpEffects[ sites[i].factorIdx ] * sum;
        if( actIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = bindingWts[ i ] * par.txpEffects[ sites[ i ].factorIdx ] * sum;
            //cout << "1: " << par.txpEffects[ sites[ i ].factorIdx ] << endl;
        }
        if( repIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = bindingWts[ i ] * par.repEffects[ sites[ i ].factorIdx ] * sum;
            //cout << "2: " << par.repEffects[ sites[ i ].factorIdx ] << endl;
        }
        //cout << "DEBUG 0: " << sum << "\t" << Zt[ i - 1] << endl;
        Zt[i] = Z[i] + Zt[i - 1];
        /*if( actIndicators[ sites[ i ].factorIdx ] )
            cout << "DEBUG 1: " << Zt[i] << "\t" << bindingWts[i]*par.txpEffects[sites[i].factorIdx]*(Zt[ i - 1] + 1) << endl;
        if( repIndicators[ sites[ i ].factorIdx ] )
            cout << "DEBUG 2: " << Zt[i] << "\t" << bindingWts[i]*par.repEffects[sites[i].factorIdx]*(Zt[ i - 1] + 1) << endl;*/
    }

    return Zt[n];
}


double Quenching_ExprFunc::compPartFuncOn() const
{
    int n = sites.size() - 1;
    int N0 = maxContact;
    Matrix Z1(n+1, N0+1);
    Matrix Z0(n+1, N0+1);

    // k = 0
    for ( int i = 0; i <= n; i++ )
    {
        double sum1 = 1, sum0 = 0;
        for ( int j = 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            bool R = testRepression( sites[j], sites[i] );
            double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1.getElement(j,0) + Z0.getElement(j,0) );
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
            double sum1 = 0, sum0 = 0;
            for ( int j = 1; j < i; j++ )
            {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                bool R = testRepression( sites[j], sites[i] );
                double effect = actIndicators[sites[j].factorIdx] * ( 1 - testRepression( sites[i], sites[j] ) ) * Z1.getElement(j,k-1) * par.txpEffects[sites[j].factorIdx];
                double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1.getElement(j,k) + Z0.getElement(j,k) + effect );
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
    double Z_on = 1;
    for ( int i = 1; i <= n; i++ )
    {
        for ( int k = 0; k <= N0; k++ )
        {
            double term = Z1.getElement(i,k) + Z0.getElement(i,k);
            Z_on += term;
        }
        for ( int k = 0; k <= N0 - 1; k++ )
        {
            Z_on += actIndicators[sites[i].factorIdx] * Z1.getElement(i,k) * par.txpEffects[sites[i].factorIdx];
        }
    }
    return Z_on;
}


double ChrModUnlimited_ExprFunc::compPartFuncOn() const
{
    int n = sites.size()- 1;

    // initialization
    vector< double > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< double > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        double sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            double dist = sites[i].start - sites[j].start;
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;

            // sum for Z0
            sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j];

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] )
            {
                sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * par.txpEffects[ sites[i].factorIdx ] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];
}


double ChrModLimited_ExprFunc::compPartFuncOn() const
{
    int n = sites.size()- 1;

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
            double sum0 = Zt.getElement(boundaries[i],k);
	    double sum0A = k > 0 ? Zt.getElement(boundaries[i],k-1) : 0.0;
	    double sum1 = sum0;

            for ( int j = boundaries[i] + 1; j < i; j++ )
            {
                double dist = sites[i].start - sites[j].start;
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;

                // sum for Z0
                sum0 += compFactorInt( sites[i], sites[j] ) * Z0.getElement(j,k);
                sum0A += k > 0 ? compFactorInt( sites[i], sites[j] ) * Z0.getElement(j,k-1) : 0;
                if ( dist > repressionDistThr )
                {
                    sum0 += Z1.getElement(j,k);
                    sum0A += k > 0 ? Z1.getElement(j,k-1) : 0;
                }

                // sum for Z1
                if ( repIndicators[ sites[i].factorIdx ] )
                {
                    sum1 += compFactorInt( sites[i], sites[j] ) * Z1.getElement(j,k);
                    if ( dist > repressionDistThr ) sum1 += Z0.getElement(j,k);
                }
            }
            Z0.setElement(i,k,bindingWts[i] * sum0);
            if ( actIndicators[sites[i].factorIdx] ) Z0(i,k) += k > 0 ? bindingWts[i] * par.txpEffects[sites[i].factorIdx] * sum0A : 0;
            if ( repIndicators[ sites[i].factorIdx ] ) Z1(i,k) = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1;
            else Z1.setElement(i,k,0.0);
            Zt.setElement(i,k,Z0.getElement(i,k) + Z1.getElement(i,k) + Zt.getElement(i - 1,k));
            //             cout << "i = " << i << " k = " << k << " Z0 = " << Z0[i][k] << " Z1 = " << Z1[i][k] << " Zt = " << Zt[i][k] << endl;
        }
    }

    // the partition function
    //     cout << "Zt[n] = " << Zt[n] << endl;
    return sum( Zt.getRow(n) );//And we end up with a vector anyway. See about fixing this.
}


double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );
    double maxInt = par.factorIntMat( a.factorIdx, b.factorIdx );
    double dist = abs( a.start - b.start );
    bool orientation = ( a.strand == b.strand );
    return intFunc->compFactorInt( maxInt, dist, orientation );
}


bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );

    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}
