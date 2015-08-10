#include "ExprPar.h"
#include "ExprPredictor.h"
#include "conf/ExprParConf.hpp"

ParFactory::ParFactory( const ExprModel& in_model, int in_nSeqs) : expr_model(in_model), nSeqs(in_nSeqs)
{
  //maximums = ExprPar( expr_model.motifs.size(), nSeqs );
  maximums = createDefaultMinMax(true);
  //minimums = ExprPar( expr_model.motifs.size(), nSeqs );
  minimums = createDefaultMinMax(false);
}

void ParFactory::separateParams(const ExprPar& input, vector<double>& free_output, vector<double>& fixed_output, const vector<bool>& indicator_bool)
{
  //Hassan's code for separating parameters
  vector<double> pars;

  input.getFreePars( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators );
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

void ParFactory::joinParams(const vector<double>& free_pars, const vector<double>& fix_pars, vector<double>& output, const vector<bool>& indicator_bool)
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

ExprPar ParFactory::createDefaultMinMax(bool min_or_max)
{
  //TODO: the model really should know about the number of factors in the model.
  int _nFactors = expr_model.motifs.size();
  assert( _nFactors > 0 );

  ExprPar tmp_par = ExprPar( _nFactors, nSeqs );

  tmp_par.maxBindingWts.assign(tmp_par.maxBindingWts.size(), min_or_max ? ExprPar::max_weight : ExprPar::min_weight); // ExprPar::min_weight
  //set the interaction maximums
  tmp_par.factorIntMat.setAll(min_or_max ? ExprPar::max_interaction : ExprPar::min_interaction);
  tmp_par.txpEffects.assign(tmp_par.txpEffects.size(), min_or_max ? ExprPar::max_effect_Thermo : ExprPar::min_effect_Thermo);//TODO: Doesn't handle Logistic model
  tmp_par.repEffects.assign(tmp_par.repEffects.size(), min_or_max ? ExprPar::max_repression : ExprPar::min_repression);
  tmp_par.basalTxps.assign(tmp_par.basalTxps.size(), min_or_max ? ExprPar::max_basal_Thermo : ExprPar::min_basal_Thermo);//TODO: Doesn't handle Logistic model
  tmp_par.pis.assign(tmp_par.pis.size(), min_or_max ? ExprPar::max_pi : ExprPar::min_pi);
  tmp_par.betas.assign(tmp_par.betas.size(), min_or_max ? ExprPar::max_beta : ExprPar::min_beta);
  tmp_par.energyThrFactors.assign(tmp_par.energyThrFactors.size(), min_or_max ? ExprPar::max_energyThrFactors : ExprPar::min_energyThrFactors);
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
    pis.assign( numBTMS, ExprPar::default_pi );

    //for the beta parameters
    betas.assign( numBTMS, ExprPar::default_beta );

    energyThrFactors.assign( numBTMS, ExprPar::default_energyThrFactors );
}


ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector <double>& _energyThrFactors ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxps( _basalTxps ), pis( _pis), betas( _betas ), nSeqs( _nSeqs ), energyThrFactors( _energyThrFactors )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() );
    assert( txpEffects.size() == maxBindingWts.size() && repEffects.size() == maxBindingWts.size() );
    assert( basalTxps.size() == one_qbtm_per_crm ? nSeqs : 1);
}


ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs ) : factorIntMat()
{

    int _nFactors = actIndicators.size();
    nSeqs = _nSeqs;
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
    //     assert( pars.size() == ( _nFactors * ( _nFactors + 1 ) / 2 + 2 * _nFactors + 2 );
    int counter = 0;

    // set maxBindingWts
    if ( estBindingOption )
    {
        for ( int i = 0; i < _nFactors; i++ )
        {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) : exp( pars[counter++] );
            maxBindingWts.push_back( weight );
        }
    }
    else
    {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
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
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT )
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


void ExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();

    // write maxBindingWts
    if ( estBindingOption )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            if( weight != weight )
            {
                cout << "DEBUG samee: getFreePars() " << i << endl;
                exit(1);
            }
            pars.push_back( weight );
        }
    }

    // write the interaction matrix
    if ( modelOption != LOGISTIC )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            for ( int j = 0; j <= i; j++ )
            {
                if ( coopMat( i, j ) )
                {
                    double interaction = searchOption == CONSTRAINED ? infty_transform( log( factorIntMat( i, j ) ), log( min_interaction ), log( max_interaction ) ) : log( factorIntMat( i, j ) );
                    pars.push_back( interaction );
                }
            }
        }
    }

    // write the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( modelOption == LOGISTIC )
        {
            double effect = searchOption == CONSTRAINED ? infty_transform( txpEffects[i], min_effect_Logistic, max_effect_Logistic ) : txpEffects[i];
            pars.push_back( effect );
        }
        /*else if ( modelOption == DIRECT ) {
                    double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
                    pars.push_back( effect );
                }*/
        else
        {
            if ( actIndicators[i] )
            {
                double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
                pars.push_back( effect );
            }
        }
    }

    // write the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( repIndicators[i] )
            {
                double repression = searchOption == CONSTRAINED ? infty_transform( log( repEffects[i] ), log( min_repression ), log( max_repression ) ) : log( repEffects[i] );
                pars.push_back( repression );
            }
        }
    }

    for( int i = 0; i < basalTxps.size(); i++ )
    {
        // write the basal transcription
        if ( modelOption == LOGISTIC )
        {
            double basal = searchOption == CONSTRAINED ? infty_transform( basalTxps[ i ], min_basal_Logistic, max_basal_Logistic ) : basalTxps[ i ];
            pars.push_back( basal );
        }
        else
        {
            double basal = searchOption == CONSTRAINED ? infty_transform( log( basalTxps[ i ] ), log( min_basal_Thermo ), log( max_basal_Thermo ) ) : log( basalTxps[ i ] );
            pars.push_back( basal );
        }
    }

    //write the pis
    for( int i = 0; i < pis.size(); i++ )
    {
        double pi_val = searchOption == CONSTRAINED? infty_transform( log( pis[ i ] ), log( min_pi ), log( max_pi ) ) : log( pis[ i ] );
        pars.push_back( pi_val );
    }

    //write the betas
    for( int i = 0; i < betas.size(); i++ )
    {
        double beta_val = searchOption == CONSTRAINED ? infty_transform( log( betas[ i ] ), log( min_beta ), log( max_beta ) ) : log( betas[ i ] );
        pars.push_back( beta_val );
    }
    for( int i = 0; i < energyThrFactors.size(); i++ )
    {
        double energyThrFactor_val = searchOption == CONSTRAINED ? infty_transform( log( energyThrFactors[ i ]),log(min_energyThrFactors), log(max_energyThrFactors) ) : log( energyThrFactors[ i ] );
        pars.push_back( energyThrFactor_val );
    }
}


void ExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const
{
//    os.setf( ios::fixed );
//    os.precision( 50 );

    // print the factor information
    for ( int i = 0; i < nFactors(); i++ )
    {
        os << motifNames[i] << "\t" << maxBindingWts[i] << "\t" << txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT ) os << "\t" << repEffects[i];
        os << endl;
    }

    // print the basal transcription
    os << "basal_transcription = " << basalTxps[ 0 ] << endl;
    for( int _i = 1; _i < basalTxps.size(); _i++ )
    {
        os << basalTxps[ _i ] << endl;
    }

    //print the pi vals
    for( int _i = 0; _i < pis.size(); _i++ )
    {
        os << pis[ _i ] << endl;
    }
    //print the beta values
    for( int i = 0; i < betas.size(); i++ )
    {
        os << betas[ i ] << endl;
    }
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


int ExprPar::load( const string& file, const int num_of_coop_pairs )
{
    // open the file
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open parameter file " << file << endl; exit( 1 ); }

    // read the factor information
    vector< string > motifNames( nFactors() );
    for ( int i = 0; i < nFactors(); i++ )
    {
        fin >> motifNames[i] >> maxBindingWts[i] >> txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT ) fin >> repEffects[i];
    }

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ )
    {
        factorIdxMap[motifNames[i]] = i;
    }

    // read the basal transcription
    string symbol, eqSign, value;
    fin >> symbol >> eqSign >> value;
    if ( symbol != "basal_transcription" || eqSign != "=" ) return RET_ERROR;
    double basalTxp_val = atof( value.c_str() );
    basalTxps[ 0 ] =  basalTxp_val ;
    if( one_qbtm_per_crm )
    {
        for( int _i = 1; _i < nSeqs; _i++ )
        {
            fin >> value;
            double basalTxp_val = atof( value.c_str() );
            basalTxps[ _i ] = basalTxp_val;
        }
    }
    //read the pi values
    for( int i = 0; i < nSeqs; i++ )
    {
        fin >> value;
        pis[ i ] = atof( value.c_str() );
    }

    //read the beta values
    for( int i = 0; i < nSeqs; i++ )
    {
        fin >> value;
        betas[ i ] = atof( value.c_str() );
    }

    // read the cooperative interactions
    string factor1, factor2;
    double coopVal;
    for( int i = 0; i < num_of_coop_pairs; i++ )
    {
        fin >> factor1 >> factor2 >> coopVal;
        if( !factorIdxMap.count( factor1 ) || !factorIdxMap.count( factor2 ) ) return RET_ERROR;
        int idx1 = factorIdxMap[factor1];
        int idx2 = factorIdxMap[factor2];
        factorIntMat( idx1, idx2 ) = coopVal;
        factorIntMat( idx2, idx1 ) = coopVal;
        //cout << factor1 << "\t" << factor2 << "\t" << idx1 << "\t" << idx2 << endl;
    }
    double factor_thr_val;
    energyThrFactors.clear();
    while( fin >> factor_thr_val )
    {
        energyThrFactors.push_back( factor_thr_val );
    }

    return 0;
}


void ExprPar::adjust( const IntMatrix& coopMat  )
{
    // adjust binding paramters
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) ) { maxBindingWts[i] *= 2.0;}
        if ( maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) ) { maxBindingWts[i] /= 2.0;}
    }

    // adjust the interaction matrix
    for ( int i = 0; i < nFactors(); i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
            if ( coopMat( i, j ) &&   factorIntMat( i, j ) < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) )
            {
                factorIntMat( i, j ) *= 2.0;
                factorIntMat( j, i ) = factorIntMat( i, j );
            }
            if ( coopMat( i, j ) &&  factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) )
            {
                factorIntMat( i, j ) /= 2.0;
                factorIntMat( j, i ) = factorIntMat( i, j );
            }
        }
    }
    // adjust transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( modelOption == LOGISTIC )
        {
            if ( txpEffects[i] < ExprPar::min_effect_Logistic + ExprPar::delta ) txpEffects[i] /= 2.0;
            if ( txpEffects[i] > ExprPar::max_effect_Logistic - ExprPar::delta ) txpEffects[i] /= 2.0;
        }
        else
        {
            if ( txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) ) { txpEffects[i] *= 2.0;}
            if ( txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) ) { txpEffects[i] /= 2.0;}
        }

    }

    // adjust the repression effects
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED || modelOption == DIRECT )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) )
            {
                if( repEffects[ i ] > 0 ) repEffects[i] *= 2.0;
                else repEffects[ i ] = ExprPar::min_repression * ( 1 + 2 * ExprPar::delta );
            }
            if ( repEffects[i] > ExprPar::max_repression * ( 1.0 - ExprPar::delta ) ) {repEffects[i] /= 2.0;}
        }
    }

    // adjust the basl transcription
    for( int _i = 0; _i < basalTxps.size(); _i ++ )
    {
        if ( modelOption == LOGISTIC )
        {
            if ( basalTxps[ _i ] < ExprPar::min_basal_Logistic + ExprPar::delta ) basalTxps[ _i ] /= 2.0;
            if ( basalTxps[ _i ] > ExprPar::max_basal_Logistic - ExprPar::delta ) basalTxps[ _i ] *= 2.0;
        }
        else
        {
            if ( basalTxps[ _i ] < ExprPar::min_basal_Thermo * ( 1.0 + ExprPar::delta ) ) basalTxps[ _i ] *= 2.0;
            if ( basalTxps[ _i ] > ExprPar::max_basal_Thermo * ( 1.0 - ExprPar::delta ) ) basalTxps[ _i ] /= 2.0;
        }
    }

    //adjust the pis
    for( int _i = 0; _i < pis.size(); _i++ )
    {
        if( pis[ _i ] < 0 ) pis[ _i] = 0;
        if( pis[ _i ] > ExprPar::max_pi * ( 1.0 - ExprPar::delta ) ) pis[ _i ] /= 2.0;
    }
    //adjust the betas
    for( int i = 0; i < betas.size(); i++ )
    {

        if( betas[ i ] < ExprPar::min_beta * ( 1.0 + ExprPar::delta ) ) betas[ i ] *= 2.0;
        if( betas[ i ] > ExprPar::max_beta )
        {
            betas[ i ] = ExprPar::max_beta * ( 1.0 - ExprPar::delta );
        }
        else if( betas[ i ] > ExprPar::max_beta * ( 1.0 - ExprPar::delta ) )
        {
            betas[ i ] /= 2.0;
        }
    }
    //adjust the energyThrFactors
    for( int i = 0; i < nFactors(); i++ )
    {
        if( energyThrFactors[ i ] < ExprPar::min_energyThrFactors * ( 1.0 + ExprPar::delta ) ) energyThrFactors[ i ] *= 2.0;
        if( energyThrFactors[ i ] > ExprPar::max_energyThrFactors * ( 1.0 - ExprPar::delta ) ) energyThrFactors[ i ] /= 2.0;
    }
}
