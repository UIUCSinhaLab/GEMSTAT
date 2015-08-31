#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include "ExprPredictor.h"
#include "conf/ExprParConf.hpp"

ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    if ( toupperStr( objOptionStr ) == "PGP" ) return PGP;

    cerr << "objOptionStr is not a valid option of objective function" << endl;
    exit(1);
}


string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    if ( objOption == PGP ) return "PGP";

    return "Invalid";
}


string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}


ExprPar::ExprPar( int _nFactors, int _nSeqs ) : factorIntMat()
{
    assert( _nFactors > 0 );

    for ( int i = 0; i < _nFactors; i++ )
    {
        maxBindingWts.push_back( ExprPar::default_weight );
    }

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );

    for ( int i = 0; i < _nFactors; i++ )
    {
        double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        txpEffects.push_back( defaultEffect );
        repEffects.push_back( ExprPar::default_repression );
    }

    nSeqs = _nSeqs;
    if( one_qbtm_per_crm  )
    {
        for( int i = 0; i < nSeqs; i++ )
        {
            double basalTxp_val = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo;
            basalTxps.push_back( basalTxp_val );
        }
    }
    else
    {
        double basalTxp_val = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo;
        basalTxps.push_back( basalTxp_val );
    }
    //for the pausing parameters

    for( int i = 0; i < nSeqs; i++ )
    {
        pis.push_back( ExprPar::default_pi );
    }
    //for the beta parameters
    for( int i = 0; i < nSeqs; i++ )
    {
        betas.push_back( ExprPar::default_beta );
    }
    for( int i = 0; i < _nFactors; i++ )
    {
        energyThrFactors.push_back( ExprPar::default_energyThrFactors );
    }
}


ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector <double>& _energyThrFactors ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxps( _basalTxps ), pis( _pis), betas( _betas ), nSeqs( _nSeqs ), energyThrFactors( _energyThrFactors )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() );
    assert( txpEffects.size() == maxBindingWts.size() && repEffects.size() == maxBindingWts.size() );
    if ( one_qbtm_per_crm )
    {
        assert( basalTxps.size() == nSeqs );
    }
    else
    {
        assert( basalTxps.size() == 1 );
    }
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

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : motifs( _motifs ), intFunc( _intFunc ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), par( _par )
{
    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
}


double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
    bindingWts.clear(); boundaries.clear();

    if( !one_qbtm_per_crm )
        seq_num = 0;
    // store the sequence
    int n = _sites.size();
    sites = _sites;
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

    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ )
    {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio );
        double samee = par.maxBindingWts[sites[i].factorIdx]*factorConcs[sites[i].factorIdx]*sites[i].prior_probability*sites[i].wtRatio;
        if(samee != samee)
        {
            cout << "DEBUG: samee for " << i << "\t" << sites[i].factorIdx << "\t" <<  par.maxBindingWts[sites[i].factorIdx] <<"\t" << factorConcs[sites[i].factorIdx] <<"\t" << sites[i].prior_probability <<"\t" << sites[i].wtRatio << endl;
            exit(1);
        }
    }

    // Logistic model
    if ( modelOption == LOGISTIC )
    {
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

ModelType ExprFunc::modelOption = QUENCHING;

double ExprFunc::compPartFuncOff() const
{
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) return compPartFuncOffChrMod();

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


double ExprFunc::compPartFuncOffChrMod() const
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
    if ( modelOption == DIRECT ) return compPartFuncOnDirect();
    if ( modelOption == QUENCHING ) return compPartFuncOnQuenching();
    if ( modelOption == CHRMOD_UNLIMITED) return compPartFuncOnChrMod_Unlimited();
    if ( modelOption == CHRMOD_LIMITED ) return compPartFuncOnChrMod_Limited();
//TODO: A compiler warning is generated here. Shouldn't there be some defensive coding?
}


double ExprFunc::compPartFuncOnDirect() const
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


double ExprFunc::compPartFuncOnQuenching() const
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


double ExprFunc::compPartFuncOnChrMod_Unlimited() const
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


double ExprFunc::compPartFuncOnChrMod_Limited() const
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


ExprPredictor::ExprPredictor( const vector <Sequence>& _seqs, const vector< SiteVec >& _seqSites, const vector < SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const ExprModel& _expr_model,
		const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts ) : seqs(_seqs), seqSites( _seqSites ), r_seqSites( _r_seqSites ), seqLengths( _seqLengths ), r_seqLengths( _r_seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ),
	expr_model( _expr_model),
	indicator_bool ( _indicator_bool ), motifNames ( _motifNames ), axis_start ( _axis_start ), axis_end( _axis_end ), axis_wts( _axis_wts )
{
    //TODO: Move appropriate lines from this block to the ExprModel class.
    assert( exprData.nRows() == nSeqs() );
    assert( factorExprData.nRows() == nFactors() && factorExprData.nCols() == nConds() );
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

    // set the option of parameter estimation
    ExprPar::estBindingOption = estBindingOption;
}


double ExprPredictor::objFunc( const ExprPar& par )
{
    if ( objOption == SSE ) return compRMSE( par );
    if ( objOption == CORR ) return -compAvgCorr( par );
    if ( objOption == PGP ) return compPGP( par );
    if ( objOption == CROSS_CORR ) return -compAvgCrossCorr( par );
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

    if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ) par_model.adjust( expr_model.coopMat );
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
        par_model.adjust( expr_model.coopMat );
	#ifdef BETAOPTSEPARATE
        optimize_beta( par_model, obj_result );
        #endif
        gradient_minimize( par_result, obj_result );
        par_model = par_result;
        par_model.adjust( expr_model.coopMat);
	#ifdef BETAOPTSEPARATE
        optimize_beta( par_model, obj_result );
        #endif
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
        randSamplePar( rng, par_rand_start );
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
        randSamplePar( rng, par_curr );
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


int ExprPredictor::predict( const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const
{
    targetExprs.clear();

    // create site representation of the target sequence
    SiteVec targetSites;
    SeqAnnotator ann( motifs, par_model.energyThrFactors );
    ann.annot( seqs[ seq_num ], targetSites );

    // predict the expression
    ExprFunc* func = createExprFunc( par_model );
    for ( int j = 0; j < nConds(); j++ )
    {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr( targetSites, targetSeqLength, concs, seq_num );
        targetExprs.push_back( predicted );
    }

    delete func;
    return 0;
}


// double ExprPredictor::test( const vector< Sequence >& testSeqs, const Matrix& testExprData, int perfOption ) const
// {
// 	assert( perfOption == 0 );
// 	assert( testExprData.nRows() == testSeqs.size() && testExprData.nCols() == nConds() );
//
// 	// make predictions
// 	Matrix predicted( testSeqs.size(), nConds() );
// 	for ( int i = 0; i < testSeqs.size(); i++ ) {
//         vector< double > targetExprs;
//         predict( targetSeqs[i], targetExprs );
//         predicted.setRow( i, targetExprs );
// 	}
//
// 	// RMSE between predictions and observations
// 	if ( perfOption == 0 ) {
// 		vector< double > corrs;
// 		for ( int i = 0; i < nExps; i++ ) {
// 			corrs.push_back( correlation( predicted[ i ], testExprData[ i ] ) );
// 		}
//
// 		return mean( corrs );
// 	}
// }

int ExprPredictor::estBindingOption = 1;          // 1. estimate binding parameters; 0. not estimate binding parameters
ObjType ExprPredictor::objOption = SSE;

double ExprPredictor::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
{
    vector< int > shifts;
    for ( int s = -maxShift; s <= maxShift; s++ )
    {
        shifts.push_back( s );
    }

    vector< double > cov;
    vector< double > corr;
    cross_corr( x, y, shifts, cov, corr );
    double result = 0, weightSum = 0;
    //     result = corr[maxShift];
    result = *max_element( corr.begin(), corr.end() );
    //     for ( int i = 0; i < shifts.size(); i++ ) {
    //         double weight = pow( shiftPenalty, abs( shifts[i] ) );
    //         weightSum += weight;
    //         result += weight * corr[i];
    //     }
    //     result /= weightSum;

    return result;
}


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

int ExprPredictor::randSamplePar( const gsl_rng* rng, ExprPar& par ) const
{
    // sample binding weights
    if ( estBindingOption )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            double rand_weight = exp( gsl_ran_flat( rng, log( ExprPar::min_weight ), log( ExprPar::max_weight ) ) );
            par.maxBindingWts[i] = rand_weight;
        }
    }

    // sample the interaction matrix
    if ( expr_model.modelOption != LOGISTIC )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            for ( int j = 0; j <= i; j++ )
            {
                double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interaction ), log( ExprPar::max_interaction ) ) );
                if ( expr_model.coopMat( i, j ) ) par.factorIntMat( i, j ) = rand_interaction;
            }
        }
        for ( int i = 0; i < nFactors(); i++ )
        {
            for ( int j = i + 1; j < nFactors(); j++ )
            {
                par.factorIntMat( i, j ) = par.factorIntMat( j, i );
            }
        }
    }

    // sample the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( expr_model.modelOption == LOGISTIC )
        {
            double rand_effect = gsl_ran_flat( rng, ExprPar::min_effect_Logistic, ExprPar::max_effect_Logistic );
            par.txpEffects[i] = rand_effect;
        }
        else if ( expr_model.modelOption == DIRECT )
        {
            double rand_effect = exp( gsl_ran_flat( rng, log( ExprPar::min_effect_Thermo ), log( ExprPar::max_effect_Thermo ) ) );
            par.txpEffects[i] = rand_effect;
        }
        else
        {
            if ( expr_model.actIndicators[i] )
            {
                double rand_effect = exp( gsl_ran_flat( rng, log( ExprPar::min_effect_Thermo ), log( ExprPar::max_effect_Thermo ) ) );
                par.txpEffects[i] = rand_effect;
            }
        }
    }

    // sample the repression effects
    if ( expr_model.modelOption == CHRMOD_UNLIMITED || expr_model.modelOption == CHRMOD_LIMITED )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( expr_model.repIndicators[i] )
            {
                double rand_repression = exp( gsl_ran_flat( rng, log( ExprPar::min_repression ), log( ExprPar::max_repression ) ) );
                par.repEffects[i] = rand_repression;
            }
        }
    }

    // sample the basal transcription
    double rand_basal;
    for( int _i = 0; _i < par.basalTxps.size(); _i ++ )
    {
        if ( expr_model.modelOption == LOGISTIC )
            rand_basal = gsl_ran_flat( rng, ExprPar::min_basal_Logistic, ExprPar::max_basal_Logistic );
        else
            rand_basal = exp( gsl_ran_flat( rng, log( ExprPar::min_basal_Thermo ), log( ExprPar::max_basal_Thermo ) ) );
        par.basalTxps[ _i ] = rand_basal;
    }
    return 0;
}


bool ExprPredictor::testPar( const ExprPar& par ) const
{
    // test binding weights
    if ( estBindingOption )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( par.maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) || par.maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) )
                return false;
        }
    }

    //cout << "dbg1" << endl;

    // test the interaction matrix
    if ( expr_model.modelOption != LOGISTIC )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            for ( int j = 0; j <= i; j++ )
            {
                if ( par.factorIntMat( i, j ) < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) || par.factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    // test the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( expr_model.modelOption == LOGISTIC )
        {
            if ( par.txpEffects[i] < ExprPar::min_effect_Logistic + ExprPar::delta || par.txpEffects[i] > ExprPar::max_effect_Logistic - ExprPar::delta )
                return false;
        }
        else if ( expr_model.modelOption == DIRECT )
        {
            if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        }
        else
        {
            if ( expr_model.actIndicators[i] )
            {
                if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    //cout << "dbg2" << endl;
    // test the repression effects
    if ( expr_model.modelOption == CHRMOD_UNLIMITED || expr_model.modelOption == CHRMOD_LIMITED || expr_model.modelOption == DIRECT )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( expr_model.repIndicators[i] )
            {
                if ( par.repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) || par.repEffects[i] > ExprPar::max_repression * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    //cout << "dbg3" << endl;
    // test the basal transcription
    for( int _i = 0; _i < par.basalTxps.size(); _i++ )
    {
        if ( expr_model.modelOption == LOGISTIC )
        {
            if ( par.basalTxps[ _i ] < ExprPar::min_basal_Logistic + ExprPar::delta || par.basalTxps[ _i ] > ExprPar::max_basal_Logistic - ExprPar::delta )
                return false;
        }
        else
        {
            if ( par.basalTxps[ _i ] < ExprPar::min_basal_Thermo * ( 1.0 + ExprPar::delta ) || par.basalTxps[ _i ] > ExprPar::max_basal_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        }
    }
    //test the pi values

    for( int _i = 0; _i < par.pis.size(); _i++ )
    {
        if( par.pis[ _i ] < 0 || par.pis[ _i ] > ExprPar::max_pi * ( 1.0 - ExprPar::delta ))
        {
            return false;
        }
    }
    //cout << "dbg4" << endl;
    //test the beta values
    for( int _i = 0; _i < par.betas.size(); _i++ )
    {
        if ( par.betas[ _i ] < ExprPar::min_beta * ( 1.0 + ExprPar::delta ) || par.betas[ _i ] > ExprPar::max_beta * ( 1.0 - ExprPar::delta ) )
            return false;
    }
    //cout << "dbg5" << endl;
    //adjust the energyThrFactors
    for( int i = 0; i < nFactors(); i++ )
    {
        if( par.energyThrFactors[ i ] < ExprPar::min_energyThrFactors * ( 1.0 + ExprPar::delta ) ) return false;
        if( par.energyThrFactors[ i ] > ExprPar::max_energyThrFactors * ( 1.0 - ExprPar::delta ) ) return false;
    }
    return true;
}


void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 50 );
    //     cout.width( 8 );

    // print binding weights

    for ( int i = 0; i < nFactors(); i++ )
    {
          cout << par.maxBindingWts[i] << "\t";
    }

    // print the interaction matrix
    for ( int i = 0; i < nFactors(); i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
            if ( expr_model.coopMat( i, j ) ) cout << par.factorIntMat( i, j ) << "\t";
        }
    }

    // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        cout << par.txpEffects[i];
    }

    // print the repression effects
    if ( expr_model.modelOption == CHRMOD_UNLIMITED || expr_model.modelOption == CHRMOD_LIMITED || expr_model.modelOption == DIRECT )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            cout << par.repEffects[i] << "\t";
        }
    }

    // print the basal transcriptions
    for( int _i = 0; _i < par.basalTxps.size(); _i ++ )
    {
        cout << par.basalTxps[ _i ] << "\t";
    }
    cout << endl;
    //print the pi values
    for( int _i = 0; _i < par.pis.size(); _i++ )
    {
        cout << par.pis[ _i ] << "\t";
    }
    cout << endl;
    //print the beta values
    assert( par.betas.size() == nSeqs() );
    /*for( int i = 0; i < par.betas.size(); i++ ){
        cout << par.betas[ i ] << endl;
    }*/
    for( int _i = 0; _i < par.energyThrFactors.size(); _i++ )
    {
        cout << par.energyThrFactors[ _i ] << "\t";
    }
    cout << endl;
}


ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{
	//TODO: just make it take an expr_model as a parameter.
	//Also, it should use a factory to create the ExprFunc, but this is already close to a factory method.
    return new ExprFunc( motifs, expr_model.intFunc, expr_model.actIndicators, expr_model.maxContact, expr_model.repIndicators, expr_model.repressionMat, expr_model.repressionDistThr, par );
}


int indices_of_crm_in_gene[] =
{
};

double ExprPredictor::compRMSE( const ExprPar& par )
{
    ExprFunc* func = createExprFunc( par );
    vector< SiteVec > seqSites( seqs.size() );
    vector< int > seqLengths( seqs.size() );
    SeqAnnotator ann( motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ ) {
       	ann.annot( seqs[ i ], seqSites[ i ] );
    	seqLengths[i] = seqs[i].size();
    }
    // error of each sequence
    double squaredErr = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
            	predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i );


            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );

            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        double beta;
	#ifdef BETAOPTTOGETHER
	beta = par.betas[i];
        squaredErr += least_square( predictedExprs, observedExprs, beta, true );
	#else
        squaredErr += least_square( predictedExprs, observedExprs, beta );
	#endif
	//(func -> getPar()).betas[ i ] = beta;
    }

    double rmse = sqrt( squaredErr / ( nSeqs() * nConds() ) );//use this line if least_square() is used
	double penalty = 0;

	//TODO: R_SEQ Either remove this dead code, or make this a conditional option.
	//for the random sequences: start
        /*double sum_max_rand = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
    	double max = 0;
        for ( int j = 0; j < nConds(); j++ ) {
		//cout << "inside for cond: " << j + 1 << endl;
            	vector< double > concs = factorExprData.getCol( i * nConds() +  j );
            	double predicted = func->predictExpr( r_seqSites[ i ], r_seqLengths[i], concs, i );
		if( predicted > max ){
			max = predicted;
		}
        }
	sum_max_rand += max ;
    }
	double avg_max_rand = sum_max_rand/nSeqs();
	if( avg_max_rand > 0.01 )
		penalty = -1000;

//cout << "end looking at RSs:" << endl;
*/
	//for the random sequences: end

    return rmse - penalty;
}

double ExprPredictor::compAvgCorr( const ExprPar& par )
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );
    vector< SiteVec > seqSites( seqs.size() );
    vector< int > seqLengths( seqs.size() );
    SeqAnnotator ann( motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ )
    {
        ann.annot( seqs[ i ], seqSites[ i ] );
        seqLengths[i] = seqs[i].size();
    }

    // Pearson correlation of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ )
    {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ )
        {
            double predicted = -1;
            vector< double > concs = factorExprData.getCol( j );
            predicted = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs, i );

            if( predicted != predicted )
            {
                cout << "DEBUG: nan in predicted value" << endl;
                exit(1);
            }

            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
	//BETAOPTOGETHER has no meaning for correlation
        totalSim += corr( predictedExprs, observedExprs );
        //         cout << "Sequence " << i << "\t" << corr( predictedExprs, observedExprs ) << endl;
    }
    delete func;
    return totalSim / nSeqs();
}


double ExprPredictor::compPGP( const ExprPar& par )
{

    /*cout << "start samee debug" << endl;
    printPar( par );
    cout << "end samee debug" << endl;	*/
    ExprFunc* func = createExprFunc( par );
    vector< SiteVec > seqSites( seqs.size() );
    vector< int > seqLengths( seqs.size() );
    SeqAnnotator ann( motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ )
    {
        ann.annot( seqs[ i ], seqSites[ i ] );
        seqLengths[i] = seqs[i].size();
    }
    // PGP of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ )
    {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ )
        {
            double predicted = -1;
            vector< double > concs = factorExprData.getCol( j );
            predicted = func->predictExpr( seqSites[ i ], seqLengths[ i ], concs, i );

            if( predicted != predicted )
            {
                cout << "DEBUG: nan in predicted value" << endl;
                exit(1);
            }
            predictedExprs.push_back( predicted );
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        double beta;
	#ifdef BETAOPTTOGETHER
	beta = par.betas[i];
        totalSim += pgp( predictedExprs, observedExprs, beta, true);
	#else
	totalSim += pgp( predictedExprs, observedExprs, beta );
	#endif
    }
    delete func;
    return totalSim / nSeqs();
}


double ExprPredictor::compAvgCrossCorr( const ExprPar& par )
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );

    // cross correlation similarity of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ )
    {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ )
        {
            double predicted = -1;
            vector< double > concs = factorExprData.getCol( j );
            predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i );

            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );

            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += exprSimCrossCorr( predictedExprs, observedExprs );
    }

    delete func;
    return totalSim / nSeqs();
}


int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector < double > pars;
    /*debug*/
    //for(int i = 0; i < actIndicators.size(); i++) cout << "act: " << i << " " << (actIndicators[i] > 0? 1 : 0) << endl;
    //for(int i = 0; i < repIndicators.size(); i++) cout << "rep: " << i << " " << (repIndicators[i] > 0? 1 : 0) << endl;
    //cout << "DEBUG: in getFreePars()" << endl;
    par_model.getFreePars( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators );
    //cout << "pars.size() = " << pars.size() << endl;
    //cout << "DEBUG: out getFreePars()" << endl;
    //Hassan start:
    int pars_size = pars.size();
    /*DEBUG*///cout << "DEBUG: " << pars_size << endl;
    /*DEBUG*///exit(1);
    fix_pars.clear();
    free_pars.clear();
    for( int index = 0; index < pars_size; index++ )
    {
        if( indicator_bool[ index ] )
        {
            free_pars.push_back( pars[ index ]);
        }
        else
        {
            fix_pars.push_back( pars[ index ] );
        }
    }
    //cout << "DEBUG: free_pars.size() = " << free_pars.size() << "\t" << "fix_pars.size() = " << fix_pars.size() << endl;
    pars.clear();
    pars = free_pars;

    //Hassan end
    // set the objective function
    gsl_multimin_function my_func;
    my_func.f = &gsl_obj_f;
    my_func.n = pars.size();
    my_func.params = (void*)this;

    // set the initial values to be searched
    gsl_vector* x = vector2gsl( pars );
    //     for ( int i = 0; i < v.size(); i++ ) gsl_vector_set( x, i, v[ i ] );

    // CHECK POINT: evaluate gsl_obj_f() function
    // 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl;

    // choose the method of optimization and set its parameters
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_vector* ss = gsl_vector_alloc( my_func.n );
    gsl_vector_set_all( ss, 1.0 );

    // create the minimizer
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );
    gsl_multimin_fminimizer_set( s, &my_func, x, ss );

    // iteration
    size_t iter = 0;
    int status;
    double size;
    do
    {
        double f_prev = iter ? s->fval : 1.0E6;   // the function starts with some very large number

        iter++;
        status = gsl_multimin_fminimizer_iterate( s );

        // check for error
        if ( status ) break;
        //Hassan start:
        free_pars = gsl2vector( s-> x);
        pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < pars_size; index ++ )
        {
            if( indicator_bool[ index ] )
            {
                pars.push_back( free_pars[ free_par_counter ++ ]);
            }
            else
            {
                pars.push_back( fix_pars[ fix_par_counter ++ ]);
            }
        }
        //cout << "DEBUG: pars.size() = " << pars.size() << endl;
        //cout << "DEBUG: free_pars.size() = " << free_pars.size() << "\t" << "fix_pars.size() = " << fix_pars.size() << "\t";
        //cout << "free_par_counter = " << free_par_counter << "\t" << "fix_par_counter = " << fix_par_counter << endl;
        //cout << "samee: " << fix_par_counter << endl;
        //cout << "DEBUG: init ExprPar start" << endl;
        ExprPar par_curr = ExprPar ( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators, nSeqs() );
        //cout << "pars.size() = " << pars.size() << "\tpars_size = " << pars_size << endl;
        //printPar( par_curr );
        //cout << "DEBUG: init ExprPar end" << endl;
        //Hassan end
        // check if the current values of parameters are valid
        //the following line should be uncommented if you remove all the changes by Hassan
        //ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );

        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) break;

        // check for stopping condition
        //         double f_curr = s->fval;
        //         double delta_f = abs( f_curr - f_prev );
        //         if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
        //         if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
        //         if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;

        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size, 1e-1 );
        // 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

        // print the current parameter and function values
        /*cout << "======================================" << endl;
        cout << "======================================" << endl;*/
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
        //par_curr.print( cout, motifNames, coopMat);
        //cout << "##" << endl << endl;
    } while ( status == GSL_CONTINUE && iter < nSimplexIters );

    // get the results
    //     vector< double > expv;
    //     for ( int i = 0; i < ( s->x )->size; i++ ) expv.push_back( exp( gsl_vector_get( s->x, i ) ) );
    //Hassan start:
    free_pars = gsl2vector( s-> x);
    pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < pars_size; index ++ )
    {
        if( indicator_bool[ index ] )
        {
            pars.push_back( free_pars[ free_par_counter ++ ]);
        }
        else
        {
            pars.push_back( fix_pars[ fix_par_counter ++ ]);
        }
    }
    //cout << "DEBUG: init par_result start." << endl;
    par_result = ExprPar ( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators, nSeqs() );
    printPar( par_result );
    //cout << "DEBUG: init par_result end." << endl;
    //Hassan end
    //uncomment the following line if you remove all the changes by Hassan
    //par_result = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
    obj_result = s->fval;

    // free the minimizer
    gsl_vector_free( x );
    gsl_vector_free( ss );
    gsl_multimin_fminimizer_free( s );

    return 0;
}


int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    //cout << "DEBUG: in getFreePars()" << endl;
    par_model.getFreePars( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators );
    //cout << "DEBUG: out getFreePars()" << endl;

    //Hassan start:
    int pars_size = pars.size();
    fix_pars.clear();
    free_pars.clear();
    for( int index = 0; index < pars_size; index++ )
    {
        if( indicator_bool[ index ] )
        {
            //cout << "testing 1: " << pars[ index ] << endl;
            free_pars.push_back( pars[ index ]);
        }
        else
        {
            //cout << "testing 2: " << pars[ index ] << endl;
            fix_pars.push_back( pars[ index ] );
        }
    }

    pars.clear();
    pars = free_pars;
    //Hassan end
    // set the objective function and its gradient
    gsl_multimin_function_fdf my_func;
    my_func.f = &gsl_obj_f;
    my_func.df = &gsl_obj_df;
    my_func.fdf = &gsl_obj_fdf;
    my_func.n = pars.size();
    my_func.params = (void*)this;

    // set the initial values to be searched
    gsl_vector* x = vector2gsl( pars );

    // CHECK POINT: evaluate gsl_obj_f() function
    // 	cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl;

    // choose the method of optimization and set its parameters
    // 	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;
    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs;

    // create the minimizer
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
    double init_step = 0.02, tol = 0.1;
    gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );

    // iteration
    size_t iter = 0;
    int status;
    do
    {
        double f_prev = iter ? s->f : 1.0E6;      // the function starts with some very large number

        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
        // 	    if ( prev_f - curr_f < 0.001 ) break;

        // check for error
        if ( status ) break;

        // check if the current values of parameters are valid
        //Hassan start:
        free_pars = gsl2vector( s-> x);
        pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < pars_size; index ++ )
        {
            if( indicator_bool[ index ] )
            {
                pars.push_back( free_pars[ free_par_counter ++ ]);
            }
            else
            {
                pars.push_back( fix_pars[ fix_par_counter ++ ]);
            }
        }

        ExprPar par_curr = ExprPar ( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators, nSeqs() );
        //Hassan end
        //ExprPar par_curr = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ){
		/* TODO: Fix the way exception handlign happens here. See issue #28 on github.
		 * If this exception happens, the best parameter vector seen so far is not the one that gets saved, maybe it reverts to before the whole run of optimization?
		 * However, the objective function value of the rejected parameter vector _does_ get saved.
		 */
	       	break;
	}

        // check for stopping condition
        double f_curr = s->f;
        double delta_f = abs( f_curr - f_prev );
        if ( objOption == SSE && delta_f < min_delta_f_SSE ) break;
        if ( objOption == CORR && delta_f < min_delta_f_Corr ) break;
        if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
        if ( objOption == PGP && delta_f < min_delta_f_PGP ) break;

        status = gsl_multimin_test_gradient( s->gradient, 5e-4 );
        // 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

        // print the current parameter and function values
        /*cout << "========================================" << endl;
        cout << "========================================" << endl;*/
        cout << iter << "\t";
        printPar( par_curr );
        printf( "\tf() = %8.5f\n", s->f );
        //par_curr.print( cout, motifNames, coopMat);
        //cout << "##" << endl << endl;
    } while ( status == GSL_CONTINUE && iter < nGradientIters );

    // get the results
    //     vector< double > expv;
    //     for ( int i = 0; i < ( s->x )->size; i++ ) expv.push_back( exp( gsl_vector_get( s->x, i ) ) );
    //Hassan start:
    free_pars = gsl2vector( s-> x);
    pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < pars_size; index ++ )
    {
        if( indicator_bool[ index ] )
        {
            pars.push_back( free_pars[ free_par_counter ++ ]);
        }
        else
        {
            pars.push_back( fix_pars[ fix_par_counter ++ ]);
        }
    }

    par_result = ExprPar ( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators, nSeqs() );
    //Hassan end
    //par_result = ExprPar( gsl2vector( s->x ), coopMat, actIndicators, repIndicators );
    obj_result = s->f;

    // free the minimizer
    gsl_vector_free( x );
    gsl_multimin_fdfminimizer_free( s );

    return 0;
}

int ExprPredictor::optimize_beta( ExprPar& par_result, double &obj_result ){
        SiteVec unuzedSV = SiteVec();

        for ( int i = 0; i < nSeqs(); i++ ){
                vector< double > targetExprs;
                vector< double > observedExprs = exprData.getRow( i );

                predict( unuzedSV, seqs[i].size(), targetExprs, i );

                double beta = par_model.betas[ i ];
                double error_or_score;

                if( ExprPredictor::objOption == SSE )
                {
                        error_or_score = sqrt( least_square( targetExprs, observedExprs, beta , false) / nConds() );
                }
                if( ExprPredictor::objOption == PGP ){
                        error_or_score = pgp( targetExprs, observedExprs, beta, false );
                }
                if( ExprPredictor::objOption == CORR ){
                        error_or_score = corr( targetExprs, observedExprs, beta, false);
                }
                if( ExprPredictor::objOption == CROSS_CORR ){
                        error_or_score = ExprPredictor::exprSimCrossCorr( observedExprs, targetExprs );
                }

                //TODO: Check for the appropriate indicator value to see if this beta is an optimization parameter or not.
                par_result.betas[i] = beta;

        }

	obj_result = objFunc( par_result );
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
    all_pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < predictor -> indicator_bool.size(); index ++ )
    {
        if( predictor ->  indicator_bool[ index ]  )
        {
            all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
        }
        else
        {
            all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
        }
    }

    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(), predictor -> nSeqs() );
    //ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );

    // call the ExprPredictor object to evaluate the objective function
    double obj = predictor->objFunc( par );
    return obj;
}

void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    double step = 1.0E-3;
    numeric_deriv( grad, gsl_obj_f, v, params, step );
}


void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params );
    gsl_obj_df( v, params, grad );
}
