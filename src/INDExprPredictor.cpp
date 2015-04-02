#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include "INDExprPredictor.h"

INDExprPar::INDExprPar( int _nFactors, int _nSeqs ) : ExprPar( _nFactors, _nSeqs)
{	
	cic_att = INDExprPar::default_cic_att;
}

INDExprPar::INDExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, const vector < double >& _basalTxps, const vector <double>& _pis, const vector <double>& _betas, int _nSeqs, const vector <double>& _energyThrFactors , double _cic_att) : ExprPar( _maxBindingWts, _factorIntMat, _txpEffects, _repEffects, _basalTxps, _pis, _betas, _nSeqs, _energyThrFactors)
{
	cic_att = _cic_att;
}


INDExprPar::INDExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators, int _nSeqs ) : ExprPar( pars, coopMat, actIndicators, repIndicators, _nSeqs)
{	

	cic_att = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[ pars.size() - 1 ], log( min_cic_att ), log( max_cic_att ) ) ) : exp ( pars[ pars.size() - 1 ] );

}

void INDExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    ExprPar::getFreePars( pars, coopMat, actIndicators, repIndicators );
    double cic_att_val = searchOption == CONSTRAINED ? infty_transform( log( cic_att ), log( min_cic_att ), log( max_cic_att ) ) : log( cic_att );
    	pars.push_back( cic_att_val );
}

void INDExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const
{
    ExprPar::print( os, motifNames, coopMat );
    os << cic_att;
    os << endl;

}

int INDExprPar::load( const string& file, const int num_of_coop_pairs )
{
	int first_load_results = ExprPar::load( file, num_of_coop_pairs );
	if( 0 != first_load_results )
		return first_load_results;
	//TOOD: This is just to make it compile. WE NEED THIS. I guess load should call a separate parser helper function. Load will open the filestream, and the parser functions will read from that.
	//Sub-classes should have their parser call their parent-classes parser.
	//fin >> cic_att;

    return 0;
}

void INDExprPar::adjust( const IntMatrix& coopMat  )
{
    ExprPar::adjust( coopMat );
    if( cic_att < INDExprPar::min_cic_att * ( 1.0 + ExprPar::delta ) ) cic_att *= 2.0;
    if( cic_att > INDExprPar::max_cic_att * ( 1.0 - ExprPar::delta ) ) cic_att /= 2.0;
}

double INDExprPar::min_cic_att = 0.99;
double INDExprPar::max_cic_att = 32.01;
double INDExprPar::default_cic_att = 10;

INDExprFunc::INDExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const INDExprPar& _par ) : ExprFunc( _motifs, _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, (ExprPar)_par ), par( _par)
{
	dperk_ExprData = Matrix();
	//dperk_ExprData = *(new Matrix());//TODO: I bet this leaks memory
}

double INDExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
	cerr << "This function is being invoked from objective functions for CC, w-PGP, or CrossCC. Change their code so that the function with cic_att is invoked.\n";
	exit( 1 );
	return 0;
}
double INDExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num, double _dperk_conc )
{
    double _cic_att = par.cic_att;	
    double att_mult = exp( -_cic_att * _dperk_conc );

    bindingWts.clear(); boundaries.clear();

    if( !one_qbtm_per_crm )
        seq_num = 0;
    // store the sequence
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
    	double stat_weight = par.maxBindingWts[ sites[i].factorIdx ] * 
				factorConcs[sites[i].factorIdx] * 
				sites[i].prior_probability * 
				sites[i].wtRatio;
	if( sites[i].factorIdx == 2 ){
		stat_weight *= att_mult;
	}
	assert( !(sites[i].prior_probability != 1) );
	/*
	if( sites[i].factorIdx == 3 ){
		stat_weight = 0;
	}
	if( sites[i].factorIdx == 4 ){
		double attenuating_factor = exp( -factorConcs[ 3 ] * pow( 2, 1.35 ) );
		stat_weight *= attenuating_factor;
	}*/
        bindingWts.push_back( stat_weight );
    }

    // Logistic model
    if ( modelOption == LOGISTIC ) {
        vector< double > factorOcc( motifs.size(), 0 ); // total occupancy of each factor
        for ( int i = 1; i < sites.size(); i++ ) {
            factorOcc[ sites[i].factorIdx ] += bindingWts[i] / ( 1.0 + bindingWts[i] );
        }
        double totalEffect = 0;
//         cout << "factor\toccupancy\ttxp_effect" << endl;
        for ( int i = 0; i < motifs.size(); i++ ) {
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
     //cout << "Z_off = " << Z_off << "\t";       
    double Z_on = compPartFuncOn();
     //cout << "Z_on = " << Z_on << endl;

    // compute the expression (promoter occupancy)
    double efficiency = Z_on / Z_off;
    double promoterOcc = efficiency * par.basalTxps[ seq_num ] / ( 1.0 + efficiency * par.basalTxps[ seq_num ] );
    return promoterOcc;
}

void INDExprFunc::set_dperk_expr( const Matrix& _dperk_ExprData ){
	dperk_ExprData = _dperk_ExprData;
}

//TODO : Take a chainsaw to this function (make it just call the superclass constructor.)
INDExprPredictor::INDExprPredictor( const vector <Sequence>& _seqs, const vector< SiteVec >& _seqSites, const vector < SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const Matrix& _dperk_ExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts ) :
ExprPredictor( _seqs, _seqSites, _r_seqSites, _seqLengths, _r_seqLengths, _exprData, _motifs, _factorExprData, _intFunc, _coopMat, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _indicator_bool, _motifNames, _axis_start, _axis_end, _axis_wts), dperk_ExprData( _dperk_ExprData )
{}

//TODO: Replace with a virtual member function of par.
bool INDExprPredictor::testPar( const ExprPar& par ) const
{
    if( ! ExprPredictor::testPar( par )) return false;
    if( ((INDExprPar)par).cic_att < INDExprPar::min_cic_att * ( 1.0 + INDExprPar::delta ) ) return false;
    if( ((INDExprPar)par).cic_att > INDExprPar::max_cic_att * ( 1.0 - INDExprPar::delta ) ) return false;
   
   return true;    
}

//TODO: Replace with a virtual member function of par.
void INDExprPredictor::printPar( const ExprPar& par ) const
{
    ExprPredictor::printPar( par ) ;
    cout << endl << ((INDExprPar)par).cic_att;
    cout << endl;
}

ExprFunc* INDExprPredictor::createExprFunc( const ExprPar& par ) const
{
    cerr << "DEBUG: virtual dispatch createExprFunc " << endl << flush;
    //TODO: //This here is where you want to create a factory pattern
	// Well, I guess this can _be_ the factory pattern for now.	
    INDExprFunc* retval = new INDExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par );
    retval->set_dperk_expr( dperk_ExprData );
    return retval;
}

//TODO: Use decent software engineering to make this function go away entirely.
int INDExprPredictor::predict( const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num , const ExprPar* _in_pars /* = NULL */) const
{
    cerr << "DEBUG: virtual dispatch of predict " << endl << flush;
    targetExprs.clear();
   
    const INDExprPar *my_pars = ((NULL == _in_pars) ? &(par_model) : (INDExprPar*)_in_pars);

    // create site representation of the target sequence
     SiteVec targetSites;
     SeqAnnotator ann( motifs, my_pars->energyThrFactors );    
     ann.annot( seqs[ seq_num ], targetSites );
            
    // predict the expression
    INDExprFunc* func = (INDExprFunc*)createExprFunc( *(my_pars) );       
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        vector <double> t_dperk_conc = dperk_ExprData.getCol( j );
        double _dperk_conc = t_dperk_conc[0];
        double predicted = func->predictExpr( targetSites, targetSeqLength, concs, seq_num, _dperk_conc);
        targetExprs.push_back( predicted );
    }

    return 0;
}
