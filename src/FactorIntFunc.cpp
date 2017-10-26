#include "FactorIntFunc.h"

#include <iostream>

#include "Tools.h"

FactorIntType getIntOption( const string& int_option_str )
{
	if ( toupperStr( int_option_str ) == "BINARY" ) return BINARY;
    if ( toupperStr( int_option_str ) == "GAUSSIAN" ) return GAUSSIAN;
    if ( toupperStr( int_option_str ) == "HELICAL" ) return HELICAL;

    cerr << "int_option_str is not a interaction option" << endl;
    exit(1);
}

string getIntOptionStr( FactorIntType intOption )
{
		if ( intOption == BINARY ) return "Binary";
		if ( intOption == GAUSSIAN ) return "Gaussian";
		if ( intOption == HELICAL ) return "Helical";

		return "Invalid";
}

double FactorIntFuncBinary::compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const
{
    assert( dist >= 0 );
	bool orientation = ( a_strand == b_strand );

    double spacingTerm = ( dist < distThr ? normalInt : 1.0 );
    double orientationTerm = orientation ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
}


double FactorIntFuncGaussian::compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const
{
    assert( dist >= 0 );
	//bool orientation = ( a_strand == b_strand );

    double GaussianInt = dist < distThr ? normalInt * exp( - ( dist * dist ) / ( 2.0 * sigma * sigma ) ) : 1.0;
    return max( 1.0, GaussianInt );
}


double FactorIntFuncGeometric::compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const
{
    assert( dist >= 0 );
	bool orientation = ( a_strand == b_strand );

    double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) );
    double orientationTerm = orientation ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
}

double FactorIntFuncHelical::compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const
{
    assert( dist >= 0 );
	//bool orientation = ( a_strand == b_strand );

    	double spacingTerm = (dist < distThr ? normalInt : 1.0);
	if(dist >= distThr) return 1.0;
	double coeff = M_PI*32.7/180.0;
	if(dist <= 5.0) return 1.0;

	double phasing = 0.5*(cos(coeff * dist) + 1.0)*(spacingTerm - 1.0) + 1.0;


    //double orientationTerm = orientation ? 1.0 : orientationEffect;
    return phasing;
}



double Dimer_FactorIntFunc::compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const
{
    assert( dist >= 0 );
	if( (a_strand != expected_a_strand) || (b_strand != expected_b_strand) || (dist > distThr) ){ return 1.0; }
	return normalInt;
}


double HalfDirectional_FactorIntFunc::compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const
{
    assert( dist >= 0 );
	if( ( enforce_a && (a_strand != expected_a_strand) ) || (enforce_b && (b_strand != expected_b_strand)) || (dist > distThr) ){ return 1.0; }
	return normalInt;
}
