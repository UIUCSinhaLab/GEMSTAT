#include "FactorIntFunc.h"

string getIntOptionStr( FactorIntType intOption )
{
		if ( intOption == BINARY ) return "Binary";
		if ( intOption == GAUSSIAN ) return "Gaussian";
		if ( intOption == HELICAL ) return "Helical";

		return "Invalid";
}

double FactorIntFuncBinary::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );

    double spacingTerm = ( dist < distThr ? normalInt : 1.0 );
    double orientationTerm = orientation ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
}


double FactorIntFuncGaussian::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );

    double GaussianInt = dist < distThr ? normalInt * exp( - ( dist * dist ) / ( 2.0 * sigma * sigma ) ) : 1.0;
    return max( 1.0, GaussianInt );
}


double FactorIntFuncGeometric::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );

    double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) );
    double orientationTerm = orientation ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
}

double FactorIntFuncHelical::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
		if(dist >= distThr) return 1.0;
		double half_dna_pitch = 19.0/2.0;
		double recip_half_dna_pitch = 2.0/19.0;
		if(dist <= half_dna_pitch) return 0.0;

		double phasing = 0.5*(cos(M_PI * recip_half_dna_pitch * dist) + 1.0);


    //double orientationTerm = orientation ? 1.0 : orientationEffect;
    return phasing;
}
