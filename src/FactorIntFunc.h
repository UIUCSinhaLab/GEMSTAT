#ifndef FACTOR_INT_FUNC_H
#define FACTOR_INT_FUNC_H

#include <cassert>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

enum FactorIntType
{
    BINARY,                                       // Binary model of interaction
    GAUSSIAN,                                      // Gaussian model of interaction
    HELICAL                                       //Helical phasing model of interaction
};

FactorIntType getIntOption( const string& int_option_str );
string getIntOptionStr( FactorIntType intOption );

/*****************************************************
 * Factor-Factor Interactions
 ******************************************************/

/* FactorIntFunc class: distance-dependent function of TF-TF interaction  */
class FactorIntFunc
{
    public:
        FactorIntFunc( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr), orientationEffect( _orientationEffect) { assert( distThr >= 0 ); }
        virtual ~FactorIntFunc(){};
        // compute the factor interaction, given the normal interaction (when they are close enough)
        virtual double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const = 0;

        double getMaxDist() const
        {
            return distThr;
        }
    protected:
        double distThr;                           // if distance < thr, the "normal" value; otherwise 1 (no interaction)
        double orientationEffect;                 // the effect of orientation: if at different strands, the effect should be multiplied this value
};

class Null_FactorIntFunc : public FactorIntFunc
{
    public:
        Null_FactorIntFunc ( ) : FactorIntFunc( 0, 1.0){};
        double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const { return 1.0;}
};

/* FactorIntFuncBinary class: binary distance function */
class FactorIntFuncBinary : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncBinary( double _distThr, double _orientationEffect = 1.0 ) : FactorIntFunc( _distThr, _orientationEffect){};

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const;


};

/* FactorIntFuncGaussian class: Gaussian distance function*/
class FactorIntFuncGaussian : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncGaussian( double _distThr, double _sigma, double _orientationEffect = 1.0) : FactorIntFunc(_distThr, _orientationEffect), sigma( _sigma )
        {
            assert( sigma > 0 );
        }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const;
    private:
        double sigma;                             // standard deviation of
};

/* FactorIntFuncGeometric class: distance function decays geometrically (but never less than 1) */
class FactorIntFuncGeometric : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncGeometric( double _distThr, double _spacingEffect, double _orientationEffect ) : FactorIntFunc(_distThr, _orientationEffect), spacingEffect( _spacingEffect )
        {}

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const;

    private:
        double spacingEffect;                     // the effect of spacing
};

/* FactorIntFuncHelical class: binary distance function */
class FactorIntFuncHelical : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncHelical( double _distThr, double _orientationEffect = 1.0 ) : FactorIntFunc(_distThr, _orientationEffect)
        {}

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const;
};

/* Dimer interaction class: binary distance function */
class Dimer_FactorIntFunc : public FactorIntFunc
{
    public:
        // constructors
        Dimer_FactorIntFunc( double _distThr, bool _a_strand, bool _b_strand ) : FactorIntFunc(_distThr, 1.0)
        {
            expected_a_strand = _a_strand;
            expected_b_strand = _b_strand;
        }
        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool a_strand, bool b_strand ) const;

        bool expected_a_strand;
        bool expected_b_strand;
};

#endif
