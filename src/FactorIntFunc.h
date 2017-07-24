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

string getIntOptionStr( FactorIntType intOption );

/*****************************************************
 * Factor-Factor Interactions
 ******************************************************/

/* FactorIntFunc class: distance-dependent function of TF-TF interaction  */
class FactorIntFunc
{
    public:
        // compute the factor interaction, given the normal interaction (when they are close enough)
        virtual double compFactorInt( double normalInt, double dist, bool orientation ) const = 0;

        // the maximum distance beyond which there is no interaction
        virtual double getMaxDist() const = 0;
};

/* FactorIntFuncBinary class: binary distance function */
class FactorIntFuncBinary : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncBinary( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyond which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // if distance < thr, the "normal" value; otherwise 1 (no interaction)
        double orientationEffect;                 // the effect of orientation: if at different strands, the effect should be multiplied this value
};

/* FactorIntFuncGaussian class: Gaussian distance function*/
class FactorIntFuncGaussian : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncGaussian( double _distThr, double _sigma ) : distThr( _distThr ), sigma( _sigma )
        {
            assert( distThr > 0 && sigma > 0 );
        }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyone which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // no interaction if distance is greater than thr.
        double sigma;                             // standard deviation of
};

/* FactorIntFuncGeometric class: distance function decays geometrically (but never less than 1) */
class FactorIntFuncGeometric : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncGeometric( double _distThr, double _spacingEffect, double _orientationEffect ) : distThr( _distThr ), spacingEffect( _spacingEffect ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyond which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // if distance < thr, the "normal" value; otherwise decay with distance (by parameter spacingEffect)
        double spacingEffect;                     // the effect of spacing
        double orientationEffect;                 // the effect of orientation: if at different strands, the effect should be multiplied this value
};

/* FactorIntFuncHelical class: binary distance function */
class FactorIntFuncHelical : public FactorIntFunc
{
    public:
        // constructors
        FactorIntFuncHelical( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

        // compute the factor interaction
        double compFactorInt( double normalInt, double dist, bool orientation ) const;

        // the maximum distance beyond which there is no interaction
        double getMaxDist() const
        {
            return distThr;
        }
    private:
        double distThr;                           // if distance < thr, the "normal" value; otherwise 1 (no interaction)
        double orientationEffect;                 // the effect of orientation: if at different strands, the effect should be multiplied this value
};

#endif
