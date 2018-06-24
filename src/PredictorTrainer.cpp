/*
 * PredictorTrainer.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: lunt
 */

#include "PredictorTrainer.h"

ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    if ( toupperStr( objOptionStr ) == "PGP" ) return PGP;
    if ( toupperStr( objOptionStr ) == "LOGISTIC_REGRESSION") return LOGISTIC_REGRESSION;
    if ( toupperStr( objOptionStr ) == "PEAK_WEIGHTED") return PEAK_WEIGHTED;
    if ( toupperStr( objOptionStr ) == "WEIGHTED_SSE") return WEIGHTED_SSE;


    cerr << "objOptionStr is not a valid option of objective function" << endl;
    exit(1);
}


string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    if ( objOption == PGP ) return "PGP";
    if ( objOption == LOGISTIC_REGRESSION ) return "LOGISTIC_REGRESSION";
    if ( objOption == PEAK_WEIGHTED ) return "PEAK_WEIGHTED";
    if ( objOption == WEIGHTED_SSE ) return "WEIGHTED_SSE";

    return "Invalid";
}


string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}
