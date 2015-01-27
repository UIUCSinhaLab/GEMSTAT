#ifndef GEMSTAT_IO_H
#define GEMSTAT_IO_H

#include <string>
#include <map>

#include "Tools.h"

int readEdgelistGraph( const string& filename, const map<string, int>& factorIdxMap, IntMatrix& destination, bool directed);

int readFactorThresholdFile( const string& filename, vector< double > destination, int nFactors);

int readFactorRoleFile(const string& filename, const map<string, int>& factorIdxMap,  vector< bool>& actIndicators, vector<bool>& repIndicators);

#endif
