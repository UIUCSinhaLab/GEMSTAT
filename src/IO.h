#ifndef GEMSTAT_IO_H
#define GEMSTAT_IO_H

#include <string>
#include <map>

#include "Tools.h"

int readEdgelistGraph( const string& filename, const map<string, int>& factorIdxMap, IntMatrix& destination, bool directed);

#endif
