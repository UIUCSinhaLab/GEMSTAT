#include "IO.h"

int readEdgelistGraph( const string& filename, const map<string, int>& factorIdxMap, IntMatrix& destination, bool directed){
	string factor1, factor2;

	ifstream matrix_infile(filename.c_str());
	if( !matrix_infile ){
		return 1;
	}

	while( matrix_infile >> factor1 >> factor2 ){
		assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
		int idx1 = factorIdxMap.at(factor1);
		int idx2 = factorIdxMap.at(factor2);

		destination( idx1, idx2 ) = true;
		if(!directed){
			destination( idx2, idx1 ) = false;
		}
	}
	
	return 0;
}
