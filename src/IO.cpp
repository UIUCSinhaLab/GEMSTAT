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
			destination( idx2, idx1 ) = true;
		}
	}
	
	return 0;
}

int readFactorThresholdFile( const string& filename, vector< double >& destination, int nFactors){
	
	ifstream factor_thr_input( filename.c_str() );

	if(!factor_thr_input.is_open() ){
		cerr << "Cannot open the factor threshold file " << filename << endl;
		return RET_ERROR;
	}
	
	destination.clear();

	for( int index = 0; index < nFactors; index++ )
        {
            double temp;
            factor_thr_input >> temp;
            destination.push_back( temp );
        }
	//Ensure there is only whitespace remaining in the file!
        factor_thr_input.close();

	return 0;
}

int readFactorRoleFile(const string& filename, const map<string, int>& factorIdxMap,  vector< bool>& actIndicators, vector<bool>& repIndicators){

	ifstream finfo( filename.c_str());
	if(!finfo.is_open()){
		cerr << "Cannot open the factor information file " << filename << endl;
		return RET_ERROR;
	}

	string name;
	int i = 0, actRole, repRole;
	while( finfo >> name >> actRole >> repRole )
	{
		if( factorIdxMap[name] != i){
			cerr << "An entry in the factor information file was out of order or otherwise invalid: " << filename << ":" << i+1 << endl;
			return RET_ERROR;
		}
		
		if( (actRole != 0 && actRole != 1) || (repRole != 0 && repRole !=1) ){
			cerr << "An invalid role setting was provided in the factor information file: " << filename << ":" << i+1 << end;
		}

		actIndicators[i] = (1 == actRole);//This does not rely on outside initialization of the actIndicators or repIndicators, other than instantiation.
		repIndicators[i] = (1 == repRole);
		i++;
	}
	
	return 0;
}	
