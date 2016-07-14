#include <set>

#include "IO.h"

int readEdgelistGraph( const string& filename, const map<string, int>& factorIdxMap, IntMatrix& destination, bool directed){
	string factor1, factor2;

	ifstream matrix_infile(filename.c_str());
	if( !matrix_infile.is_open() ){
		return RET_ERROR;
	}

	while( matrix_infile >> factor1 >> factor2 ){
		assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
		int idx1 = factorIdxMap.at(factor1);
		int idx2 = factorIdxMap.at(factor2);

		destination.setElement( idx1, idx2, true);
		if(!directed){
			destination.setElement( idx2, idx1, true);
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
		if( factorIdxMap.at(name) != i){
			cerr << "An entry in the factor information file was out of order or otherwise invalid: " << filename << ":" << i+1 << endl;
			return RET_ERROR;
		}

		if( (actRole != 0 && actRole != 1) || (repRole != 0 && repRole !=1) ){
			cerr << "An invalid role setting was provided in the factor information file: " << filename << ":" << i+1 << endl;
			return RET_ERROR;
		}

		actIndicators[i] = (1 == actRole);//This does not rely on outside initialization of the actIndicators or repIndicators, other than instantiation.
		repIndicators[i] = (1 == repRole);
		i++;
	}

	finfo.close();

	return 0;
}

int readAxisWeights(const string& filename, vector< int >& axis_start, vector< int >& axis_end, vector< double >& axis_wts){
	//Seems to load a file of ranges and weights. All weights should total 100.
	//Ranges start at 0
	ifstream axis_wtInfo ( filename.c_str() );
        if( !axis_wtInfo )
        {
            cerr << "Cannot open the axis weight information file " << filename << endl;
            exit( 1 );
        }
        int temp1, temp2;
        double temp3;
        double temp3_sum = 0;
        while( axis_wtInfo >> temp1 >> temp2 >> temp3 ) //TODO: Change to readline to enforce end of lines
        {
            axis_start.push_back( temp1 );
            axis_end.push_back( temp2 );
            axis_wts.push_back( temp3 );
            temp3_sum += temp3;
        }
        assert( !( temp3_sum > 100 ) && !( temp3_sum < 100 ));

	return 0;
}

int writePredictions(const string& filename, ExprPredictor& predictor, Matrix& exprData, vector< string >& expr_condNames, vector< int > enhancer_to_promoter_mapping, vector< string > promoter_names, bool fix_beta /*= false*/ ){
	// re figure out the enhancer weights
	vector< double> enhancer_weights = vector< double >(enhancer_to_promoter_mapping.size(), 0.0);

	set< int > tmp_promoter_IDs;
	for(int i = 0;i<enhancer_to_promoter_mapping.size();i++){
		tmp_promoter_IDs.insert(enhancer_to_promoter_mapping[i]);
	}
	int num_promoters = tmp_promoter_IDs.size();

	vector<int> tmp_enhancer_per_promoter_count = vector<int>(num_promoters,0);
	for(int i = 0;i < enhancer_to_promoter_mapping.size();i++){
		tmp_enhancer_per_promoter_count[enhancer_to_promoter_mapping[i]] += 1;
	}

	vector< double > inverse_number_of_enhancers_per_promoter = vector<double>(num_promoters,1.0);
	for(int i = 0;i < num_promoters;i++) {
		inverse_number_of_enhancers_per_promoter[i] = 1.0/tmp_enhancer_per_promoter_count[i];
	}

	for(int i_enhancer = 0;i_enhancer < enhancer_to_promoter_mapping.size();i_enhancer++){
		enhancer_weights[i_enhancer] = inverse_number_of_enhancers_per_promoter[enhancer_to_promoter_mapping[i_enhancer]];
	}


	// print the predictions
	ofstream fout( filename.c_str() );
	if ( !fout )
	{
		cerr << "Cannot open file " << filename << endl;
		exit( 1 );
	}

	ExprPar par = predictor.getPar();
	const vector< SiteVec >& allSeqSites = predictor.getSeqSites();

	fout << "Rows\t" << expr_condNames << endl;

	vector< vector< double > > per_enhancer_predictions;

	for( int i = 0; i < predictor.nSeqs();i++){
		vector< double > targetExprs;
		predictor.predict( allSeqSites[i], predictor.seqs[i].size(), targetExprs, i );
		per_enhancer_predictions.push_back(targetExprs);
	}

	vector< vector< double > > summed_predictions(exprData.nRows(),vector< double >(exprData.nCols(),0.0));
/*    for(int i = 0;i < ground_truth.size();i++){
		summed_predictions.push_back(vector<double>(ground_truth[i].size(),0.0));
	}
	*/

	for(int i = 0;i < per_enhancer_predictions.size();i++){
		int promoter_index = enhancer_to_promoter_mapping[i];
		for(int j = 0; j < per_enhancer_predictions[i].size();j++){
			summed_predictions[promoter_index][j] += enhancer_weights[i]*per_enhancer_predictions[i][j];
		}
	}

	for ( int i = 0; i < exprData.nRows(); i++ )
    {
        vector< double > targetExprs = summed_predictions[i];
        vector< double > observedExprs = exprData.getRow( i );

        // error
        // print the results
	// observations
        fout << promoter_names[i] << "\t" << observedExprs << endl;
        fout << promoter_names[i];

        double beta = par.getBetaForSeq(i);

				double error_or_score;

				//this ensures that the right scaling parameter is applied.
				vector<vector<double> > multiple_predictions(exprData.nRows(), vector< double >(exprData.nCols(),0.0));
				vector<vector<double> > multiple_observations(exprData.nRows(), vector< double >(exprData.nCols(),0.0));
				multiple_predictions[i] = targetExprs;
				multiple_observations[i] = observedExprs;

				error_or_score = ((MultiEnhancerObjFunc*)predictor.trainingObjective)->wrapped_obj_func->eval(multiple_observations, multiple_predictions, &par);

        for ( int j = 0; j < predictor.nConds(); j++ ){
						fout << "\t" << ( beta * targetExprs[j] );
				}
        fout << endl;

        // print the agreement bewtween predictions and observations
        cout << promoter_names[i] << "\t" << beta << "\t" << error_or_score << endl;
    }
		fout.flush();
		fout.close();


		ofstream individual_fout( (filename + string(".individual")).c_str() );
		if ( !individual_fout )
		{
			cerr << "Cannot open file " << filename << ".individual" << endl;
			exit( 1 );
		}

		individual_fout << "Rows\t" << expr_condNames << endl;
		for(int i = 0;i < per_enhancer_predictions.size();i++){

				vector< double > targetExprs = per_enhancer_predictions[i];
				vector< double > observedExprs = exprData.getRow( enhancer_to_promoter_mapping[i] );

				double beta_for_enhancer = par.getBetaForSeq(enhancer_to_promoter_mapping[i]);
				individual_fout << predictor.seqs[i].getName() << "\t" << observedExprs << endl;
				individual_fout << predictor.seqs[i].getName();


				for ( int j = 0; j < predictor.nConds(); j++ ){
					individual_fout << "\t" << ( beta_for_enhancer * enhancer_weights[i]*targetExprs[j] );
				}
				individual_fout << endl;

		}
		individual_fout.flush();
		individual_fout.close();
	return 0;
}
