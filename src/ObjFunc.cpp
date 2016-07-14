#include <set>

#include "ObjFunc.h"

double RMSEObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double squaredErr = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
      double beta = 1.0;
      #ifdef BETAOPTTOGETHER
        if(NULL != par)
          beta = par->getBetaForSeq(i);
        squaredErr += least_square( prediction[i], ground_truth[i], beta, true );
      #else
        squaredErr += least_square( prediction[i], ground_truth[i], beta );
      #endif
  }

    double rmse = sqrt( squaredErr / ( nSeqs * nConds ) );
    return rmse;
}

double AvgCorrObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double totalSim = 0.0;

    for(int i = 0;i<ground_truth.size();i++){

      totalSim += corr(  prediction[i], ground_truth[i] );
  }

    return -totalSim/nSeqs;
  }

double PGPObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

        assert(ground_truth.size() == prediction.size());
        int nSeqs = ground_truth.size();
        int nConds = ground_truth[0].size();
        double totalPGP = 0.0;

        for(int i = 0;i<ground_truth.size();i++){
          double beta = 1.0;
          #ifdef BETAOPTTOGETHER
        	beta = par->getBetaForSeq(i);
                totalPGP += pgp(  prediction[i], ground_truth[i], beta, true);
        	#else
        	totalPGP += pgp(  prediction[i], ground_truth[i], beta );
        	#endif
      }

      return totalPGP / nSeqs;
  }

double AvgCrossCorrObjFunc::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
  {
      vector< int > shifts;
      for ( int s = -maxShift; s <= maxShift; s++ )
      {
          shifts.push_back( s );
      }

      vector< double > cov;
      vector< double > corr;
      cross_corr( x, y, shifts, cov, corr );
      double result = 0, weightSum = 0;
      //     result = corr[maxShift];
      result = *max_element( corr.begin(), corr.end() );
      //     for ( int i = 0; i < shifts.size(); i++ ) {
      //         double weight = pow( shiftPenalty, abs( shifts[i] ) );
      //         weightSum += weight;
      //         result += weight * corr[i];
      //     }
      //     result /= weightSum;

      return result;
  }

double AvgCrossCorrObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double totalSim = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
        totalSim += exprSimCrossCorr( prediction[i], ground_truth[i] );
    }

  return -totalSim / nSeqs;
  }


double LogisticRegressionObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){
    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double totalLL = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
      vector<double> Y = ground_truth[i];
      vector<double> Ypred = prediction[i];

      double one_sequence_LL = 0.0;

      for(int j = 0;j<Y.size();j++){
        double one_gt = Y[i];
        double pred_prob = logistic(w*(Ypred[j] - bias));
        double singleLL = one_gt*log(pred_prob) + (1.0 - one_gt)*log(1.0 - pred_prob);
        one_sequence_LL += singleLL;
      }

      totalLL += one_sequence_LL;
    }
    return -totalLL;
}

MultiEnhancerObjFunc::MultiEnhancerObjFunc(ObjFunc *to_wrap, vector<int> in_enhancer_promoter_mapping){
    wrapped_obj_func = to_wrap;
    enhancer_to_promoter_mapping = in_enhancer_promoter_mapping;

    enhancer_weights = vector< double >(enhancer_to_promoter_mapping.size(), 0.0);

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
}


double MultiEnhancerObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){
    //sum the predictions, creat a new prediction vector vector, and pass that to the wrapped objective function.

    vector< vector< double > > summed_predictions(ground_truth.size(),vector< double >(ground_truth[0].size(),0.0));
/*    for(int i = 0;i < ground_truth.size();i++){
      summed_predictions.push_back(vector<double>(ground_truth[i].size(),0.0));
    }
    */

    for(int i = 0;i < prediction.size();i++){
      int promoter_index = enhancer_to_promoter_mapping[i];
      for(int j = 0; j < prediction[i].size();j++){
        summed_predictions[promoter_index][j] += enhancer_weights[i]*prediction[i][j];
      }
    }

    /*beta can just handle this. need this code if we go to learning he relative weights.
    for(int i = 0; i < summed_predictions.size();i++){
      for(int j = 0; j < summed_predictions[i].size();j++){
        summed_predictions[i][j] *= inverse_number_of_enhancers_per_promoter;
      }
    }
    */

    return wrapped_obj_func->eval(ground_truth,summed_predictions,par);
}
