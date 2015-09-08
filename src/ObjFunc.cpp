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
          beta = par->betas[i];
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
        	beta = par->betas[i];
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
        double pred_prob = logistic(w*Ypred[j] - bias);
        double singleLL = one_gt*log(pred_prob) + (1.0 - one_gt)*log(1.0 - pred_prob);
        one_sequence_LL += singleLL;
      }

      totalLL += one_sequence_LL;
    }
    return -totalLL;
}
