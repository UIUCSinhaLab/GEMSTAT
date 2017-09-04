#ifndef OBJFUNC_H
#define OBJFUNC_H

#include "ExprPar.h"

class ObjFunc {
public:
  ObjFunc(){}
  virtual ~ObjFunc(){}

public:
  //pure virtual function
  virtual double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par){return -1;}
};


class RMSEObjFunc: public ObjFunc {
public:
  ~RMSEObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
};

class AvgCorrObjFunc: public ObjFunc {
public:
  ~AvgCorrObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
};

class PGPObjFunc: public ObjFunc {
public:
  ~PGPObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
};

class AvgCrossCorrObjFunc: public ObjFunc {
public:
  AvgCrossCorrObjFunc(int _maxShift, double _shiftPenalty){  maxShift = _maxShift; shiftPenalty = _shiftPenalty; }
  ~AvgCrossCorrObjFunc(){}
  double exprSimCrossCorr( const vector< double >& x, const vector< double >& y );
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
private:
  int maxShift;
  double shiftPenalty;
};

class LogisticRegressionObjFunc: public ObjFunc {
public:
  LogisticRegressionObjFunc(){ w = 1.0; bias = 0.5;}
  LogisticRegressionObjFunc(double _w, double _bias){w = _w; bias = _bias;}
  ~LogisticRegressionObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
private:
  double w;
  double bias;
};

class RegularizedObjFunc: public ObjFunc {
public:
  RegularizedObjFunc(ObjFunc* wrapped_obj_func, const ExprPar& centers, const ExprPar& l1, const ExprPar& l2);
  ~RegularizedObjFunc(){delete my_wrapped_obj_func;}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
  ObjFunc* my_wrapped_obj_func;
private:
  vector<double> my_centers;
  vector<double> lambda1;
  vector<double> lambda2;

  vector<double> cache_pars;
  //vector<double> cache_diffs;
  //vector<double> cache_sq_diffs;
};

class PeakWeightedObjFunc: public ObjFunc {
public:
  ~PeakWeightedObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
};

class Weighted_RMSEObjFunc: public RMSEObjFunc {
public:
    Weighted_RMSEObjFunc() : RMSEObjFunc(), weights(NULL) {}
  ~Weighted_RMSEObjFunc(){if(NULL != weights){delete weights;}}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);

  void set_weights(Matrix *in_weights);
private:
    Matrix *weights;
    double total_weight;
};

#endif
