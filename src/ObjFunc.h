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

class MultiEnhancerObjFunc: public ObjFunc {
public:
  MultiEnhancerObjFunc(ObjFunc *to_wrap, vector<int> in_enhancer_promoter_mapping);
  ~MultiEnhancerObjFunc(){delete wrapped_obj_func;}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);

  ObjFunc *wrapped_obj_func;
  vector<int> enhancer_to_promoter_mapping;
private:
  /* won't need this quite yet. Beta will handle this.
  vector<double> inverse_number_of_enhancers_per_promoter;
  */
};

#endif
