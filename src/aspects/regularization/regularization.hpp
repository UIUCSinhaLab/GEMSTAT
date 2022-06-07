#ifndef REGULARIZATION_HPP
#define REGULARIZATION_HPP

#include "ObjFunc.h"

#include "ExprPredictor.h"

/**
Decorator/container objective function, it wraps other objective functions.
*/
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

void regularization_cmndline_init(ExprPredictor *predictor, int argc, char* argv[] );

#endif //REGULARIZATION_HPP
