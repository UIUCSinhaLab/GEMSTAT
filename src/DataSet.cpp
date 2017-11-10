
#include "DataSet.h"

DataSet::DataSet(const Matrix& tf_concentrations, const Matrix& output_values) : exprData(output_values), factorExprData(tf_concentrations)
{
  assert(exprData.nCols() == factorExprData.nCols());
}

int DataSet::nConds() const{
  return factorExprData.nCols();
}

Condition DataSet::getCondition(int i, ExprPar signalling_params) const{
  vector< double > values = factorExprData.getCol(i);
  return Condition(values);

}

Condition::Condition(vector< double > _concs) : concs(_concs){

}
