
#include "DataSet.h"

Condition::Condition(vector< double > _concs) : concs(_concs){
}

DataSet::DataSet(const Matrix& tf_concentrations, const Matrix& output_values) : exprData(output_values), factorExprData(tf_concentrations)
{
  assert(exprData.nCols() == factorExprData.nCols());
}

int DataSet::nConds() const{
  return factorExprData.nCols();
}

int DataSet::n_rows_output() const{
	return exprData.nRows();
}
int DataSet::n_rows_input() const{
	return factorExprData.nRows();
}

Condition DataSet::getCondition(int i, const ExprPar &signalling_params) const{
  vector< double > values = factorExprData.getCol(i);
  return Condition(values);

}

vector< double > DataSet::get_output_row(int i) const{return this->exprData.getRow(i);}
vector< double > DataSet::get_output_col(int j) const{return this->exprData.getCol(j);}
const Matrix& DataSet::get_output_matrix()const{return this->exprData;}
