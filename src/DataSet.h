/*
 * DataSet.h
 *
 *  Created on: April 26, 2017
 *      Author: lunt
 */

#ifndef SRC_DATASET_H_
#define SRC_DATASET_H_

#include "Tools.h"
#include "ExprPar.h"

class Condition{
  public:
    Condition(vector< double > _concs);
    vector< double > concs;
};

class DataSet{
public:
  DataSet(const Matrix& tf_concentrations, const Matrix& output_values);
  DataSet(const DataSet &other) : exprData(other.exprData), factorExprData(other.factorExprData){};
  ~DataSet(){};

  virtual int nConds() const;//number of columns
  virtual int n_rows_output() const;
  virtual int n_rows_input() const;

  virtual Condition getCondition(int i , const ExprPar &signalling_params) const;
  virtual vector< double > get_output_row(int i) const;
  virtual vector< double > get_output_col(int j) const;
  virtual const Matrix& get_output_matrix() const;

protected:
  const Matrix& exprData;                   // expressions of the corresponding sequences across multiple conditions
  const Matrix& factorExprData;             // [TF] of all factors over multiple conditions
};

#endif
