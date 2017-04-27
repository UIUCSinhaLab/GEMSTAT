/*
 * DataSet.h
 *
 *  Created on: April 26, 2017
 *      Author: lunt
 */

#ifndef SRC_DATASET_H_
#define SRC_DATASET_H_

#include "Tools.h"

class Condition{
  public:
    Condition(vector< double > _concs);
    vector< double > concs;
};

class DataSet{
public:
  DataSet(const Matrix& tf_concentrations, const Matrix& output_values);
  ~DataSet(){};

  int nConds() const;

  Condition getCondition(int i) const;

  const Matrix& exprData;                   // expressions of the corresponding sequences across multiple conditions
  const Matrix& factorExprData;             // [TF] of all factors over multiple conditions
};

#endif
