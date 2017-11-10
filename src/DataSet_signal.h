/*
 * DataSet.h
 *
 *  Created on: April 26, 2017
 *      Author: lunt
 */

#ifndef SRC_DATASET_SIGNAL_H_
#define SRC_DATASET_SIGNAL_H_

#include "DataSet.h"

#include "Tools.h"
#include "ExprPar.h"

class DataSet_Signal : public DataSet {
public:
  DataSet_Signal(const Matrix& tf_concentrations, const Matrix& output_values, const Matrix& signaling_values);
  ~DataSet_Signal(){};

  void set_row_names(const vector<string> &in_names);
  void set_signal_row_names(const vector<string> &in_names);

  int nConds() const;

  virtual Condition getCondition(int i , const ExprPar &signaling_params) const;

  const Matrix& signaling_data;
  vector<string> row_names;
  map<string, int> row_names_to_row;
  vector<string> signal_row_names;
  map<string, int> signal_row_names_to_row;
};

#endif
