
#include "DataSet_signal.h"

DataSet_Signal::DataSet_Signal(const Matrix& tf_concentrations, const Matrix& output_values, const Matrix& signaling_values) : TrainingDataset(tf_concentrations,output_values),
signaling_data(signaling_values), row_names(), row_names_to_row()
{
  assert(exprData.nCols() == factorExprData.nCols());
}

int DataSet_Signal::nConds() const{
  return factorExprData.nCols();
}

void DataSet_Signal::set_row_names(const vector<string> &in_names){
    row_names = in_names;
    for(int i = 0;i<row_names.size();i++){
        row_names_to_row[row_names[i]] = i;
    }
}

void DataSet_Signal::set_signal_row_names(const vector<string> &in_names){
    signal_row_names = in_names;
    for(int i = 0;i<signal_row_names.size();i++){
        signal_row_names_to_row[signal_row_names[i]] = i;
    }
}

Condition DataSet_Signal::getCondition(int i, const ExprPar &signaling_params) const{


    vector< double > values = this->factorExprData.getCol(i);

    //We're just going to put CIC attenuation right here. blah. It should be a subclass.

    int cic_i = this->row_names_to_row.at("cic");

    int dperk_i = this->signal_row_names_to_row.at("dperk");


    //cerr << "TYpE : " << ((gsparams::DictList)signaling_params.my_pars).at("signaling").at("cic_att").my_type << endl;
    double att_param = 0.0;
    try{
        att_param = 0.01;
        att_param = ((gsparams::DictList)(signaling_params.my_pars)).at("signaling").at("cic_att").v();
    }catch(std::out_of_range e){
        //cerr << signaling_params.my_pars << endl;
        throw e;
    }

    double foobar = signaling_data.getElement(dperk_i,i);
    //double foobar = 0.0;
    values[cic_i] = values[cic_i]*pow(1.0 + foobar*exp(att_param),-1.0);

    //pMAD (signalled by dpp)
    if( this->row_names_to_row.find("mad") != this->row_names_to_row.end() && this->signal_row_names_to_row.find("dpp") != this->signal_row_names_to_row.end() ){
    	int mad_i = this->row_names_to_row.at("mad");
	int dpp_i = this->signal_row_names_to_row.at("dpp");
	double dpp_mad_to_pmad = 1.0;
	try{
		 dpp_mad_to_pmad = ((gsparams::DictList)(signaling_params.my_pars)).at("signaling").at("dpp_mad_pmad").v();
	}catch(std::out_of_range e){
		cerr << "could not find 'signaling::dpp_mad_pmad' in parameters" << endl;
	}
	values[mad_i] = values[mad_i]*dpp_mad_to_pmad*signaling_data.getElement(dpp_i,i);
    }

    return Condition(values);
}
