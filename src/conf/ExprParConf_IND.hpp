#ifndef EXPRPARCONF_H
#define EXPRPARCONF_H

ModelType ExprPar::modelOption = CHRMOD_UNLIMITED;
SearchType ExprPar::searchOption = UNCONSTRAINED;

double ExprPar::default_weight = 1.0;
double ExprPar::default_interaction = 1.0;
double ExprPar::default_effect_Logistic = 0.0;
double ExprPar::default_effect_Thermo = 1.0;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -5.0;
double ExprPar::default_basal_Thermo = 0.001;
double ExprPar::default_pi = 1E-50;
double ExprPar::min_weight = 0.0099;
double ExprPar::max_weight = 10000.0001;
double ExprPar::min_interaction = 0.99;
double ExprPar::max_interaction = 100.0001;
double ExprPar::min_effect_Logistic = -5;
double ExprPar::max_effect_Logistic = 5;
double ExprPar::min_effect_Thermo = 0.99;
double ExprPar::max_effect_Thermo = 10.0001;
double ExprPar::min_repression = 9.9E-6;
double ExprPar::max_repression = 0.990001;
double ExprPar::min_basal_Logistic = -9.0;
double ExprPar::max_basal_Logistic = -1.0;
double ExprPar::min_basal_Thermo = 9E-4;
double ExprPar::max_basal_Thermo = 0.010001;
double ExprPar::delta = -1.0E-15;
double ExprPar::default_beta = 5;
double ExprPar::min_beta = 1.0E-4;
double ExprPar::max_beta = 500;

double ExprPar::min_energyThrFactors = 0.005;
double ExprPar::max_energyThrFactors = 0.99;
double ExprPar::default_energyThrFactors = 0.9;

double ExprPar::min_pi = 1.0E-75;
double ExprPar::max_pi = 1E10;

bool ExprPar::one_qbtm_per_crm = false;
bool ExprFunc::one_qbtm_per_crm = false;

double SeqAnnotator::alpha = 6.008;
double SeqAnnotator::beta = 0.207;

#endif
