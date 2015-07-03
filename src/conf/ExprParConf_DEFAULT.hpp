#ifndef EXPRPARCONF_H
#define EXPRPARCONF_H

ModelType ExprPar::modelOption = CHRMOD_UNLIMITED;
SearchType ExprPar::searchOption = UNCONSTRAINED;
int ExprPar::estBindingOption = 1;                // 1. estimate binding parameters; 0. not estimate binding parameters

double ExprPar::default_weight = 1.0;
double ExprPar::default_interaction = 1.0;
double ExprPar::default_effect_Logistic = 0.0;
double ExprPar::default_effect_Thermo = 1.0;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -5.0;
double ExprPar::default_basal_Thermo = 0.01;
double ExprPar::default_pi = 1E-50;
double ExprPar::min_weight =0.01;
double ExprPar::max_weight = 250.5;
double ExprPar::min_interaction = 0.01;
double ExprPar::max_interaction = 100.5;
double ExprPar::min_effect_Logistic = -5;
double ExprPar::max_effect_Logistic = 5;
// double ExprPar::min_effect_Direct = 0.01;
double ExprPar::min_effect_Thermo = 1;
double ExprPar::max_effect_Thermo = 4.55;
double ExprPar::min_repression = 1E-40;
double ExprPar::max_repression = 1;
double ExprPar::min_basal_Logistic = -9.0;
double ExprPar::max_basal_Logistic = -1.0;
double ExprPar::min_basal_Thermo = 1.0E-5;
double ExprPar::max_basal_Thermo = 0.105;
double ExprPar::delta = -1.0E-75;       //We accept parameter vectors up to this far (per parameter) INSIDE the constrained space. Negative values mean slack outside the constrained interval. Positive values expose another bug that cannot be fixed immediately.
double ExprPar::default_beta = 5;
double ExprPar::min_beta = 1.0E-4;
double ExprPar::max_beta = 500;

double ExprPar::min_energyThrFactors = 0.1;
double ExprPar::max_energyThrFactors = 0.99;
double ExprPar::default_energyThrFactors = 0.9;

double ExprPar::min_pi = 1.0E-75;
double ExprPar::max_pi = 1E10;

bool ExprPar::one_qbtm_per_crm = false;
bool ExprFunc::one_qbtm_per_crm = false;

double SeqAnnotator::alpha = 6.008;
double SeqAnnotator::beta = 0.207;

#endif
