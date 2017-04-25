#ifndef EXPR_FUNC_H
#define EXPR_FUNC_H

#include "ExprModel.h"
#include "FactorIntFunc.h"
#include "ExprPar.h"


/*****************************************************
 * Expression Model and Parameters
 ******************************************************/


/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class ExprFunc
{
    public:
        // constructors
        ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par );

        // access methods
        const vector< Motif >& getMotifs() const
        {
            return motifs;
        }

        // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
        double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num );
        const ExprPar& getPar() const { return par; }

        static ModelType modelOption;             // model option
        static bool one_qbtm_per_crm;
    private:
        // TF binding motifs
        const vector< Motif >& motifs;

        // control parameters
        const FactorIntFunc* intFunc;             // function to compute distance-dependent TF-TF interactions
        const vector< bool >& actIndicators;      // 1 if the TF is in the activator set
        int maxContact;                           // the maximum contact
        const vector< bool >& repIndicators;      // 1 if the TF is in the repressor set
        const IntMatrix& repressionMat;           // repression matrix: R(f,f') = 1 if f can repress f'
        double repressionDistThr;                 // distance threshold for repression: d_R

        // model parameters
        ExprPar par;//NOTE: Removing "const" here caused the copy constructor to be called. Thus the ExprFunc gets its own copy that will not have problems when the original par is changed.

        // the sequence whose expression is to be predicted
        SiteVec sites;
        vector< int > boundaries;                 // left boundary of each site beyond which there is no interaction

        // intermediate computational results
        vector< double > bindingWts;

        // compute the partition function when the basal transcriptional machinery (BTM) is not bound
        double compPartFuncOff() const;

        // compute the partition function when the BTM is not bound: ChrMod model
        double compPartFuncOffChrMod() const;

        // compute the partition function when the BTM is bound
        double compPartFuncOn() const;

        // compute the paritition function when the BTM is bound: Direct model
        double compPartFuncOnDirect() const;

        // compute the paritition function when the BTM is bound: Quenching model
        double compPartFuncOnQuenching() const;

        // compute the paritition function when the BTM is bound: ChrMod_Unlimited model
        double compPartFuncOnChrMod_Unlimited() const;

        // compute the paritition function when the BTM is bound: ChrMod_Limited model
        double compPartFuncOnChrMod_Limited() const;

        // compute the TF-TF interaction between two occupied sites
        double compFactorInt( const Site& a, const Site& b ) const;

        double compDen() const;
        double compNum( int TFid ) const;

        // test if one site represses another site
        bool testRepression( const Site& a, const Site& b ) const;

        //setup functions that may be useful to subclasses
        inline void setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num);
        inline void setupBindingWeights(const vector< double >& factorConcs);
};

#endif
