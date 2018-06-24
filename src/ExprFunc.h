#ifndef EXPR_FUNC_H
#define EXPR_FUNC_H

#include "ExprModel.h"
#include "FactorIntFunc.h"
#include "ExprPar.h"
#include "DataSet.h"

/*****************************************************
 * Expression Model and Parameters
 ******************************************************/

typedef long double gemstat_dp_t;

/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class ExprFunc
{
    public:
        // constructors
        ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num);

        // access methods
        const vector< Motif >& getMotifs() const
        {
            return motifs;
        }

        // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
        virtual double predictExpr( const vector< double >& factorConcs );
        virtual double predictExpr( const Condition& in_condition );
        const ExprPar& getPar() const { return par; }

        //static ModelType modelOption;             // model option
        static bool one_qbtm_per_crm;
    protected:
        //setup functions that may be useful to subclasses
        virtual void setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num);
        void setupBindingWeights(const vector< double >& factorConcs);
        // TF binding motifs

        const ExprModel* expr_model;
        const vector< Motif >& motifs;

        // control parameters
//        const FactorIntFunc* intFunc;             // function to compute distance-dependent TF-TF interactions
        const vector< bool >& actIndicators;      // 1 if the TF is in the activator set
        int maxContact;                           // the maximum contact
        const vector< bool >& repIndicators;      // 1 if the TF is in the repressor set
        const IntMatrix& repressionMat;           // repression matrix: R(f,f') = 1 if f can repress f'
        double repressionDistThr;                 // distance threshold for repression: d_R
        int seq_number;
        int seq_length;

        // model parameters
        ExprPar par;//NOTE: Removing "const" here caused the copy constructor to be called. Thus the ExprFunc gets its own copy that will not have problems when the original par is changed.
                          //NOTE: (Additional) put const back, copying is slow, and par should be constant for an ExprFunc, this also fixed a memory leak.

        // the sequence whose expression is to be predicted
        SiteVec sites;
        int n_sites;  //Useful because the sitevec with pseudosites etc. might change, this should be the number of true sites.
        vector< int > boundaries;                 // left boundary of each site beyond which there is no interaction

        // intermediate computational results
        vector< gemstat_dp_t > bindingWts;


        // compute the TF-TF interaction between two occupied sites
        double compFactorInt( const Site& a, const Site& b ) const;

        double compDen() const;
        double compNum( int TFid ) const;

        // test if one site represses another site
        bool testRepression( const Site& a, const Site& b ) const;

        // compute the partition function when the BTM is bound
        virtual gemstat_dp_t compPartFuncOn() const;
        // compute the partition function when the basal transcriptional machinery (BTM) is not bound
        virtual gemstat_dp_t compPartFuncOff() const;


        /*
        At setup time, this will get populated from the SNOT object.
        */
        // parameters
        vector < GEMSTAT_PAR_FLOAT_T > maxBindingWts;          // binding weight of the strongest site for each TF: K(S_max) [TF_max]
        Matrix factorIntMat;                      // (maximum) interactions between pairs of factors: omega(f,f')
        vector < GEMSTAT_PAR_FLOAT_T > txpEffects;             // transcriptional effects: alpha for Direct and Quenching model, exp(alpha) for Logistic model (so that the same default values can be used). Equal to 1 if a TF is not an activator under the Quenching model
        vector < GEMSTAT_PAR_FLOAT_T > repEffects;             // repression effects: beta under ChrMod models (the equlibrium constant of nucleosome association with chromatin). Equal to 0 if a TF is not a repressor.
        //vector < GEMSTAT_PAR_FLOAT_T > pis;
        //     double expRatio; 		// constant factor of measurement to prediction

        vector < GEMSTAT_PAR_FLOAT_T > betas;
        //vector < GEMSTAT_PAR_FLOAT_T > energyThrFactors;

    private:
};

class Logistic_ExprFunc : public ExprFunc {
  public:
      // constructors
      Logistic_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;

      double predictExpr( const vector< double >& factorConcs );
};

class Markov_ExprFunc : public ExprFunc {
  public:
      Markov_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num);// : ExprFunc( _model, _par , sites_, seq_len, seq_num);
      double predictExpr( const vector< double >& factorConcs );
  protected:
      //void setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num);

      virtual double expr_from_config( const vector< double >& marginals);
      vector<int> rev_bounds;
};

class Direct_ExprFunc : public ExprFunc {
  public:
      // constructors
      Direct_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;

};

class Quenching_ExprFunc : public ExprFunc {
  public:
      // constructors
      Quenching_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;

};

class ChrMod_ExprFunc : public ExprFunc {
  public:
      // constructors
      ChrMod_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    virtual gemstat_dp_t compPartFuncOn() const = 0;
    // compute the partition function when the basal transcriptional machinery (BTM) is not bound
    gemstat_dp_t compPartFuncOff() const;

};

class ChrModUnlimited_ExprFunc : public ChrMod_ExprFunc {
  public:
      // constructors
      ChrModUnlimited_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ChrMod_ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;
};

class ChrModLimited_ExprFunc : public ChrMod_ExprFunc {
  public:
      // constructors
      ChrModLimited_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ChrMod_ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;
};


#endif
