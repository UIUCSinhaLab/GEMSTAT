#ifndef EXPR_FUNC_H
#define EXPR_FUNC_H

#include "ExprModel.h"
#include "FactorIntFunc.h"
#include "ExprPar.h"
#include "DataSet.h"

/*****************************************************
 * Expression Model and Parameters
 ******************************************************/


/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class ExprFunc
{
    public:
        // constructors
        ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par );

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
        const vector< Motif >& motifs;

        // control parameters
        const FactorIntFunc* intFunc;             // function to compute distance-dependent TF-TF interactions
        const vector< bool >& actIndicators;      // 1 if the TF is in the activator set
        int maxContact;                           // the maximum contact
        const vector< bool >& repIndicators;      // 1 if the TF is in the repressor set
        const IntMatrix& repressionMat;           // repression matrix: R(f,f') = 1 if f can repress f'
        double repressionDistThr;                 // distance threshold for repression: d_R
        int seq_number;
        int seq_length;

        // model parameters
        const ExprPar par;//NOTE: Removing "const" here caused the copy constructor to be called. Thus the ExprFunc gets its own copy that will not have problems when the original par is changed.
                          //NOTE: (Additional) put const back, copying is slow, and par should be constant for an ExprFunc, this also fixed a memory leak.

        // the sequence whose expression is to be predicted
        SiteVec sites;
        int n_sites;  //Useful because the sitevec with pseudosites etc. might change, this should be the number of true sites.
        vector< int > boundaries;                 // left boundary of each site beyond which there is no interaction

        // intermediate computational results
        vector< double > bindingWts;


        // compute the TF-TF interaction between two occupied sites
        double compFactorInt( const Site& a, const Site& b ) const;

        double compDen() const;
        double compNum( int TFid ) const;

        // test if one site represses another site
        bool testRepression( const Site& a, const Site& b ) const;

        // compute the partition function when the BTM is bound
        virtual double compPartFuncOn() const;
        // compute the partition function when the basal transcriptional machinery (BTM) is not bound
        virtual double compPartFuncOff() const;

    private:
};

class Logistic_ExprFunc : public ExprFunc {
  public:
      // constructors
      Logistic_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;

      double predictExpr( const vector< double >& factorConcs );
};

class Markov_ExprFunc : public ExprFunc {
  public:
      Markov_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par );// : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par);
      double predictExpr( const vector< double >& factorConcs );
  protected:
      //void setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num);

      virtual double expr_from_config( const vector< double >& marginals);
      vector<int> rev_bounds;
};

class Direct_ExprFunc : public ExprFunc {
  public:
      // constructors
      Direct_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;
  protected:
    // compute the partition function when the BTM is bound
    double compPartFuncOn() const;

};

class Quenching_ExprFunc : public ExprFunc {
  public:
      // constructors
      Quenching_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;
  protected:
    // compute the partition function when the BTM is bound
    double compPartFuncOn() const;

};

class ChrMod_ExprFunc : public ExprFunc {
  public:
      // constructors
      ChrMod_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;
  protected:
    // compute the partition function when the BTM is bound
    virtual double compPartFuncOn() const = 0;
    // compute the partition function when the basal transcriptional machinery (BTM) is not bound
    double compPartFuncOff() const;

};

class ChrModUnlimited_ExprFunc : public ChrMod_ExprFunc {
  public:
      // constructors
      ChrModUnlimited_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ChrMod_ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;
  protected:
    // compute the partition function when the BTM is bound
    double compPartFuncOn() const;
};

class ChrModLimited_ExprFunc : public ChrMod_ExprFunc {
  public:
      // constructors
      ChrModLimited_ExprFunc( const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ChrMod_ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;
  protected:
    // compute the partition function when the BTM is bound
    double compPartFuncOn() const;
};


#endif
