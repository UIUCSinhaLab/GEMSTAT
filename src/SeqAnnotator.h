#ifndef SEQ_ANNOTATOR_H
#define SEQ_ANNOTATOR_H

#include <cctype>
#include <cstring>

#include "Tools.h"

using namespace std;

/*****************************************************
 * DNA Sequences
 ******************************************************/

// alphabet
const char ALPHABET[] = { 'A', 'C', 'G', 'T', 'N', '-' };
const int NBASES = 4;                             // number of bases
const int ALPHABET_SIZE = 6;
const int MISSING = 4;                            // missing nt.
const int GAP = 5;                                // '-'

/* Utility functions */
// test if a is a nucleotide {A,G,C,T}
bool isNt( int a );
// the complement of a nt.
int complement( int a );
// conversion of a symbol (nt. or gap or N) to its index in ALPHABET
int symbolToInt( char c );
// character of the strand (+ or -)
char strand2char( bool strand );
// strand (+ or -) from the character
bool char2strand( char c );
// nucleotide distribution from GC content
vector< double > createNtDistr( double gcContent );

// file formats for Sequence
enum FileFormats { FASTA, PHYLIP };

/* Sequence class */
class Sequence
{
    public:
        // constructors
        Sequence() : nts() {}
        Sequence( const vector< int >& _nts ) : nts( _nts ) {}
        Sequence( const string& str );
        Sequence( const Sequence& other, int start, int length, bool strand = true );
        void copy( const Sequence& other ) { nts = other.nts; }
        Sequence( const Sequence& other ) { copy( other ); }

        // assignment
        Sequence& operator=( const Sequence& other ) { copy( other ); return *this; }

        // sequence comparison
        bool operator==( const Sequence& other )
        {
            if ( nts == other.nts ) return true;
            else return false;
        }

        // access methods
        int size() const { return nts.size(); }
        const int& operator[]( int pos ) const
        {
            assert( pos >= 0 && pos < size() );
            return nts[ pos ];
        }
        int& operator[]( int pos )
        {
            assert( pos >= 0 && pos < size() );
            return nts[ pos ];
        }

        // insert a nt. or a sequence segment at the end of sequence
        int push_back( int nt );
        int push_back( const Sequence& elem );

        // the reverse complement of a sequence
        Sequence compRevCompl() const;

        // information of the sequence
        void getNtCounts( vector< int >& counts ) const;
        bool containsMissing() const;

        // clear the sequence
        void clear() { nts.clear(); }

        // load Sequence from a file
        int load( const string& file, string& name, int format = FASTA );
        int load( const string& file, int format = FASTA );

        // output
        friend ostream& operator<<( ostream& os, const Sequence& seq );
    private:
        vector< int > nts;                        // nts[ i ]: the i-th nt.
};

// read sequences from a file
int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format = FASTA );
int readSequences( const string& file, vector< Sequence >& seqs, int format = FASTA );

// write sequences to a file
int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format = FASTA );
int writeSequences( const string& file, const vector< Sequence >& seqs, int format = FASTA );

/*****************************************************
 * DNA Sequence Motifs and Transcription Factors
 ******************************************************/

const double PSEUDO_COUNT = 0.25;

// construct the position weight matrix from the count matrix
Matrix compWtmx( const Matrix& countMatrix, double pseudoCount );

/* Motif class: transcription factor binding site motif */
class Motif
{
    public:
        // constructors
        Motif() : pwm(), background() {}
        Motif( const Matrix& _pwm, const vector< double >& _background );
                                                  // countMatrix in Transfac format
        Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background );
        void copy( const Motif& other ) { pwm = other.pwm; background = other.background; LLRMat = other.LLRMat; maxSite = other.maxSite; maxLLR = other.maxLLR; }
        Motif( const Motif& other ) { copy( other ); }

        // assignment
        Motif& operator=( const Motif& other ) { copy( other ); return *this; }

        // access methods
        int length() const { return pwm.nRows(); }
        const Matrix& getPwm() const { return pwm; }
        const vector< double >& getBackground() const { return background; }
        const Matrix& getLLRMat() const { return LLRMat; }
        const Sequence& getMaxSite() const { return maxSite; }
        double getMaxLLR() const { return maxLLR; }

        // compute the log-likelihood ratio of a sequence element
        double LLR( const Sequence& elem ) const;

        // compute the energy of a sequence element, relative to the strongest site (thus always >= 0)
        double energy( const Sequence& elem ) const;

        // sample a site from PWM
        void sample( const gsl_rng* rng, Sequence& elem, bool strand = true ) const;

        // load a Motif
        int load( const string& file, const vector< double >& background, string& name );
        int load( const string& file, const vector< double >& background );

        // output
        friend ostream& operator<<( ostream& os, const Motif& motif );
    private:
        Matrix pwm;                               // the position weight matrix: f_i(b), the frequency of b at position i
        vector< double > background;              // background distribution
        Matrix LLRMat;                            // LLR matrix, M(i,b) = log( f_i(b) / p(b) ), where f_i(b) is the frequency of b at position i of PWM, and p(b) the frequency of b at the background
        Sequence maxSite;                         // the sequence of the strongest site
        double maxLLR;                            // LLR of the strongest site

        // initialization: compute the LLR matrix and maxLLR
        void init();
};

// read the motifs (PWMs) in a FASTA-like format: pseudocounts within the file
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names );
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs );

/*****************************************************
 * Annotation of Sequences
 ******************************************************/

/* Site class: a TFBS in a sequence with its relevant information */
class Site
{
    public:
        // constructors
        Site() : start( 0 ), strand( true ), factorIdx( 0 ), energy( 0 ), wtRatio( 1 ), prior_probability( 1 ) {}
        Site( int _start, bool _strand, int _factorIdx, double _prior_probability ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( 0 ), wtRatio( 1 ), prior_probability( _prior_probability ) {}
        Site( int _start, bool _strand, int _factorIdx, double _energy, double _prior_probability ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( _energy ), prior_probability( _prior_probability ) { wtRatio = exp( -energy ); }
        void copy( const Site& other ) { start = other.start; strand = other.strand; factorIdx = other.factorIdx; energy = other.energy; wtRatio = other.wtRatio; prior_probability = other.prior_probability; }
        Site( const Site& other ) { copy( other ); }

        // assignment
        Site& operator=( const Site& other ) { copy( other ); return *this; }

        friend ostream& operator<<( ostream& os, const Site& site );

        int start;                                // start position: 0-based
        bool strand;                              // 1: positive; 0: negative
        int factorIdx;                            // the index of the associated TF, starting from 0
        double energy;                            // the energy relative to the strongest site (nonnegative)
        double wtRatio;                           // the binding weight ratio (<= 1) of site vs the strongest site: K(S) / K(S_max) = q(S) / q(S_max)
        double prior_probability;
};

// test if two sites overlap
bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs );

// representation of Sequence as a Site vector
typedef vector< Site > SiteVec;

// read sites (of potentially multiple sequences)
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy = false );
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy = false );

/* SeqAnnotator class: annotate a given sequence by extracting its TFBSs */
class SeqAnnotator
{
    public:
        // constructors
        SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrFactors ) : motifs( _motifs ), energyThrFactors( _energyThrFactors ) { assert( motifs.size() == energyThrFactors.size() ); }

        // annotate a sequence
        int annot( const Sequence& seq, SiteVec& sites ) const;
        //    int annot( const Sequence& seq, SiteVec& sites, const vector < double >& dnase_start, const vector < double >& dnase_end ) const;
        int annot( const Sequence& seq, SiteVec& sites, const vector < double >& dnase_start, const vector < double >& dnase_end, const vector < double >& scores, const double seq_start ) const;
        // compute the energy of sites (update the input sites)
        int compEnergy( const Sequence& seq, SiteVec& sites ) const;
        double sigmoidal( const double score ) const;
        static double alpha;
        static double beta;
    private:
        vector< Motif > motifs;                   // all the TF binding motifs
        vector< double > energyThrFactors;        // energy thresholds for all the motifs
};
#endif
