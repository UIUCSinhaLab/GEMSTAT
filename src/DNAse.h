#ifndef DNASE_H
#define DNASE_H

#include <string>
#include <map>
#include <vector>

#include "SeqAnnotator.h"

typedef struct {
	double start;
	double end;
	double score;
} DNAse_region;

map< string, vector< DNAse_region > > read_DNAse_file(const string& filename);

class DNAseAwareSeqAnnotator : public SeqAnnotator {
	public:
		DNAseAwareSeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrFactors ) : SeqAnnotator(_motifs,_energyThrFactors) {}
		DNAseAwareSeqAnnotator( const SeqAnnotator& other) : SeqAnnotator(other.motifs, other.energyThrFactors) {}
		~DNAseAwareSeqAnnotator();
		int annot( const Sequence& seq, SiteVec& sites, const vector< DNAse_region >& regions, const double seq_start ) const;
};

#endif
