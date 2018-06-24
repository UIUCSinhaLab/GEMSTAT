/*****************************************************
 * Train and test the expression model
 * Input: sequence, expression, motif, factor expr, cooperativity rules,
 *   activator list, repression rules
 * Output: trained model and expression of training data
 * File formats:
 * (1) Sequence: Fasta format
 * (2) Expression: one line per sequence
 *       <seq_name expr_1 ... expr_m>
 * (3) Motif: Stubb format
 * (4) Factor expression: one line per factor
 *       <factor expr_1 ... expr_m>
 * (5) Cooperativity: list of cooperative factor pairs
 * (6) Factor roles: the role field is either 1 (yes) or 0 (no)
 *       <factor activator_role repressor_role>
 * (7) Repression: list of <R A> where R represses A
 * (8) Parameters: the format is:
 *       <factor binding activation repression>
 *       <factor1 factor2 coop>
 *       <basal_transcription = x>
 *     where repression is optional, and the coop. lines are optional.
 * Note that (5), (6), (7) and (8) may be empty
 ******************************************************/

#include <stdexcept>

#include "utils/gs_errors.h"


#include "SeqAnnotator.h"




double SeqAnnotator::alpha = 6.008;
double SeqAnnotator::beta = 0.207;

/**
Copied from IO.h so we don't have to include that and a bunch of other stuff.
*/
int readFactorThresholdFile( const string& filename, vector< double >& destination, int nFactors){

	ifstream factor_thr_input( filename.c_str() );

	if(!factor_thr_input.is_open() ){
		cerr << "Cannot open the factor threshold file " << filename << endl;
		return RET_ERROR;
	}

	destination.clear();

	for( int index = 0; index < nFactors; index++ )
        {
            double temp;
            factor_thr_input >> temp;
            destination.push_back( temp );
        }
	//Ensure there is only whitespace remaining in the file!
        factor_thr_input.close();

	return 0;
}


void altered_site_print( ostream& os, const Site& site , const vector< Motif >& motifs )
{
	char strandChar = site.strand ? '+' : '-';
    os << site.start + 1;
	if(site.end != -1){
		os << ".." << site.end+1;
	}
	os << "\t" << strandChar << "\t" << motifs[site.factorIdx].get_name() << "\t" << site.energy << "\t" << site.wtRatio;
}

int main( int argc, char* argv[] )
{
    // command line processing
    string seqFile, r_seqFile, annFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, axis_wtFile;
    string outFile;                               // output file
    string dnase_file;
    string factor_thr_file;
    double eTF = 0.60;
    bool traditional_format = false;

    for ( int i = 1; i < argc; i++ )
    {
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];
        else if ( !strcmp( "-m", argv[ i ] ) )
            motifFile = argv[ ++i ];
        else if ( !strcmp( "-et", argv[i] ) )
            eTF = atof( argv[ ++i ] );
        else if ( !strcmp( "-df", argv[ i ]))
            dnase_file = argv[ ++i ];
        else if ( !strcmp( "-ft", argv[ i ]))
            factor_thr_file = argv[ ++i ];
	else if ( !strcmp( "-tf", argv[ i ]))
	    traditional_format = true;
    }

    if ( seqFile.empty() || motifFile.empty() )
    {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -m motifFile [ -et default_factor_energy=<" << eTF << "> -ft <factorThresholdsFile -df <dnase_file> -tf (for traditional output format)]" << endl;
        exit( 1 );
    }

    //     bool readSites = false;     // whether read sites (if true) or read sequences

    // additional control parameters
    double gcContent = 0.5;

    // read the sequences
    vector< Sequence > seqs;
    vector< string > seqNames;
    int rval = readSequences( seqFile, seqs, seqNames );
    ASSERT_MESSAGE(rval != RET_ERROR, "Could not read the sequence file.");
    int nSeqs = seqs.size();

    // read the motifs
    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    rval = readMotifs( motifFile, background, motifs, motifNames );
    ASSERT_MESSAGE( rval != RET_ERROR , "Could not read the motifs.");
    int nFactors = motifs.size();

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < motifNames.size(); i++ )
    {
        factorIdxMap[motifNames[i]] = i;
    }

    //initialize the energy threshold factors
    vector < double > energyThrFactors(nFactors, eTF);

    if( ! factor_thr_file.empty() )
    {
	int readFactorRet = readFactorThresholdFile(factor_thr_file, energyThrFactors, nFactors);
	ASSERT_MESSAGE( 0==readFactorRet , "Difficulty opening the factor_thr_input file.");
    }

    // site representation of the sequences
    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrFactors );
    if ( annFile.empty() )                        // construct site representation
    {

        if( dnase_file.empty() )
        {
            for ( int i = 0; i < nSeqs; i++ )
            {
                //cout << "Annotated sites for CRM: " << seqNames[i] << endl;
                ann.annot( seqs[ i ], seqSites[ i ] );
                seqLengths[i] = seqs[i].size();
            }
        }
        else
        {
            for ( int i = 0; i < nSeqs; i++ )
            {
				throw std::runtime_error("DNAse aware GEMSTAT is not currently implmemented. Sorry.");
				//TODO: Read the dnaase file. The code below is a start. But this feature is not used.
				//See DNAse.h and DNAse.cpp
				/*
				map< string, vector< DNAse_region > > dnase_data = read_DNAse_file(dnase_file);
				DNAseAwareSeqAnnotator dnann(ann);
	            for ( int i = 0; i < nSeqs; i++ )
	            {
					string this_seq_name;//TODO: set this
	                dnann.annot( seqs[ i ], seqSites[ i ], dnase_data[this_seq_name], temp_start );
	                seqLengths[i] = seqs[i].size();
	            }
				*/
            }
        }
    }                                             // read the site representation and compute the energy of sites
    else
    {
        rval = ann.readSites( annFile, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ )
        {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }

    // CHECK POINT
    //     cout << "Sequences:" << endl;
    //     for ( int i = 0; i < seqs.size(); i++ ) cout << seqNames[i] << endl << seqs[i] << endl;
    //     cout << "Expression: " << endl << exprData << endl;
    //     cout << "Factor motifs:" << endl;
    //     for ( int i = 0; i < motifs.size(); i++ ) cout << motifNames[i] << endl << motifs[i] << endl;
    //     cout << "Factor expression:" << endl << factorExprData << endl;
    //     cout << "Cooperativity matrix:" << endl << coopMat << endl;
    //     cout << "Activators:" << endl << actIndicators << endl;
    //     cout << "Repressors:" << endl << repIndicators << endl;
    //     cout << "Repression matrix:" << endl << repressionMat << endl;
         cerr << "Site representation of sequences:" << endl;
         for ( int i = 0; i < nSeqs; i++ ) {
             cout << ">" << seqNames[i] << endl;
             for ( int j = 0; j < seqSites[i].size(); j++ ){
		    if( traditional_format){
			    cout << seqSites[i][j] << endl;
		    }else{
		    	altered_site_print(cout, seqSites[i][j], motifs);
		    	cout << endl;
		    }
		}
         }

    return 0;
}
