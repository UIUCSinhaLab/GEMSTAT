#include <string>
#include <map>
#include <vector>

#include <iostream>
#include <fstream>

#include "DNAse.h"

using namespace std;

map< string, vector< DNAse_region > > read_DNAse_file(const string& filename){

	//cout << "Annotated sites for CRM: " << seqNames[i] << endl;
	ifstream dnase_input( filename.c_str() );
	//TODO: Make it an exception, not an assertion.
	assert( dnase_input.is_open());

	string temp_s;
	string temp_gen;
	string chr;
	double temp_start, temp_end;
	vector < double > dnase_start;
	vector < double > dnase_end;
	vector < double > scores;

	map< string, vector< DNAse_region > > gathered_dnase_data();

	while( dnase_input >> temp_s )
	{
		dnase_input >> temp_gen >> temp_gen >> temp_gen >> temp_gen >> temp_gen;
		dnase_input >> chr;
		dnase_input >> temp_gen >> temp_start >> temp_end >> temp_gen;
		if( temp_s == seqNames[ i ] )
		{
			//cout << "Processing for:\t" << temp_s << endl;
			dnase_start.clear();
			dnase_end.clear();
			scores.clear();
			ifstream chr_input( ("chr" + chr + ".bed").c_str() );
			//cout << "File: " << "chr" + chr + ".bed" << "opened" << endl;
			//TODO: Make an exception, not an assertion
			assert( chr_input.is_open() );
			double chr_start, chr_end, chr_score;
			//cout << "Starting location on chromosome: " << (long long int)temp_start << endl;
			//cout << "Ending location on chromosome: " << (long long int)temp_end << endl;
			while( chr_input >> temp_gen >> chr_start >> chr_end >> temp_gen >> chr_score )
			{
				if( ( chr_start < temp_start && chr_end < temp_start ) || ( chr_start > temp_end && chr_end > temp_end ) )
				{
					;
				}
				else
				{
					dnase_start.push_back( chr_start );
					dnase_end.push_back( chr_end );
					scores.push_back( chr_score );
					//cout << "Inserting: " << (long long int)chr_start << "\t" << (long long int)chr_end << "\t" << chr_score << endl;
				}
			}
			chr_input.close();
			break;
		}
	}

	dnase_input.close();

}

int DNAseAwareSeqAnnotator::annot( const Sequence& seq, SiteVec& sites, const vector< DNAse_region >& regions, const double seq_start ) const
{
    //cout << "start annotation:" << endl;
    sites.clear();

    // scan the sequence for the sites of all motifs
    int dnase_data_size = regions.size();
    for ( int i = 0; i < seq.size(); i++ )
    {
        // test for each motif
        for ( int k = 0; k < motifs.size(); k++ )
        {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) continue;
            double energy;

            //cout << "For motif: " << k << ", having len: " << l << ", at position: " << i << endl;

            // positive strand
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrFactors[ k ] * motifs[ k ].getMaxLLR() )
            {
                double win_start = seq_start + i;
                double win_end = seq_start + i + l - 1;
                //cout << "window start: " << (long long int)win_start << ", window end: " << (long long int)win_end << endl;
                int win_1 = 0;
                //cout << "accessibility data: " << (long long int) dnase_start[win_1] << "\t" << (long long int)dnase_end[win_1] << endl;
                while( win_end > regions[win_1].end )
                {
                    win_1++;
                }
                int win_2;
                if( regions[ win_1 ].start > win_start )
                {
                    win_2 = win_1 - 1;
                }
                else
                {
                    win_2 = win_1;
                }
                //cout << "window 1: " << win_1 << ", window 2: " << win_2 << endl;
                //cout << "window 1 score: " << scores[ win_1 ] << ", window 2 score: " << scores[ win_2 ] << endl;

                double score = ( regions[ win_1 ].score + regions[ win_2 ].score ) / 2;
                //cout << "Score: " << score << endl;
                double prior_prob = sigmoidal( score );
                //cout << "prior_prob: " << prior_prob << endl;
                sites.push_back( Site( i, i+l-1, 1, k, energy, prior_prob ) );
            }

            // negative strand
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrFactors[ k ]  * motifs[k].getMaxLLR() )
            {
                double win_start = seq_start + i;
                double win_end = seq_start + i + l - 1;
                //cout << "window start: " << (long long int)win_start << ", window end: " << (long long int)win_end << endl;
                int win_1 = 0;
                //cout << "accessibility data: " << (long long int) dnase_start[win_1] << "\t" << (long long int)dnase_end[win_1] << endl;
                while( win_end > regions[ win_1 ].end )
                {
                    win_1++;
                }
                int win_2;
                if( regions[ win_1 ].start > win_start )
                {
                    win_2 = win_1 - 1;
                }
                else
                {
                    win_2 = win_1;
                }
                //cout << "window 1: " << win_1 << ", window 2: " << win_2 << endl;
                //cout << "window 1 score: " << scores[ win_1 ] << ", window 2 score: " << scores[ win_2 ] << endl;
                double score = ( regions[ win_1 ].score + regions[ win_2 ].score ) / 2;
                //cout << "Score: " << score << endl;
                double prior_prob = sigmoidal( score );
                //cout << "prior_prob: " << prior_prob << endl;
                sites.push_back( Site( i, i+l-1, 0, k, energy, prior_prob ) );
            }
        }
    }

    //cout << "end annotation" << endl;
    return sites.size();
}
