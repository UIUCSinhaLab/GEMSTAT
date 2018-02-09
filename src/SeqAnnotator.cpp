#include <sstream>
#include <stdexcept>

#include "SeqAnnotator.h"

#include "utils/gs_parsing.h"

bool isNt( int a )
{
    if ( a < 0 || a > 3 ) return false;
    else return true;
}


int complement( int a )
{
    assert( a >= 0 && a < ALPHABET_SIZE );

    if ( a == 0 ) return 3;
    if ( a == 1 ) return 2;
    if ( a == 2 ) return 1;
    if ( a == 3 ) return 0;
    if ( a == MISSING ) return MISSING;
    if ( a == GAP ) return GAP;
}


int symbolToInt( char c )
{
    char upper = toupper( c );
    for ( int i = 0; i < ALPHABET_SIZE; i++ )
    {
        if ( ALPHABET[ i ] == upper ) return i;
    }

    return -1;
}


char strand2char( bool strand )
{
    if ( strand ) return '+';
    else return '-';
}


bool char2strand( char c )
{
    assert( c == '+' || c == '-' );

    if ( c == '+' ) return true;
    else return false;
}


vector< double > createNtDistr( double gcContent )
{
    assert( gcContent >= 0 && gcContent <= 1.0 );

    vector< double > freqs( 4 );
    freqs[0] = ( 1.0 - gcContent ) / 2.0;
    freqs[1] = gcContent / 2.0;
    freqs[2] = freqs[1];
    freqs[3] = freqs[0];

    return freqs;
}


Sequence::Sequence( const string& str )
{
    for ( int i = 0; i < str.size(); i++ )
    {
        int nt = symbolToInt( str[ i ] );         // could be a NNN or gap
        if ( nt >= 0 && nt < ALPHABET_SIZE )
        {
            nts.push_back( nt );
        }
        else
        {
            cerr << "Illegal symbol: " << nt << " in " << str << endl;
            exit( 1 );
        }
    }
}


Sequence::Sequence( const Sequence& other, int start, int length, bool strand )
{
    assert( start >= 0 && length >= 0 && ( start + length ) <= other.size() );

    for ( int i = 0; i < length; i++ )
    {
        if ( strand ) { nts.push_back( other[ start + i ] ); }
        else { nts.push_back( complement( other[ start + length - 1 - i ] ) ); }
    }
}


int Sequence::push_back( int nt )
{
    assert( nt >= 0 && nt < ALPHABET_SIZE );
    nts.push_back( nt );

    return 0;
}


int Sequence::push_back( const Sequence& elem )
{
    for ( int i = 0; i < elem.size(); i++ ) push_back( elem[ i ] );
    return 0;
}


Sequence Sequence::compRevCompl() const
{
    return Sequence( *this, 0, size(), false );
}

const string Sequence::getName() const
{
    return name;
}

void Sequence::setName(const string iname)
{
    name = iname;
}

void Sequence::getNtCounts( vector< int >& counts ) const
{
    counts.clear();
    for ( int i = 0; i < NBASES; i++ )
    {
        counts.push_back( 0 );
    }

    for ( int i = 0; i < nts.size(); i++ )
    {
        if ( nts[ i ] != GAP ) counts[ nts[ i ] ]++;
    }
}


bool Sequence::containsMissing() const
{
    for ( int i = 0; i < nts.size(); i++ )
    {
        if ( nts[ i ] == MISSING ) return true;
    }

    return false;
}


int Sequence::load( const string& file, string& name, int format )
{
    vector< Sequence > seqs;
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );
    if ( rval == RET_ERROR ) return RET_ERROR; //TODO: throw an exception

    copy( seqs[ 0 ] );
    name = names[ 0 ];
    return rval;
}


int Sequence::load( const string& file, int format )
{
    string name;
    int rval = load( file, name, format );

    return rval;
}


ostream& operator<<( ostream& os, const Sequence& seq )
{
    // output the nts
    for ( int i = 0; i < seq.size(); i++ )
    {
        os << ALPHABET[ seq[ i ] ];
    }

    return os;
}


int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format )
{
    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }//TODO: throw an exception
    seqs.clear();
    names.clear();

    // 	open the file
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; return RET_ERROR; }//TODO: throw an exception

    string line;
    Sequence seq;

    // read sequences: FASTA format
    if ( format == FASTA )
    {
        while ( getline( fin, line ) )
        {
            // add the sequence and start a new sequence if the line starts with >
            //cout << line << endl;
            if ( line[ 0 ] == '>' )
            {
                if ( seq.size() )
                {
                    seqs.push_back( seq );
                    seq.clear();
                }

                stringstream ss( line.substr( 1 ) );
                string name;
                ss >> name;
                names.push_back( name );
                seq.setName( string(name) );
            }
            else
            {
                // check if the line contains content
                int start = line.find_first_not_of( " \t\r" );
                int last = line.find_last_not_of( " \t\r" );
                if ( start == string::npos || last == string::npos ) continue;

                // append the sequence
                for ( int i = start; i <= last; i++ )
                {
                                                  // could be a NNN or gap
                    int nt = symbolToInt( line[ i ] );
                    if ( nt >= 0 && nt < ALPHABET_SIZE )
                    {
                        seq.push_back( nt );
                    }
                    else
                    {
                        cerr << "Illegal symbol: " << nt << " in " << file << endl;
                        return RET_ERROR; //TODO: throw an exception
                    }
                }
            }
        }

        // add the last sequence
        if( seq.size() ) seqs.push_back( seq );

        return 0;
    }
}


int readSequences( const string& file, vector< Sequence >& seqs, int format )
{
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );
    return rval;
}


int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format )
{
    assert( seqs.size() == names.size() );

    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }//TODO: throw an exception

    ofstream fout( file.c_str() );

    if ( format == FASTA )
    {
        for ( int i = 0; i < seqs.size(); i++ )
        {
            fout << ">" << names[ i ] << endl;
            fout << seqs[ i ] << endl;
        }
    }

    return 0;
}


int writeSequences( const string& file, const vector< Sequence >& seqs, int format )
{
    // default name: integer starting from 1
    vector< string > names;
    for ( int i = 0; i < seqs.size(); i++ )
    {
        char buffer[ 10 ];
        sprintf( buffer, "%i", i );
        names.push_back( string( buffer ) );
    }

    // print
    return writeSequences( file, seqs, names, format );
}


Matrix compWtmx( const Matrix& countMatrix, double pseudoCount )
{
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );

    int l = countMatrix.nRows();                  // l: the length of motif
    Matrix pwm( l, 4 );

    //     // the sum of each position/column should be a const. (number of sequences)
    //     double n = 0;		// number of sites used in the count matrix
    //     for ( int j = 0; j < 4; j++ ) {
    //         n += countMatrix( 0, j );
    //     }
    //     for ( int i = 1; i < l; i++ ) {
    //         double count = 0;
    //         for ( int j = 0; j < 4; j++ ) {
    //             count += countMatrix( i, j );
    //         }
    //         if ( count != n ) { cout << "count matrix incorrect" << endl; exit( 1 ); }
    //     }

    // the multinomial distribution at each column
    for ( int i = 0; i < l; i++ )
    {
        double n = 0;                             // total counts at this position
        for ( int j = 0; j < 4; j++ )
        {
            n += countMatrix( i, j );
        }
        for ( int j = 0; j < 4; j++ )
        {
            pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }
    }

    return pwm;
}


Motif::Motif( const Matrix& _pwm, const vector< double >& _background ) : pwm( _pwm ), background( _background ), LLRMat( pwm.nRows(), 4 )
{
    assert( background.size() == 4 );

    init();
}


Motif::Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background ) : background( _background ), LLRMat( countMatrix.nRows(), 4 )
{
    assert( background.size() == 4 );

    pwm = compWtmx( countMatrix, pseudoCount );
    init();
}


double Motif::LLR( const Sequence& elem ) const
{
    int l = pwm.nRows();
    if ( elem.size() != l ) return GSL_NEGINF;
    if ( elem.containsMissing() ) return GSL_NEGINF;

    double result = 0;
    for ( int i = 0; i < l; i++ )
    {
        result += LLRMat( i, elem[ i ] );
    }

    return result;
}


double Motif::energy( const Sequence& elem ) const
{
    return ( -LLR( elem ) + maxLLR );
}


void Motif::sample( const gsl_rng* rng, Sequence& elem, bool strand ) const
{
    assert( rng != NULL );

    int l = pwm.nRows();
    Sequence sampleElem;
    for ( int i = 0; i < l; i++ )
    {
        // nt. distribution at position i
        vector< double > distr = pwm.getRow( i );

        // sample nt. from this distribution
        int nt = sampleMul( rng, distr );
        sampleElem.push_back( nt );
    }

    if ( strand == 0 ) elem = sampleElem.compRevCompl();
    else elem = sampleElem;
}


int Motif::load( const string& file, const vector< double >& background, string& name )
{
    vector< Motif > motifs;
    vector< string > names;
    int rval = readMotifs( file, background, motifs, names );
    if ( rval == RET_ERROR ) return RET_ERROR;//TODO: throw an exception

    copy( motifs[ 0 ] );
    name = names[ 0 ];
    return rval;
}


int Motif::load( const string& file, const vector< double >& background )
{
    string name;
    int rval = load( file, background, name );

    return rval;
}


ostream& operator<<( ostream& os, const Motif& motif )
{
    os << motif.pwm;

    return os;
}


void Motif::init()
{
    int l = pwm.nRows();

    // compute the LLR matrix
    for ( int i = 0; i < l; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            LLRMat( i, j ) = log( pwm( i, j ) / background[ j ] );
        }
    }

    // the strongest site
    for ( int i = 0; i < l; i++ )
    {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite.push_back( b_max );
    }

    // compute the LLR of the strongest site
    maxLLR = 0;
    for ( int i = 0; i < l; i++ )
    {
        maxLLR += LLRMat( i, maxSite[ i ] );
    }
}


int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names )
{
    // 	open the file
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; return RET_ERROR; } //TODO: throw an exception
    motifs.clear();
    names.clear();

    string line;
    SimpleStringTokenizer my_tok;

    // read the motifs
    do
    {
        getline( fin, line );
        int comment_starts_at = line.find("#");
        if( 0 == comment_starts_at ) continue;//whole line is a comment line
        if( 1 <= comment_starts_at){
            line = line.substr(0,comment_starts_at);//Truncate line to remove comment
        }

        if(fin.eof()){//TODO: maybe this is not good?
            break;
        }

        if ( line[ 0 ] != '>' ) throw std::runtime_error("Bad format for PWM (Expected start of PWM)");

        // read the names, length and pseudocount

        my_tok.tokenize(line);
        string name = my_tok[0].substr(1);//Remember to throw away the ">".
        int length = atoi( my_tok[1].c_str() );
        double pseudoCount = my_tok.size() >= 3 ? atof(my_tok[2].c_str()) : PSEUDO_COUNT;

        // read the count matrix
        Matrix countMat( length, NBASES );
        for ( int i = 0; i < length; ++i )
        {
            getline( fin, line );//No comment handling...
            my_tok.tokenize(line);

            if(my_tok.size() != NBASES ) throw std::runtime_error("Encountered a line with the wrong number of bases while reading PWM.");

            //read one line out of the count matrix.
            for ( int j = 0; j < NBASES; ++j )
            {
                countMat( i, j ) = atof(my_tok[j].c_str());
            }
        }

        //END OF MOTIF
        getline( fin, line );
        if (0 != line.compare("<")) throw std::runtime_error("Bad format for PWM (unclosed PWM)");

        // create the motif
        names.push_back( string( name ) );
        motifs.push_back( Motif( countMat, pseudoCount, background ) );
    } while ( !fin.eof() );

    return 0;
}


int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs )
{
    vector< string > names;
    return readMotifs( file, background, motifs, names );
}


ostream& operator<<( ostream& os, const Site& site )
{
    char strandChar = site.strand ? '+' : '-';
    os << site.start + 1 << "\t" << strandChar << "\t" << site.factorIdx << "\t" << site.energy << "\t" << site.wtRatio;

    return os;
}


bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;

    return true;
}


int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy )
{
    ifstream fin( file.c_str() );
    if ( !fin )
    {
        return RET_ERROR; //TODO: throw an exception
    }
    sites.clear();
    names.clear();

    SiteVec currVec;
    int nrecords = 0;                             // number of ">" read so far
    while ( !fin.eof() )
    {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;

        if ( line.substr( 0, 1 ) == ">" )
        {
            stringstream ss( line.substr( 1 ) );
            string name;
            ss >> name;
            names.push_back( name );
            nrecords++;
            if ( nrecords > 1 )
            {
                sites.push_back( currVec );
                currVec.clear();
            }
        }
        else
        {
            int start;
            char strandChar;
            string factor;
            double energy = 0;
            stringstream ss( line );
            ss >> start >> strandChar >> factor;
            if ( readEnergy ) ss >> energy;
            bool strand = strandChar == '+' ? 1 : 0;
            map<string, int>::const_iterator iter = factorIdxMap.find( factor );
            if(iter == factorIdxMap.end()){
              cerr << "The site annotation file reffered to a factor that doesn't exist. \n(Did you use one with factor numbers instead of names? The third column must be textual names.)" << endl;
              exit(1);
            } //TODO: Throw an exception if the factor couldn't be found!
            currVec.push_back( Site( start - 1, strand, iter->second , energy, 1 ) );
        }
    }

    sites.push_back( currVec );

    return 0;
}


int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy )
{
    vector< string > names;
    return readSites( file, factorIdxMap, sites, names, readEnergy );
}


int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites, const vector < double >& dnase_start, const vector < double >& dnase_end, const vector < double >& scores, const double seq_start ) const
{
    //cout << "start annotation:" << endl;
    sites.clear();

    // scan the sequence for the sites of all motifs
    int dnase_data_size = dnase_start.size();
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
                while( win_end > dnase_end[ win_1 ] )
                {
                    win_1++;
                }
                int win_2;
                if( dnase_start[ win_1 ] > win_start )
                {
                    win_2 = win_1 - 1;
                }
                else
                {
                    win_2 = win_1;
                }
                //cout << "window 1: " << win_1 << ", window 2: " << win_2 << endl;
                //cout << "window 1 score: " << scores[ win_1 ] << ", window 2 score: " << scores[ win_2 ] << endl;

                double score = ( scores[ win_1 ] + scores[ win_2 ] ) / 2;
                //cout << "Score: " << score << endl;
                double prior_prob = sigmoidal( score );
                //cout << "prior_prob: " << prior_prob << endl;
                sites.push_back( Site( i, 1, k, energy, prior_prob ) );
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
                while( win_end > dnase_end[ win_1 ] )
                {
                    win_1++;
                }
                int win_2;
                if( dnase_start[ win_1 ] > win_start )
                {
                    win_2 = win_1 - 1;
                }
                else
                {
                    win_2 = win_1;
                }
                //cout << "window 1: " << win_1 << ", window 2: " << win_2 << endl;
                //cout << "window 1 score: " << scores[ win_1 ] << ", window 2 score: " << scores[ win_2 ] << endl;
                double score = ( scores[ win_1 ] + scores[ win_2 ] ) / 2;
                //cout << "Score: " << score << endl;
                double prior_prob = sigmoidal( score );
                //cout << "prior_prob: " << prior_prob << endl;
                sites.push_back( Site( i, 0, k, energy, prior_prob ) );
            }
        }
    }

    //cout << "end annotation" << endl;
    return sites.size();
}


int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites ) const
{
    //cout << "start annotation:" << endl;
    sites.clear();

    // scan the sequence for the sites of all motifs
    for ( int i = 0; i < seq.size(); i++ )
    {
        // test for each motif
        for ( int k = 0; k < motifs.size(); k++ )
        {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) continue;
            double energy;

            // positive strand
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrFactors[ k ] * motifs[ k ].getMaxLLR() )
            {
                sites.push_back( Site( i, 1, k, energy, 1 ) );
            }

            // negative strand
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrFactors[ k ]  * motifs[k].getMaxLLR() )
            {
                sites.push_back( Site( i, 0, k, energy, 1 ) );
            }
        }
    }

    //cout << "end annotation" << endl;
    return sites.size();
}


int SeqAnnotator::compEnergy( const Sequence& seq, SiteVec& sites ) const
{
    for ( int i = 0; i < sites.size(); i++ )
    {
        Sequence elem( seq, sites[i].start, motifs[sites[i].factorIdx].length(), sites[i].strand );
        sites[i].energy = motifs[sites[i].factorIdx].energy( elem );
        sites[i].wtRatio = exp( -sites[i].energy );
    }

    return 0;
}


double SeqAnnotator::sigmoidal( const double score ) const
{
    //cout << "In sigmoidal: " << endl;
    //cout << "score: " << score << endl;
    //cout << "beta: " << beta << endl;
    //cout << "alpha: " << alpha << endl;
    //cout << "exp(-beta * score + alpha): " << exp( -beta * score + alpha ) << endl;
    //cout << "prior_prob: " << 1 / ( 1 + exp( -beta * score + alpha ) ) << endl;
    double prior_prob = 1 / ( 1 + exp( -beta * score + alpha ) );
    //cout << "Out sigmoidal: " << endl;
    return prior_prob;
}
