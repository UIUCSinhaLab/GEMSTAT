/*
* the library of general tools
*/
#ifndef TOOLS_H
#define TOOLS_H
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_int.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

using namespace std;

/* constants */
const int RET_ERROR = -1;		// return error of a function
const int INT_INF = numeric_limits< int >::max() / 2;

/*****************************************************
* Vectors and Matrices 
******************************************************/

/* functions for gsl_vector and gsl_vector_int */
vector< double > gsl2vector( const gsl_vector* v );
gsl_vector* vector2gsl( const vector< double >& v );
vector< int > gsl_int2vector( const gsl_vector_int* v );
gsl_vector_int* vector2gsl_int( const vector< int >& v );

/* Matrix class */
class Matrix {
    gsl_matrix* data;	// data
public:
    // constructors
    Matrix() : data( NULL ) {}
    Matrix( int nRows, int nCols );
    Matrix( int nRows, int nCols, double init );
    Matrix( const gsl_matrix* _data );
    Matrix( const double **_data, int nRows, int nCols );
    Matrix( const vector< vector< double > >& _data );
    void copy( const Matrix& other ) 
    { 
        if ( other.isEmpty() ) { data = NULL; return; }
        setDimensions( other.nRows(), other.nCols() ); 
        gsl_matrix_memcpy( data, other.data );
    }
    Matrix( const Matrix& other ) : data( NULL ) { copy( other ); }

    // destructor
    ~Matrix();
    
    // assignment
    Matrix& operator=( const Matrix& other ) { copy(other); return *this; }
			
    // access methods
    int nRows() const { if ( !isEmpty() ) return data->size1; else return 0; }
    int nCols() const { if ( !isEmpty() ) return data->size2; else return 0; }
    double getElement( int row, int col ) const { return gsl_matrix_get( data, row, col ); } 
    void setElement( int row, int col, double x ) { gsl_matrix_set( data, row, col, x ); }
    const double& operator()( int row, int col) const {
        assert( row >= 0 && row < data->size1  && col >= 0 && col < data->size2 );
        return *gsl_matrix_ptr( data, row, col );		
    }
    double& operator()( int row, int col ) {
        assert( row >= 0 && row < data->size1  && col >= 0 && col < data->size2 );
        return *gsl_matrix_ptr( data, row, col );
    }
    gsl_matrix* getData() const { return data; }
    void setDimensions( int nRows, int nCols);
    vector< double > getRow( int row ) const;
    vector< double > getCol( int col ) const;
    
    // information
    bool operator==( const Matrix& other ) const;
    bool isSquare() const;
    bool isEmpty() const { if ( data ) return false; else return true; }
    bool isSymmetric() const; 

    // set matrix elements
    void setZero() { gsl_matrix_set_zero( data ); }
    void setAll( double x ) { gsl_matrix_set_all( data, x ); }
    void setRow( int row, const vector< double >& v );
    void setCol( int col, const vector< double >& v );
    void setRows( const vector< double >& v );
    void setCols( const vector< double >& v );
    
    // set special matrices 
    void setIdentityMatrix( int n );		// identity matrix of dimension n 
    void setDiagonalMatrix( const vector< double >& diag );	// diagonal matrix

    // transpose
    Matrix transpose() const; 
    
    // addition
    Matrix operator+( const Matrix& other ) const;
    Matrix& operator+=( const Matrix& other );
    
    // multiplication of a constant
    Matrix operator*( const double c ) const;
    Matrix& operator*=( const double c ); 
    
    // matrix multiplication
    Matrix operator*( const Matrix& other ) const;	
    
    // output operator
    friend ostream& operator<<( ostream& os, const Matrix& m );
    
    // load Matrix from a file
    // If readDims = true, read dimensions (first row of the file); o/w, only read data, assuming the correct dimensions are known
    int load( const string& fileName, bool readDims = false );
    
    // save Matrix to a file
    void save( const string& fileName ) const;
};

/* IntMatrix class: integer matrix */
class IntMatrix {
    gsl_matrix_int* data;	// data
public:
    // constructors
    IntMatrix() : data( NULL ) {}
    IntMatrix( int nRows, int nCols );
    IntMatrix( int nRows, int nCols, int init );
    IntMatrix( const gsl_matrix_int* _data );
    IntMatrix( const int **_data, int nRows, int nCols );
    IntMatrix( const vector< vector< int > >& _data );
    void copy( const IntMatrix& other ) 
    { 
        if ( other.isEmpty() ) { data = NULL; return; }
        setDimensions( other.nRows(), other.nCols() ); 
        gsl_matrix_int_memcpy( data, other.data );
    }
    IntMatrix( const IntMatrix& other ) : data( NULL ) { copy( other ); }

    // destructor
    ~IntMatrix();
    
    // assignment
    IntMatrix& operator=( const IntMatrix& other ) { copy(other); return *this; }
                    
    // access methods
    int nRows() const { if ( !isEmpty() ) return data->size1; else return 0; }
    int nCols() const { if ( !isEmpty() ) return data->size2; else return 0; }
    int getElement( int row, int col ) const { return gsl_matrix_int_get( data, row, col ); } 
    void setElement( int row, int col, int v ) { gsl_matrix_int_set( data, row, col, v ); }
    const int& operator()( int row, int col) const {
        assert( row >= 0 && row < data->size1  && col >= 0 && col < data->size2 );
        return *gsl_matrix_int_ptr( data, row, col );		
    }
    int& operator()( int row, int col ) {
        assert( row >= 0 && row < data->size1  && col >= 0 && col < data->size2 );
        return *gsl_matrix_int_ptr( data, row, col );
    }
    gsl_matrix_int* getData() const { return data; }
    void setDimensions( int nRows, int nCols);
    vector< int > getRow( int row ) const;
    vector< int > getCol( int col ) const;
    
    // information
    bool operator==( const IntMatrix& other ) const;
    bool isSquare() const;
    bool isEmpty() const { if ( data ) return false; else return true; }
    bool isSymmetric() const; 

    // set matrix elements
    void setZero() { gsl_matrix_int_set_zero( data ); }
    void setAll( int x ) { gsl_matrix_int_set_all( data, x ); }
    void setRow( int row, const vector< int >& v );
    void setCol( int col, const vector< int >& v );
    void setRows( const vector< int >& v );
    void setCols( const vector< int >& v );
    
    // set special matrices 
    void setIdentityMatrix( int n );		// identity matrix of dimension n 
    void setDiagonalMatrix( const vector< int >& diag );	// diagonal matrix

    // transpose
    IntMatrix transpose() const; 
    
    // addition
    IntMatrix operator+( const IntMatrix& other ) const;
    IntMatrix& operator+=( const IntMatrix& other );
    
    // multiplication of a constant
    IntMatrix operator*( const int c ) const;
    IntMatrix& operator*=( const int c ); 
    
    // output operator
    friend ostream& operator<<( ostream& os, const IntMatrix& m );
    
    // load IntMatrix from a file
    // If readDims = true, read dimensions (first row of the file); o/w, only read data, assuming the correct dimensions are known
    int load( const string& fileName, bool readDims = false );
    
    // save IntMatrix to a file
    void save( const string& fileName ) const;
};

/*****************************************************
* Mathematical Functions 
******************************************************/

/* table used for log_add */
class log_add_table {
public:
    double delta;				// interval size
    vector< double > x_array;	// x[0], ..., x[N]
    vector< double > T;			// T[0], ..., T[N] where T[k] = log( 1 + exp( x[k] ) )
    
    // constructor
    log_add_table( double a, double b, int N );
};

/* log-add algorithm 
 * log_add( x, y ) = log( exp(x) + exp(y) ), etc. */
double log_add( double p, double q );
double log_add( const vector< double >& p );
double log_add( double p, double q, double r );
 
// log. of a vector 
vector< double > log( const vector< double >& v );

// exp. of a vector
vector< double > exp( const vector< double >& v );

// log. of a matrix
Matrix log( const Matrix& M );

// Eucledian distance bewteen two vectors
double Eucledian_dist( const vector< double >& x, const vector< double >& y );

// transformation of a number between (a,b) and (-inf, +inf)
double infty_transform( double x, double a, double b );     // (a,b) to (-inf, +inf)
double inverse_infty_transform( double z, double a, double b );     // (-inf, +inf) to (a,b)

// transformation involving weight parameters (w: 0 <= w_i <= 1) s.t. the new parameters in (-inf,+inf)
vector< double > weight_transform( const vector< double > w );		// w -> u
vector< double > inverse_weight_transform( const vector< double > u );		// u -> w

// finite difference method of gradient of a function
void numeric_deriv( gsl_vector* grad, double (*f)( const gsl_vector*, void* ), const gsl_vector* v, void* params, double step );

// sum
template< class T >
inline T sum( const vector< T >& v )
{
    T result = 0;
    for ( int i = 0; i < v.size(); i++ ) 
        result += v[ i ];
            
    return result;	
}

// max and min
double max( const vector< double >& v, int &arg );

// the threshold value at a given percentage in a vector of numbers
// if highest is true, then the value at the highest p; otherwise, the value at the lowest p
double elementAt( const vector< double >& v, double p, bool highest = true );

// rounding function: round x towards zero to the nearest integer (e.g. trunc (1.5) is 1.0 and trunc (-1.5) is -1.0 )
double trunc( double x );

// logit function & inverse logit function
double logit( double x );
double inv_logit( double x );

// logistic function
double logistic( double x );

/*****************************************************
* Statistics
******************************************************/
// mean of a vector
double mean( const vector< double >& x );

// median of a vector
double median( const vector< double >& x );

// standar deviation 
double std_dev( const vector< double >& x );

// Pearson correlation of v1 and v2, they must have equal sizes
double corr( const vector< double >& x, const vector< double >& y );

// Cross correlation of two time series (vectors): x and y. The vector t stores the values of the time lag to be evaluted. 
int cross_corr( const vector< double >& x, const vector< double >& y, const vector< int >& lag, vector< double >& cov, vector< double >& corr );

// least-square fit of two vectors: beta - coefficient, return RSS
double least_square( const vector< double >& x, const vector< double >& y, double& beta );

// check if a vector of real numbers is probablity mass function
bool isPmf( const vector< double > &p );

// sample one outcome from a multinomial distribution
int sampleMul( const gsl_rng * r, const vector< double > &p );

// sample a positive integer from a truncated geometric series: 1, r, r^2, ... r^(n - 1)
int sampleTruncatedGeometric( const gsl_rng* rng, double r, int n );

// mixture of two multinomial distributions, where the weight of the first one is "weight"
vector< double > multinomialMixture( const vector< double >& distr1, const vector< double >& distr2, double weight );

/*****************************************************
* String & I/O Functions 
******************************************************/

// tokenize/split a string
void tokenize( const string& str, vector< string >& tokens, const string& delimiters = " \t\n" );

// string to upper case
string toupperStr( const string& str );

// string to lower case
string tolowerStr( const string& str );

// convert a string in the format of "{bcd,Kr,hb}" (separated by , or space) into a vector of strings
void stringToVector( vector< string >& result, const string& str, const string& leftBoundary = "{", const string& rightBoundary = "}", const string& sep = " ," );

// convert a string in the format of "{0.3,0.2,0.2,0.3}" (separated by , or space) into a vector of numbers
void stringToVector( vector< double >& result, const string& str, const string& leftBoundary = "{", const string& rightBoundary = "}", const string& sep = " ," );

// read parameter file into a <field, value> table. Return 0 if successful, -1 otherwise
// Format: field = value [# ...] where # indicates the start of comments (ignored)
// Ex. kappa = 2.0 # kappa: the transition/transversion bias
int readParams( const string& file, map< string, string >& params, const char* comment = "#" );

// control the output format
class IO_Ctrl {
public:
    static string vector_delimiter;
    static string map_delimiter;	
    static string field_separator;
};

// output map< T1, T2 >
template< class T1, class T2 >
ostream& operator<<( ostream& os, const map< T1, T2 >& table )
{
    class map< T1, T2 >::const_iterator it;
    for ( it = table.begin(); it != table.end(); it++ ) { 
        os << it->first << IO_Ctrl::field_separator << it->second << IO_Ctrl::map_delimiter;    
    }

    return os;
}

// output vector< T >
template< class T >
ostream& operator<<( ostream& os, const vector< T >& v )
{
    if ( !v.empty() ) os << v[ 0 ];
    for ( int i = 1; i < v.size(); i++ ) {
        os << IO_Ctrl::vector_delimiter << v[ i ];
    }
    
    return os;	
}

// print matrix
template< class T >
void printMatrix( ostream& os, T** &x, int m, int n )
{	
    for ( int i = 0; i < m; i++ ) {
        for ( int j = 0; j < n; j++ ) {
                os << x[ i ][ j ] << " ";	
        }	
        os << endl;
    }	
}

// read vector
template< class T >
int readVector( const string& file, vector< T >& v )
{
    v.clear();
    
    ifstream fin( file.c_str() );
    if ( !fin ) {
        cerr << "Cannot open " << file << endl;
        exit( 1 );
    }	
    T elem;
    while ( fin >> elem ) v.push_back( elem );
    
    return 0;
} 	

// read matrix (2D array) data with the first column as labels
template< class T > 
int readMatrix( const string& file, vector< string >& labels, vector< vector< T > >& data )
{
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        labels.push_back( name );
        T val;
        vector< T > vals;
        while ( ss >> val ) vals.push_back( val );
        data.push_back( vals );
    }

    return 0;
}

// read matrix (2D array) data with the first row as column lables and the first column as row labels
template< class T > 
int readMatrix( const string& file, vector< string >& rowLabels, vector< string >& colLabels, vector< vector< T > >& data )
{
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }

    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    // read the first row: row labels (ignore the first field)
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    
    // read the data
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        T val;
        vector< T > vals;
        while ( ss >> val ) vals.push_back( val );
        data.push_back( vals );
    }

    return 0;
}

/*****************************************************
* Data Structure
******************************************************/

// transformation of an array to a vector
template< class T >
void array2vector( const T* arr, int n, vector< T >& v )
{
    v.clear();
    
    for ( int i = 0; i < n; i++ ) v.push_back( arr[ i ] );
}

// subset of a vector
template< class T >
int vectSubset( const vector< T >& v, const vector< int >& indices, vector< T >& results )
{
    results.clear();
    
    for ( int i = 0; i < indices.size(); i++ ) {
        if ( indices[ i ] >= v.size() ) return RET_ERROR;
        results.push_back( v[ indices[ i ] ] );	
    }	
    
    return 0;
}

// complement of an index set
int indexCompl( const vector< int >& indices, int n, vector< int >& complIndices );

#endif
