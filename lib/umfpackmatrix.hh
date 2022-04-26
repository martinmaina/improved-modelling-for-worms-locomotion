#ifndef UMFPACKMATRIX_HH
#define UMFPACKMATRIX_HH

#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>

struct UmfpackMatrixHolder
{
  UmfpackMatrixHolder( const int n_ = 0 )
    : n(n_), nz(0), Ap( n_+1, 0 )
  {}

  void clear()
  {
    n = 0;
    nz = 0;
    Ap.clear();
    Ai.clear();
    Ax.clear();
  }

  template< class Vector >
  void operator()( const Vector& in, Vector& out ) const
  {
    assert( in.size() == n );
    assert( out.size() == n );

    // clear out
    for( int i = 0; i < n; ++i )
      out[ i ] = 0;

    for( int i = 0; i < n; ++i )
      {
	for( int j = Ap.at( i ); j < Ap.at( i+1 ); ++j )
	  {
	    out[ Ai.at( j ) ] += Ax.at( j ) * in.at( i );
	  }
      }
  }


  int n;
  int nz;
  std::vector<int> Ap; //colptr
  std::vector<int> Ai; //row ind
  std::vector<double> Ax;  //value.d
};

class UmfpackMatrix
{
  typedef std::pair< unsigned int, double > ColValPairType;
  typedef std::vector< ColValPairType > RowType;
  typedef std::vector< RowType > DataType;

public:
  UmfpackMatrix( const unsigned int n )
    : data_( n ), n_( n )
  {}

  void add( const unsigned int rowIdx, const unsigned int colIdx, const double val )
  {
    RowType& row = data_.at( rowIdx );
    ColValPairType p( colIdx, val );

    for( auto& r : row )
      {
	if( r.first == p.first )
	  {
	    r.second += p.second;
	    return;
	  }
      }

    // not found
    row.push_back( p );
  }

  double get( const unsigned int rowIdx, const unsigned int colIdx ) const
  {
    const RowType& row = data_.at( rowIdx );

    for( const auto& r : row )
      {
	if( r.first == colIdx )
	  {
	    return r.second;
	  }
      }

    return 0;
  }

  void clear()
  {
    for( auto& row : data_ )
      {
	row.clear();
      }
  }

  void setDirichletRow( const unsigned int rowIdx )
  {
    bool doneDiagonal = false;
    for( unsigned int i = 0; i < n_; ++i )
      {
	for( auto& p : data_.at( i ) )
	  {
	    if( p.first == rowIdx )
	      {
		if( i == rowIdx )
		  {
		    p.second = 1.0;
		    doneDiagonal = true;
		  }
		else
		  {
		    p.second = 0.0;
		  }
	      }
	  }
      }

    if( not doneDiagonal )
      add( rowIdx, rowIdx, 1.0 );
  }

  template< class Vector >
  void operator()( const Vector& x, Vector& b ) const
  {
    // clear b
    for( unsigned int i = 0; i < n_; ++i )
      b.at( i ) = 0;

    for( unsigned int i = 0; i < n_; ++i )
      {
	const RowType& rowi = data_.at( i );
	for( const auto& p : rowi )
	  {
	    b.at( p.first ) += x.at( i ) * p.second;
	  }
      }
  }

  void print( std::ostream& s ) const
  {
    s << "***** WARNING: THIS IS THE TRANSPOSE *****" << std::endl;

    for( auto row : data_ )
      {
	if( row.empty() )
	  std::cout << "[empty row]";

	// first sort each row
	std::sort( row.begin(), row.end() );

	for( auto p : row )
	  s << " ( " << p.first << ", " << p.second << " )";
	s << std::endl;
      }
  }

  void matrix( UmfpackMatrixHolder& m ) const
  {
    m.clear();

    // set matrix size
    m.n = n_;

    int rowStartIdx = 0;
    m.Ap.push_back( rowStartIdx );
    for( auto row : data_ )
      {
	// first sort each row
	std::sort( row.begin(), row.end() );

	// add data to matrix holder
	for( const ColValPairType& pair : row )
	  {
	    m.Ai.push_back( pair.first );
	    m.Ax.push_back( pair.second );
	    rowStartIdx++;
	    m.nz++;
	  }

	// add end of row counter
	m.Ap.push_back( rowStartIdx );
      }
  }

  unsigned int n() { return n_; }

private:
  DataType data_;
  const unsigned int n_;
};

#endif
