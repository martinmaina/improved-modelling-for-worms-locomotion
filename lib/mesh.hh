#ifndef MESH_HH
#define MESH_HH

#include <cassert>

struct Mesh
{
  Mesh( const double a, const double b, const unsigned int N )
    : a_( a ), b_( b ), N_( N )
  {}

  double h( unsigned int j ) const
  {
    assert( j > 0 and j < N_ );
    return u( j ) - u( j - 1 );
  }

  double u( unsigned int j ) const
  {
    return a_ + j * ( b_ - a_ ) / ( N_ - 1);
  }

  double a() const { return a_; }
  double b() const { return b_; }
  unsigned int N() const { return N_; }

private:
  const double a_, b_;
  const unsigned int N_;
};

#endif
