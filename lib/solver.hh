#ifndef SOLVER_HH
#define SOLVER_HH

#include <cassert>
#include <exception>
#include <sstream>

extern "C" {
#include <umfpack.h>
}

#include "umfpackmatrix.hh"

struct UmfpackError : public std::exception
{
  UmfpackError( const int& e )
    : e_( e )
  {}

  virtual const char* what() const throw()
  {
    switch( e_ )
      {
      case UMFPACK_OK:
	return "UMFPACK was successful";
      case UMFPACK_WARNING_singular_matrix:
	return "Matrix is singular. There are exact zeros on the diagonal of U";
      case UMFPACK_WARNING_determinant_underflow:
	return "The determinant is nonzero, but smaller in magnitude than the smallest positive floating-point number.";


      case UMFPACK_ERROR_n_nonpositive:
      case UMFPACK_ERROR_out_of_memory:
	return "\
Insufficient memory to perform the symbolic analysis. If the \
analysis requires more than 2GB of memory and you are using \
the 32-bit (\"int\") version of UMFPACK, then you are guaranteed \
to run out of memory. Try using the 64-bit version of UMFPACK.";
      case UMFPACK_ERROR_argument_missing:
	return "One or more required arguments is missing.";
      case UMFPACK_ERROR_internal_error:
	return "\
Something very serious went wrong. This is a bug. \
Please contact the author (davis@cise.ufl.edu).";
      default:
	std::stringstream ss;
	ss << "unknown error: " << e_;
	return ss.str().c_str();
      }
  }
private:
  const int e_;
};

template< class Vector >
void umfpackSolve( const UmfpackMatrix& A, Vector&x, const Vector& b )
{
  // get matrix data
  UmfpackMatrixHolder m;
  A.matrix( m );

  // check sizes
  assert( x.size() == m.n );
  assert( b.size() == m.n );

  // helper variables
  int status;
  double Control [UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric ;

  // set defaults
  umfpack_di_defaults(Control);

  // symbolic factorisation
  status = umfpack_di_symbolic ( m.n, m.n, &*(m.Ap.begin()), &*(m.Ai.begin()), &*(m.Ax.begin()),
				 &Symbolic, Control, Info) ;
  if( status != UMFPACK_OK )
    throw UmfpackError( status );

  // numeric factorisation
  status = umfpack_di_numeric ( &*(m.Ap.begin()), &*(m.Ai.begin()), &*(m.Ax.begin()),
				Symbolic, &Numeric, Control, Info) ;
  if( status != UMFPACK_OK )
    throw UmfpackError( status );

  // free symbolic factorisation
  umfpack_di_free_symbolic (&Symbolic) ;

  // solve
  status = umfpack_di_solve (UMFPACK_A, &*(m.Ap.begin()), &*(m.Ai.begin()), &*(m.Ax.begin()),
			     &*(x.begin()), &*(b.begin()), Numeric, Control, Info) ;
  if( status != UMFPACK_OK )
    throw UmfpackError( status );

  // free numerical factorisation
  umfpack_di_free_numeric (&Numeric) ;

  // check solve
  Vector res( x );
  A( x, res );
  res.axpy(-1,b);
 

//for( auto bb : b )
//std::cout << " " << bb;
//std::cout << std::endl;

  if( not ( res.norm2() < 1.0e-10 ) )
    {
      std::cerr << "residual = " << res.norm2() << std::endl;
      throw "solver error";
    }

}

#endif
