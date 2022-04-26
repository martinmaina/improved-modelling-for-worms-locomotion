#ifndef MATRIX_HH
#define MATRIX_HH

#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

#include "mesh.hh"
#include "parameter.hh"
#include "timeprovider.hh"
#include "model.hh"
#include "umfpackmatrix.hh"


class MeshMatrix
  : public UmfpackMatrix
{
  typedef UmfpackMatrix BaseType;

public:
  MeshMatrix( const Mesh& mesh )
    : BaseType( mesh.N() ), mesh_( mesh )
  {}

  MeshMatrix( const Mesh& mesh, const unsigned int n )
    : BaseType( n ), mesh_( mesh )
  {}

  void setDirichletBoundaries()
  {
    const int dim = n() / mesh().N();

    for( int d = 0; d < dim; ++d )
      {
	setDirichletRow( 0*dim + d );
	setDirichletRow( (mesh().N()-1)*dim + d );
      }
  }

  virtual void assemble() = 0;
  void assemble( const bool ) { std::cerr << "not me!" << std::endl; }

protected:
  const Mesh& mesh() const
  {
    return mesh_;
  }

private:
  const Mesh& mesh_;
};

template< class PositionMeshFunction >
class MassMatrix
  : public MeshMatrix
{
  typedef typename PositionMeshFunction :: RangeVectorType RangeVectorType;

public:
  MassMatrix( const Mesh& mesh, const PositionMeshFunction& X )
    : MeshMatrix( mesh, X.dim * mesh.N() ), X_( X )
  {
    assemble();
  }

  void assemble()
  {
    // clear first
    clear();

    // element loop
    for( unsigned int e = 0; e < mesh().N()-1; ++e )
      {
	// construct area element
	const RangeVectorType xe = X_.evaluate( e );
	const RangeVectorType xep = X_.evaluate( e+1 );
	const double q = ( xe - xep ).norm();

	// find value
	const double value = 0.5 * q;

	// add to vertex in this element
	for( unsigned int d = 0; d < X_.dim; ++d )
	  {
	    for( unsigned int ie = 0; ie < 2; ++ie )
	      {
		add( X_.dim*(ie+e) + d, X_.dim*(ie+e) + d, value );
	      }
	  }
      }
  }

  void setDirichletBoundaries()
  {
    for( unsigned int d = 0; d < X_.dim; ++d )
      {
	setDirichletRow( d );
	setDirichletRow( mesh().N()-1+d );
      }
  }


private:
  const PositionMeshFunction& X_;
};

template< class PositionMeshFunction >
class InverseMassMatrix
  : public MeshMatrix
{
  typedef typename PositionMeshFunction :: RangeVectorType RangeVectorType;

public:
  InverseMassMatrix( const Mesh& mesh, const PositionMeshFunction& X )
    : MeshMatrix( mesh, X.dim * mesh.N() ), X_( X )
  {
    assemble();
  }

  void assemble()
  {
    // clear first
    clear();

    // element loop
    for( unsigned int e = 0; e < mesh().N()-1; ++e )
      {
	// construct area element
	const RangeVectorType xe = X_.evaluate( e );
	const RangeVectorType xep = X_.evaluate( e+1 );
	const double q = ( xe - xep ).norm();

	// find value
	double value = 0.5 * q;

	if( e > 0 )
	  {
	    // construct previous area element
	    const RangeVectorType xem = X_.evaluate( e-1 );
	    const double qm = ( xe - xem ).norm();

	    // find value
	    value += 0.5 * qm;
	  }


	// add to vertex in this element
	for( unsigned int d = 0; d < X_.dim; ++d )
	  {
	    add( X_.dim*e + d, X_.dim*e + d, 1.0/ value );
	  }
      }
  }

private:
  const PositionMeshFunction& X_;
};

template< class PositionMeshFunction >
class StiffnessMatrix
  : public MeshMatrix
{
  typedef typename PositionMeshFunction :: RangeVectorType RangeVectorType;

public:
  StiffnessMatrix( const Mesh& mesh, const PositionMeshFunction& X )
    : MeshMatrix( mesh, X.dim * mesh.N() ), X_( X )
  {
    assemble();
  }

  void assemble()
  {
    // clear first
    clear();

    // element loop
    for( unsigned int e = 0; e < mesh().N()-1; ++e )
      {
	// construct area element
	const RangeVectorType xe = X_.evaluate( e );
	const RangeVectorType xep = X_.evaluate( e+1 );
	const double q = ( xe - xep ).norm();

	// find value
	const double value = 1.0 / q;

	// add to vertex in this element
	for( unsigned int d = 0; d < X_.dim; ++d )
	  {
	    // diagonal terms
	    for( unsigned int ie = 0; ie < 2; ++ie )
	      {
		add( X_.dim*(ie+e) + d, X_.dim*(ie+e) + d, value );
	      }

	    // off diagonal terms
	    add( X_.dim*(e) + d, X_.dim*(e+1) + d, -value );
	    add( X_.dim*(e+1) + d, X_.dim*(e) + d, -value );
	  }
      }
  }

private:
  const PositionMeshFunction& X_;
};

template < class PositionMeshFunction, class Model >
class SystemMatrix
  : public MeshMatrix
{
public:
  typedef typename PositionMeshFunction :: RangeVectorType RangeVectorType;
  static const int dim = RangeVectorType::dim;
  using ModelType = Model;

  SystemMatrix( const Mesh& mesh, const PositionMeshFunction& X,
		const PositionMeshFunction& Y,
		const ModelType& model, const bool implicit )
    : MeshMatrix( mesh, 2*dim*mesh.N() + mesh.N()-1 ),
      X_( X ), Y_( Y ),
      model_( model ),
      implicit_( implicit ),
      assembled_( false )
  {}


  template< class Vector >
  void operator()( const Vector& x, Vector& b ) const
  {
    if( not assembled_ )
      throw "system matrix not assembled";

    MeshMatrix::operator()( x, b );

    if( not implicit_ )
      {
	// extract sub-vectors
	const auto& X = x.position();
	auto& rhsX = b.position();
	auto& rhsY = b.curvature();
	auto& rhsP = b.pressure();

	// element loop
	for( unsigned int e = 0; e < mesh().N()-1; ++e )
	  {
	    // extract geometry
	    const double s = 0.5 * ( mesh().u( e ) + mesh().u( e+1 ) );
	    const auto xe = X.evaluate( e );
	    const auto xep = X.evaluate( e+1 );
	    const double q = ( xe - xep ).norm();
	    RangeVectorType nu = ( xe - xep ).perp();
	    nu /= q;
//std::cout << "  " << atan(nu[1]/nu[0]) << std::endl;
	    // integrate position and curvature forcing
	    for( unsigned int ie = 0; ie < 2; ++ie )
	      {
		const auto xie = X.evaluate( e + ie );
		const auto sie = mesh().u( e + ie );  // mesh().u(j) is jth element down midline

		RangeVectorType f = model().f( xie );

		RangeVectorType valueX;
		for( unsigned int d = 0; d < X.dim; ++d )
		  {
		    valueX[ d ] = 0.5 * q * f[ d ];
		  }

	//	assert(Y_.evaluate( e )[ 0 ] - b.curvature().evaluate( e )[ 0 ]+ model().beta( sie ) * nu[ 0 ] < 1e-5); // kappa = b.curvature().evaluate( e )[ 0 ]          model().scalarCurvature( 0 )
		rhsX.add( e+ie, valueX );

		RangeVectorType valueY;
		for( unsigned int d = 0; d < X.dim; ++d )
		  {
		    valueY[ d ] = 0.5 * q * model().beta( sie ) * nu[ d ];
		  }

		rhsY.add( e+ie, valueY );
	      }

	    // integrate pressure forcing
	    if( model().constraint() == ModelType :: Constraint :: measure )
	      rhsP[e] = mesh().h( e+1 ) * -model().gamma( s );
	  }

	// extra forcing from model
	model().extraForcing( X, rhsX );
      }
  }
/*
  void TestScale()
  {
int e = 0;
	    // extract geometry
	    const double s = 0.5 * ( mesh().u( e ) + mesh().u( e+1 ) );
	    const auto xe = X.evaluate( e );
	    const auto xep = X.evaluate( e+1 );
	    const double q = ( xe - xep ).norm();
	    RangeVectorType nu = ( xe - xep ).perp();
	    nu /= q;

	
	auto equal = Y_.evaluate( e )[0] - curvature_.evaluate( e )[0] + nu[0]*model().beta( e );
  }
*/
  void assemble()
  {
    // clear first
    clear();

    // time step
    const double deltaT = model().timeProvider().deltaT();

    // element loop
    for( unsigned int e = 0; e < mesh().N()-1; ++e )
      {
	// model parameters
	const double s = 0.5 * ( mesh().u( e ) + mesh().u( e+1 ) );
	const double K = model().K( s );
	const double E = model().e( s );

	// construct area element and tangent
	const RangeVectorType xe = X_.evaluate( e );
	const RangeVectorType xep = X_.evaluate( e+1 );
	const double q = ( xep - xe ).norm();
	RangeVectorType tau = ( xep - xe );
	tau /= q;
	assert(q==q);
	assert(q);

//std::cout << "e: " << e << "      tau: " << tau << std::endl;
	const RangeVectorType ye = Y_.evaluate( e );
	const RangeVectorType yep = Y_.evaluate( e+1 );

	const double energyDensity = model().energyDensity( (xe + xep)*0.5,
							    (ye + yep)*0.5 );

	// add to vertex in this element
	for( unsigned int k = 0; k < X_.dim; ++k )
	  {
	    for( unsigned int l = 0; l < X_.dim; ++l )
	      {
		// projection matrix
		const double Ptau = ( (double) (k == l) - tau[l]*tau[k] );

		// drag term
		// if( k == l )
		  {
		    const double id = (double) (k == l);
		    const double hilf = ( 1.0 - K ) * tau[l]*tau[k];
		    const double drag = 0.5 * q * ( K * id +  hilf ) / deltaT;

		    // add to matrix
		    add( dim*(e) + l, dim*(e) + k, drag );
		    add( dim*(e+1) + l, dim*(e+1) + k, drag );
		  }

		// implicit terms
		if( implicit_ )
		  {
		    // projected stiffness term
		    const double pStiff = -Ptau * E / q;

		    // add to matrix
		    const unsigned int pStiffOffset = dim*mesh().N();
		    add( pStiffOffset + dim*(e) + l, dim*(e) + k, pStiff );
		    add( pStiffOffset + dim*(e+1) + l, dim*(e+1) + k, pStiff );
		    add( pStiffOffset + dim*(e) + l, dim*(e+1) + k, -pStiff );
		    add( pStiffOffset + dim*(e+1) + l, dim*(e) + k, -pStiff );

		    // projected stiffness term
		    const double wStiff = energyDensity / q;
		    // add to matrix
		    add( dim*(e) + l, dim*(e) + k, wStiff );
		    add( dim*(e+1) + l, dim*(e+1) + k, wStiff );
		    add( dim*(e) + l, dim*(e+1) + k, -wStiff );
		    add( dim*(e+1) + l, dim*(e) + k, -wStiff );

		    // mass matrix
		    if( k == l )
		      {
			const double mass = 0.5 * q;
			// add to matrix
			const unsigned int massOffset = dim*mesh().N();
			add( massOffset + dim*(e) + l, massOffset + dim*(e) + k, mass );
			add( massOffset + dim*(e+1) + l, massOffset + dim*(e+1) + k, mass );
		      }

		    // stiffness matrix
		    if( k == l )
		      {
			const double stiff = 1.0 / q;
			const unsigned int stiffOffset = dim*mesh().N();
			add( dim*(e) + l, stiffOffset + dim*(e) + k, stiff );
			add( dim*(e+1) + l, stiffOffset + dim*(e+1) + k, stiff );
			add( dim*(e) + l, stiffOffset + dim*(e+1) + k, -stiff );
			add( dim*(e+1) + l, stiffOffset + dim*(e) + k, -stiff );

			if( e == 0 )
			  {
			    add( dim*(e) + l, stiffOffset + dim*(e) + k, -1/q );
			    add( dim*(e+1) + l, stiffOffset + dim*(e) + k, 1/q );
			  }
			if( e == mesh().N()-2 )
			  {
			    add( dim*(e) + l, stiffOffset + dim*(e+1) + k, 1/q );
			    add( dim*(e+1) + l, stiffOffset + dim*(e+1) + k, -1/q );
			  }
		      }
		  }
	      }

	    // non symmetric terms
	    const double pres = tau[k];
	    // add to matrix
	    const unsigned int presOffset = 2 * dim * mesh().N();
	    if( implicit_ or model().constraint() == ModelType :: Constraint :: velocity )
	      {
		add( presOffset + e, dim*e + k, pres );
		add( presOffset + e, dim*(e+1) + k, -pres );

		add( dim*e + k, presOffset + e, pres );
		add( dim*(e+1) + k, presOffset + e, -pres );
	      }
	  }
      }

    // set flag
    assembled_ = true;
  }

protected:
  const ModelType& model() const { return model_; }

private:
  const PositionMeshFunction& X_;
  const PositionMeshFunction& Y_;
  const ModelType& model_;
  const bool implicit_;

  bool assembled_;
};


#endif // #ifndef MATRIX_H
