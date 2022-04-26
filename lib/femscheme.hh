#ifndef FEMSCHEME_HH
#define FEMSCHEME_HH

#include <iomanip>
#include "mesh.hh"
#include "discretefunction.hh"
#include "matrix.hh"
#include "solver.hh"
#include "timeprovider.hh"

template< class Model, const int wd >
struct FemScheme
{
  using ModelType = Model;
  static const int worlddim = wd;

  using SolutionTupleType = PositionCurvaturePressureTuple< worlddim >;
  using PositionFunctionType = typename SolutionTupleType :: PositionFunctionType;
  using CurvatureFunctionType = typename SolutionTupleType :: CurvatureFunctionType;
  using PressureFunctionType = typename SolutionTupleType :: PressureFunctionType;

  using RangeVectorType = typename PositionFunctionType :: RangeVectorType;

  using SystemOperatorType = SystemMatrix< PositionFunctionType, ModelType >;

  FemScheme( const Mesh& mesh, const ModelType& model )
    : mesh_( mesh ),
      model_( model ),
      // solution variables
      solutionTuple_( mesh ),
      oldSolutionTuple_( mesh ),
      rhsTuple_( mesh ),
      // operators
      implicitOperator_( mesh, oldSolutionTuple_.position(), oldSolutionTuple_.curvature(), model_, true ),
      explicitOperator_( mesh, oldSolutionTuple_.position(), oldSolutionTuple_.curvature(), model_, false )
      //Y_out(NULL)  
  {
      /*  Y_out = new std::ofstream("Y_out.txt");
	if( !Y_out->is_open() )
        {
          std::cerr << "Error: Y_out.txt didn't open" << std::endl;
        }*/
  }

/*~FemScheme() 
{
  Y_out->close();
  if( Y_out )
    delete Y_out;
}*/

  const PositionFunctionType& position() const { return solutionTuple_.position(); }
  const CurvatureFunctionType& curvature() const { return solutionTuple_.curvature(); }  //    const auto& Y = solutionTuple_.curvature(); from "energy" function
  const PressureFunctionType& pressure() const { return solutionTuple_.pressure(); }
  const PositionFunctionType& oldPosition() const { return oldSolutionTuple_.position(); }

  void initialiseXY()
  {
    RangeVectorType Xu(0);
    for( unsigned int j = 0; j < mesh_.N(); ++j )
      {
	const double uj = mesh_.u(j);
	const auto Xuj = model_.X0( uj );
	solutionTuple_.position().assign( j, Xuj );
      }

    // compute curvature
    computeCurvature();
  }

  void initialiseXY( const RangeVectorType& middle )
  {
    RangeVectorType Xu(0);
    for( unsigned int j = 0; j < mesh_.N(); ++j )
      {
	for( unsigned int d = 0; d < worlddim; ++d )
	  {
	    Xu[d] = middle[d] + model_.gamma( mesh_.u(j) ) * ( mesh_.u(j) - 0.5 ) / std::sqrt(3.0);
	  }

	solutionTuple_.position().assign( j, Xu );
      }


    // compute curvature
    computeCurvature();
  }

  void initialiseXY( const PositionFunctionType& old )
  {
    // assign from old
    solutionTuple_.position().assign( old );

    // compute curvature
    computeCurvature();
  }

  void computeCurvature()
  {
    // construct mass and stiffness matrices of size worlddim*N
    typedef InverseMassMatrix< PositionFunctionType > InverseMassMatrixType;
    InverseMassMatrixType inverseMass( mesh_, position() );

    typedef StiffnessMatrix< PositionFunctionType > StiffnessMatrixType;
    StiffnessMatrixType stiff( mesh_, position() );

    // find curvature
    Vector rhsData( worlddim * mesh_.N() );
    PositionFunctionType rhs( mesh_, SubArray< Vector >( rhsData.begin(), rhsData.end() ) ); 
    stiff( position(), rhs );
    rhs *= -1;

    // set dirichlet boundary conditions
    rhs.setDirichletBoundaries( 0 );
    inverseMass.setDirichletBoundaries();
    inverseMass( rhs, solutionTuple_.curvature() );
  }


  void prepare()
  {
    // update old solution
    oldSolutionTuple_.assign( solutionTuple_ );

    // construct rhs
    explicitOperator_.assemble();
    explicitOperator_( oldSolutionTuple_, rhsTuple_ );

    // set Dirichlet conditions
    rhsTuple_.curvature().setDirichletBoundaries( 0 );
  }

  void solve()
  {
    // assemble system operator
    implicitOperator_.assemble();

    // solve
    umfpackSolve( implicitOperator_, solutionTuple_, rhsTuple_ );
  }

  /*void output_Y()
  {
	const auto& X = solutionTuple_.position();
	const auto& Y = solutionTuple_.curvature();
	
	for(unsigned int e = 0; e < X.mesh().N() - 1; ++e)
	{
	*Y_out << " " << Y.evaluate( e )[0] << " " << Y.evaluate( e )[1];
 	}

	*Y_out << "\n";
  }*/

  double energy( const bool verbose = false ) const
  {
    const auto& X = solutionTuple_.position();
    const auto& Y = solutionTuple_.curvature();

    const unsigned int N = X.mesh().N();
    double ret = 0;

    for( unsigned int e = 0; e < N-1; ++e )
      {
	// find element geometry
	const auto xe = X.evaluate( e );
	const auto xep = X.evaluate( e+1 );
	const auto xd = xe - xep;
	const double q = xd.norm();

	// find curvature at point
	const auto ye = Y.evaluate( e );

	// find energy density
	const double ed = model_.energyDensity( xe, ye );

	ret += ed * q;
      }

    // extra energy from model
    ret += model_.extraEnergy( X );

    if( verbose )
      std::cout << "energy: " << ret << std::endl;

    return ret;
  }

  double energy( std::ofstream& file, const bool verbose = false ) const
  {
    file << model_.timeProvider().time() << ",";

    const auto& X = solutionTuple_.position();
    const auto& Y = solutionTuple_.curvature();

    const unsigned int N = X.mesh().N();
    double ret = 0;

    std::vector<double> energyTuple( model_.nCam+1,0.0);
    std::vector<double> list;

    for( unsigned int e = 0; e < N-1; ++e )
      {
	// find element geometry
	const auto xe = X.evaluate( e );
	const auto xep = X.evaluate( e+1 );
	const auto xd = xe - xep;
	const double q = xd.norm();

	// find curvature at point
	const auto ye = Y.evaluate( e );

	// find energy density
	const double ed = model_.energyDensity( xe, ye, list );

	ret += ed * q;

	for( unsigned int d = 0; d < list.size(); ++d )
	  energyTuple.at(d) += list.at(d) * q;
      }

    for( auto e : energyTuple )
      file << e << ",";

    // extra energy from model
    ret += model_.extraEnergy( X );

    file << model_.extraEnergy( X );

    if( verbose )
      std::cout << "energy: " << ret << std::endl;

    file << "," << updateSize();
    file << std::endl;

    return ret;
  }

  double updateSize() const
  {
    const auto& X = solutionTuple_.position();
    const auto& Xold = oldSolutionTuple_.position();

    const unsigned int N = X.mesh().N();

    double ret = 0;
    for( unsigned int e = 0; e < N-1; ++e )
      {
	// find element geometry
	const auto xe = X.evaluate( e );
	const auto xep = X.evaluate( e+1 );
	const auto xd = xe - xep;
	const double q = xd.norm();

	const auto xeOld = Xold.evaluate( e );
	const auto xepOld = Xold.evaluate( e+1 );

	ret += ( ( xeOld - xe ).norm2() + ( xepOld - xep ).norm2() ) * q;
      }

    return std::sqrt( ret );
  }

  std::pair<double,double> constraintError() const
  {
    const auto& X = solutionTuple_.position();
    const unsigned int N = X.mesh().N();

    std::pair< double, double > minMax = { 1e10, 0.0 };
    for( unsigned int e = 0; e < N-1; ++e )
      {
	// find element geometry
	const auto xe = X.evaluate( e );
	const auto xep = X.evaluate( e+1 );
	const auto xd = xe - xep;
	const double q = xd.norm() / X.mesh().h( e+1 );

	minMax.first = std::min( minMax.first, q );
	minMax.second = std::max( minMax.second, q );
      }

    return minMax;
  }


private:
  const Mesh& mesh_;
  const ModelType& model_;

  SolutionTupleType solutionTuple_;
  SolutionTupleType oldSolutionTuple_;
  SolutionTupleType rhsTuple_;

  SystemOperatorType implicitOperator_;
  SystemOperatorType explicitOperator_;

 // std::ofstream* Y_out;

};

template< class Model, const int wd >
struct BaseFemScheme
{
  using ModelType = Model;
  static const int worlddim = wd;

  using SolutionTupleType = PositionCurvaturePressureTuple< worlddim >;
  using PositionFunctionType = typename SolutionTupleType :: PositionFunctionType;
  using CurvatureFunctionType = typename SolutionTupleType :: CurvatureFunctionType;
  using PressureFunctionType = typename SolutionTupleType :: PressureFunctionType;

  using RangeVectorType = typename PositionFunctionType :: RangeVectorType;

  using SystemOperatorType = SystemMatrix< PositionFunctionType, ModelType >;

  BaseFemScheme( const Mesh& mesh, const ModelType& model )
    : mesh_( mesh ),
      model_( model ),
      // solution variables
      solutionTuple_( mesh ),
      oldSolutionTuple_( mesh ),
      rhsTuple_( mesh ),
      // operators
      implicitOperator_( mesh, oldSolutionTuple_.position(), oldSolutionTuple_.curvature(), model_, true ),
      explicitOperator_( mesh, oldSolutionTuple_.position(), oldSolutionTuple_.curvature(), model_, false )
  {}

  const PositionFunctionType& position() const { return solutionTuple_.position(); }
  const CurvatureFunctionType& curvature() const { return solutionTuple_.curvature(); }
  const PressureFunctionType& pressure() const { return solutionTuple_.pressure(); }
  const PositionFunctionType& oldPosition() const { return oldSolutionTuple_.position(); }

  void initialiseXY()
  {
    RangeVectorType Xu(0);
    for( unsigned int j = 0; j < mesh_.N(); ++j )
      {
	const double uj = mesh_.u(j);
	const auto Xuj = model_.X0( uj );
	solutionTuple_.position().assign( j, Xuj );
      }

    // compute curvature
    computeCurvature();
  }

  void computeCurvature()
  {
    // construct mass and stiffness matrices of size worlddim*N
    typedef InverseMassMatrix< PositionFunctionType > InverseMassMatrixType;
    InverseMassMatrixType inverseMass( mesh_, position() );

    typedef StiffnessMatrix< PositionFunctionType > StiffnessMatrixType;
    StiffnessMatrixType stiff( mesh_, position() );

    // find curvature
    Vector rhsData( worlddim * mesh_.N() );
    PositionFunctionType rhs( mesh_, SubArray< Vector >( rhsData.begin(), rhsData.end() ) );
    stiff( position(), rhs );
    rhs *= -1;

    // set dirichlet boundary conditions
    rhs.setDirichletBoundaries( 0 );
    inverseMass.setDirichletBoundaries();
    inverseMass( rhs, solutionTuple_.curvature() );
  }

  void prepare()
  {
    // update old solution
    oldSolutionTuple_.assign( solutionTuple_ );

    // construct rhs
    explicitOperator_.assemble();
    explicitOperator_( oldSolutionTuple_, rhsTuple_ );

    // set Dirichlet conditions
    rhsTuple_.curvature().setDirichletBoundaries( 0 );
  }

  void solve()
  {
    // assemble system operator
    implicitOperator_.assemble();

    // solve
    umfpackSolve( implicitOperator_, solutionTuple_, rhsTuple_ );
  }

private:
  const Mesh& mesh_;
  const ModelType& model_;

  SolutionTupleType solutionTuple_;
  SolutionTupleType oldSolutionTuple_;
  SolutionTupleType rhsTuple_;

  SystemOperatorType implicitOperator_;
  SystemOperatorType explicitOperator_;
};

#endif // #ifndef FEMSCHEME_HH
