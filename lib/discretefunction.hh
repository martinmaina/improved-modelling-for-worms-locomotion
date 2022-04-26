#ifndef DISCRETEFUNCTION_HH
#define DISCRETEFUNCTION_HH

#include <fstream>

#include "vector.hh"
#include "mesh.hh"

/**
 *  \class DiscreteFunctionInterface
 *
 *  \brief interface class for discrete functions
 *
 *  This class defines what it means to be a discrete function over a
 *  mesh.
 */
template< unsigned int mydim >
class DiscreteFunctionInterface
{
public:
  static const int dim = mydim;

  typedef RangeVector< dim > RangeVectorType;

  /**
   *  \brief constructor
   */
  explicit DiscreteFunctionInterface( const Mesh& mesh )
    : mesh_( mesh )
  {}

  /**
   *  \brief copy constructor
   */
  DiscreteFunctionInterface( const DiscreteFunctionInterface<dim>& other )
    : mesh_( other.mesh_ )
  {}

  /**
   *  \brief evaluation operator
   *
   *  evaluates discrete function at entry j
   *  \param[in]  j index
   *  \return range vector
   */
  virtual const RangeVectorType evaluate( const unsigned int j ) const = 0;

  /**
   *  \brief assignment operator
   *
   *  assigns discrete function at entry j
   *  \param[in]  j index
   *  \param[in]  val range vector
   */
  virtual void assign( const unsigned int j, const RangeVectorType& val ) = 0;

  /**
   *  \brief assignment operator
   *
   *  assigns discrete function at entry j
   *  \param[in]  j index
   *  \param[in]  val all values at set to val
   */
  virtual void assign( const unsigned int j, const double val ) = 0;

  /**
   *  \brief addition operator
   *
   *  addition operator for discrete function at entry j.
   *  the result should be data[ j ][ d ] = val [ d ]
   *  \param[in]  j index
   *  \param[in]  val range vector
   */
  virtual void add( const unsigned int j, const RangeVectorType& val ) = 0;

  /**
   *  \brief addition operator
   *
   *  addition operator for discrete function at entry j.
   *  the result should be data[ j ][ d ] = val
   *  \param[in]  j index
   *  \param[in]  val all dimension use val
   */
  virtual void add( const unsigned int j, const double val ) = 0;

  /**
   *  \brief scaling operator
   *
   *  scales discrete function by scalar
   *  \param val scale factor
   */
  virtual void operator*=( const double val ) = 0;

  /**
   *  \brief number of range vectors stored
   */
  virtual unsigned int N() const = 0;

  /**
   *  \brief number of doubles stored
   */
  virtual unsigned int size() const = 0;

  /**
   *  \brief access to individual entries
   */
  virtual double& operator[]( const unsigned int j ) = 0;

  /**
   *  \brief const access to individual entries
   */
  virtual double operator[]( const unsigned int j ) const = 0;

  /**
   *  \brief access to individual entries
   */
  virtual double& at( const unsigned int j )
  {
    return operator[]( j );
  }

  /**
   *  \brief const access to individual entries
   */
  virtual double at( const unsigned int j ) const
  {
    return operator[]( j );
  }

protected:
  /**
   *  \brief const access to the mesh
   */
  const Mesh& mesh() const
  {
    return mesh_;
  }

private:
  const Mesh& mesh_;
};

/**
 *  \class PiecewiseLinearFunction
 *
 *  \brief implementation of piecewise linear function
 *
 *  data is stored elsewhere
 */
template< unsigned int mydim >
class PiecewiseLinearFunction
  : public DiscreteFunctionInterface< mydim >
{
  typedef DiscreteFunctionInterface<mydim > BaseType;
  typedef SubArray< Vector > DataType;

public:
  static const unsigned int dim = BaseType :: dim;
  typedef typename BaseType :: RangeVectorType RangeVectorType;


/**
* \brief constructor
*
* \param[in] mesh, vector of type Vector
*/

  PiecewiseLinearFunction( const Mesh& mesh, DataType data )
    : BaseType( mesh ), data_( data )
  {}

  PiecewiseLinearFunction( const Mesh& mesh, Vector vector )
    : BaseType( mesh ), data_( vector )
  {}
 
/**
* \brief copy constructor
*
* \param[in] name of other function of type PiecewiseLinearFunction
*/

  PiecewiseLinearFunction( const PiecewiseLinearFunction& other )
    : BaseType( other ), data_( other.data_ )
  {}

/**
* \brief return vector at index j
*
* \param[in] index j
* \param[out] vector of type RangeVectorType
*
* \detailed returns range vector of length dim (e.g. dim=2) 
*  for index j. j represents length down the worm body(?).
*/

  virtual const RangeVectorType evaluate( const unsigned int j ) const
  {
    RangeVectorType ret;
    for( unsigned int d = 0; d < dim; ++d )
    {
      ret[ d ] = data_[ j * dim + d ];
    }
    return ret;
  }

/**
* \brief assign vector at index j
*
* \param[in] index j
* \param[in] val range vector
*/

  virtual void assign( const unsigned int j, const RangeVectorType& val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      const double v = val.at( d );
      data_[ j * dim + d ] = v;
    }
  }


/**
* \brief assign vector at index j
*
* \param[in] index j
* \param[in] val double input to data_ vector at j
*/

  virtual void assign( const unsigned int j, const double val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      data_[ j * dim + d ] = val;
    }
  }

  virtual void assign( const PiecewiseLinearFunction<mydim>& other )
  {
    for( unsigned int n = 0; n < size(); ++n )
    {
      data_[ n ] = other.data()[ n ];
    }
  }

/**
* \brief add double 'v' to range vector 'val' at index(s) j
*
* \param[in] index j
* \param[in] val range vector
*/

  virtual void add( const unsigned int j, const RangeVectorType& val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      const double v = val.at( d );
      data_[ j * dim + d ] += v;
    }
  }

/**
* \brief add double val to data_ at index(s) j
*
* \param[in] index j
* \param[in] val double
*/

  virtual void add( const unsigned int j, const double val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      data_[ j * dim + d ] += val;
    }
  }
/**
* \brief multiply data by double val
*
* \param[in] val double
*/
  virtual void operator*=( const double val )
  {
    for( auto& d : data_ )
      d *= val;
  }

  virtual unsigned int N() const
  {
    return mesh().N();
  }

/**
* \brief access value of data_ at index j
*
* \param[in] index j
* \param[out] double data_ point at j 
*/

  virtual double& operator[]( const unsigned int j )
  {
    return data_[ j ];
  }

  virtual double operator[]( const unsigned int j ) const
  {
    return data_[ j ];
  }

  typename DataType :: iterator begin()
  {
    return data_.begin();
  }

  typename DataType :: iterator end()
  {
    return data_.end();
  }

  typename DataType :: const_iterator begin() const
  {
    return data_.begin();
  }

  typename DataType :: const_iterator end() const
  {
    return data_.end();
  }

/**
* \detailed takes vector at meshpoint and outputs 
* according to output stream e.g. cout
*
* \param[in] output stream s
& \param[out] depends upon input s
*/
  void print( std::ostream& s ) const
  {
    for( unsigned int i = 0; i < N(); ++i )
    {
      RangeVectorType v = evaluate( i );
      for( auto vc : v )
      {
        s << vc << " ";
      }
      s << std::endl;
    }
  }

  unsigned int size() const
  {
    return data_.size();
  }

  void setDirichletBoundaries( const double val )
  {
    assign( 0, val );
    assign( N()-1, val );
  }

  using BaseType::mesh;

  DataType& data() { return data_; }
  const DataType& data() const { return data_; }

private:
  DataType data_;
};

/**
 *  \class PiecewiseConstantFunction
 *
 *  \brief implementation of piecewise constant function
 *
 *  data is stored elsewhere
 */
template< unsigned int mydim >
class PiecewiseConstantFunction
  : public DiscreteFunctionInterface< mydim >
{
  typedef DiscreteFunctionInterface<mydim > BaseType;
  typedef SubArray< Vector > DataType;

public:
  static const unsigned int dim = BaseType :: dim;
  typedef typename BaseType :: RangeVectorType RangeVectorType;

  PiecewiseConstantFunction( const Mesh& mesh, DataType data )
    : BaseType( mesh ), data_( data )
  {}

  PiecewiseConstantFunction( const PiecewiseConstantFunction& other )
    : BaseType( other ), data_( other.data_ )
  {}

/**
* \brief returns data at index j
*
* \param[in] index j
* \param[out] vector of data at j
*
* data_ is a one dimensional array
*/
  virtual const RangeVectorType evaluate( const unsigned int j ) const
  {
    RangeVectorType ret;
    for( unsigned int d = 0; d < dim; ++d )
    {
      ret[ d ] = data_[ j * dim + d ];
    }
    return ret;
  }

  virtual void assign( const unsigned int j, const RangeVectorType& val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      data_[ j * dim + d ] = val.at( j );
    }
  }

  virtual void assign( const unsigned int j, const double val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      data_[ j * dim + d ] = val;
    }
  }

  virtual void add( const unsigned int j, const RangeVectorType& val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      data_[ j * dim + d ] += val.at( j );
    }
  }

  virtual void add( const unsigned int j, const double val )
  {
    for( unsigned int d = 0; d < dim; ++d )
    {
      data_[ j * dim + d ] += val;
    }
  }

  virtual void operator*=( const double val )
  {
    for( auto& d : data_ )
      d *= val;
  }

  virtual unsigned int N() const
  {
    return mesh().N()-1;
  }

  virtual double& operator[]( const unsigned int j )
  {
    return data_[ j ];
  }

  virtual double operator[]( const unsigned int j ) const
  {
    return data_[ j ];
  }

  typename DataType :: iterator begin()
  {
    return data_.begin();
  }

  typename DataType :: iterator end()
  {
    return data_.end();
  }

  typename DataType :: const_iterator begin() const
  {
    return data_.begin();
  }

  typename DataType :: const_iterator end() const
  {
    return data_.end();
  }

  unsigned int size() const
  {
    return data_.size();
  }

  void print( std::ostream& s ) const
  {
    for( unsigned int i = 0; i < N(); ++i )
    {
      RangeVectorType v = evaluate( i );
      for( auto vc : v )
      {
        s << vc << " ";
      }
      s << std::endl;
    }
  }

protected:
  using BaseType::mesh;

private:
  DataType data_;
};

template< unsigned int mydim >
class PositionCurvaturePressureTuple
  : public Vector
{
  typedef Vector BaseType;

public:
  static const unsigned int dim = mydim;

  typedef PiecewiseLinearFunction< dim > PositionFunctionType;    // position and curvature are piecewise linear functions
  typedef PiecewiseLinearFunction< dim > CurvatureFunctionType;
  typedef PiecewiseConstantFunction< 1 > PressureFunctionType;    // pressure is a piecewise constant function

  typedef RangeVector< dim > RangeVectorType;

  typedef SubArray< Vector > SubArrayType;

  PositionCurvaturePressureTuple( const Mesh& mesh )
    : BaseType( ( 2*dim * mesh.N() + (mesh.N()-1) ) ), mesh_( mesh ),
      position_( mesh, SubArrayType( begin(), begin()+dim*mesh.N() ) ),
      curvature_( mesh, SubArrayType( begin()+dim*mesh.N(), begin()+2*dim*mesh.N() ) ),
      pressure_( mesh, SubArrayType( begin()+2*dim*mesh.N(), end() ) )
  {}

/**
* \brief Getters: return private members position_, curvature_, pressure_
*
* \param[out] private members of type piecewiselinear, or piecewise constant, function
*/
  const PositionFunctionType& position() const { return position_; }
  PositionFunctionType& position() { return position_; }
  const CurvatureFunctionType& curvature() const { return curvature_; }
  CurvatureFunctionType& curvature() { return curvature_; }
  const PressureFunctionType& pressure() const { return pressure_; }
  PressureFunctionType& pressure() { return pressure_; }

  void assign( const PositionCurvaturePressureTuple<mydim> &other )
  {
    for( unsigned int i = 0; i < size(); ++i )
    {
      at( i ) = other.at( i );
    }
  }

private:
  const Mesh& mesh_;

  PositionFunctionType position_;
  CurvatureFunctionType curvature_;
  PressureFunctionType pressure_;
};

#endif // #ifndef DISCRETEFUNCTION_HH
