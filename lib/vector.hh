#ifndef VECTOR_HH
#define VECTOR_HH

#ifdef OpenCV_FOUND
#include "opencv2/opencv.hpp"
#endif

#include <array>
#include <vector>
#include <cmath>

struct Vector : public std::vector<double>
{
  // typedef the Base type and the type of this class 
  // since they are used below
  typedef std::vector<double> Base;
  typedef Vector This;
  
  // construction of a vector providing the size,
  // note that the vector values are not initialized
  explicit Vector(const unsigned int N) : Base(N) {}
  // a second constructor is provided which also initializes 
  // the values
  Vector(const unsigned int N, const double a) : Base(N,a) {}
  // also add a copy constructor (although the defaul would do the same)
  Vector(const Vector &other) : Base(other) {}

  // compute the scalar product of this vector with a second one
  double operator*(const This &y) const
  {
    double ret = 0;
    for (unsigned int i=0;i<size();++i)
      ret += (*this)[i]*y[i];
    return ret;
  }
  // axpy operation this -> a * x + this 
  void axpy(double a, const This &x)
  {
    for (unsigned int i=0;i<size();++i)
      (*this)[i] += a*x[i];
  }
  // this -> x + a * this 
  void xpay(double a, const This &x)
  {
    for (unsigned int i=0;i<size();++i)
      (*this)[i] = x[i] + a*(*this)[i];
  }
  // compute the squared l2 norm of this vector
  double norm2() const
  {
    return (*this)*(*this);
  }

  void clear()
  {
    for( unsigned int i = 0; i < size(); ++i )
      (*this)[ i ] = 0;
  }
};

template< unsigned int mydim >
struct RangeVector : public std::array< double, mydim >
{
  typedef std::array< double, mydim > Base;
  typedef RangeVector< mydim > This;
  static const int dim = mydim;

  explicit RangeVector() : Base() {}
  explicit RangeVector( const double a )
    : Base()
  {
    this->fill( a );
  }

#ifdef OpenCV_FOUND
  RangeVector( const cv::Point2d& pt )
    : Base()
  {
    if( dim != 2 )
      {
	std::stringstream ss;
	ss << "unable to convert cv::Point2d to RangeVector<" << dim << ">";
	throw ss.str();
      }
    this->at( 0 ) = pt.x;
    this->at( 1 ) = pt.y;
  }
#endif

#ifdef OpenCV_FOUND
  RangeVector( const cv::Point3d& pt )
    : Base()
  {
    if( dim != 3 )
      {
	std::stringstream ss;
	ss << "unable to convert cv::Point3d to RangeVector<" << dim << ">";
	throw ss.str();
      }
    this->at( 0 ) = pt.x;
    this->at( 1 ) = pt.y;
    this->at( 2 ) = pt.z;
  }
#endif

  This operator+( const This &other ) const
  {
    This out;
    for( unsigned int d = 0; d < dim; ++d )
      {
	out[ d ] = (*this)[ d ] + other[ d ];
      }
    return out;
  }

  This operator-( const This &other ) const
  {
    This out;
    for( unsigned int d = 0; d < dim; ++d )
      {
	out[ d ] = (*this)[ d ] - other[ d ];
      }
    return out;
  }

  This operator*( const double l ) const
  {
    This out;
    for( unsigned int d = 0; d < dim; ++d )
      {
	out[ d ] = (*this)[ d ] * l;
      }
    return out;
  }

  void operator+=( const This &other )
  {
    for( unsigned int d = 0; d < dim; ++d )
      {
	(*this)[ d ] += other[ d ];
      }
  }

  void operator-=( const This &other )
  {
    for( unsigned int d = 0; d < dim; ++d )
      {
	(*this)[ d ] -= other[ d ];
      }
  }

  void operator*=( const double l )
  {
    for( unsigned int d = 0; d < dim; ++d )
      (*this)[ d ] *= l;
  }

  void operator/=( const double l )
  {
    for( unsigned int d = 0; d < dim; ++d )
      (*this)[ d ] /= l;
  }

  // compute the scalar product of this vector with a second one
  double operator*(const This &y) const
  {
    double ret = 0;
    for (unsigned int i=0;i<dim;++i)
      ret += (*this)[i]*y[i];
    return ret;
  }

  double norm2() const
  {
    return (*this)*(*this);
  }

  double norm() const
  {
    return sqrt( norm2() );
  }

  double dist( const This& other ) const
  {
    double sum = 0;
    for (unsigned int d = 0; d < dim; ++d )
      {
	sum += ( (*this)[d] - other[d] ) * ( (*this)[d] - other[d] );
      }
    return std::sqrt( sum );
  }

  void clear()
  {
    this->fill( 0.0 );
  }

  This perp() const
  {
    if( dim != 2 )
      {
	static bool warn = true;
	if( warn )
	  {
	    std::cerr << __FILE__ << ":" << __LINE__ << ":WARNING: perp() not implemented for dim != 2" << std::endl;
	    warn = false;
	  }

	This ret;
	ret.clear();
	return ret;
      }

    This ret;
    ret[0] = (*this)[1];
    ret[1] = -(*this)[0];

    return ret;
  }

#ifdef OpenCV_FOUND
  // implicit conversions to opencv point types
  operator cv::Point2d() const
  {
    if( dim != 2 )
      {
	std::stringstream ss;
	ss << "unable to convert RangeVector<" << dim << "> to cv::Point2d";
	throw ss.str();
      }

    cv::Point2d pt;
    pt.x = this->at( 0 );
    pt.y = this->at( 1 );
    return pt;
  }

  operator cv::Point3d() const
  {
    if( dim != 3 )
      {
	std::stringstream ss;
	ss << "unable to convert RangeVector<" << dim << "> to cv::Point3d";
	throw ss.str();
      }

    cv::Point3d pt;
    pt.x = this->at( 0 );
    pt.y = this->at( 1 );
    pt.z = this->at( 2 );
    return pt;
  }
#endif

  template< const unsigned int thedim >
  friend std::ostream& operator<<(std::ostream& os, const RangeVector<thedim>& v);
};

template< const unsigned int dim >
std::ostream& operator<<(std::ostream& os, const RangeVector<dim>& v)
{
  os << "(";
  for( const auto vv : v )
    {
      os << " " << vv;
    }
  os << " )";
  return os;
}

template< class T >
struct SubArray
{
  typedef typename T :: value_type value_type;
  typedef typename T :: iterator iterator;
  typedef typename T :: const_iterator const_iterator;

  SubArray( T& array )
    : begin_( array.begin() ), end_( array.end() )
  {}

  SubArray( const iterator& begin, const iterator& end )
    : begin_( begin ), end_( end )
  {}

  SubArray( const SubArray<T>& other )
    : begin_( other.begin_ ), end_( other.end_ )
  {}

  value_type& operator[]( const unsigned int i )
  {
    assert( i < size() );
    iterator it = begin() + i;
    value_type& v = *it;
    return v;
  }

  value_type operator[]( const unsigned int i ) const
  {
    assert( i < size() );
    const_iterator it = begin() + i;
    const value_type& v = *it;
    return v;
  }

  iterator begin()
  {
    return begin_;
  }

  iterator end()
  {
    return end_;
  }

  const_iterator begin() const
  {
    return begin_;
  }

  const_iterator end() const
  {
    return end_;
  }

  void operator/=( const double a )
  {
    for( auto& d : *this )
      d /= a;
  }

  unsigned int N() const
  {
    return std::distance( begin_, end_ );
  }

  unsigned int size() const
  {
    return N();
  }

private:
  const iterator begin_;
  const iterator end_;
};

#endif
