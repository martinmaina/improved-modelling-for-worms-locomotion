#ifndef TIMEPROVIDER_HH
#define TIMEPROVIDER_HH

#include "parameter.hh"

struct TimeProvider
{
  TimeProvider( double startTime = 0.0 )
    : time_( startTime ),
      it_( 0 ),
      dt_( 0.0 ),
      dtValid_( false )
  {}

  TimeProvider( const TimeProvider& other )
    : time_( other.time_ ),
      it_( other.it_ ),
      dt_( other.dt_ ),
      dtValid_( other.dtValid_ )
  {}

  double time() const
  {
    return time_;
  }

  unsigned int iteration() const
  {
    return it_;
  }

  void next( const double dt )
  {
    setDeltaT( dt );
    time_ += deltaT();
    it_ +=1 ;
  }

  double deltaT() const
  {
    if( not dtValid_ )
      throw "invalid time step";
    return dt_;
  }

  void setDeltaT( const double dt )
  {
    dt_ = dt;
    dtValid_ = true;
  }


private:
  double time_;
  unsigned int it_;

  double dt_;
  bool dtValid_;
};

#endif // #ifndef TIMEPROVIDER_HH
