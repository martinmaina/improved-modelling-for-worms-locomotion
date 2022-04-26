#ifndef MODEL_HH
#define MODEL_HH

#include <memory>
#include <iostream>      
#include <functional>
#include "vector.hh"
#include "timeprovider.hh"
#include "parameter.hh"
#include "matrix.hh"
#include "discretefunction.hh"
#include<vector>
#include<math.h>
#ifdef OpenCV_FOUND
#include "distancefunction.hh"
#endif


template< unsigned int mydim >
struct ModelInterface
{
  static const unsigned int dim = mydim;
  using RangeVectorType = RangeVector< dim >;

  ModelInterface( const TimeProvider& timeProvider ) // custom constructor for struct "ModelInterface"
    : timeProvider_( timeProvider )
  {}

  virtual RangeVectorType X0( const double s ) const = 0;
  virtual RangeVectorType f( const RangeVectorType& x ) const = 0;
  virtual double gamma( const double s ) const = 0;

  virtual double K( const double s ) const = 0;
  virtual double e( const double s ) const = 0;

  virtual double beta( const double s ) const = 0;
  enum Constraint { velocity, measure };
  virtual Constraint constraint() const = 0;

  const TimeProvider& timeProvider() const
  {
    return timeProvider_;
  }

private:
  const TimeProvider& timeProvider_;
};

template< const unsigned int mydim >
struct ModelDefault
  : public ModelInterface< mydim >
{
private:
  using BaseType = ModelInterface< mydim >;

public:
  static const unsigned int dim = BaseType :: dim;
  using RangeVectorType = typename BaseType :: RangeVectorType;
  using Constraint = typename BaseType :: Constraint;

  ModelDefault( ConfigParameters& parameters, const TimeProvider& timeProvider ) //constructor
    : BaseType( timeProvider ),
      gamma_( parameters.get<double>( "model.gamma" ) ),
      K_( parameters.get<double>( "model.K" ) ),
      e_( parameters.get<double>( "model.e" ) )
  {}

  virtual RangeVectorType X0( const double s ) const
  {
    RangeVectorType ret;
    ret[ 0 ] = s;
    for( unsigned int d = 1; d < dim; ++d )
      {
	ret[ d ] = 0.0;
      }

    return ret;
  }

  virtual RangeVectorType f( const RangeVectorType& /*x*/ ) const
  {
    RangeVectorType ret(0);
    return ret;
  }

  virtual double gamma( const double /*s*/ ) const
  {
    return gamma_;
  }

  virtual double K( const double /*s*/ ) const
  {
    return K_;
  }

  virtual double e( const double /*s*/ ) const
  {
    return e_;
  }

  virtual double energyDensity( const RangeVectorType& X, const RangeVectorType& Y ) const
  {
    return 0;
  }

  virtual double beta( const double s ) const
  {
    return 0;
  }

  virtual Constraint constraint() const
  {
    return Constraint::measure;
  }

  template< class DF >
  void extraForcing( const DF& X, DF& rhs ) const
  {}

  template< class DF >
  double extraEnergy( const DF& /*X*/ ) const
  {
    return 0;
  }

  using BaseType :: timeProvider;

private:
  const double gamma_;
  const double K_;
  const double e_;
};

#ifdef OpenCV_FOUND
template< const unsigned int mydim, const unsigned int nc >
struct DistanceModel
  : public ModelDefault< mydim >
{
private:
  using BaseType = ModelDefault< mydim >;

public:
  static const unsigned int dim = BaseType :: dim;
  static const unsigned int nCam = nc;
  using RangeVectorType = typename BaseType :: RangeVectorType;

  using DistanceFunctionType = DistanceFunction<dim>;

  DistanceModel( ConfigParameters const& parameters, DistanceFunctionType& distanceFunction,
		 const TimeProvider& timeProvider, const TimeProvider& imageTimeProvider )
    : BaseType( const_cast<ConfigParameters&>(parameters) // BUG
			 , timeProvider ),
      distance_( distanceFunction ),
      imageTimeProvider_( imageTimeProvider ),
      firstStepMaxPenality_( parameters.get<double>( "model.e.firststep.maxpenality") ),
      firstStepRelaxationTime_( parameters.get<double>( "model.e.firststep.relaxationtime") )
  {}

  virtual RangeVectorType X0( const double s ) const
  {
    RangeVectorType ret = distance_.middle();
    ret[ 1 ] += gamma(s)*(s-0.5);

    return ret;
  }

  virtual RangeVectorType f( const RangeVectorType& x ) const
  {
    // distance forcing
    RangeVectorType ret = distance_.normalDistance( x );
    ret *= -1.0;

    return ret;
  }

  double energyDensity( const RangeVectorType& X, const RangeVectorType& Y ) const
  {
    // fidelity term
    const double fid = 0.5 * distance_.distance( X );

    // curvature term
    const double curvature = 0.5 * this->e(0) * Y.norm();

    return fid + curvature;
  }

  double energyDensity( const RangeVectorType& X, const RangeVectorType& Y,
			std::vector<double>& list ) const
  {
    // fidelity term
    const double fid = 0.5 * distance_.distance( X, list );

    // curvature term
    const double curvature = 0.5 * BaseType::e(0) * Y.norm();
//   std::cout << "curv: " << curvature << std::endl;

    // separate into list
    list.push_back( curvature );

    return fid + curvature;
  }

  template< class DF >
  double extraEnergy( const DF& X ) const
  {
    return distance_.extraEnergy( X );
  }

  template< class DF >
  void extraForcing( const DF& X, DF& rhs ) const
  {
    distance_.extraForcing( X, rhs );
  }

  virtual double beta( const double s ) const
  {
    return 0;
  }

  virtual double e( const double s ) const
  {
    const double olde = BaseType::e(s);

    if( imageTimeProvider_.iteration() > 0 )
      return olde;

    // rescale first step regularisation
    const double time = this->timeProvider().time();
    return olde * ( 1.0 + firstStepMaxPenality_ * exp( - time / firstStepRelaxationTime_ ) );
  }

  using BaseType::gamma;

private:
  DistanceFunctionType& distance_;
  const TimeProvider& imageTimeProvider_;

  const double firstStepMaxPenality_;
  const double firstStepRelaxationTime_;
};
#endif // #ifndef OpenCV_FOUND

template< const unsigned int mydim >
struct LambdaMuscleModel
  : public ModelDefault< mydim >
{
private:
  using BaseType = ModelDefault< mydim >;

public:
  static const unsigned int dim = BaseType :: dim;
  using RangeVectorType = typename BaseType :: RangeVectorType;
  using Constraint = typename BaseType :: Constraint;

  LambdaMuscleModel( ConfigParameters& parameters,
		     const std::function< RangeVectorType( double ) >& X0,
		     const std::function< double( double, double ) >& beta,
		     const TimeProvider& timeProvider )
    : BaseType( parameters, timeProvider ),
      X0_( X0 ),
      beta_( beta ),
      constraint_( parameters.getFromList( "model.constraint", { "velocity", "measure" } ) ),
      eps_( parameters.get<double>( "model.I2.eps" ) )
  {}

  virtual RangeVectorType X0( const double s ) const
  {
    return X0_( s );
  }

  using BaseType::f;
  using BaseType::gamma;
  using BaseType::K;

  virtual double e( const double s ) const
  {
    const double eps = eps_;
    const double num = 2.0 * std::sqrt( (eps + s)*(eps + 1.0 - s) );
    const double den = 1.0 + 2.0 * eps;
std::cout << "hellooooo" << std::endl;
    return BaseType::e(s) * 10.0; // ( num*num*num ) / ( den*den*den ); // constant radius?? <<<<<<<<<<<<<<<<

  }


  virtual double beta( const double s ) const
  {
    const double t = timeProvider().time();
    return beta_( s, t );
  }

  virtual double energyDensity( const RangeVectorType& X, const RangeVectorType& Y ) const
  {
    return 0;
  }

  virtual Constraint constraint() const
  {
    switch( constraint_ ) {
    case 0:
      return Constraint::velocity;
    case 1:
      return Constraint::measure;
    default:
      return Constraint::measure;
    }

  }

  using BaseType::timeProvider;

private:
  const std::function< RangeVectorType( double ) >& X0_;
  const std::function< double( double, double ) >& beta_;
  const int constraint_;
  const double eps_;
};

template< const unsigned int mydim >
struct DrivenMuscleModel
  : public ModelDefault< mydim >
{
private:
  using BaseType = ModelDefault< mydim >;

public:
  static const unsigned int dim = BaseType :: dim;
  using RangeVectorType = typename BaseType :: RangeVectorType;


  DrivenMuscleModel(const Mesh& mesh, ConfigParameters& parameters, const TimeProvider& timeProvider )
    : BaseType( parameters, timeProvider ),
      beta0_( parameters.get<double>( "model.muscle.beta0", M_PI ) ),
      lambda_( parameters.get<double>( "model.muscle.lambda", 2.0 /1.5 ) ),
      omega_( parameters.get<double>( "model.muscle.omega", 2.0 * M_PI ) ),
      curvatureData_( mesh.N()*dim ),
      curvature_( mesh, SubArray< Vector >( curvatureData_.begin(), curvatureData_.end() ) )
	
  {}

  virtual RangeVectorType X0( const double s ) const
  {
    RangeVectorType ret;
    ret[ 0 ] = s;
    for( unsigned int d = 1; d < dim; ++d )
      {
	ret[ d ] = 0.0;
      }

    return ret;
  }

  using BaseType::f;
  using BaseType::gamma;
  using BaseType::K;
  using BaseType::e;

  using CurvatureType = PiecewiseLinearFunction<dim>;
  Vector curvatureData_;
  CurvatureType curvature_;
  void storeCurvature( const CurvatureType& curvature )
  {
    curvature_.assign( curvature );
  }

  virtual double beta( const double s ) const
  {
/*
    // forcing is:    beta = beta0 sin( 2 pi s / lambda - omega t )
    const double t = timeProvider().time();
    const double q = 2.0 * M_PI / lambda_ * s - omega_ * t;

    return beta0_ * sin( q );
*/
    const double t = timeProvider().time();
    const double q = 2.0 * M_PI * s/0.3 - 2.0* M_PI *0.5*3.3*t;

    return 10.0 * sin( q );

  }

using BaseType::timeProvider;

private:
  const double beta0_, lambda_, omega_;
};





template< const unsigned int mydim >
struct ReadMuscleModel
  : public ModelDefault< mydim >
{
private:
  using BaseType = ModelDefault< mydim >;

public:
  static const unsigned int dim = BaseType :: dim;
  using RangeVectorType = typename BaseType :: RangeVectorType;

  ReadMuscleModel( const Mesh& mesh, ConfigParameters& parameters, const TimeProvider& timeProvider )
    : BaseType( parameters, timeProvider ),
      curvatureData_( mesh.N()*dim ),
      curvature_( mesh, SubArray< Vector >( curvatureData_.begin(), curvatureData_.end() ) ),
      positionData_( mesh.N()*dim ),
      position_( mesh, SubArray< Vector >(positionData_.begin(), positionData_.end() ) ),
      N(parameters.get<double>("mesh.N")),
      beta0_( parameters.get<double>( "model.muscle.beta0", M_PI ) ),
      lambda_( parameters.get<double>( "model.muscle.lambda", 2.0 /1.5 ) ),
      omega_( parameters.get<double>( "model.muscle.omega", 2.0 * M_PI ) ),
      eps_( parameters.get<double>( "model.I2.eps" ) ), 
      timeToUpdate(0.0),
      freq(parameters.get<double>("model.beta.freq")*3.3),
      amp(parameters.get<double>("model.beta.amp")),
      wave(parameters.get<double>("model.beta.wave")),
      propFrom(parameters.get<double>("model.read.propFrom")),
      propTo(parameters.get<double>("model.read.propTo")),
      propStrength(parameters.get<double>("model.read.propStrength")*3.3),
      a_(parameters.get<double>("model.torque.a")*3.3),
      b_(parameters.get<double>("model.torque.b")*3.3),
      tau(N),
      dtau(N),
      scalarCurvature(N),
      propFeedback(N),
      //nu_out(NULL),
      Activation(NULL),
      Beta(NULL),
      Kappa(NULL),
      pointx(NULL),
      pointy(NULL)
 {
       /* nu_out = new std::ofstream("nu_out.txt");
	if( !nu_out->is_open() )
        {
          std::cerr << "Error: nu_out.txt didn't open" << std::endl;
        }*/
        Activation = new std::ofstream("Activation.txt");
	if( !Activation->is_open() )
        {
          std::cerr << "Error: Activation.txt didn't open" << std::endl;
        }
        Beta = new std::ofstream("Beta.txt");
	if( !Beta->is_open() )
        {
          std::cerr << "Error: Beta didn't open" << std::endl;
        }
	Kappa = new std::ofstream("Kappa.txt");
        if( !Kappa->is_open() )
        {
          std::cerr << "Error: Kappa didn't open" << std::endl;
        }
	pointx = new std::ofstream("pointx");
        if( !pointx->is_open() )
        {
          std::cerr << "Error: pointx file didn't open" << std::endl;
        }
	pointy = new std::ofstream("pointy");
        if( !pointy->is_open() )
        {
          std::cerr << "Error: pointy file didn't open" << std::endl;
        }
 }


~ReadMuscleModel() // destructor causes core dump for some reason :s 
{
/*  nu_out->close();
  if( nu_out )
    delete nu_out;*/
  Activation->close();
  //if( Activation )
  //  delete Activation;
  Beta->close();
  //if( Beta )
  //  delete Beta;
  Kappa->close();
  //if( Kappa )
  //  delete Kappa;
  pointx->close();
  pointy->close();
}

  virtual double e( const double s ) const
  {
    const double eps = eps_;
    const double num = 2.0 * std::sqrt( (eps + s)*(eps + 1.0 - s) );
    const double den = 1.0 + 2.0 * eps;

    return BaseType::e(s) * ( num*num*num ) / ( den*den*den ); // keep in mind Ma = -E*I2*beta still <<<<<<<<<<<<<<<<

  }

  virtual double rad( const double s ) const
  {
    const double eps = 1e-2;
    const double num = 2.0 * std::sqrt( (eps + s)*(eps + 1.0 - s) );
    const double den = 1.0 + 2.0 * eps;

    return ( num*num*num ) / ( den*den*den ); // keep in mind Ma = -E*I2*beta still <<<<<<<<<<<<<<<<

  }

  virtual RangeVectorType X0( const double s ) const
  {
    RangeVectorType ret;
	ret[0] = s; 		//cos(-M_PI*s)/M_PI;
	ret[1] = 0.0; //0.0; //s*(1.-s)		//sin(M_PI*s)/M_PI;

    for( unsigned int d = 2; d < dim; ++d )
      {
	ret[ d ] = 0.0;
      }

    return ret;
  }

  using BaseType::f;
  using BaseType::gamma;
  using BaseType::K;
  using BaseType::e;

  using CurvatureType = PiecewiseLinearFunction<dim>;
  Vector curvatureData_;
  CurvatureType curvature_;

  void storeCurvature( const CurvatureType& curvature )
  {
    curvature_.assign( curvature );
  }

  using PositionType = PiecewiseLinearFunction<dim>;
  Vector positionData_;
  PositionType position_;

  void storePosition( const PositionType& position )
  {
    position_.assign( position );
  }  

  double boundary(const double s) const
  {

	double out;
	if(std::abs(s) < 1e-7) // 1e-7
	{
		out = 1.0;
	}
	else if(std::abs(s) > 0.5 ) // 1e-7
	{
		out = 0.0;
	}
	else
	{
		out = 1.0;
		//std::cout << out << std::endl;
	}

	return out; //rad(s);

  }


mutable double N;

  void write_Activation()
  {

	*Activation << " " <<  timeProvider().time();

	for(int j=0; j < N; j++) // technically j < M, but for now (continuous torque), M=N effectively
	{
		double pos = double(j)/double(N);
		*Activation << " " << a_*activation(pos); //
	}
	
	*Activation << "\n";

  }

  void write_beta()
  {

	*Beta << " " <<  timeProvider().time();

	for(int j=0; j < N; j++)
	{
		double pos = double(j)/double(N);
		*Beta << " " << b_*beta(pos); //
	}
	
	*Beta << "\n";

  }


  void write_kappa()
  {

        *Kappa << " " <<  timeProvider().time();

        for(int j=0; j < N; j++)
        {
	double pos = double(j)/double(N);
		if(j < 1e-8)
		{
			*Kappa << " " << -beta(pos);
		}
		else if(abs(j-(N-1)) < 1e-8 )
		{
			*Kappa << " " << -beta(pos);
		}
		else
		{
	                *Kappa << " " << scalarCurvature[j];
	        } 
	}

        *Kappa << "\n"; 
  }



  void write_pointx()
  {

        *pointx << " " <<  timeProvider().time();

        for(int j=0; j < N; j++)
        {
                *pointx << " " << position_.evaluate(j)[0] ;
        }

        *pointx << "\n";
  }

  void write_pointy()
  {

        *pointy << " " <<  timeProvider().time();

        for(int j=0; j < N; j++)
        {
                *pointy << " " << position_.evaluate(j)[1] ;
        }

        *pointy << "\n";
  }


/*  void write_nu()
  {

	for(int j=0; j < N-1; j++)
	{
	const auto xe  = position_.evaluate( j  );
	const auto xep = position_.evaluate( j+1 );
	RangeVectorType nu = (xe - xep).perp();
	const auto a = nu.norm();
	nu[0] = nu[0]/a;
	nu[1] = nu[1]/a;
	*nu_out << " " << nu[0] << " " << nu[1];
	}
  *nu_out << "\n";
  }*/

 void updateScalarCurvature() const // BUG: acts on internal meshpoints only!! (j != 0,N-1)
  {
	assert(scalarCurvature.size() == N);
	scalarCurvature[0] = -beta(0); // to fit BCs 
	scalarCurvature[N-1] = -beta(1); // possible bug: beta same timestep??

	for(unsigned int j=1; j < N-1; j++ )
	{

		const auto xem = position_.evaluate( j-1 );
		const auto xe  = position_.evaluate( j  );
		const auto xep = position_.evaluate( j+1 );
	
		RangeVectorType num = (xem - xe).perp();
		RangeVectorType nup = (xe - xep).perp();
		const auto a = num.norm();
		const auto b = nup.norm();
		num[0] = num[0]/a;
		num[1] = num[1]/a;
		nup[0] = nup[0]/b;
		nup[1] = nup[1]/b;

		RangeVectorType nubar = num + nup;
		const auto c = nubar.norm();
		nubar[0] = nubar[0]/c;
		nubar[1] = nubar[1]/c;
		scalarCurvature[j] = curvature_.evaluate(j)[0]*nubar[0] + curvature_.evaluate(j)[1]*nubar[1];
	}
  }

  void updatePropFeedback() const
   {
/*	double sum;
	int S;
	double a;
	double b;
	int ia;
	int ib;

	for(double s = 0.0 ; s <= 1.0 ; s = s + 1/N)
	{
		
		S = round(s*(N-1));
		a = std::max(s + propFrom, 0.0);
		a = std::min(a, 1.0);
		b = std::min(s + propTo,   1.0 + 1.0/N); // so ia-ib != 0
		ia = a*N;
		ib = b*N;
		sum = 0.0;

		for(int i = ia ; i < ib ; ++i)
		{
			if( i < 1 || i > N-1 )
			{
				sum = sum +  0;
			}
			else
			{
				sum = sum + scalarCurvature[i];
			} 
		}
		propFeedback[S] = sum/double(ib - ia);
	}*/
   }


void printScalarCurvature() const
  {

	std::cout << timeProvider().time() << " , " <<  scalarCurvature.size()  << " curve: ";
	for(unsigned int k = 0; k < N ; k++  )
	{
	 std::cout << scalarCurvature[ k ] << " ,  " ;	
	}
	std::cout << std::endl;
  }

void printPropFeedback() const
  {

	std::cout << timeProvider().time() << " , (tau[k]) ";
	for(unsigned int k = 0; k < N ; k++  )
	{
	 std::cout << propFeedback[ k ] << " ,  " << tau[ k ];	
	}
	std::cout << std::endl;
  }

  double activation( const double s ) const
   {

	// Output travelling sine wave
	const double t = timeProvider().time();
	const double q = 2.0 * M_PI * s/wave - 2.0* M_PI *freq*t;
	double output = amp * sin( q ); // boundary(s)*

	return output;
   }


  virtual double beta( const double s ) const
  {

	int S = round(s*(N-1));

	dtau[S] = a_*activation(s) - b_*tau[S] ;// + propStrength * propFeedback[S]; // +ve curvature feedback (due to weird -ve beta to kappa relationship - it's actually inhibitory feedback)

	if(timeToUpdate <= timeProvider().time())
	{
		for(int k=0; k<N; k++)
		{
			tau[k] = tau[k] + timeProvider().deltaT()*dtau[k]; 
		}
	timeToUpdate = timeToUpdate + timeProvider().deltaT();
	}

	return activation(s); //tau[S];   //tau[S]; //activation(s);
  }

using BaseType::timeProvider;

private:
  const double beta0_, lambda_, omega_;
  const double eps_;
  mutable double timeToUpdate;
  const double freq, amp, wave;
  const double propFrom, propTo, propStrength;
  const double a_, b_;
  mutable std::vector<double> tau;
  mutable std::vector<double> dtau;
  mutable std::vector<double> scalarCurvature;
  mutable std::vector<double> propFeedback;
	//std::ofstream* nu_out;
	std::ofstream* Activation;
	std::ofstream* Beta;
	std::ofstream* Kappa;
	std::ofstream* pointx;
	std::ofstream* pointy;
}; // ReadMuscleModel

#endif // #ifndef MODEL_HH
