#ifndef MODELTHRESHOLD_HH
#define MODELTHRESHOLD_HH

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

template< const unsigned int mydim >
struct ThresholdModel
  : public ModelDefault< mydim >
{
private:
  using BaseType = ModelDefault< mydim >;

public:
  static const unsigned int dim = BaseType :: dim;
  using RangeVectorType = typename BaseType :: RangeVectorType;

  ThresholdModel( const Mesh& mesh, ConfigParameters& parameters, const TimeProvider& timeProvider )
    : BaseType( parameters, timeProvider ),
      curvatureData_( mesh.N()*dim ),
      curvature_( mesh, SubArray< Vector >( curvatureData_.begin(), curvatureData_.end() ) ),
      positionData_( mesh.N()*dim ),
      position_( mesh, SubArray< Vector >(positionData_.begin(), positionData_.end() ) ),
      N(parameters.get<double>("mesh.N")),
      M(parameters.get<double>("model.beta.M")),
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
      a_(parameters.get<double>("model.torque.a")*3.3),
      b_(parameters.get<double>("model.torque.b")*3.3),
      tau(N),
      dtau(N),
      Switch(2, std::vector<double>(N)), // [0] Dorsal  [1] Ventral
      thd_0(parameters.get<double>("model.threshold.thd_0")),
      thd_1(parameters.get<double>("model.threshold.thd_1")),
      thv_00(parameters.get<double>("model.threshold.thv_00")),
      thv_01(parameters.get<double>("model.threshold.thv_01")),
      thv_10(parameters.get<double>("model.threshold.thv_10")),
      thv_11(parameters.get<double>("model.threshold.thv_11")),
      scalarCurvature(N),
      propFeedback(N),
      lastFB(N),
      //nu_out(NULL),
      reset(parameters.get<bool>("model.threshold.reset")),
      Activation(NULL),
      Beta(NULL),
      Kappa(NULL),
      feedback(NULL),
      pointx(NULL),
      pointy(NULL)
 {
	for(int j=0;j<N;j++)
	{
		Switch[0][j] = 0;
		Switch[1][j] = 0; //0
	}
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
	feedback = new std::ofstream("feedback.txt");
	if( !feedback->is_open() )
	{
  	  std::cerr << "Error: feedback file didn't open" << std::endl;
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


~ThresholdModel() // destructor causes core dump for some reason :s 
{
  Activation->close();
  //if( Activation )
  //  delete Activation;
  Beta->close();
  //if( Beta )
  //  delete Beta;
  Kappa->close();
  //if( Kappa )
  //  delete Kappa;
  feedback->close();
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
    const double eps = eps_;
    const double num = 2.0 * std::sqrt( (eps + s)*(eps + 1.0 - s) );
    const double den = 1.0 + 2.0 * eps;

    return ( num*num*num ) / ( den*den*den ); // keep in mind Ma = -E*I2*beta still <<<<<<<<<<<<<<<<

  }

  virtual RangeVectorType X0( const double s ) const
  {
    RangeVectorType ret;
	ret[0] = s; 		  //cos(-M_PI*s)/M_PI;
	ret[1] = 0.0; //0.0; //sin(M_PI*s)/M_PI; //s*(1.-s)

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

mutable double N;
mutable double M;

  void write_Activation()
  {

	*Activation << " " <<  timeProvider().time();

	for(int j=0; j < N; j++) // technically j < M, but for now (continuous torque), M=N effectively
	{
		double pos = double(j)/double(N);
		*Activation << " " << Switch[0][j] - Switch[1][j]; //
	}
	
	*Activation << "\n";

  }

  void write_beta()
  {

	*Beta << " " <<  timeProvider().time();

	for(int j=0; j < N; j++)
	{
		*Beta << " " << tau[j];
	}
	
	*Beta << "\n";

  }


  void write_kappa()
  {

        *Kappa << " " <<  timeProvider().time();

        for(int j=0; j < N; j++)
        {
        	*Kappa << " " << scalarCurvature[j];
	}

        *Kappa << "\n"; 
  }

  void write_feedback()
  {

        *feedback << " " <<  timeProvider().time();

        for(int j=0; j < N; j++)
        {
		double pos = double(j)/double(N);
		*feedback << " " << propFeedback[j] ; 
	}

        *feedback << "\n"; 
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



 void updateScalarCurvature() const // BUG: acts on internal meshpoints only!! (j != 0,N-1)
  {
	assert(scalarCurvature.size() == N);

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
	
	scalarCurvature[0]   = scalarCurvature[1];    // -beta(0)
	scalarCurvature[N-1] = scalarCurvature[N-2];  // -beta(1)
  }


  void updatePropFeedback() const
   {
	double sum;
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
		propFeedback[S] = sum/std::abs(double(ib - ia)); // *5/N
	}
   }


/*
// conserve prop range along all body -- no scaling required
  void updatePropFeedback() const
   {
        double sum;
        int S;
        double a;
        double b;
        int ia;
        int ib;
        for(double s = 0.0 ; s <= 1.0 ; s = s + 1/N)
        {
                S = round(s*(N-1));
		double propMin = 0.25 ;// propTo - propFrom; //  >> if propRange > propMin, range shrinks slightly towards tail, but not completely
                a = std::max(s + propFrom, 0.0);
                a = std::min(a, 1.0 - (propMin)); // prop range fixed towards tail
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
                        else if(i < S) // range becomes anterior to s
			{
				sum = sum - scalarCurvature[i];
			}
			else
                        {
                                sum = sum + scalarCurvature[i];
                        }
                }
                propFeedback[S] = sum/std::abs(double(N*(b - a))); // *N
        }
   }
*/
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

void printSwitch() const
  {
	std::cout << "time: " << timeProvider().time() << std::endl;
	std::cout << "Dorsal:   ";	
	for(int j=90; j<N; j++)
	{
		std::cout << Switch[0][j] << " ";
	}

	std::cout << std::endl;
	std::cout << "Ventral:  ";	
	for(int j=90; j<N; j++)
	{
		std::cout << Switch[1][j] << " ";
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
/*
  void threshold_FFD(const double s) const // When driving by noramlised sinewave
   {
	int S = round(s*(N-1));
	
	int thd_on  = .7 ; // +ve
	int thd_off = .2 ;
	int thv_on  = -.6 ; // -ve
	int thv_off = -.25 ;

	// Dorsal
	if( Switch[0][S] < 1e-10 && activation(s) >= thd_on )
		{
			Switch[0][S] = 1;
		}
	else if ( std::abs(Switch[0][S] - 1 ) < 1e-10 && activation(s) < thd_off )
		{
			Switch[0][S] = 0;
		}

	// Ventral
	if( Switch[1][S] < 1e-10 && activation(s) <= thv_on ) // thresh has opposite signs than dorsal, since curvature has flipped sign
		{
			Switch[1][S] = 1;
		}
	else if ( std::abs(Switch[1][S] - 1 ) < 1e-10 && activation(s) > thv_off )
		{
			Switch[1][S] = 0;
		}
	
   }
*/

  void threshold_prop() const// When proprioception drives locomotion
   {

       for(int j=0; j<N; ++j)
        {
        	if( propFeedback[j] < thd_0 ) /**/ { Switch[0][j] = 1;}
        	else if( propFeedback[j] > thd_1  ) /**/ { Switch[0][j] = 0;}
        	Switch[1][j] = 1-Switch[0][j];
	}

   }

  virtual double beta( const double s ) const
  {

	int S = round(s*(N-1));
//	std::cout  << S << " tau: " << -b_*tau[S] <<  " , curve: " << a_*scalarCurvature[ S ]  << " , prop: " << propStrength*propFeedback[ S ]  << std::endl
	dtau[S] =  -a_*(Switch[0][S] - Switch[1][S]) - b_*tau[S]; // +ve curvature feedback (due to weird -ve beta to kappa relationship - it's actually inhibitory feedback)

	if(timeToUpdate <= timeProvider().time())
	{
		for(int k=0; k<N; k++)
		{
			tau[k] = tau[k] + timeProvider().deltaT()*dtau[k]; 
		}
	timeToUpdate = timeToUpdate + timeProvider().deltaT();
	}

	return tau[S];   //tau[S]; //activation(s);
  }

using BaseType::timeProvider;

private:
  const double beta0_, lambda_, omega_;
  const double eps_;
  mutable double timeToUpdate;
  const double freq, amp, wave;
  const double propFrom, propTo;
  const double a_, b_;
  mutable std::vector<double> tau;
  mutable std::vector<double> dtau;
  mutable std::vector<std::vector<double>> Switch;
  const double thd_0, thd_1, thv_00, thv_01, thv_10, thv_11;
  mutable std::vector<double> scalarCurvature;
  mutable std::vector<double> propFeedback;
  mutable std::vector<double> lastFB;
  mutable bool reset;
	//std::ofstream* nu_out;
	mutable std::ofstream* Activation;
	mutable std::ofstream* Beta;
	mutable std::ofstream* Kappa;
	mutable std::ofstream* feedback;
	mutable std::ofstream* pointx;
	mutable std::ofstream* pointy;

}; // thresholdModel

#endif // #ifndef MODELTHRESHOLD_HH
