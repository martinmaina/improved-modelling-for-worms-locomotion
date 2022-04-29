#include "config.hh"

#include <iostream>
#include <vector>  // added by jack
#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <exception>
#include <limits>

// input parameters
#include "lib/parameter.hh"

// mesh and time helpers
#include "lib/mesh.hh"
#include "lib/timeprovider.hh"

// problem definitions
#include "lib/model.hh"
#include "lib/modelThreshold.hh"
#include "lib/modelOneThreshold.hh"
//#include "lib/model_neural.hh"  // think about doing this
#include "lib/femscheme.hh"

// output helpers
#include "lib/vtuwriter.hh"

#include "lib/vector.hh"
#include "lib/discretefunction.hh"


void algorithm( ConfigParameters& parameters, const Mesh &mesh)
{

  static const unsigned int dim = WORLD_DIM;

  // time holder
  TimeProvider timeProvider;
  const double endTime = parameters.get< double >( "time.end" )/3.3;
  double deltaT = parameters.get< double >( "time.deltaT" );
  double transient = parameters.get<double>("io.transient")/3.3;
  timeProvider.setDeltaT( deltaT );

  // find model 
  using ModelType = ThresholdModel< dim > ; //ThresholdModel<dim> ;// ReadMuscleModel< dim >;

  ModelType model( mesh, parameters, timeProvider );

  // initialise scheme
  using FemSchemeType = FemScheme< ModelType, WORLD_DIM >;
  FemSchemeType scheme( mesh, model );
  scheme.initialiseXY();

  // io time helpers
  double nextWriteTime = 0;
  const double nextWriteTimeStep = parameters.get< double >( "io.timestep" );

  VTUWriter writer( parameters, timeProvider );

  //model.printSwitch();
  model.storePosition( scheme.position() ); // Jack  
  scheme.computeCurvature(); //Jack
  model.storeCurvature( scheme.curvature() ); // Jack
  model.updateScalarCurvature();
  model.updatePropFeedback();

  // output important data
  if( timeProvider.time() >= nextWriteTime ) 
  {
    //model.storePosition( scheme.position() );
    //scheme.computeCurvature();
    writer.writeVtu( scheme.position(), scheme.curvature(), scheme.pressure() );
    nextWriteTime += nextWriteTimeStep;
  }

  // update time
  timeProvider.next( deltaT );
  std::ofstream file;
  file.open( parameters.getString( "io.prefix" ) + "/energy.txt" );
  file.precision( std::numeric_limits<double>::max_digits10 );

  for( ; timeProvider.time() < endTime; timeProvider.next( deltaT ) )
  {

    // prepare
    scheme.prepare();

    // solve
    scheme.solve();

    const double e = scheme.energy();
    const auto p = scheme.constraintError();
    file << timeProvider.time() << " " <<  e << " "
   << p.first << " " << p.second << std::endl;

   //model.printSwitch();
   scheme.computeCurvature(); 			//Jack
   model.storeCurvature( scheme.curvature() );  //Jack
   model.storePosition( scheme.position() );    //Jack
   model.updateScalarCurvature(); 		//Jack
   model.updatePropFeedback();
   model.threshold_prop();

//std::cout << "============== ONE TIME LOOP COMPLETED ================" << std::endl;

    // output important data
    if( timeProvider.time() - transient>= nextWriteTime )
    {

//      const double e = scheme.energy();
//      const auto p = scheme.constraintError();
//      file << timeProvider.time() << "," <<  e << ","
//	   << p.first << "," << p.second << std::endl;

      // output important data
//      if( timeProvider.time() >= nextWriteTime ) 
//	{

          // output beta (Jack)
	  model.write_Activation();
	  model.write_beta();
	  model.write_kappa();
	  model.write_feedback();
	  model.write_pointx();
	  model.write_pointy();
	
	  //model.printSwitch();
	  //model.printTorque();

	  //scheme.computeCurvature();
	  writer.writeVtu( scheme.position(), scheme.curvature(), scheme.pressure() );
	  nextWriteTime += nextWriteTimeStep;


//	}
    }
  }

  file.close();
}

int main( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "usage: " << argv[0] << " {parameter file}" << std::endl;
    return 1;
  }

  try
  {
   
     ConfigParameters parameters( argv[1] );
    for( int i = 2; i < argc; ++i )
    {
      parameters.add( argv[ i ] );
    }
  
    std::cout << "compile time variables:" << std::endl;
    std::cout << " - WORLD_DIM: " << WORLD_DIM << std::endl;

    // make mesh
    const unsigned int N = parameters.get< unsigned int >( "mesh.N" );
    const double a = parameters.get< double >( "mesh.a" );
    const double b = parameters.get< double >( "mesh.b" );
    Mesh mesh( a, b, N );

    // run algorithm
    algorithm( parameters, mesh);

    return 0;
  }
  catch( const std::string& e )
  {
    std::cout << "Exception: " << e << std::endl;
    return 1;
  }
  catch( const char* e )
  {
    std::cout << "Exception: " << e << std::endl;
    return 1;
  }
  catch( std::exception& e )
  {
    std::cerr << "exception caught: " << e.what() << std::endl;
    return 1;
  }
}
