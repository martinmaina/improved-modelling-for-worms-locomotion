#ifndef VTUWRITER_HH
#define VTUWRITER_HH

#include <iomanip>
#include <sstream>
#include <limits>

#include "parameter.hh"
#include "mesh.hh"
#include "timeprovider.hh"

struct VTUWriter
{
  VTUWriter( ConfigParameters const& parameters, const TimeProvider& timeProvider,
	     const std::string fnPrefix = "" )
    : parameters_( parameters ), timeProvider_( timeProvider ), fnPrefix_( fnPrefix )
  {}

  VTUWriter( const VTUWriter& other )
    : parameters_( other.parameters_), timeProvider_( other.timeProvider_ ), fnPrefix_( other.fnPrefix_ )
  {}

  ~VTUWriter()
  {
    if( pvd_ )
      closePvd();
  }

  template< class DV1, class DV2, class DV3 >
  void writeVtu( const DV1 &X, const DV2 &Y, const DV3& P )
  {
    std::stringstream ss_head;
    ss_head << fnPrefix_ << "worm_" << std::setfill('0') << std::setw(6) <<  timeProvider_.iteration() << ".vtu";
    writeVtu( X, Y, P, ss_head.str() );
  }

  template< class DV1, class DV2, class DV3 >
  void writeVtu( const DV1 &X, const DV2 &Y, const DV3& P, const std::string& head )
  {
    using RangeVectorType = typename DV1::RangeVectorType;

    const unsigned int N = X.N();

    std::stringstream ss;
    const std::string outputDir = parameters_.get< std::string >( "io.prefix" );
    ss << outputDir << "/" << head;

    std::ofstream file;
    file.open( ss.str().c_str() );
    if( not file.is_open() )
      throw "unable to open file " + ss.str();

    // set to full precision
    file.precision( std::numeric_limits<double>::max_digits10 );

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << N << "\" NumberOfCells=\"" << N-1 << "\">\n";
    file << "      <PointData Scalars=\"W\">\n";
    file << "        <DataArray type=\"Float32\" Name=\"W\" NumberOfComponents=\"" << Y.dim << "\">\n";
    file << "         ";
    for( const double Yj : Y )
      {
	file << " " << Yj;
      }
    file << "\n";
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float32\" Name=\"u\" NumberOfComponents=\"" << 1 << "\">\n";
    file << "         ";
    for( unsigned int j = 0; j < X.mesh().N(); ++j )
      {
	file << " " << X.mesh().u(j);
      }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </PointData>\n";

    file << "      <CellData Scalars=\"P\">\n";
    file << "        <DataArray type=\"Float32\" Name=\"P\" NumberOfComponents=\"" << P.dim << "\">\n";
    file << "         ";
    for( const double Pj : P )
      {
	file << " " << Pj;
      }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for( unsigned int j = 0 ; j < N; ++j )
      {
	file << "         ";
	RangeVectorType Xj = X.evaluate( j );
	for( unsigned int d = 0; d < X.dim; ++d )
	  {
	    file <<  " " << Xj[ d ];
	  }
	for( unsigned int d = X.dim; d < 3; ++d )
	  {
	    file << " 0";
	  }
	file << "\n";
      }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    file << "          ";
    for( unsigned int j = 0; j < N-1; ++j )
      {
	file << j << " " << j+1 << " ";
      }
    file << "\n";
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    file << "          ";
    for( unsigned int j = 0; j < N-1; ++j )
      {
	file << 2*(j+1) << " ";
      }
    file << "\n";
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    file << "          ";
    for( unsigned int j = 0; j < N-1; ++j )
      {
	file << 3 << " ";
      }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();

    if( not pvd_.is_open() )
      openPvd();

    pvd_ << "    <DataSet timestep=\"" << timeProvider_.time() << "\" " << "group=\"\" part=\"0\"\n"
	 << "             file=\"" << head <<"\"/>" << std::endl;

    std::cout << ss.str() << " written" << std::endl;
  }

protected:
  void openPvd()
  {
    if( pvd_.is_open() )
      return;

    const std::string outputDir = parameters_.get< std::string >( "io.prefix" );
    const std::string pvdFilename = outputDir + "/" + fnPrefix_ + "worm.pvd";
    pvd_.open( pvdFilename.c_str() );

    if( not pvd_.is_open() )
      {
        throw "unable to open " + pvdFilename;
      }
    else
      {
	pvd_ << "<?xml version=\"1.0\"?>" << std::endl;
	pvd_ << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
	pvd_ << "  <Collection>" << std::endl;
      }
  }

  void closePvd()
  {
    if( pvd_.is_open() )
      {
	pvd_ << "  </Collection>" << std::endl;
	pvd_ << "</VTKFile>" << std::endl;

	pvd_.close();
	std::cout << "data collated in " << fnPrefix_ << "worm.pvd" << std::endl;
      }
  }

private:
  ConfigParameters parameters_;
  std::ofstream pvd_;
  const TimeProvider& timeProvider_;
  const std::string fnPrefix_;
};

#endif
