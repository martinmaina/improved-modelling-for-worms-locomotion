#ifndef PARAMETER_HH
#define PARAMETER_HH

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include "trace.hh"

/**
 *  \class ConfigParameters
 *
 *  \brief class for accessing parameters from files and command line
 *
 *  This class contains method to help with the use of parameters in
 *  the algorithms we want to use. The parameters are specified in a
 *  file with one key-value pair per line for example:
 \code{.unparsed}

  key1:value1
  key2:value2
  key3:value3

 \endcode
 *  Note that there are no spaces either side of the colon. In
 *  addition to a parameter file, extra key-value pairs can be added

 *  using the add method.
 */
class ConfigParameters
{
  typedef std::map< std::string, std::string > DataType;


public:
  /**
   * \brief constructor with parameter filename input
   */
  explicit ConfigParameters( const std::string& filename )
  {
    std::ifstream file;
    file.open( filename.c_str() );
    if( !file.is_open() )
      throw "invalid parameter file";

    std::cout << "reading parameters from " << filename << std::endl;
    readConfigFile( file );

    file.close();
  }

  /**
   * \brief set parameter with key
   */
  void set( const std::string& key, const std::string& value){
	  data_[key]=value;
  }

  /**
   * \brief get parameter with key
   */
  template< class T >
  T get( const std::string& key ) const
  {
    // get value from data
    const auto valueIt = data_.find( key );

    if( valueIt == data_.end() )
      {
	throw "unable to find key: " + key;
      }

    // extract data
    const std::string value = valueIt->second;

    // assign to string stream
    std::stringstream ss;
    ss << value;

    // extract from string stream
    T retValue;
    ss >> retValue;

    return retValue;
  }


  /**
   * \brief get parameter with key
   */
  template< class T >
  T get( const std::string& key, const T& def )
  {
    // get value from data
    const auto valueIt = data_.find( key );

    if( valueIt == data_.end() )
      {
	std::cout << "unable to find key: " + key << "\n"
		  << "using default " << key << " = " << def << std::endl;
	return def;
      }

    // extract data
    const std::string value = valueIt->second;

    // assign to string stream
    std::stringstream ss;
    ss << value;

    // extract from sting stream
    T retValue;
    ss >> retValue;

    return retValue;
  }

  std::string getString( const std::string& key ) const
  {
    // get value from data
    const auto valueIt = data_.find( key );

    if( valueIt == data_.end() )
      {
	throw "unable to find key: " + key;
      }

    // extract data
    const std::string value = valueIt->second;
    return value;
  }

  /**
   * \brief get parameter with key plus camera number
   * the format is original_key.camN
   */
  template< class T >
  T get( const std::string& key, const unsigned int cam )
  {
    std::stringstream ss;
    ss << key << ".cam" << cam;
    return get< T >( ss.str() );
  }

  std::string getString( const std::string& key, const unsigned int cam ) const
  {
    std::stringstream ss;
    ss << key << ".cam" << cam;
    return getString( ss.str() );
  }

  /**
   * \brief get vector with key base
   */
  template< class T >
  std::vector< T > getVector( const std::string& key, const unsigned int nCams )
  {
    std::vector< T > ret;
    for( unsigned int cam = 0; cam < nCams; ++cam )
      ret.push_back( get<T>( key, cam ) );

    return ret;
  }

  std::vector< std::string > getStringVector( const std::string& key, const unsigned int nCams )
  {
    std::vector< std::string > ret;
    for( unsigned int cam = 0; cam < nCams; ++cam )
      ret.push_back( getString( key, cam ) );

    return ret;
  }

  int getFromList( const std::string& key, const std::vector< std::string >& names )
  {
    const std::string value = getString( key );
    const auto it = std::find( names.begin(), names.end(), value );

    if( it != names.end() )
      {
	return std::distance( names.begin(), it );
      }
    else
      {
	std::stringstream ss;
	ss << "unable to find value \"" << value << "\" in names:";
	for( const std::string& n : names )
	  ss << " " << n;
	throw ss.str();
      }
  }

  /**
   * \brief add parameter (from command line arguments)
   */
  void add( const char* in )
  {
    std::istringstream is_line( in );
    std::string key;
    if( std::getline( is_line, key, ':' ) )
      {
	std::string value;
	if( std::getline( is_line, value ) )
	  store_line( key, value );
      }
  }

protected:
  /**
   * \brief extract key-value pairs for store_line to read
   */
  void readConfigFile( std::ifstream& stream )
  {
    std::string line;
    while( std::getline( stream, line) )
      {
	std::istringstream is_line(line);
	std::string key;
	if( std::getline(is_line, key, ':') )
	  {
	    std::string value;
	    if( std::getline(is_line, value) )
	      store_line(key, value);
	  }
      }
  }

  /**
   * \brief stores key-value pairs, over writes where necessary
   */
  void store_line( const std::string &key, const std::string& value )
  {
    auto it = data_.insert( std::pair< std::string, std::string >( key, value ) );

    if( it.second == true )
      {
	std::cout << " - " << key << " = " << value << std::endl;
	return;
      }

    // otherwise over writing
    const std::string oldValue = it.first->second;
    it.first->second = value;
    std::cout << " * " << key << " = " << value
	      << " ( was " << oldValue << " )" << std::endl;
  }

private:
  DataType data_;
};

template<typename >
void set( ConfigParameters& c, const std::string& key, const bool& value )
{
	std::string vs=value?"true":"false";
	c.set(vs, vs);
}

template<class T>
T get( ConfigParameters const& c, const std::string& key )
{
	return c.get<T>(key);
}

template<>
bool get( ConfigParameters const& c, const std::string& key )
{
	std::string s=c.get<std::string>(key);

	if(s=="true"){ untested();
		return true;
	}else if(s=="1"){ untested();
		return true;
	}else if(s=="yes"){ untested();
		return true;
	} // ....
	else{ untested();
		return false;
	}
}

template<class T>
T get( ConfigParameters const& c, const std::string& key, T def )
{
	incomplete();
	return c.get<T>(key);
}

template<>
bool get( ConfigParameters const& c, const std::string& key, bool def )
{ untested();
	bool ret;
	try{ untested();
		ret = get<bool>(c, key);
	}catch( std::string ){ untested();
	  	// eek, no BaseException!
		ret = def;
	}

	return ret;
}


#endif // #ifndef PARAMETER_HH
