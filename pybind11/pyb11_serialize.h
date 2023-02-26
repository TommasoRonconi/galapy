#include <pybind11/pybind11.h>
#include <memory>
#include <string>
#include <exception>
#include <iostream>

namespace py = pybind11;

namespace utl {
  
  template < class T >
  py::bytes __getstate__ ( const T & o ) {
  
    std::size_t size = o.serialize_size();
    char * data = new char[ size ];
    o.serialize( data );

    auto out = py::bytes( data, size );
    delete [] data;

    return out;

  }

  template< class T >
  T __setstate__ ( const py::bytes & b ) {

    std::string strtmp { b.cast< std::string >() };
    T out {};
    
    out.deserialize( strtmp.c_str() );
    return out;

  }

} // endnamespace utl
