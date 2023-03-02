// External includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// STL includes
#include <exception>
#include <vector>
#include <memory>
#include <string>
// Internal includes
#include <pyb11_serialize.h>
#include <syn.h>

namespace py = pybind11;

// ==========================================================================================

PYBIND11_MODULE( SYN_core, m ) {
  
  py::class_< sed::syn >( m, "CSYN" )
    .def( py::init< const std::vector< double >& >() )
    .def( "set_params",
	  [] ( sed::syn & o,
	       const py::array_t< double > & params ) {
	    o.set_params( static_cast< const double * >( params.data() ) );
	  } )
    .def( "opt_depth", &sed::syn::opt_depth )
    .def( "energy", py::vectorize(&sed::syn::energy) )
    .def(py::pickle(
    		    []( const sed::syn &o ) { //__getstate__
		      return utl::__getstate__< sed::syn >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< sed::syn >( b ); }
    		    ) );

} // end PYBIND11_MODULE( SYN_core, m )

// ==========================================================================================
