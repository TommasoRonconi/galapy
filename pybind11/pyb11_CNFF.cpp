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
#include <nff.h>

namespace py = pybind11;

// ==========================================================================================

PYBIND11_MODULE( NFF_core, m ) {
  
  py::class_< sed::nff >( m, "CNFF" )
    .def( py::init< const std::vector< double >& >() )
    .def( "set_params",
	  [] ( sed::nff & o,
	       const py::array_t< double > & params ) {
	    o.set_params( static_cast< const double * >( params.data() ) );
	  } )
    .def( "Te", [] ( const sed::nff & o, const double Zgas ){
		  return std::exp( o.lTe( std::log10( 50 * Zgas ) ) );
		} )
    .def( "gff", [] ( const sed::nff & o,
		      const std::size_t il,
		      const double Te ) {
		   return o.gff( il, std::log( Te ) - 4 * sed::cnst::ln_10 );
		 } )
    .def( "emission", py::vectorize(&sed::nff::emission) )
    .def(py::pickle(
    		    []( const sed::nff &o ) { //__getstate__
		      return utl::__getstate__< sed::nff >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< sed::nff >( b ); }
    		    ) );

} // end PYBIND11_MODULE( NFF_core, m )

// ==========================================================================================
