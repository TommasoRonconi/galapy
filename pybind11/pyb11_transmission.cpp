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
#include <transmission.h>

namespace py = pybind11;

// ==========================================================================================

PYBIND11_MODULE( BandpassTransmission, m ) {
  
  py::class_< utl::transmission >( m, "BPT" )
    .def( py::init<
	  const std::vector< double >&,
	  const std::vector< double >& >() )
    .def( "get_bandpass_flux", []( const utl::transmission & o,
				   const py::array_t< double > & LL,
				   const py::array_t< double > & FL ) {

				 auto ll = static_cast< const double * >( LL.data() );
				 auto fl = static_cast< const double * >( FL.data() );
				 return o.get_bandpass_flux( ll, fl, LL.shape(0) );
				 
			       } )
    .def_readonly( "lmin", &utl::transmission::lmin )
    .def_readonly( "lmax", &utl::transmission::lmax )
    .def_readonly( "lpiv", &utl::transmission::lpiv )
    .def_readonly( "norm", &utl::transmission::norm )
    .def( "get_lmin", [] ( const utl::transmission & o ) { return o.lmin; } )
    .def( "get_lmax", [] ( const utl::transmission & o ) { return o.lmax; } )
    .def( "get_lpiv", [] ( const utl::transmission & o ) { return o.lpiv; } )
    .def( "get_norm", [] ( const utl::transmission & o ) { return o.norm; } )
    .def( "get_xaxis", &utl::transmission::get_xaxis )
    .def( "get_yaxis", &utl::transmission::get_yaxis )
    .def(py::pickle(
    		    []( const utl::transmission &o ) { //__getstate__
		      return utl::__getstate__< utl::transmission >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< utl::transmission >( b ); }
    		    ) );

} // end PYBIND11_MODULE( BandpassTransmission, m )

// ==========================================================================================
