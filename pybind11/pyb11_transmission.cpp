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
  // py::options options;
  // options.disable_user_defined_docstrings();
  
  py::class_< utl::transmission >( m, "BPT" )
    .def( py::init<
	  const std::vector< double >&,
	  const std::vector< double >& >(), R"(
Class handling Bandpass Transmission filters.
  
This is defined on an input wavelength grid, the 2 inputs should have the same size
  
Parameters
----------
ll : array-like
   input wavelength grid where the transmission is defined
tt : array-like
   transmission, in photons, of the given filter  
)"
	  )
    .def( "get_bandpass_flux", []( const utl::transmission & o,
				   const py::array_t< double > & LL,
				   const py::array_t< double > & FL ) {

				 auto ll = static_cast< const double * >( LL.data() );
				 auto fl = static_cast< const double * >( FL.data() );
				 return o.get_bandpass_flux( ll, fl, LL.shape(0) ); },
	  R"(
Function returning the transmitted integrated flux within the transmittion band

Parameters
----------
LL : array-like
   wavelength grid over which to compute the transmitted flux
FL : array-like
   flux grid, corresponding to the grid LL, for which to compute the transmission

Returns
-------
: float
   transmitted flux within the filter"
)", py::arg("LL"), py::arg("FL") )
    .def_readonly( "lmin", &utl::transmission::lmin, R"(minimum wavelength of the grid)" )
    .def_readonly( "lmax", &utl::transmission::lmax, R"(maximum wavelength of the grid)" )
    .def_readonly( "lpiv", &utl::transmission::lpiv, R"(pivot wavelength of the filter)" )
    .def_readonly( "norm", &utl::transmission::norm, R"(normalisation of the transmission)" )
    .def( "get_lmin", [] ( const utl::transmission & o ) { return o.lmin; },
	  R"(function returning the minimum wavelength of the grid)" )
    .def( "get_lmax", [] ( const utl::transmission & o ) { return o.lmax; },
	  R"(function returning the maximum wavelength of the grid)" )
    .def( "get_lpiv", [] ( const utl::transmission & o ) { return o.lpiv; },
	  R"(function returning pivot wavelength of the filter)" )
    .def( "get_norm", [] ( const utl::transmission & o ) { return o.norm; },
	  R"(functino returning the normalisation of the transmission)" )
    .def( "get_xaxis", &utl::transmission::get_xaxis,
	  R"(function returning the x-axis of the interpolation grid (i.e. wavelength))" )
    .def( "get_yaxis", &utl::transmission::get_yaxis,
	  R"(function returning the y-axis of the interpolation grid (i.e. normalized transmission))" )
    .def(py::pickle(
    		    []( const utl::transmission &o ) { //__getstate__
		      return utl::__getstate__< utl::transmission >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< utl::transmission >( b ); }
    		    ) );

} // end PYBIND11_MODULE( BandpassTransmission, m )

// ==========================================================================================
