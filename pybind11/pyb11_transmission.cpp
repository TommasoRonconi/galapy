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
	  const std::vector< double >& >(),
	  "Class handling Bandpass Transmission filters."
	  "\n"
	  "\nThis is defined on an input wavelength grid, the 2 inputs should have the same size"
	  "\n"
	  "\nParameters"
	  "\n----------"
	  "\nll : array-like"
	  "\n\tinput wavelength grid where the transmission is defined"
	  "\ntt : array-like"
	  "\n\ttransmission, in photons, of the given filter"
	  )
    .def( "get_bandpass_flux", []( const utl::transmission & o,
				   const py::array_t< double > & LL,
				   const py::array_t< double > & FL ) {

				 auto ll = static_cast< const double * >( LL.data() );
				 auto fl = static_cast< const double * >( FL.data() );
				 return o.get_bandpass_flux( ll, fl, LL.shape(0) ); },
	  "function returning the transmitted integrated flux within the transmittion band\n"
	  "\nParameters"
	  "\n----------"
	  "\nLL : array-like"
	  "\n\twavelength grid over which to compute the transmitted flux"
	  "\nFL : array-like"
	  "\n\tflux grid, corresponding to the grid LL, for which to compute the transmission"
	  "\n"
	  "\nReturns"
	  "\n-------"
	  ": float"
	  "\ttransmitted flux within the filter", py::arg("LL"), py::arg("FL")
	  )
    .def_readonly( "lmin", &utl::transmission::lmin, "minimum wavelength of the grid" )
    .def_readonly( "lmax", &utl::transmission::lmax, "maximum wavelength of the grid" )
    .def_readonly( "lpiv", &utl::transmission::lpiv, "pivot wavelength of the filter" )
    .def_readonly( "norm", &utl::transmission::norm, "normalisation of the transmission" )
    .def( "get_lmin", [] ( const utl::transmission & o ) { return o.lmin; },
	  "function returning the minimum wavelength of the grid" )
    .def( "get_lmax", [] ( const utl::transmission & o ) { return o.lmax; },
	  "function returning the maximum wavelength of the grid" )
    .def( "get_lpiv", [] ( const utl::transmission & o ) { return o.lpiv; },
	  "function returning pivot wavelength of the filter" )
    .def( "get_norm", [] ( const utl::transmission & o ) { return o.norm; },
	  "functino returning the normalisation of the transmission" )
    .def( "get_xaxis", &utl::transmission::get_xaxis,
	  "function returning the x-axis of the interpolation grid (i.e. wavelength)" )
    .def( "get_yaxis", &utl::transmission::get_yaxis,
	  "function returning the y-axis of the interpolation grid "
	  "(i.e. normalized transmission)" )
    .def(py::pickle(
    		    []( const utl::transmission &o ) { //__getstate__
		      return utl::__getstate__< utl::transmission >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< utl::transmission >( b ); }
    		    ) );

} // end PYBIND11_MODULE( BandpassTransmission, m )

// ==========================================================================================
