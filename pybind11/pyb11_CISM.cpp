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
#include <ism.h>

namespace py = pybind11;

// ==========================================================================================

PYBIND11_MODULE( ISM_core, m ) {
  
  py::class_< sed::diffuse >( m, "CDD" )
    .def( py::init<>() )
    .def( "set_params",
	  [] ( sed::diffuse & o,
	       const py::array_t< double > & params ) {
	    o.set_params( static_cast< const double * >( params.data() ) );
	  } )
    .def( "set_temperature", &sed::diffuse::set_temperature,
	  "Sets internally the temperature" )
    .def( "set_slopes",  &sed::diffuse::set_slopes )
    .def( "temperature", []( sed::diffuse & o, const double Etot ){
			   return o.temperature( Etot ); },
	  "Computes temperature necessary for "
	  "energy balance and returns it",
	  py::arg("Etot") )
    .def( "extinction",  py::vectorize(&sed::diffuse::extinction) )
    .def( "attenuation", py::vectorize(&sed::diffuse::attenuation) )
    .def( "emission",    py::vectorize(&sed::diffuse::emission) )
    .def( "A_V", []( const sed::diffuse & o ){ return o.get_params()[0]; }  )
    .def(py::pickle(
    		    []( const sed::diffuse &o ) { //__getstate__
		      return utl::__getstate__< sed::diffuse >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< sed::diffuse >( b ); }
    		    ) )
    .def( "T", &sed::diffuse::get_temperature );
  
  // ========================================================================================
  
  py::class_< sed::cloud >( m, "CMC" )
    .def( py::init<>() )
    .def( "set_params",
	  [] ( sed::cloud & o,
	       const py::array_t< double > & params ) {
	    o.set_params( static_cast< const double * >( params.data() ) );
	  } )
    .def( "set_temperature", &sed::cloud::set_temperature,
	  "Sets internally the temperature" )
    .def( "set_slopes",  &sed::cloud::set_slopes )
    .def( "temperature", []( sed::cloud & o, const double Etot ){
			   return o.temperature( Etot ); },
	  "Computes temperature necessary for "
	  "energy balance and returns it",
	  py::arg("Etot") )
    .def( "extinction",  py::vectorize(&sed::cloud::extinction) )
    .def( "attenuation", py::vectorize(&sed::cloud::attenuation) )
    .def( "emission",    py::vectorize(&sed::cloud::emission) )
    .def( "eta", py::vectorize(&sed::cloud::eta) )
    .def( "A_V", []( const sed::cloud & o ){ return o.get_params()[0]; }  )
    .def(py::pickle(
    		    []( const sed::cloud &o ) { //__getstate__
		      return utl::__getstate__< sed::cloud >( o ); },
    		    []( const py::bytes b ) {     //__setstate__
		      return utl::__setstate__< sed::cloud >( b ); }
    		    ) )
    .def( "T", &sed::cloud::get_temperature );

  m.def( "total_attenuation", py::vectorize(&sed::total_attenuation) );

} // end PYBIND11_MODULE( ISM_core, m )

// ==========================================================================================
