// External includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// STL includes
#include <exception>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
// Internal includes
#include <pyb11_serialize.h>
#include <csp.h>

namespace py = pybind11;

// ==========================================================================================

namespace sed {

  py::array_t< float >
  get_csp_emission ( const sed::csp & ccsp,
		     const py::array_t<
		     std::size_t,
		     py::array::c_style | py::array::forcecast > & il,
		     const py::array_t<
		     double,
		     py::array::c_style | py::array::forcecast > & Tfact ) {

    auto res = py::array_t<float>(il.shape(0));

    auto ptrR = static_cast<float*>(res.mutable_data());
    auto ptrL = static_cast<const std::size_t*>(il.data());

    if ( Tfact.shape(0) > 0 ) {
      auto ptrT = static_cast<const double*>(Tfact.data());
      for ( py::ssize_t ii = 0; ii < il.shape(0); ++ii )
	ptrR[ii] = ccsp.emission( ptrL[ii], ptrT + ii * ccsp.tau_size() );    
    }
    else {
      auto ptrT = new double [ ccsp.tau_size() ];
      std::fill( ptrT, ptrT + ccsp.tau_size(), 1.0 );
      for ( py::ssize_t ii = 0; ii < il.shape(0); ++ii )
	ptrR[ii] = ccsp.emission( ptrL[ii], ptrT );
      delete [] ptrT;
    }

    return res;
  
  }

// ==========================================================================================

  py::array_t< float >
  get_csp_kernel_emission ( const sed::csp & ccsp,
			    const py::array_t<
			    std::size_t,
			    py::array::c_style | py::array::forcecast > & il,
			    const py::array_t<
			    double,
			    py::array::c_style | py::array::forcecast > & MCfact,
			    const py::array_t<
			    double,
			    py::array::c_style | py::array::forcecast > & DDfact ) {

    auto res0 = py::array_t<float>(il.shape(0));
    auto res1 = py::array_t<float>(il.shape(0));
    auto res2 = py::array_t<float>(il.shape(0));

    auto ptrR0 = static_cast<float*>(res0.mutable_data());
    auto ptrR1 = static_cast<float*>(res1.mutable_data());
    auto ptrR2 = static_cast<float*>(res2.mutable_data());
    auto ptrIL = static_cast<const std::size_t*>(il.data());
    auto ptrMC = static_cast<const double*>(MCfact.data());
    auto ptrDD = static_cast<const double*>(DDfact.data());

    std::vector< double > out;
    for ( py::ssize_t ii = 0; ii < il.shape(0); ++ii ) {
      out = ccsp.kernel_emission( ptrIL[ii],
				  ptrMC + ii * ccsp.tau_size(),
				  ptrDD + ii * ccsp.tau_size() );
      ptrR0[ii] = out[0];
      ptrR1[ii] = out[1];
      ptrR2[ii] = out[2];
    }

    return py::make_tuple( res0, res1, res2 );
  
  }

  // ========================================================================================

  py::tuple loadSSP ( const std::string & file_name ) {

    // Open binary input stream      
    std::ifstream ifs ( file_name, std::ios_base::in | std::ios_base::binary );
    if ( ifs.fail() )
      throw std::ios_base::failure( "Unable to open input file." );
      
    // Read from binary file
    std::size_t Nl;
    ifs.read( ( char * ) &Nl, sizeof( std::size_t ) );
    auto lambda_a = py::array_t<double>( Nl );
    auto lambda   = static_cast<double*>(lambda_a.mutable_data());
    ifs.read( reinterpret_cast< char* >( lambda ), sizeof( double ) * Nl );
    
    std::size_t Nt;
    ifs.read( ( char * ) &Nt, sizeof( std::size_t ) );
    auto tau_a = py::array_t<double>( Nt );
    auto tau   = static_cast<double*>(tau_a.mutable_data());
    ifs.read( reinterpret_cast< char* >( tau ), sizeof( double ) * Nt );
    
    std::size_t NZ;
    ifs.read( ( char * ) &NZ, sizeof( std::size_t ) );
    auto Z_a = py::array_t<double>( NZ );
    auto Z   = static_cast<double*>(Z_a.mutable_data());
    ifs.read( reinterpret_cast< char* >( Z ), sizeof( double ) * NZ );
    
    std::size_t NSSP;
    ifs.read( ( char * ) &NSSP, sizeof( std::size_t ) );
    auto SSP_a = py::array_t<double>( NSSP );
    auto SSP   = static_cast<double*>(SSP_a.mutable_data());
    ifs.read( reinterpret_cast< char* >( SSP ), sizeof( double ) * NSSP );

    // Raise RuntimeError if EoF is not reached:
    if ( !ifs ) {
      ifs.clear(); ifs.close();
      throw std::ios_base::failure( "Error reading input file: "
				    "reading ended without reaching EoF" );
    }
    ifs.clear(); ifs.close();
    
    return py::make_tuple( lambda_a, tau_a, Z_a, SSP_a );
    
  } // endfunc loadSSP

} // endnamespace sed

// ==========================================================================================

PYBIND11_MODULE( CSP_core, m ) {

  py::class_< sed::csp >( m, "CCSP" )
    .def( py::init<
	  const std::vector< double >&,
	  const std::vector< double >&,
	  const std::vector< double >&,
	  const std::vector< double >&,
	  const int
	  >(),
	  py::arg("lambda"), py::arg("tau"), py::arg("Z"),
	  py::arg("LSSP"), py::arg("do_CCSN") )
    .def( "set_params", &sed::csp::set_params )
    .def( "SSP",        &sed::csp::luminosity )
    .def( "emission",   &sed::get_csp_emission,
	  "computes the composite stellar population emission",
	  py::arg("il"), py::arg("Tfact") = py::array_t<double>() )
    .def( "_kernel_emission",   &sed::get_csp_kernel_emission,
	  "computes the composite stellar population emission "
	  "along with the two phases of attenuation",
	  py::arg("il"), py::arg("MCfact"), py::arg("DDfact") )
    .def( "RCCSN",      &sed::csp::RCCSN )
    .def( py::pickle(
		     []( const sed::csp &o ) { //__getstate__
		       return utl::__getstate__< sed::csp >( o ); },
		     []( const py::bytes b ) { //__setstate__
		       return utl::__setstate__< sed::csp >( b ); }
		     ) );

  m.def( "loadSSP", &sed::loadSSP );

} // end PYBIND11_MODULE( CSP_core, m )

// ==========================================================================================
