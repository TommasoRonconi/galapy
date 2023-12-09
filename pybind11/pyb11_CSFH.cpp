// External includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// STL includes
#include <exception>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
// Internal includes
#include <pyb11_serialize.h>
#include <sfh.h>

namespace py = pybind11;

// ==========================================================================================

std::unique_ptr<sed::sfh_base> set_sfh_model ( const std::string & model ) {

  if ( model == "insitu" )
    return std::make_unique<sed::sfh_insitu>();

  if ( model == "constant" )
    return std::make_unique<sed::sfh_constant>();

  if ( model == "delayedexp" )
    return std::make_unique<sed::sfh_delayedexp>();

  if ( model == "lognormal" )
    return std::make_unique<sed::sfh_lognorm>();

  if ( model == "interpolated" )
    return std::make_unique<sed::sfh_interpolated>();

  // if ( model == "burst" )
  //   return new sed::sfh_burst{};

  throw std::invalid_argument( "SFH model not valid. "
			       "Valid models are: "
			       "'insitu', 'constant', "
			       "'delayedexp', 'lognormal', "
			       "'interpolated'" );
  
}

// ==========================================================================================

namespace sed {
  struct CSFH {
    
    std::unique_ptr< sed::sfh_base > ptr = nullptr;
    std::string model = "";

    CSFH ( const double tau_quench = 2.e+10,
	   const std::string & model = "insitu" ) :
      ptr{ set_sfh_model( model ) }, model { model } {

      ptr->set_tau_quench( tau_quench );

    }
    CSFH ( const CSFH & other ) {

      model = other.model;
      ptr = set_sfh_model( model );
      if ( other.ptr )
	*ptr = *other.ptr;

    }
    ~CSFH () = default;
    double operator() ( const double xx ) const noexcept { return ( *ptr )( xx ); }
    void set_params ( const std::vector< double > & param ) noexcept {
      ptr->set_params( param.data() );
      return;
    }
    void set_tau_quench ( const double tau_q ) { ptr->set_tau_quench( tau_q ); }
    void set_interpolator ( const std::vector< double > & x,
			    const std::vector< double > & y ) { ptr->set_interpolator( x, y ); }
    double Mstar ( const double tau, const std::size_t npoints = 1000 ) const noexcept {
      return ptr->get_Mstar( tau, npoints );
    }
    double Mdust ( const double tau ) const noexcept { return ptr-> get_Mdust( tau ); }
    double Mgas ( const double tau ) const noexcept { return ptr-> get_Mgas( tau ); }
    double Zgas ( const double tau ) const noexcept { return ptr-> get_Zgas( tau ); }
    double Zstar ( const double tau ) const noexcept { return ptr-> get_Zstar( tau ); }
    py::array_t< double > dMstar () const noexcept {
      auto dM = ptr->get_dMstar();
      return py::array_t< double >( dM.size(), dM.data() );
    }

    std::size_t serialize_size () const { return ptr->serialize_size(); }
    char * serialize ( char * data ) const { return ptr->serialize( data ); }
    const char * deserialize ( const char * data ) { return ptr->deserialize(data); }
    
  }; // endstruct CSFH
} // endnamespace sed

// ==========================================================================================

PYBIND11_MODULE( SFH_core, m ) {
  
  py::class_< sed::CSFH >( m, "CSFH" )
    .def( py::init< const double, const std::string & >(),
	  py::arg("tau_quench") = 2.e+10,
	  py::arg("model") = "insitu" )
    .def( "__call__", py::vectorize(&sed::CSFH::operator()),
	  "Compute SFR at given time",
	  py::arg("tau"))
    .def( "set_params",       &sed::CSFH::set_params )
    .def( "set_tau_quench",   &sed::CSFH::set_tau_quench )
    .def( "set_interpolator", &sed::CSFH::set_interpolator )
    .def( "Mstar", py::vectorize(&sed::CSFH::Mstar) )
    .def( "Mdust", py::vectorize(&sed::CSFH::Mdust) )
    .def( "Mgas",  py::vectorize(&sed::CSFH::Mgas) )
    .def( "Zgas",  py::vectorize(&sed::CSFH::Zgas) )
    .def( "Zstar", py::vectorize(&sed::CSFH::Zstar) )
    .def( "dMstar", &sed::CSFH::dMstar )
    .def( "time_grid",
	  []( const sed::CSFH & csfh,
	      const double age,
	      const std::vector< double > & tgrid,
	      const std::vector< double > & Zgrid ) {

	    std::vector< double > out_dMgrid (tgrid.size()), out_Zgrid(tgrid.size());
	    std::vector< std::size_t > out_Zidx(tgrid.size());
	    std::size_t out_last_idx;

	    csfh.ptr->time_grid( age, tgrid, Zgrid,
				 out_dMgrid, out_Zgrid, out_Zidx, out_last_idx );
	    
	    return py::make_tuple( py::array_t< double >( out_dMgrid.size(),
							  out_dMgrid.data() ),
				   py::array_t< double >( out_Zgrid.size(),
							  out_Zgrid.data() ),
				   py::array_t< std::size_t >( out_Zidx.size(),
							       out_Zidx.data() ),
				   out_last_idx );
	    
	  } )
    .def(py::pickle(
    		    []( const sed::CSFH &o ) { //__getstate__

		      py::bytes b = utl::__getstate__< sed::CSFH >( o );
		      return py::make_tuple( py::float_( o.ptr->get_tau_quench() ),
					     py::str( o.model ), b );

		    },
    		    []( const py::tuple t ) { //__setstate__

		      sed::CSFH o { t[0].cast< double >(), t[1].cast< std::string >() };
		      std::string strtmp { t[2].cast< std::string >() };
		      o.deserialize( strtmp.c_str() );
		      return o;
		      
		    }
    		    ) );

} // end PYBIND11_MODULE( SFH_core, m )

// ==========================================================================================
