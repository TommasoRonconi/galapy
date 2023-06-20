/**
 *  @file utl/include/transmission.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __TRANSMISSION_H__
#define __TRANSMISSION_H__

// external includes
#include <cmath>
#include <iostream>

// internal includes
#include <utilities.h>
#include <serialize.h>
#include <interpolation.h>

namespace utl {

  class transmission : public Serializable {

  private :

    // private variables of the class
    utl::interpolator< utl::lin_interp > _ftran {};

    // private functions of the class
    // [ ... ]

  public :

    double lmin = -1.;
    double lmax = -1.;
    double lpiv = 0.;
    double norm = 0.;

    transmission () = default;
    
    transmission ( const std::vector< double > & ll,
		   const std::vector< double > & fl ) {

      // compute the normalization factor and the pivot wavelength
      // by integrating with the trapezoid rule
      for ( std::size_t ii = 1; ii < ll.size(); ++ii ) {
        norm +=
	  ( fl[ ii ]/ll[ ii ] + fl[ ii-1 ]/ll[ ii-1 ] ) *
	  ( ll[ ii ] - ll[ ii-1 ] );
	lpiv +=
	  ( fl[ ii ]*ll[ ii ] + fl[ ii-1 ]*ll[ ii-1 ] ) *
	  ( ll[ ii ] - ll[ ii-1 ] );
      }
      norm = 2. / norm;
      lpiv = std::sqrt( 0.5 * lpiv * norm );

      // compute the normalized transmission function
      std::vector< double > ff ( fl.size() );
      for ( std::size_t ii = 0; ii < ll.size(); ++ii )
	ff[ ii ] = norm * fl[ ii ] / ll[ ii ];

      // build the interpolator
      _ftran = utl::interpolator< utl::lin_interp > { ll, ff, "linear" };
      
      // store the minimum and the maximum wavelength
      lmin = _ftran.get_xmin();
      lmax = _ftran.get_xmax();
      
    }
    virtual ~transmission () = default;
    
    inline double operator() ( const double ll ) const noexcept { return _ftran( ll ); }

    double get_bandpass_flux ( const double * const ll,
			       const double * const fl,
			       const std::size_t len ) const noexcept{

      double flux = 0.;
      for ( std::size_t ii = 1; ii < len; ++ii )
	flux +=
	  ( _ftran( ll[ ii ] )*fl[ ii ] + _ftran( ll[ ii-1 ] )*fl[ ii-1 ] ) *
	  ( ll[ ii ] - ll[ ii-1 ] );
      return flux * 0.5;
	
    }

    std::vector< double > get_xaxis () { return _ftran.get_xv(); }
    std::vector< double > get_yaxis () { return _ftran.get_fv(); }

    // =============================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	_ftran.serialize_size() +
	SerialPOD< double >::serialize_size( lmin ) +
	SerialPOD< double >::serialize_size( lmax ) +
	SerialPOD< double >::serialize_size( lpiv ) +
	SerialPOD< double >::serialize_size( norm );

    }

    virtual char * serialize ( char * data ) const {

      data = _ftran.serialize( data );
      data = SerialPOD< double >::serialize( data, lmin );
      data = SerialPOD< double >::serialize( data, lmax );
      data = SerialPOD< double >::serialize( data, lpiv );
      data = SerialPOD< double >::serialize( data, norm );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = _ftran.deserialize( data );
      data = SerialPOD< double >::deserialize( data, lmin );
      data = SerialPOD< double >::deserialize( data, lmax );
      data = SerialPOD< double >::deserialize( data, lpiv );
      data = SerialPOD< double >::deserialize( data, norm );
      return data;
      
    }

    // =============================================================================

  }; // endclass transmission

} // endnamespace utl

#endif //__TRANSMISSION_H__
