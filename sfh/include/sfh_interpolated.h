/**
 *  @file sfh/include/sfh_interpolated.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __SFH_INTERPOLATED_H__
#define __SFH_INTERPOLATED_H__

// internal includes
#include <serialize.h>
#include <sfh_empirical.h>
#include <interpolation.h>

namespace sed {

  // ================================================================
  // Interpolated SFH model
  
  class sfh_interpolated : public sfh_empirical {

  private :

    utl::interpolator< utl::lin_interp > _fsfh {};

    // private functions of the class
    double _sfr ( const double tau ) const noexcept override {

      return _fsfh( tau );

    }

  protected :

  public :

    sfh_interpolated () : sfh_empirical{} {}

    sfh_interpolated ( const double tau_quench,
		       const std::vector< double > & tau,
		       const std::vector< double > & sfr ) noexcept
      : sfh_empirical { tau_quench } {

      _fsfh = utl::interpolator< utl::lin_interp > { tau, sfr };

    }

    void set_interpolator ( const std::vector< double > & tau,
			    const std::vector< double > & sfr ) noexcept override {

      _fsfh = utl::interpolator< utl::lin_interp > { tau, sfr };
      return;

    }

    // In:
    // param[ 0 ] = M_dust
    // param[ 1 ] = Z_gxy
    void set_params ( const double * const param = nullptr ) noexcept override {

      if ( param ) {
	this->set_Mdust( param[ 0 ] );
	this->set_Zgxy( param[ 1 ] );
      }
      else {
	this->set_Mdust( _Mdust );
	this->set_Zgxy( _Zgxy );
      }
      return;
      
    }


    // =============================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
    	sfh_empirical::serialize_size() +
        _fsfh.serialize_size();
      
    }

    virtual char * serialize ( char * data ) const {
      
      data = sfh_empirical::serialize( data );
      data = _fsfh.serialize( data );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = sfh_empirical::deserialize( data );
      data = _fsfh.deserialize( data );
      return data;
      
    }

    // =============================================================

  }; // endclass sfh_interpolated


} // endnamespace sed

#endif //__SFH_INTERPOLATED_H__
