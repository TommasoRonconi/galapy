/**
 *  @file ism/include/ism.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __ISM_H__
#define __ISM_H__

// external includes
#include <cmath>

// internal includes
#include <utilities.h>
#include <serialize.h>
#include <interpolation.h>
#include <pah.h>

namespace sed {

  // in units of \f$\AA\f$ ( 100 micron = 10^6 angstrom )
  inline double
  delta ( const double lambda,
	  const float low,
	  const float high ) {
    
    return low * ( lambda <= 1.e+6 ) + high * ( lambda > 1.e+6 );

  }
  
  // ================================================================
  // ISM Base class

  class ism : public Serializable {

  private :

    const double _logTmin = -1, _logTmax = 4;
    const double _logLmin = 0, _logLmax = 10;

  protected :

    double _Temp;
    std::vector< double > _paramsrc;
    double _low, _upp, _ext_norm_cont;
    virtual double _A_Vband ( const double * const param = nullptr ) const noexcept = 0;
    inline double _delta ( const double lambda ) const noexcept {

      return delta( lambda, _low, _upp );
      
    }
      
    virtual double _fact_greybody ( const double * const param = nullptr ) const noexcept = 0;

    inline double _emission ( const double lambda,
			      const double temp,
			      const double fact ) const noexcept {

      return fact * ( 1 - attenuation( lambda ) ) * utl::blackbody( lambda * 1.e-8, temp );

    }
    
  public :

    ism () = default;
    ism ( const double low, const double upp ) { set_slopes( low, upp ); }
    // ism ( const double low, const double upp ) : _low{ low }, _upp{ upp } {
    //   _ext_norm_cont = std::pow( 1.81818181818181818e+2, _upp - _low );
    // }
      
    virtual ~ism () = default;

    virtual void set_params ( const double * const param ) noexcept = 0;
    void set_temperature ( const double temp ) { _Temp = temp; }

    void set_slopes ( const double low, const double upp ) noexcept {
      
      _low = low;
      _upp = upp;

      // ( 100 micron / 5500 Angstrom ) = 1.818181[...]*10^2 Angstrom
      _ext_norm_cont = std::pow( 1.81818181818181818e+2, _upp - _low );
      
    }
    
    std::vector< double > get_params () { return _paramsrc; }

    double extinction ( const double lambda ) const noexcept {

      return _paramsrc[ 0 ] *
	std::pow( lambda * 1.81818181818181818e-4, -_delta( lambda ) ) *
	delta( lambda, 1, _ext_norm_cont );

    }
    
    double attenuation ( const double lambda ) const noexcept {

      // double ret = std::pow( 10, -0.4 * this->extinction( lambda ) );
      // return ( ret < 1.0 ) ? ret : 0.99999;
      return std::pow( 10, -0.4 * this->extinction( lambda ) );
      
    }

    virtual double emission ( const double lambda ) const noexcept = 0;
    
    virtual double temperature ( const double Etot ) noexcept;

    // ================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	SerialVecPOD< double >::serialize_size( _paramsrc ) +
	SerialPOD< double >::serialize_size( _Temp ) +
	SerialPOD< double >::serialize_size( _low ) +
	SerialPOD< double >::serialize_size( _upp ) +
	SerialPOD< double >::serialize_size( _ext_norm_cont );

    }

    virtual char * serialize ( char * data ) const {

      data = SerialVecPOD< double >::serialize( data, _paramsrc );
      data = SerialPOD< double >::serialize( data, _Temp );
      data = SerialPOD< double >::serialize( data, _low );
      data = SerialPOD< double >::serialize( data, _upp );
      data = SerialPOD< double >::serialize( data, _ext_norm_cont );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = SerialVecPOD< double >::deserialize( data, _paramsrc );
      data = SerialPOD< double >::deserialize( data, _Temp );
      data = SerialPOD< double >::deserialize( data, _low );
      data = SerialPOD< double >::deserialize( data, _upp );
      data = SerialPOD< double >::deserialize( data, _ext_norm_cont );
      return data;

    }

    // ================================================================

  }; // endclass ism
  
  // ================================================================
  // diffuse class

  class diffuse : public ism {

  private :

    double _Labs = 0.;
    const utl::interpolator< utl::lin_interp > _fpah { pah::lpah, pah::fpah, "linear" };
    
  protected :

    // In:
    // param[ 0 ] = f_MC
    // param[ 1 ] = N_V_band
    // param[ 2 ] = M_dust
    // param[ 3 ] = R_dust [ pc ]
    virtual double _A_Vband ( const double * const param = nullptr ) const noexcept override;
    virtual double _fact_greybody ( const double * const param = nullptr ) const noexcept override;
    
  public :
    
    diffuse () : ism{ 0.7, 2 } {}
    ~diffuse () = default;
        
    // In:
    // param[ 0 ] = f_MC
    // param[ 1 ] = N_V^diff
    // param[ 2 ] = M_dust
    // param[ 3 ] = R_dust [ pc ]
    // param[ 4 ] = f_PAH
    // param[ 5 ] = delta_low
    // param[ 6 ] = delta_upp
    // Out:
    // idx_0 = A_V_diff
    // idx_1 = fact_greybody
    // idx_2 = f_PAH
    void set_params ( const double * const param ) noexcept override {

      _paramsrc = { _A_Vband( param ),
		    _fact_greybody( param ),
		    param[ 4 ] };
      set_slopes( param[ 5 ], param[ 6 ] );
      return;
      
    }

    double emission ( const double lambda ) const noexcept override;

    double temperature ( const double Etot ) noexcept override {

      _Labs = Etot;
      return ism::temperature( ( 1 - _paramsrc[ 2 ] ) * Etot );

    }

    // ================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	ism::serialize_size() +
	SerialPOD< double >::serialize_size( _Labs );

    }

    virtual char * serialize ( char * data ) const {

      data = ism::serialize( data );
      data = SerialPOD< double >::serialize( data, _Labs );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = ism::deserialize( data );
      data = SerialPOD< double >::deserialize( data, _Labs );
      return data;

    }

    // ================================================================

  }; // endclass diffuse

  // ================================================================
  // cloud class
  
  class cloud : public ism {

  private :

  protected :

    virtual double _A_Vband ( const double * const param = nullptr ) const noexcept override;
    virtual double _fact_greybody ( const double * const param = nullptr ) const noexcept override;
    
  public :

    cloud () : ism{ 1.3, 1.6 } {}
    ~cloud () = default;
    
    // In:
    // param[ 0 ] = f_MC
    // param[ 1 ] = N_V^MC
    // param[ 2 ] = N_MC
    // param[ 3 ] = R_MC [ pc ]
    // param[ 4 ] = Z_gas / Z_sol
    // param[ 5 ] = tau_esc
    // param[ 6 ] = M_gas
    // param[ 7 ] = delta_low
    // param[ 8 ] = delta_upp
    // Out:
    // idx_0 = A_V_MC
    // idx_1 = fact_greybody
    // idx_2 = tau_esc
    // idx_3 = 1. / tau_esc
    void set_params ( const double * const param ) noexcept override {

      _paramsrc = { _A_Vband( param ),
		    _fact_greybody( param ),
		    param[ 5 ],
		    1 / param[ 5 ] };
      set_slopes( param[ 7 ], param[ 8 ] );
      return;
      
    }
    
    double eta ( const double tau ) const noexcept {

      return
        ( tau <= _paramsrc[ 2 ] ) +
    	( 2. - tau * _paramsrc[ 3 ] ) *
    	( ( _paramsrc[ 2 ] < tau ) &
    	  ( tau <= 2 * _paramsrc[ 2 ] ) );

    }

    double emission ( const double lambda ) const noexcept override;

    double time_attenuation ( const double lambda, const double tau ) const noexcept;

    double average_attenuation ( const double lambda, const double tau ) const noexcept;
    
  }; // endclass cloud

  // ================================================================
  // Non-member functions

  double total_attenuation ( const double wavelenght,
			     const double att_dd_ll,
			     const double att_mc_ll,
			     const double eta_mc_tt );

} // endnamespace sed

#endif //__ISM_H__
