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
#include <iostream>

// internal includes
#include <utilities.h>
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

  class ism {

  private :

    const double _logTmin = -1, _logTmax = 4;
    const double _logLmin = 0, _logLmax = 10;

  protected :

    double _Temp;
    std::vector< double > _current_params;
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
    ism ( const double low, const double upp ) : _low{ low }, _upp{ upp } {

      _ext_norm_cont = std::pow( 1.81818181818181818e+2, _upp - _low );

    }
      
    virtual ~ism () = default;

    virtual void set_params ( const double * const param ) noexcept = 0;
    void set_temperature ( const double temp ) { _Temp = temp; }
    
    std::vector< double > get_params () { return _current_params; }

    double extinction ( const double lambda ) const noexcept {

      return _current_params[ 0 ] *
	std::pow( lambda * 1.81818181818181818e-4, -_delta( lambda ) ) *
	delta( lambda, 1, _ext_norm_cont );

    }
    
    double attenuation ( const double lambda ) const noexcept {

      return std::pow( 10, -0.4 * this->extinction( lambda ) );
      
    }

    virtual double emission ( const double lambda ) const noexcept = 0;
    
    virtual double temperature ( const double Etot ) noexcept;

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
    // Out:
    // idx_0 = A_V_diff
    // idx_1 = fact_greybody
    // idx_2 = f_PAH
    void set_params ( const double * const param ) noexcept override {

      _current_params = { _A_Vband( param ),
			  _fact_greybody( param ),
			  param[ 4 ] };
      return;
      
    }

    double emission ( const double lambda ) const noexcept override;

    double temperature ( const double Etot ) noexcept override {

      _Labs = Etot;
      return ism::temperature( ( 1 - _current_params[ 2 ] ) * Etot );

    }

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
    // Out:
    // idx_0 = A_V_MC
    // idx_1 = fact_greybody
    // idx_2 = tau_esc
    // idx_3 = 1. / tau_esc
    void set_params ( const double * const param ) noexcept override {

      _current_params = { _A_Vband( param ),
			  _fact_greybody( param ),
			  param[ 5 ],
			  1 / param[ 5 ] };
      return;
      
    }
    
    double eta ( const double tau ) const noexcept {

      return
        ( tau <= _current_params[ 2 ] ) +
    	( 2. - tau * _current_params[ 3 ] ) *
    	( ( _current_params[ 2 ] < tau ) &
    	  ( tau <= 2 * _current_params[ 2 ] ) );

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
