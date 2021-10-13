/**
 *  @file sfh/include/insitu_sfh.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __SFH_INSITU_H__
#define __SFH_INSITU_H__

// internal includes
#include <sfh_base.h>

namespace sed {

  class sfh_insitu : public sfh_base {

  private :

    // constant parameters
    const double _s = 3.;
    const double _i_s = 1. / _s;
    const double _kSN = 10;
    const double _eps_acc = 1.e+6;
    const double _yZ = 0.04, _yD = 3.8e-4;

    // tunable parameters
    double _Psi_max = 100.;   // [ Msol * yr^-1 ]
    double _tau_star = 1.e+8; // [ yr ]

    // private functions of the class
    double _sfr ( const double tau ) const noexcept override;
    double _dust_to_gas ( const double tau ) const noexcept;

  public :

    sfh_insitu () : sfh_base{} {}

    sfh_insitu ( const double tau_quench ) noexcept
      : sfh_base { tau_quench } {

      _paramsrc = std::vector< double >( 4 );
      set_params();

    }
    
    // In:
    // param[ 0 ] = Psi_max
    // param[ 1 ] = tau_star
    // Out:
    // idx_0 = Psi_max
    // idx_1 = 1 + 3 * std::pow( Psi_max, -0.3 )
    // idx_2 = tau_star
    // idx_3 = 1. / tau_star
    void set_params ( const double * const param = nullptr ) noexcept override {

      if ( param ) {
	this->set_Psi_max( param[ 0 ] );
	this->set_tau_star( param[ 1 ] );
      }
      else {
	this->set_Psi_max( _Psi_max );
	this->set_tau_star( _tau_star );
      }
      return;
      
    }

    // Overriding galactic content functions
    double get_Mdust ( const double tau ) const noexcept override;
    double get_Mgas  ( const double tau ) const noexcept override;
    double get_Zgas  ( const double tau ) const noexcept override;
    double get_Zstar ( const double tau ) const noexcept override;

    // Tunable parameters Getters/Setters 
    double get_Psi_max () const noexcept { return _Psi_max; }
    double get_tau_star () const noexcept { return _tau_star; }
    void set_Psi_max ( const double Psi_max ) noexcept {
      _Psi_max = Psi_max;
      _paramsrc[ 0 ] = _Psi_max;
      _paramsrc[ 1 ] = 1 + 3 * std::pow( _Psi_max, -0.3 );
      return;
    }
    void set_tau_star ( const double tau_star ) noexcept {
      _tau_star = tau_star;
      _paramsrc[ 2 ] = _tau_star;
      _paramsrc[ 3 ] = 1. / _tau_star;
      return;
    }

  }; // endclass sfh_insitu

} // endnamespace sed

#endif //__SFH_INSITU_H__
