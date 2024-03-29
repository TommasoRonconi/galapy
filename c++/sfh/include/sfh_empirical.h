/**
 *  @file sfh/include/sfh_empirical.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __SFH_EMPIRICAL_H__
#define __SFH_EMPIRICAL_H__

// internal includes
#include <serialize.h>
#include <sfh_base.h>

namespace sed {

  // ================================================================
  // Base class of empirical SFH models
  
  class sfh_empirical : public sfh_base {

  private :
    
    // private functions of the class
    virtual double _gas_to_dust ( const double tau ) const noexcept {

      return 100 * std::pow( get_Zgas( tau ) * sed::cnst::solZ, 0.85 );

    }

  protected :

    // tunable parameters
    double _Mdust = 1.e+8; // [ Msol ]
    double _Zgxy = 0.01;    // [ percent / 100 ]

  public :

    sfh_empirical () : sfh_base{} {}

    sfh_empirical ( const double tau_quench ) noexcept
      : sfh_base { tau_quench } {}

    // Overriding Getters for tunable parameters
    double get_Mdust ( __attribute__((unused)) const double tau )
      const noexcept override { return _Mdust; }
    double get_Mgas  ( const double tau )
      const noexcept override { return _Mdust * _gas_to_dust( tau ); }
    double get_Zgas  ( __attribute__((unused)) const double tau )
      const noexcept override { return _Zgxy; }
    double get_Zstar ( __attribute__((unused)) const double tau )
      const noexcept override { return _Zgxy; }

    // Tunable parameter setters
    void set_Mdust ( const double Mdust ) { _Mdust = Mdust; return; }
    void set_Zgxy ( const double Zgxy ) { _Zgxy = Zgxy; return; }

    // =============================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	sfh_base::serialize_size() +
	SerialPOD< double >::serialize_size( _Mdust ) +
	SerialPOD< double >::serialize_size( _Zgxy );

    }

    virtual char * serialize ( char * data ) const {

      data = sfh_base::serialize( data );
      data = SerialPOD< double >::serialize( data, _Mdust );
      data = SerialPOD< double >::serialize( data, _Zgxy );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = sfh_base::deserialize( data );
      data = SerialPOD< double >::deserialize( data, _Mdust );
      data = SerialPOD< double >::deserialize( data, _Zgxy );
      return data;

    }

    // =============================================================

  }; // endclass sfh_empirical

  // ================================================================
  // Constant SFH

  class sfh_constant : public sfh_empirical {

  private :

    // tunable parameters
    double _Psi = 1.; // [ Msol * yr^-1 ]
    
    // private functions of the class
    double _sfr ( __attribute__((unused)) const double tau ) const noexcept override {

      return _paramsrc[ 0 ];

    }

  public :

    // sfh_constant () : sfh_empirical{} { _paramsrc = std::vector< double >( 1 ); }

    sfh_constant ( const double tau_quench = 2.e+10 ) noexcept
      : sfh_empirical { tau_quench } {

      _paramsrc = std::vector< double >( 1 );
      set_params();
      
    }

    // In:
    // param[ 0 ] = Psi
    // param[ 1 ] = M_dust
    // param[ 2 ] = Z_gxy
    void set_params ( const double * const param = nullptr ) noexcept override {

      if ( param ) {
	this->set_Psi( param[ 0 ] );
	this->set_Mdust( param[ 1 ] );
	this->set_Zgxy( param[ 2 ] );
      }
      else {
	this->set_Psi( _Psi );
	this->set_Mdust( _Mdust );
	this->set_Zgxy( _Zgxy );
      }
      return;
      
    }

    // Tunable parameters setters/getters
    double get_Psi () const noexcept { return _Psi; }
    void set_Psi ( const double Psi ) noexcept {
      _Psi = Psi;
      _paramsrc[ 0 ] = _Psi;
      return;
    }

    // =============================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	sfh_empirical::serialize_size() +
	SerialPOD< double >::serialize_size( _Psi );

    }

    virtual char * serialize ( char * data ) const {

      data = sfh_empirical::serialize( data );
      data = SerialPOD< double >::serialize( data, _Psi );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = sfh_empirical::deserialize( data );
      data = SerialPOD< double >::deserialize( data, _Psi );
      return data;

    }

    // =============================================================
    
  }; // endclass sfh_constant

  // ================================================================
  // Delayed Exponential SFH

  class sfh_delayedexp : public sfh_empirical {

  private :

    // tunable parameters
    double _Psi_norm = 1.;    // [ Msol * yr^-1 ]
    double _k_shape = 0.2;    // []
    double _tau_star = 1.e+8; // [ yr ]
    
    // private functions of the class
    double _sfr ( const double tau ) const noexcept override {

      return _paramsrc[ 0 ] * std::pow( tau, _paramsrc[ 1 ] ) * \
	std::exp( - tau * _paramsrc[ 2 ] );

    }

  public :

    // sfh_delayedexp () : sfh_empirical{} { _paramsrc = std::vector< double >( 3 ); }

    sfh_delayedexp ( const double tau_quench = 2.e+10 ) noexcept
      : sfh_empirical { tau_quench } {

      _paramsrc = std::vector< double >( 3 );
      set_params();
      
    }
    
    // In:
    // param[ 0 ] = Psi_norm
    // param[ 1 ] = k
    // param[ 2 ] = tau_star
    // param[ 3 ] = M_dust
    // param[ 4 ] = Z_gxy
    // Out:
    // idx_0 = Psi_norm
    // idx_1 = k
    // idx_2 = 1. / tau_star
    void set_params ( const double * const param = nullptr ) noexcept override {

      if ( param ) {
	this->set_Psi_norm( param[ 0 ] );
	this->set_k_shape( param[ 1 ] );
	this->set_tau_star( param[ 2 ] );
	this->set_Mdust( param[ 3 ] );
	this->set_Zgxy( param[ 4 ] );
      }
      else {
	this->set_Psi_norm( _Psi_norm );
	this->set_k_shape( _k_shape );
	this->set_tau_star( _tau_star );
	this->set_Mdust( _Mdust );
	this->set_Zgxy( _Zgxy );
      }
      return;
      
    }

    // Tunable parameters setters/getters
    double get_Psi_norm () const noexcept { return _Psi_norm; }
    double get_k_shape () const noexcept { return _k_shape; }
    double get_tau_star () const noexcept { return _tau_star; }
    void set_Psi_norm ( const double Psi_norm ) noexcept {
      _Psi_norm = Psi_norm;
      _paramsrc[ 0 ] = _Psi_norm;
      return;
    }
    void set_k_shape ( const double k_shape ) noexcept {
      _k_shape = k_shape;
      _paramsrc[ 1 ] = _k_shape;
      return;
    }
    void set_tau_star ( const double tau_star ) noexcept {
      _tau_star = tau_star;
      _paramsrc[ 2 ] = 1. / _tau_star;
      return;
    }

    // =============================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	sfh_empirical::serialize_size() +
	SerialPOD< double >::serialize_size( _Psi_norm ) +
	SerialPOD< double >::serialize_size( _k_shape ) +
	SerialPOD< double >::serialize_size( _tau_star );

    }

    virtual char * serialize ( char * data ) const {

      data = sfh_empirical::serialize( data );
      data = SerialPOD< double >::serialize( data, _Psi_norm );
      data = SerialPOD< double >::serialize( data, _k_shape );
      data = SerialPOD< double >::serialize( data, _tau_star );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = sfh_empirical::deserialize( data );
      data = SerialPOD< double >::deserialize( data, _Psi_norm );
      data = SerialPOD< double >::deserialize( data, _k_shape );
      data = SerialPOD< double >::deserialize( data, _tau_star );
      return data;

    }

    // =============================================================
    
  }; // endclass sfh_delayedexp

  // ================================================================
  // Log-Normal SFH

  class sfh_lognorm : public sfh_empirical {

  private :

    // tunable parameters
    double _Psi_norm = 100.;  // [ Msol * yr^-1 ]
    double _sigma_star = 2.;  // []
    double _tau_star = 1.e+8; // [ yr ]
    
    // private functions of the class
    double _sfr ( const double tau ) const noexcept override {

      double lntau = std::log( tau * _paramsrc[ 2 ] );

      return _paramsrc[ 0 ] * std::sqrt( sed::cnst::ip * _paramsrc[ 1 ] ) * \
	std::exp( - lntau * lntau * _paramsrc[ 1 ] );

    }

  public :

    // sfh_lognorm () : sfh_empirical{} { _paramsrc = std::vector< double >( 3 ); }

    sfh_lognorm ( const double tau_quench = 2.e+10 ) noexcept
      : sfh_empirical { tau_quench } {

      _paramsrc = std::vector< double >( 3 );
      set_params();
      
    }

    // In:
    // param[ 0 ] = Psi_norm
    // param[ 1 ] = sigma_star
    // param[ 2 ] = tau_star
    // param[ 3 ] = M_dust
    // param[ 4 ] = Z_gxy
    // Out:
    // idx_0 = Psi_norm
    // idx_1 = 1 / ( 2 * sigma_star^2 )
    // idx_2 = 1. / tau_star    
    void set_params ( const double * const param = nullptr ) noexcept override {

      if ( param ) {
	this->set_Psi_norm( param[ 0 ] );
	this->set_sigma_star( param[ 1 ] );
	this->set_tau_star( param[ 2 ] );
	this->set_Mdust( param[ 3 ] );
	this->set_Zgxy( param[ 4 ] );
      }
      else {
	this->set_Psi_norm( _Psi_norm );
	this->set_sigma_star( _sigma_star );
	this->set_tau_star( _tau_star );
	this->set_Mdust( _Mdust );
	this->set_Zgxy( _Zgxy );
      }
      return;
      
    }

    // Tunable parameters setters/getters
    double get_Psi_norm () const noexcept { return _Psi_norm; }
    double get_sigma_star () const noexcept { return _sigma_star; }
    double get_tau_star () const noexcept { return _tau_star; }
    void set_Psi_norm ( const double Psi_norm ) noexcept {
      _Psi_norm = Psi_norm;
      _paramsrc[ 0 ] = _Psi_norm;
      return;
    }
    void set_sigma_star ( const double sigma_star ) noexcept {
      _sigma_star = sigma_star;
      _paramsrc[ 1 ] = 1. / ( 2 * _sigma_star * _sigma_star );
      return;
    }
    void set_tau_star ( const double tau_star ) noexcept {
      _tau_star = tau_star;
      _paramsrc[ 2 ] = 1. / _tau_star;
      return;
    }

    // =============================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	sfh_empirical::serialize_size() +
	SerialPOD< double >::serialize_size( _Psi_norm ) +
	SerialPOD< double >::serialize_size( _sigma_star ) +
	SerialPOD< double >::serialize_size( _tau_star );

    }

    virtual char * serialize ( char * data ) const {

      data = sfh_empirical::serialize( data );
      data = SerialPOD< double >::serialize( data, _Psi_norm );
      data = SerialPOD< double >::serialize( data, _sigma_star );
      data = SerialPOD< double >::serialize( data, _tau_star );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = sfh_empirical::deserialize( data );
      data = SerialPOD< double >::deserialize( data, _Psi_norm );
      data = SerialPOD< double >::deserialize( data, _sigma_star );
      data = SerialPOD< double >::deserialize( data, _tau_star );
      return data;

    }

    // =============================================================
    
  }; // endclass sfh_lognorm

} // endnamespace sed

#endif //__SFH_EMPIRICAL_H__
