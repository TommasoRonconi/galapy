/**
 *  @file sfh/include/sfh.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __SFH_BASE_H__
#define __SFH_BASE_H__

// external includes
#include <vector>
#include <cmath>

// internal includes
#include <utilities.h>
#include <serialize.h>
#include <imf.h>

namespace sed {

  // =============================================================================
  // SFH base class

  class sfh_base : public Serializable {

  private :

    // private variables
    double _tau_quench;
    std::vector< double > _dM = std::vector< double >( 2e5, 0.0 );

  protected :
    
    std::vector< double > _paramsrc;

    virtual double _sfr ( const double tau ) const noexcept = 0;
    
  public :

    sfh_base () = default;

    sfh_base ( const double tau_q ) noexcept
      : _tau_quench { tau_q } {}

    virtual ~sfh_base () = default;

    virtual void set_params ( const double * const param ) noexcept = 0;
    void set_tau_quench ( const double tau_q ) { _tau_quench = tau_q; };
    virtual void set_interpolator ( __attribute__((unused)) const std::vector< double > &,
				    __attribute__((unused)) const std::vector< double > & )
      noexcept { return; };
    
    double get_tau_quench () { return _tau_quench; };
    
    std::vector< double > get_params () { return _paramsrc; }
    
    double operator() ( const double xx ) const noexcept;

    void model ( const double * const tau,
		 double * const psi,
		 const std::size_t size ) const noexcept;

    void eval ( const double * const tau,
	        double * const psi,
		const std::size_t size,
		const double * const param = nullptr ) noexcept;

    double get_Mstar ( const double tau,
    		       const std::size_t npoints = 100 ) const noexcept;

    std::vector< double > get_dMstar () const noexcept { return _dM; }
    
    // Virtual galactic content functions
    virtual double get_Mdust ( const double tau ) const noexcept = 0;
    virtual double get_Mgas  ( const double tau ) const noexcept = 0;
    virtual double get_Zgas  ( const double tau ) const noexcept = 0;
    virtual double get_Zstar ( const double tau ) const noexcept = 0;

    // Compute on time grid
    // on PyBind11
    void time_grid ( const double age,
		     const std::vector< double > & tgrid,
		     const std::vector< double > & Zgrid,
		     std::vector< double > & out_dMgrid,
		     std::vector< double > & out_Zgrid,
		     std::vector< std::size_t > & out_Zidx,
		     std::size_t & out_last_idx );
    // on CPython
    // void time_grid ( const double age,
    // 		     const double * const tgrid,
    // 		     const std::size_t tgrid_size,
    // 		     const double * const Zgrid,
    // 		     const std::size_t Zgrid_size,
    // 		     double * const * const out_psigrid,
    // 		     double * const * const out_Zgrid,
    // 		     std::size_t * const * const out_Zidx,
    // 		     std::size_t * const out_last_idx );

    // =============================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	SerialVecPOD< double >::serialize_size( _paramsrc ) +
	SerialPOD< double >::serialize_size( _tau_quench );

    }

    virtual char * serialize ( char * data ) const {

      data = SerialVecPOD< double >::serialize( data, _paramsrc );
      data = SerialPOD< double >::serialize( data, _tau_quench );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = SerialVecPOD< double >::deserialize( data, _paramsrc );
      data = SerialPOD< double >::deserialize( data, _tau_quench );
      return data;

    }

    // =============================================================

  }; // endclass sfh_base
  
} // endnamespace sed

#endif //__SFH_BASE_H__
