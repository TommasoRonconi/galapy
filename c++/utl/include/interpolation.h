#ifndef __INTERPOLATOR__
#define __INTERPOLATOR__

/// internal includes
#include <utilities.h>
#include <interp/base_interface.h>
#include <interp/ibstree_interface.h>

namespace utl {

  template< class T = lin_interp >
  class interpolator {

  private:

    T _interface;

  public:

    interpolator () = default;

    interpolator ( const T & interface )
      : _interface{ interface } {}

    // generic forwarding constructor
    template< class ... Args >
    interpolator ( Args && ... args )
      : _interface{ std::forward< Args >( args )... } {}

    // destructor
    virtual ~interpolator () = default;

    double operator() ( const double xx ) const noexcept {
  
      return _interface.eval( xx );
  
    }

    // move constructor
    interpolator ( interpolator && ii ) = default;

    // copy contructor
    interpolator ( const interpolator & ii ) = default;

    // move assignment
    interpolator & operator= ( interpolator && ii ) = default;

    // copy assignment
    interpolator & operator= ( interpolator & ii ) = default;

    double integrate ( const double aa, const double bb ) const noexcept {

      return _interface.integrate( aa, bb );

    }

    size_t get_thinness () const noexcept { return _interface.get_thinness(); }
      
    double get_xmin () const noexcept { return _interface.get_xmin(); }
      
    double get_xmax () const noexcept { return _interface.get_xmax(); }
      
    std::vector< double > get_xv () const noexcept { return _interface.get_xv(); }

    std::vector< double > get_fv () const noexcept { return _interface.get_fv(); }

    size_t size () const noexcept { return _interface.size(); }

    interpolator & operator+= ( const interpolator & rhs ) {

      _interface += rhs._interface;
	
      return * this;

    }
      
    friend interpolator operator+ ( interpolator lhs,
				    const interpolator & rhs ) {

      lhs += rhs;
      return lhs;

    }

    interpolator & operator-= ( const interpolator & rhs ) {

      _interface -= rhs._interface;
	
      return * this;

    }
      
    friend interpolator operator- ( interpolator lhs,
				    const interpolator & rhs ) {

      lhs -= rhs;
      return lhs;

    }

    interpolator & operator*= ( const interpolator & rhs ) {

      _interface *= rhs._interface;
	
      return * this;

    }
      
    friend interpolator operator* ( interpolator lhs,
				    const interpolator & rhs ) {

      lhs *= rhs;
      return lhs;

    }

    interpolator & operator/= ( const interpolator & rhs ) {

      _interface /= rhs._interface;
	
      return * this;

    }
      
    friend interpolator operator/ ( interpolator lhs,
				    const interpolator & rhs ) {

      lhs /= rhs;
      return lhs;

    }

    interpolator & operator+= ( const double & rhs ) {

      _interface += rhs;
	
      return * this;

    }
      
    friend interpolator operator+ ( interpolator lhs,
				    const double & rhs ) {

      lhs += rhs;
      return lhs;

    }
      
    friend interpolator operator+ ( const double & lhs,
				    interpolator rhs ) {

      rhs += lhs;
      return rhs;

    }

    interpolator & operator-= ( const double & rhs ) {

      _interface += -rhs;
	
      return * this;

    }
      
    friend interpolator operator- ( interpolator lhs,
				    const double & rhs ) {

      lhs -= rhs;
      return lhs;

    }
      
    friend interpolator operator- ( const double & lhs,
				    interpolator rhs ) {

      rhs -= lhs;
      rhs *= -1;
      return lhs;

    }

    interpolator & operator*= ( const double & rhs ) {

      _interface *= rhs;
	
      return * this;

    }
      
    friend interpolator operator* ( interpolator lhs,
				    const double & rhs ) {

      lhs *= rhs;
      return lhs;

    }
      
    friend interpolator operator* ( const double & lhs,
				    interpolator rhs ) {

      rhs *= lhs;
      return lhs;

    }

    interpolator & operator/= ( const double & rhs ) {

      _interface *= 1. / rhs;
	
      return * this;

    }
      
    friend interpolator operator/ ( interpolator lhs,
				    const double & rhs ) {

      lhs *= 1. / rhs;
      return lhs;

    }
      
    friend interpolator operator/ ( const double & lhs,
				    interpolator rhs ) {

      return interpolator{ operator/( lhs, rhs._interface ) };

    }

    // =============================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return _interface.serialize_size();

    }

    virtual char * serialize ( char * data ) const {

      return _interface.serialize( data );

    }

    virtual const char * deserialize ( const char * data ) {

      return _interface.deserialize( data );

    }

    // =============================================================================

  }; //endclass interpolator

} //endnamespace utl

#endif //__INTERPOLATOR__
