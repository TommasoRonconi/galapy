#ifndef __BASE_INTERFACE__
#define __BASE_INTERFACE__

#include <utilities.h>

namespace utl {

  class base_interface {

  private:

    double _x_min;

    double _x_max;

  protected:

    size_t _thinness;
      
    std::vector< double > _xv;
      
    std::vector< double > _fv;

  public:

    base_interface () = default;

    base_interface ( const double x_min, const double x_max,
		     const size_t thinness )
      : _x_min{ x_min }, _x_max{ x_max }, _thinness{ thinness } {

	_xv.reserve( _thinness );
	_fv.reserve( _thinness );

      }
    
    /// move constructor
    base_interface ( base_interface && ii )
      : _x_min{ std::move( ii._x_min ) }, _x_max{ std::move( ii._x_max ) },
	_thinness{ std::move( ii._thinness ) },
	_xv{ std::move( ii._xv ) }, _fv{ std::move( ii._fv ) } {}

    /// copy constructor
    base_interface ( const base_interface & ii )
      : _x_min{ ii._x_min }, _x_max{ ii._x_max }, _thinness{ ii._thinness } {

	_xv = ii._xv;
	_fv = ii._fv;

      }

    void swap ( base_interface & ii ) noexcept {
	
      // enable ADL (not always necessary, but good practice)
      using std::swap;

      swap( this->_x_min, ii._x_min );
      swap( this->_x_max, ii._x_max );
      swap( this->_thinness, ii._thinness );
      swap( this->_xv, ii._xv );
      swap( this->_fv, ii._fv );

      return;

    }

    virtual ~base_interface () = default;

    virtual double eval ( const double xx ) const = 0;

    virtual double integrate ( const double aa, const double bb ) const = 0;

    virtual size_t get_thinness () const { return _thinness; }

    virtual double get_xmin () const { return _x_min; }

    virtual double get_xmax () const { return _x_max; }

    virtual std::vector< double > get_xv () const { return _xv; }

    virtual std::vector< double > get_fv () const { return _fv; }

    virtual size_t size () const { return _thinness; }
      
  }; // endclass base_interface

} // endnamespace utl

#endif //__BASE_INTERFACE__
