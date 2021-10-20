#ifndef __INTERVAL__
#define __INTERVAL__

#include <utility>
#include <iostream>

namespace utl {

  template< class T >
  class interval {

  private:
    std::pair< T, T > _limits;
  
  public:

    /// operator put-stream overload
    friend std::ostream& operator<< ( std::ostream& os, const interval< T >& in ) {

      os << "[ " << in.low() << ", " << in.upp() << " )";
      return os;

    }
  
    interval () = default;
    interval ( const T & lower, const T & upper ) : _limits{lower,upper} {}
    ~interval () = default;
    T low () const { return _limits.first; }
    T upp () const { return _limits.second; }

  }; // end of class interval

  // lhs is lower than the interval if lhs < rhs.low()
  template< class T, class U >
  inline bool operator< ( const T& lhs, const interval<U>& rhs )
  { return lhs < rhs.low(); }

  // lhs is greater than the interval if lhs >= rhs.upp()
  template< class T, class U >
  inline bool operator> ( const T& lhs, const interval<U>& rhs )
  { return !(lhs < rhs.upp()); }

  // rhs is larger than the interval if rhs >= lhs.upp()
  template< class T, class U >
  inline bool operator< ( const interval<T>& lhs, const U& rhs )
  { return rhs > lhs; }

  // rhs is smaller than the interval if rhs < rhs.low()
  template< class T, class U >
  inline bool operator> ( const interval<T>& lhs, const U& rhs )
  { return rhs < lhs; }

  // lhs is in interval if rhs.low() <= lhs < rhs.upp()
  template< class T, class U >
  inline bool operator== ( const T& lhs, const interval<U>& rhs )
  { return !( lhs < rhs.low() ) && ( lhs < rhs.upp() ); }

  // rhs is in interval if lhs.low() <= rhs < lhs.upp()
  template< class T, class U >
  inline bool operator== ( const interval<T>& lhs, const U& rhs )
  { return rhs == lhs; }

  // lhs is an interval smaller than rhs if lhs.upp() <= rhs.low()
  template< class T, class U >
  inline bool operator< ( const interval<T>& lhs,
			  const interval<U>& rhs )
  { return !( lhs.upp() > rhs.low() ); }

  // lhs is an interval bigger than rhs if lhs.low() >= rhs.upp()
  template< class T, class U >
  inline bool operator> ( const interval<T>& lhs,
			  const interval<U>& rhs )
  { return rhs < lhs; }

} // endnamespace utl

#endif //__INTERVAL__
