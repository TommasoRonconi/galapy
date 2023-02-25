/**
 *  @file ibstree/iterator.h
 *
 *  @brief The classes iterator and const_iterator
 *
 *  This file defines the interface of the class iterator and 
 *  of class const_iterator
 *
 *  @author Herbert Nguruwe, Tommaso Ronconi
 *
 *  @author hknguruwe@gmail.com, tronconi@sissa.it
 */


#ifndef __ITERATOR__
#define __ITERATOR__

// STL includes
#include <memory>

// internal includes
#include "node.h"

namespace utl {

  /**
   *  @class iterator iterator.h "ibstree/iterator.h"
   *
   *  @brief The class iterator
   *
   *  This class is used to handle objects of type <EM> iterator
   *  </EM>. It is templated on two types T and U for the  
   *  <EM> key </EM> and <EM> value </EM>, respectively, of the node
   *  contained by the iterator is .
   */
  template < class T, class U >
  class iterator {

    /// Type node keyword definition
    using node = struct node<T, U>;

    /// Actual content of the iterator: a pointer to some node
    node * current;
  
  public:

    /// Default constructor
    iterator () noexcept = default;

    /**
     *  @brief Custom constructor of iterator to some node n
     *
     *  @param n raw pointer to the node to which the iterator will iterator
     *           (iterators gonna iterate..lol)
     */ 
    iterator( node * n ) : current{ n } {}

    /// De-reference operator overload, it returns a reference to the de-referenced object of
    /// type node contained in the iterator
    node& operator*() { return *current; }

    /// Member-access operator overload, it returns the raw pointer to the object of type node
    /// contained in the iterator
    node* operator->() const { return current;  }
  
    /**
     *  @brief Pre-increment operator overload, it copies into <b>*this</b> a pointer to 
     *         the next element in the Tree hierarchy.
     *
     *  @return A reference to the resulting iterator
     */
    iterator& operator++() {

      if ( current ) {
	if ( current->right )
	  current = current->right->leftmost() ;
	else
	  current = current->parent;
      }
      
      return *this;
    
    }
  
    /**
     *  @brief Post-increment operator overload, implementation calls overloaded pre-increment 
     *         operator 
     *
     *  @return The resulting iterator
     */
    iterator operator++( int ) {
    
      iterator it{*this};
    
      ++(*this);
    
      return it;
    
    }
  
    /**
     *  @brief Logical-equality operator overload
     *
     *  @param other the r-value to be compared
     *
     *  @return bool, true if the two iterator contain the same pointer, false if not
     */
    bool operator==(const iterator& other) { return current == other.current; }
  
    /**
     *  @brief Logical-inequality operator overload
     *
     *  @param other the r-value to be compared
     *
     *  @return Implementation is done in terms of overloaded equality operator
     */
    bool operator!=(const iterator& other) { return !( *this == other ); }
  
  }; // end of class iterator


  // ===========================================================================


  /**
   *  @class const_iterator const_iterator.h "ibstree/iterator.h"
   *
   *  @brief The class const_iterator
   *
   *  This class is used to handle objects of type <EM> const_iterator </EM>.
   *  It inherits from class iterator and is templated on two types T and U for  
   *  the <EM> key </EM> and <EM> value </EM>, respectively, of the node
   *  contained by the iterator is .
   */
  template < class T, class U >
  class const_iterator : public iterator<T, U> {
  
    /// Type node keyword definition
    using node = struct node<T, U>;
    
  public:
  
    /// 'parent' keyword definition, aliases the iterator class 
    using parent = iterator<T, U>;
  
    using parent::iterator;

    /// De-reference operator overload (constant version), it returns a reference
    /// to the de-referenced object of type node contained in the iterator
    const node* operator->() const { return parent::operator->();  }
  
    /// Member-access operator overload (constant version), it returns the
    /// raw pointer to the object of type node contained in the iterator
    const node& operator*() const { return parent::operator*(); }
  
    using parent::operator==;
  
    using parent::operator!=;
  
  }; // end of class const_iterator

} // endnamespace utl
  
#endif //__ITERATOR__
