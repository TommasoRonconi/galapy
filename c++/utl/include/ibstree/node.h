/**
 *  @file ibstree/node.h
 *
 *  @brief The class node
 *
 *  This file defines the interface of the class node
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */


#ifndef __NODE__
#define __NODE__

// STL includes
#include <memory>
#include <utility>

// internal includes
#include "interval.h"

namespace utl {

  /**
   *  @class node node.h "ibstree/node.h"
   *
   *  @brief The class node
   *
   *  This class is used to handle objects of type <EM> node
   *  </EM>. It is templated on two types T and U for the  
   *  <EM> key </EM> and <EM> value </EM> retained by the node, respectively.
   */
  template< class T, class U >
  struct node {

    using IT = interval< T >;

    /// Content of the node, a std::pair templated on the key (first
    /// argument of pair) and value (second argument of std::pair) types
    std::pair< IT, U > content;
    
    /// Raw pointer to parent node
    node * parent = nullptr;

    /// Pointer to left node (std::unique_ptr)
    std::unique_ptr<node> left = nullptr;

    /// Pointer to right node (std::unique_ptr)
    std::unique_ptr<node> right = nullptr;

    /**
     *  @name Constructors/Destructor
     */
    ///@{

    /// default constructor
    node () : content{ std::pair< IT, U>{} }, parent{ nullptr } {}

    /**
     *  @brief Constructor that builds a new node setting the content and position wrt parent
     *
     *  @param key the key of the node
     *
     *  @param value the value contained in the node
     *
     *  @param par parent node to which the node is attached (default = nullptr)
     */
    node ( IT key, U value, node * par = nullptr )
      : content{ std::pair<IT,U>( key, value ) }, parent{ par } {}

  
    // ==================================================
  

    /// default destructor of node
    ~node() = default;

    ///@}

    /**
     *  @name Public functions of the struct node
     */
    ///@{

    /// returns the key of the node
    const IT& key() const { return content.first; }

    /// returns the value of the node
    const U& value() const { return content.second; }

    /// returns the value of the node
    U& value() { return content.second; }

    node * find ( const T key );

    /**
     *  @brief Recursive function to insert a new node lower in hierarchy with respect
     *         to current node. If node has no childs it inserts the new one as the 
     *         correct new child, otherwise calls itself again from the right child of
     *         current node.
     *
     *  @param key key of the new node to be generated
     *  
     *  @param value value of the new node to be generated
     *
     *  @param substitute whether to substitute or not the value if a node with same 
     *                    key  exists
     *
     *  @return raw pointer to last node inserted
     */
    node * insert ( const IT key, const U value );


    /**
     *  @brief Recursive function to extract the current node and all the nodes 
     *         lower in hierarchy with respect to current.
     *         First it stores a raw pointer to the current node than calls itself again 
     *         from the childs, first left and then right.
     *         Finally, it returns void.
     *
     *  @param (out) a vector of raw constant pointers of type node<T,U>
     *
     *  @return void
     */
    void extract ( std::vector< const node * > & store ) const noexcept;

    /**
     *  @brief Recursive function to find the leftmost node in hierarchy from current
     *
     *  @return Raw pointer to leftmost node, if current has no left child returns <b>this</b>
     */
    node * leftmost () {

      if ( left )
	return left->leftmost();
      else
	return this;
   
    }

    /**
     *  @brief Recursive function for removing all nodes lower in hierarchy than current.
     *         If current has left/right child calls itself again from left/right 
     *         than resets left/right.
     *
     *  @return void
     */
    void clear () {

      // reset everything on the left
      left.reset();

      // reset everything on the right
      right.reset();
    
    }


    ///@}
    
  }; // end of class node

#include "node.tpp"

} // endnamespace utl

#endif //__NODE__
