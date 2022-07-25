/**
 *  @file ibstree/ibstree.h
 *
 *  @brief The class ibstree
 *
 *  This file defines the interface of the class ibstree
 *
 *  @author Herbert Nguruwe, Tommaso Ronconi
 *
 *  @author hknguruwe@gmail.com, tronconi@sissa.it
 */


#ifndef __TREE__
#define __TREE__

// STL includes
#include <iostream>
#include <memory>
#include <vector>
#include <functional>
#include <cmath>

// internal includes
#include "node.h"
#include "iterator.h"

namespace utl {

  /**
   *  @class ibstree ibstree.h "ibstree/ibstree.h"
   *
   *  @brief The class ibstree
   *
   *  This class is used to handle objects of type <EM> ibstree
   *  </EM>. It is templated on two types T and U for the  
   *  <EM> key </EM> and <EM> value </EM> of the tree, respectively.
   */
  template < class T, class U >
  class ibstree {
  
    /// structure 'node<T,U>' to 'node'
    using node = struct node< T, U >;

    /**
     *  @name Private variables of the class
     */
    ///@{

    /// The root node of the BST
    std::unique_ptr<node> root = nullptr;

    /// The tail of the BST
    node * tail = nullptr;

    ///@}

    /**
     *  @name Private functions of the class
     */
    ///@{


    /**
     *  @brief Templeted private function for in-place balance a section of the ibstree. 
     *         It receives a vector of raw pointers to node and recursiverly calls itself
     *         setting each branch of the ibstree with a top-down strategy.
     *         The implemented algorithm recursively splits the vector in 2 parts, inserts the
     *         mid-point into the current position and calls itself again once for the first half
     *         and once for the second half of the remaining nodes.
     *         When the vector is empty it sets left and right to nullptr and returns.
     *
     *  @param here iterator to current position
     *
     *  @param nodes vector of raw pointers to nodes
     *
     *  @return void
     */
    void kernel_balance ( iterator< T, U > here, const std::vector<node*>& nodes ); 

    ///@}
  
  public:

    /// class 'iterator< T, U >' to 'iterator'
    using iterator = class iterator< T, U >;

    /// class 'const_iterator< T, U >' to 'const_iterator'
    using const_iterator = class const_iterator< T, U >;

    /**
     *  @name Friends of the class ibstree
     */
    ///@{

    /// operator<< overload
    template < class ot, class ou >
    friend std::ostream& operator<< ( std::ostream&, const ibstree< ot, ou >& );
  
    ///@}

    /**
     *  @name Constructors/Destructor
     */
    ///@{


    /// default constructor
    ibstree () = default;
  
    /// Copy-constructor builds new ibstree by deep-coping argument ibstree 
    ibstree ( const ibstree & T_other );

    /// Copy-assignment operator overload, assign to existing ibstree
    /// a deep copy of some other ibstree
    ibstree& operator= ( const ibstree & T_other ) {

      this->clear();
      auto tmp = T_other;
      (*this) = std::move( tmp );

      return *this;
    
    }

    /// move-constructor
    ibstree ( ibstree&& T_other ) : root{ std::move( T_other.root ) },
				    tail{ std::move( T_other.tail ) } {}

    /// move-assignment operator
    ibstree& operator= ( ibstree&& T_other ) {

      root = std::move( T_other.root );
      tail = std::move( T_other.tail );

      return *this;

    }

    // default destructor
    ~ibstree() noexcept = default;

    ///@}

    /**
     *  @name iterators declaration
     */
    ///@{

    /// iterator to the tail of the tree (aka the node with smallest key)
    iterator begin() { return iterator{ tail }; }

    /// iterator to nullptr (aka the parent of the node with the largest key)
    iterator end() { return iterator{ nullptr }; }

    /// iterator to the root of the tree (aka the highest level node)
    iterator top() { return iterator{ root.get() }; }

    /// const_iterator to the tail of the tree (aka the node with smallest key)
    const_iterator cbegin() const { return const_iterator{ tail }; }

    /// const_iterator to nullptr (aka the parent of the node with the largest key)
    const_iterator cend() const { return const_iterator{ nullptr }; }

    /// const_iterator to the root of the tree (aka the highest level node)
    const_iterator ctop() const { return const_iterator{ root.get() }; }

    ///@}

    /**
     *  @name Public functions of the class ibstree
     */
    ///@{

    bool planted () const noexcept {
      
      if ( root ) return true;
      return false;
      
    }

    /**
     *  @brief Templated function to insert a new node in the ibstree.
     *         If the root already exists it calls function node::insert, otherwise
     *         it inserts root.
     *
     *  @param key the key that will identify the new node
     *
     *  @param value the value contained by the new node
     *
     *  @param substitute <b>true</b>: substitute the value in node if key is present<br/>
     *                    <b>false</b>: do not substitute and exit function
     *
     *  @return iterator to inserted node
     */  
    iterator insert ( const interval<T> key, const U value );


    /**
     *  @brief Templated function to extract a vector of raw pointers to the nodes in the ibstree.
     *         If the root exists it calls the recursive function node::insert. 
     *
     *  @param (out) reference to the allocated vector used to store the raw pointers
     *
     *  @return void
     */  
    void extract ( std::vector< const node * > & store ) const noexcept {

      if ( root )
	root->extract( store );
      return;
      
    }

    /**
     *  @brief Function to remove all nodes from ibstree. 
     *         Resets <b>root</b> unique pointer.
     *
     *  @return void
     */
    void clear () {
    
      root.reset();
    
    }

    /**
     *  @brief Function to balance the ibstree.
     *         <ol>
     *         <li> Fills an ordered vector with raw pointers to the nodes composing the ibstree;</li>
     *         <li> Resets root node with midpoint of said vector;</li>
     *         <li> Call to private function ibstree::kernel_balance() with arguments:
     *                   <ul>
     *                   <li><EM>iterator{ root.get() }</EM>: iterator to root</li>
     *                   <li><EM>nodes</EM>: ordered vector of raw pointers to node</li>
     *                   </ul></li>
     *         </ol>
     *         If nodes vector is empty resets root to nullptr.
     *
     *  @return void
     */
    void balance ();

    /**
     *  @brief Function to find an element with given key.
     *         It calls the recursive function find() of the root node
     *
     *  @param key constant key value to be searched
     *
     *  @return iterator to position in ibstree where key is found, if key is not found
     *  returns iterator to end of ibstree (i.e. nullptr)
     *
     */
    iterator find ( const T key ) { return iterator { root->find( key ) }; }

    /**
     *  @brief Function to find an element with given key.
     *         It calls the recursive function find() of the root node
     *
     *  @param key constant key value to be searched
     *
     *  @return const_iterator to position in ibstree where key is found, if key is not found
     *  returns const_iterator to end of ibstree (i.e. nullptr)
     *
     */
    const_iterator find ( const T key ) const { return const_iterator { root->find( key ) }; }

    /** 
     *  @brief Exception handler for key not found
     */
    struct key_not_found {
    
      std::string message;
      key_not_found( const std::string &s ) : message{s} {}
    
    };

    /**
     *  @brief Array subscript operator overload, const version
     *
     *  @param key the key of the referred node
     *
     *  @return returns a reference to the value contained in the node.
     *  If the node does not exist it throws an exception key_not_found.
     */
    const U & operator[]( const T key ) const {
    
      const_iterator n = find( key );
    
      if ( n == cend() )
	throw key_not_found{ "\nI'm sorry dave, I'm afraid I can't do that\n(key " +
	    std::to_string(key) +
	    " is out of bounds.)"};
    
      return n->value();
    
    }

    ///@}
  
  
  }; // end of class ibstree declaration

#include "ibstree.tpp"

} // endnamespace utl

#endif //__TREE__
