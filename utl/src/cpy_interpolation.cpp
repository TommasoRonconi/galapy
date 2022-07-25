#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cpy_utilities.h>
#include <cpy_serialize.h>
#include <interpolation.h>
#include <iostream>

template< class T >
struct _PyInterp {
  PyObject_HEAD /* macro for fixed-size structure initialization */
  utl::interpolator< T > * ptrObj;
};

// initialize _PyInterp Object
template< class T, class U >
int _PyInterp_init ( T * self, PyObject * args, PyObject * kwds ) {

  // input arrays
  PyArrayObject * NPyXV = NULL;
  PyArrayObject * NPyFV = NULL;
  std::vector< double > xv, fv;
  
  // interpolation type
  char const * type = "linear";
    
  // work-around to silence compile warning (const char* -> char* forbidden in C++):
  char ** kwlist = new char * [ 4 ] { (char*)"xv", (char*)"fv", (char*)"kind", NULL };
    
  // Get the passed input arrays and interpolation type
  if ( !PyArg_ParseTupleAndKeywords( args, kwds, "O!O!|s", kwlist,
				     &PyArray_Type, &NPyXV,
				     &PyArray_Type, &NPyFV,
				     &type ) ) {
    return -1;
  }

  // check for errors in passed interpolation type
  // if ( !( ( strcoll( type, "linear" )  == 0 ) ||
  // 	  ( strcoll( type, "cspline" ) == 0 ) ) ) {
  if ( !( ( strcoll( type, "linear" )  == 0 ) ) ) {
    // PyErr_SetString( PyExc_TypeError,
    // 		     "Interpolation type not valid. "
    // 		     "Valid types are: 'linear', 'cspline'" );
    PyErr_SetString( PyExc_TypeError,
		     "Interpolation type not valid. "
		     "Valid types are: 'linear'" );
    return -1;
  }

  /* Convert NumPy-arrays to C++ vectors */
  if ( NPyArrayToCxxVector1D< double >( NPyXV, xv ) == -1 ) return -1;
  if ( NPyArrayToCxxVector1D< double >( NPyFV, fv ) == -1 ) return -1;
  
  // Allocate c++ templated interpolator object
  self->ptrObj = new utl::interpolator< U >{ xv, fv, std::string{ type } };

  return 0;
    
}

template< class T >
void _PyInterp_dealloc ( T * self ) {

  delete self->ptrObj;
  Py_TYPE( self )->tp_free( self );
  return;
    
}

template< class T >
PyObject *  _PyInterp_call ( T * self, PyObject * args ) {

  PyObject * NPyBuf;

  /* Parse arguments */
  if ( !PyArg_ParseTuple( args, "O", &NPyBuf ) ) return NULL;

  /* Scalar input */
  if ( PyObject_TypeCheck( NPyBuf, &PyFloat_Type ) ) {
    double xx = PyFloat_AS_DOUBLE( NPyBuf );
    return PyFloat_FromDouble( ( *( self->ptrObj ) )( xx ) );
  }
  /* Array-like input */
  else {
    
    /* Convert Numpy-array to C-array */
    double * xx;
    std::size_t size;
    if ( NPyArrayToCArray1D< double >( (PyArrayObject*)NPyBuf, &xx, &size ) == -1 ) return NULL;

    /* Call C++ member function */
    double * outarr = new double [ size ];
    for ( unsigned int ii = 0; ii < size; ++ii )
      outarr[ ii ] = ( *( self->ptrObj ) )( xx[ ii ] );

    /* Clear heap */
    delete [] xx;

    PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
						reinterpret_cast< void * >( outarr ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    return ret;

  }
     
}

template< class T >
PyObject *  _PyInterp_integrate ( T * self, PyObject * args ) {

  double aa, bb;
  
  if ( !PyArg_ParseTuple( args, "dd", &aa, &bb ) ) {
    return NULL;
  }

  return PyFloat_FromDouble( ( self->ptrObj )->integrate( aa, bb ) );

}

extern "C" {

  // ========================================================================================
  // ========================================= LIN ==========================================
  // ========================================================================================

  // define the linear interpolator type
  typedef struct _PyInterp< utl::lin_interp > PyLinInterp;
  
  // ========================================================================================
  
  // initialize PyLinInterp Object
  static int PyLinInterp_init ( PyLinInterp * self, PyObject * args, PyObject * kwds ) {
    
    return _PyInterp_init< PyLinInterp, utl::lin_interp >( self, args, kwds );

  }
  
  // ========================================================================================
  
  // deallocate PyLinInterp Object
  static void PyLinInterp_dealloc ( PyLinInterp * self ) {

    _PyInterp_dealloc< PyLinInterp >( self );
    return;

  }
  
  // ========================================================================================

  static PyObject * PyLinInterp_call ( PyLinInterp * self, PyObject * args ) {

    return _PyInterp_call< PyLinInterp >( self, args );

  }
  
  // ========================================================================================
  
  static const char DocString_integrate[] =
    "Function for integrating the function in a range.\n"
    "\nParameters"
    "\n----------"
    "\naa : float\n"
    "\tlower limit for integration\n"
    "\nbb : float\n"
    "\tupper limit for integration\n"
    "\nReturns"
    "\n-------"
    "\n: float\n"
    "\tAn approximated integral of the given function\n"; 
  static PyObject * PyLinInterp_integrate ( PyLinInterp * self, PyObject * args ) {

    return _PyInterp_integrate< PyLinInterp >( self, args );

  }
  
  // ========================================================================================

  /* Pickle the object */
  static PyObject * PyLinInterp___getstate__ ( PyLinInterp * self, PyObject * Py_UNUSED(ignored) ) {

    PyObject * ret = CPy___getstate__< PyLinInterp >( self, NULL );
    if ( !ret ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to get state from CPySYN object" );
      return NULL;
    }

    return ret;
    
  }
  
  // ========================================================================================

  /* Un-Pickle the object */
  static PyObject * PyLinInterp___setstate__ ( PyLinInterp * self, PyObject * state ) {

    if ( !CPy___setstate__< PyLinInterp, utl::interpolator< utl::lin_interp > >( self, state ) ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to set state from PyLinInterp object" );
      return NULL;
    }
    Py_RETURN_NONE;
    
  }
  
  // ========================================================================================
    
  static PyMethodDef PyLinInterp_Methods[] =
    {
     { "integrate",
       (PyCFunction)PyLinInterp_integrate,
       METH_VARARGS,
       DocString_integrate },
     { "__getstate__",
       (PyCFunction) PyLinInterp___getstate__,
       METH_NOARGS,
       "Pickle the Custom object" },
     { "__setstate__",
       (PyCFunction) PyLinInterp___setstate__,
       METH_O,
       "Un-pickle the Custom object" },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject PyLinInterp_t = { PyVarObject_HEAD_INIT( NULL, 0 )
					"galapy.internal.interp.lin_interp"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef interp_module = {
					     PyModuleDef_HEAD_INIT,
					     "galapy.internal.interp",
					     "",
					     -1,
					     NULL, NULL, NULL, NULL, NULL
  }; /* endPyModuleDef interp_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_interp( void ) {

    /* -------------------- */
    /* Create interp module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &interp_module );
    if ( m == NULL )
        return NULL;

    /* --------------------------------------------- */
    /* Adding new object lin_interp to interp module */
    /* --------------------------------------------- */
    PyLinInterp_t.tp_new       = PyType_GenericNew;
    PyLinInterp_t.tp_basicsize = sizeof( PyLinInterp );
    PyLinInterp_t.tp_dealloc   = (destructor) PyLinInterp_dealloc;
    PyLinInterp_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    PyLinInterp_t.tp_doc       = "Linear interpolator object";
    PyLinInterp_t.tp_call      = (ternaryfunc)PyLinInterp_call;
    PyLinInterp_t.tp_methods   = PyLinInterp_Methods;
    PyLinInterp_t.tp_init      = (initproc) PyLinInterp_init;
    if ( PyType_Ready( &PyLinInterp_t ) < 0 )
      return NULL;
    Py_INCREF( &PyLinInterp_t );
    PyModule_AddObject( m, "lin_interp", (PyObject *)&PyLinInterp_t );

    /* -------------------------------- */
    /* Initialize NumPy Array interface */
    /* -------------------------------- */
    import_array();

    /* return new module */
    return m;
    
  }
  
} // endextern "C"

// ========================================================================================
