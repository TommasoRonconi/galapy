#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cpy_utilities.h>
#include <interpolation.h>

template< class T >
struct _PyInterp {
  PyObject_HEAD /* macro for fixed-size structure initialization */
  utl::interpolator< T > * ptrObj;
};

// initialize _PyInterp Object
template< class T, class U >
int _PyInterp_init ( T * self, PyObject * args, PyObject * kwds ) {

  // input arrays
  PyObject *xv_buf, *fv_buf;
  Py_buffer xv, fv;

  // interpolation type
  char const * type = "linear";
    
  // work-around to silence compile warning (const char* -> char* forbidden in C++):
  char ** kwlist = new char * [ 4 ] { (char*)"xv", (char*)"fv", (char*)"kind", NULL };
    
  // Get the passed input arrays and interpolation type
  if ( !PyArg_ParseTupleAndKeywords( args, kwds, "OO|s", kwlist,
				     &xv_buf, &fv_buf, &type ) ) {
    return -1;
  }
  
  // check for errors in passed interpolation type
  if ( !( ( strcoll( type, "linear" )  == 0 ) ||
	  ( strcoll( type, "cspline" ) == 0 ) ) ) {
    PyErr_SetString( PyExc_TypeError,
		     "Interpolation type not valid. "
		     "Valid types are: 'linear', 'cspline'" );
    PyBuffer_Release( &xv );
    PyBuffer_Release( &fv );
    return -1;
  }

  // convert input x-axis ndarray to buffer
  if ( return_1D_array( xv_buf, &xv, "d") == -1 ) {
    return -1;
  }

  // convert input y-axis ndarray to buffer
  if ( return_1D_array( fv_buf, &fv, "d") == -1 ) {
    return -1;
  }

  // Allocate c++ templated interpolator object
  self->ptrObj =
    new utl::interpolator< U >
    { std::vector< double >( (double*)xv.buf,
			     (double*)xv.buf + xv.shape[ 0 ] ),
      std::vector< double >( (double*)fv.buf,
			     (double*)fv.buf + fv.shape[ 0 ] ),
      std::string{ type } };

  // Indicate we're done working with the buffer
  PyBuffer_Release( &xv );
  PyBuffer_Release( &fv );
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

  PyObject * xx_buf;
  
  if ( !PyArg_ParseTuple( args, "O", &xx_buf ) ) {
    return NULL;
  }
  if ( PyObject_TypeCheck( xx_buf, &PyFloat_Type ) ) {
    double xx = PyFloat_AS_DOUBLE( xx_buf );
    return PyFloat_FromDouble( ( *( self->ptrObj ) )( xx ) );
  }
  else {
    Py_buffer xx;
    if ( return_1D_array( xx_buf, &xx, "d" ) == -1 ) {
      return NULL;
    }
    PyObject * pyList = PyList_New( xx.shape[ 0 ] );
    std::vector< double > xv ( (double*)xx.buf, (double*)xx.buf + xx.shape[ 0 ] );
    PyBuffer_Release( &xx );
    for ( std::size_t ii = 0; ii < xv.size(); ++ii )
      if ( PyList_SetItem( pyList, ii,
			   PyFloat_FromDouble( ( *( self->ptrObj ) )( xv[ ii ] ) )
			   ) != 0 ) return NULL;
    return pyList;
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
  
  // initialize PyLinInterp Object
  static int PyLinInterp_init ( PyLinInterp * self, PyObject * args, PyObject * kwds ) {
    
    return _PyInterp_init< PyLinInterp, utl::lin_interp >( self, args, kwds );

  }
  
  // deallocate PyLinInterp Object
  static void PyLinInterp_dealloc ( PyLinInterp * self ) {

    _PyInterp_dealloc< PyLinInterp >( self );
    return;

  }

  static PyObject * PyLinInterp_call ( PyLinInterp * self, PyObject * args ) {

    return _PyInterp_call< PyLinInterp >( self, args );

  }
  
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
    
  static PyMethodDef PyLinInterp_Methods[] =
    {
     { "integrate",
       (PyCFunction)PyLinInterp_integrate,
       METH_VARARGS,
       DocString_integrate },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject PyLinInterp_t = { PyVarObject_HEAD_INIT( NULL, 0 )
					"interpy.lin_interp"   /* tp_name */
  };

  // ========================================================================================
  // ========================================= LOG ==========================================
  // ========================================================================================

  // [ ... ]
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef interpy_module = {
					      PyModuleDef_HEAD_INIT,
					      "interpy",
					      "",
					      -1,
					      NULL, NULL, NULL, NULL, NULL
  }; /* endPyModuleDef interpy_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_interpy( void ) {

    /* --------------------- */
    /* Create interpy module */
    /* --------------------- */
    PyObject * m;
    m = PyModule_Create( &interpy_module );
    if ( m == NULL )
        return NULL;

    /* ---------------------------------------------- */
    /* Adding new object lin_interp to interpy module */
    /* ---------------------------------------------- */
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

    /* ---------------------------------------------- */
    /* Adding new object log_interp to interpy module */
    /* ---------------------------------------------- */
    // [ ... ]

    /* return new module */
    return m;
    
  }
  
} // endextern "C"

// ========================================================================================
