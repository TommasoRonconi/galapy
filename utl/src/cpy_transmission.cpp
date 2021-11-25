// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <transmission.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

extern "C" {

  // ========================================================================================
  // ======================================== CPyBPT ========================================
  // ========================================================================================

  typedef struct {
    PyObject_HEAD
    utl::transmission * ptrObj;
  } CPyBPT;
  
  // ========================================================================================
  
  // initialize CPyBPT Object
  static const char DocString_BPTinit[] =
    "BPT( self, tau_quench, model = 'insitu')\n"
    "--\n\n"
    "Class defining ...\n\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\n : \n"
    "\t\n";
  static int CPyBPT_init ( CPyBPT *self, PyObject *args, PyObject *kwds ) {
    
    PyArrayObject * NPyLL = NULL;
    PyArrayObject * NPyFL = NULL;
    std::vector< double > ll, fl;
    
    /* Get the passed Python object */
    if ( !PyArg_ParseTuple( args, "O!O!", &PyArray_Type, &NPyLL, &PyArray_Type, &NPyFL ) )
      return -1;
    
    /* Convert NumPy-arrays to C++ vectors */
    if ( NPyArrayToCxxVector1D< double >( NPyLL, ll ) == -1 ) return -1;
    if ( NPyArrayToCxxVector1D< double >( NPyFL, fl ) == -1 ) return -1;

    self->ptrObj = new utl::transmission{ ll, fl };

    return 0;
    
  }

  // ========================================================================================

  static void CPyBPT_dealloc ( CPyBPT * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }

  // ========================================================================================

  static PyObject * CPyBPT_call ( CPyBPT * self, PyObject * args ) {

    PyArrayObject * NPyBuf = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &NPyBuf ) ) return NULL;
    
    npy_intp size = PyArray_DIM( NPyBuf, 0 );
    double * bpt = new double [size];
    for ( std::size_t ii = 0; ii < size; ++ii )
      bpt[ ii ] = ( *( self->ptrObj ) )
	( *reinterpret_cast< double * >( PyArray_GETPTR1( NPyBuf, ii ) ) );

    return PyArray_SimpleNewFromData( 1, &size, NPY_DOUBLE, reinterpret_cast< void * >( bpt ) );

  }

  // ========================================================================================
  
  // Get the bandpass flux
  static const char DocString_get_bandpass_flux[] =
    "get_bandpass_flux( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_bandpass_flux ( CPyBPT * self, PyObject * args ) {

    PyArrayObject * NPyLL = NULL;
    PyArrayObject * NPyFL = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!O!",
			    &PyArray_Type, &NPyLL,
			    &PyArray_Type, &NPyFL ) ) return NULL;

    npy_intp Lsize = PyArray_DIM( NPyLL, 0 );
    npy_intp Fsize = PyArray_DIM( NPyFL, 0 );

    /* Here better raise an exception */
    if ( Lsize != Fsize ) return NULL;
    
    return PyFloat_FromDouble( self->ptrObj->get_bandpass_flux
			      ( reinterpret_cast< double * >( PyArray_DATA( NPyLL ) ),
				reinterpret_cast< double * >( PyArray_DATA( NPyFL ) ),
				Lsize )
			      );

  }

  // ========================================================================================
  
  static PyMethodDef CPyBPT_Methods[] =
    {
     { "get_bandpass_flux",
       (PyCFunction) CPyBPT_get_bandpass_flux,
       METH_VARARGS,
       DocString_get_bandpass_flux },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject CPyBPT_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "BandpassTransmission.BPT"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef bpt_module = {
					  PyModuleDef_HEAD_INIT,
					  "BandpassTransmission",
					  "C-Python implementation of "
					  "bandpass transmission (BPT).\n",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef bpt_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_BandpassTransmission( void ) {

    /* ---------------------------------- */
    /* Create BandpassTransmission module */
    /* ---------------------------------- */
    PyObject * m;
    m = PyModule_Create( &bpt_module );
    if ( m == NULL )
        return NULL;

    /* ---------------------------------------------------- */
    /* Adding new object BPT to BandpassTransmission module */
    /* ---------------------------------------------------- */
    // Describe Transmission object
    CPyBPT_t.tp_new       = PyType_GenericNew;
    CPyBPT_t.tp_basicsize = sizeof( CPyBPT );
    CPyBPT_t.tp_dealloc   = (destructor) CPyBPT_dealloc;
    CPyBPT_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPyBPT_t.tp_doc       = DocString_BPTinit;
    CPyBPT_t.tp_call      = (ternaryfunc) CPyBPT_call;
    CPyBPT_t.tp_methods   = CPyBPT_Methods;
    //~ CPyBPT_t.tp_members=Noddy_members;
    CPyBPT_t.tp_init      = (initproc) CPyBPT_init;
    if ( PyType_Ready( &CPyBPT_t ) < 0 )
      return NULL;
    Py_INCREF( &CPyBPT_t );
    // Add BPT object to the module
    PyModule_AddObject( m, "BPT", (PyObject *)&CPyBPT_t );

    /* -------------------------------- */
    /* Initialize NumPy Array interface */
    /* -------------------------------- */
    import_array();

    /* ----------------- */
    /* return new module */
    /* ----------------- */
    return m;
    
  }
  
  // ========================================================================================

} // endextern "C"

