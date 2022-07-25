// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <nff.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

extern "C" {
  
  // ========================================================================================
  // ======================================== CPyNFF ========================================
  // ========================================================================================

  // Base class NFF
  typedef struct {
    PyObject_HEAD
    sed::nff * ptrObj;
  } CPyNFF;

  // ========================================================================================

  // Initialize object Bremsstrahlung
  static int CPyNFF_init ( CPyNFF * self, PyObject *args ) {

    PyArrayObject * Lbuf = NULL;
    std::vector< double > Lvec;
        
    /* Get the passed Python object */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &Lbuf ) ) return -1;
  
    /* Convert NumPy-array to C++ vector */
    if ( NPyArrayToCxxVector1D< double >( Lbuf, Lvec ) == -1 ) return -1;
    
    self->ptrObj = new sed::nff{ Lvec };
    return 0;
    
  }
  
  // ========================================================================================

  static void CPyNFF_dealloc ( CPyNFF * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }
  
  // ========================================================================================
  
  static const char DocString_set_params[] =
    "Function for setting the parameters of the model.\n"
    "\nParameters"
    "\n----------"
    "\nparams : array-like\n"
    "\nReturns"
    "\n-------"
    "\n: None"; 
  static PyObject * CPyNFF_set_params ( CPyNFF * self, PyObject * args ) {

    PyArrayObject * NPyBuf = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &NPyBuf ) ) return NULL;

    /* Convert NumPy-array to C-array */
    double * params;
    std::size_t len;
    if ( NPyArrayToCArray1D< double >( NPyBuf, &params, &len ) == -1 ) return NULL;

    /* Call C++ member function */
    self->ptrObj->set_params( params );

    // equivalent to return None
    Py_RETURN_NONE;

  }
  
  // ========================================================================================
  
  static const char DocString_Te[] =
    "Computes the electron temperature of Nebulae for given gas metallicity.\n"
    "\nParameters"
    "\n----------"
    "\nZgas : scalar float\n"
    "\tAbsolute value of the gas metallicity\n"
    "\nReturns"
    "\n-------"
    "\nT : scalar float\n"
    "\tElectron temperature of nebulae.\n"; 
  static PyObject * CPyNFF_Te ( CPyNFF * self, PyObject * args ) {

    double Zgas;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "d", &Zgas ) ) return NULL;

    /* Call C++ member function */
    return PyFloat_FromDouble( std::exp( self->ptrObj->lTe( std::log10( 50 * Zgas ) ) ) );

  }
  
  // ========================================================================================
  
  static const char DocString_gff[] =
    "Computes the gaunt factor of Nebulae for given "
    "wavelenght element in grid and at given electron temperature.\n"
    "\nParameters"
    "\n----------"
    "\nil : scalar int\n"
    "\tIndex in the wavelenght grid\n"
    "\nTe : scalar float\n"
    "\tElectron temperature in units of Kelvin\n"
    "\nReturns"
    "\n-------"
    "\nT : scalar float\n"
    "\tGaunt factor of nebulae.\n"; 
  static PyObject * CPyNFF_gff ( CPyNFF * self, PyObject * args ) {

    int il;
    double Te;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "id", &il, &Te ) ) return NULL;

    /* Call C++ member function */
    return PyFloat_FromDouble( self->ptrObj->gff( il, std::log( Te ) - 4 * sed::cnst::ln_10 ) );

  }
  
  // ========================================================================================

    static const char DocString_emission[] =
    "emission( self, il, Q_H )\n"
    "--\n\n"
    "Computes the NFF emission at given index in the wavelenght-grid.\n"
    "\nParameters"
    "\n----------"
    "\nil : array of int\n"
    "\tarray of indexes of the positions in the wavelenght-grid"
    " for which to compute the emission.\n"
    "\nQ_H : float\n"
    "\tIntrinsic photo-ionization rate in units of photons per second [s^-1]\n"
    "\nReturns"
    "\n-------"
    "\nL_FF : array or scalar float\n"
    "\tthe unattenuated Bremsstrahlung emission from Nebular regions\n";
  static PyObject * CPyNFF_emission ( CPyNFF * self, PyObject * args ) {
    
    PyArrayObject * il_buf;
    double Q_H; 

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!d", &PyArray_Type, &il_buf, &Q_H ) ) return NULL;

    /* Convert Numpy-array to C-array */
    std::size_t * il_arr;
    std::size_t il_size;
    if ( NPyArrayToCArray1D< std::size_t >( (PyArrayObject*)il_buf, &il_arr, &il_size ) == -1 )
      return NULL;

    double * outarr = new double [ il_size ];
    for ( unsigned int ii = 0; ii < il_size; ++ii ) 
      outarr[ ii ] = self->ptrObj->emission( il_arr[ ii ], Q_H );

    /* Clear heap */
    delete [] il_arr;

    PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&il_size, NPY_DOUBLE,
						reinterpret_cast< void * >( outarr ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    return ret;
    
  }
  
  // ========================================================================================
  
  static PyMethodDef CPyNFF_Methods[] = {
					 { "set_params",
					   (PyCFunction) CPyNFF_set_params,
					   METH_VARARGS,
					   DocString_set_params },
					 { "Te",
					   (PyCFunction) CPyNFF_Te,
					   METH_VARARGS,
					   DocString_Te },
					 { "gff",
					   (PyCFunction) CPyNFF_gff,
					   METH_VARARGS,
					   DocString_gff },
					 { "emission",
					   (PyCFunction) CPyNFF_emission,
					   METH_VARARGS,
					   DocString_emission },
					 {NULL, NULL, 0, NULL}
  };

  static PyTypeObject CPyNFF_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "NFF_core.CNFF"   /* tp_name */
  };
  
  // ========================================================================================

  static struct PyModuleDef nff_module = {
					  PyModuleDef_HEAD_INIT,
					  "NFF_core",
					  "Python wrap of c++ NFF component implementation.\n"
					  "Build an object of type nff as:\n"
					  ">>> import galapy.NFF_core as cnff\n"
					  ">>> import numpy as np\n"
					  ">>> ll = np.logspace(1,8,100)\n"
					  ">>> nff = cnff.CNFF( ll )\n",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef nff_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_NFF_core( void ) {

    /* -------------------- */
    /* Create CPyNFF module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &nff_module );
    if ( m == NULL )
        return NULL;

    /* -------------------------------------- */
    /* Adding new object CNFF to CPyNFF module */
    /* -------------------------------------- */
    // Describe CNFF object
    CPyNFF_t.tp_new       = PyType_GenericNew;
    CPyNFF_t.tp_basicsize = sizeof( CPyNFF );
    CPyNFF_t.tp_dealloc   = (destructor) CPyNFF_dealloc;
    CPyNFF_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPyNFF_t.tp_doc       = "Nebular Free-Free object";
    CPyNFF_t.tp_methods   = CPyNFF_Methods;
    CPyNFF_t.tp_init      = (initproc) CPyNFF_init;
    if ( PyType_Ready( &CPyNFF_t ) < 0 )
      return NULL;
    Py_INCREF( &CPyNFF_t );
    // Add CNFF object to the module
    PyModule_AddObject( m, "CNFF", (PyObject *)&CPyNFF_t );

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
