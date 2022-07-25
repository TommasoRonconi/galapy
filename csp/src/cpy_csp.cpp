// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
// STL headers
#include <vector>
#include <string>
#include <csp.h>
#include <cstring>
#include <fstream>
#include <algorithm>

extern "C" {

  // ========================================================================================
  // =================================== Non-Type Methods ===================================
  // ========================================================================================
  
  static const char DocString_loadSSP[] =
    "loadSSP( file_name )\n"
    "--\n\n"
    "Function for ...\n"
    "\nParameters"
    "\n----------"
    "\nfile_name : string\n"
    "\tposition in the file-system of the binary file "
    "storing the SSP table."
    "\nReturns"
    "\n-------"
    "\n : "; 
  static PyObject * method_loadSSP ( PyObject * self, PyObject * args ) {

    PyArrayObject * Pyl = NULL, * Pyt = NULL, * PyZ = NULL, * PySSP = NULL;
    const char * file_name;

    if ( !PyArg_ParseTuple( args, "s", &file_name ) ) return NULL;

    // Open binary input stream
    std::ifstream ifs ( file_name, std::ios_base::in | std::ios_base::binary );
    if ( ifs.fail() ) {
      PyErr_SetString( PyExc_FileNotFoundError,
  		       "Unable to open input file." );
      return NULL;
    }

    // [ read here ]
    npy_intp Nl, Nt, NZ, NSSP;
    double *lambda, *tau, *Z, *SSP;

    // read the wavelenght sampling vector from file
    ifs.read( ( char * ) & Nl, sizeof( std::size_t ) );
    lambda = new double [ Nl ];
    ifs.read( reinterpret_cast< char* >( lambda ), sizeof( double ) * Nl );
    Pyl = ( PyArrayObject * )PyArray_SimpleNewFromData( 1, &Nl, NPY_DOUBLE,
							reinterpret_cast< void * >( lambda ) );
    PyArray_ENABLEFLAGS(Pyl, NPY_ARRAY_OWNDATA);

    // read the time sampling vector from file
    ifs.read( ( char * ) & Nt, sizeof( std::size_t ) );
    tau = new double [ Nt ];
    ifs.read( reinterpret_cast< char* >( tau ), sizeof( double ) * Nt );
    Pyt = ( PyArrayObject * )PyArray_SimpleNewFromData( 1, &Nt, NPY_DOUBLE,
							reinterpret_cast< void * >( tau ) );
    PyArray_ENABLEFLAGS(Pyt, NPY_ARRAY_OWNDATA);

    // read the time sampling vector from file
    ifs.read( ( char * ) & NZ, sizeof( std::size_t ) );
    Z = new double [ NZ ];
    ifs.read( reinterpret_cast< char* >( Z ), sizeof( double ) * NZ );
    PyZ = ( PyArrayObject * )PyArray_SimpleNewFromData( 1, &NZ, NPY_DOUBLE,
							reinterpret_cast< void * >( Z ) );
    PyArray_ENABLEFLAGS(PyZ, NPY_ARRAY_OWNDATA);

    // read the time sampling vector from file
    ifs.read( ( char * ) & NSSP, sizeof( std::size_t ) );
    SSP = new double [ NSSP ];
    ifs.read( reinterpret_cast< char* >( SSP ), sizeof( double ) * NSSP );
    PySSP = ( PyArrayObject * )PyArray_SimpleNewFromData( 1, &NSSP, NPY_DOUBLE,
							  reinterpret_cast< void * >( SSP ) );
    PyArray_ENABLEFLAGS(PySSP, NPY_ARRAY_OWNDATA);
    
    // Raise RuntimeError if EoF is not reached:
    if ( !ifs ) {
      ifs.clear(); ifs.close();
      PyErr_SetString( PyExc_RuntimeError,
  		       "Error reading input file: "
  		       "reading ended without reaching EoF" );
      return NULL;
    }
    ifs.clear(); ifs.close();

    // here I want to return a tuple with ( l, t, Z, LltZ )
    return PyTuple_Pack( 4, Pyl, Pyt, PyZ, PySSP );
    
  }
    
  // ========================================================================================

  static PyMethodDef CPyCSP_NT_Methods[] =
    {
     { "loadSSP",
       (PyCFunction) method_loadSSP,
       METH_VARARGS,
       DocString_loadSSP },
     { NULL, NULL, 0, NULL }
    };
  
  // ========================================================================================
  // ======================================== CPyCSP ========================================
  // ========================================================================================

  typedef struct {
    PyObject_HEAD
    sed::csp * ptrObj;
  } CPyCSP;

  // ========================================================================================

  static const char DocString_CCSPinit[] =
    "CCSP( self, lambda, tau, Z, SSPtable )\n"
    "--\n\n"
    "Composite Stellar Population (CSP) object.\n"
    "\nParameters"
    "\n----------"
    "\nlambda   : array-like\n"
    "\ntau      : array-like\n"
    "\nZ        : array-like\n"
    "\nSSPtable : array-like\n"; 
  static int CPyCSP_init ( CPyCSP *self, PyObject *args, PyObject *kwds ) {

    PyArrayObject * Lbuf = NULL, * Tbuf = NULL, * Zbuf = NULL, * SSPbuf = NULL;
    std::vector< double > Lvec, Tvec, Zvec, SSPvec;
    int do_CCSN_rate;
        
    /* Get the passed Python object */
    if ( !PyArg_ParseTuple( args, "O!O!O!O!i",
			    &PyArray_Type, &Lbuf,
			    &PyArray_Type, &Tbuf,
			    &PyArray_Type, &Zbuf,
			    &PyArray_Type, &SSPbuf,
			    &do_CCSN_rate ) ) return -1;
  
    /* Convert NumPy-array to C++ vector */
    if ( NPyArrayToCxxVector1D< double >( Lbuf, Lvec ) == -1 ) return -1;
    if ( NPyArrayToCxxVector1D< double >( Tbuf, Tvec ) == -1 ) return -1;
    if ( NPyArrayToCxxVector1D< double >( Zbuf, Zvec ) == -1 ) return -1;
    if ( NPyArrayToCxxVector1D< double >( SSPbuf, SSPvec ) == -1 ) return -1;

    /* Call member constructor */
    self->ptrObj = new sed::csp{ Lvec, Tvec, Zvec, SSPvec, do_CCSN_rate };

    /* Return 0 on success */
    return 0;
    
  }

  // ========================================================================================

  static void CPyCSP_dealloc ( CPyCSP * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }

  // ========================================================================================
  
  static const char DocString_set_params[] =
    "set_params( self, psi, Zstar, iZl, itl )\n"
    "--\n\n"
    "Function for setting the parameters of the model.\n"
    "\nParameters"
    "\n----------"
    "\npsi   : array-like\n"
    "\nZstar : array-like\n"
    "\niZl   : array-like\n"
    "\nitl   : integer\n"
    "\nReturns"
    "\n-------"
    "\n: None"; 
  static PyObject * CPyCSP_set_params ( CPyCSP * self, PyObject * args ) {

    PyArrayObject * psi_buf = NULL, * Zstar_buf = NULL, * izl_buf = NULL;
    unsigned long itl;
    std::vector< double > psi_vec, Zstar_vec;
    std::vector< std::size_t > izl_vec;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!O!O!k",
			    &PyArray_Type,  &psi_buf,
			    &PyArray_Type,  &Zstar_buf,
			    &PyArray_Type,  &izl_buf,
			    &itl ) ) return NULL;

  
    /* Convert NumPy arrays to C++ vectors */
    if ( NPyArrayToCxxVector1D< double >( psi_buf, psi_vec ) == -1 ) return NULL;
    if ( NPyArrayToCxxVector1D< double >( Zstar_buf, Zstar_vec ) == -1 ) return NULL;
    if ( NPyArrayToCxxVector1D< std::size_t >( izl_buf, izl_vec ) == -1 ) return NULL;

    /* Call C++ member function */
    self->ptrObj->set_params( psi_vec, Zstar_vec, izl_vec, (std::size_t)itl );

    // equivalent to return None
    Py_RETURN_NONE;

  }
  
  // ========================================================================================
  
  static const char DocString_getitem[] =
    "SSP(self, il, it, iz)\n"
    "--\n\n"
    "Function for ...\n"
    "\nParameters"
    "\n----------"
    "\nil : uint\n"
    "\nit : uint\n"
    "\niz : uint\n"
    "\nReturns"
    "\n-------"
    "\n : "; 
  static PyObject * CPyCSP_getitem ( CPyCSP * self, PyObject * args ) {

    int il, it, iz;
    
    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "iii", &il, &it, &iz ) ) {
      return NULL;
    }

    /* Call C++ member function */
    return PyFloat_FromDouble( self->ptrObj->luminosity( il, it, iz ) );

  }

  // ========================================================================================

  static const char DocString_emission[] =
    "emission( self, il, Ftau = 1. )\n"
    "--\n\n"
    "Computes the CSP emission at given index in the wavelenght-grid.\n"
    "It approximates the integral:\n"
    "\n.. math::\n"
    "\n\tL_\\lambda^\\text{CSP,unatt}(\\tau') = "
    "\\int_0^{\\tau'}\\text{d}\\tau F(\\tau)\\cdot "
    "L_\\lambda^\\text{SSP}\\bigl[\\tau, Z_\\ast(\\tau'-\\tau)\\bigr]\\psi(\\tau'-\\tau)\n"
    "\nParameters"
    "\n----------"
    "\nil : array of int\n"
    "\tarray of indexes of the positions in the wavelenght-grid"
    " for which to compute the emission.\n"
    "\nFtau : array\n"
    "\tarray containing a function of time to convolve the integrand with "
    "(must have same dimension of the SSP's time-grid times the lenght of the `il` array)\n"
    "\nReturns"
    "\n-------"
    "\nL_CSP : array or scalar float\n"
    "\tthe emission of the galaxy's CSP filtered by the function of time :math:`F(\\tau)`\n";
  static PyObject * CPyCSP_emission ( CPyCSP * self, PyObject * args, PyObject * kwds ) {
    
    PyArrayObject * il_buf;
    PyArrayObject * Tfact_buf = NULL;
    double * Tfact_arr; 
    
    // work-around to silence compile warning (const char* -> char* forbidden in C++):
    char ** kwlist = new char * [ 3 ] { (char*)"", (char*)"Tfact", NULL };

    /* Parse arguments */
    if ( !PyArg_ParseTupleAndKeywords( args, kwds, "O!|O!", kwlist,
				       &PyArray_Type, &il_buf,
				       &PyArray_Type, &Tfact_buf ) ) return NULL;
    
    /* Clear heap */
    delete [] kwlist;

    /* Convert Numpy-array to C-array */
    std::size_t * il_arr;
    std::size_t il_size;
    if ( NPyArrayToCArray1D< std::size_t >( (PyArrayObject*)il_buf, &il_arr, &il_size ) == -1 ) return NULL;

    /* If kwarg 'Tfact' is provided, copy its content into C-style array, 
       otherwise allocate C-style array filled with 1s of the correct size */
    std::size_t tdim;
    if ( Tfact_buf ) {
      if ( NPyArrayToCArray1D< double >( Tfact_buf, &Tfact_arr, &tdim ) == -1 ) return NULL;
    }
    else {
      tdim = il_size * self->ptrObj->tau_size();
      Tfact_arr = new double [ tdim ];
      std::fill_n( Tfact_arr, tdim, 1. );
    }

    double * outarr = new double [ il_size ];
    for ( unsigned int ii = 0; ii < il_size; ++ii ) 
      outarr[ ii ] = self->ptrObj->emission( il_arr[ ii ], Tfact_arr + self->ptrObj->tau_size() * ii );

    /* Clear heap */
    delete [] Tfact_arr;
    delete [] il_arr;

    PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&il_size, NPY_DOUBLE,
						reinterpret_cast< void * >( outarr ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    return ret;
    
  } 
  
  // ========================================================================================

  static const char DocString_RCCSN[] =
    "Returns the extinction value in the visible band.\n";
  static PyObject * CPyCSP_RCCSN ( CPyCSP * self ) {

    return PyFloat_FromDouble( self->ptrObj->RCCSN() );
    
  }  
  
  // ========================================================================================
  
  static PyMethodDef CPyCSP_Methods[] = {
					 { "set_params",
					   (PyCFunction) CPyCSP_set_params,
					   METH_VARARGS,
					   DocString_set_params },
					 { "SSP",
					   (PyCFunction) CPyCSP_getitem,
					   METH_VARARGS,
					   DocString_getitem },
					 { "emission",
					   (PyCFunction) CPyCSP_emission,
					   METH_VARARGS | METH_KEYWORDS,
					   DocString_emission },
					{ "RCCSN",
					  (PyCFunction) CPyCSP_RCCSN,
					  METH_NOARGS,
					  DocString_RCCSN },
					 {NULL, NULL, 0, NULL}
  };

  static PyTypeObject CPyCSP_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "CSP_core.CCSP"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef csp_module = {
					  PyModuleDef_HEAD_INIT,
					  "CSP_core",
					  "Python wrap of c++ CSP component implementation.\n"
					  "Build an object of type csp as:\n"
					  "\t>>> import galapy.CSP_core as ccsp\n"
					  "\t>>> csp = ccsp.CCSP()",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef csp_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_CSP_core( void ) {

    /* -------------------- */
    /* Create CPyCSP module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &csp_module );
    if ( m == NULL )
        return NULL;

    /* --------------------------------------- */
    /* Adding new object CCSP to CPyCSP module */
    /* --------------------------------------- */
    // Describe CCSP object
    CPyCSP_t.tp_new       = PyType_GenericNew;
    CPyCSP_t.tp_basicsize = sizeof( CPyCSP );
    CPyCSP_t.tp_dealloc   = (destructor) CPyCSP_dealloc;
    CPyCSP_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPyCSP_t.tp_doc       = DocString_CCSPinit; //"CSP objects";
    // CPyCSP_t.tp_call      = (ternaryfunc) CPyCSP_call;
    CPyCSP_t.tp_methods   = CPyCSP_Methods;
    //~ CPyCSP_t.tp_members=Noddy_members;
    CPyCSP_t.tp_init      = (initproc) CPyCSP_init;
    if ( PyType_Ready( &CPyCSP_t ) < 0 )
      return NULL;
    Py_INCREF( &CPyCSP_t );
    // Add CCSP object to the module
    PyModule_AddObject( m, "CCSP", (PyObject *)&CPyCSP_t );

    /* ---------------------------------------- */
    /* Adding Non-Type methods to CPyCSP module */
    /* ---------------------------------------- */
    PyModule_AddFunctions( m, CPyCSP_NT_Methods );

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
