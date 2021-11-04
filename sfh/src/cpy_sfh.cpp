// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <sfh.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

PyObject * _Create_SFHmodelDict () {

  PyObject * dict = PyDict_New();
  if ( PyDict_SetItemString( dict, "insitu",     PyLong_FromLong( 1 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "constant",   PyLong_FromLong( 2 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "delayedexp", PyLong_FromLong( 3 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "lognormal",  PyLong_FromLong( 4 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "burst",      PyLong_FromLong( 5 ) ) < 0 ) return NULL;
  return dict;
  
}

extern "C" {

  // ========================================================================================
  // ======================================== CPySFH ========================================
  // ========================================================================================

  typedef struct {
    PyObject_HEAD
    sed::sfh_base * ptrObj;
  } CPySFH;

  // ========================================================================================

  static PyObject * sfh_model = _Create_SFHmodelDict();

  // ========================================================================================

  // initialize CPySFH Object
  static const char DocString_CSFHinit[] =
    "CSFH( self, tau_quench, model = 'insitu')\n"
    "--\n\n"
    "Class defining the model of Star Formation History (SFH).\n\n"
    "The possible models to choose are\n\n"
    "#. 'insitu'\n"
    "#. 'constant'\n"
    "#. 'delayedexp'\n"
    "#. 'lognormal'\n"
    "#. 'burst'\n"
    "\nParameters"
    "\n----------"
    "\ntau_quench : float\n"
    "\tEventual abrupt quenching time for star formation.\n"
    "\tShould be expressed in years. Refers to the age of the galaxy.\n"
    "\nmodel : string\n"
    "\tOne among ('insitu', 'constant', 'delayedexp', 'lognormal', 'burst').\n"
    "\tDefault is 'insitu'.\n";
  static int CPySFH_init ( CPySFH *self, PyObject *args, PyObject *kwds ) {
    
    double tau_quench;
    char const * model = "insitu";
    
    // work-around to silence compile warning (const char* -> char* forbidden in C++):
    char ** kwlist = new char * [ 3 ] { (char*)"", (char*)"model", NULL };
    
    /* Get the passed Python object */
    if ( !PyArg_ParseTupleAndKeywords( args, kwds, "d|s", kwlist, &tau_quench, &model ) ) {
      return -1;
    }

    switch ( (int) PyLong_AsLong( PyDict_GetItemString( sfh_model, model ) ) ) {

    case 1 :
      self->ptrObj = new sed::sfh_insitu{ tau_quench };
      break;

    case 2 :
      self->ptrObj = new sed::sfh_constant{ tau_quench };
      break;

    case 3 :
      self->ptrObj = new sed::sfh_delayedexp{ tau_quench };
      break;

    case 4 :
      self->ptrObj = new sed::sfh_lognorm{ tau_quench };
      break;

    // case 5 :
    //   self->ptrObj = new sed::sfh_burst{};
    //   break;

    default :
      PyErr_SetString( PyExc_TypeError,
    		       "SFH model not valid. "
    		       "Valid models are: "
    		       "'insitu', 'constant', "
    		       "'delayedexp', 'lognormal', "
    		       "'burst'" );
      return -1;

    } // endswitch ( model )

    /* Clear heap */
    delete [] kwlist;

    return 0;
    
  }

  // ========================================================================================

  static void CPySFH_dealloc ( CPySFH * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }

  // ========================================================================================

  static PyObject * CPySFH_call ( CPySFH * self, PyObject * args ) {

    return Py_call_fromScalarOrArray< CPySFH >( self, args );

  }

  // ========================================================================================
  
  static const char DocString_set_params[] =
    "set_params( self, params )\n"
    "--\n\n"
    "Function for setting the parameters of the model.\n"
    "\nParameters"
    "\n----------"
    "\nparams : 1d-array\n"
    "\tarray containing the list of parameters of the model. "
    "Parameters have to be passed in the correct order:\n\n"
    "\t#. :code:`'insitu': [psi_max, tau_star]`\n"
    "\t#. :code:`'constant': [psi]`\n"
    "\t#. :code:`'delayedexp': [psi_norm, k_shape, tau_star]`\n"
    "\t#. :code:`'lognormal': [psi_norm, sigma_star, tau_star]`\n"
    "\t#. :code:`'burst': [...]`\n"
    "\nReturns"
    "\n-------"
    "\n: None"; 
  static PyObject * CPySFH_set_params ( CPySFH * self, PyObject * args ) {

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
  
  static const char DocString_Mstar[] =
    "Mstar( self, tau, npoints = 100 )\n"
    "--\n\n"
    "computes the stellar mass at a given age of the galaxy.\n"
    "It approximates the integral:\n"
    "\n.. math::\n"
    "\n\tM_\\ast(\\tau') = \\int_0^{\\tau'}\\text{d}\\tau "
    "\\bigl[1 - \\mathcal{R}_\\text{IMF}(\\tau)\\bigr]\\psi(\\tau)\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tgalaxy age in years.\n"
    "\nnpoints : int\n"
    "\tthinness for approximated integral computation (default is 100)\n"
    "\nReturns"
    "\n-------"
    "\n: float\n"
    "\tthe stellar mass of the galaxy at time :math:`\\tau`\n"; 
  static PyObject * CPySFH_Mstar ( CPySFH * self, PyObject * args, PyObject * kwds ) {

    double tau;
    Py_ssize_t npoints = 100;
    
    // work-around to silence compile warning (const char* -> char* forbidden in C++):
    char ** kwlist = new char * [ 3 ] { (char*)"", (char*)"npoints", NULL };

    // Parse arguments
    if ( !PyArg_ParseTupleAndKeywords( args, kwds, "d|n", kwlist, &tau, &npoints ) ) {
      return NULL;
    }

    /* Clear heap */
    delete [] kwlist;

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Mstar( tau, (std::size_t)npoints ) );

  }

  // ========================================================================================
  
  static const char DocString_Mdust[] =
    "Mdust( self, tau )\n"
    "--\n\n"
    "Returns the dust mass at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is a free parameter.\n"
    "For the :code:`insitu` model, the dust mass is given by "
    ":math:`M_\\text{gas}(\\tau)D(\\tau)` where :math:`M_\\text{gas}=\\psi(\\tau)\\tau_\\ast` "
    "and where :math:`D(\\tau)` is the gas mass ratio "
    "(for an analytic expression of this quantity see Pantoni et al. 2019 and Lapi et al. 2020).\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tDust content in solar masses (:math:`M_\\odot`) at give time :math:`\\tau`.\n"; 
  static PyObject * CPySFH_Mdust ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Mdust( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_Mgas[] =
    "Mgas( self, tau )\n"
    "--\n\n"
    "Returns the gas mass at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is given by :math:`M_\\text{gas} = M_\\text{dust}/D` "
    "where :math:`D \\sim 0.01 (Z_\\text{gas}/Z_\\odot)^{-0.85}` is the "
    "dust-to-gas mass ratio, derived from observations."
    "For the :code:`insitu` model, the gas mass is given by"
    " :math:`M_\\text{gas}=\\psi(\\tau)\\tau_\\ast`\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tGas content in solar masses (:math:`M_\\odot`) at give time :math:`\\tau`.\n";
  static PyObject * CPySFH_Mgas ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Mgas( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_Zgas[] =
    "Zgas( self, tau )\n"
    "--\n\n"
    "Returns the gas metallicity at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is a free parameter with "
    ":math:`Z_\\text{gas} = Z_\\ast`.\n"
    "For the :code:`insitu` model, it is instead given by\n"
    "\n.. math::\n"
    "\n\tZ_\\text{gas}=\\dfrac{s y_Z}{s\\gamma-1}"
    "\\biggl[1 - \\dfrac{(s\\gamma-1)x}{e^{(s\\gamma-1)x}-1}\\biggr]\n"
    "\nwhere :math:`x\\equiv\\tau/s\\tau_\\ast` and :math:`y_Z\\approx0.04` "
    "is the metal production yield (including recycling) for a Chabrier IMF.\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tGas absolute metallicity at give time :math:`\\tau`.\n";
  static PyObject * CPySFH_Zgas ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Zgas( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_Zstar[] =
    "Zstar( self, tau )\n"
    "--\n\n"
    "Returns the stellar metallicity at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is a free parameter with "
    ":math:`Z_\\ast = Z_\\text{gas}`.\n"
    "For the :code:`insitu` model, it is instead given by\n"
    "\n.. math::\n"
    "\n\tZ_\\ast=\\dfrac{y_Z}{\\gamma-1}"
    "\\biggl[1 - \\dfrac{s\\gamma}{(s\\gamma-1}"
    "\\dfrac{e^{-x}-e^{-s\\gamma x}[1 + (s\\gamma -1)x]}"
    "{s\\gamma -1 + e^{-s\\gamma x}- s\\gamma e^{-x}}\\biggr]\n"
    "\nwhere :math:`x\\equiv\\tau/s\\tau_\\ast` and :math:`y_Z\\approx0.04` "
    "is the metal production yield (including recycling) for a Chabrier IMF.\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tStellar absolute metallicity at give time :math:`\\tau`.\n";
  static PyObject * CPySFH_Zstar ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Zstar( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_time_grid[] =
    "Function computing ...\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tgalactic age\n"
    "\nnpoints : int\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n"; 
  static PyObject * CPySFH_time_grid ( CPySFH * self, PyObject * args, PyObject * kwds ) {

    double age;
    PyArrayObject * tgrid_NPyBuf = NULL, * Zgrid_NPyBuf = NULL;
    PyArrayObject * NPyPsiGrid = NULL, * NPyZGrid = NULL, * NPyZidxGrid = NULL;
    std::vector< double > tgrid_CxxVec, Zgrid_CxxVec;
    
    if ( !PyArg_ParseTuple( args, "dO!O!",
			    &age,
			    &PyArray_Type, &tgrid_NPyBuf,
			    &PyArray_Type, &Zgrid_NPyBuf ) ) return NULL;

    /* Convert NumPy-array to C-array */
    if ( NPyArrayToCxxVector1D< double >( tgrid_NPyBuf, tgrid_CxxVec ) == -1 ) return NULL;
    if ( NPyArrayToCxxVector1D< double >( Zgrid_NPyBuf, Zgrid_CxxVec ) == -1 ) return NULL;

    /* Call member functions and convert to NumPy-arrays */
    self->ptrObj->time_grid( age, tgrid_CxxVec, Zgrid_CxxVec );
    NPyPsiGrid  =
      ( PyArrayObject * )
      CxxVectorToNPyArray1D< double, NPY_DOUBLE >( self->ptrObj->get_psi_grid() );
    NPyZGrid    =
      ( PyArrayObject * )
      CxxVectorToNPyArray1D< double, NPY_DOUBLE >( self->ptrObj->get_Z_grid() );
    NPyZidxGrid =
      ( PyArrayObject * )
      CxxVectorToNPyArray1D< std::size_t, NPY_UINT64 >( self->ptrObj->get_Zidx_grid() );
    
    // return tuple
    return PyTuple_Pack( 4, NPyPsiGrid, NPyZGrid, NPyZidxGrid,
  			 PyLong_FromSize_t( self->ptrObj->get_last_grid_idx() ) );
    
  }

  // ========================================================================================
  
  static PyMethodDef CPySFH_Methods[] =
    {
     { "set_params",
       (PyCFunction) CPySFH_set_params,
       METH_VARARGS,
       DocString_set_params },
     { "Mstar",
       (PyCFunction) CPySFH_Mstar,
       METH_VARARGS | METH_KEYWORDS,
       DocString_Mstar },
     { "Mdust",
       (PyCFunction) CPySFH_Mdust,
       METH_VARARGS,
       DocString_Mdust },
     { "Mgas",
       (PyCFunction) CPySFH_Mgas,
       METH_VARARGS,
       DocString_Mgas },
     { "Zgas",
       (PyCFunction) CPySFH_Zgas,
       METH_VARARGS,
       DocString_Zgas },
     { "Zstar",
       (PyCFunction) CPySFH_Zstar,
       METH_VARARGS,
       DocString_Zstar },
     { "time_grid",
       (PyCFunction) CPySFH_time_grid,
       METH_VARARGS,
       DocString_time_grid },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject CPySFH_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "SFH_core.CSFH"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef sfh_module = {
					  PyModuleDef_HEAD_INIT,
					  "SFH_core",
					  "Python wrap of c++ SFH component implementation.\n"
					  "Build an object of type sfh as:\n"
					  "\t>>> import galapy.SFH_core as csfh\n"
					  "\t>>> sfh = SFH_core.CSFH()",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef sfh_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_SFH_core( void ) {

    /* -------------------- */
    /* Create CPySFH module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &sfh_module );
    if ( m == NULL )
        return NULL;

    /* --------------------------------------- */
    /* Adding new object CSFH to CPySFH module */
    /* --------------------------------------- */
    // Describe CSFH object
    CPySFH_t.tp_new       = PyType_GenericNew;
    CPySFH_t.tp_basicsize = sizeof( CPySFH );
    CPySFH_t.tp_dealloc   = (destructor) CPySFH_dealloc;
    CPySFH_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPySFH_t.tp_doc       = DocString_CSFHinit;
    CPySFH_t.tp_call      = (ternaryfunc) CPySFH_call;
    CPySFH_t.tp_methods   = CPySFH_Methods;
    //~ CPySFH_t.tp_members=Noddy_members;
    CPySFH_t.tp_init      = (initproc) CPySFH_init;
    if ( PyType_Ready( &CPySFH_t ) < 0 )
      return NULL;
    Py_INCREF( &CPySFH_t );
    // Add CSFH object to the module
    PyModule_AddObject( m, "CSFH", (PyObject *)&CPySFH_t );

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

