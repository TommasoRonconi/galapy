import setuptools
import os
import sysconfig
import numpy as np

from setuptools import setup, find_packages, Extension

################################################################################
# Functions for reading the version 

def read( rel_path ):
    here = os.path.abspath( os.path.dirname( __file__ ) )
    with open( os.path.join( here, rel_path ), 'r' ) as fp :
        return fp.read()

def get_version( rel_path ):
    for line in read( rel_path ).splitlines() :
        if line.startswith( '__version__' ) :
            delim = '"' if '"' in line else "'"
            return line.split( delim )[ 1 ]
    else:
        raise RuntimeError( "Unable to find version string." )


################################################################################
# Global variables for compiler

# extra_compile_args = []
# # extra_compile_args = [ el
# #                        for el
# #                        in sysconfig.get_config_var('CFLAGS').split()
# #                        if ( el != '-Wstrict-prototypes' ) and ( el != '-O2' ) ]
# extra_compile_args += ["-DNDEBUG", "-O3", "-std=c++14", "-fPIC", "-shared"]
extra_compile_args = [ "-DNDEBUG", "-O3" ]

# extra_link_args = []
# extra_link_args += [ el
#                      #if el != '-Wl,--as-needed'
#                      #else '-Wl,--no-as-needed'
#                      for el
#                      in sysconfig.get_config_var('LDFLAGS').split() ]
# extra_link_args += [ '-Wl,--no-undefined' ]

def main():

    os.environ["CC"] = "g++"
    os.environ["CXX"] = "g++"    

    #############################################################################
    # C++ implementation of the interpolation class
    
    ext_intp = Extension( "galapy.internal.interp",
                          [ os.path.join( 'utl', 'src', 'cpy_interpolation.cpp' )
                          ],
                          include_dirs = [ os.path.join( 'utl', 'include' ),
                                           np.get_include()
                          ],
                          extra_compile_args=extra_compile_args,
                          language="c++14",
                          libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of SFH functions and types
    
    ext_sfh = Extension( "galapy.SFH_core",
                         [ os.path.join( 'sfh', 'src', 'cpy_sfh.cpp' ),
                           os.path.join( 'sfh', 'src', 'sfh_base.cpp' ),
                           os.path.join( 'sfh', 'src', 'sfh_insitu.cpp')
                         ],
                         include_dirs = [ os.path.join( 'sfh', 'include'),
                                          os.path.join( 'imf', 'include'),
                                          os.path.join( 'utl', 'include'),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++14",
                         libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of CSP functions and types
    
    ext_csp = Extension( "galapy.CSP_core",
                         [ os.path.join( 'csp', 'src', 'cpy_csp.cpp' ),
                           os.path.join( 'csp', 'src', 'csp.cpp' ),
                         ],
                         include_dirs = [ os.path.join( 'csp', 'include' ),
                                          os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++14",
                         libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of ISM functions and types
    
    ext_ism = Extension( "galapy.ISM_core",
                         [ os.path.join( 'ism', 'src', 'cpy_ism.cpp' ),
                           os.path.join( 'ism', 'src', 'ism.cpp' ),
                         ],
                         include_dirs = [ os.path.join( 'ism', 'include' ),
                                          os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++14",
                         libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of the Nebular Free-Free emission (Bremsstrahlung)
    
    ext_nff = Extension( "galapy.NFF_core",
                         [ os.path.join( 'rad', 'src', 'cpy_nff.cpp' ),
                         ],
                         include_dirs = [ os.path.join( 'rad', 'include' ),
                                          os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++14",
                         libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of the Synchrotron emission (generic parametric)
    
    ext_syn = Extension( "galapy.SYN_core",
                         [ os.path.join( 'rad', 'src', 'cpy_syn.cpp' ),
                         ],
                         include_dirs = [ os.path.join( 'rad', 'include' ),
                                          os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++14",
                         libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of BPT functions and types
    
    ext_bpt = Extension( "galapy.BandpassTransmission",
                         [ os.path.join( 'utl', 'src', 'cpy_transmission.cpp' )
                         ],
                         include_dirs = [ os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++14",
                         libraries = [ "m", "stdc++" ]
    )

    #############################################################################
    # C++ implementation of IGM functions and types
            
    # ext_igm = Extension( "galapy.internal.PyIGM",
    #                      [ "cpython/cpython_igm.cpp",
    #                        "cpp_modules/igm/src/igm.cpp" ],
    #                      include_dirs = [ "cpp_modules/igm/include",
    #                                       "cpp_modules/utilities/include",
    #                                       "cpython/include" ],
    #                      extra_compile_args=extra_compile_args,
    #                      language="c++14",
    #                      libraries = [ "gsl", "gslcblas" ],
    # )

    #############################################################################
    # Call setup
    
    setup( name        = "galapy",
           version     = get_version( os.path.join('galapy', '__init__.py') ),
           description = "GalaPy - Galactic spectral analysis tools in Python",
           package_dir = {
               'galapy' : 'galapy',
               'galapy.sampling' : os.path.join( 'galapy', 'sampling' ),
               'galapy.configuration' : os.path.join( 'galapy', 'configuration' ),
               'galapy.internal' : os.path.join( 'galapy', 'internal' ),
               'galapy.analysis' : os.path.join( 'galapy', 'analysis' )
           },
           packages = [ 'galapy',
                        'galapy.configuration',
                        'galapy.internal',
                        'galapy.sampling',
                        'galapy.analysis' ],
           ext_modules = [
               ext_intp,
               ext_sfh,
               ext_csp,
               ext_ism,
               ext_nff,
               ext_syn,
               ext_bpt
           ],
           include_package_data = True,
           entry_points = {
               'console_scripts' : [
                   'galapy-fit = galapy.sampling.Run:run',
                   'galapy-genparams = galapy.sampling.Run:generate_parameter_file',
                   'galapy-download-database = galapy.internal.data:_entrypoint_download_database',
               ]
           },
           install_requires = [
               'numpy',
               'scipy',
               'emcee',
               'dynesty',
               'matplotlib',
               'getdist',
               'requests',
               'pytest'
           ]
    )

if __name__ == "__main__":
    main()
