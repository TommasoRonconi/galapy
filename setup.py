import setuptools
import os
import sysconfig
import numpy as np

from setuptools import setup, find_packages

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
# Functions for building the library

def return_extensions () :
    from pybind11.setup_helpers import Pybind11Extension
    
    #############################################################################
    # C++ implementation of the interpolation class

    ext_intp = Pybind11Extension(
        "galapy.internal.interp",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_interpolation.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )
    
    #############################################################################
    # C++ implementation of SFH functions and types

    ext_sfh = Pybind11Extension(
        "galapy.SFH_core",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_CSFH.cpp' ),
              os.path.join( 'c++', 'sfh', 'src', 'sfh_base.cpp' ),
              os.path.join( 'c++', 'sfh', 'src', 'sfh_insitu.cpp') ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'sfh', 'include'),
                                 os.path.join( 'c++', 'imf', 'include'),
                                 os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )

    #############################################################################
    # C++ implementation of CSP functions and types

    ext_csp = Pybind11Extension(
        "galapy.CSP_core",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_CCSP.cpp' ),
              os.path.join( 'c++', 'csp', 'src', 'csp.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'csp', 'include' ),
                                 os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )

    #############################################################################
    # C++ implementation of ISM functions and types

    ext_ism = Pybind11Extension(
        "galapy.ISM_core",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_CISM.cpp' ),
              os.path.join( 'c++', 'ism', 'src', 'ism.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'ism', 'include' ),
                                 os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )

    #############################################################################
    # C++ implementation of the Nebular Free-Free emission (Bremsstrahlung)

    ext_nff = Pybind11Extension(
        "galapy.NFF_core",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_CNFF.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'rad', 'include' ),
                                 os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )

    #############################################################################
    # C++ implementation of the Synchrotron emission (generic parametric)

    ext_syn = Pybind11Extension(
        "galapy.SYN_core",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_CSYN.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'rad', 'include' ),
                                 os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )    

    #############################################################################
    # C++ implementation of BPT functions and types

    ext_bpt = Pybind11Extension(
        "galapy.BandpassTransmission",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_transmission.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'utl', 'include' ),
                                 os.path.join( 'pybind11' ) ] ),
        libraries = [ "m" ],
        extra_compile_args=[ '-std=c++17' ]
    )

    return [
        ext_intp,
        ext_sfh,
        ext_csp,
        ext_ism,
        ext_nff,
        ext_syn,
        ext_bpt
    ]

def main():
    
    #############################################################################
    # Call setup
    
    setup( name        = "galapy-fit",
           version     = get_version( os.path.join('galapy', '__init__.py') ),
           description = "GalaPy - Spectral modelling tool for galaxies in Python",
           package_dir = {
               'galapy' : 'galapy',
               'galapy.sampling' : os.path.join( 'galapy', 'sampling' ),
               'galapy.configuration' : os.path.join( 'galapy', 'configuration' ),
               'galapy.internal' : os.path.join( 'galapy', 'internal' ),
               'galapy.analysis' : os.path.join( 'galapy', 'analysis' ),
               'galapy.io' : os.path.join( 'galapy', 'io' ),
           },
           packages = [ 'galapy',
                        'galapy.configuration',
                        'galapy.internal',
                        'galapy.sampling',
                        'galapy.analysis',
                        'galapy.io' ],
           ext_modules = return_extensions(),
           include_package_data = True,
           entry_points = {
               'console_scripts' : [
                   'galapy-fit = galapy.sampling.Run:_run',
                   'galapy-genparams = galapy.sampling.Run:_generate_parameter_file',
                   'galapy-download-database = galapy.internal.data:_entrypoint_download_database',
               ]
           },
           install_requires = [
               'pybind11',
               'numpy',
               'scipy',
               'tqdm',
               'emcee',
               'dynesty',
               'matplotlib',
               'getdist',
               'requests',
               'h5py',
               'pytest'
           ]
    )

if __name__ == "__main__":
    main()
