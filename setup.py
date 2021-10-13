import setuptools
import os
import sysconfig
import numpy as np

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
    
extra_compile_args = []
extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ["-std=c++14", "-fPIC", "-shared"]
extra_compile_args.remove( "-Wstrict-prototypes" )

extra_link_args = []
extra_link_args += [ el
                     #if el != '-Wl,--as-needed'
                     #else '-Wl,--no-as-needed'
                     for el
                     in sysconfig.get_config_var('LDFLAGS').split() ]
extra_link_args += [ '-Wl,--no-undefined' ]

class custom_build_ext(build_ext):
    def build_extensions(self):
        # Override the compiler executables. Importantly, this
        # removes the "default" compiler flags that would
        # otherwise get passed on to to the compiler, i.e.,
        # distutils.sysconfig.get_var("CFLAGS").
        self.compiler.set_executable("compiler", " ".join([ "g++" ] + extra_compile_args ) )
        #self.compiler.set_executable("compiler_so", " ".join([ "g++" ] + extra_compile_args ) )
        #self.compiler.set_executable("compiler_cxx", " ".join([ "g++" ] + extra_compile_args ) )
        self.compiler.set_executable("linker_so", "g++")
        build_ext.build_extensions(self)

def main():

    os.environ["CC"] = "g++"
    os.environ["CXX"] = "g++"
    os.environ["LDSHARED"] = ' '.join([ 'g++'
                                        if el == 'gcc'
                                        else el
                                        for el
                                        in sysconfig.get_config_var('BLDSHARED').split() ])
    

    #############################################################################
    # C++ implementation of SFH functions and types
    
    ext_sfh = Extension( "galapy.internal.CPySFH",
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
                         language="c++",
                         libraries = [ "m" ]
    )

    #############################################################################
    # C++ implementation of CSP functions and types
    
    ext_csp = Extension( "galapy.internal.CPyCSP",
                         [ os.path.join( 'csp', 'src', 'cpy_csp.cpp' ),
                           os.path.join( 'csp', 'src', 'csp.cpp' ),
                         ],
                         include_dirs = [ os.path.join( 'csp', 'include' ),
                                          os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++",
                         libraries = [ "m" ]
    )

    #############################################################################
    # C++ implementation of ISM functions and types
    
    ext_ism = Extension( "galapy.internal.CPyISM",
                         [ os.path.join( 'ism', 'src', 'cpy_ism.cpp' ),
                           os.path.join( 'ism', 'src', 'ism.cpp' ),
                         ],
                         include_dirs = [ os.path.join( 'ism', 'include' ),
                                          os.path.join( 'utl', 'include' ),
                                          np.get_include()
                         ],
                         extra_compile_args=extra_compile_args,
                         language="c++",
                         libraries = [ "m" ]
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
           version     = "0.0.1",
           description = "GalaPy - Galactic spectral analysis tools in Python",
           package_dir = {
               'galapy' : 'galapy',
               'galapy.internal' : os.path.join( 'galapy', 'internal' )
           },
           packages = [ 'galapy', 'galapy.internal' ],
           ext_modules = [
               ext_sfh,
               ext_csp,
               ext_ism
           ],
           include_package_data = True,
    )

if __name__ == "__main__":
    main()
