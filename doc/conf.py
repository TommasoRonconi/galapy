# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import subprocess
import glob
import os
import sys
sys.path.insert( 0, os.path.abspath( '/home/tomi/.venvs/test_setups/lib/python3.8/site-packages' ) )
sys.path.insert( 0, os.path.abspath( '/home/tomi/.venvs/test_setups/lib/python3.8/site-packages/galapy' ) )
import galapy
import galapy.internal
#sys.path.insert( 0, os.path.abspath( '/opt/sadfit/sadfit-0.0.1test/lib/python3.8/site-packages/sadfit' ) )
#sys.path.insert( 0, os.path.abspath( '/opt/sadfit/sadfit-0.0.1test/lib/python3.8/site-packages/' ) )
#sys.path.insert( 0, os.path.abspath( '../python_sector' ) )
#sys.path.insert( 0, os.path.abspath( '../python_sector/sadfit' ) )
#sys.path.insert( 0, os.path.abspath( '../python_sector/sadfit/internal' ) )

# -- Configuration for ReadTheDocs setup -------------------------------------
 
# def configureDoxyfile(input_dir, output_dir):
#     with open('../../doxyfile.in', 'r') as file :
#         filedata = file.read()
 
#     filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
#     filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)
 
#     with open('doxyfile', 'w') as file:
#         file.write(filedata)
 
# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
 
# breathe_projects = {}
 
# if read_the_docs_build:
#     input_dir = '../../../mbh'
#     output_dir = 'build'
#     configureDoxyfile(input_dir, output_dir)
#     subprocess.call('doxygen', shell=True)
#     breathe_projects['MBH'] = output_dir + '/xml'

# -- Convert the tutorials ----------------------------------------------------

# for fn in glob.glob("../examples/*.ipynb"):
#     name = os.path.splitext(os.path.split(fn)[1])[0]
#     outfn = os.path.join("tutorials", name + ".rst")
#     print("Building {0}...".format(name))
#     subprocess.check_call(
#         "jupyter nbconvert --to rst "
#         + fn
#         + " --output-dir tutorials",
#         shell=True,
#     ) #  --template tutorials/tutorial_rst
#     subprocess.check_call("python fix_internal_links.py " + outfn, shell=True)
    

# -- Project information -----------------------------------------------------

project = 'GaLapy'
copyright = '2021, Tommaso Ronconi'
author = 'Tommaso Ronconi'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    #'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon'
]

# Breathe configuration
# breathe_default_project = 'ScamPy'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = ['.rst', '.md', '.txt', '.ipynb']

master_doc = 'index'

# -- Autodoc configuration ---------------------------------------------------

autodoc_mock_imports = [ 'numpy', 'scipy', 'galapy.internal' ]#, 'galapy', 'galapy.internal' ]
autodoc_default_options = {
    'members': True
    }

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# import sphinx_pdj_theme
# html_theme = 'sphinx_pdj_theme'
# html_theme_path = [sphinx_pdj_theme.get_html_theme_path()]

# html_theme = 'sphinx_rtd_theme'
# html_theme = 'karma_sphinx_theme'
# html_theme = 'insegel'
html_theme = 'sphinx_material'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# Define position of the logo
# html_logo = '../images/SAD.png'
