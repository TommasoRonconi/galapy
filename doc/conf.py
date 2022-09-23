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
# sys.path.insert( 0, os.path.abspath( '/home/tomi/.venvs/test_setups/lib/python3.8/site-packages' ) )
# sys.path.insert( 0, os.path.abspath( '/home/tomi/.venvs/test_setups/lib/python3.8/site-packages/galapy' ) )
# sys.path.insert( 0, os.path.abspath( '/home/tomi/.venvs/test_setups/lib/python3.8/site-packages/galapy/internal' ) )
import galapy
import galapy.sampling
import galapy.configuration
import galapy.internal

# -- Configuration for ReadTheDocs setup -------------------------------------

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
 
# -- Convert the tutorials ----------------------------------------------------

for fn in glob.glob("../notebooks/*.ipynb"):
    name = os.path.splitext(os.path.split(fn)[1])[0]
    outfn = os.path.join("tutorials", name + ".rst")
    print("Building {0}...".format(name))
    subprocess.check_call(
        "jupyter nbconvert --to rst "
        + fn
        + " --output-dir tutorials",
        shell=True,
    ) #  --template tutorials/tutorial_rst
    subprocess.check_call("python fix_internal_links.py " + outfn, shell=True)    

# -- Project information -----------------------------------------------------

project = 'GaLapy'
copyright = '2022, Tommaso Ronconi'
author = 'Tommaso Ronconi'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    #'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [ 'requirements.txt' ]

source_suffix = ['.rst', '.md', '.txt', '.ipynb']

master_doc = 'index'

# -- Autodoc configuration ---------------------------------------------------

autodoc_mock_imports = [ 'numpy', 'scipy' ]
autodoc_default_options = {
    'members': True
    }

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# import sphinx_pdj_theme
# html_theme = 'sphinx_pdj_theme'
# html_theme_path = [sphinx_pdj_theme.get_html_theme_path()]

html_theme = 'sphinx_rtd_theme'
# html_theme = 'karma_sphinx_theme'
# html_theme = 'insegel'
# html_theme = 'sphinx_material'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [ "_static" ]

def setup(app):
    # overrides for wide tables in RTD theme
    app.add_css_file("theme_overrides.css")  # path relative to static

# Define position of the logo
# html_logo = '../images/SAD.png'

# -- Options for LaTeX output ---------------------------------------------

latex_engine = 'pdflatex'
latex_elements = {}

# latex_logo = '_static/logo.jpg'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'main.tex', 'GaLapy Documentation',
     'Tommaso Ronconi', 'manual')
]
