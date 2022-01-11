SOURCEDIR=./

LOC=/tmp

OUTPUTDIR=$(LOC)/galapy_doc/

FILEDIR=/home/tomi/.venvs/test_setups/lib/python3.8/site-packages/galapy

all: docs_html docs_latex


docs_html: conf.py index.rst python_doc/*
	sphinx-build -b html $(SOURCEDIR) $(OUTPUTDIR)/html/

docs_latex: conf.py index.rst python_doc/*
	sphinx-build -b latex $(SOURCEDIR) $(OUTPUTDIR)/latex/

clean_html:
	rm -rf $(OUTPUTDIR)/html

clean_latex:
	rm -rf $(OUTPUTDIR)/latex

clean:
	rm -rf $(OUTPUTDIR)

purge_tutorials: 
	rm -rf $(SOURCEDIR)tutorials/*

purge_python_doc:
	rm -rf $(SOURCEDIR)python_doc/_generated/*

purge_all: clean purge_tutorials purge_python_doc
