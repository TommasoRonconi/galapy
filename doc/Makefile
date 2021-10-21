SOURCEDIR=./

OUTPUTDIR=/tmp/galapy_doc/

FILEDIR=/home/tomi/.venvs/test_setups/lib/python3.8/site-packages/galapy

all: docs


docs: conf.py index.rst python_doc/*
	sphinx-build -b html $(SOURCEDIR) $(OUTPUTDIR)
