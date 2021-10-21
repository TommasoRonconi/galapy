SOURCEDIR=./

OUTPUTDIR=../

FILEDIR=../galapy

all: docs


docs: conf.py index.rst python_doc/*
	sphinx-build -b html $(SOURCEDIR) $(OUTPUTDIR)
