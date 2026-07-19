# Minimal Makefile for Sphinx docs
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: html clean tutorials-sync figures

html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)

tutorials-sync:
	jupytext --sync tutorials/*.py

figures:
	python make_docs_figures.py
	python make_showcase_figures.py

clean:
	rm -rf "$(BUILDDIR)"
	rm -rf api/generated
