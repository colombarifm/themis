SPHINXBUILD = sphinx-build
SOURCEDIR = source
BUILDDIR = build

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html

pdf:
	$(SPHINXBUILD) -b latex $(SOURCEDIR) $(BUILDDIR)/latex
	cd $(BUILDDIR)/latex && make

clean:
	rm -rf $(BUILDDIR)

.PHONY: html pdf clean

