LDADD       = @top_builddir@/ncl/libncl.la
AM_CPPFLAGS = -I@top_srcdir@/ncl
noinst_PROGRAMS = ncltest
ncltest_SOURCES = ncltest.cpp ncltest.h
check-local: ncltest
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py $(top_builddir)/example/ncltest/ncltest $(top_srcdir)/test/OldValidIn $(top_srcdir)/test/OldValidOut

