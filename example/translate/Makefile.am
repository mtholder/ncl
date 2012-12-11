LDADD       = @top_builddir@/ncl/libncl.la
AM_CPPFLAGS = -I@top_srcdir@/ncl
noinst_PROGRAMS = NEXUStranslate
NEXUStranslate_SOURCES = translateNuc.cpp
NEXUStranslate_CPPFLAGS = $(AM_CPPFLAGS)

check-local: NEXUStranslate
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py -o $(top_builddir)/example/translate/NEXUStranslate  $(top_srcdir)/test/TranslateIn $(top_srcdir)/test/TranslateOut
