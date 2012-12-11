LDADD       = @top_builddir@/ncl/libncl.la
AM_CPPFLAGS = -I@top_srcdir@/ncl
bin_PROGRAMS = NEXUSnormalizer NEXUSvalidator NCLconverter
noinst_PROGRAMS = NEXUSinspector NEX_us2ml NEXUSunion

EXTERNAL_FLAG=-e

check-local: NEXUSnormalizer
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py $(top_builddir)/example/normalizer/NEXUSnormalizer $(top_srcdir)/data/characters.nex $(top_srcdir)/test/data
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py $(top_builddir)/example/normalizer/NEXUSnormalizer $(top_srcdir)/data/sample.tre $(top_srcdir)/test/data
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py $(top_builddir)/example/normalizer//NEXUSnormalizer $(top_srcdir)/test/NTSValidIn $(top_srcdir)/test/NTSValidOut
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py -x $(top_builddir)/example/normalizer/NEXUSnormalizer $(top_srcdir)/test/funkyValidIn $(top_srcdir)/test/funkyValidOut
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py $(EXTERNAL_FLAG) $(top_builddir)/example/normalizer/NEXUSnormalizer $(top_srcdir)/test/ExternalValidIn $(top_srcdir)/test/ExternalValidOut
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py -i -e $(top_builddir)/example/normalizer/NEXUSnormalizer $(top_srcdir)/test/ExternalInvalidIn
	$(PYTHON) $(top_srcdir)/test/roundTripNCLTest.py -a $(top_builddir)/example/normalizer/NEXUSnormalizer
	

NEXUSnormalizer_SOURCES = normalizer.cpp normalizer.h

NEXUSvalidator_SOURCES = normalizer.cpp normalizer.h
NEXUSvalidator_CPPFLAGS = $(AM_CPPFLAGS) -DJUST_VALIDATE_NEXUS

NEXUSinspector_SOURCES = normalizer.cpp normalizer.h
NEXUSinspector_CPPFLAGS = $(AM_CPPFLAGS) -DJUST_REPORT_NEXUS

NEX_us2ml_SOURCES = normalizer.cpp normalizer.h us2ml.cpp
NEX_us2ml_CPPFLAGS = $(AM_CPPFLAGS) -DTO_NEXML_CONVERTER

NEXUSunion_SOURCES = normalizer.cpp normalizer.h
NEXUSunion_CPPFLAGS = $(AM_CPPFLAGS) -DMULTIFILE_NEXUS_READER

NCLconverter_SOURCES = normalizer.cpp normalizer.h us2ml.cpp converter.cpp
NCLconverter_CPPFLAGS = $(AM_CPPFLAGS) -DNCL_CONVERTER_APP
