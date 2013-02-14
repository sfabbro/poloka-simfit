CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-*

lib_LTLIBRARIES = libpoloka-lc.la

src_includedir = $(includedir)/poloka
src_include_HEADERS = \
	fiducial.h \
	gausspsf.h \
	lcio.h \
	lightcurve.h \
	lightcurvepoint.h \
	photstar.h \
	refstar.h \
	simfit.h \
	simfitphot.h \
	simfitvignet.h \
	vignet.h \
	vignetphot.h \
	vignetserver.h


libpoloka_lc_la_SOURCES = \
	$(src_include_HEADERS) \
	gausspsf.cc \
	lcio.cc \
	lightcurve.cc \
	lightcurvepoint.cc \
	photstar.cc \
	refstar.cc \
	simfit.cc \
	simfitphot.cc \
	simfitvignet.cc \
	vignet.cc \
	vignetphot.cc \
	vignetserver.cc

libpoloka_lc_la_CPPFLAGS = @POLOKA_CORE_CFLAGS@ @POLOKA_PSF_CFLAGS@ @POLOKA_SUB_CFLAGS@

libpoloka_lc_la_LIBADD = @POLOKA_CORE_LIBS@ @POLOKA_PSF_LIBS@ @POLOKA_SUB_LIBS@
