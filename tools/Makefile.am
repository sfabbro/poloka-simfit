CLEANFILES = *~
MAINTAINERCLEANFILES = \
        Makefile.in \
        stamp-*

AM_CPPFLAGS = @POLOKA_CORE_CFLAGS@

AM_DEFAULT_SOURCE_EXT = .cc

bin_PROGRAMS = \
	calibrate calibratefixedpos fitnight \
	make_lightcurve modelvignet

LDADD = $(top_builddir)/poloka/libpoloka-lc.la
