## Initialization
AC_PREREQ(2.59)
AC_INIT([poloka-lc], 0.1.0, [poloka team])
AC_CONFIG_SRCDIR([poloka/simfit.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE([foreign])

## Keep this one before libtool check
AM_DISABLE_STATIC

## Checks for programs to build poloka-core library
AC_PROG_CXX
AC_PROG_LIBTOOL

## Check for mandatory poloka-core
PKG_CHECK_MODULES([POLOKA_CORE],
		  [poloka-core],,
		  [AC_MSG_ERROR([Cound not find the poloka-core library])])

## Check for mandatory poloka-sub
PKG_CHECK_MODULES([POLOKA_SUB],
		  [poloka-sub],,
		  [AC_MSG_ERROR([Cound not find the poloka-sub library])])

## Check for mandatory poloka-psf
PKG_CHECK_MODULES([POLOKA_PSF],
		  [poloka-psf],,
		  [AC_MSG_ERROR([Cound not find the poloka-psf library])])

AC_CONFIG_FILES([Makefile
		 poloka-lc.pc
                 doc/Doxyfile
                 doc/Makefile
                 poloka/Makefile
		 tools/Makefile
                 ])
AC_OUTPUT
