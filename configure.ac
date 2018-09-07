 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(metaphysicl, 0.3.0.1, roystgnr@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(metaphysicl_config.h.tmp)
dnl Prefixes everything in the temporary file with METAPHYSICL_
AX_PREFIX_CONFIG_H(metaphysicl_config.h,METAPHYSICL,metaphysicl_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.11 no-define color-tests parallel-tests subdir-objects])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

dnl Turn off AM_MAINTAINER_MODE for bootstrapped release by
dnl uncommenting this:
dnl AM_MAINTAINER_MODE # git trashes timestamps, forcing autoreconf needlessly

LT_INIT

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=metaphysicl

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])



dnl------------------------------
dnl Checks for compilers
dnl------------------------------

AC_PROG_CC
AC_PROG_CXX

dnl Added for AX_CODE_COVERAGE macro to work correctly, even though there's no Fortran here.
AC_PROG_FC

dnl We're going to need at least C++11 for some of our test cases
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--disable-cxx11],
                             [build without C++11 support]),
              [AS_CASE("${enableval}",
                       [yes], [enablecxx11=yes],
                       [no],  [enablecxx11=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --disable-cxx11)])],
              [enablecxx11=auto])

dnl We're going to need at least C++14 for some of our test cases
AC_ARG_ENABLE(cxx14,
              AC_HELP_STRING([--disable-cxx14],
                             [build without C++14 support]),
              [AS_CASE("${enableval}",
                       [yes], [enablecxx14=yes],
                       [no],  [enablecxx14=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --disable-cxx14)])],
              [enablecxx14=auto])

AS_IF([test x$enablecxx14 = xyes -a x$enablecxx11 = xno],
      [AC_MSG_ERROR(can't enable C++14 while disabling C++11)])

AS_IF([test x$enablecxx14 != xno -a x$enablecxx11 != xno],
      [AX_CXX_COMPILE_STDCXX(14,noext,optional)],
      [HAVE_CXX14=0])

AM_CONDITIONAL(CXX14_ENABLED,test x$HAVE_CXX14 = x1)

AS_IF([test x$HAVE_CXX14 = x1],
      [HAVE_CXX11=1],
      [AS_IF([test x$enablecxx11 != xno],
             [AX_CXX_COMPILE_STDCXX(11,noext,optional)],
             [HAVE_CXX11=0])])
AM_CONDITIONAL(CXX11_ENABLED,test x$HAVE_CXX11 = x1)

dnl -Wall warnings, -Wall the time.
AX_CXXFLAGS_WARN_ALL

dnl---------------------------------------------------------
dnl Checks for  library prerequisites for other libraries...
dnl---------------------------------------------------------

dnl Optional check for MASA
AX_PATH_MASA(0.20,no)

dnl Optional check for VexCL
AS_IF([test x$HAVE_CXX11 = x1],
      [AX_PATH_VEXCL(0.0.0, no)],
      [AM_CONDITIONAL(VEXCL_ENABLED,false)])
if (test x$HAVE_VEXCL = x1); then
  antioch_optional_test_INCLUDES="$VEXCL_CPPFLAGS $antioch_optional_test_INCLUDES"
  antioch_optional_test_LDFLAGS="$VEXCL_LDFLAGS $antioch_optional_test_LDFLAGS"
  antioch_optional_test_LIBS="$VEXCL_LIBS $antioch_optional_test_LIBS"
fi

dnl--------------------------
dnl Checks for code coverage
dnl--------------------------
AX_CODE_COVERAGE

dnl---------------------------------
dnl Query configuration environment
dnl Must be done *before* Doxygen to
dnl get some happy variables defined.
dnl---------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------------------
dnl            Doxygen support
dnl------------------------------------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(metaphysicl, doxygen/metaphysicl.dox, docs/doxygen)

AC_CONFIG_FILES(doxygen/txt_common/about_vpath.page)

dnl------------------------------------------
dnl License stamping for non-dist builds
dnl------------------------------------------

# PB: Although the documentation says that ac_top_srcdir should be defined,
#     it doesn't get activated until _AC_SRCDIRS gets called
#     (used to be _AC_SRCPATHS), according to this thread:
#     http://lists.gnu.org/archive/html/autoconf-patches/2003-02/msg00003.html
#     My own hacking revealed that we could also do the following:
#     TOP_SEARCH_DIR=$ac_pwd/$srcdir
#     I'm not sure which will be more "future proof", but this is the alternative
#     to using $(shell) calls in the Makefile.am which 1. requires
#     GNU Make and 2. Causes automake to issue a warning
_AC_SRCDIRS(.)
ABS_TOP_SRC_DIR=$ac_abs_top_srcdir
AM_CONDITIONAL(METAPHYSICL_LICENSESTAMPEXISTS, [test -f $ABS_TOP_SRC_DIR/src/common/lic_utils/update_license.pl])

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/metaphysicl/metaphysicl_version.h)
            
dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([
  doxygen/Makefile
  doxygen/metaphysicl.dox
  Makefile
  metaphysicl.pc
  src/Makefile
  test/Makefile
  ])

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_ENV

AX_SUMMARIZE_CONFIG
