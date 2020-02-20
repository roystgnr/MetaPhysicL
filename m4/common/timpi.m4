# SYNOPSIS
#
#   Test for TIMPI Library
#
#   AM_PATH_TIMPI( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-timpi=DIR option. Searches --with-timpi,
#   $TIMPI_DIR, and the usual places for TIMPI headers and libraries.
#
#   On success, sets TIMPI_CXXFLAGS and #defines HAVE_TIMPI.
#   Assumes package is optional unless overridden with $2=yes
#

AC_DEFUN([AX_PATH_TIMPI],
[

AC_ARG_VAR(TIMPI_DIR,[root directory of TIMPI installation])

AC_ARG_WITH(timpi,
  [AS_HELP_STRING([--with-timpi[=DIR]],[root directory of TIMPI installation (default = TIMPI_DIR)])],
  [with_timpi=$withval
if test "${with_timpi}" != yes; then
    TIMPI_PREFIX=$withval
fi
],[
# assume a sensible default of --with-timpi=yes
with_timpi=yes
if test "x${TIMPI_DIR}" != "x"; then
   TIMPI_PREFIX=${TIMPI_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_TIMPI=0
HAVE_TIMPI_LIB=0

if test "${with_timpi}" != no ; then

    if test -f "${TIMPI_PREFIX}/bin/timpi-config" ; then
        dnl According to the automake variable AM_CPPFLAGS and INCLUDES
        dnl provide the same functionality, but the latter is deprecated
        TIMPI_CPPFLAGS=`${TIMPI_PREFIX}/bin/timpi-config --include`
        dnl LIBS is a variable that automake inherits from autoconf
        TIMPI_LIBS=`${TIMPI_PREFIX}/bin/timpi-config --libs`
    fi

    ac_TIMPI_save_CPPFLAGS="$CPPFLAGS"
    ac_TIMPI_save_LIBS="$LIBS"

    CPPFLAGS="${TIMPI_CPPFLAGS} ${CPPFLAGS}"
    LIBS="${TIMPI_LIBS} ${LIBS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([timpi/timpi.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_timpi_version=ifelse([$1], ,0.10, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_timpi_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_timpi_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_timpi_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    dnl begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for timpi - version >= $min_timpi_version)
        version_succeeded=no

    	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
       	@%:@include <timpi/timpi.h>
            ]], [[
            #if TIMPI_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (TIMPI_MAJOR_VERSION >= $MAJOR_VER) && (TIMPI_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (TIMPI_MAJOR_VERSION >= $MAJOR_VER) && (TIMPI_MINOR_VERSION >= $MINOR_VER) && (TIMPI_MICRO_VERSION >= $MICRO_VER)
            /* Winner winner, chicken dinner */
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])

      if test "$version_succeeded" != "yes";then
         if test "$is_package_required" = yes; then
         	  AC_MSG_ERROR([

                                 Your TIMPI library version does not meet the minimum versioning
                                 requirements ($min_timpi_version).  Please use --with-timpi to specify the location
                                 of an updated installation or consider upgrading the system version.

            ])
         fi
      fi

      # Library availability

      AC_MSG_CHECKING([for -ltimpi linkage])

      AC_LINK_IFELSE(
      [AC_LANG_PROGRAM([#include <timpi/timpi.h>],
      [
        TIMPI::TIMPIInit init();
      ])],
      [AC_MSG_RESULT(yes)
      found_library=yes ],[AC_MSG_RESULT(no)])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CPPFLAGS="$ac_TIMPI_save_CPPFLAGS"
    LIBS="$ac_TIMPI_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
             succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
       	  AC_MSG_ERROR([TIMPI not found.  Try either --with-timpi or setting TIMPI_DIR.])
       else
          AC_MSG_NOTICE([optional TIMPI library not found])
       fi
    else
        HAVE_TIMPI=1
        dnl Write these variables into the Makefile
        AC_SUBST(TIMPI_PREFIX)
        AC_SUBST(TIMPI_CPPFLAGS)
        dnl Define C processor macro; this will show up in metaphysicl_config
        AC_DEFINE(HAVE_TIMPI,1,[Define if TIMPI is available])
        if test "$found_library" = yes; then
          dnl This is for building the test suite only
          HAVE_TIMPI_LIB=1
          AC_DEFINE(HAVE_TIMPI_LIB,1,[Define if a TIMPI library is available])
          AC_SUBST(TIMPI_LIBS)
        fi
    fi
fi

AM_CONDITIONAL(TIMPI_ENABLED,test x$HAVE_TIMPI = x1)
AM_CONDITIONAL(TIMPI_LIB_ENABLED,test x$HAVE_TIMPI_LIB = x1)

])
