# SYNOPSIS
#
#   Test for HelloWorld
#
#   AX_PATH_HELLOWORLD( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-helloworld=DIR option. Searches
#   --with-helloworld, $HELLOWORLD_DIR, and the usual places for
#   HELLOWORLD headers and libraries.
#
#   On success, sets HELLOWORLD_CPPFLAGS, HELLOWORLD_LIBS, and
#   #defines HAVE_HELLOWORLD.
#   Also defines automake conditional HELLOWORLD_ENABLED.  Assumes
#   package is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: helloworld.m4 37037 2013-02-16 01:03:09Z pbauman $
#
# COPYLEFT
#
#   Copyright (c) 2013 Paul T. Bauman <pbauman@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_HELLOWORLD],
[

AC_ARG_VAR(HELLOWORLD_DIR,[root directory of HelloWorld installation])

AC_ARG_WITH(helloworld,
  [AS_HELP_STRING([--with-helloworld[=DIR]],[root directory of HelloWorld installation (default = HELLOWORLD_DIR)])],
  [with_helloworld=$withval
if test "${with_helloworld}" != yes; then
    HELLOWORLD_PREFIX=$withval
fi
],[
with_helloworld=$withval
if test "x${HELLOWORLD_DIR}" != "x"; then
   HELLOWORLD_PREFIX=${HELLOWORLD_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_HELLOWORLD=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_helloworld}" != no ; then

    if test -d "${HELLOWORLD_PREFIX}/lib" ; then
       HELLOWORLD_LIBS="-L${HELLOWORLD_PREFIX}/lib -lhelloworld -Wl,-rpath,${HELLOWORLD_PREFIX}/lib"
    fi

    if test -d "${HELLOWORLD_PREFIX}/include" ; then
       HELLOWORLD_CPPFLAGS="-I${HELLOWORLD_PREFIX}/include -I${HELLOWORLD_PREFIX}/src"
    fi

    ac_HELLOWORLD_save_CPPFLAGS="$CPPFLAGS"
    ac_HELLOWORLD_save_LDFLAGS="$LDFLAGS"
    ac_HELLOWORLD_save_LIBS="$LIBS"

    CPPFLAGS="${HELLOWORLD_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${HELLOWORLD_LIBS} ${LDFLAGS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([helloworld/helloworld_version.h],[found_header=yes],[found_header=no])
    AC_LANG_POP([C++])

    #-----------------------
    # Minimum version check
    #----------------------

    min_helloworld_version=ifelse([$1], ,0.0.0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_helloworld_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_helloworld_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_helloworld_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for helloworld - version >= $min_helloworld_version)
        version_succeeded=no

	AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include "helloworld/helloworld_version.h"
            ]], [[
            #if HELLOWORLD_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (HELLOWORLD_MAJOR_VERSION >= $MAJOR_VER) && (HELLOWORLD_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (HELLOWORLD_MAJOR_VERSION >= $MAJOR_VER) && (HELLOWORLD_MINOR_VERSION >= $MINOR_VER) && (HELLOWORLD_MICRO_VERSION >= $MICRO_VER)
            /* I feel like chicken tonight, like chicken tonight? */
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])
	AC_LANG_POP([C++])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your HELLOWORLD version does not meet the minimum versioning
   requirements ($min_helloworld_version).  Please use --with-helloworld to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lhelloworld linkage])

    AC_LANG_PUSH([C++])

    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "helloworld/helloworld_version.h"],
                                   [HELLOWORLD::get_helloworld_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CPPFLAGS="$ac_HELLOWORLD_save_CPPFLAGS"
    LDFLAGS="$ac_HELLOWORLD_save_LDFLAGS"
    LIBS="$ac_HELLOWORLD_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
           if test "$found_library" = yes; then
              succeeded=yes
           fi
        fi
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([HELLOWORLD not found.  Try either --with-helloworld or setting HELLOWORLD_DIR.])
       else
          AC_MSG_NOTICE([optional HELLOWORLD library not found])
          HELLOWORLD_CPPFLAGS=""   # HELLOWORLD_CFLAGS empty on failure
          HELLOWORLD_LIBS=""     # HELLOWORLD_LIBS empty on failure
       fi
    else
        HAVE_HELLOWORLD=1
        AC_DEFINE(HAVE_HELLOWORLD,1,[Define if HELLOWORLD is available])
        AC_SUBST(HELLOWORLD_CPPFLAGS)
        AC_SUBST(HELLOWORLD_LIBS)
        AC_SUBST(HELLOWORLD_PREFIX)
    fi

    AC_SUBST(HAVE_HELLOWORLD)

# fi

AM_CONDITIONAL(HELLOWORLD_ENABLED,test x$HAVE_HELLOWORLD = x1)

])
