AC_DEFUN([AC_HAVE_FEEXCEPT],
[
dnl Try to build a program with feenableexcept
dnl
dnl It's not safe to just AC_COMPILE_IFELSE, because that gives me a
dnl false positive on a conda environment with clang 12.0.1 - C
dnl allowed implicit function declarations until C99, and for some
dnl reason this error according to a 1999 standard is only a warning
dnl according to a 2021 compiler.

ac_FEEXCEPT_save_LIBS="$LIBS"

LIBS="$LIBS -lm"

AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([@%:@define _GNU_SOURCE
                              @%:@include <fenv.h>],
               [feenableexcept(FE_DIVBYZERO | FE_INVALID);])],
            [ac_cv_have_feenableexcept=yes],
            [ac_cv_have_feenableexcept=no])

AS_IF([test "x$ac_cv_have_feenableexcept" = "xyes"],
      [AC_DEFINE([HAVE_FEENABLEEXCEPT],[1],[define if the compiler supports feenableexcept])])

dnl Try to compile a program with fedisableexcept
AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([@%:@define _GNU_SOURCE
                              @%:@include <fenv.h>],
               [fedisableexcept(FE_DIVBYZERO | FE_INVALID);])],
            [ac_cv_have_fedisableexcept=yes],
            [ac_cv_have_fedisableexcept=no])

AS_IF([test "x$ac_cv_have_fedisableexcept" = "xyes"],
      [AC_DEFINE([HAVE_FEDISABLEEXCEPT],[1],[define if the compiler supports fedisableexcept])])

LIBS="$ac_FEEXCEPT_save_LIBS"
])
