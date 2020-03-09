AC_DEFUN([TIMPI_CONTROL_ARGS],
[
AC_ARG_VAR(TIMPI_DIR,[root directory of TIMPI installation])

dnl We dont default future_timpi_dir to an environment TIMPI_DIR because we want the user
dnl to be very explicit if they're using this configure capability
enable_future_timpi=no
AC_ARG_WITH(future-timpi-dir,
            [AS_HELP_STRING([--with-future-timpi-dir=DIR], [root directory of future TIMPI installation. Forces METAPHYSICL_HAVE_TIMPI if MPI installation found])],
            [
              future_timpi_dir=$withval
              if test "$future_timpi_dir" = yes; then
                AC_MSG_ERROR([You must specify a location for the future timpi directory])
              fi
              if test "$future_timpi_dir" != no; then
                enable_future_timpi=yes
                TIMPI_PREFIX="$future_timpi_dir"
              fi
            ])

dnl After this block enable_installed_timpi should be either yes or no
enable_installed_timpi=no
AC_ARG_WITH(timpi,
  [AS_HELP_STRING([--with-timpi@<:@=DIR@:>@],[root directory of TIMPI installation (default = TIMPI_DIR)])],
  [
    if test "$withval" != yes; then
        TIMPI_PREFIX=$withval
        dnl At this point with_timpi is either no or a directory path
        if test "$withval" != no; then
          enable_installed_timpi=yes
        fi
    fi
    if test "enable_installed_timpi" = yes && test "$enable_future_timpi" = yes; then
      AC_MSG_ERROR([You cannot pass both the --with-timpi and --with-future-timpi-dir configure options])
    fi
  ],
  [
    dnl assume a sensible default of --with-timpi=yes unless enable_future_timpi=yes
    AS_IF([test "$enable_future_timpi" = no],
          [
            if test "x${TIMPI_DIR}" != "x"; then
               enable_installed_timpi=yes
               TIMPI_PREFIX=${TIMPI_DIR}
            fi
          ],
          [enable_installed_timpi=no])
  ])
])
