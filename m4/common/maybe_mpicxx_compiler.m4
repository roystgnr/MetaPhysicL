# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([MAYBE_MPICXX_COMPILER],
[
  AC_REQUIRE([ACSM_COMPILER_CONTROL_ARGS])
  AC_REQUIRE([ACSM_SCRAPE_PETSC_CONFIGURE])

  # --------------------------------------------------------------
  # look for a decent C++ compiler or honor --with-cxx=...
  CXX_TRY_LIST="g++ icpc icc pgCC c++"

  AC_ARG_WITH([cxx],
              AS_HELP_STRING([--with-cxx=CXX], [C++ compiler to use]),
              [CXX="$withval"],
              [])

  dnl If we already have CXX either from --with-cxx or from the environment, then there's no sense going
  dnl any further. Moreover if we are not enabling mpi then we don't have to query for mpi compilers
  dnl or for a compiler from PETSc
  AS_IF([test -z "$CXX" && test "$enablempi" != no],
        [
          dnl Did we get --with-mpi=DIR or was there a MPI_HOME or MPIHOME set?
          AS_IF([test x"$MPI" != x],
                [
                  dnl Inspect $MPI/bin
                  AS_IF([test -d "$MPI/bin"],
                        [
                          AC_CHECK_PROGS(LOCAL_CXX, [mpicxx mpiCC mpicc], [], ["$MPI/bin"])
                          AS_IF([test -z "$LOCAL_CXX"],
                                [AS_ECHO(["None of the wrappers we look for exist in $MPI/bin. We will not try to use mpi compiler wrappers"])],
                                [MPI_USING_WRAPPERS=1;CXX="$MPI/bin/$LOCAL_CXX"])
                        ],
                        [AS_ECHO(["An MPI directory was specified, but $MPI/bin does not exist. We will not try to use mpi compiler wrappers"])])
                ],
                [
                  dnl No MPI directory specified. If we have PETSc, let's try to snoop some
                  dnl information from there. We'll use this information further below and in
                  dnl mpi.m4
                  AS_IF([test x"$PETSC_HAVE_MPI" = x1 && test x"$PETSC_CXX" != x],
                        [],
                        dnl PETSc doesn't define a CXX so we'll just try to pull one from the environment
                        [CXX_TRY_LIST="mpicxx mpiCC mpicc $CXX_TRY_LIST"])
                ]) dnl AS_IF([test x"$MPI" != x])
        ])

  dnl See whether we are using PETSC_CXX. Unfortunately PETSC_CXX may
  dnl be prefixed with a PATH, and its not straightforward to strip it off.
  dnl If CXX is not set AC_PROG_CXX will call AC_CHECK_TOOLS which will prefix
  dnl every argument in CXX_TRY_LIST with values in $PATH, so we will
  dnl not find something like PATH/PETSC_PREFIX/mpicxx. The solution
  dnl then is just to set CXX to PETSC_CXX so that AC_CHECK_TOOLS
  dnl never gets called
  AS_IF([test -z "$CXX" && test x"$PETSC_HAVE_MPI" = x1 && test x"$PETSC_CXX" != x],
        [CXX="$PETSC_CXX"])

  dnl If we still don't have a CXX set then we will try to pick one up from CXX_TRY_LIST
  AC_PROG_CXX([$CXX_TRY_LIST])
])
