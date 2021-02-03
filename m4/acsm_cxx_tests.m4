dnl ----------------------------------------------------------------
dnl Test callback for use with ACSM.
dnl ----------------------------------------------------------------

dnl Currently we don't bother with extra tests ourselves.
dnl
dnl This is in the ACSM namespace because it is a "callback" from an
dnl ACSM macro.

AC_DEFUN([ACSM_TEST_CXX_ALL],
  [
    # Roll the dice!
    have_cxx_all=yes
  ])
