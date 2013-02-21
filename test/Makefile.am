BUILT_SOURCES   = .license.stamp

check_PROGRAMS  =
check_PROGRAMS += core_name_unit
check_PROGRAMS += input_file_unit

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #helloworld_version.h

AM_LDFLAGS = $(top_builddir)/src/libhelloworld.la

# Sources for these tests
core_name_unit_SOURCES = core_name_unit.C
input_file_unit_SOURCES = input_file_unit.C

#Define tests to actually be run
TESTS  =
TESTS += core_name_unit
TESTS += input_file_unit.sh

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
