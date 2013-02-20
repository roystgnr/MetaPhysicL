BUILT_SOURCES   = .license.stamp

check_PROGRAMS  =
check_PROGRAMS += chem_mixture_unit
check_PROGRAMS += kinetics_unit

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h

AM_LDFLAGS = $(top_builddir)/src/libantioch.la

# Sources for these tests
chem_mixture_unit_SOURCES = chem_mixture_unit.C
kinetics_unit_SOURCES = kinetics_unit.C

#Define tests to actually be run
TESTS  =
TESTS += chem_mixture_unit
TESTS += kinetics_unit_air_5sp.sh

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
