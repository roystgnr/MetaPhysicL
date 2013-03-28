BUILT_SOURCES   = .license.stamp

check_PROGRAMS  =
check_PROGRAMS += derivs_unit
check_PROGRAMS += identities_unit
check_PROGRAMS += main_unit
check_PROGRAMS += navier_unit
check_PROGRAMS += pde_unit
check_PROGRAMS += testheaders_unit
check_PROGRAMS += testopt_unit

if CXX11_ENABLED
  check_PROGRAMS += physics_unit
endif

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/graphs/include
AM_CPPFLAGS += -I$(top_srcdir)/src/numerics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #metaphysicl_version.h

AM_LDFLAGS = $(top_builddir)/src/libmetaphysicl.la

# Sources for these tests
derivs_unit_SOURCES = derivs_unit.C
identities_unit_SOURCES = identities_unit.C
main_unit_SOURCES = main_unit.C
navier_unit_SOURCES = navier_unit.C
pde_unit_SOURCES = pde_unit.C
physics_unit_SOURCES = physics_unit.C
testheaders_unit_SOURCES = testheaders_unit.C
testopt_unit_SOURCES = testopt_unit.C

#Define tests to actually be run
TESTS  =
TESTS += derivs_unit
TESTS += identities_unit
TESTS += main_unit
TESTS += navier_unit
TESTS += pde_unit
TESTS += testheaders_unit
TESTS += testopt_unit

if CXX11_ENABLED
  TESTS += physics_unit
endif

#-------------
# MASA support
#-------------
if MASA_ENABLED
  AM_CPPFLAGS += $(MASA_CXXFLAGS)
  AM_LDFLAGS  += $(MASA_LIBS)
endif

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
