#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = metaphysicl_version

lib_LTLIBRARIES = libmetaphysicl.la

libmetaphysicl_la_LDFLAGS = -version-info 10:0:0

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
libmetaphysicl_la_SOURCES =

# utilities
libmetaphysicl_la_SOURCES += utilities/src/metaphysicl_version.C
libmetaphysicl_la_SOURCES += utilities/src/metaphysicl_exceptions.C

#----------------------------
# INCLUDES we want distributed
#----------------------------
includedir = $(prefix)/include/metaphysicl
include_HEADERS =

# core
include_HEADERS += core/include/metaphysicl/compare_types.h
include_HEADERS += core/include/metaphysicl/ct_set.h
include_HEADERS += core/include/metaphysicl/ct_types.h
include_HEADERS += core/include/metaphysicl/metaphysicl_common.h

# graphs
include_HEADERS += graphs/include/metaphysicl/physics.h

# numerics
include_HEADERS += numerics/include/metaphysicl/dualderivatives.h
include_HEADERS += numerics/include/metaphysicl/dualdynamicsparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/dualdynamicsparsenumberarray_decl.h
include_HEADERS += numerics/include/metaphysicl/dualdynamicsparsenumbervector.h
include_HEADERS += numerics/include/metaphysicl/dualdynamicsparsenumbervector_decl.h
include_HEADERS += numerics/include/metaphysicl/dualexpression.h
include_HEADERS += numerics/include/metaphysicl/dualnamedarray.h
include_HEADERS += numerics/include/metaphysicl/dualnumber_forward.h
include_HEADERS += numerics/include/metaphysicl/dualnumber.h
include_HEADERS += numerics/include/metaphysicl/dualnumber_decl.h
include_HEADERS += numerics/include/metaphysicl/nddualnumber.h
include_HEADERS += numerics/include/metaphysicl/dualnumber_surrogate.h
include_HEADERS += numerics/include/metaphysicl/dualnumber_surrogate_decl.h
include_HEADERS += numerics/include/metaphysicl/dualnumberarray.h
include_HEADERS += numerics/include/metaphysicl/dualnumbervector.h
include_HEADERS += numerics/include/metaphysicl/dualshadow.h
include_HEADERS += numerics/include/metaphysicl/dualshadowdynamicsparsearray.h
include_HEADERS += numerics/include/metaphysicl/dualshadowdynamicsparsevector.h
include_HEADERS += numerics/include/metaphysicl/dualshadowsparsestruct.h
include_HEADERS += numerics/include/metaphysicl/dualshadowsparsevector.h
include_HEADERS += numerics/include/metaphysicl/dualshadowvector.h
include_HEADERS += numerics/include/metaphysicl/dualsparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/dualsparsenumberstruct.h
include_HEADERS += numerics/include/metaphysicl/dualsparsenumbervector.h
include_HEADERS += numerics/include/metaphysicl/dynamicsparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/dynamicsparsenumberarray_decl.h
include_HEADERS += numerics/include/metaphysicl/dynamicsparsenumberbase.h
include_HEADERS += numerics/include/metaphysicl/dynamicsparsenumberbase_decl.h
include_HEADERS += numerics/include/metaphysicl/dynamicsparsenumbervector.h
include_HEADERS += numerics/include/metaphysicl/dynamicsparsenumbervector_decl.h
include_HEADERS += numerics/include/metaphysicl/namedindexarray.h
include_HEADERS += numerics/include/metaphysicl/numberarray.h
include_HEADERS += numerics/include/metaphysicl/numbervector.h
include_HEADERS += numerics/include/metaphysicl/raw_type.h
include_HEADERS += numerics/include/metaphysicl/shadownumber.h
include_HEADERS += numerics/include/metaphysicl/sparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/sparsenumberstruct.h
include_HEADERS += numerics/include/metaphysicl/sparsenumberutils.h
include_HEADERS += numerics/include/metaphysicl/sparsenumbervector.h
include_HEADERS += numerics/include/metaphysicl/dualsemidynamicsparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/dualsemidynamicsparsenumberarray_decl.h
include_HEADERS += numerics/include/metaphysicl/semidynamicsparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/semidynamicsparsenumberarray_decl.h
include_HEADERS += numerics/include/metaphysicl/parallel_dualnumber.h
include_HEADERS += numerics/include/metaphysicl/parallel_dynamicsparsenumberarray.h
include_HEADERS += numerics/include/metaphysicl/parallel_numberarray.h
include_HEADERS += numerics/include/metaphysicl/parallel_semidynamicsparsenumberarray.h

# utilities
include_HEADERS += utilities/include/metaphysicl/metaphysicl_asserts.h
include_HEADERS += utilities/include/metaphysicl/metaphysicl_cast.h
include_HEADERS += utilities/include/metaphysicl/metaphysicl_exceptions.h
include_HEADERS += utilities/include/metaphysicl/metaprogramming.h
include_HEADERS += utilities/include/metaphysicl/testable.h
include_HEADERS += utilities/include/metaphysicl/dynamic_std_array_wrapper.h
include_HEADERS += utilities/include/metaphysicl/parallel_dynamic_std_array_wrapper.h

# Needs to be builddir since this is generated by configure
include_HEADERS += $(top_builddir)/src/utilities/include/metaphysicl/metaphysicl_version.h
include_HEADERS += $(top_builddir)/src/utilities/include/metaphysicl/metaphysicl_config.h

# Version app
metaphysicl_version_SOURCES = apps/version.C
metaphysicl_version_LDADD = libmetaphysicl.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS  =
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/graphs/include
AM_CPPFLAGS += -I$(top_srcdir)/src/numerics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #metaphysicl_version.h

AM_CXXFLAGS = $(ACSM_ANY_WERROR_FLAG) $(ACSM_ANY_PARANOID_FLAGS)

CLEANFILES =

#--------------
# Code coverage
#--------------
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

if METAPHYSICL_LICENSESTAMPEXISTS
STAMPED_FILES = $(libmetaphysicl_la_SOURCES) $(include_HEADERS) $(metaphysicl_version_SOURCES)

BUILT_SOURCES = .license.stamp

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl -S $(top_srcdir)/src $(top_srcdir)/LICENSE $(STAMPED_FILES)
	echo 'updated source license headers' >$@

CLEANFILES += .license.stamp
endif
