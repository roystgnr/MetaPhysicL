#include <metaphysicl/metaphysicl_config.h>

#if defined(METAPHYSICL_HAVE_TIMPI) && defined(METAPHYSICL_HAVE_TIMPI_LIB)

#include <metaphysicl/parallel_dynamicsparsenumberarray.h>
#include <metaphysicl/parallel_dualnumber.h>
#include <metaphysicl/parallel_semidynamicsparsenumberarray.h>

#include <timpi/communicator.h>
#include <timpi/parallel_implementation.h>
#include <timpi/timpi_init.h>

#include <vector>

#define METAPHYSICL_UNIT_ASSERT(expr)                                                              \
  if (!(expr))                                                                                     \
  metaphysicl_error()

#define METAPHYSICL_UNIT_FP_ASSERT(test, true_value, tolerance)                                    \
  if (!(test < true_value + tolerance && test > true_value - tolerance))                           \
  metaphysicl_error()

#define TOLERANCE 1e-12

using namespace TIMPI;
using namespace MetaPhysicL;

Communicator * TestCommWorld;

template <typename D, bool asd>
void
testContainerAllGather()
{
  typedef DualNumber<double, D, asd> DualReal;

  std::vector<DualReal> vals;
  const unsigned int my_rank = TestCommWorld->rank();

  // Initialize value
  DualReal in = my_rank;
  // Initialize derivative
  in.derivatives().insert(my_rank) = 1.;

  TestCommWorld->allgather(in, vals);

  const std::size_t comm_size = TestCommWorld->size();
  const std::size_t vec_size = vals.size();

  METAPHYSICL_UNIT_ASSERT(comm_size == vec_size);

  for (std::size_t i = 0; i < vec_size; ++i)
  {
    const auto & dn = vals[i];
    METAPHYSICL_UNIT_ASSERT(dn.derivatives().size() == 1);
    METAPHYSICL_UNIT_FP_ASSERT(dn.value(), double(i), TOLERANCE);
    METAPHYSICL_UNIT_FP_ASSERT(dn.derivatives()[i], double(1), TOLERANCE);
  }
}

template <typename NonBuiltin>
void
testStandardTypeAssignment()
{
  NonBuiltin ex{};
  StandardType<NonBuiltin> a(&ex);
  StandardType<NonBuiltin> b(&ex);
  a = b;
}

int
main(int argc, const char * const * argv)
{
  TIMPI::TIMPIInit init(argc, argv);
  TestCommWorld = &init.comm();

  testContainerAllGather<DynamicSparseNumberArray<double, unsigned int>, true>();
  testContainerAllGather<DynamicSparseNumberArray<double, unsigned int>, false>();
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<50>>, true>();
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<50>>, false>();

  testStandardTypeAssignment<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<50>>>();
  testStandardTypeAssignment<DynamicStdArrayWrapper<double, NWrapper<50>>>();
  testStandardTypeAssignment<DualNumber<double>>();

  return 0;
}

#else

int
main()
{
  return 0;
}

#endif // METAPHYSICL_HAVE_TIMPI && METAPHYSICL_HAVE_TIMPI_LIB
