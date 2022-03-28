#include <metaphysicl/metaphysicl_config.h>

#if defined(METAPHYSICL_HAVE_TIMPI) && defined(METAPHYSICL_HAVE_TIMPI_LIB)

#include <metaphysicl/parallel_dynamicsparsenumberarray.h>
#include <metaphysicl/parallel_dualnumber.h>
#include <metaphysicl/parallel_numberarray.h>
#include <metaphysicl/parallel_semidynamicsparsenumberarray.h>

#include <timpi/communicator.h>
#include <timpi/parallel_implementation.h>
#include <timpi/timpi_init.h>

#include <vector>

#define METAPHYSICL_UNIT_ASSERT(expr)                                                              \
  do {                                                                                             \
    if (!(expr))                                                                                   \
      metaphysicl_error();                                                                         \
  } while (0)

#define METAPHYSICL_UNIT_FP_ASSERT(test, true_value, tolerance)                                    \
  do {                                                                                             \
    if (!(test < true_value + tolerance && test > true_value - tolerance))                         \
      metaphysicl_error();                                                                         \
  } while (0)

#define TOLERANCE 1e-12

using namespace TIMPI;
using namespace MetaPhysicL;

Communicator * TestCommWorld;

constexpr unsigned int maxarraysize = 50;

template <typename D, bool asd>
void
testContainerAllGather(bool fixed_size = false)
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
    if (fixed_size)
      METAPHYSICL_UNIT_ASSERT(dn.derivatives().size() == maxarraysize);
    else
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


template <typename D, bool asd>
void
testBroadcast(bool fixed_size = false)
{
  typedef DualNumber<double, D, asd> DualReal;

  const unsigned int my_rank = TestCommWorld->rank();

  // Initialize value
  DualReal dr = my_rank+4;
  // Initialize derivative
  dr.derivatives().insert(my_rank) = (my_rank+1);

  TestCommWorld->broadcast(dr);

  METAPHYSICL_UNIT_ASSERT(dr.value() == 4.0);
  if (fixed_size)
    METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == maxarraysize);
  else
    METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == 1);
  METAPHYSICL_UNIT_ASSERT(dr.derivatives()[0] == 1.0);
}


template <typename D, bool asd>
void
testDualSum()
{
  typedef DualNumber<double, D, asd> DualReal;

  const unsigned int my_rank = TestCommWorld->rank();
  const unsigned int comm_size = TestCommWorld->size();

  // Initialize value
  DualReal dr = my_rank+4;
  // Initialize derivative
  dr.derivatives() = my_rank+2;

  TestCommWorld->sum(dr);

  METAPHYSICL_UNIT_ASSERT(dr.value() == 4.0*comm_size + comm_size*(comm_size-1)/2);
  METAPHYSICL_UNIT_ASSERT(dr.derivatives() == 2.0*comm_size + comm_size*(comm_size-1)/2);
}


template <typename C>
void
testContainerSum(bool fixed_size = false)
{
  const unsigned int my_rank = TestCommWorld->rank();
  const unsigned int comm_size = TestCommWorld->size();
  const unsigned int full_size = std::min(maxarraysize-1, comm_size);

  // Initialize values
  C c;
  if (my_rank < full_size)
    {
      c.insert(my_rank) = (my_rank+1);
      c.insert(my_rank+1) = (my_rank+2);
    }

  TestCommWorld->sum(c);

  if (fixed_size)
    METAPHYSICL_UNIT_ASSERT(c.size() == maxarraysize);
  else
    METAPHYSICL_UNIT_ASSERT(c.size() == full_size+1);

  METAPHYSICL_UNIT_ASSERT(c[0] == 1.0);
  for (unsigned int p = 1; p != full_size; ++p)
    METAPHYSICL_UNIT_ASSERT(c[p] == 2*p+2);
  METAPHYSICL_UNIT_ASSERT(c[full_size] == full_size+1);
}


template <typename D, bool asd>
void
testDualContainerSum(bool fixed_size = false)
{
  typedef DualNumber<double, D, asd> DualReal;

  const unsigned int my_rank = TestCommWorld->rank();
  const unsigned int comm_size = TestCommWorld->size();
  const unsigned int full_size = std::min(maxarraysize-1, comm_size);

  // Initialize value
  DualReal dr = my_rank+4;
  // Initialize derivative
  if (my_rank < full_size)
    {
      dr.derivatives().insert(my_rank) = (my_rank+1);
      dr.derivatives().insert(my_rank+1) = (my_rank+2);
    }

  TestCommWorld->sum(dr);

  METAPHYSICL_UNIT_ASSERT(dr.value() == 4.0*full_size + full_size*(full_size-1)/2);
  if (fixed_size)
    METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == maxarraysize);
  else
    METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == full_size+1);

  METAPHYSICL_UNIT_ASSERT(dr.derivatives()[0] == 1.0);
  for (unsigned int p = 1; p != full_size; ++p)
    METAPHYSICL_UNIT_ASSERT(dr.derivatives()[p] == 2*p+2);
  METAPHYSICL_UNIT_ASSERT(dr.derivatives()[full_size] == full_size+1);
}


int
main(int argc, const char * const * argv)
{
  TIMPI::TIMPIInit init(argc, argv);
  TestCommWorld = &init.comm();

  testDualSum<double, true>();
  testDualSum<double, false>();

  testBroadcast<DynamicSparseNumberArray<double, unsigned int>, true>();
  testBroadcast<DynamicSparseNumberArray<double, unsigned int>, false>();
  testBroadcast<NumberArray<maxarraysize, double>, true>(true);
  testBroadcast<NumberArray<maxarraysize, double>, false>(true);
  testBroadcast<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testBroadcast<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();

  testContainerAllGather<DynamicSparseNumberArray<double, unsigned int>, true>();
  testContainerAllGather<DynamicSparseNumberArray<double, unsigned int>, false>();
  testContainerAllGather<NumberArray<maxarraysize, double>, true>(true);
  testContainerAllGather<NumberArray<maxarraysize, double>, false>(true);
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();

  testContainerSum<NumberArray<maxarraysize, double>>(true);
  testContainerSum<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>>();

/*
  // These rely on reduction support for packed-range types!
  testDualContainerSum<DynamicSparseNumberArray<double, unsigned int>, true>();
  testDualContainerSum<DynamicSparseNumberArray<double, unsigned int>, false>();
*/
  testDualContainerSum<NumberArray<maxarraysize, double>, true>(true);
  testDualContainerSum<NumberArray<maxarraysize, double>, false>(true);
  testDualContainerSum<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testDualContainerSum<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();

  testStandardTypeAssignment<NumberArray<maxarraysize, double>>();
  testStandardTypeAssignment<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>>();
  testStandardTypeAssignment<DynamicStdArrayWrapper<double, NWrapper<maxarraysize>>>();
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
