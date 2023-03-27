#include <metaphysicl/metaphysicl_config.h>

#if defined(METAPHYSICL_HAVE_TIMPI) && defined(METAPHYSICL_HAVE_TIMPI_LIB)

#include "metaphysicl/metaphysicl_exceptions.h"
#include <metaphysicl/parallel_dynamicsparsenumberarray.h>
#include <metaphysicl/parallel_dualnumber.h>
#include <metaphysicl/parallel_numberarray.h>
#include <metaphysicl/parallel_semidynamicsparsenumberarray.h>

#include <timpi/communicator.h>
#include <timpi/parallel_implementation.h>
#include <timpi/parallel_sync.h>
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
  if (fixed_size)
    in.derivatives().insert(std::min(my_rank,maxarraysize-1)) = 1.;
  else
    in.derivatives().insert(my_rank) = 1.;

  TestCommWorld->allgather(in, vals);

  const unsigned int comm_size = TestCommWorld->size();
  const unsigned int vec_size = vals.size();

  METAPHYSICL_UNIT_ASSERT(comm_size == vec_size);

  for (unsigned int i = 0; i < vec_size; ++i)
  {
    const auto & dn = vals[i];
    if (fixed_size)
      METAPHYSICL_UNIT_ASSERT(dn.derivatives().size() == maxarraysize);
    else
      METAPHYSICL_UNIT_ASSERT(dn.derivatives().size() == 1);
    METAPHYSICL_UNIT_FP_ASSERT(dn.value(), double(i), TOLERANCE);

    if (fixed_size)
      {
        const unsigned int one_j = std::min(i,maxarraysize-1);
        METAPHYSICL_UNIT_FP_ASSERT(dn.derivatives()[one_j], double(1), TOLERANCE);
        for (unsigned int j = 0; j < maxarraysize; ++j)
          if (j != one_j)
            METAPHYSICL_UNIT_FP_ASSERT(dn.derivatives()[j], double(0), TOLERANCE);
      }
    else
      METAPHYSICL_UNIT_FP_ASSERT(dn.derivatives()[i], double(1), TOLERANCE);
  }
}


template <typename D, bool asd>
void
testPackedAllGather()
{
  typedef DualNumber<double, D, asd> DualReal;

  typedef std::map<int, std::vector<DualReal>> Container;
  Container vals;
  const unsigned int my_rank = TestCommWorld->rank();

  // Initialize values
  if (my_rank == 0)
    {
      vals[0] = { DualReal(1), DualReal(2), DualReal(3) };
      vals[1] = { DualReal(4), DualReal(5) };
    }
  else if (my_rank == 1)
    {
      vals[2] = { DualReal(6) };
      vals[3] = { DualReal(7), DualReal(8) };
    }

  std::vector<Container> all_vals;

  TestCommWorld->allgather(vals, all_vals);

  const std::size_t comm_size = TestCommWorld->size();
  const std::size_t vec_size = all_vals.size();

  METAPHYSICL_UNIT_ASSERT(comm_size == vec_size);

  METAPHYSICL_UNIT_ASSERT(all_vals[0][0].size() == 3);
  METAPHYSICL_UNIT_ASSERT(all_vals[0][1].size() == 2);

  METAPHYSICL_UNIT_FP_ASSERT(all_vals[0][0][0], double(1), TOLERANCE);
  METAPHYSICL_UNIT_FP_ASSERT(all_vals[0][0][1], double(2), TOLERANCE);
  METAPHYSICL_UNIT_FP_ASSERT(all_vals[0][0][2], double(3), TOLERANCE);
  METAPHYSICL_UNIT_FP_ASSERT(all_vals[0][1][0], double(4), TOLERANCE);
  METAPHYSICL_UNIT_FP_ASSERT(all_vals[0][1][1], double(5), TOLERANCE);

  if (vec_size > 1)
    {
      METAPHYSICL_UNIT_ASSERT(all_vals[1][2].size() == 1);
      METAPHYSICL_UNIT_ASSERT(all_vals[1][3].size() == 2);

      METAPHYSICL_UNIT_FP_ASSERT(all_vals[1][2][0], double(6), TOLERANCE);
      METAPHYSICL_UNIT_FP_ASSERT(all_vals[1][3][0], double(7), TOLERANCE);
      METAPHYSICL_UNIT_FP_ASSERT(all_vals[1][3][1], double(8), TOLERANCE);
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
  if (!fixed_size || my_rank < maxarraysize)
    dr.derivatives().insert(my_rank) = (my_rank+1);

  TestCommWorld->broadcast(dr);

  METAPHYSICL_UNIT_ASSERT(dr.value() == 4.0);
  if (fixed_size)
    {
      METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == maxarraysize);
      for (std::size_t i = 1; i < maxarraysize; ++i)
        METAPHYSICL_UNIT_ASSERT(dr.derivatives()[i] == 0.0);
    }
  else
    {
      METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == 1);
    }
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
  C c {0};
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

  METAPHYSICL_UNIT_ASSERT(dr.value() == 4.0*comm_size + comm_size*(comm_size-1)/2);
  if (fixed_size)
    METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == maxarraysize);
  else
    METAPHYSICL_UNIT_ASSERT(dr.derivatives().size() == full_size+1);

  METAPHYSICL_UNIT_ASSERT(dr.derivatives()[0] == 1.0);
  for (unsigned int p = 1; p != full_size; ++p)
    METAPHYSICL_UNIT_ASSERT(dr.derivatives()[p] == 2*p+2);
  METAPHYSICL_UNIT_ASSERT(dr.derivatives()[full_size] == full_size+1);
}


template <typename D, bool asd>
void
testDualContainerSync()
{
  typedef DualNumber<double, D, asd> DualReal;
  typedef DualReal Datum;

  const unsigned int my_rank = TestCommWorld->rank();
  const unsigned int comm_size = TestCommWorld->size();
  const unsigned int full_size = std::min(maxarraysize-1, comm_size);

  std::unordered_map<processor_id_type, std::vector<Datum>> push_data;

  for (unsigned int i=0; i!=5; ++i)
    {
      const processor_id_type p = (my_rank+(i*i))%comm_size;
      DualReal dr = my_rank;
      if (my_rank < full_size)
        {
          dr.derivatives().insert(my_rank) = (my_rank+1);
          dr.derivatives().insert(my_rank+1) = (p+1);
        }
      push_data[p].push_back(dr);
    }

  auto action_functor = [my_rank, full_size]
    (const processor_id_type src_p,
     const std::vector<Datum> & received)
    {
      for (const auto & r : received)
        {
          METAPHYSICL_UNIT_ASSERT(r.value() == src_p);
          for (unsigned int p = 0; p != full_size; ++p)
            {
              const auto & d = r.derivatives()[p];
              if (p == src_p)
                METAPHYSICL_UNIT_ASSERT(d == src_p+1);
              else if (p == src_p+1)
                METAPHYSICL_UNIT_ASSERT(d == my_rank+1);
              else
                METAPHYSICL_UNIT_ASSERT(d == 0);
            }
        }
    };

  TIMPI::push_parallel_vector_data(*TestCommWorld, push_data, action_functor);
}


template <typename D, bool asd>
void
testDualContainerTupleSync()
{
  typedef DualNumber<double, D, asd> DualReal;
  typedef std::tuple<unsigned int, DualReal, double> Datum;

  const unsigned int my_rank = TestCommWorld->rank();
  const unsigned int comm_size = TestCommWorld->size();
  const unsigned int full_size = std::min(maxarraysize-1, comm_size);

  std::unordered_map<processor_id_type, std::vector<Datum>> push_data;

  for (unsigned int i=0; i!=5; ++i)
    {
      const processor_id_type p = (my_rank+(i*i))%comm_size;
      DualReal dr = my_rank;
      if (my_rank < full_size)
        {
          dr.derivatives().insert(my_rank) = (my_rank+1);
          dr.derivatives().insert(my_rank+1) = (p+1);
        }
      push_data[p].emplace_back(my_rank+2, dr, p);
    }

  auto action_functor = [my_rank, full_size]
    (const processor_id_type src_p,
     const std::vector<Datum> & received)
    {
      for (const auto & t : received)
        {
          METAPHYSICL_UNIT_ASSERT(std::get<0>(t) == src_p+2);
          const DualReal & r = std::get<1>(t);
          METAPHYSICL_UNIT_ASSERT(r.value() == src_p);
          for (unsigned int p = 0; p != full_size; ++p)
            {
              const auto & d = r.derivatives()[p];
              if (p == src_p)
                METAPHYSICL_UNIT_ASSERT(d == src_p+1);
              else if (p == src_p+1)
                METAPHYSICL_UNIT_ASSERT(d == my_rank+1);
              else
                METAPHYSICL_UNIT_ASSERT(d == 0);
            }
          METAPHYSICL_UNIT_ASSERT(std::get<2>(t) == my_rank);
        }
    };

  TIMPI::push_parallel_vector_data(*TestCommWorld, push_data, action_functor);
}


int
main(int argc, const char * const * argv)
{
  MetaPhysicL::enableFPE(true);

  TIMPI::TIMPIInit init(argc, argv);
  TestCommWorld = &init.comm();

  testDualSum<double, true>();
  testDualSum<double, false>();

  testBroadcast<DynamicSparseNumberArray<double, unsigned int>, true>();
  testBroadcast<DynamicSparseNumberArray<double, unsigned int>, false>();
  testBroadcast<DynamicSparseNumberArray<double, unsigned long long>, true>();
  testBroadcast<DynamicSparseNumberArray<double, unsigned long long>, false>();
  testBroadcast<NumberArray<maxarraysize, double>, true>(true);
  testBroadcast<NumberArray<maxarraysize, double>, false>(true);
  testBroadcast<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testBroadcast<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();
  testBroadcast<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, true>();
  testBroadcast<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, false>();

  testContainerAllGather<DynamicSparseNumberArray<double, unsigned int>, true>();
  testContainerAllGather<DynamicSparseNumberArray<double, unsigned int>, false>();
  testContainerAllGather<DynamicSparseNumberArray<double, unsigned long long>, true>();
  testContainerAllGather<DynamicSparseNumberArray<double, unsigned long long>, false>();
  testContainerAllGather<NumberArray<maxarraysize, double>, true>(true);
  testContainerAllGather<NumberArray<maxarraysize, double>, false>(true);
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, true>();
  testContainerAllGather<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, false>();

  testPackedAllGather<DynamicSparseNumberArray<double, unsigned int>, true>();
  testPackedAllGather<DynamicSparseNumberArray<double, unsigned int>, false>();
  testPackedAllGather<DynamicSparseNumberArray<double, unsigned long long>, true>();
  testPackedAllGather<DynamicSparseNumberArray<double, unsigned long long>, false>();
  testPackedAllGather<NumberArray<maxarraysize, double>, true>();
  testPackedAllGather<NumberArray<maxarraysize, double>, false>();
  testPackedAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testPackedAllGather<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();
  testPackedAllGather<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, true>();
  testPackedAllGather<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, false>();

  testContainerSum<NumberArray<maxarraysize, double>>(true);
  testContainerSum<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>>();
  testContainerSum<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>>();

/*
  // These rely on reduction support for packed-range types!
  testDualContainerSum<DynamicSparseNumberArray<double, unsigned int>, true>();
  testDualContainerSum<DynamicSparseNumberArray<double, unsigned int>, false>();
*/
  testDualContainerSum<NumberArray<maxarraysize, double>, true>(true);
  testDualContainerSum<NumberArray<maxarraysize, double>, false>(true);
  testDualContainerSum<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testDualContainerSum<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();
  testDualContainerSum<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, true>();
  testDualContainerSum<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, false>();

  testDualContainerSync<DynamicSparseNumberArray<double, unsigned int>, true>();
  testDualContainerSync<DynamicSparseNumberArray<double, unsigned int>, false>();
  testDualContainerSync<DynamicSparseNumberArray<double, unsigned long long>, true>();
  testDualContainerSync<DynamicSparseNumberArray<double, unsigned long long>, false>();
  testDualContainerSync<NumberArray<maxarraysize, double>, true>();
  testDualContainerSync<NumberArray<maxarraysize, double>, false>();
  testDualContainerSync<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testDualContainerSync<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();
  testDualContainerSync<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, true>();
  testDualContainerSync<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, false>();

  testDualContainerTupleSync<DynamicSparseNumberArray<double, unsigned int>, true>();
  testDualContainerTupleSync<DynamicSparseNumberArray<double, unsigned int>, false>();
  testDualContainerTupleSync<DynamicSparseNumberArray<double, unsigned long long>, true>();
  testDualContainerTupleSync<DynamicSparseNumberArray<double, unsigned long long>, false>();
  testDualContainerTupleSync<NumberArray<maxarraysize, double>, true>();
  testDualContainerTupleSync<NumberArray<maxarraysize, double>, false>();
  testDualContainerTupleSync<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, true>();
  testDualContainerTupleSync<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>, false>();
  testDualContainerTupleSync<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, true>();
  testDualContainerTupleSync<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>, false>();

  testStandardTypeAssignment<NumberArray<maxarraysize, double>>();
  testStandardTypeAssignment<SemiDynamicSparseNumberArray<double, unsigned int, NWrapper<maxarraysize>>>();
  testStandardTypeAssignment<SemiDynamicSparseNumberArray<double, unsigned long long, NWrapper<maxarraysize>>>();
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
