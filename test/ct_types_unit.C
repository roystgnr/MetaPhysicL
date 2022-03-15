#include <metaphysicl/dynamicsparsenumberarray.h>
#include <metaphysicl/dynamicsparsenumbervector.h>
#include <metaphysicl/semidynamicsparsenumberarray.h>
#include <metaphysicl/numberarray.h>
#include <metaphysicl/numbervector.h>
#include <metaphysicl/dualnumber.h>
#include <metaphysicl/ct_types.h>

#define METAPHYSICL_UNIT_ASSERT(expr)                                                              \
  if (!(expr))                                                                                     \
  metaphysicl_error()

#ifndef METAPHYSICL_COMMA
#define METAPHYSICL_COMMA ,
#endif

using namespace MetaPhysicL;

int
main()
{
  NumberArray<1, float> f_na;
  DynamicSparseNumberArray<float, unsigned int> f_dsna;
  SemiDynamicSparseNumberArray<float, unsigned int, NWrapper<1>> f_sdsna;
  std::array<float, 1> f_a;
  std::vector<float> f_v;

  NumberArray<1, float> d_na;
  DynamicSparseNumberArray<float, unsigned int> d_dsna;
  SemiDynamicSparseNumberArray<float, unsigned int, NWrapper<1>> d_sdsna;
  std::array<double, 1> d_a;
  std::vector<double> d_v;

  // Test container replacement
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<decltype(f_na) METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA decltype(d_na)>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<decltype(f_dsna) METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA decltype(d_dsna)>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<decltype(f_sdsna) METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA decltype(d_sdsna)>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<decltype(f_a) METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA decltype(d_a)>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<decltype(f_v) METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA decltype(d_v)>::value);

  // Test atomic/nuclear algebraic type replacement
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<char METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<signed char METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<unsigned char METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<short int METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<unsigned short int METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<int METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<unsigned int METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<long METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<long long METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<unsigned long METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<unsigned long long METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<float METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<long double METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<DualNumber<double> METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(std::is_same<typename ReplaceAlgebraicType<
                              DynamicSparseNumberVector<double METAPHYSICL_COMMA unsigned int>
                                  METAPHYSICL_COMMA double>::type METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<
          NumberVector<1 METAPHYSICL_COMMA double> METAPHYSICL_COMMA double>::type
                       METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ReplaceAlgebraicType<double METAPHYSICL_COMMA float>::type
                       METAPHYSICL_COMMA float>::value);

  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ValueType<std::vector<double>>::type METAPHYSICL_COMMA double>::value);
  METAPHYSICL_UNIT_ASSERT(
      std::is_same<typename ValueType<std::array<double METAPHYSICL_COMMA 1>>::type
                       METAPHYSICL_COMMA double>::value);
}
