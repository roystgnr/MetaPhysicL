
// MetaPhysicL
#include "metaphysicl/sparsenumbervector.h"
#include "metaphysicl/dualnamedarray.h"

#include "metaphysicl_config.h"

// VexCL
#ifdef METAPHYSICL_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// C++
#include <iostream>

using namespace MetaPhysicL;

int main(void)
{
  typedef
    NamedIndexArray
      <double,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<3>::type> >
    indexed_by_three;

  typedef
    DualExpression<indexed_by_three, indexed_by_three> dual_three;

  dual_three test_val;
  test_val.value().raw_data() = 0.5;
  test_val.derivatives().raw_data() = 1;
  test_val.value().raw_sizes().get<3>() = 1;
  test_val.derivatives().raw_sizes().get<3>() = 1;

#ifdef METAPHYSICL_HAVE_VEXCL
  vex::Context ctx (vex::Filter::Env && vex::Filter::Count(1));
  std::cout << ctx << std::endl;
      

  typedef
    NamedIndexArray
      <vex::vector<double>,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<2>::type> >
    vex_indexed_by_two;

  typedef
    DualExpression<vex_indexed_by_two, vex_indexed_by_two> dual_two;

  typedef
    NamedIndexArray
      <vex::vector<double>,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<1>::type> >
    vex_indexed_by_one;

  typedef
    DualExpression<vex_indexed_by_one, vex_indexed_by_one> dual_one;

  vex_indexed_by_one test_one_val(vex::vector<double>(ctx, 5),0);
  test_one_val.raw_sizes().template get<1>() = 5;

  vex_indexed_by_one test_one_deriv(vex::vector<double>(ctx, 5),0);
  test_one_deriv.raw_sizes().template get<1>() = 5;

  dual_one test_one(test_one_val, test_one_deriv);

  vex_indexed_by_two test_two_val(vex::vector<double>(ctx, 3),0);
  test_two_val.raw_sizes().template get<2>() = 3;

  vex_indexed_by_two test_two_deriv(vex::vector<double>(ctx, 3),0);
  test_two_deriv.raw_sizes().template get<2>() = 3;

  dual_two test_two(test_two_val, test_two_deriv);

  test_one.value().raw_data()[2] = 7;
  test_two.value().raw_data()[1] = 2;

  auto test_three_val_val = test_one_val * test_two_val;
  auto test_three_deriv_val = test_one_deriv * test_two_val;

  // Compilation Failure here:
  // auto test_three_a = test_one * test_two_val;

  // auto test_three = test_one * test_two;

/*
  if (test_three.value().raw_sizes().template get<1>() != 5)
    return 1;

  if (test_three.value().raw_sizes().template get<2>() != 3)
    return 1;

  vex::vector<double> test_output = test_three.value().raw_data();

  if (test_output[7] != 14)
    return 1;
*/
#endif

  return 0;
}
