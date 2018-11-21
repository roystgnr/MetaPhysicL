#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"

#include "math_structs.h"

#include <iostream>

using namespace MetaPhysicL;

#define nd_derivs_expect_near(double1, double2, tolerance)                                         \
  {                                                                                                \
    int new_returnval = std::abs(double1 - double2) > tolerance;                                   \
    if (new_returnval)                                                                             \
      std::cerr << "Failed test at line " << __LINE__ << std::endl;                                \
    returnval = returnval || new_returnval;                                                        \
  }

int
main()
{
  NDDualNumber<double, NumberArray<2, double>> scalar_ad_prop;
  NDDualNumber<VectorValue<double>, NumberArray<2, VectorValue<double>>> vector_ad_prop;
  NDDualNumber<TensorValue<double>, NumberArray<2, TensorValue<double>>> tensor_ad_prop;

  double scalar_reg_prop;
  VectorValue<double> vector_reg_prop;
  TensorValue<double> tensor_reg_prop;

  int returnval = 0;
  double tol = 1e-8;

  scalar_ad_prop = 2.;
  scalar_ad_prop.derivatives()[0] = 1.;
  scalar_ad_prop.derivatives()[1] = 2.;
  scalar_ad_prop *= scalar_ad_prop;
  nd_derivs_expect_near(scalar_ad_prop.value(), 4., tol);
  nd_derivs_expect_near(scalar_ad_prop.derivatives()[0], 4., tol);
  nd_derivs_expect_near(scalar_ad_prop.derivatives()[1], 8., tol);

  scalar_ad_prop = scalar_ad_prop * scalar_ad_prop;
  nd_derivs_expect_near(scalar_ad_prop.value(), 16., tol);
  nd_derivs_expect_near(scalar_ad_prop.derivatives()[0], 32., tol);
  nd_derivs_expect_near(scalar_ad_prop.derivatives()[1], 64., tol);

  vector_ad_prop = scalar_ad_prop * VectorValue<double>(1., 1., 1.);
  for (decltype(3) i = 0; i != 3; ++i)
  {
    nd_derivs_expect_near(vector_ad_prop.value()(i), 16., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[0](i), 32., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[1](i), 64., tol);
  }
  scalar_ad_prop = vector_ad_prop * vector_ad_prop;

  nd_derivs_expect_near(scalar_ad_prop.value(), 768., tol);
  nd_derivs_expect_near(scalar_ad_prop.derivatives()[0], 3072., tol);
  nd_derivs_expect_near(scalar_ad_prop.derivatives()[1], 6144., tol);

  vector_ad_prop *= 3.;
  for (decltype(3) i = 0; i != 3; ++i)
  {
    nd_derivs_expect_near(vector_ad_prop.value()(i), 48., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[0](i), 96., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[1](i), 192., tol);
  }
  vector_ad_prop = 4. * vector_ad_prop + vector_ad_prop * 5.;
  for (decltype(3) i = 0; i != 3; ++i)
  {
    nd_derivs_expect_near(vector_ad_prop.value()(i), 432., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[0](i), 864., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[1](i), 1728., tol);
  }

  tensor_ad_prop = scalar_ad_prop * TensorValue<double>(1., 1., 1., 1., 1., 1., 1., 1., 1.);
  for (decltype(3) i = 0; i != 3; ++i)
  {
    for (decltype(3) j = 0; j != 3; ++j)
    {
      nd_derivs_expect_near(tensor_ad_prop.value()(i, j), 768., tol);
      nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](i, j), 3072., tol);
      nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](i, j), 6144., tol);
    }
  }
  tensor_ad_prop *= 1. / 768.;
  for (decltype(3) i = 0; i != 3; ++i)
  {
    for (decltype(3) j = 0; j != 3; ++j)
    {
      nd_derivs_expect_near(tensor_ad_prop.value()(i, j), 1., tol);
      nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](i, j), 4., tol);
      nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](i, j), 8., tol);
      tensor_ad_prop.value()(i, j) = i + j;
      tensor_ad_prop.derivatives()[0](i, j) = 2. * i + 4. * j;
      tensor_ad_prop.derivatives()[1](i, j) = 3. * i + 5. * j;
    }
  }
  vector_ad_prop *= 1. / 432.;
  for (decltype(3) i = 0; i != 3; ++i)
  {
    nd_derivs_expect_near(vector_ad_prop.value()(i), 1., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[0](i), 2., tol);
    nd_derivs_expect_near(vector_ad_prop.derivatives()[1](i), 4., tol);
    vector_ad_prop.value()(i) = 10. - i;
    vector_ad_prop.derivatives()[0](i) = 10. - 2. * i;
    vector_ad_prop.derivatives()[1](i) = 10. - 3. * i;
  }

  vector_ad_prop = tensor_ad_prop * vector_ad_prop;
  nd_derivs_expect_near(vector_ad_prop.value()(0), 25, tol);
  nd_derivs_expect_near(vector_ad_prop.value()(1), 52, tol);
  nd_derivs_expect_near(vector_ad_prop.value()(2), 79, tol);
  nd_derivs_expect_near(vector_ad_prop.derivatives()[0](0), 120, tol);
  nd_derivs_expect_near(vector_ad_prop.derivatives()[0](1), 198, tol);
  nd_derivs_expect_near(vector_ad_prop.derivatives()[0](2), 276, tol);
  nd_derivs_expect_near(vector_ad_prop.derivatives()[1](0), 140, tol);
  nd_derivs_expect_near(vector_ad_prop.derivatives()[1](1), 242, tol);
  nd_derivs_expect_near(vector_ad_prop.derivatives()[1](2), 344, tol);

  tensor_ad_prop *= tensor_ad_prop;
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 0), 5, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 0), 30, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 0), 40, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 1), 8, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 1), 54, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 1), 70, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 2), 11, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 2), 78, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 2), 100, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 0), 8, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 0), 42, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 0), 58, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 1), 14, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 1), 84, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 1), 112, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 2), 20, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 2), 126, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 2), 166, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 0), 11, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 0), 54, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 0), 76, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 1), 20, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 1), 114, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 1), 154, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 2), 29, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 2), 174, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 2), 232, tol);

  tensor_ad_prop = tensor_ad_prop * tensor_ad_prop;
  tensor_ad_prop = 7. * tensor_ad_prop + tensor_ad_prop * 8.;
  tensor_ad_prop = 9. * tensor_ad_prop - tensor_ad_prop * 10.;

  nd_derivs_expect_near(tensor_ad_prop.value()(0, 0), -3150, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 0), -37800, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 0), -50400, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 1), -5580, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 1), -71280, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 1), -93600, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 2), -8010, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 2), -104760, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 2), -136800, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 0), -5580, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 0), -62640, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 0), -84960, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 1), -9900, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 1), -118800, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 1), -158400, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 2), -14220, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 2), -174960, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 2), -231840, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 0), -8010, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 0), -87480, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 0), -119520, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 1), -14220, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 1), -166320, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 1), -223200, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 2), -20430, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 2), -245160, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 2), -326880, tol);

  scalar_reg_prop = scalar_ad_prop.value();
  vector_reg_prop = vector_ad_prop.value();
  tensor_reg_prop = tensor_ad_prop.value();
  nd_derivs_expect_near(scalar_reg_prop, 768, tol);
  nd_derivs_expect_near(vector_reg_prop(0), 25, tol);
  nd_derivs_expect_near(vector_reg_prop(1), 52, tol);
  nd_derivs_expect_near(vector_reg_prop(2), 79, tol);
  nd_derivs_expect_near(tensor_reg_prop(0, 0), -3150, tol);
  nd_derivs_expect_near(tensor_reg_prop(0, 1), -5580, tol);
  nd_derivs_expect_near(tensor_reg_prop(0, 2), -8010, tol);
  nd_derivs_expect_near(tensor_reg_prop(1, 0), -5580, tol);
  nd_derivs_expect_near(tensor_reg_prop(1, 1), -9900, tol);
  nd_derivs_expect_near(tensor_reg_prop(1, 2), -14220, tol);
  nd_derivs_expect_near(tensor_reg_prop(2, 0), -8010, tol);
  nd_derivs_expect_near(tensor_reg_prop(2, 1), -14220, tol);
  nd_derivs_expect_near(tensor_reg_prop(2, 2), -20430, tol);

  TensorValue<DualNumber<double, NumberArray<2, double>>> inner_template;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
    {
      inner_template(i, j).value() = tensor_ad_prop.value()(i, j);
      for (unsigned di = 0; di < 2; ++di)
        inner_template(i, j).derivatives()[di] = tensor_ad_prop.derivatives()[di](i, j);
    }

  tensor_ad_prop = (tensor_ad_prop + tensor_ad_prop.transpose()) / 2.;
  auto strain = tensor_ad_prop.tr();
  tensor_ad_prop(0, 0) += (strain - tensor_ad_prop.tr() + tensor_ad_prop(2, 2)) / 3.;

  inner_template = (inner_template + inner_template.transpose()) / 2.;
  auto strain_inner = inner_template.tr();
  inner_template(0, 0) += (strain_inner - inner_template.tr() + inner_template(2, 2)) / 3.;

  nd_derivs_expect_near(tensor_ad_prop.value()(0, 0), -9960, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 1), -5580, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(0, 2), -8010, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 0), -5580, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 1), -9900, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(1, 2), -14220, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 0), -8010, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 1), -14220, tol);
  nd_derivs_expect_near(tensor_ad_prop.value()(2, 2), -20430, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 0), -119520, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 1), -66960, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](0, 2), -96120, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 0), -66960, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 1), -118800, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](1, 2), -170640, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 0), -96120, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 1), -170640, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[0](2, 2), -245160, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 0), -159360, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 1), -89280, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](0, 2), -128160, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 0), -89280, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 1), -158400, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](1, 2), -227520, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 0), -128160, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 1), -227520, tol);
  nd_derivs_expect_near(tensor_ad_prop.derivatives()[1](2, 2), -326880, tol);

  nd_derivs_expect_near(inner_template(0, 0).value(), -9960, tol);
  nd_derivs_expect_near(inner_template(0, 1).value(), -5580, tol);
  nd_derivs_expect_near(inner_template(0, 2).value(), -8010, tol);
  nd_derivs_expect_near(inner_template(1, 0).value(), -5580, tol);
  nd_derivs_expect_near(inner_template(1, 1).value(), -9900, tol);
  nd_derivs_expect_near(inner_template(1, 2).value(), -14220, tol);
  nd_derivs_expect_near(inner_template(2, 0).value(), -8010, tol);
  nd_derivs_expect_near(inner_template(2, 1).value(), -14220, tol);
  nd_derivs_expect_near(inner_template(2, 2).value(), -20430, tol);
  nd_derivs_expect_near(inner_template(0, 0).derivatives()[0], -119520, tol);
  nd_derivs_expect_near(inner_template(0, 1).derivatives()[0], -66960, tol);
  nd_derivs_expect_near(inner_template(0, 2).derivatives()[0], -96120, tol);
  nd_derivs_expect_near(inner_template(1, 0).derivatives()[0], -66960, tol);
  nd_derivs_expect_near(inner_template(1, 1).derivatives()[0], -118800, tol);
  nd_derivs_expect_near(inner_template(1, 2).derivatives()[0], -170640, tol);
  nd_derivs_expect_near(inner_template(2, 0).derivatives()[0], -96120, tol);
  nd_derivs_expect_near(inner_template(2, 1).derivatives()[0], -170640, tol);
  nd_derivs_expect_near(inner_template(2, 2).derivatives()[0], -245160, tol);
  nd_derivs_expect_near(inner_template(0, 0).derivatives()[1], -159360, tol);
  nd_derivs_expect_near(inner_template(0, 1).derivatives()[1], -89280, tol);
  nd_derivs_expect_near(inner_template(0, 2).derivatives()[1], -128160, tol);
  nd_derivs_expect_near(inner_template(1, 0).derivatives()[1], -89280, tol);
  nd_derivs_expect_near(inner_template(1, 1).derivatives()[1], -158400, tol);
  nd_derivs_expect_near(inner_template(1, 2).derivatives()[1], -227520, tol);
  nd_derivs_expect_near(inner_template(2, 0).derivatives()[1], -128160, tol);
  nd_derivs_expect_near(inner_template(2, 1).derivatives()[1], -227520, tol);
  nd_derivs_expect_near(inner_template(2, 2).derivatives()[1], -326880, tol);

  VectorValue<DualNumber<double, NumberArray<2, double>>> vector_inner_template;
  for (unsigned i = 0; i < 3; ++i)
  {
    vector_inner_template(i).value() = vector_ad_prop.value()(i);
    for (unsigned di = 0; di < 2; ++di)
      vector_inner_template(i).derivatives()[di] = vector_ad_prop.derivatives()[di](i);
  }

  auto outer_norm = vector_ad_prop.norm();
  auto inner_norm = vector_inner_template.norm();
  nd_derivs_expect_near(outer_norm.value(), inner_norm.value(), tol);
  nd_derivs_expect_near(outer_norm.derivatives()[0], inner_norm.derivatives()[0], tol);
  nd_derivs_expect_near(outer_norm.derivatives()[1], inner_norm.derivatives()[1], tol);

  vector_ad_prop.value() = VectorValue<double>(.05, 0, 0);
  vector_ad_prop.derivatives()[0] = VectorValue<double>(-0.5, 0, 0);
  vector_ad_prop.derivatives()[1] = VectorValue<double>(0.5, 0, 0);
  double dx[2] = {-0.5, 0.5};
  double dy[2] = {0, 0};
  double dz[2] = {0, 0};
  vector_inner_template(0) = DualNumber<double, NumberArray<2, double>>(.05, dx);
  vector_inner_template(1) = DualNumber<double, NumberArray<2, double>>(0, dy);
  vector_inner_template(2) = DualNumber<double, NumberArray<2, double>>(0, dz);
  outer_norm = vector_ad_prop.norm();
  inner_norm = vector_inner_template.norm();
  nd_derivs_expect_near(outer_norm.value(), inner_norm.value(), tol);
  nd_derivs_expect_near(outer_norm.derivatives()[0], inner_norm.derivatives()[0], tol);
  nd_derivs_expect_near(outer_norm.derivatives()[1], inner_norm.derivatives()[1], tol);

  DualNumberSurrogate<double, NumberArray<2, double*>> dns(0);
  double zero(0);
  dns.derivatives()[0] = &zero;
  dns.derivatives()[1] = &zero;
  DualNumberSurrogate<double, NumberArray<2, double*>> dns2(dns);
  DualNumberSurrogate<double, NumberArray<2, double*>> dns3(scalar_ad_prop);

  DualNumber<double, NumberArray<2, double>> new_scalar_ad_prop(dns);
  nd_derivs_expect_near(new_scalar_ad_prop.value(), 0, tol);
  nd_derivs_expect_near(new_scalar_ad_prop.derivatives()[0], 0, tol);
  nd_derivs_expect_near(new_scalar_ad_prop.derivatives()[1], 0, tol);

  new_scalar_ad_prop = dns3;
  nd_derivs_expect_near(new_scalar_ad_prop.value(), scalar_ad_prop.value(), tol);
  nd_derivs_expect_near(new_scalar_ad_prop.derivatives()[0], scalar_ad_prop.derivatives()[0], tol);
  nd_derivs_expect_near(new_scalar_ad_prop.derivatives()[1], scalar_ad_prop.derivatives()[1], tol);

  return returnval;
}
