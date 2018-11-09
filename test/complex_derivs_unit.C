#include <complex>
#include <iostream>

#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"

using namespace MetaPhysicL;

#define expect_near(double1, double2, tolerance)                                                   \
  {                                                                                                \
    int new_returnval = fabs(double1 - double2) > tolerance;                                       \
    if (new_returnval)                                                                             \
      std::cerr << "Failed test at line " << __LINE__ << std::endl;                                \
    returnval = returnval || new_returnval;                                                        \
  }

#define expect_nan(dualnumber)                                                                     \
  {                                                                                                \
    int new_returnval =                                                                            \
        !(std::isnan(dualnumber.derivatives()[0]) && std::isnan(dualnumber.derivatives()[1]));     \
    if (new_returnval)                                                                             \
      std::cerr << "Failed test at line " << __LINE__ << std::endl;                                \
    returnval = returnval || new_returnval;                                                        \
  }

#define expect_complex_nan(dualnumber)                                                             \
  {                                                                                                \
    int new_returnval = !(std::isnan(dualnumber.derivatives()[0].real()) &&                        \
                          std::isnan(dualnumber.derivatives()[1].real()) &&                        \
                          std::isnan(dualnumber.derivatives()[0].imag()) &&                        \
                          std::isnan(dualnumber.derivatives()[1].imag()));                         \
    if (new_returnval)                                                                             \
      std::cerr << "Failed test at line " << __LINE__ << std::endl;                                \
    returnval = returnval || new_returnval;                                                        \
  }

int
main()
{
  auto sqrt2 = std::sqrt(2.0);
  double tol = 1e-8;
  int returnval = 0;

  DualNumber<std::complex<double>, NumberArray<2, std::complex<double>>> cdn{
      std::complex<double>{1., 1.}};

  DualNumber<double, NumberArray<2, double>> dn_real = std::real(cdn);
  expect_near(dn_real.value(), 1, tol);
  expect_nan(dn_real);

  DualNumber<double, NumberArray<2, double>> dn_imag = std::imag(cdn);
  expect_near(dn_imag.value(), 1, tol);
  expect_nan(dn_imag);

  DualNumber<double, NumberArray<2, double>> dn_norm = std::norm(cdn);
  expect_near(dn_norm.value(), 2, tol);
  expect_nan(dn_norm);

  DualNumber<double, NumberArray<2, double>> dn_abs = std::abs(cdn);
  expect_near(dn_abs.value(), sqrt2, tol);
  expect_nan(dn_abs);

  DualNumber<std::complex<double>, NumberArray<2, std::complex<double>>> dn_conj = std::conj(cdn);
  expect_near(dn_conj.value().real(), 1, tol);
  expect_near(dn_conj.value().imag(), -1, tol);
  expect_complex_nan(dn_conj);

  return returnval;
}
