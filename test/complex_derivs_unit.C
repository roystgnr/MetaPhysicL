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

#define expect_nan(number)                                                                         \
  {                                                                                                \
    int new_returnval = !std::isnan(number);                                                       \
    if (new_returnval)                                                                             \
      std::cerr << "Failed test at line " << __LINE__ << std::endl;                                \
    returnval = returnval || new_returnval;                                                        \
  }

#define expect_nan_dualnumber(dualnumber)                                                          \
  expect_nan(dualnumber.derivatives()[0]);                                                         \
  expect_nan(dualnumber.derivatives()[1])

#define expect_complex_nan_dualnumber(dualnumber)                                                  \
  expect_nan(dualnumber.derivatives()[0].real());                                                  \
  expect_nan(dualnumber.derivatives()[1].real());                                                  \
  expect_nan(dualnumber.derivatives()[0].imag());                                                  \
  expect_nan(dualnumber.derivatives()[1].imag())

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
  expect_nan_dualnumber(dn_real);

  DualNumber<double, NumberArray<2, double>> dn_imag = std::imag(cdn);
  expect_near(dn_imag.value(), 1, tol);
  expect_nan_dualnumber(dn_imag);

  DualNumber<double, NumberArray<2, double>> dn_norm = std::norm(cdn);
  expect_near(dn_norm.value(), 2, tol);
  expect_nan_dualnumber(dn_norm);

  DualNumber<double, NumberArray<2, double>> dn_abs = std::abs(cdn);
  expect_near(dn_abs.value(), sqrt2, tol);
  expect_nan_dualnumber(dn_abs);

  DualNumber<std::complex<double>, NumberArray<2, std::complex<double>>> dn_conj = std::conj(cdn);
  expect_near(dn_conj.value().real(), 1, tol);
  expect_near(dn_conj.value().imag(), -1, tol);
  expect_complex_nan_dualnumber(dn_conj);

  DualNumber<std::complex<double>> cdn2{std::complex<double>{1., 1.}};

  DualNumber<double> dn_real2 = std::real(cdn2);
  expect_near(dn_real2.value(), 1, tol);
  expect_nan(dn_real2.derivatives());

  DualNumber<double> dn_imag2 = std::imag(cdn2);
  expect_near(dn_imag2.value(), 1, tol);
  expect_nan(dn_imag2.derivatives());

  DualNumber<double> dn_norm2 = std::norm(cdn2);
  expect_near(dn_norm2.value(), 2, tol);
  expect_nan(dn_norm2.derivatives());

  DualNumber<double> dn_abs2 = std::abs(cdn2);
  expect_near(dn_abs2.value(), sqrt2, tol);
  expect_nan(dn_abs2.derivatives());

  DualNumber<std::complex<double>> dn_conj2 = std::conj(cdn2);
  expect_near(dn_conj2.value().real(), 1, tol);
  expect_near(dn_conj2.value().imag(), -1, tol);
  expect_nan(dn_conj2.derivatives().real());
  expect_nan(dn_conj2.derivatives().imag());

  return returnval;
}
