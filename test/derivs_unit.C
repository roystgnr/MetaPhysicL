#include <cstdlib> // rand()
#include <iostream>
#include <limits>

#include "metaphysicl_config.h"

#include "metaphysicl/dualnumberarray.h"
#include "metaphysicl/dualnumbervector.h"

static const unsigned int N = 10; // test pts.

using namespace MetaPhysicL;

#define one_test(test_func, error_quant) \
  error_quant = raw_value(test_func); \
  { int new_returnval = test_error_quant(random_quant, error_quant); \
  if (new_returnval) \
    std::cerr << "Failed test: " << #test_func << std::endl; \
  returnval = returnval || new_returnval; }

template <typename DualScalar,
          typename Scalar,
          typename std::enable_if<ScalarTraits<Scalar>::value
                                  && ScalarTraits<DualScalar>::value, int>::type = 0>
int test_error_quant (const DualScalar& random_scalar,
                      const Scalar& error_scalar)
{
  using std::max;
  using std::fabs;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  if (fabs(error_scalar) > tol || error_scalar != error_scalar)
  {
    std::cerr << "Value " << random_scalar <<
      "\nError " << error_scalar <<
      "\nTol   " << tol << std::endl;
    return 1;
  }

  return 0;
}

template <typename DualVector,
          typename Vector,
          typename std::enable_if<!ScalarTraits<Vector>::value
                                  && !ScalarTraits<DualVector>::value, int>::type = 0>
int test_error_quant (const DualVector& random_vec,
                      const Vector& error_vec)
{
  using std::max;
  using std::fabs;

  typedef typename ValueType<Vector>::type Scalar;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  Scalar max_abs_error = 0;

  for (unsigned int i=0; i != error_vec.size(); ++i)
    {
      max_abs_error = max(max_abs_error, fabs(error_vec[i]));

      // Handle NaNs properly.  Testing max_abs_error for NaN is
      // impossible because IEEE sucks:
      // https://en.wikipedia.org/wiki/IEEE_754_revision#min_and_max
      if (max_abs_error > tol || error_vec[i] != error_vec[i])
        {
	  std::cerr << "Value " << random_vec[i] <<
		       "\nError " << error_vec[i] <<
		       "\nTol   " << tol << std::endl;
	  return 1;
        }
    }

  return 0;
}

template <typename DualScalar,
          typename Vector,
          typename std::enable_if<!ScalarTraits<Vector>::value
                                  && ScalarTraits<DualScalar>::value, int>::type = 0>
int test_error_quant (const DualScalar& random_scalar,
                      const Vector& error_vec)
{
  using std::max;
  using std::fabs;

  typedef typename ValueType<Vector>::type Scalar;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  Scalar max_abs_error = 0;

  for (unsigned int i=0; i != error_vec.size(); ++i)
    {
      max_abs_error = max(max_abs_error, fabs(error_vec[i]));

      // Handle NaNs properly.  Testing max_abs_error for NaN is
      // impossible because IEEE sucks:
      // https://en.wikipedia.org/wiki/IEEE_754_revision#min_and_max
      if (max_abs_error > tol || error_vec[i] != error_vec[i])
        {
	  std::cerr << "Value " << random_scalar.derivatives()[i] <<
		       "\nError " << error_vec[i] <<
		       "\nTol   " << tol << std::endl;
	  return 1;
        }
    }

  return 0;
}

template <typename T, typename T2>
int test_func_values(const T & random_quant, T2 & error_quant)
{
  using std::abs;
  using std::acos;
  using std::asin;
  using std::atan;
  using std::ceil;
  using std::cos;
  using std::cosh;
  using std::exp;
  using std::fabs;
  using std::floor;
  using std::log;
  using std::log10;
  using std::pow;
  using std::sin;
  using std::sinh;
  using std::sqrt;
  using std::tan;
  using std::tanh;

#if __cplusplus >= 201103L
  using std::exp2;
  using std::log2;
  using std::expm1;
  using std::log1p;
  using std::cbrt;
  using std::asinh;
  using std::acosh;
  using std::atanh;
  using std::erf;
  using std::erfc;
  using std::trunc;
  using std::round;
#endif // __cplusplus >= 201103L

  int returnval = 0;

  typedef typename ValueType<T2>::type Scalar;
  Scalar pi = acos(Scalar(-1));

  one_test(2*random_quant - random_quant - random_quant, error_quant);

  one_test(3*random_quant - random_quant*3, error_quant);

  one_test((random_quant + random_quant)/2 - random_quant, error_quant);

  one_test(sqrt(random_quant) * sqrt(random_quant) - random_quant, error_quant);
  one_test(random_quant*random_quant - pow(random_quant,2), error_quant);
  one_test(sqrt(random_quant) - pow(random_quant,Scalar(.5)), error_quant);

  one_test(log(exp(random_quant)) - random_quant, error_quant);
  one_test(exp(log(random_quant)) - random_quant, error_quant);
  one_test(exp(random_quant) - pow(exp(Scalar(1)), random_quant), error_quant);

  one_test(tan(random_quant) - sin(random_quant)/cos(random_quant), error_quant);
  one_test(random_quant - sin(asin(random_quant)), error_quant);
  one_test(random_quant - cos(acos(random_quant)), error_quant);
  one_test(random_quant - tan(atan(random_quant)), error_quant);
  one_test(1 - pow(sin(random_quant), 2) - pow(cos(random_quant), 2), error_quant);
  one_test(cos(random_quant) - sin(random_quant + pi/2), error_quant);

  one_test(tanh(random_quant) - sinh(random_quant)/cosh(random_quant), error_quant);
  one_test(1 + pow(sinh(random_quant), 2) - pow(cosh(random_quant), 2), error_quant);

  one_test(log10(random_quant) - log(random_quant)/log(Scalar(10)), error_quant);

  one_test(floor(random_quant / 2), error_quant);
  one_test(ceil(random_quant / 2 - 1), error_quant);

  one_test(abs(random_quant) - random_quant, error_quant);
  one_test(fabs(random_quant-.75) - abs(random_quant-.75), error_quant);

#if __cplusplus >= 201103L
  one_test(log2(exp2(random_quant)) - random_quant, error_quant);
  one_test(exp2(log2(random_quant)) - random_quant, error_quant);
  one_test(expm1(random_quant) - exp(random_quant) + 1, error_quant);
  one_test(log1p(random_quant) - log(random_quant + 1), error_quant);
  one_test(cbrt(random_quant) - pow(random_quant, Scalar(1)/3), error_quant);
  one_test(asinh(sinh(random_quant)) - random_quant, error_quant);
  one_test(acosh(cosh(random_quant)) - random_quant, error_quant);
  one_test(atanh(tanh(random_quant)) - random_quant, error_quant);
  one_test(1 - erf(random_quant) - erfc(random_quant), error_quant);
  one_test(trunc(3 * random_quant - 1.5), error_quant);
  one_test(round(2 * random_quant - 1), error_quant);
#endif // __cplusplus >= 201103L

  return returnval;
}

template <typename T, typename T2, typename T3, typename T4>
int test_func_derivatives(const T & random_quant,
                          const T & zero_quant,
                          T2 & error_quant,
                          const T3 & function_value,
                          const T4 & analytic_multiplier)
{
  using std::abs;
  using std::acos;
  using std::asin;
  using std::atan;
  using std::ceil;
  using std::cos;
  using std::cosh;
  using std::exp;
  using std::fabs;
  using std::floor;
  using std::log;
  using std::log10;
  using std::pow;
  using std::sin;
  using std::sinh;
  using std::sqrt;
  using std::tan;
  using std::tanh;

#if __cplusplus >= 201103L
  using std::exp2;
  using std::log2;
  using std::expm1;
  using std::log1p;
  using std::cbrt;
  using std::asinh;
  using std::acosh;
  using std::atanh;
  using std::trunc;
  using std::round;
  using std::nearbyint;
  using std::rint;
#endif // __cplusplus >= 201103L

  int returnval = 0;

  typedef typename ValueType<T2>::type Scalar;

  one_test(derivatives(pow(sin(random_quant-2),2)) -
	   2*sin(function_value-2)*cos(function_value-2)*analytic_multiplier, error_quant);

  one_test(derivatives(cos(2*random_quant)) + 2*sin(2*function_value)*analytic_multiplier, error_quant);
  one_test(derivatives(tan(.5*random_quant)) - .5/pow(cos(.5*function_value),2)*analytic_multiplier, error_quant);

  one_test(derivatives(sqrt(random_quant+1)) - 1/sqrt(function_value+1)/2*analytic_multiplier, error_quant);

  one_test(derivatives((random_quant-1)*(random_quant-1)) - 2*(function_value-1)*analytic_multiplier, error_quant);

  one_test(derivatives(pow(random_quant,1.5)) -
           1.5*pow(function_value,.5)*analytic_multiplier, error_quant);

  one_test(derivatives(exp(pow(random_quant,3))) -
           exp(pow(function_value,3))*3*pow(function_value,2)*analytic_multiplier, error_quant);

  one_test(derivatives(exp(random_quant)) -
           exp(function_value)*analytic_multiplier, error_quant);

  one_test(derivatives(pow(2,random_quant)) -
	   pow(2,function_value)*log(Scalar(2))*analytic_multiplier, error_quant);

  one_test(derivatives(asin(random_quant)) -
           1/sqrt(1-function_value*function_value)*analytic_multiplier, error_quant);

  one_test(derivatives(sinh(random_quant)) - cosh(function_value)*analytic_multiplier, error_quant);
  one_test(derivatives(cosh(random_quant)) - sinh(function_value)*analytic_multiplier, error_quant);

  one_test(derivatives(tanh(random_quant)) -
	   derivatives(sinh(random_quant)/cosh(random_quant)), error_quant);

#if __cplusplus >= 201103L
  one_test(derivatives(exp2(random_quant)) - exp2(random_quant) *
           analytic_multiplier * std::log(Scalar(2)), error_quant);
  one_test(derivatives(log2(random_quant)) - 1 / (random_quant * log(Scalar(2))), error_quant);
  one_test(derivatives(expm1(random_quant)) - (expm1(random_quant) + 1) * analytic_multiplier, error_quant);
  one_test(derivatives(log1p(random_quant)) - 1 / (1 + random_quant) * analytic_multiplier, error_quant);
  one_test(derivatives(cbrt(random_quant)) - pow(random_quant, -Scalar(2)/3) / 3 * analytic_multiplier, error_quant);
  one_test(pow(derivatives(asinh(random_quant)), 2) - 1 / (random_quant * random_quant + 1) * analytic_multiplier, error_quant);
  one_test(pow(derivatives(acosh(random_quant+1)), 2) * random_quant - 1 / (random_quant + 2) * analytic_multiplier, error_quant);
  one_test(derivatives(atanh(random_quant)) - 1 / (1 - random_quant * random_quant), error_quant);
  one_test(derivatives(trunc(random_quant*10)), error_quant);
  one_test(derivatives(round(random_quant*10)), error_quant);
  one_test(derivatives(nearbyint(random_quant*10)), error_quant);
  one_test(derivatives(rint(random_quant*10)), error_quant);

#endif // __cplusplus >= 201103L

  // Some non-random tests, too:
  one_test(derivatives(pow(zero_quant,2)), error_quant);
  one_test(derivatives(pow(zero_quant,1)) - 1, error_quant);

  return returnval;
}

template <typename Vector>
int scalartester (void)
{
  typedef typename Vector::value_type Scalar;
  typedef DualNumber<Scalar, Vector> DN;

  DN random_quant;

  Vector error_vec = 0;
  Vector unity_vec = 1;
  Scalar error_scalar = 0;

  std::srand(12345);

  random_quant.value() = .25 + (static_cast<Scalar>(std::rand())/static_cast<Scalar>(RAND_MAX)/2);
  DN zero_quant = 0;
  for (unsigned int i=0; i != N; ++i)
    {
      random_quant.derivatives()[i] = 1;
      zero_quant.derivatives() = 1;
    }


  int returnval = test_func_values(random_quant, error_scalar);
  returnval = returnval ||
    test_func_derivatives(random_quant, zero_quant, error_vec,
                          random_quant.value(), unity_vec);

  return returnval;
}

template <typename Vector>
int vectester (void)
{
  typedef typename ValueType<Vector>::type DualScalar;
  typedef typename DualScalar::value_type Scalar;

  Vector random_quant, zero_quant;

  typename DerivativeType<Vector>::type error_vec = 0;

  std::srand(12345); // Fixed seed for reproduceability of failures

  // Avoid divide by zero errors or acos(x>1) NaNs later
  for (unsigned int i=0; i != N; ++i)
    {
      random_quant[i] = .25 + (static_cast<Scalar>(std::rand())/static_cast<Scalar>(RAND_MAX)/2);
      random_quant[i].derivatives() = 1;
      zero_quant[i] = 0;
      zero_quant[i].derivatives() = 1;
    }

  int returnval =  test_func_values(random_quant, error_vec);
  returnval = returnval || test_func_derivatives(random_quant, zero_quant, error_vec, random_quant, 1.);

  return returnval;
}

int main(void)
{
  int returnval = 0;
  returnval = returnval || vectester<NumberArray<N, DualNumber<float> > >();
  returnval = returnval || vectester<NumberArray<N, DualNumber<double> > >();
  returnval = returnval || vectester<NumberArray<N, DualNumber<long double> > >();
  returnval = returnval || scalartester<NumberArray<N, float> >();
  returnval = returnval || scalartester<NumberArray<N, double> >();
  returnval = returnval || scalartester<NumberArray<N, long double> >();

  // We no longer treat vectors like arrays for built-in functions, so
  // most of the identities above make no sense.
  /*
  returnval = returnval || vectester<NumberVector<N, DualNumber<float> > >();
  returnval = returnval || vectester<NumberVector<N, DualNumber<double> > >();
  returnval = returnval || vectester<NumberVector<N, DualNumber<long double> > >();
  */

  return returnval;
}
