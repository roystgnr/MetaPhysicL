#include <cstdlib> // rand()
#include <iostream>
#include <limits>

#include "metaphysicl_config.h"

#include "metaphysicl/numberarray.h"

static const unsigned int N = 10; // test pts.

using namespace MetaPhysicL;

#define one_test(test_func) \
  error_vec = test_func; \
  { int new_returnval = test_error_vec(random_vec, error_vec); \
  if (new_returnval) \
    std::cerr << "Failed test: " << #test_func << std::endl; \
  returnval = returnval || new_returnval; }

template <typename Scalar>
int test_error_vec (const NumberArray<N,Scalar>& random_vec,
                    const NumberArray<N,Scalar>& error_vec)
{
  using std::max;
  using std::fabs;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  Scalar max_abs_error = 0;
  for (unsigned int i=0; i != N; ++i)
    {
      max_abs_error = max(max_abs_error, fabs(error_vec[i]));

      if (max_abs_error > tol)
        {
	  std::cerr << "Value " << random_vec[i] <<
		       "\nError " << error_vec[i] <<
		       "\nTol   " << tol << std::endl;
	  return 1;
        }
    }

  return 0;
}

template <typename Scalar>
int vectester (void)
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

  NumberArray<N, Scalar> random_vec;

  NumberArray<N, Scalar> error_vec = 0;

  std::srand(12345); // Fixed seed for reproduceability of failures

  // Avoid divide by zero errors later
  for (unsigned int i=0; i != N; ++i)
    random_vec[i] = .25 + (static_cast<Scalar>(std::rand())/RAND_MAX);

  Scalar pi = acos(Scalar(-1));

  int returnval = 0;

  one_test(2*random_vec - random_vec - random_vec);

  one_test(3*random_vec - random_vec*3);

  one_test((random_vec + random_vec)/2 - random_vec);

  one_test(sqrt(random_vec) * sqrt(random_vec) - random_vec);
  one_test(random_vec*random_vec - pow(random_vec,2));
  one_test(sqrt(random_vec) - pow(random_vec,Scalar(.5)));

  one_test(log(exp(random_vec)) - random_vec);
  one_test(exp(log(random_vec)) - random_vec);
  one_test(exp(random_vec) - pow(exp(Scalar(1)), random_vec));

  one_test(tan(random_vec) - sin(random_vec)/cos(random_vec));
  one_test(random_vec - sin(asin(random_vec)));
  one_test(random_vec - cos(acos(random_vec)));
  one_test(random_vec - tan(atan(random_vec)));
  one_test(1 - pow(sin(random_vec), 2) - pow(cos(random_vec), 2));
  one_test(cos(random_vec) - sin(random_vec + pi/2));

  one_test(tanh(random_vec) - sinh(random_vec)/cosh(random_vec));
  one_test(1 + pow(sinh(random_vec), 2) - pow(cosh(random_vec), 2));

  one_test(log10(random_vec) - log(random_vec)/log(Scalar(10)));

  one_test(floor(random_vec / 2));
  one_test(ceil(random_vec / 2 - 1));

  one_test(abs(random_vec) - random_vec);
  one_test(fabs(random_vec-.75) - abs(random_vec-.75));

  return returnval;
}

int main(void)
{
  int returnval = 0;
  returnval = returnval || vectester<float>();
  returnval = returnval || vectester<double>();
  returnval = returnval || vectester<long double>();

  return returnval;
}
