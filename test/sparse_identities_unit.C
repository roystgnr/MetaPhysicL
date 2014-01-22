#include <cstdlib> // rand()
#include <iostream>
#include <limits>

#include "metaphysicl_config.h"

#include "metaphysicl/numberarray.h"
#include "metaphysicl/numbervector.h"
#include "metaphysicl/sparsenumberarray.h"

static const unsigned int N = 10; // test pts.

using namespace MetaPhysicL;

#define one_test(test_func) \
  error_vec = test_func; \
  { int new_returnval = test_error_vec(random_vec, error_vec); \
  if (new_returnval) \
    std::cerr << "Failed test: " << #test_func << std::endl; \
  returnval = returnval || new_returnval; }

template <typename Vector>
int test_error_vec (const Vector& random_vec,
                    const Vector& error_vec)
{
  using std::max;
  using std::fabs;

  typedef typename Vector::value_type Scalar;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  Scalar max_abs_error = 0;

  for (unsigned int i=0; i != error_vec.size(); ++i)
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

template <typename Vector>
int vectester (void)
{
  using std::abs;
  using std::acos;
  using std::asin;
  using std::atan;
  using std::ceil;
  using std::fabs;
  using std::floor;
  using std::pow;
  using std::sin;
  using std::sqrt;
  using std::tan;

  typedef typename Vector::value_type Scalar;

  Vector random_vec;

  Vector error_vec = 0;

  std::srand(12345); // Fixed seed for reproduceability of failures

  for (unsigned int i=0; i != N; ++i)
    random_vec[i] = .25 + (static_cast<Scalar>(std::rand())/RAND_MAX);

  int returnval = 0;

  one_test(2*random_vec - random_vec - random_vec);

  one_test(3*random_vec - random_vec*3);

  one_test((random_vec + random_vec)/2 - random_vec);

  one_test(sqrt(random_vec) * sqrt(random_vec) - random_vec);
  one_test(random_vec*random_vec - pow(random_vec,2));
  one_test(sqrt(random_vec) - pow(random_vec,Scalar(.5)));

  one_test(random_vec - sin(asin(random_vec)));
  one_test(random_vec - tan(atan(random_vec)));

  one_test(floor(random_vec / 2));
  one_test(ceil(random_vec / 2 - 1));

  one_test(abs(random_vec) - random_vec);
  one_test(fabs(random_vec-.75) - abs(random_vec-.75));


  return returnval;
}

int main(void)
{
  int returnval = 0;
  returnval = returnval || vectester<NumberArray<N, float> >();
  returnval = returnval || vectester<NumberArray<N, double> >();
  returnval = returnval || vectester<NumberArray<N, long double> >();

  returnval = returnval || vectester<NumberVector<N, float> >();
  returnval = returnval || vectester<NumberVector<N, double> >();
  returnval = returnval || vectester<NumberVector<N, long double> >();

/*
  returnval = returnval || vectester<SparseNumberArrayOf<4,
                                                         0, float,
                                                         1, float,
                                                         2, float,
                                                         3, float>::type >();
  returnval = returnval || vectester<SparseNumberArrayOf<4,
                                                         0, double,
                                                         1, double,
                                                         2, double,
                                                         3, double>::type >();
  returnval = returnval || vectester<SparseNumberArrayOf<4,
                                                         0, long double,
                                                         1, long double,
                                                         2, long double,
                                                         3, long double>::type >();
*/

  return returnval;
}
