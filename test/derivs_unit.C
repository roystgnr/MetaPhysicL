#include <algorithm>
#include <fenv.h>
#include <iostream>

// List builtin classes, useful for disambiguation
template <typename T, typename Enable=void>
struct BuiltinTraits {
      static const bool value = false;
};

#define ScalarBuiltin_true(type) \
template<> \
struct BuiltinTraits<type> { static const bool value = true; }

ScalarBuiltin_true(int);
ScalarBuiltin_true(float);
ScalarBuiltin_true(double);

template<typename S, typename T, typename Enable=void>
struct CompareTypes {};

template<typename T, typename Enable>
struct CompareTypes<T, T, Enable> {
  typedef T supertype;
};


#define CompareTypes_super(a,b,super) \
	template<typename Enable> \
	struct CompareTypes<a, b, Enable> { \
	  typedef super supertype; \
	}

#define CompareTypes_all(mysub,mysuper) \
        CompareTypes_super(mysub, mysuper, mysuper); \
        CompareTypes_super(mysuper, mysub, mysuper)

CompareTypes_all(int, float);
CompareTypes_all(int, double);
CompareTypes_all(float, double);


template <typename T, typename D>
class DualNumber
{
public:
  T& value() { return _val; }

  const T& value() const { return _val; }

  D& derivatives() { return _deriv; }

  const D& derivatives() const { return _deriv; }

  DualNumber () = default;
  DualNumber (const DualNumber & a) = default;

  template <typename T2, typename D2>
  DualNumber (const DualNumber<T2,D2>& a) :
    _val(a.value()), _deriv(a.derivatives()) {}

  template <typename T2>
  DualNumber<T, D> & operator/= (const T2& a) {
    _val /= a;
    _deriv /= a;
    return *this;
  }

private:
  T _val;
  D _deriv;
};


template<typename T, typename T2, typename D>
struct CompareTypes<DualNumber<T, D>, T2,
                    typename std::enable_if<BuiltinTraits<T2>::value >::type> {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<
                       typename CompareTypes<D, T2>::supertype,
                       T
                     >::supertype> supertype;
};

template<typename T, typename D, typename T2, typename D2>
struct CompareTypes<DualNumber<T, D>, DualNumber<T2, D2> > {
  typedef DualNumber<typename CompareTypes<T,T2>::supertype,
                     typename CompareTypes<
                       typename CompareTypes<T, T2>::supertype,
                       typename CompareTypes<D, D2>::supertype
                     >::supertype> supertype;
};


#define roy_assert_less(expr1,expr2)  do { if (!(expr1 < expr2)) { std::cerr << "Assertion failed." << std::endl; } } while(0)

template <std::size_t N, typename T>
class NumberArray
{
public:

  NumberArray () = default;
  NumberArray (const NumberArray & a) = default;

  template <typename T2>
  NumberArray (const NumberArray<N,T2>& a)
    { for (std::size_t i=0; i != N; ++i) _data[i] = a[i]; }

  T& operator[](std::size_t i)
    { roy_assert_less(i, N); // Remove this assert and no more SIGFPE!
      return _data[i]; }

  const T& operator[](std::size_t i) const
    { roy_assert_less(i, N); // Remove this assert and no more SIGFPE!
      return _data[i]; }

  template <typename T2>
  NumberArray<N,T>& operator/= (const T2& a)
    { for (std::size_t i=0; i != N; ++i) _data[i] /= a; return *this; }

private:
  T _data[N];
};


template<std::size_t N, typename T>
struct CompareTypes<NumberArray<N,T>, NumberArray<N,T>> {
  typedef NumberArray<N, T> supertype; 
};

template<std::size_t N, typename T, typename T2>
struct CompareTypes<NumberArray<N,T>, NumberArray<N,T2>> {
  typedef NumberArray<N, typename CompareTypes<T, T2>::supertype> supertype;
};

template<std::size_t N, typename T, typename T2>
struct CompareTypes<T2, NumberArray<N, T>> {
  typedef NumberArray<N, typename CompareTypes<T, T2>::supertype> supertype;
};

template<std::size_t N, typename T, typename T2>
struct CompareTypes<NumberArray<N, T>, T2> {
  typedef NumberArray<N, typename CompareTypes<T, T2>::supertype> supertype;
};


template <typename T, typename D, typename T2>
inline
typename CompareTypes<DualNumber<T,D>,T2>::supertype
operator / (const DualNumber<T,D>& a, const T2& b)
{
  typename CompareTypes<DualNumber<T,D>,T2>::supertype
    returnval = a;
  returnval /= b;
  return returnval;
}


static const unsigned int N = 10; // test pts.

int main(int, char *[])
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID);

  typedef float Scalar;
  DualNumber<Scalar, NumberArray<N, Scalar>> dn;

  dn.value() = .25;
  for (unsigned int i=0; i != N; ++i)
    dn.derivatives()[i] = 1;

  auto working_divide = dn/3;
  // 3.0 doesn't FPE! 2 doesn't!
  DualNumber<double, NumberArray<N, double>> broken_manually{dn/3};

  (void)working_divide;
  (void)broken_manually;

  return 0;
}
