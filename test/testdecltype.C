#include <complex>

template <typename T>
struct testtypes
{
  typedef float type1;
  typedef double type2;
};

auto testfunc(float f, double d)
  -> decltype(std::complex<double>(f) + d)
{
  auto c = std::complex<double>(f);
  auto r = c + d;
  return r;
}

template <typename T>
struct testdecl
{

  typedef decltype(::testfunc(testtypes<int>::type1(), testtypes<int>::type2())) type3;
  typedef decltype(::testfunc(testtypes<T>::type1(), testtypes<T>::type2())) type4;
};

int main(void)
{
  testfunc(1, 2);

  decltype(testfunc(testtypes<int>::type1(), testtypes<int>::type2())) a;

  return 0.;
}
