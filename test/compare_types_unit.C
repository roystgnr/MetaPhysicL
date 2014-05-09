
#include "metaphysicl/compare_types.h"

using namespace MetaPhysicL;

template <typename T1, typename T2>
struct SubInstantiator {
  typename CompareTypes<T1 ,       T2 >::supertype st1;
  typename CompareTypes<T1 , const T2 >::supertype st3;
  typename CompareTypes<T1 , const T2&>::supertype st4;
};

template <typename T1, typename T2>
struct Instantiator {
  SubInstantiator<      T1 , T2> si1;
  SubInstantiator<const T1 , T2> si2;
  SubInstantiator<const T1&, T2> si4;
};


int main (void)
{
  Instantiator<float, float> i1;
  Instantiator<float, double> i2;
  Instantiator<double, double> i3;
  Instantiator<int, int> i4;
  return 0;
}
