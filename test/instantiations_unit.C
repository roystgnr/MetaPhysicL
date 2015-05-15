
#include "metaphysicl/dualnamedarray.h"
#include "metaphysicl/dualnumberarray.h"
#include "metaphysicl/dualshadowsparsestruct.h"
#include "metaphysicl/dualshadowsparsevector.h"
#include "metaphysicl/dualshadowvector.h"
#include "metaphysicl/namedindexarray.h"
#include "metaphysicl/sparsenumberarray.h"

using namespace MetaPhysicL;

template <typename T1, typename T2>
struct Instantiator {
  DualExpression<T1, T2> test_de;
  DualNumber<T1, T2> test_dn;
  ShadowNumber<T1, T2> test_shadow;
  NumberArray<5, T1> test_na;
  NumberVector<5, T1> test_nv;
  typename SparseNumberArrayOf<4, 2, T1, 3, T2, 5, T1, 7, T2>::type
          test_sna;
  typename SparseNumberVectorOf<4, 2, T1, 3, T2, 5, T1, 7, T2>::type
          test_snv;
  typename SparseNumberStructOf<4, 2, T1, 3, T2, 5, T1, 7, T2>::type
          test_sns;
};

int main (void)
{
  Instantiator<float, float> i1;
  Instantiator<float, double> i2;
  Instantiator<double, double> i3;
  Instantiator<int, int> i4;
  return 0;
}
