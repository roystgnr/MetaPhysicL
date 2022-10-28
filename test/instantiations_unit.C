
#if __cplusplus >= 201402L
#  include "metaphysicl/dualnamedarray.h"
#endif

#include "metaphysicl/dualnumberarray.h"
#include "metaphysicl/dualshadowsparsestruct.h"
#include "metaphysicl/dualshadowsparsevector.h"
#include "metaphysicl/dualshadowvector.h"
#include "metaphysicl/dynamicsparsenumberarray.h"
#include "metaphysicl/dynamicsparsenumbervector.h"
#include "metaphysicl/metaphysicl_exceptions.h"
#include "metaphysicl/semidynamicsparsenumberarray.h"
#include "metaphysicl/sparsenumberarray.h"

#include <ostream>

using namespace MetaPhysicL;

template <typename T1, typename T2>
struct Instantiator {

  Instantiator() {
    test_sdsna.insert(13)=6;
    test_dsna.insert(14)=7;
    test_dsnv.insert(15)=8;
  }

#if __cplusplus >= 201402L
  DualExpression<T1, T2> test_de;
#endif

  DualNumber<T1, T2> test_dn {1};
  ShadowNumber<T1, T2> test_shadow {2};
  NumberArray<5, T1> test_na {3};
  NumberVector<5, T1> test_nv {4};
  SemiDynamicSparseNumberArray<T1,int,NWrapper<5>> test_sdsna;
  DynamicSparseNumberArray<T1,int> test_dsna;
  DynamicSparseNumberVector<T1,int> test_dsnv;
  typename SparseNumberArrayOf<4, 2, T1, 3, T2, 5, T1, 7, T2>::type
          test_sna {0};
  typename SparseNumberVectorOf<4, 2, T1, 3, T2, 5, T1, 7, T2>::type
          test_snv {0};
  typename SparseNumberStructOf<4, 2, T1, 3, T2, 5, T1, 7, T2>::type
          test_sns {0};

#if __cplusplus >= 201402L
  NamedIndexArray
    <double,
     SparseNumberVector
       <long unsigned int,
        ULongSetConstructor<2>::type> > indexed_by_two {5, 0};
#endif
};

template <typename T1, typename T2>
void test_out(std::ostream & output, const Instantiator<T1,T2> & i)
{
#if __cplusplus >= 201402L
  output << "DualExpression: " << i.test_de << '\n';
  output << "NamedIndexArray: " << i.indexed_by_two << '\n';
#endif
  output << "DualNumber: " << i.test_dn << '\n';
  output << "ShadowNumber: " << i.test_shadow << '\n';
  output << "NumberArray: " << i.test_na << '\n';
  output << "NumberVector: " << i.test_nv << '\n';
  output << "SemiDynamicSparseNumberArray: " << i.test_sdsna << '\n';
  output << "DynamicSparseNumberArray: " << i.test_dsna << '\n';
  output << "DynamicSparseNumberVector: " << i.test_dsnv << '\n';
  output << "SparseNumberArray: " << i.test_sna << '\n';
  output << "SparseNumberVector: " << i.test_snv << '\n';
  output << "SparseNumberStruct: " << i.test_sns << '\n';
}

int main (void)
{
  MetaPhysicL::enableFPE(true);

  Instantiator<float, float> i1;
  Instantiator<float, double> i2;
  Instantiator<double, double> i3;
  Instantiator<int, int> i4;

  std::ostream output(std::cout.rdbuf());
  test_out(output, i1);
  test_out(output, i2);
  test_out(output, i3);
  test_out(output, i4);
  return 0;
}
