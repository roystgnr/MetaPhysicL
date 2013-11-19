
#include <iostream>

#include "metaphysicl/sparsenumbervector.h"
#include "metaphysicl/namedindexarray.h"

using namespace MetaPhysicL;

int main(void)
{
  typedef
    NamedIndexArray
      <double,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<3>::type> >
    indexed_by_three;

  indexed_by_three test_val;
  test_val.raw_data() = 1;
  test_val.raw_sizes().get<3>() = 1;

  return 0;
}
