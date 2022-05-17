#include "metaphysicl/dynamicsparsenumberarray.h"
#include "metaphysicl/metaphysicl_exceptions.h"

#define METAPHYSICL_UNIT_ASSERT(expr)                                                              \
  if (!(expr))                                                                                     \
  metaphysicl_error()

using namespace MetaPhysicL;

int main()
{
  MetaPhysicL::enableFPE(true);

  DynamicSparseNumberArray<double, unsigned int> test;
  test.insert(0) = 1e-14;
  test.sparsity_trim();
  METAPHYSICL_UNIT_ASSERT(test.size() == 1);
  test.sparsity_trim(1e-12);
  METAPHYSICL_UNIT_ASSERT(test.size() == 0);
}
