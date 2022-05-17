#include <iostream>
#include <exception>

#include "metaphysicl_config.h"

#include "metaphysicl/dynamic_std_array_wrapper.h"
#include "metaphysicl/metaphysicl_exceptions.h"

using namespace MetaPhysicL;

#define EXPECT_EQ(v1, v2, name)                                                                    \
  {                                                                                                \
    if (v1 != v2) {                                                                                \
      std::cerr << "Failed test '" << name << "' at line " << __LINE__ << std::endl;               \
      returnval++;                                                                                 \
    }                                                                                              \
  }

int main(int, char * [])
{
  MetaPhysicL::enableFPE(true);

  int returnval = 0;

  DynamicStdArrayWrapper<int, NWrapper<7>> test;
  test.resize(5);
  for (unsigned int i = 0; i < 5; ++i)
    test[i] = i + 1;
  EXPECT_EQ(test.size(), 5, "constructed size");

  // test forward iterator
  {
    unsigned int r = 0;
    unsigned int n = 0;
    for (auto it = test.begin(); it != test.end(); ++it)
      r += ++n * *it;
    EXPECT_EQ(r, 1*1 + 2*2 + 3*3 + 4*4 + 5*5, "forward iterator");
  }

  // test reverse iterator
  {
    unsigned int r = 0;
    unsigned int n = 0;
    for (auto it = test.rbegin(); it != test.rend(); ++it)
      r += ++n * *it;
    EXPECT_EQ(r, 5*1 + 4*2 + 3*3 + 2*4 + 1*5, "reverse iterator");
  }

  // test resize
  test.resize(3);
  EXPECT_EQ(test.size(), 3, "resize");
  {
    unsigned int r = 0;
    for (auto & i : test)
      r += i;
    EXPECT_EQ(r, 6, "resized data");
  }

  // // test const_iterator
  // const auto & ctest = test;
  // count = 0;
  // for (auto & i : ctest)
  //   count += i;
  // EXPECT_EQ(count, 1 + 2 + 3 + 4 + 5 + 6);

  return returnval;
}
