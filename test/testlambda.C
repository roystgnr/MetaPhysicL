#include <iostream>

const auto mylambda1 = [](const float& s, const float& t){return s + t;};

struct WackyFunctor
{
/*
  template <typename S, typename T>
  void MyLambda(auto mylambda = [](const S& s, const T& t){return s + t;})
  {
  };
*/

/*
  template <typename S, typename T>
  auto operator()(const S& sin, const T& tin)
  -> decltype(MyLambda<S,T>::mylambda(sin, tin)) {
    return MyLambda<S,T>::mylambda(sin, tin);
  }
*/
};

int main(void)
{
//  std::cout << WackyFunctor()(1, 2) << std::endl;

  return 0.;
}
