#ifndef METAPHYSICL_MATH_STRUCTS_H
#define METAPHYSICL_MATH_STRUCTS_H

#include <tuple>
#include <map>
#include <utility>
#include <cmath>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/nddualnumber.h"

namespace MetaPhysicL
{

template <typename T>
class VectorValue
{
public:
  VectorValue() : _coords() {}

  VectorValue(const T x, const T y = 0, const T z = 0)
  {
    _coords[0] = x;
    _coords[1] = y;
    _coords[2] = z;
  }

  typedef T value_type;
  typedef std::tuple<unsigned int> index_type;

private:
  T _coords[3];

public:
  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  auto operator*(const T2 & scalar) const -> VectorValue<decltype(this->_coords[0] * scalar)>
  {
    return {_coords[0] * scalar, _coords[1] * scalar, _coords[2] * scalar};
  }

  template <typename T2>
  auto operator*(const VectorValue<T2> & vector) const -> decltype(this->_coords[0] * vector(0))
  {
    decltype(_coords[0] * vector(0)) ret = 0;
    for (unsigned i = 0; i < 3; ++i)
      ret += _coords[i] * vector(i);
    return ret;
  }

  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  VectorValue<T> & operator*=(const T2 & scalar)
  {
    for (unsigned i = 0; i < 3; ++i)
      _coords[i] *= scalar;
    return *this;
  }

  template <typename T2>
  auto operator+(const VectorValue<T2> & vector) const
      -> VectorValue<decltype(this->_coords[0] + vector(0))>
  {
    return {_coords[0] + vector(0), _coords[1] + vector(1), _coords[2] + vector(2)};
  }

  template <typename T2>
  auto operator-(const VectorValue<T2> & vector) const
      -> VectorValue<decltype(this->_coords[0] - vector(0))>
  {
    return {_coords[0] - vector(0), _coords[1] - vector(1), _coords[2] - vector(2)};
  }

  VectorValue<T> operator-() const { return {-_coords[0], -_coords[1], -_coords[2]}; }

  const T & operator()(unsigned i) const { return _coords[i]; }
  T & operator()(unsigned i) { return _coords[i]; }

  T norm() const
    {
      return std::sqrt(_coords[0] * _coords[0] + _coords[1] * _coords[1] + _coords[2] * _coords[2]);
    }
};

template <typename T, typename T2, typename std::enable_if<ScalarTraits<T>::value, int>::type = 0>
auto operator*(const T & scalar, const VectorValue<T2> & vector)
    -> VectorValue<decltype(scalar * vector(0))>
{
  return {scalar * vector(0), scalar * vector(1), scalar * vector(2)};
}

template <bool reverseorder>
struct PlusType<VectorValue<double>, VectorValue<double>, reverseorder>
{
  typedef VectorValue<double> supertype;
};

template <bool reverseorder>
struct MinusType<VectorValue<double>, VectorValue<double>, reverseorder>
{
  typedef VectorValue<double> supertype;
};

template <bool reverseorder>
struct MultipliesType<VectorValue<double>, VectorValue<double>, reverseorder>
{
  typedef double supertype;
};

template <bool reverseorder>
struct MultipliesType<double, VectorValue<double>, reverseorder>
{
  typedef VectorValue<double> supertype;
};

template <>
struct DividesType<VectorValue<double>, double>
{
  typedef VectorValue<double> supertype;
};

template <typename>
class TensorValue;

template <typename T>
struct FunctionReturnTensorValue
{
  typedef TensorValue<T> type;
};

template <typename T, template <size_t, typename> class DerivativeContainer, size_t N>
struct FunctionReturnTensorValue<DualNumber<T, DerivativeContainer<N, T>>>
{
  typedef NotADuckDualNumber<TensorValue<T>, DerivativeContainer<N, TensorValue<T>>> type;
};

template <typename T>
class TensorValue
{
public:
  TensorValue() : _coords() {}

  TensorValue(const T xx,
              const T xy = 0,
              const T xz = 0,
              const T yx = 0,
              const T yy = 0,
              const T yz = 0,
              const T zx = 0,
              const T zy = 0,
              const T zz = 0)
  {
    _coords[0] = xx;
    _coords[1] = xy;
    _coords[2] = xz;
    _coords[3] = yx;
    _coords[4] = yy;
    _coords[5] = yz;
    _coords[6] = zx;
    _coords[7] = zy;
    _coords[8] = zz;
  }

  TensorValue(const VectorValue<T> & row1, const VectorValue<T> & row2, const VectorValue<T> & row3)
  {
    _coords[0] = row1(0);
    _coords[1] = row1(1);
    _coords[2] = row1(2);
    _coords[3] = row2(0);
    _coords[4] = row2(1);
    _coords[5] = row2(2);
    _coords[6] = row3(0);
    _coords[7] = row3(1);
    _coords[8] = row3(2);
  }

  typedef T value_type;
  typedef std::tuple<unsigned int, unsigned int> index_type;

private:
  T _coords[9];

public:
  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  auto operator*(const T2 & scalar) const -> TensorValue<decltype(scalar * this->_coords[0])>
  {
    TensorValue<decltype(scalar * _coords[0])> ret;
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        ret(i, j) = (*this)(i, j) * scalar;
    return ret;
  }

  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  auto operator/(const T2 & scalar) const -> TensorValue<decltype(this->_coords[0] / scalar)>
  {
    TensorValue<decltype(scalar * _coords[0])> ret;
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        ret(i, j) = (*this)(i, j) / scalar;
    return ret;
  }

  template <typename T2>
  auto operator*(const VectorValue<T2> & vector) const
      -> VectorValue<decltype(this->_coords[0] * vector(0))>
  {
    VectorValue<decltype(_coords[0] * vector(0))> ret;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        ret(i) += (*this)(i, j) * vector(j);

    return ret;
  }

  template <typename T2>
  auto operator*(const TensorValue<T2> & tensor) const
      -> TensorValue<decltype(this->_coords[0] * tensor(0, 0))>
  {
    TensorValue<decltype(_coords[0] * tensor(0, 0))> ret;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        for (unsigned int k = 0; k < 3; ++k)
          ret(i, j) += (*this)(i, k) * tensor(k, j);

    return ret;
  }

  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  TensorValue<T> & operator*=(const T2 & scalar)
  {
    for (unsigned i = 0; i < 9; ++i)
      _coords[i] *= scalar;
    return *this;
  }

  template <typename T2>
  TensorValue<T> & operator*=(const TensorValue<T2> & tensor)
  {
    TensorValue<T> ret;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        for (unsigned int k = 0; k < 3; ++k)
          ret(i, j) += (*this)(i, k) * tensor(k, j);

    *this = ret;
    return *this;
  }

  template <typename T2>
  auto operator+(const TensorValue<T2> & tensor) const
      -> TensorValue<decltype(this->_coords[0] + tensor(0, 0))>
  {
    TensorValue<decltype(_coords[0] + tensor(0, 0))> ret;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        ret(i, j) = (*this)(i, j) + tensor(i, j);
    return ret;
  }

  template <typename T2>
  TensorValue<T> & operator+=(const TensorValue<T2> & tensor)
  {
    for (unsigned int i = 0; i < 9; ++i)
      _coords[i] += tensor(i);
    return *this;
  }

  template <typename T2>
  auto operator-(const TensorValue<T2> & tensor) const
      -> TensorValue<decltype(this->_coords[0] - tensor(0, 0))>
  {
    TensorValue<decltype(_coords[0] - tensor(0, 0))> ret;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        ret(i, j) = (*this)(i, j) - tensor(i, j);
    return ret;
  }

  TensorValue<T> operator-() const
  {
    TensorValue<T> ret;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        ret(i, j) = -(*this)(i, j);
    return ret;
  }

  const T & operator()(unsigned i, unsigned j) const { return _coords[i * 3 + j]; }

  T & operator()(unsigned i, unsigned j) { return _coords[i * 3 + j]; }

  const T & operator()(unsigned i) const { return _coords[i]; }
  T & operator()(unsigned i) { return _coords[i]; }

  TensorValue<T> transpose() const
  {
    return {_coords[0],
            _coords[3],
            _coords[6],
            _coords[1],
            _coords[4],
            _coords[7],
            _coords[2],
            _coords[5],
            _coords[8]};
  }

  T tr() const { return _coords[0] + _coords[4] + _coords[8]; }
};

template <typename T, typename T2, typename std::enable_if<ScalarTraits<T>::value, int>::type = 0>
auto operator*(const T & scalar, const TensorValue<T2> & tensor)
    -> TensorValue<decltype(scalar * tensor(0, 0))>
{
  TensorValue<decltype(scalar * tensor(0, 0))> ret;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      ret(i, j) = scalar * tensor(i, j);

  return ret;
}

template <bool reverseorder>
struct PlusType<TensorValue<double>, TensorValue<double>, reverseorder>
{
  typedef TensorValue<double> supertype;
};

template <bool reverseorder>
struct MinusType<TensorValue<double>, TensorValue<double>, reverseorder>
{
  typedef TensorValue<double> supertype;
};

template <bool reverseorder>
struct MultipliesType<TensorValue<double>, TensorValue<double>, reverseorder>
{
  typedef TensorValue<double> supertype;
};

template <bool reverseorder>
struct MultipliesType<double, TensorValue<double>, reverseorder>
{
  typedef TensorValue<double> supertype;
};

template <>
struct MultipliesType<TensorValue<double>, VectorValue<double>>
{
  typedef VectorValue<double> supertype;
};

template <>
struct DividesType<TensorValue<double>, double>
{
  typedef TensorValue<double> supertype;
};

template <typename D>
class NotADuckDualNumber<VectorValue<double>, D> : public DualNumber<VectorValue<double>, D>
{
public:
  NotADuckDualNumber() : DualNumber<VectorValue<double>, D>() {}

  using DualNumber<VectorValue<double>, D>::DualNumber;

  typedef typename D::template resize<3>::other ResizedDerivatives;
  typedef typename D::template rebind<double>::other::template resize<3>::other DuckNumberDerivatives;

  template <typename KeyType, typename ValueType>
  using Map = std::map<KeyType, ValueType>;

  NotADuckDualNumber<VectorValue<double>, D> operator-() const
  {
    return NotADuckDualNumber<VectorValue<double>, D>(-this->_val, -this->_deriv);
  }
  NotADuckDualNumber<VectorValue<double>, D> operator!() const
  {
    return NotADuckDualNumber<VectorValue<double>, D>(!this->_val, !this->_deriv);
  }

  VectorValue<DualNumber<double, DuckNumberDerivatives>> inner_template() const
  {
    DualNumber<double, DuckNumberDerivatives> x{this->value()(0), _e0};
    DualNumber<double, DuckNumberDerivatives> y{this->value()(1), _e1};
    DualNumber<double, DuckNumberDerivatives> z{this->value()(2), _e2};
    return {x, y, z};
  }

  NotADuckDualNumber<VectorValue<double>, ResizedDerivatives>
  convert_to_outer(const VectorValue<DualNumber<double, DuckNumberDerivatives>> & inner) const
  {
    NotADuckDualNumber<VectorValue<double>, ResizedDerivatives> nd = VectorValue<double>(inner(0).value(), inner(1).value(), inner(2).value());

    for (unsigned di = 0; di < 3; ++di)
      for (unsigned i = 0; i < 3; ++i)
        nd.derivatives()[di](i) = inner(i).derivatives()[di];
    return nd;
  }

  const DualNumber<double, DuckNumberDerivatives> &
  convert_to_outer(const DualNumber<double, DuckNumberDerivatives> & inner) const
    {
      return inner;
    }

  static constexpr double _e0[3] = {1, 0, 0};
  static constexpr double _e1[3] = {0, 1, 0};
  static constexpr double _e2[3] = {0, 0, 1};

  metaphysicl_const_return_decl(operator());
  metaphysicl_nonconst_return_decl(operator());
  metaphysicl_const_return_decl(norm);

  metaphysicl_map_api_decl(VectorValue<double>);

private:
  metaphysicl_map_decl(VectorValue<double>);
};

metaphysicl_const_return_def(operator(), VectorValue<double>)
metaphysicl_nonconst_return_def(operator(), VectorValue<double>)
metaphysicl_const_return_def(norm, VectorValue<double>)

metaphysicl_map_api_def(VectorValue<double>)

template <typename D>
constexpr double NotADuckDualNumber<VectorValue<double>, D>::_e0[3];
template <typename D>
constexpr double NotADuckDualNumber<VectorValue<double>, D>::_e1[3];
template <typename D>
constexpr double NotADuckDualNumber<VectorValue<double>, D>::_e2[3];


template <typename D>
class NotADuckDualNumber<TensorValue<double>, D> : public DualNumber<TensorValue<double>, D>
{
public:
  NotADuckDualNumber() : DualNumber<TensorValue<double>, D>() {}

  template <typename D2>
  NotADuckDualNumber(const NotADuckDualNumber<VectorValue<double>, D2> & row1,
                     const NotADuckDualNumber<VectorValue<double>, D2> & row2,
                     const NotADuckDualNumber<VectorValue<double>, D2> & row3)
    : DualNumber<TensorValue<double>, D>()
  {
    this->value() = TensorValue<double>(row1.value(), row2.value(), row3.value());
    auto size = this->derivatives().size();
    for (decltype(size) i = 0; i < size; ++i)
      this->derivatives()[i] =
          TensorValue<double>(row1.derivatives()[i], row2.derivatives()[i], row3.derivatives()[i]);
  }

  using DualNumber<TensorValue<double>, D>::DualNumber;

  typedef typename D::template resize<9>::other ResizedDerivatives;
  typedef typename D::template rebind<double>::other::template resize<9>::other DuckNumberDerivatives;

  template <typename KeyType, typename ValueType>
  using Map = std::map<KeyType, ValueType>;

  NotADuckDualNumber<TensorValue<double>, D> operator-() const
  {
    return NotADuckDualNumber<TensorValue<double>, D>(-this->_val, -this->_deriv);
  }
  NotADuckDualNumber<TensorValue<double>, D> operator!() const
  {
    return NotADuckDualNumber<TensorValue<double>, D>(!this->_val, !this->_deriv);
  }

  TensorValue<DualNumber<double, DuckNumberDerivatives>> inner_template() const
  {
    DualNumber<double, DuckNumberDerivatives> xx{this->value()(0, 0), _e0};
    DualNumber<double, DuckNumberDerivatives> xy{this->value()(0, 1), _e1};
    DualNumber<double, DuckNumberDerivatives> xz{this->value()(0, 2), _e2};
    DualNumber<double, DuckNumberDerivatives> yx{this->value()(1, 0), _e3};
    DualNumber<double, DuckNumberDerivatives> yy{this->value()(1, 1), _e4};
    DualNumber<double, DuckNumberDerivatives> yz{this->value()(2, 1), _e5};
    DualNumber<double, DuckNumberDerivatives> zx{this->value()(2, 0), _e6};
    DualNumber<double, DuckNumberDerivatives> zy{this->value()(2, 1), _e7};
    DualNumber<double, DuckNumberDerivatives> zz{this->value()(2, 2), _e8};
    return {xx, xy, xz, yx, yy, yz, zx, zy, zz};
  }

  NotADuckDualNumber<TensorValue<double>, ResizedDerivatives>
  convert_to_outer(const TensorValue<DualNumber<double, DuckNumberDerivatives>> & inner) const
  {
    NotADuckDualNumber<TensorValue<double>, ResizedDerivatives> nd = TensorValue<double>(inner(0, 0).value(),
                                                                               inner(0, 1).value(),
                                                                               inner(0, 2).value(),
                                                                               inner(1, 0).value(),
                                                                               inner(1, 1).value(),
                                                                               inner(1, 2).value(),
                                                                               inner(2, 0).value(),
                                                                               inner(2, 1).value(),
                                                                               inner(2, 2).value());
    for (unsigned di = 0; di < 9; ++di)
      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
          nd.derivatives()[di](i, j) = inner(i, j).derivatives()[di];
    return nd;
  }

  const DualNumber<double, DuckNumberDerivatives> &
  convert_to_outer(const DualNumber<double, DuckNumberDerivatives> & inner) const
    {
      return inner;
    }

  static constexpr double _e0[9] = {1, 0, 0, 0, 0, 0, 0, 0, 0};
  static constexpr double _e1[9] = {0, 1, 0, 0, 0, 0, 0, 0, 0};
  static constexpr double _e2[9] = {0, 0, 1, 0, 0, 0, 0, 0, 0};
  static constexpr double _e3[9] = {0, 0, 0, 1, 0, 0, 0, 0, 0};
  static constexpr double _e4[9] = {0, 0, 0, 0, 1, 0, 0, 0, 0};
  static constexpr double _e5[9] = {0, 0, 0, 0, 0, 1, 0, 0, 0};
  static constexpr double _e6[9] = {0, 0, 0, 0, 0, 0, 1, 0, 0};
  static constexpr double _e7[9] = {0, 0, 0, 0, 0, 0, 0, 1, 0};
  static constexpr double _e8[9] = {0, 0, 0, 0, 0, 0, 0, 0, 1};

  metaphysicl_const_return_decl(operator());
  metaphysicl_nonconst_return_decl(operator());
  metaphysicl_const_return_decl(tr);
  metaphysicl_const_return_decl(transpose);

  metaphysicl_map_api_decl(TensorValue<double>);

private:
  metaphysicl_map_decl(TensorValue<double>);
};

metaphysicl_const_return_def(operator(), TensorValue<double>)
metaphysicl_nonconst_return_def(operator(), TensorValue<double>)
metaphysicl_const_return_def(tr, TensorValue<double>)
metaphysicl_const_return_def(transpose, TensorValue<double>)

metaphysicl_map_api_def(TensorValue<double>)

template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e0[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e1[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e2[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e3[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e4[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e5[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e6[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e7[9];
template <typename D>
constexpr double NotADuckDualNumber<TensorValue<double>, D>::_e8[9];


} // namespace MetaPhysicL
#endif
