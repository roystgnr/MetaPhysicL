#ifndef METAPHYSICL_MATH_STRUCTS_H
#define METAPHYSICL_MATH_STRUCTS_H

#include <tuple>
#include <map>
#include "metaphysicl/compare_types.h"
#include "metaphysicl/dualnumber.h"

namespace MetaPhysicL
{

class RealVector
{
public:
  RealVector();

  RealVector(const double x, const double y = 0, const double z = 0);

  typedef double value_type;
  typedef std::tuple<unsigned int> index_type;

  RealVector operator*(const double &)const;
  double operator*(const RealVector &)const;
  RealVector & operator*=(const double &);

  RealVector operator+(const RealVector &) const;
  RealVector operator-(const RealVector &) const;
  RealVector operator-() const;

  const double & operator()(unsigned) const;
  double & operator()(unsigned);

private:
  double _coords[3];
};

RealVector operator*(const double &, const RealVector &);

template <bool reverseorder>
struct PlusType<RealVector, RealVector, reverseorder>
{
  typedef RealVector supertype;
};

template <bool reverseorder>
struct MinusType<RealVector, RealVector, reverseorder>
{
  typedef RealVector supertype;
};

template <bool reverseorder>
struct MultipliesType<RealVector, RealVector, reverseorder>
{
  typedef double supertype;
};

template <bool reverseorder>
struct MultipliesType<double, RealVector, reverseorder>
{
  typedef RealVector supertype;
};

template <>
struct DividesType<RealVector, double>
{
  typedef RealVector supertype;
};

class RealTensor
{
public:
  RealTensor();

  RealTensor(const double xx,
             const double xy = 0,
             const double xz = 0,
             const double yx = 0,
             const double yy = 0,
             const double yz = 0,
             const double zx = 0,
             const double zy = 0,
             const double zz = 0);

  RealTensor(const RealVector & row1, const RealVector & row2, const RealVector & row3);

  typedef double value_type;
  typedef std::tuple<unsigned int, unsigned int> index_type;

  RealTensor operator*(const double &)const;
  RealTensor operator/(const double &) const;
  RealVector operator*(const RealVector &)const;
  RealTensor operator*(const RealTensor &)const;
  RealTensor & operator*=(const double &);
  RealTensor & operator*=(const RealTensor &);

  RealTensor operator+(const RealTensor &) const;
  RealTensor operator-(const RealTensor &) const;
  RealTensor operator-() const;

  const double & operator()(unsigned, unsigned) const;
  double & operator()(unsigned, unsigned);

  RealTensor transpose() const;
  double tr() const;

private:
  double _coords[9];
};

RealTensor operator*(const double &, const RealTensor &);

template <bool reverseorder>
struct PlusType<RealTensor, RealTensor, reverseorder>
{
  typedef RealTensor supertype;
};

template <bool reverseorder>
struct MinusType<RealTensor, RealTensor, reverseorder>
{
  typedef RealTensor supertype;
};

template <bool reverseorder>
struct MultipliesType<RealTensor, RealTensor, reverseorder>
{
  typedef RealTensor supertype;
};

template <bool reverseorder>
struct MultipliesType<double, RealTensor, reverseorder>
{
  typedef RealTensor supertype;
};

template <>
struct MultipliesType<RealTensor, RealVector>
{
  typedef RealVector supertype;
};

template <>
struct DividesType<RealTensor, double>
{
  typedef RealTensor supertype;
};

template <typename D>
class NotADuckDualNumber<RealTensor, D> : public DualNumber<RealTensor, D>
{
public:
  NotADuckDualNumber() : DualNumber<RealTensor, D>() {}

  NotADuckDualNumber(
      const NotADuckDualNumber<RealVector, typename D::template rebind<RealVector>::other> & row1,
      const NotADuckDualNumber<RealVector, typename D::template rebind<RealVector>::other> & row2,
      const NotADuckDualNumber<RealVector, typename D::template rebind<RealVector>::other> & row3)
    : DualNumber<RealTensor, D>()
  {
    this->value() = RealTensor(row1.value(), row2.value(), row3.value());
    auto size = this->derivatives().size();
    for (decltype(size) i = 0; i < size; ++i)
      this->derivatives()[i] =
          RealTensor(row1.derivatives()[i], row2.derivatives()[i], row3.derivatives()[i]);
  }

  using DualNumber<RealTensor, D>::DualNumber;

  template <typename KeyType, typename ValueType>
  using Map = std::map<KeyType, ValueType>;

  NotADuckDualNumber<RealTensor, D> operator-() const
  {
    return NotADuckDualNumber<RealTensor, D>(-this->_val, -this->_deriv);
  }
  NotADuckDualNumber<RealTensor, D> operator!() const
  {
    return NotADuckDualNumber<RealTensor, D>(!this->_val, !this->_deriv);
  }

  metaphysicl_const_return_decl(operator());
  metaphysicl_nonconst_return_decl(operator());
  metaphysicl_const_return_decl(tr);
  metaphysicl_const_return_decl(transpose);

  metaphysicl_map_api_decl(RealTensor);

private:
  metaphysicl_map_decl(RealTensor);
};

metaphysicl_const_return_def(operator(), RealTensor);
metaphysicl_nonconst_return_def(operator(), RealTensor);
metaphysicl_const_return_def(tr, RealTensor);
metaphysicl_const_return_def(transpose, RealTensor);

metaphysicl_map_api_def(RealTensor);
} // namespace MetaPhysicL
#endif
