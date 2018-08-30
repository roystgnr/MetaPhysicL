#include "math_structs.h"

namespace MetaPhysicL
{

RealVector::RealVector() : _coords() {}

RealVector::RealVector(const double x, const double y, const double z)
{
  _coords[0] = x;
  _coords[1] = y;
  _coords[2] = z;
}

const double &
RealVector::operator()(unsigned i) const
{
  return _coords[i];
}

double &
RealVector::operator()(unsigned i)
{
  return _coords[i];
}

RealVector RealVector::operator*(const double & scalar) const
{
  return {_coords[0] * scalar, _coords[1] * scalar, _coords[2] * scalar};
}

double RealVector::operator*(const RealVector & vector) const
{
  double ret = 0;
  for (unsigned i = 0; i < 3; ++i)
    ret += _coords[i] * vector(i);
  return ret;
}

RealVector &
RealVector::operator*=(const double & scalar)
{
  for (unsigned i = 0; i < 3; ++i)
    _coords[i] *= scalar;
  return *this;
}

RealVector
RealVector::operator+(const RealVector & vector) const
{
  return {_coords[0] + vector(0), _coords[1] + vector(1), _coords[2] + vector(2)};
}

RealVector
RealVector::operator-(const RealVector & vector) const
{
  return {_coords[0] - vector(0), _coords[1] - vector(1), _coords[2] - vector(2)};
}

RealVector
RealVector::operator-() const
{
  return {-_coords[0], -_coords[1], -_coords[2]};
}

RealTensor::RealTensor() : _coords() {}

RealTensor::RealTensor(const double xx,
                       const double xy,
                       const double xz,
                       const double yx,
                       const double yy,
                       const double yz,
                       const double zx,
                       const double zy,
                       const double zz)
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

RealTensor::RealTensor(const RealVector & row1, const RealVector & row2, const RealVector & row3)
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

const double &
RealTensor::operator()(unsigned i, unsigned j) const
{
  return _coords[i * 3 + j];
}

double &
RealTensor::operator()(unsigned i, unsigned j)
{
  return _coords[i * 3 + j];
}

RealTensor RealTensor::operator*(const double & scalar) const
{
  RealTensor ret;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      ret(i, j) = (*this)(i, j) * scalar;
  return ret;
}

RealTensor
RealTensor::operator/(const double & scalar) const
{
  RealTensor ret;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      ret(i, j) = (*this)(i, j) / scalar;
  return ret;
}

RealVector RealTensor::operator*(const RealVector & vec) const
{
  RealVector ret;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      ret(i) += (*this)(i, j) * vec(j);

  return ret;
}

RealTensor RealTensor::operator*(const RealTensor & tensor) const
{
  RealTensor ret;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        ret(i, j) += (*this)(i, k) * tensor(k, j);

  return ret;
}

RealTensor &
RealTensor::operator*=(const double & scalar)
{
  for (unsigned i = 0; i < 9; ++i)
    _coords[i] *= scalar;
  return *this;
}

RealTensor &
RealTensor::operator*=(const RealTensor & tensor)
{
  RealTensor ret;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        ret(i, j) += (*this)(i, k) * tensor(k, j);

  *this = ret;
  return *this;
}

RealTensor
RealTensor::operator+(const RealTensor & tensor) const
{
  RealTensor ret;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      ret(i, j) = (*this)(i, j) + tensor(i, j);
  return ret;
}

RealTensor
RealTensor::operator-(const RealTensor & tensor) const
{
  RealTensor ret;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      ret(i, j) = (*this)(i, j) - tensor(i, j);
  return ret;
}

RealTensor
RealTensor::operator-() const
{
  RealTensor ret;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      ret(i, j) = -(*this)(i, j);
  return ret;
}

RealTensor
RealTensor::transpose() const
{
  return RealTensor(_coords[0],
                    _coords[3],
                    _coords[6],
                    _coords[1],
                    _coords[4],
                    _coords[7],
                    _coords[2],
                    _coords[5],
                    _coords[8]);
}

double
RealTensor::tr() const
{
  return _coords[0] + _coords[4] + _coords[8];
}

RealVector operator*(const double & scalar, const RealVector & vector)
{
  return {scalar * vector(0), scalar * vector(1), scalar * vector(2)};
}

RealTensor operator*(const double & scalar, const RealTensor & tensor)
{
  RealTensor ret;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      ret(i, j) = scalar * tensor(i, j);

  return ret;
}
}
