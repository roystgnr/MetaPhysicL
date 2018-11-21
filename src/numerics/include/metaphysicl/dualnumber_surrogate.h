#ifndef METAPHYSICL_DUALNUMBER_SURROGATE_H
#define METAPHYSICL_DUALNUMBER_SURROGATE_H

#include "metaphysicl/dualnumber_surrogate_decl.h"

namespace MetaPhysicL
{

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(
    DualNumber<T, typename D::template rebind<T>::other> & dn)
  : _value(dn.value())
{
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    _derivatives[i] = &dn.derivatives()[i];
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(
    DualNumber<T, typename D::template rebind<T>::other> && dn)
  : _value(dn.value())
{
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    _derivatives[i] = &dn.derivatives()[i];
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(T & n) : _value(n)
{
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(T && n) : _value(n)
{
}

template <typename T, typename D>
template <typename T2, typename D2, class... Args>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumber<T2, D2> & dn, Args &&... args)
  : _value(dn.value()(std::forward<Args>(args)...))
{
  auto size = _derivatives.size();
  for (decltype(size) di = 0; di < size; ++di)
    _derivatives[di] = &dn.derivatives()[di](std::forward<Args>(args)...);
}

template <typename T, typename D>
template <typename T2, typename D2, class... Args>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumber<T2, D2> && dn, Args &&... args)
  : _value(dn.value()(std::forward<Args>(args)...))
{
  auto size = _derivatives.size();
  for (decltype(size) di = 0; di < size; ++di)
    _derivatives[di] = &dn.derivatives()[di](std::forward<Args>(args)...);
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumberSurrogate<T, D> & dns)
  : _value(dns.value()), _derivatives(dns.derivatives())
{
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumberSurrogate<T, D> && dns)
  : _value(dns.value()), _derivatives(dns.derivatives())
{
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const DualNumberSurrogate<T, D> & dns)
{
  _value = dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives()[i])
      *_derivatives[i] = *dns.derivatives()[i];
    else
      *_derivatives[i] = 0;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const DualNumberSurrogate<T2, D2> & dns)
{
  _value = dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives()[i])
      *_derivatives[i] = *dns.derivatives()[i];
    else
      *_derivatives[i] = 0;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(DualNumberSurrogate<T2, D2> && dns)
{
  _value = dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives()[i])
      *_derivatives[i] = *dns.derivatives()[i];
    else
      *_derivatives[i] = 0;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const T2 & in_value)
{
  _value = in_value;
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] = 0;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const DualNumber<T2, D2> & in_dn)
{
  _value = in_dn.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] = in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator+=(const DualNumberSurrogate<T2, D2> & dns)
{
  _value += dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    if (dns.derivatives()[i])
      *_derivatives[i] += *dns.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator-=(const DualNumberSurrogate<T2, D2> & dns)
{
  _value -= dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    if (dns.derivatives()[i])
      *_derivatives[i] -= *dns.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator*=(const DualNumberSurrogate<T2, D2> & dns)
{
  _value *= dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives()[i])
      *_derivatives[i] = dns.value() * *_derivatives[i] + _value * *dns.derivatives()[i];
    else
      *_derivatives[i] *= dns.value();
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator/=(const DualNumberSurrogate<T2, D2> & dns)
{
  _value /= dns.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives()[i])
      *_derivatives[i] = (dns.value() * *_derivatives[i] - _value * *dns.derivatives()[i]) /
                         (dns.value() * dns.value());
    else
      *_derivatives[i] /= dns.value();
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator+=(const T2 & in_value)
{
  _value += in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator-=(const T2 & in_value)
{
  _value -= in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator*=(const T2 & in_value)
{
  _value *= in_value;
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] *= in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator/=(const T2 & in_value)
{
  _value /= in_value;
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] /= in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator+=(const DualNumber<T2, D2> & in_dn)
{
  _value += in_dn.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] += in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator-=(const DualNumber<T2, D2> & in_dn)
{
  _value -= in_dn.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] -= in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator*=(const DualNumber<T2, D2> & in_dn)
{
  _value *= in_dn.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] = in_dn.value() * *_derivatives[i] + _value * in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator/=(const DualNumber<T2, D2> & in_dn)
{
  _value /= in_dn.value();
  auto size = _derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *_derivatives[i] = (in_dn.value() * *_derivatives[i] - _value * in_dn.derivatives()[i]) /
                       (in_dn.value() * in_dn.value);
  return *this;
}

template <typename T, typename D>
inline const T &
DualNumberSurrogate<T, D>::value() const
{
  return _value;
}

template <typename T, typename D>
inline T &
DualNumberSurrogate<T, D>::value()
{
  return _value;
}

template <typename T, typename D>
inline const D &
DualNumberSurrogate<T, D>::derivatives() const
{
  return _derivatives;
}

template <typename T, typename D>
inline D &
DualNumberSurrogate<T, D>::derivatives()
{
  return _derivatives;
}

#define metaphysicl_DNS_compares(comparator)                                                       \
  template <typename T, typename D, typename T2>                                                   \
  bool operator comparator(const DualNumberSurrogate<T, D> & dns, const T2 & solo)                 \
  {                                                                                                \
    return dns.value() comparator solo;                                                            \
  }                                                                                                \
  template <typename T2, typename T, typename D>                                                   \
  bool operator comparator(const T2 & solo, const DualNumberSurrogate<T, D> & dns)                 \
  {                                                                                                \
    return solo comparator dns.value();                                                            \
  }                                                                                                \
  template <typename T, typename D, typename T2, typename D2>                                      \
  bool operator comparator(const DualNumberSurrogate<T, D> & dns1,                                 \
                           const DualNumberSurrogate<T2, D2> & dns2)                               \
  {                                                                                                \
    return dns1.value() comparator dns2.value();                                                   \
  }

metaphysicl_DNS_compares(>)
metaphysicl_DNS_compares(<)
metaphysicl_DNS_compares(==)
metaphysicl_DNS_compares(!=)
metaphysicl_DNS_compares(>=)
metaphysicl_DNS_compares(<=)
}
#endif
