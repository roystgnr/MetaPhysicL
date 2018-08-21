#ifndef METAPHYSICL_DUALNUMBER_SURROGATE_H
#define METAPHYSICL_DUALNUMBER_SURROGATE_H

#include "metaphysicl/dualnumber_surrogate_decl.h"
#include "metaphysicl/dualnumber_decl.h"

namespace MetaPhysicL
{

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumber<T, D> & dn) : value(dn.value())
{
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = &dn.derivatives()[i];
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumber<T, D> && dn) : value(dn.value())
{
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = &dn.derivatives()[i];
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(const T & n) : value(n)
{
}

template <typename T, typename D>
template <typename T2, typename D2, class... Args>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumber<T2, D2> & dn, Args &&... args)
  : value(dn.value()(std::forward<Args>(args)...))
{
  auto size = derivatives.size();
  for (decltype(size) di = 0; di < size; ++di)
    derivatives[di] = &dn.derivatives()[di](std::forward<Args>(args)...);
}

template <typename T, typename D>
template <typename T2, typename D2, class... Args>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumber<T2, D2> && dn, Args &&... args)
  : value(dn.value()(std::forward<Args>(args)...))
{
  auto size = derivatives.size();
  for (decltype(size) di = 0; di < size; ++di)
    derivatives[di] = &dn.derivatives()[di](std::forward<Args>(args)...);
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumberSurrogate<T, D> & dns)
  : value(dns.value), derivatives(dns.derivatives)
{
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(const DualNumberSurrogate<T, D> & dns)
  : value(dns.value), derivatives(dns.derivatives)
{
}

template <typename T, typename D>
inline DualNumberSurrogate<T, D>::DualNumberSurrogate(DualNumberSurrogate<T, D> && dns)
  : value(dns.value), derivatives(dns.derivatives)
{
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(DualNumberSurrogate<T2, D2> & dns)
{
  value = dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives[i])
      *derivatives[i] = *dns.derivatives[i];
    else
      *derivatives[i] = 0;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const DualNumberSurrogate<T2, D2> & dns)
{
  value = dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives[i])
      *derivatives[i] = *dns.derivatives[i];
    else
      *derivatives[i] = 0;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(DualNumberSurrogate<T2, D2> && dns)
{
  value = dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives[i])
      *derivatives[i] = *dns.derivatives[i];
    else
      *derivatives[i] = 0;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const T2 & in_value)
{
  value = in_value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] = 0;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator=(const DualNumber<T2, D2> & in_dn)
{
  value = in_dn.value();
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] = in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator+=(const DualNumberSurrogate<T2, D2> & dns)
{
  value += dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    if (dns.derivatives[i])
      *derivatives[i] += *dns.derivatives[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator-=(const DualNumberSurrogate<T2, D2> & dns)
{
  value -= dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    if (dns.derivatives[i])
      *derivatives[i] -= *dns.derivatives[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator*=(const DualNumberSurrogate<T2, D2> & dns)
{
  value *= dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives[i])
      *derivatives[i] = dns.value * *derivatives[i] + value * *dns.derivatives[i];
    else
      *derivatives[i] *= dns.value;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator/=(const DualNumberSurrogate<T2, D2> & dns)
{
  value /= dns.value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
  {
    if (dns.derivatives[i])
      *derivatives[i] =
          (dns.value * *derivatives[i] - value * *dns.derivatives[i]) / (dns.value * dns.value);
    else
      *derivatives[i] /= dns.value;
  }
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator+=(const T2 & in_value)
{
  value += in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator-=(const T2 & in_value)
{
  value -= in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator*=(const T2 & in_value)
{
  value *= in_value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] *= in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator/=(const T2 & in_value)
{
  value /= in_value;
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] /= in_value;
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator+=(const DualNumber<T2, D2> & in_dn)
{
  value += in_dn.value();
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] += in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator-=(const DualNumber<T2, D2> & in_dn)
{
  value -= in_dn.value();
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] -= in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator*=(const DualNumber<T2, D2> & in_dn)
{
  value *= in_dn.value();
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] = in_dn.value() * *derivatives[i] + value * in_dn.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberSurrogate<T, D> &
DualNumberSurrogate<T, D>::operator/=(const DualNumber<T2, D2> & in_dn)
{
  value /= in_dn.value();
  auto size = derivatives.size();
  for (decltype(size) i = 0; i < size; ++i)
    *derivatives[i] = (in_dn.value() * *derivatives[i] - value * in_dn.derivatives()[i]) /
                      (in_dn.value() * in_dn.value);
  return *this;
}

#define metaphysicl_DNS_compares(comparator)                                                       \
  template <typename T, typename D, typename T2>                                                   \
  bool operator comparator(const DualNumberSurrogate<T, D> & dns, const T2 & solo)                 \
  {                                                                                                \
    return dns.value comparator solo;                                                              \
  }                                                                                                \
  template <typename T2, typename T, typename D>                                                   \
  bool operator comparator(const T2 & solo, const DualNumberSurrogate<T, D> & dns)                 \
  {                                                                                                \
    return solo comparator dns.value;                                                              \
  }                                                                                                \
  template <typename T, typename D, typename T2, typename D2>                                      \
  bool operator comparator(const DualNumberSurrogate<T, D> & dns1,                                 \
                           const DualNumberSurrogate<T2, D2> & dns2)                               \
  {                                                                                                \
    return dns1.value comparator dns2.value;                                                       \
  }                                                                                                \
  void macro_syntax_function()

metaphysicl_DNS_compares(>);
metaphysicl_DNS_compares(<);
metaphysicl_DNS_compares(==);
metaphysicl_DNS_compares(!=);
metaphysicl_DNS_compares(>=);
metaphysicl_DNS_compares(<=);
}
#endif
