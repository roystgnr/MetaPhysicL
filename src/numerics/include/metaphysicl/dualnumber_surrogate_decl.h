#ifndef METAPHYSICL_DUALNUMBER_SURROGATE_DECL_H
#define METAPHYSICL_DUALNUMBER_SURROGATE_DECL_H

namespace MetaPhysicL
{

template <typename T, typename D, typename Enable>
class NotADuckDualNumber;
template <typename T, typename D = T, typename Enable = void>
using NDDualNumber = NotADuckDualNumber<T, D, Enable>;

// Surrogate structure that refers/points to a component of a DualNumber, e.g. it might correspond
// to the ij'th component of a DualNumber<TypeTensor<T>, N>

template <typename T, typename D>
struct DualNumberSurrogate
{
  DualNumberSurrogate(DualNumber<T, D> & dn);

  DualNumberSurrogate(DualNumber<T, D> && dn);

  DualNumberSurrogate(const T & n);

  template <typename T2, typename D2, class... Args>
  DualNumberSurrogate(DualNumber<T2, D2> & dn, Args &&... args);

  template <typename T2, typename D2, class... Args>
  DualNumberSurrogate(DualNumber<T2, D2> && dn, Args &&... args);

  DualNumberSurrogate(DualNumberSurrogate<T, D> & dns);

  DualNumberSurrogate(const DualNumberSurrogate<T, D> & dns);

  DualNumberSurrogate(DualNumberSurrogate<T, D> && dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator=(DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator=(const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator=(DualNumberSurrogate<T2, D2> && dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator=(const T2 & in_value);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator=(const DualNumber<T2, D2> & in_dn);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator+=(const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator-=(const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator*=(const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator/=(const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator+=(const T2 & in_value);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator-=(const T2 & in_value);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator*=(const T2 & in_value);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator/=(const T2 & in_value);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator+=(const DualNumber<T2, D2> & in_dn);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator-=(const DualNumber<T2, D2> & in_dn);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator*=(const DualNumber<T2, D2> & in_dn);

  template <typename T2, typename D2>
  DualNumberSurrogate<T, D> & operator/=(const DualNumber<T2, D2> & in_dn);

  T & value;
  D derivatives;
};

#define metaphysicl_DNS_compares_decl(comparator)                                                  \
  template <typename T, typename D, typename T2>                                                   \
  bool operator comparator(const DualNumberSurrogate<T, D> & dns, const T2 & solo);                \
                                                                                                   \
  template <typename T2, typename T, typename D>                                                   \
  bool operator comparator(const T2 & solo, const DualNumberSurrogate<T, D> & dns);                \
                                                                                                   \
  template <typename T, typename D, typename T2, typename D2>                                      \
  bool operator comparator(const DualNumberSurrogate<T, D> & dns1,                                 \
                           const DualNumberSurrogate<T2, D2> & dns2)

metaphysicl_DNS_compares_decl(>);
metaphysicl_DNS_compares_decl(<);
metaphysicl_DNS_compares_decl(==);
metaphysicl_DNS_compares_decl(!=);
metaphysicl_DNS_compares_decl(>=);
metaphysicl_DNS_compares_decl(<=);
}
#endif
