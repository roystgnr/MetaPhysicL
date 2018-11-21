#ifndef METAPHYSICL_NDDUALNUMBER_H
#define METAPHYSICL_NDDUALNUMBER_H

#include "metaphysicl/dualnumber.h"
#include "metaphysicl/dualnumber_surrogate.h"

#include <memory>

namespace MetaPhysicL
{

template <typename T, typename D>
class NotADuckDualNumber : public DualNumber<T, D>
{
public:
  using DualNumber<T, D>::DualNumber;

  NotADuckDualNumber<T, D> operator-() const
  {
    return NotADuckDualNumber<T, D>(-this->value(), -this->derivatives());
  }
  NotADuckDualNumber<T, D> operator!() const
  {
    return NotADuckDualNumber<T, D>(!this->value(), !this->derivatives());
  }
};

template <typename T, typename D>
using NDDualNumber = NotADuckDualNumber<T, D>;

#define NDDualNumber_op(opname, functorname, dn_first_calc, dn_second_calc, dualcalc)              \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline auto operator opname(const NDDualNumber<T, D> & a, const NDDualNumber<T2, D2> & b)        \
      ->NDDualNumber<decltype(a.value() opname b.value()), decltype(dualcalc)>                     \
  {                                                                                                \
    auto value = a.value() opname b.value();                                                       \
    auto derivatives = dualcalc;                                                                   \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline auto operator opname(const DualNumber<T, D> & a, const NDDualNumber<T2, D2> & b)          \
      ->NDDualNumber<decltype(a.value() opname b.value()), decltype(dualcalc)>                     \
  {                                                                                                \
    auto value = a.value() opname b.value();                                                       \
    auto derivatives = dualcalc;                                                                   \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline auto operator opname(const NDDualNumber<T, D> & a, const DualNumber<T2, D2> & b)          \
      ->NDDualNumber<decltype(a.value() opname b.value()), decltype(dualcalc)>                     \
  {                                                                                                \
    auto value = a.value() opname b.value();                                                       \
    auto derivatives = dualcalc;                                                                   \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline auto operator opname(const T & a, const NDDualNumber<T2, D> & b)                          \
      ->NDDualNumber<decltype(a opname b.value()), decltype(dn_second_calc)>                       \
  {                                                                                                \
    auto value = a opname b.value();                                                               \
    auto derivatives = dn_second_calc;                                                             \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D, typename T2>                                                   \
  inline auto operator opname(const NDDualNumber<T, D> & a, const T2 & b)                          \
      ->NDDualNumber<decltype(a.value() opname b), decltype(dn_first_calc)>                        \
  {                                                                                                \
    auto value = a.value() opname b;                                                               \
    auto derivatives = dn_first_calc;                                                              \
    return {value, derivatives};                                                                   \
  }

NDDualNumber_op(+, Plus, a.derivatives(), b.derivatives(), a.derivatives() + b.derivatives())

NDDualNumber_op(-, Minus, a.derivatives(), -b.derivatives(), a.derivatives() - b.derivatives())

NDDualNumber_op(*,
                Multiplies,
                a.derivatives() * b,
                a * b.derivatives(),
                a.value() * b.derivatives() + a.derivatives() * b.value())

NDDualNumber_op(/,
                Divides,
                a.derivatives() / b,
                -a * b.derivatives() / (b.value() * b.value()),
                (b.value() * a.derivatives() - b.derivatives() * a.value()) /
                    (b.value() * b.value()))

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_plus(const NDDualNumber<T, D> & nd, const DualNumberSurrogate<T2, DPtrs> & dns) ->
    typename D::template rebind<decltype(nd.value() + dns.value())>::other
{
  typename D::template rebind<decltype(nd.value() + dns.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = nd.derivatives()[i] + *dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_plus(const DualNumberSurrogate<T2, DPtrs> & dns, const NDDualNumber<T, D> & nd) ->
    typename D::template rebind<decltype(nd.value() + dns.value())>::other
{
  typename D::template rebind<decltype(nd.value() + dns.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = nd.derivatives()[i] + *dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_minus(const NDDualNumber<T, D> & nd, const DualNumberSurrogate<T2, DPtrs> & dns) ->
    typename D::template rebind<decltype(nd.value() - dns.value())>::other
{
  typename D::template rebind<decltype(nd.value() - dns.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = nd.derivatives()[i] - *dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_minus(const DualNumberSurrogate<T2, DPtrs> & dns, const NDDualNumber<T, D> & nd) ->
    typename D::template rebind<decltype(dns.value() - nd.value())>::other
{
  typename D::template rebind<decltype(dns.value() - nd.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i] - nd.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_multiply(const NDDualNumber<T, D> & nd, const DualNumberSurrogate<T2, DPtrs> & dns) ->
    typename D::template rebind<decltype(dns.value() * nd.value())>::other
{
  typename D::template rebind<decltype(dns.value() * nd.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i] * nd.value() + nd.derivatives()[i] * dns.value();
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_multiply(const DualNumberSurrogate<T2, DPtrs> & dns, const NDDualNumber<T, D> & nd) ->
    typename D::template rebind<decltype(dns.value() * nd.value())>::other
{
  typename D::template rebind<decltype(dns.value() * nd.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i] * nd.value() + nd.derivatives()[i] * dns.value();
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_divide(const NDDualNumber<T, D> & nd, const DualNumberSurrogate<T2, DPtrs> & dns) ->
    typename D::template rebind<decltype(nd.value() / dns.value())>::other
{
  typename D::template rebind<decltype(nd.value() / dns.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = (nd.deriatives()[i] * dns.value() - nd.value() * *dns.derivatives()[i]) /
                     (dns.value() * dns.value());
  return derivatives;
}

template <typename T, typename D, typename T2, typename DPtrs>
auto
nd_dns_divide(const DualNumberSurrogate<T2, DPtrs> & dns, const NDDualNumber<T, D> & nd) ->
    typename D::template rebind<decltype(dns.value() / nd.value())>::other
{
  typename D::template rebind<decltype(dns.value() / nd.value())>::other derivatives;
  auto size = nd.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = (*dns.deriatives()[i] * nd.value() - dns.value() * nd.derivatives()[i]) /
                     (nd.value() * nd.value());
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_plus(const T2 & t, const DualNumberSurrogate<T, D> & dns) ->
    typename D::template rebind<decltype(t + dns.value())>::other
{
  typename D::template rebind<decltype(t + dns.value())>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_plus(const DualNumberSurrogate<T, D> & dns, const T2 & t) ->
    typename D::template rebind<decltype(t + dns.value())>::other
{
  typename D::template rebind<decltype(t + dns.value())>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_minus(const T2 & t, const DualNumberSurrogate<T, D> & dns) ->
    typename D::template rebind<decltype(t - dns.value())>::other
{
  typename D::template rebind<decltype(t - dns.value())>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = -*dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_minus(const DualNumberSurrogate<T, D> & dns, const T2 & t) ->
    typename D::template rebind<decltype(dns.value() - t)>::other
{
  typename D::template rebind<decltype(dns.value() - t)>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i];
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_multiply(const T2 & t, const DualNumberSurrogate<T, D> & dns) ->
    typename D::template rebind<decltype(dns.value() * t)>::other
{
  typename D::template rebind<decltype(dns.value() * t)>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i] * t;
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_multiply(const DualNumberSurrogate<T, D> & dns, const T2 & t) ->
    typename D::template rebind<decltype(dns.value() * t)>::other
{
  typename D::template rebind<decltype(dns.value() * t)>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i] * t;
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_divide(const T2 & t, const DualNumberSurrogate<T, D> & dns) ->
    typename D::template rebind<decltype(t / dns.value())>::other
{
  typename D::template rebind<decltype(t / dns.value())>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = -t * *dns.derivatives()[i] / (dns.value() * dns.value());
  return derivatives;
}

template <typename T, typename D, typename T2>
auto
dns_t_divide(const DualNumberSurrogate<T, D> & dns, const T2 & t) ->
    typename D::template rebind<decltype(dns.value() / t)>::other
{
  typename D::template rebind<decltype(dns.value() / t)>::other derivatives;
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    derivatives[i] = *dns.derivatives()[i] / t;
  return derivatives;
}

#define ND_DNS_op(opname, optype)                                                                  \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline auto operator opname(const DualNumberSurrogate<T, D> & a, const NDDualNumber<T2, D2> & b) \
      ->NDDualNumber<decltype(a.value() opname b.value()),                                         \
                     typename D::template rebind<decltype(a.value() opname b.value())>::other>     \
  {                                                                                                \
    auto value = a.value() opname b.value();                                                       \
    auto derivatives = std::move(nd_dns_##optype(a, b));                                           \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline auto operator opname(const NDDualNumber<T, D> & a, const DualNumberSurrogate<T2, D2> & b) \
      ->NDDualNumber<decltype(a.value() opname b.value()),                                         \
                     typename D::template rebind<decltype(a.value() opname b.value())>::other>     \
  {                                                                                                \
    auto value = a.value() opname b.value();                                                       \
    auto derivatives = std::move(nd_dns_##optype(a, b));                                           \
    return {value, derivatives};                                                                   \
  }                                                                                                \
  template <typename T, typename D, typename T2>                                                   \
  inline auto operator opname(const T2 & a, const DualNumberSurrogate<T, D> & b)                   \
      ->NDDualNumber<decltype(a opname b.value()),                                                 \
                     typename D::template rebind<decltype(a opname b.value())>::other>             \
  {                                                                                                \
    auto value = a opname b.value();                                                               \
    auto derivatives = std::move(dns_t_##optype(a, b));                                            \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D, typename T2>                                                   \
  inline auto operator opname(const DualNumberSurrogate<T, D> & a, const T2 & b)                   \
      ->NDDualNumber<decltype(a.value() opname b),                                                 \
                     typename D::template rebind<decltype(a.value() opname b)>::other>             \
  {                                                                                                \
    auto value = a.value() opname b;                                                               \
    auto derivatives = std::move(dns_t_##optype(a, b));                                            \
    return {value, derivatives};                                                                   \
  }

ND_DNS_op(+, plus)
ND_DNS_op(-, minus)
ND_DNS_op(*, multiply)
ND_DNS_op(/, divide)

//
// NotADuck helper functions, method declarations, and method definitions
//

template <typename T, typename TBase, typename D, class... PtrArgs, class... ParamArgs>
void
const_void_helper(void (TBase::*fn)(PtrArgs...) const,
                  const NDDualNumber<T, D> & calling_dn,
                  ParamArgs &&... args)
{
  (calling_dn.value().*fn)(std::forward<ParamArgs>(args)...);
  auto size = calling_dn.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    (calling_dn.derivatives()[i].*fn)(std::forward<ParamArgs>(args)...);
}

template <typename T, typename TBase, typename D, class... PtrArgs, class... ParamArgs>
void
void_helper(void (TBase::*fn)(PtrArgs...), NDDualNumber<T, D> & calling_dn, ParamArgs &&... args)
{
  (calling_dn.value().*fn)(std::forward<ParamArgs>(args)...);
  auto size = calling_dn.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    (calling_dn.derivatives()[i].*fn)(std::forward<ParamArgs>(args)...);
}

#define metaphysicl_const_return_decl(method_name)                                                 \
  template <class... Args>                                                                         \
  auto method_name(Args &&... args) const->NDDualNumber<                                           \
      typename std::remove_const<typename std::remove_reference<decltype(                          \
          this->value().method_name(std::forward<Args>(args)...))>::type>::type,                   \
      typename D::template rebind<typename std::remove_const<typename std::remove_reference<       \
          decltype(this->value().method_name(std::forward<Args>(args)...))>::type>::type>::other>

#define metaphysicl_nonconst_return_decl(method_name)                                              \
  template <class... Args>                                                                         \
  auto method_name(Args &&... args)                                                                \
      ->DualNumberSurrogate<                                                                       \
          typename std::remove_reference<decltype(                                                 \
              this->value().method_name(std::forward<Args>(args)...))>::type,                      \
          typename D::template rebind<typename std::remove_reference<decltype(                     \
              this->value().method_name(std::forward<Args>(args)...))>::type *>::other> &

#define metaphysicl_const_void_decl(method_name)                                                   \
  template <class... Args>                                                                         \
  void method_name(Args &&... args) const

#define metaphysicl_nonconst_void_decl(method_name)                                                \
  template <class... Args>                                                                         \
  void method_name(Args &&... args)

#define metaphysicl_try_emplace_decl(SpecialType)                                                  \
  template <class... ConstructionArgs>                                                             \
  void dns_try_emplace(typename SpecialType::index_type, ConstructionArgs &&... args)

#define metaphysicl_at_decl(SpecialType)                                                           \
  DualNumberSurrogate<typename SpecialType::value_type,                                            \
                      typename D::template rebind<typename SpecialType::value_type *>::other> &    \
      dns_at(typename SpecialType::index_type)

#define metaphysicl_map_decl(SpecialType)                                                          \
  Map<typename SpecialType::index_type,                                                            \
      std::shared_ptr<DualNumberSurrogate<                                                         \
          typename SpecialType::value_type,                                                        \
          typename D::template rebind<typename SpecialType::value_type *>::other>>>                \
      _dual_number_surrogates

#define metaphysicl_map_api_decl(SpecialType)                                                      \
  metaphysicl_try_emplace_decl(SpecialType);                                                       \
  metaphysicl_at_decl(SpecialType)

#define metaphysicl_const_return_def(method_name, SpecialType)                                     \
  template <typename D>                                                                            \
  template <class... Args>                                                                         \
  auto NDDualNumber<SpecialType, D>::method_name(Args &&... args) const->NDDualNumber<             \
      typename std::remove_const<typename std::remove_reference<decltype(                          \
          this->value().method_name(std::forward<Args>(args)...))>::type>::type,                   \
      typename D::template rebind<typename std::remove_const<typename std::remove_reference<       \
          decltype(this->value().method_name(std::forward<Args>(args)...))>::type>::type>::other>  \
  {                                                                                                \
    typedef typename D::template rebind<typename std::remove_const<typename std::remove_reference< \
        decltype(this->value().method_name(std::forward<Args>(args)...))>::type>::type>::other     \
        DerivativeType;                                                                            \
    DerivativeType deriv;                                                                          \
    auto outer_size = this->derivatives().size();                                                  \
    auto inner_template_dn =                                                                       \
        this->convert_to_outer(this->inner_template().method_name(std::forward<Args>(args)...));   \
    auto & inner_template_value = inner_template_dn.value();                                       \
    auto & inner_template_derivatives = inner_template_dn.derivatives();                           \
    auto inner_size = inner_template_derivatives.size();                                           \
    for (decltype(outer_size) outer_di = 0; outer_di < outer_size; ++outer_di)                     \
    {                                                                                              \
      deriv[outer_di] = 0;                                                                         \
      for (decltype(inner_size) inner_di = 0; inner_di < inner_size; ++inner_di)                   \
        deriv[outer_di] +=                                                                         \
            inner_template_derivatives[inner_di] * this->derivatives()[outer_di](inner_di);        \
    }                                                                                              \
    return {inner_template_value, deriv};                                                          \
  }

#define metaphysicl_nonconst_return_def(method_name, SpecialType)                                  \
  template <typename D>                                                                            \
  template <class... Args>                                                                         \
  auto NDDualNumber<SpecialType, D>::method_name(Args &&... args)                                  \
      ->DualNumberSurrogate<                                                                       \
          typename std::remove_reference<decltype(                                                 \
              this->value().method_name(std::forward<Args>(args)...))>::type,                      \
          typename D::template rebind<typename std::remove_reference<decltype(                     \
              this->value().method_name(std::forward<Args>(args)...))>::type *>::other> &          \
  {                                                                                                \
    typename SpecialType::index_type key(std::forward<Args>(args)...);                             \
    dns_try_emplace(key, std::forward<Args>(args)...);                                             \
    return dns_at(key);                                                                            \
  }

#define metaphysicl_const_void_def(method_name, SpecialType)                                       \
  template <typename D>                                                                            \
  template <class... Args>                                                                         \
  void NDDualNumber<SpecialType, D>::method_name(Args &&... args) const                            \
  {                                                                                                \
    const_void_helper(&SpecialType::method_name, *this, std::forward<Args>(args)...);              \
  }

#define metaphysicl_nonconst_void_def(method_name, SpecialType)                                    \
  template <typename D>                                                                            \
  template <class... Args>                                                                         \
  void NDDualNumber<SpecialType, D>::method_name(Args &&... args)                                  \
                                                                                                   \
  {                                                                                                \
    void_helper(&SpecialType::method_name, *this, std::forward<Args>(args)...);                    \
  }

#define metaphysicl_try_emplace_def(SpecialType)                                                   \
  template <typename D>                                                                            \
  template <class... ConstructionArgs>                                                             \
  inline void NDDualNumber<SpecialType, D>::dns_try_emplace(typename SpecialType::index_type key,  \
                                                            ConstructionArgs &&... args)           \
  {                                                                                                \
    typedef DualNumberSurrogate<                                                                   \
        typename SpecialType::value_type,                                                          \
        typename D::template rebind<typename SpecialType::value_type *>::other>                    \
        dns_type;                                                                                  \
    typename Map<typename SpecialType::index_type, std::shared_ptr<dns_type>>::iterator it =       \
        _dual_number_surrogates.find(key);                                                         \
    if (it == _dual_number_surrogates.end())                                                       \
    {                                                                                              \
      std::shared_ptr<dns_type> dns =                                                              \
          std::make_shared<dns_type>(*this, std::forward<ConstructionArgs>(args)...);              \
      _dual_number_surrogates.emplace(key, std::move(dns));                                        \
    }                                                                                              \
  }

#define metaphysicl_at_def(SpecialType)                                                            \
  template <typename D>                                                                            \
  inline DualNumberSurrogate<                                                                      \
      typename SpecialType::value_type,                                                            \
      typename D::template rebind<typename SpecialType::value_type *>::other> &                    \
  NDDualNumber<SpecialType, D>::dns_at(typename SpecialType::index_type key)                       \
  {                                                                                                \
    return *_dual_number_surrogates.at(key);                                                       \
  }


#define metaphysicl_map_api_def(SpecialType)                                                       \
  metaphysicl_try_emplace_def(SpecialType)                                                         \
  metaphysicl_at_def(SpecialType)
}

#endif
