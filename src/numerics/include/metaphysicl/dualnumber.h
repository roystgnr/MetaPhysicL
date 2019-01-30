//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MetaPhysicL - A metaprogramming library for physics calculations
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------


#ifndef METAPHYSICL_DUALNUMBER_H
#define METAPHYSICL_DUALNUMBER_H

#include "metaphysicl/dualnumber_decl.h"
#include "metaphysicl/dualnumber_surrogate.h"

namespace MetaPhysicL {

template <typename T, typename D>
class NotADuckDualNumber;

template <typename T, typename D>
inline
T&
DualNumber<T,D>::value() { return _val; }

template <typename T, typename D>
inline
const T&
DualNumber<T,D>::value() const { return _val; }

template <typename T, typename D>
inline
D&
DualNumber<T,D>::derivatives() { return _deriv; }

template <typename T, typename D>
inline
const D&
DualNumber<T,D>::derivatives() const { return _deriv; }

template <typename T, typename D>
inline
bool
DualNumber<T,D>::boolean_test() const { return _val; }

template <typename T, typename D>
inline
DualNumber<T,D>
DualNumber<T,D>::operator- () const { return DualNumber<T,D>(-_val, -_deriv); }

template <typename T, typename D>
inline
DualNumber<T,D>
DualNumber<T,D>::operator! () const { return DualNumber<T,D>(!_val, !_deriv); }

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D> &
DualNumber<T,D>::operator=(const DualNumber<T2,D2> & dn)
{
  _val = dn.value();
  _deriv = dn.derivatives();
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D> &
DualNumber<T,D>::operator=(const NotADuckDualNumber<T2,D2> & nd_dn)
{
  _val = nd_dn.value();
  _deriv = nd_dn.derivatives();
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>::DualNumber(const DualNumberSurrogate<T2,D2> & dns) : _val(dns.value())
{
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    _deriv[i] = *dns.derivatives()[i];
}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D> &
DualNumber<T,D>::operator=(const DualNumberSurrogate<T2,D2> & dns)
{
  _val = dns.value();
  auto size = dns.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    _deriv[i] = *dns.derivatives()[i];
  return *this;
}

template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D> &
DualNumber<T,D>::operator=(const T2 & scalar)
{
  _val = scalar;
  _deriv = 0;
  return *this;
}

//
// Member function definitions
//
template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D>::DualNumber(const T2& val) :
  _val  (DualNumberConstructor<T,D>::value(val)),
  _deriv(DualNumberConstructor<T,D>::deriv(val)) {}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>::DualNumber(const T2& val,
                            const D2& deriv) :
  _val  (DualNumberConstructor<T,D>::value(val,deriv)),
  _deriv(DualNumberConstructor<T,D>::deriv(val,deriv)) {}

// FIXME: these operators currently do automatic type promotion when
// encountering DualNumbers of differing levels of recursion and
// differentiability.  But what we really want is automatic type
// *demotion*, to avoid pretending we have accurate derivatives which
// we don't have.  If we could do that right then some potential
// subtle run-time user errors would turn into compile-time user
// errors.

#define DualNumber_preop(opname, functorname, simplecalc, dualcalc) \
template <typename T, typename D> \
template <typename T2> \
inline \
DualNumber<T,D>& \
DualNumber<T,D>::operator opname##= (const T2& in) \
{ \
  simplecalc; \
  this->value() opname##= in; \
  return *this; \
} \
 \
template <typename T, typename D> \
template <typename T2, typename D2> \
inline \
DualNumber<T,D>& \
DualNumber<T,D>::operator opname##= (const DualNumber<T2,D2>& in) \
{ \
  dualcalc; \
  this->value() opname##= in.value(); \
  return *this; \
} \
 \
template <typename T, typename D> \
template <typename T2, typename D2> \
inline \
DualNumber<T,D> & \
DualNumber<T,D>::operator opname##= (const NotADuckDualNumber<T2,D2>& in) \
{ \
  dualcalc; \
  this->value() opname##= in.value(); \
  return *this; \
} \
 \
template <typename T, typename D, typename T2, typename D2> \
inline \
typename functorname##Type<DualNumber<T,D>,DualNumber<T2,D2> >::supertype \
operator opname (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  typedef typename \
    functorname##Type<DualNumber<T,D>,DualNumber<T2,D2> >::supertype \
    DS; \
  DS returnval = a; \
  returnval opname##= b; \
  return returnval; \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename functorname##Type<DualNumber<T2,D>,T,true>::supertype \
operator opname (const T& a, const DualNumber<T2,D>& b) \
{ \
  typedef typename \
    functorname##Type<DualNumber<T2,D>,T,true>::supertype DS; \
  DS returnval = a; \
  returnval opname##= b; \
  return returnval; \
} \
 \
template <typename T, typename D, typename T2> \
inline \
typename functorname##Type<DualNumber<T,D>,T2,false>::supertype \
operator opname (const DualNumber<T,D>& a, const T2& b) \
{ \
  typedef typename \
    functorname##Type<DualNumber<T,D>,T2,false>::supertype DS; \
  DS returnval = a; \
  returnval opname##= b; \
  return returnval; \
}



// With C++11, define "move operations" where possible.  We should be
// more complete and define the move-from-b alternatives as well, but
// those would require additional support to correctly handle
// division, subtraction, or non-commutative addition/multiplication
#if __cplusplus >= 201103L
#define DualNumber_op(opname, functorname, simplecalc, dualcalc) \
        DualNumber_preop(opname, functorname, simplecalc, dualcalc) \
 \
template <typename T, typename D, typename T2, typename D2> \
inline \
typename functorname##Type<DualNumber<T,D>,DualNumber<T2,D2> >::supertype \
operator opname (DualNumber<T,D>&& a, const DualNumber<T2,D2>& b) \
{ \
  typedef typename \
    functorname##Type<DualNumber<T,D>,DualNumber<T2,D2> >::supertype \
    DS; \
  DS returnval = std::move(a); \
  returnval opname##= b; \
  return returnval; \
} \
 \
template <typename T, typename D, typename T2> \
inline \
typename functorname##Type<DualNumber<T,D>,T2,false>::supertype \
operator opname (DualNumber<T,D>&& a, const T2& b) \
{ \
  typedef typename \
    functorname##Type<DualNumber<T,D>,T2,false>::supertype DS; \
  DS returnval = std::move(a); \
  returnval opname##= b; \
  return returnval; \
}

#else
#define DualNumber_op(opname, functorname, simplecalc, dualcalc) \
        DualNumber_preop(opname, functorname, simplecalc, dualcalc)
#endif

DualNumber_op(+, Plus, , this->derivatives() += in.derivatives())

DualNumber_op(-, Minus, , this->derivatives() -= in.derivatives())

DualNumber_op(*,
              Multiplies,
              this->derivatives() *= in,
              this->derivatives() = this->derivatives() * in.value() +
                                    this->value() * in.derivatives())

DualNumber_op(/,
              Divides,
              this->derivatives() /= in,
              this->derivatives() = this->derivatives() / in.value() -
                                    in.derivatives() * this->value() /
                                        (in.value() * in.value()))


#define DualNumber_compare(opname)                          \
template <typename T, typename D, typename T2, typename D2> \
inline \
bool \
operator opname  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  return (a.value() opname b.value()); \
} \
 \
template <typename T, typename T2, typename D2> \
inline \
typename boostcopy::enable_if_class< \
  typename CompareTypes<DualNumber<T2,D2>,T>::supertype, \
  bool \
>::type \
operator opname  (const T& a, const DualNumber<T2,D2>& b) \
{ \
  return (a opname b.value()); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename boostcopy::enable_if_class< \
  typename CompareTypes<DualNumber<T,D>,T2>::supertype, \
  bool \
>::type \
operator opname  (const DualNumber<T,D>& a, const T2& b) \
{ \
  return (a.value() opname b); \
}

    DualNumber_compare(>) DualNumber_compare(>=) DualNumber_compare(<) DualNumber_compare(<=)
        DualNumber_compare(==) DualNumber_compare(!=) DualNumber_compare(&&) DualNumber_compare(||)

            template <typename T, typename D>
            inline std::ostream &
            operator<<(std::ostream & output, const DualNumber<T, D> & a)
{
  return output << '(' << a.value() << ',' << a.derivatives() << ')';
}


template <typename T, typename D>
inline
D gradient(const DualNumber<T, D>& a)
{
  return a.derivatives();
}
} // namespace MetaPhysicL


namespace std {

using MetaPhysicL::DualNumber;
using MetaPhysicL::CompareTypes;

template <typename T, typename D>
inline bool isnan (const DualNumber<T,D> & a)
{
  using std::isnan;
  return isnan(a.value());
}


#if __cplusplus >= 201103L
#define DualNumber_std_unary(funcname, derivative, precalc) \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname (const DualNumber<T,D> & in) \
{ \
  DualNumber<T,D> returnval = in; \
  T funcval = std::funcname(in.value()); \
  precalc; \
  returnval.derivatives() *= derivative; \
  returnval.value() = funcval; \
  return returnval; \
} \
 \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname (DualNumber<T,D> && in) \
{ \
  T funcval = std::funcname(in.value()); \
  precalc; \
  in.derivatives() *= derivative; \
  in.value() = funcval; \
  return std::move(in); \
}

#define DualNumber_equiv_unary(funcname, equivalent) \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname (const DualNumber<T,D> & in) \
{ \
  return std::equivalent(in); \
} \
 \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname (DualNumber<T,D> && in) \
{ \
  return std::equivalent(in); \
}

#else

#define DualNumber_std_unary(funcname, derivative, precalc) \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname (DualNumber<T,D> in) \
{ \
  T funcval = std::funcname(in.value()); \
  precalc; \
  in.derivatives() *= derivative; \
  in.value() = funcval; \
  return std::move(in); \
}

#define DualNumber_equiv_unary(funcname, equivalent) \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname (DualNumber<T,D> in) \
{ \
  return std::equivalent(in); \
}
 
#endif

#define DualNumber_equivfl_unary(funcname) \
DualNumber_equiv_unary(funcname##f, funcname) \
DualNumber_equiv_unary(funcname##l, funcname)


DualNumber_std_unary(sqrt, 1 / (2 * funcval),)
DualNumber_std_unary(exp, funcval,)
DualNumber_std_unary(log, 1 / in.value(),)
DualNumber_std_unary(log10, 1 / in.value() * (1/std::log(T(10.))),)
DualNumber_std_unary(sin, std::cos(in.value()),)
DualNumber_std_unary(cos, -std::sin(in.value()),)
DualNumber_std_unary(tan, sec_in * sec_in, T sec_in = 1 / std::cos(in.value()))
DualNumber_std_unary(asin, 1 / std::sqrt(1 - in.value()*in.value()),)
DualNumber_std_unary(acos, -1 / std::sqrt(1 - in.value()*in.value()),)
DualNumber_std_unary(atan, 1 / (1 + in.value()*in.value()),)
DualNumber_std_unary(sinh, std::cosh(in.value()),)
DualNumber_std_unary(cosh, std::sinh(in.value()),)
DualNumber_std_unary(tanh, sech_in * sech_in, T sech_in = 1 / std::cosh(in.value()))
DualNumber_std_unary(abs, (in.value() > 0) - (in.value() < 0),) // std < and > return 0 or 1
DualNumber_equiv_unary(fabs, abs)
DualNumber_std_unary(norm, 2. * in.value(),)
DualNumber_std_unary(ceil, 0,)
DualNumber_std_unary(floor, 0,)

#if __cplusplus >= 201103L
DualNumber_equiv_unary(llabs, abs)
DualNumber_equiv_unary(imaxabs, abs)
DualNumber_equivfl_unary(fabs)
DualNumber_equivfl_unary(exp)
DualNumber_std_unary(exp2, std::log(T(2))*funcval,)
DualNumber_equivfl_unary(exp2)
DualNumber_std_unary(expm1, std::exp(in.value()),)
DualNumber_equivfl_unary(expm1)
DualNumber_equivfl_unary(log)
DualNumber_equivfl_unary(log10)
DualNumber_std_unary(log2, 1 / in.value() * (1/std::log(T(2))),)
DualNumber_equivfl_unary(log2)
DualNumber_std_unary(log1p, 1 / (in.value() + 1),)
DualNumber_equivfl_unary(log1p)
DualNumber_equivfl_unary(sqrt)
DualNumber_std_unary(cbrt, 1 / (3 * funcval * funcval),)
DualNumber_equivfl_unary(cbrt)
DualNumber_equivfl_unary(sin)
DualNumber_equivfl_unary(cos)
DualNumber_equivfl_unary(tan)
DualNumber_equivfl_unary(asin)
DualNumber_equivfl_unary(acos)
DualNumber_equivfl_unary(atan)
DualNumber_equivfl_unary(sinh)
DualNumber_equivfl_unary(cosh)
DualNumber_equivfl_unary(tanh)
DualNumber_std_unary(asinh, 1 / sqrt(1 + in.value()*in.value()),)
DualNumber_equivfl_unary(asinh)
DualNumber_std_unary(acosh, 1 / sqrt(in.value()*in.value() - 1),)
DualNumber_equivfl_unary(acosh)
DualNumber_std_unary(atanh, 1 / (1 - in.value()*in.value()),)
DualNumber_equivfl_unary(atanh)
// 2/sqrt(pi) = 1/sqrt(atan(1.0))
DualNumber_std_unary(erf, 1/sqrt(atan(T(1)))*exp(-in.value()*in.value()),)
DualNumber_equivfl_unary(erf)
DualNumber_std_unary(erfc, -1/sqrt(atan(T(1)))*exp(-in.value()*in.value()),)
DualNumber_equivfl_unary(erfc)
// FIXME: how do we differentiate tgamma and lgamma without boost?
DualNumber_equivfl_unary(ceil)
DualNumber_equivfl_unary(floor)
DualNumber_std_unary(trunc, 0,)
DualNumber_equivfl_unary(trunc)
DualNumber_std_unary(round, 0,)
DualNumber_equivfl_unary(round)
DualNumber_std_unary(nearbyint, 0,)
DualNumber_equivfl_unary(nearbyint)
DualNumber_std_unary(rint, 0,)
DualNumber_equivfl_unary(rint)
#endif // __cplusplus >= 201103L

#define DualNumber_complex_std_unary_real(funcname) \
template <typename T, typename D> \
inline DualNumber<T, typename D::template rebind<T>::other> \
funcname(const DualNumber<std::complex<T>, D> & in) \
{ \
  return {std::funcname(in.value()), std::numeric_limits<double>::quiet_NaN()}; \
} \
template <typename T> \
inline DualNumber<T> \
funcname(const DualNumber<std::complex<T>> & in)    \
{ \
  return {std::funcname(in.value()), std::numeric_limits<double>::quiet_NaN()}; \
}

DualNumber_complex_std_unary_real(real)
DualNumber_complex_std_unary_real(imag)
DualNumber_complex_std_unary_real(norm)
DualNumber_complex_std_unary_real(abs)

#define DualNumber_complex_std_unary_complex_pre(funcname) \
template <typename T, typename D> \
inline DualNumber<std::complex<T>, D> \
funcname(const DualNumber<std::complex<T>, D> & in) \
{ \
  return {std::funcname(in.value()), std::complex<T>{std::numeric_limits<double>::quiet_NaN(), \
                                                     std::numeric_limits<double>::quiet_NaN()}}; \
} \
template <typename T> \
inline DualNumber<std::complex<T>> \
funcname(const DualNumber<std::complex<T>> & in) \
{ \
  return {std::funcname(in.value()), std::complex<T>{std::numeric_limits<double>::quiet_NaN(), \
                                                     std::numeric_limits<double>::quiet_NaN()}}; \
}

#if __cplusplus >= 201103L
#define DualNumber_complex_std_unary_complex(funcname) \
DualNumber_complex_std_unary_complex_pre(funcname) \
template <typename T, typename D> \
inline DualNumber<std::complex<T>, D> \
funcname(DualNumber<std::complex<T>, D> && in) \
{ \
  in.value() = std::funcname(in.value()); \
  in.derivatives() = std::complex<T>(std::numeric_limits<double>::quiet_NaN(), \
                                     std::numeric_limits<double>::quiet_NaN()); \
  return in; \
} \
template <typename T> \
inline DualNumber<std::complex<T>> \
funcname(DualNumber<std::complex<T>> && in) \
{ \
  in.value() = std::funcname(in.value()); \
  in.derivatives() = std::complex<T>(std::numeric_limits<double>::quiet_NaN(), \
                                     std::numeric_limits<double>::quiet_NaN()); \
  return in; \
}
#else
#define DualNumber_complex_std_unary_complex(funcname) \
DualNumber_complex_std_unary_complex_pre(funcname)
#endif

DualNumber_complex_std_unary_complex(conj)

#define DualNumber_std_binary(funcname, derivative) \
template <typename T, typename D, typename T2, typename D2> \
inline \
typename CompareTypes<DualNumber<T,D>,DualNumber<T2,D2> >::supertype \
funcname (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<DualNumber<T,D>,DualNumber<T2,D2> >::supertype type; \
 \
  TS funcval = std::funcname(a.value(), b.value()); \
  return type(funcval, derivative); \
} \
 \
template <typename T, typename D> \
inline \
DualNumber<T,D> \
funcname (const DualNumber<T,D>& a, const DualNumber<T,D>& b) \
{ \
  T funcval = std::funcname(a.value(), b.value()); \
  return DualNumber<T,D>(funcval, derivative); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename CompareTypes<DualNumber<T2,D>,T,true>::supertype \
funcname (const T& a, const DualNumber<T2,D>& b) \
{ \
  typedef typename CompareTypes<DualNumber<T2,D>,T,true>::supertype type; \
  type newa(a); \
  return std::funcname(newa, b); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename CompareTypes<DualNumber<T,D>,T2>::supertype \
funcname (const DualNumber<T,D>& a, const T2& b) \
{ \
  typedef typename CompareTypes<DualNumber<T,D>,T2>::supertype type; \
  type newb(b); \
  return std::funcname(a, newb); \
}

#define DualNumber_equiv_binary(funcname, equivalent) \
template <typename T, typename D, typename T2, typename D2> \
inline \
typename CompareTypes<DualNumber<T,D>,DualNumber<T2,D2> >::supertype \
funcname (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  return std::equivalent(a,b); \
} \
 \
template <typename T, typename D> \
inline \
DualNumber<T,D> \
funcname (const DualNumber<T,D>& a, const DualNumber<T,D>& b) \
{ \
  return std::equivalent(a,b); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename CompareTypes<DualNumber<T2,D>,T,true>::supertype \
funcname (const T& a, const DualNumber<T2,D>& b) \
{ \
  return std::equivalent(a,b); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename CompareTypes<DualNumber<T,D>,T2>::supertype \
funcname (const DualNumber<T,D>& a, const T2& b) \
{ \
  return std::equivalent(a,b); \
}

#define DualNumber_equivfl_binary(funcname) \
DualNumber_equiv_binary(funcname##f, funcname) \
DualNumber_equiv_binary(funcname##l, funcname)


// if_else is necessary here to handle cases where a is negative but b
// is 0; we should have a contribution of 0 from those, not NaN.
DualNumber_std_binary(pow,
  funcval * (b.value() * a.derivatives() / a.value() +
  MetaPhysicL::if_else(b.derivatives(), b.derivatives() * std::log(a.value()), b.derivatives())))
DualNumber_std_binary(atan2,
  (b.value() * a.derivatives() - a.value() * b.derivatives()) /
  (b.value() * b.value() + a.value() * a.value()))
DualNumber_std_binary(max,
  (a.value() > b.value()) ? a.derivatives() : b.derivatives())
DualNumber_std_binary(min,
  (a.value() > b.value()) ? b.derivatives() : a.derivatives())
DualNumber_std_binary(fmod, a.derivatives())

#if __cplusplus >= 201103L
DualNumber_equivfl_binary(pow)
DualNumber_equivfl_binary(fmod)
DualNumber_std_binary(remainder, a.derivatives())
DualNumber_equivfl_binary(remainder)
DualNumber_equiv_binary(fmax, max)
DualNumber_equivfl_binary(fmax)
DualNumber_equiv_binary(fmin, min)
DualNumber_equivfl_binary(fmin)
DualNumber_std_binary(fdim, if_else(a.value() > b.value(),
                                    a.derivatives() - b.derivatives(), 0))
DualNumber_equivfl_binary(fdim)
DualNumber_std_binary(hypot, (a.value()*a.derivatives() +
                              b.value()*b.derivatives()) /
                              hypot(a.value(), b.value()))
DualNumber_equivfl_binary(hypot)
DualNumber_equivfl_binary(atan2)
#endif // __cplusplus >= 201103L

} // namespace std


#endif // METAPHYSICL_DUALNUMBER_H
