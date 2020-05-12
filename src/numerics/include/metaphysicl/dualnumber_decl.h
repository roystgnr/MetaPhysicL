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


#ifndef METAPHYSICL_DUALNUMBER_DECL_H
#define METAPHYSICL_DUALNUMBER_DECL_H

#include <ostream>
#include <limits>
#include <vector>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/dualderivatives.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/testable.h"
#include "metaphysicl/dualnumber_forward.h"

namespace MetaPhysicL {

template <typename T, typename D=T>
class NotADuckDualNumber;
template <typename T, typename D>
class DualNumberSurrogate;

/**
 * \p T denotes the type of the value
 * \p D denotes the derivative type. This may be a builtin or it could be a container
 * \p asd is short for allow_skipping_derivatives, e.g. at compile time a user may
 *    choose whether the do_derivatives static member can influence the computation
 *    of derivatives
 */
template <typename T, typename D, bool asd>
class DualNumber : public safe_bool<DualNumber<T,D,asd> >
{
public:
  typedef T value_type;

  typedef D derivatives_type;

  // allow clearer code
  static constexpr bool allow_skipping_derivatives = asd;

  DualNumber() = default;

  template <typename T2>
  DualNumber(const T2& val);

  template <typename T2, typename D2>
  DualNumber(const T2& val, const D2& deriv);

#if __cplusplus >= 201103L
  // Move constructors are useful when all your data is on the heap
  DualNumber(DualNumber<T, D, asd> && /*src*/);

  // Move assignment avoids heap operations too
  DualNumber& operator= (DualNumber<T, D, asd> && /*src*/);

  // Standard copy operations get implicitly deleted upon move
  // constructor definition, so we redefine them.
  DualNumber(const DualNumber<T, D, asd> & /*src*/);

  DualNumber& operator= (const DualNumber<T, D, asd> & /*src*/);
#endif

  template <typename T2, typename D2>
  DualNumber & operator=(const DualNumber<T2,D2,asd> & dn);

  template <typename T2, typename D2>
  DualNumber & operator=(const NotADuckDualNumber<T2,D2> & nd_dn);

  template <typename T2, typename D2>
  explicit DualNumber(const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2, typename D2>
  DualNumber & operator= (const DualNumberSurrogate<T2, D2> & dns);

  template <typename T2>
  DualNumber & operator= (const T2 & scalar);

  T& value();

  const T& value() const;

  D& derivatives();

  const D& derivatives() const;

  bool boolean_test() const;

  DualNumber<T,D,asd> operator- () const;

  DualNumber<T,D,asd> operator! () const;

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator+= (const DualNumber<T2,D2,asd>& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator+= (const NotADuckDualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T, D, asd> & operator+= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator-= (const DualNumber<T2,D2,asd>& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator-= (const NotADuckDualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T, D, asd> & operator-= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator*= (const DualNumber<T2,D2,asd>& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator*= (const NotADuckDualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T, D, asd> & operator*= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator/= (const DualNumber<T2,D2,asd>& a);

  template <typename T2, typename D2>
  DualNumber<T, D, asd> & operator/= (const NotADuckDualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T, D, asd> & operator/= (const T2& a);


  static bool do_derivatives;

private:
  T _val;
  D _deriv;
};

// Helper class to handle partial specialization for DualNumber
// constructors

template <typename T, typename D, bool asd = false>
struct DualNumberConstructor
{
  static T value(const DualNumber<T,D,asd>& v) { return v.value(); }

  template <typename T2>
  static T value(const T2& v) { return v; }

  template <typename T2, typename D2>
  static T value(const T2& v, const D2&) { return v; }

  template <typename T2, typename D2>
  static T value(const DualNumber<T2,D2,asd>& v) {
    return DualNumberConstructor<T,D,asd>::value(v.value());
  }

  template <typename T2>
  static D deriv(const T2&) { return 0.; }

  template <typename T2, typename D2>
  static D deriv(const DualNumber<T2,D2,asd>& v) { return v.derivatives(); }

  template <typename T2, typename D2>
  static D deriv(const T2&, const D2& d) { return d; }
};

template <typename T, typename D, bool asd, typename DD>
struct DualNumberConstructor<DualNumber<T,D,asd>, DD, asd>
{
  template <typename T2, typename D2, typename D3>
  static DualNumber<T,D,asd> value(const DualNumber<DualNumber<T2,D2,asd>, D3,asd>& v) { return v.value(); }

  template <typename T2>
  static DualNumber<T,D,asd> value(const T2& v) { return v; }

  template <typename T2, typename D2>
  static DualNumber<T,D,asd> value(const T2& v, const D2& d) { return DualNumber<T,D,asd>(v,d); }

  template <typename D2>
  static DualNumber<T,D,asd> value(const DualNumber<T,D,asd>& v, const D2&) { return v; }

  template <typename T2>
  static DD deriv(const T2&) { return 0; }

  template <typename T2, typename D2>
  static DD deriv(const DualNumber<T2,D2,asd>& v) { return v.derivatives(); }

  template <typename T2, typename D2>
  static DD deriv(const T2&, const D2& d) { return d; }
};

// FIXME: these operators currently do automatic type promotion when
// encountering DualNumbers of differing levels of recursion and
// differentiability.  But what we really want is automatic type
// *demotion*, to avoid pretending we have accurate derivatives which
// we don't have.  If we could do that right then some potential
// subtle run-time user errors would turn into compile-time user
// errors.

#define DualNumber_decl_preop(opname, functorname) \
template <typename T, typename D, typename T2, typename D2, bool asd>  \
inline \
typename functorname##Type<DualNumber<T,D,asd>,DualNumber<T2,D2,asd> >::supertype \
operator opname (const DualNumber<T,D,asd>& a, const DualNumber<T2,D2,asd>& b); \
 \
 \
template <typename T, typename T2, typename D, bool asd> \
inline \
typename functorname##Type<DualNumber<T2,D,asd>,T,true>::supertype \
operator opname (const T& a, const DualNumber<T2,D,asd>& b); \
 \
 \
template <typename T, typename D, typename T2, bool asd> \
inline \
typename functorname##Type<DualNumber<T,D,asd>,T2,false>::supertype \
operator opname (const DualNumber<T,D,asd>& a, const T2& b);



// With C++11, define "move operations" where possible.  We should be
// more complete and define the move-from-b alternatives as well, but
// those would require additional support to correctly handle
// division, subtraction, or non-commutative addition/multiplication
#if __cplusplus >= 201103L
#define DualNumber_decl_op(opname, functorname) \
        DualNumber_decl_preop(opname, functorname) \
 \
template <typename T, typename D, typename T2, typename D2, bool asd> \
inline \
typename functorname##Type<DualNumber<T,D,asd>,DualNumber<T2,D2,asd> >::supertype \
operator opname (DualNumber<T,D,asd>&& a, const DualNumber<T2,D2,asd>& b); \
 \
 \
template <typename T, typename D, typename T2, bool asd> \
inline \
typename functorname##Type<DualNumber<T,D,asd>,T2,false>::supertype \
operator opname (DualNumber<T,D,asd>&& a, const T2& b); \

#else
#define DualNumber_decl_op(opname, functorname) \
        DualNumber_decl_preop(opname, functorname)
#endif

DualNumber_decl_op(+, Plus)
DualNumber_decl_op(-, Minus)
DualNumber_decl_op(*, Multiplies)
DualNumber_decl_op(/, Divides)

#define DualNumber_decl_compare(opname)                     \
template <typename T, typename D, typename T2, typename D2, bool asd> \
inline \
bool \
operator opname  (const DualNumber<T,D,asd>& a, const DualNumber<T2,D2,asd>& b); \
 \
 \
template <typename T, typename T2, typename D2, bool asd> \
inline \
typename boostcopy::enable_if_class< \
  typename CompareTypes<DualNumber<T2,D2,asd>,T>::supertype, \
  bool \
>::type \
operator opname  (const T& a, const DualNumber<T2,D2,asd>& b); \
 \
 \
template <typename T, typename T2, typename D, bool asd> \
inline \
typename boostcopy::enable_if_class< \
  typename CompareTypes<DualNumber<T,D,asd>,T2>::supertype, \
  bool \
>::type \
operator opname  (const DualNumber<T,D,asd>& a, const T2& b);

DualNumber_decl_compare(>)
DualNumber_decl_compare(>=)
DualNumber_decl_compare(<)
DualNumber_decl_compare(<=)
DualNumber_decl_compare(==)
DualNumber_decl_compare(!=)
DualNumber_decl_compare(&&)
DualNumber_decl_compare(||)

template <typename T, typename D, bool asd>
inline
std::ostream&
operator<< (std::ostream& output, const DualNumber<T,D,asd>& a);


// ScalarTraits, RawType, CompareTypes specializations

template <typename T, typename D, bool asd>
struct ScalarTraits<DualNumber<T, D, asd> >
{
  static const bool value = ScalarTraits<T>::value;
};

template <typename T, typename D, bool asd>
struct RawType<DualNumber<T, D, asd> >
{
  typedef typename RawType<T>::value_type value_type;

  static value_type value(const DualNumber<T, D, asd>& a) { return raw_value(a.value()); }
};

// vector specialization
template <typename T, typename D, bool asd>
struct RawType<std::vector<DualNumber<T,D,asd>>>
{
  typedef std::vector<typename RawType<T>::value_type> value_type;

  static value_type value(const std::vector<DualNumber<T,D,asd>> & in)
  {
    value_type ret_val(in.size());
    for (std::size_t i = 0; i < in.size(); ++i)
      ret_val[i] = raw_value(in[i]);
  }
};

// vector of vectors specialization
template <typename T, typename D, bool asd>
struct RawType<std::vector<std::vector<DualNumber<T,D,asd>>>>
{
  typedef std::vector<std::vector<typename RawType<T>::value_type>> value_type;

  static value_type value(const std::vector<std::vector<DualNumber<T,D,asd>>> & in)
  {
    value_type ret_val(in.size());
    for (std::size_t i = 0; i < in.size(); ++i)
    {
      ret_val[i].resize(in[i].size());
      for (std::size_t j = 0; j < in[i].size(); ++j)
        ret_val[i][j] = raw_value(in[i][j]);
    }
  }
};

template<typename T, typename T2, typename D, bool asd, bool reverseorder>
struct PlusType<DualNumber<T, D, asd>, T2, reverseorder,
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricPlusType<T, T2, reverseorder>::supertype, D, asd> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd, bool reverseorder>
struct PlusType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricPlusType<T, T2, reverseorder>::supertype,
                     typename SymmetricPlusType<D, D2, reverseorder>::supertype,
                     asd> supertype;
};

template<typename T, typename D, bool asd>
struct PlusType<DualNumber<T, D, asd>, DualNumber<T, D, asd> > {
  typedef DualNumber<typename SymmetricPlusType<T,T>::supertype,
                     typename SymmetricPlusType<D,D>::supertype,
                     asd> supertype;
};


template<typename T, typename T2, typename D, bool asd, bool reverseorder>
struct MinusType<DualNumber<T, D, asd>, T2, reverseorder,
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricMinusType<T, T2, reverseorder>::supertype, D, asd> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd, bool reverseorder>
struct MinusType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricMinusType<T, T2, reverseorder>::supertype,
                     typename SymmetricMinusType<D, D2, reverseorder>::supertype,
                     asd> supertype;
};

template<typename T, typename D, bool asd, bool reverseorder>
struct MinusType<DualNumber<T, D, asd>, DualNumber<T, D, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricMinusType<T,T>::supertype,
                     typename SymmetricMinusType<D,D>::supertype,
                     asd> supertype;
};


template<typename T, typename T2, typename D, bool asd, bool reverseorder>
struct MultipliesType<DualNumber<T, D, asd>, T2, reverseorder,
                      typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricMultipliesType<T, T2, reverseorder>::supertype,
                     typename SymmetricMultipliesType<D, T2, reverseorder>::supertype,
                     asd> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd, bool reverseorder>
struct MultipliesType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricMultipliesType<T, T2, reverseorder>::supertype,
                     typename SymmetricPlusType<
                       typename SymmetricMultipliesType<T, D2, reverseorder>::supertype,
                       typename SymmetricMultipliesType<D, T2, reverseorder>::supertype>::supertype,
                     asd> supertype;
};

template<typename T, typename D, bool asd, bool reverseorder>
struct MultipliesType<DualNumber<T, D, asd>, DualNumber<T, D, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricMultipliesType<T, T, reverseorder>::supertype,
                     typename SymmetricMultipliesType<T, D, reverseorder>::supertype,
                     asd> supertype;
};


template<typename T, typename T2, typename D, bool asd>
struct DividesType<DualNumber<T, D, asd>, T2, false,
                      typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricDividesType<T, T2>::supertype,
                     typename SymmetricDividesType<D, T2>::supertype,
                     asd> supertype;
};

template<typename T, typename D, typename T2, bool asd>
struct DividesType<DualNumber<T, D, asd>, T2, true,
                   typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricDividesType<T2, T>::supertype,
                     typename SymmetricDividesType<
                       typename SymmetricMultipliesType<T2, D>::supertype,
                       T
                     >::supertype,
                    asd> supertype;
};


template<typename T, typename D, typename T2, typename D2, bool asd>
struct DividesType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, false> {
  typedef DualNumber<typename SymmetricDividesType<T, T2>::supertype,
                     typename SymmetricMinusType<
                       typename SymmetricDividesType<T2, D>::supertype,
                       typename SymmetricDividesType<
                         typename SymmetricMultipliesType<T, D2>::supertype,
                         T2
                       >::supertype
                     >::supertype,
                    asd> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd>
struct DividesType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, true> {
  typedef typename DividesType<DualNumber<T2, D2, asd>, DualNumber<T, D, asd>, false>::supertype supertype;
};

template<typename T, typename D, bool asd>
struct DividesType<DualNumber<T, D, asd>, DualNumber<T, D, asd>, false> {
  typedef DualNumber<T,
                     typename SymmetricMinusType<
                       typename SymmetricDividesType<T, D>::supertype,
                       typename SymmetricDividesType<
                         typename SymmetricMultipliesType<T, D>::supertype,
                         T
                       >::supertype
                     >::supertype,
                    asd> supertype;
};

template<typename T, typename D, bool asd>
struct DividesType<DualNumber<T, D, asd>, DualNumber<T, D, asd>, true> {
  typedef typename DividesType<DualNumber<T, D, asd>, DualNumber<T, D, asd>, false>::supertype supertype;
};

template<typename T, typename T2, typename D, bool asd, bool reverseorder>
struct AndType<DualNumber<T, D, asd>, T2, reverseorder,
               typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricAndType<T, T2, reverseorder>::supertype, bool> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd, bool reverseorder>
struct AndType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricAndType<T, T2, reverseorder>::supertype,
                     bool> supertype;
};

template<typename T, typename D, bool asd>
struct AndType<DualNumber<T, D, asd>, DualNumber<T, D, asd> > {
  typedef DualNumber<typename SymmetricAndType<T,T>::supertype,
                     bool> supertype;
};

template<typename T, typename T2, typename D, bool asd, bool reverseorder>
struct OrType<DualNumber<T, D, asd>, T2, reverseorder,
              typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricOrType<T, T2, reverseorder>::supertype, bool> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd, bool reverseorder>
struct OrType<DualNumber<T, D, asd>, DualNumber<T2, D2, asd>, reverseorder> {
  typedef DualNumber<typename SymmetricOrType<T, T2, reverseorder>::supertype,
                     bool> supertype;
};

template<typename T, typename D, bool asd>
struct OrType<DualNumber<T, D, asd>, DualNumber<T, D, asd> > {
  typedef DualNumber<typename SymmetricOrType<T,T>::supertype,
                     bool> supertype;
};




template<typename T, typename T2, typename D, bool asd, bool reverseorder>
struct CompareTypes<DualNumber<T, D, asd>, T2, reverseorder,
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> {
  typedef DualNumber<typename SymmetricCompareTypes<T, T2>::supertype,
                     typename SymmetricCompareTypes<
                       typename SymmetricCompareTypes<D, T2>::supertype,
                       T
                     >::supertype, asd> supertype;
};

template<typename T, typename D, typename T2, typename D2, bool asd>
struct CompareTypes<DualNumber<T, D, asd>, DualNumber<T2, D2, asd> > {
  typedef DualNumber<typename SymmetricCompareTypes<T, T2>::supertype,
                     typename SymmetricCompareTypes<
                       typename SymmetricCompareTypes<T, T2>::supertype,
                       typename SymmetricCompareTypes<D, D2>::supertype
                     >::supertype,
                    asd> supertype;
};

template<typename T, typename D, bool asd>
struct CompareTypes<DualNumber<T, D, asd>, DualNumber<T, D, asd> > {
  typedef DualNumber<T, typename SymmetricCompareTypes<T, D>::supertype, asd> supertype;
};


template <typename T, typename D, bool asd>
inline
D gradient(const DualNumber<T, D, asd>& a);

} // namespace MetaPhysicL


namespace std {

using MetaPhysicL::DualNumber;
using MetaPhysicL::CompareTypes;

template <typename T, typename D, bool asd>
inline bool isnan (const DualNumber<T,D,asd> & a);

// Some forward declarations necessary for recursive DualNumbers

#if __cplusplus >= 201103L

template <typename T, typename D, bool asd>
inline DualNumber<T,D,asd> cos  (const DualNumber<T,D,asd> & a);

template <typename T, typename D, bool asd>
inline DualNumber<T,D,asd> cos  (DualNumber<T,D,asd> && a);

template <typename T, typename D, bool asd>
inline DualNumber<T,D,asd> cosh (const DualNumber<T,D,asd> & a);

template <typename T, typename D, bool asd>
inline DualNumber<T,D,asd> cosh (DualNumber<T,D,asd> && a);

#else

template <typename T, typename D, bool asd>
inline DualNumber<T,D,asd> cos  (DualNumber<T,D,asd> a);

template <typename T, typename D, bool asd>
inline DualNumber<T,D,asd> cosh (DualNumber<T,D,asd> a);

#endif

// Now just combined declaration/definitions

#if __cplusplus >= 201103L
#define DualNumber_decl_std_unary(funcname) \
template <typename T, typename D, bool asd> \
inline \
DualNumber<T,D,asd> funcname (const DualNumber<T,D,asd> & in); \
 \
 \
template <typename T, typename D, bool asd> \
inline \
DualNumber<T,D,asd> funcname (DualNumber<T,D,asd> && in);


#else

#define DualNumber_decl_std_unary(funcname) \
template <typename T, typename D, bool asd> \
inline \
DualNumber<T,D,asd> funcname (DualNumber<T,D,asd> in);

#endif

#define DualNumber_decl_fl_unary(funcname) \
DualNumber_decl_std_unary(funcname##f) \
DualNumber_decl_std_unary(funcname##l)

#define DualNumber_decl_stdfl_unary(funcname) \
DualNumber_decl_std_unary(funcname) \
DualNumber_decl_fl_unary(funcname)

DualNumber_decl_std_unary(sqrt)
DualNumber_decl_std_unary(exp)
DualNumber_decl_std_unary(log)
DualNumber_decl_std_unary(log10)
DualNumber_decl_std_unary(sin)
DualNumber_decl_std_unary(cos)
DualNumber_decl_std_unary(tan)
DualNumber_decl_std_unary(asin)
DualNumber_decl_std_unary(acos)
DualNumber_decl_std_unary(atan)
DualNumber_decl_std_unary(sinh)
DualNumber_decl_std_unary(cosh)
DualNumber_decl_std_unary(tanh)
DualNumber_decl_std_unary(abs)
DualNumber_decl_std_unary(norm)
DualNumber_decl_std_unary(fabs)
DualNumber_decl_std_unary(ceil)
DualNumber_decl_std_unary(floor)

#if __cplusplus >= 201103L
DualNumber_decl_std_unary(llabs)
DualNumber_decl_std_unary(imaxabs)
DualNumber_decl_fl_unary(fabs)
DualNumber_decl_fl_unary(exp)
DualNumber_decl_stdfl_unary(exp2)
DualNumber_decl_stdfl_unary(expm1)
DualNumber_decl_fl_unary(log)
DualNumber_decl_fl_unary(log10)
DualNumber_decl_stdfl_unary(log2)
DualNumber_decl_stdfl_unary(log1p)
DualNumber_decl_fl_unary(sqrt)
DualNumber_decl_stdfl_unary(cbrt)
DualNumber_decl_fl_unary(sin)
DualNumber_decl_fl_unary(cos)
DualNumber_decl_fl_unary(tan)
DualNumber_decl_fl_unary(asin)
DualNumber_decl_fl_unary(acos)
DualNumber_decl_fl_unary(atan)
DualNumber_decl_fl_unary(sinh)
DualNumber_decl_fl_unary(cosh)
DualNumber_decl_fl_unary(tanh)
DualNumber_decl_stdfl_unary(asinh)
DualNumber_decl_stdfl_unary(acosh)
DualNumber_decl_stdfl_unary(atanh)
DualNumber_decl_stdfl_unary(erf)
DualNumber_decl_stdfl_unary(erfc)
DualNumber_decl_fl_unary(ceil)
DualNumber_decl_fl_unary(floor)
DualNumber_decl_stdfl_unary(trunc)
DualNumber_decl_stdfl_unary(round)
DualNumber_decl_stdfl_unary(nearbyint)
DualNumber_decl_stdfl_unary(rint)
#endif // __cplusplus >= 201103L

#define DualNumber_decl_complex_std_unary_real(funcname) \
template <typename T, typename D, bool asd> \
inline DualNumber<T, typename D::template rebind<T>::other, asd> \
funcname(const DualNumber<std::complex<T>, D, asd> & in); \
template <typename T, bool asd> \
inline DualNumber<T,T,asd> \
funcname(const DualNumber<std::complex<T>,std::complex<T>,asd> & in)

DualNumber_decl_complex_std_unary_real(real);
DualNumber_decl_complex_std_unary_real(imag);
DualNumber_decl_complex_std_unary_real(norm);
DualNumber_decl_complex_std_unary_real(abs);

#define DualNumber_decl_complex_std_unary_complex_pre(funcname) \
template <typename T, typename D, bool asd> \
inline DualNumber<std::complex<T>, D, asd> \
funcname(const DualNumber<std::complex<T>, D, asd> & in); \
template <typename T, bool asd> \
inline DualNumber<std::complex<T>,std::complex<T>,asd> \
funcname(const DualNumber<std::complex<T>,std::complex<T>,asd> & in)

#if __cplusplus >= 201103L
#define DualNumber_decl_complex_std_unary_complex(funcname) \
DualNumber_decl_complex_std_unary_complex_pre(funcname);  \
template <typename T, typename D, bool asd> \
inline DualNumber<std::complex<T>, D, asd> \
funcname(DualNumber<std::complex<T>, D, asd> && in); \
 \
template <typename T, bool asd> \
inline DualNumber<std::complex<T>,std::complex<T>,asd> \
funcname(DualNumber<std::complex<T>,std::complex<T>,asd> && in)

#else
#define DualNumber_decl_complex_std_unary_complex(funcname) \
DualNumber_complex_std_unary_complex_pre(funcname);
#endif

DualNumber_decl_complex_std_unary_complex(conj);

#define DualNumber_decl_std_binary(funcname)                \
template <typename T, typename D, typename T2, typename D2, bool asd> \
inline \
typename CompareTypes<DualNumber<T,D,asd>,DualNumber<T2,D2,asd> >::supertype \
funcname (const DualNumber<T,D,asd>& a, const DualNumber<T2,D2,asd>& b); \
 \
 \
template <typename T, typename D, bool asd> \
inline \
DualNumber<T,D,asd> \
funcname (const DualNumber<T,D,asd>& a, const DualNumber<T,D,asd>& b); \
 \
 \
template <typename T, typename T2, typename D, bool asd> \
inline \
typename CompareTypes<DualNumber<T2,D,asd>,T,true>::supertype \
funcname (const T& a, const DualNumber<T2,D,asd>& b); \
 \
 \
template <typename T, typename T2, typename D, bool asd> \
inline \
typename CompareTypes<DualNumber<T,D,asd>,T2>::supertype \
funcname (const DualNumber<T,D,asd>& a, const T2& b);

#define DualNumber_decl_fl_binary(funcname) \
DualNumber_decl_std_binary(funcname##f) \
DualNumber_decl_std_binary(funcname##l)

#define DualNumber_decl_stdfl_binary(funcname) \
DualNumber_decl_std_binary(funcname) \
DualNumber_decl_fl_binary(funcname)

DualNumber_decl_std_binary(pow)
DualNumber_decl_std_binary(atan2)
DualNumber_decl_std_binary(max)
DualNumber_decl_std_binary(min)
DualNumber_decl_std_binary(fmod)

#if __cplusplus >= 201103L
DualNumber_decl_fl_binary(pow)
DualNumber_decl_fl_binary(fmod)
DualNumber_decl_stdfl_binary(remainder)
DualNumber_decl_stdfl_binary(fmax)
DualNumber_decl_stdfl_binary(fmin)
DualNumber_decl_stdfl_binary(fdim)
DualNumber_decl_stdfl_binary(hypot)
DualNumber_decl_fl_binary(atan2)
#endif // __cplusplus >= 201103L

template <typename T, typename D, bool asd>
class numeric_limits<DualNumber<T, D, asd> > :
  public MetaPhysicL::raw_numeric_limits<DualNumber<T, D, asd>, T> {};

} // namespace std


#endif // METAPHYSICL_DUALNUMBER_DECL_H
