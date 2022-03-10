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
// $Id: core.h 37197 2013-02-21 05:49:09Z roystgnr $
//
//--------------------------------------------------------------------------


#ifndef METAPHYSICL_DYNAMICSPARSENUMBERBASE_DECL_H
#define METAPHYSICL_DYNAMICSPARSENUMBERBASE_DECL_H

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <ostream>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/ct_set.h"
#include "metaphysicl/metaphysicl_asserts.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/sparsenumberutils.h"
#include "metaphysicl/testable.h"

namespace MetaPhysicL {

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
class DynamicSparseNumberBase
{
public:
  typedef typename Data::value_type T;
  typedef typename Indices::value_type I;

  typedef T value_type;

  template <unsigned int i>
  struct entry_type {
    typedef value_type type;
  };

  typedef I index_value_type;

  std::size_t size() const;

  void resize(std::size_t s);

  DynamicSparseNumberBase();

#if __cplusplus >= 201103L
  // Move constructors are useful when all your data is on the heap
  DynamicSparseNumberBase(DynamicSparseNumberBase && src) = default;

  // Move assignment avoids heap operations too
  DynamicSparseNumberBase& operator= (DynamicSparseNumberBase && src) = default;

  // Standard copy operations get implicitly deleted upon move
  // constructor definition, so we manually enable them.
  DynamicSparseNumberBase(const DynamicSparseNumberBase & src) = default;

  DynamicSparseNumberBase& operator= (const DynamicSparseNumberBase & src) = default;
#endif

  template <typename Data2, typename Indices2, class... SubTypeArgs2>
  DynamicSparseNumberBase(
      const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & src);

  template <typename Data2, typename Indices2, class... SubTypeArgs2>
  DynamicSparseNumberBase(
      DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> && src);

  T* raw_data();

  const T* raw_data() const;

  typename Data::reference raw_at(unsigned int i);

  typename Data::const_reference raw_at(unsigned int i) const;

  I& raw_index(unsigned int i);

  const I& raw_index(unsigned int i) const;

  // FIXME: these encapsulation violations are necessary for std::pow
  // until I can figure out the right friend declaration.
  const Data& nude_data() const;

  Data& nude_data();

  const Indices& nude_indices() const;

  Indices& nude_indices();

  std::size_t runtime_index_query(index_value_type i) const;

  std::size_t runtime_index_of(index_value_type i) const;

  T& operator[](index_value_type i);

  const T& operator[](index_value_type i) const;

  template <unsigned int i>
  typename entry_type<i>::type& get();

  template <unsigned int i>
  const typename entry_type<i>::type& get() const;

  value_type& insert(unsigned int i);

  template <unsigned int i>
  typename entry_type<i>::type& insert();

  template <unsigned int i, typename T2>
  void set(const T2& val);

  bool boolean_test() const;

  SubType<SubTypeArgs...> operator- () const;

  // Since this is a dynamically allocated sparsity pattern, we can
  // increase it as needed to support e.g. operator+=
  template <typename Indices2>
  void sparsity_union (const Indices2& new_indices);

  // Since this is a dynamically allocated sparsity pattern, we can
  // decrease it when possible for efficiency
  template <typename Indices2>
  void sparsity_intersection (const Indices2& new_indices);

  // Since this is a dynamically allocated sparsity pattern, we can
  // decrease it when possible for efficiency. This method will remove
  // any index-data pairs for whom the data entry is less than or equal
  // to the prescribed tolerance
  void sparsity_trim (value_type tolerance = 0);

  // Not defineable since !0 != 0
  // SubType<SubTypeArgs...> operator! () const;

  template <class... SubTypeArgs2>
  SubType<SubTypeArgs...> & operator+=(const SubType<SubTypeArgs2...> & a);

  template <class... SubTypeArgs2>
  SubType<SubTypeArgs...> & operator-=(const SubType<SubTypeArgs2...> & a);

  template <class... SubTypeArgs2>
  SubType<SubTypeArgs...> & operator*=(const SubType<SubTypeArgs2...> & a);

  template <class... SubTypeArgs2>
  SubType<SubTypeArgs...> & operator/=(const SubType<SubTypeArgs2...> & a);

  template <typename T2>
  SubType<SubTypeArgs...> & operator*=(const T2 & a);

  template <typename T2>
  SubType<SubTypeArgs...> & operator/=(const T2 & a);

protected:

  Data _data;
  Indices _indices;
};


//
// Non-member functions
//

template <template <class...> class SubType,
          typename BoolData,
          typename BoolIndices,
          class... BoolSubTypeArgs,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename Data2,
          typename Indices2,
          class... SubTypeArgs2>
inline typename SubType<SubTypeArgs...>::template rebind<
    typename CompareTypes<typename Data::value_type, typename Data2::value_type>::supertype,
    typename CompareTypes<typename Indices::value_type,
                          typename Indices2::value_type>::supertype>::other
if_else(const DynamicSparseNumberBase<BoolData, BoolIndices, SubType, BoolSubTypeArgs...> &
            condition,
        const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & if_true,
        const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & if_false);



#define DynamicSparseNumberBase_decl_op_ab(opname, subtypename, functorname) \
template <class... AArgs, class... BArgs> \
inline typename Symmetric##functorname##Type<subtypename<AArgs...>, \
                                             subtypename<BArgs...>>::supertype \
operator opname(const subtypename<AArgs...> & a, const subtypename<BArgs...> & b);


#if __cplusplus >= 201103L

#define DynamicSparseNumberBase_decl_op(subtypename, opname, functorname) \
DynamicSparseNumberBase_decl_op_ab(opname, subtypename, functorname) \
 \
template <class... AArgs, class... BArgs> \
inline typename Symmetric##functorname##Type<subtypename<AArgs...>, \
                                             subtypename<BArgs...>>::supertype \
operator opname(subtypename<AArgs...> && a, const subtypename<BArgs...> & b);

#else

#define DynamicSparseNumberBase_decl_op(subtypename, opname, functorname) \
DynamicSparseNumberBase_decl_op_ab(opname, subtypename, functorname)

#endif

// Let's also allow scalar times vector.
// Scalar plus vector, etc. remain undefined in the sparse context.

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename MultipliesType<SubType<SubTypeArgs...>, T, true>::supertype
operator*(const T & a,
          const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & b);

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename MultipliesType<SubType<SubTypeArgs...>, T>::supertype
operator*(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a,
          const T & b);

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename DividesType<SubType<SubTypeArgs...>, T>::supertype
operator/(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a,
          const T & b);

#if __cplusplus >= 201103L

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename MultipliesType<SubType<SubTypeArgs...>, T>::supertype
operator*(DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> && a, const T & b);

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename DividesType<SubType<SubTypeArgs...>, T>::supertype
operator/(DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> && a, const T & b);

#endif

#define DynamicSparseNumberBase_decl_operator_binary(opname, functorname) \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename Data2, \
          typename Indices2, \
          class... SubTypeArgs2> \
inline typename SubType<SubTypeArgs...>::template rebind< \
    bool, \
    typename CompareTypes<typename Indices::value_type, \
                          typename Indices2::value_type>::supertype>::other \
operator opname( \
    const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, \
    const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & b); \
 \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename T> \
inline typename SubType<SubTypeArgs...>::template rebind<bool>::other operator opname( \
    const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, const T & b); \
 \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename T> \
inline typename SubType<SubTypeArgs...>::template rebind<bool>::other operator opname( \
    const T & a, const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & b);

// NOTE: unary functions for which 0-op-0 is true are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

DynamicSparseNumberBase_decl_operator_binary(<, less)
// DynamicSparseNumberBase_decl_operator_binary(<=)
DynamicSparseNumberBase_decl_operator_binary(>, greater)
// DynamicSparseNumberBase_decl_operator_binary(>=)
// DynamicSparseNumberBase_decl_operator_binary(==)
DynamicSparseNumberBase_decl_operator_binary(!=, not_equal_to)

// FIXME - make && an intersection rather than a union for efficiency
DynamicSparseNumberBase_decl_operator_binary(&&, logical_and)
DynamicSparseNumberBase_decl_operator_binary(||, logical_or)


template <template <typename, typename> class SubType,
          typename T, typename I>
inline
std::ostream&
operator<< (std::ostream& output, const DynamicSparseNumberBase<T,I,SubType>& a);


// CompareTypes, RawType, ValueType specializations

#define DynamicSparseNumberBase_comparisons(subtypename, templatename) \
template <class... SubTypeArgs, bool reverseorder> \
struct templatename<subtypename<SubTypeArgs...>, subtypename<SubTypeArgs...>, reverseorder> \
{ \
  typedef subtypename<SubTypeArgs...> supertype; \
}; \
 \
template <class... SubTypeArgs, class... SubTypeArgs2, bool reverseorder> \
struct templatename<subtypename<SubTypeArgs...>, subtypename<SubTypeArgs2...>, reverseorder> \
{ \
  typedef typename subtypename<SubTypeArgs...>::template rebind< \
      typename Symmetric##templatename<typename subtypename<SubTypeArgs...>::T, \
                                       typename subtypename<SubTypeArgs2...>::T, \
                                       reverseorder>::supertype, \
      typename CompareTypes<typename subtypename<SubTypeArgs...>::I, \
                            typename subtypename<SubTypeArgs2...>::I>::supertype>::other \
      supertype; \
}; \
 \
template <typename T, class... SubTypeArgs, bool reverseorder> \
struct templatename<subtypename<SubTypeArgs...>, \
                    T, \
                    reverseorder, \
                    typename boostcopy::enable_if<BuiltinTraits<T>>::type> \
{ \
  typedef typename subtypename<SubTypeArgs...>::template rebind< \
      typename Symmetric##templatename<typename subtypename<SubTypeArgs...>::T, \
                                       T, \
                                       reverseorder>::supertype, \
      typename subtypename<SubTypeArgs...>::I>::other supertype; \
}

} // namespace MetaPhysicL

namespace std {

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DynamicSparseNumberBase;
using MetaPhysicL::SymmetricCompareTypes;

#define DynamicSparseNumberBase_decl_std_unary(funcname) \
template <template <class...> class SubType, \
          typename Data, typename Indices, class... SubTypeArgs>  \
inline \
SubType<SubTypeArgs...> \
funcname (const DynamicSparseNumberBase<Data,Indices,SubType,SubTypeArgs...> & a);


#define DynamicSparseNumberBase_decl_std_binary_union(funcname) \
template <template <class...> class SubType, \
          typename Data, typename Indices, class... SubTypeArgs, \
          typename Data2, typename Indices2, class... SubTypeArgs2> \
inline \
typename SubType<SubTypeArgs...>::template rebind< \
      typename SymmetricCompareTypes<typename SubType<SubTypeArgs...>::T, \
                                     typename SubType<SubTypeArgs2...>::T>::supertype, \
      typename CompareTypes<typename SubType<SubTypeArgs...>::I, \
                            typename SubType<SubTypeArgs2...>::I>::supertype>::other \
funcname (const DynamicSparseNumberBase<Data,Indices,SubType,SubTypeArgs...>& a, \
          const DynamicSparseNumberBase<Data2,Indices2,SubType,SubTypeArgs2...>& b); \
 \
 \
template <template <class...> class SubType, \
          typename Data, typename Indices, class... SubTypeArgs, typename T> \
inline \
typename SubType<SubTypeArgs...>::template rebind< \
      typename SymmetricCompareTypes<typename SubType<SubTypeArgs...>::T, T>::supertype, \
      typename SubType<SubTypeArgs...>::I>::other \
funcname (const DynamicSparseNumberBase<Data,Indices,SubType,SubTypeArgs...>& a, const T& b); \
 \
 \
template <template <class...> class SubType, \
          typename Data, typename Indices, class... SubTypeArgs, typename T> \
inline \
typename SubType<SubTypeArgs...>::template rebind< \
      typename SymmetricCompareTypes<T, typename SubType<SubTypeArgs...>::T>::supertype, \
      typename SubType<SubTypeArgs...>::I>::other \
funcname (const T& a, const DynamicSparseNumberBase<Data,Indices,SubType,SubTypeArgs...>& b);


#define DynamicSparseNumberBase_decl_fl_unary(funcname) \
DynamicSparseNumberBase_decl_std_unary(funcname##f) \
DynamicSparseNumberBase_decl_std_unary(funcname##l)


#define DynamicSparseNumberBase_decl_stdfl_unary(funcname) \
DynamicSparseNumberBase_decl_std_unary(funcname) \
DynamicSparseNumberBase_decl_fl_unary(funcname)


#define DynamicSparseNumberBase_decl_fl_binary_union(funcname) \
DynamicSparseNumberBase_decl_std_binary_union(funcname##f) \
DynamicSparseNumberBase_decl_std_binary_union(funcname##l)

#define DynamicSparseNumberBase_decl_stdfl_binary_union(funcname) \
DynamicSparseNumberBase_decl_std_binary_union(funcname) \
DynamicSparseNumberBase_decl_fl_binary_union(funcname)


// Pow needs its own specialization, both to avoid being confused by
// pow<T1,T2> and because pow(x,0) isn't 0.
template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T2>
inline typename SubType<SubTypeArgs...>::template rebind<
    typename SymmetricCompareTypes<typename Data::value_type, T2>::supertype,
    typename Indices::value_type>::other
pow(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, const T2 & b);

// NOTE: unary functions for which f(0) != 0 are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

// DynamicSparseNumberBase_decl_std_binary(pow) // separate definition
// DynamicSparseNumberBase_decl_std_unary(exp)
// DynamicSparseNumberBase_decl_std_unary(log)
// DynamicSparseNumberBase_decl_std_unary(log10)
DynamicSparseNumberBase_decl_std_unary(sin)
// DynamicSparseNumberBase_decl_std_unary(cos)
DynamicSparseNumberBase_decl_std_unary(tan)
DynamicSparseNumberBase_decl_std_unary(asin)
// DynamicSparseNumberBase_decl_std_unary(acos)
DynamicSparseNumberBase_decl_std_unary(atan)
DynamicSparseNumberBase_decl_std_binary_union(atan2)
DynamicSparseNumberBase_decl_std_unary(sinh)
// DynamicSparseNumberBase_decl_std_unary(cosh)
DynamicSparseNumberBase_decl_std_unary(tanh)
DynamicSparseNumberBase_decl_std_unary(sqrt)
DynamicSparseNumberBase_decl_std_unary(abs)
DynamicSparseNumberBase_decl_std_unary(fabs)
DynamicSparseNumberBase_decl_std_binary_union(max)
DynamicSparseNumberBase_decl_std_binary_union(min)
DynamicSparseNumberBase_decl_std_unary(ceil)
DynamicSparseNumberBase_decl_std_unary(floor)
DynamicSparseNumberBase_decl_std_binary_union(fmod) // TODO: optimize this

#if __cplusplus >= 201103L
DynamicSparseNumberBase_decl_std_unary(llabs)
DynamicSparseNumberBase_decl_std_unary(imaxabs)
DynamicSparseNumberBase_decl_fl_unary(fabs)
// DynamicSparseNumberBase_decl_fl_unary(exp)
// DynamicSparseNumberBase_decl_stdfl_unary(exp2)
DynamicSparseNumberBase_decl_stdfl_unary(expm1)
// DynamicSparseNumberBase_decl_fl_unary(log)
// DynamicSparseNumberBase_decl_fl_unary(log10)
// DynamicSparseNumberBase_decl_stdfl_unary(log2)
// DynamicSparseNumberBase_decl_stdfl_unary(log1p)
DynamicSparseNumberBase_decl_fl_unary(sqrt)
DynamicSparseNumberBase_decl_stdfl_unary(cbrt)
DynamicSparseNumberBase_decl_fl_unary(sin)
// DynamicSparseNumberBase_decl_fl_unary(cos)
DynamicSparseNumberBase_decl_fl_unary(tan)
DynamicSparseNumberBase_decl_fl_unary(asin)
// DynamicSparseNumberBase_decl_fl_unary(acos)
DynamicSparseNumberBase_decl_fl_unary(atan)
DynamicSparseNumberBase_decl_fl_unary(sinh)
// DynamicSparseNumberBase_decl_fl_unary(cosh)
DynamicSparseNumberBase_decl_fl_unary(tanh)
DynamicSparseNumberBase_decl_stdfl_unary(asinh)
// DynamicSparseNumberBase_decl_stdfl_unary(acosh)
DynamicSparseNumberBase_decl_stdfl_unary(atanh)
DynamicSparseNumberBase_decl_stdfl_unary(erf)
// DynamicSparseNumberBase_decl_stdfl_unary(erfc)
DynamicSparseNumberBase_decl_fl_unary(ceil)
DynamicSparseNumberBase_decl_fl_unary(floor)
DynamicSparseNumberBase_decl_stdfl_unary(trunc)
DynamicSparseNumberBase_decl_stdfl_unary(round)
DynamicSparseNumberBase_decl_stdfl_unary(nearbyint)
DynamicSparseNumberBase_decl_stdfl_unary(rint)

DynamicSparseNumberBase_decl_fl_binary_union(fmod)
DynamicSparseNumberBase_decl_std_binary_union(remainder) // TODO: optimize this
DynamicSparseNumberBase_decl_stdfl_binary_union(fmax)
DynamicSparseNumberBase_decl_stdfl_binary_union(fmin)
DynamicSparseNumberBase_decl_stdfl_binary_union(fdim)
DynamicSparseNumberBase_decl_stdfl_binary_union(hypot)
DynamicSparseNumberBase_decl_fl_binary_union(atan2)
#endif // __cplusplus >= 201103L

#define DynamicSparseNumberBase_decl_std_unary_complex(funcname) \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs> \
inline auto funcname( \
    const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & in) \
    ->typename SubType<SubTypeArgs...>::template rebind<decltype(std::funcname( \
                                                            typename Data::value_type())), \
                                                        typename Indices::value_type>::other

DynamicSparseNumberBase_decl_std_unary_complex(real);
DynamicSparseNumberBase_decl_std_unary_complex(imag);
DynamicSparseNumberBase_decl_std_unary_complex(norm);
} // namespace std


#endif // METAPHYSICL_DYNAMICSPARSENUMBERARRAY_DECL_H
