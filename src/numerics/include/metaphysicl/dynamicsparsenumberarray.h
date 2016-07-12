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


#ifndef METAPHYSICL_DYNAMICSPARSENUMBERARRAY_H
#define METAPHYSICL_DYNAMICSPARSENUMBERARRAY_H

#include "metaphysicl/dynamicsparsenumberbase.h"

namespace MetaPhysicL {

// Forward declarations

// Data type T, index type I
template <typename T, typename I>
class DynamicSparseNumberArray;

// Helper structs

template<typename I1, typename I2, typename S, typename T, bool reverseorder>
struct DotType<DynamicSparseNumberArray<S,I1>,
               DynamicSparseNumberArray<T,I2>, reverseorder> {
  typedef
    DynamicSparseNumberArray
      <typename DotType<S,T,reverseorder>::supertype,
       typename CompareTypes<I1, I2>::supertype>
      supertype;
};

template<typename I1, typename I2, typename S, typename T, bool reverseorder>
struct OuterProductType<DynamicSparseNumberArray<S, I1>,
                        DynamicSparseNumberArray<T, I2>, reverseorder> {
  typedef
    DynamicSparseNumberArray
      <typename OuterProductType<S,T,reverseorder>::supertype,
       typename CompareTypes<I1, I2>::supertype>
      supertype;
};

template<typename S, typename I>
struct SumType<DynamicSparseNumberArray<S, I> > {
  typedef DynamicSparseNumberArray<typename SumType<S>::supertype, I> supertype;
};


template <typename T, typename I>
class DynamicSparseNumberArray :
  public DynamicSparseNumberBase<T,I,DynamicSparseNumberArray>,
  public safe_bool<DynamicSparseNumberArray<T,I> >
{
public:
  template <typename T2>
  struct rebind {
    typedef DynamicSparseNumberArray<T2, I> other;
  };

  DynamicSparseNumberArray() {}

  DynamicSparseNumberArray(const T& val) {
    // This makes no sense unless val is 0!
#ifndef NDEBUG
    if (val)
      throw std::domain_error("Cannot initialize DynamicSparseNumberArray with non-zero scalar");
#endif
  }

  template <typename T2>
  DynamicSparseNumberArray(const T2& val) {
    // This makes no sense unless val is 0!
#ifndef NDEBUG
    if (val)
      throw std::domain_error("Cannot initialize DynamicSparseNumberArray with non-zero scalar");
#endif
  }

#if __cplusplus >= 201103L
  // Move constructors are useful when all your data is on the heap
  DynamicSparseNumberArray(DynamicSparseNumberArray<T, I> && src) = default;

  // Move assignment avoids heap operations too
  DynamicSparseNumberArray& operator= (DynamicSparseNumberArray<T, I> && src) = default;

  // Standard copy operations get implicitly deleted upon move
  // constructor definition, so we manually enable them.
  DynamicSparseNumberArray(const DynamicSparseNumberArray<T, I> & src) = default;

  DynamicSparseNumberArray& operator= (const DynamicSparseNumberArray<T, I> & src) = default;
#endif

  template <typename T2, typename I2>
  DynamicSparseNumberArray(DynamicSparseNumberArray<T2, I2> src) :
    DynamicSparseNumberBase<T,I,MetaPhysicL::DynamicSparseNumberArray>(src) {}


  template <typename T2, typename I2>
  DynamicSparseNumberArray
    <typename DotType<T,T2>::supertype,
     typename CompareTypes<I, I2>::supertype>
  dot (const DynamicSparseNumberArray<T2,I2>& a) const
  {
    typedef typename DotType<T,T2>::supertype TS;
    typedef typename CompareTypes<I, I2>::supertype IS;

    DynamicSparseNumberArray<TS, IS> returnval;

    // FIXME
    metaphysicl_not_implemented();

    return returnval;
  }

  template <typename T2, typename I2>
  DynamicSparseNumberArray<
    typename OuterProductType<T,T2>::supertype,
    typename CompareTypes<I, I2>::supertype>
  outerproduct (const DynamicSparseNumberArray<T2, I2>& a) const
  {
    typedef typename OuterProductType<T,T2>::supertype TS;
    typedef typename CompareTypes<I, I2>::supertype IS;
    DynamicSparseNumberArray<TS, IS> returnval;

    // FIXME
    metaphysicl_not_implemented();

    return returnval;
  }
};


//
// Non-member functions
//

template <unsigned int N,
          unsigned int index1=0, typename Data1=void,
          unsigned int index2=0, typename Data2=void,
          unsigned int index3=0, typename Data3=void,
          unsigned int index4=0, typename Data4=void,
          unsigned int index5=0, typename Data5=void,
          unsigned int index6=0, typename Data6=void,
          unsigned int index7=0, typename Data7=void,
          unsigned int index8=0, typename Data8=void>
struct DynamicSparseNumberArrayOf
{
  typedef
  typename SymmetricCompareTypes<Data1,
    typename SymmetricCompareTypes<Data2,
      typename SymmetricCompareTypes<Data3,
        typename SymmetricCompareTypes<Data4,
          typename SymmetricCompareTypes<Data5,
            typename SymmetricCompareTypes<Data6,
              typename SymmetricCompareTypes<Data7,Data8>::supertype
            >::supertype
          >::supertype
        >::supertype
      >::supertype
    >::supertype
  >::supertype supertype;

  typedef DynamicSparseNumberArray<supertype, unsigned int> type;
};



template <std::size_t N, unsigned int index, typename T>
struct DynamicSparseNumberArrayUnitVector
{
  typedef DynamicSparseNumberArray<T, unsigned int> type;

  static type value() {
    type returnval;
    returnval.resize(1);
    returnval.raw_at(0) = 1;
    returnval.raw_index(0) = index;
    return returnval;
  }
};


template <std::size_t N, typename T>
struct DynamicSparseNumberArrayFullVector
{
  typedef DynamicSparseNumberArray<T,unsigned int> type;

  static type value() {
    type returnval;
    returnval.resize(N);
    for (unsigned int i=0; i != N; ++i)
      {
        returnval.raw_at(i) = 1;
        returnval.raw_index(i) = i;
      }
    return returnval;
  }
};



template <typename T, typename I, typename I2>
inline
DynamicSparseNumberArray<DynamicSparseNumberArray<T, I>, I2>
transpose(const DynamicSparseNumberArray<DynamicSparseNumberArray<T, I2>, I>& /*a*/)
{
  DynamicSparseNumberArray<DynamicSparseNumberArray<T, I>, I2> returnval;

  metaphysicl_not_implemented();

  return returnval;
}


template <typename T, typename I>
DynamicSparseNumberArray<typename SumType<T>::supertype, I>
sum (const DynamicSparseNumberArray<T, I> &a)
{
  std::size_t index_size = a.size();

  DynamicSparseNumberArray<typename SumType<T>::supertype, I>
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i).sum();
    returnval.raw_index(i) = a.raw_index(i);
  }

  return returnval;
}



DynamicSparseNumberBase_op(DynamicSparseNumberArray, +, Plus)       // Union)
DynamicSparseNumberBase_op(DynamicSparseNumberArray, -, Minus)      // Union)
DynamicSparseNumberBase_op(DynamicSparseNumberArray, *, Multiplies) // Intersection)
DynamicSparseNumberBase_op(DynamicSparseNumberArray, /, Divides)    // First)


// CompareTypes, RawType, ValueType specializations

DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, CompareTypes);
DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, PlusType);
DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, MinusType);
DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, MultipliesType);
DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, DividesType);
DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, AndType);
DynamicSparseNumberBase_comparisons(DynamicSparseNumberArray, OrType);


template <typename T, typename I>
struct RawType<DynamicSparseNumberArray<T, I> >
{
  typedef DynamicSparseNumberArray<typename RawType<T>::value_type, I> value_type;

  static value_type value(const DynamicSparseNumberArray<T, I>& a)
    {
      value_type returnval;
      returnval.nude_indices() = a.nude_indices();

      std::size_t index_size = a.size();
      returnval.nude_data().resize(index_size);

      for (unsigned int i=0; i != index_size; ++i)
        returnval.raw_at(i) = RawType<T>::value(a.raw_at(i));
      return returnval;
    }
};

template <typename T, typename I>
struct ValueType<DynamicSparseNumberArray<T, I> >
{
  typedef typename ValueType<T>::type type;
};

} // namespace MetaPhysicL


namespace std {

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DynamicSparseNumberArray;
using MetaPhysicL::SymmetricCompareTypes;

#define DynamicSparseNumberArray_std_unary(funcname) \
template <typename T, typename I> \
inline \
DynamicSparseNumberArray<T, I> \
funcname (DynamicSparseNumberArray<T, I> a) \
{ \
  std::size_t index_size = a.size(); \
  for (unsigned int i=0; i != index_size; ++i) \
    a.raw_at(i) = std::funcname(a.raw_at(i)); \
 \
  return a; \
}


#define DynamicSparseNumberArray_std_binary_union(funcname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, \
                         typename CompareTypes<I,I2>::supertype> \
funcname (const DynamicSparseNumberArray<T, I>& a, \
          const DynamicSparseNumberArray<T2, I2>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<I,I2>::supertype IS; \
  DynamicSparseNumberArray<TS, IS> returnval; \
 \
  std::size_t index_size = a.nude_indices.size(); \
  returnval.nude_indices = a.nude_indices; \
  returnval.nude_data.resize(index_size); \
  returnval.sparsity_union(b.nude_indices); \
 \
  typename std::vector<I>::const_iterator  index_a_it = a.nude_indices.begin(); \
  typename std::vector<I2>::const_iterator index_b_it = b.nude_indices.begin(); \
  typename std::vector<IS>::iterator     index_out_it = returnval.nude_indices.begin(); \
 \
  typename std::vector<T>::const_iterator  data_a_it = a.nude_data.begin(); \
  typename std::vector<T2>::const_iterator data_b_it = b.nude_data.begin(); \
  typename std::vector<TS>::iterator     data_out_it = returnval.nude_data.begin(); \
 \
  const IS  maxIS  = std::numeric_limits<IS>::max(); \
 \
  for (; index_out_it != returnval.nude_indices.end(); ++index_out_it, ++data_out_it) { \
    const IS index_a = (index_a_it == a.nude_indices.end()) ? maxIS : *index_a_it; \
    const IS index_b = (index_b_it == b.nude_indices.end()) ? maxIS : *index_b_it; \
    const IS index_out = *index_out_it; \
    const TS data_a  = (index_a_it == a.nude_indices.end()) ? 0: *data_a_it; \
    const TS data_b  = (index_b_it == b.nude_indices.end()) ? 0: *data_b_it; \
    typename std::vector<TS>::reference data_out = *data_out_it; \
 \
    if (index_a == index_out) { \
      if (index_b == index_out) { \
        data_out = std::funcname(data_a, data_b); \
        index_b_it++; \
        data_b_it++; \
      } else { \
        data_out = std::funcname(data_a, 0); \
      } \
      index_a_it++; \
      data_a_it++; \
    } else { \
      metaphysicl_assert_equal_to(index_b, index_out); \
      data_out = std::funcname(0, data_b); \
      index_b_it++; \
      data_b_it++; \
    } \
  } \
 \
  return returnval; \
} \
 \
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const DynamicSparseNumberArray<T, I>& a, const T2& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  DynamicSparseNumberArray<TS, I> returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices = a.nude_indices; \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = std::funcname(a.raw_at(i), b); \
 \
  return returnval; \
} \
 \
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const T& a, const DynamicSparseNumberArray<T2, I>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  DynamicSparseNumberArray<TS, I> returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices = b.nude_indices; \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = std::funcname(a, b.raw_at(i)); \
 \
  return returnval; \
}


// Pow needs its own specialization, both to avoid being confused by
// pow<T1,T2> and because pow(x,0) isn't 0.
template <typename T, typename T2, typename I>
inline
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, I>
pow (const DynamicSparseNumberArray<T, I>& a, const T2& b)
{
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS;
  DynamicSparseNumberArray<TS, I> returnval;

  std::size_t index_size = a.size();
  returnval.nude_indices() = a.nude_indices();
  returnval.nude_data().resize(index_size);

  for (unsigned int i=0; i != index_size; ++i)
    returnval.raw_at(i) = std::pow(a.raw_at(i), b);

  return returnval;
}


// NOTE: unary functions for which f(0) != 0 are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

// DynamicSparseNumberArray_std_binary(pow) // separate definition
// DynamicSparseNumberArray_std_unary(exp)
// DynamicSparseNumberArray_std_unary(log)
// DynamicSparseNumberArray_std_unary(log10)
DynamicSparseNumberArray_std_unary(sin)
// DynamicSparseNumberArray_std_unary(cos)
DynamicSparseNumberArray_std_unary(tan)
DynamicSparseNumberArray_std_unary(asin)
// DynamicSparseNumberArray_std_unary(acos)
DynamicSparseNumberArray_std_unary(atan)
DynamicSparseNumberArray_std_binary_union(atan2)
DynamicSparseNumberArray_std_unary(sinh)
// DynamicSparseNumberArray_std_unary(cosh)
DynamicSparseNumberArray_std_unary(tanh)
DynamicSparseNumberArray_std_unary(sqrt)
DynamicSparseNumberArray_std_unary(abs)
DynamicSparseNumberArray_std_unary(fabs)
DynamicSparseNumberArray_std_binary_union(max)
DynamicSparseNumberArray_std_binary_union(min)
DynamicSparseNumberArray_std_unary(ceil)
DynamicSparseNumberArray_std_unary(floor)
DynamicSparseNumberArray_std_binary_union(fmod) // TODO: optimize this


template <typename T, typename I>
class numeric_limits<DynamicSparseNumberArray<T, I> > :
  public MetaPhysicL::raw_numeric_limits<DynamicSparseNumberArray<T, I>, T> {};

} // namespace std


#endif // METAPHYSICL_DYNAMICSPARSENUMBERARRAY_H
