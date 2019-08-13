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

#ifndef METAPHYSICL_SEMIDYNAMICSPARSENUMBERARRAY_DECL_H
#define METAPHYSICL_SEMIDYNAMICSPARSENUMBERARRAY_DECL_H

#include "metaphysicl/dynamicsparsenumberbase_decl.h"
#include "metaphysicl/dynamic_std_array_wrapper.h"
#include "metaphysicl/testable.h"

namespace MetaPhysicL
{
template <typename T, typename I, typename N>
class SemiDynamicSparseNumberArray : public DynamicSparseNumberBase<DynamicStdArrayWrapper<T, N>,
                                                                    DynamicStdArrayWrapper<I, N>,
                                                                    SemiDynamicSparseNumberArray,
                                                                    T,
                                                                    I,
                                                                    N>,
                                     public safe_bool<SemiDynamicSparseNumberArray<T, I, N>>
{
public:
  template <typename T2, typename I2 = I>
  struct rebind
  {
    typedef SemiDynamicSparseNumberArray<T2, I2, N> other;
  };

  SemiDynamicSparseNumberArray() = default;

  SemiDynamicSparseNumberArray(const T & val);

  template <typename T2>
  SemiDynamicSparseNumberArray(const T2 & val);

  SemiDynamicSparseNumberArray(SemiDynamicSparseNumberArray && src) = default;

  SemiDynamicSparseNumberArray & operator=(SemiDynamicSparseNumberArray && src) = default;

  SemiDynamicSparseNumberArray(const SemiDynamicSparseNumberArray & src) = default;

  SemiDynamicSparseNumberArray & operator=(const SemiDynamicSparseNumberArray & src) = default;

  template <typename T2, typename I2>
  SemiDynamicSparseNumberArray(const SemiDynamicSparseNumberArray<T2, I2, N> & src);

  template <typename T2, typename I2>
  SemiDynamicSparseNumberArray(SemiDynamicSparseNumberArray<T2, I2, N> && src);
};

//
// Non-member functions
//

template <size_t N,
          unsigned int index1 = 0,
          typename Data1 = void,
          unsigned int index2 = 0,
          typename Data2 = void,
          unsigned int index3 = 0,
          typename Data3 = void,
          unsigned int index4 = 0,
          typename Data4 = void,
          unsigned int index5 = 0,
          typename Data5 = void,
          unsigned int index6 = 0,
          typename Data6 = void,
          unsigned int index7 = 0,
          typename Data7 = void,
          unsigned int index8 = 0,
          typename Data8 = void>
struct SemiDynamicSparseNumberArrayOf
{
  typedef typename SymmetricCompareTypes<
      Data1,
      typename SymmetricCompareTypes<
          Data2,
          typename SymmetricCompareTypes<
              Data3,
              typename SymmetricCompareTypes<
                  Data4,
                  typename SymmetricCompareTypes<
                      Data5,
                      typename SymmetricCompareTypes<
                          Data6,
                          typename SymmetricCompareTypes<Data7, Data8>::supertype>::supertype>::
                      supertype>::supertype>::supertype>::supertype>::supertype supertype;

  typedef SemiDynamicSparseNumberArray<supertype, unsigned int, NWrapper<N>> type;
};

DynamicSparseNumberBase_decl_op(SemiDynamicSparseNumberArray, +, Plus)       // Union)
DynamicSparseNumberBase_decl_op(SemiDynamicSparseNumberArray, -, Minus)      // Union)
DynamicSparseNumberBase_decl_op(SemiDynamicSparseNumberArray, *, Multiplies) // Intersection)
DynamicSparseNumberBase_decl_op(SemiDynamicSparseNumberArray, /, Divides)    // First)

// CompareTypes, RawType, ValueType specializations

DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, CompareTypes);
DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, PlusType);
DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, MinusType);
DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, MultipliesType);
DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, DividesType);
DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, AndType);
DynamicSparseNumberBase_comparisons(SemiDynamicSparseNumberArray, OrType);

template <typename T, typename I, typename N>
struct RawType<SemiDynamicSparseNumberArray<T, I, N>>
{
  typedef SemiDynamicSparseNumberArray<typename RawType<T>::value_type, I, N> value_type;

  static value_type value(const SemiDynamicSparseNumberArray<T, I, N> & a);
};

template <typename T, typename I, typename N>
struct ValueType<SemiDynamicSparseNumberArray<T, I, N>>
{
  typedef typename ValueType<T>::type type;
};

} // namespace MetaPhysicL

namespace std
{

using MetaPhysicL::SemiDynamicSparseNumberArray;

template <typename T, typename I, typename N>
class numeric_limits<SemiDynamicSparseNumberArray<T, I, N>>
  : public MetaPhysicL::raw_numeric_limits<SemiDynamicSparseNumberArray<T, I, N>, T>
{
};

} // namespace std

#endif // METAPHYSICL_SEMIDYNAMICSPARSENUMBERARRAY_DECL_H
