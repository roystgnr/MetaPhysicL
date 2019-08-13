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

#ifndef METAPHYSICL_SEMIDYNAMICSPARSENUMBERARRAY_H
#define METAPHYSICL_SEMIDYNAMICSPARSENUMBERARRAY_H

#include "metaphysicl/semidynamicsparsenumberarray_decl.h"
#include "metaphysicl/dynamicsparsenumberbase.h"

namespace MetaPhysicL
{

template <typename T, typename I, typename N>
inline SemiDynamicSparseNumberArray<T, I, N>::SemiDynamicSparseNumberArray(const T & val)
{
  // Avoid unused variable warnings in opt mode.
  (void)val;
  // This makes no sense unless val is 0!
#ifndef NDEBUG
  if (val)
    throw std::domain_error("Cannot initialize DynamicSparseNumberArray with non-zero scalar");
#endif
}

template <typename T, typename I, typename N>
template <typename T2>
inline SemiDynamicSparseNumberArray<T, I, N>::SemiDynamicSparseNumberArray(const T2 & val)
{
  // Avoid unused variable warnings in opt mode.
  (void)val;
  // This makes no sense unless val is 0!
#ifndef NDEBUG
  if (val)
    throw std::domain_error("Cannot initialize DynamicSparseNumberArray with non-zero scalar");
#endif
}

template <typename T, typename I, typename N>
template <typename T2, typename I2>
inline SemiDynamicSparseNumberArray<T, I, N>::SemiDynamicSparseNumberArray(
    const SemiDynamicSparseNumberArray<T2, I2, N> & src)
  : DynamicSparseNumberBase<DynamicStdArrayWrapper<T, N>,
                            DynamicStdArrayWrapper<I, N>,
                            MetaPhysicL::SemiDynamicSparseNumberArray,
                            T,
                            I,
                            N>(src)
{
}

template <typename T, typename I, typename N>
template <typename T2, typename I2>
inline SemiDynamicSparseNumberArray<T, I, N>::SemiDynamicSparseNumberArray(
    SemiDynamicSparseNumberArray<T2, I2, N> && src)
  : DynamicSparseNumberBase<DynamicStdArrayWrapper<T, N>,
                            DynamicStdArrayWrapper<I, N>,
                            MetaPhysicL::SemiDynamicSparseNumberArray,
                            T,
                            I,
                            N>(src)
{
}

DynamicSparseNumberBase_op(SemiDynamicSparseNumberArray, +, Plus)       // Union)
DynamicSparseNumberBase_op(SemiDynamicSparseNumberArray, -, Minus)      // Union)
DynamicSparseNumberBase_op(SemiDynamicSparseNumberArray, *, Multiplies) // Intersection)
DynamicSparseNumberBase_op(SemiDynamicSparseNumberArray, /, Divides)    // First)

template <typename T, typename I, typename N>
inline typename RawType<SemiDynamicSparseNumberArray<T, I, N>>::value_type
RawType<SemiDynamicSparseNumberArray<T, I, N>>::value(
    const SemiDynamicSparseNumberArray<T, I, N> & a)
{
  value_type returnval;
  returnval.nude_indices() = a.nude_indices();

  std::size_t index_size = a.size();
  returnval.nude_data().resize(index_size);

  for (unsigned int i = 0; i != index_size; ++i)
    returnval.raw_at(i) = RawType<T>::value(a.raw_at(i));
  return returnval;
}
}

#endif // METAPHYSICL_SEMIDYNAMICSPARSENUMBERARRAY_H
