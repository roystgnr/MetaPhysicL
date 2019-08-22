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

#ifndef METAPHYSICL_DUALSEMIDYNAMICSPARSENUMBERARRAY_DECL_H
#define METAPHYSICL_DUALSEMIDYNAMICSPARSENUMBERARRAY_DECL_H

#include "metaphysicl/dualnumber_decl.h"
#include "metaphysicl/semidynamicsparsenumberarray_decl.h"

namespace MetaPhysicL
{

template <typename T, typename I, typename N>
struct DerivativeType<SemiDynamicSparseNumberArray<T, I, N>>
{
  typedef SemiDynamicSparseNumberArray<typename DerivativeType<T>::type, I, N> type;
};

template <typename T, typename I, typename N>
struct DerivativesType<SemiDynamicSparseNumberArray<T, I, N>>
{
  typedef SemiDynamicSparseNumberArray<typename DerivativesType<T>::type, I, N> type;
};

template <typename T, typename I, typename N>
inline typename DerivativeType<SemiDynamicSparseNumberArray<T, I, N>>::type
derivative(const SemiDynamicSparseNumberArray<T, I, N> & a, unsigned int derivativeindex);

template <typename T, typename I, typename N>
inline typename DerivativesType<SemiDynamicSparseNumberArray<T, I, N>>::type
derivatives(const SemiDynamicSparseNumberArray<T, I, N> & a);

template <typename T, typename I, typename N, unsigned int derivativeindex>
struct DerivativeOf<SemiDynamicSparseNumberArray<T, I, N>, derivativeindex>
{
  static typename DerivativeType<SemiDynamicSparseNumberArray<T, I, N>>::type
  derivative(const SemiDynamicSparseNumberArray<T, I, N> & a);
};

// For a vector of values a[i] each of which has a defined gradient,
// the divergence is the sum of derivative_wrt_xi(a[i])

// For a tensor of values, we take the divergence with respect to the
// first index.
template <typename T, typename I, typename N>
inline typename DerivativeType<T>::type divergence(const SemiDynamicSparseNumberArray<T, I, N> & a);

// For a vector of values, the gradient is going to be a tensor
template <typename T, typename I, typename N>
inline SemiDynamicSparseNumberArray<typename T::derivatives_type, I, N>
gradient(const SemiDynamicSparseNumberArray<T, I, N> & a);

// DualNumber is subordinate to SemiDynamicSparseNumberArray

#define DualSemiDynamicSparseNumberArray_comparisons(templatename)                                 \
  template <typename T, typename T2, typename D, typename I, typename N, bool reverseorder>        \
  struct templatename<SemiDynamicSparseNumberArray<T2, I, N>, DualNumber<T, D>, reverseorder>      \
  {                                                                                                \
    typedef SemiDynamicSparseNumberArray<                                                          \
        typename Symmetric##templatename<T2, DualNumber<T, D>, reverseorder>::supertype,           \
        I,                                                                                         \
        N>                                                                                         \
        supertype;                                                                                 \
  }

DualSemiDynamicSparseNumberArray_comparisons(CompareTypes);
DualSemiDynamicSparseNumberArray_comparisons(PlusType);
DualSemiDynamicSparseNumberArray_comparisons(MinusType);
DualSemiDynamicSparseNumberArray_comparisons(MultipliesType);
DualSemiDynamicSparseNumberArray_comparisons(DividesType);
DualSemiDynamicSparseNumberArray_comparisons(AndType);
DualSemiDynamicSparseNumberArray_comparisons(OrType);

} // namespace MetaPhysicL

#endif // METAPHYSICL_DUALSEMIDYNAMICSPARSENUMBERARRAY_DECL_H
