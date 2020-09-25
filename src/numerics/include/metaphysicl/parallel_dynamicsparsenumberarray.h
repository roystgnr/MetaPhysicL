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

#ifndef METAPHYSICL_PARALLEL_DYNAMICSPARSENUMBERARRAY_H
#define METAPHYSICL_PARALLEL_DYNAMICSPARSENUMBERARRAY_H

#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_TIMPI

#include "metaphysicl/dynamicsparsenumberarray.h"
#include "metaphysicl/metaphysicl_cast.h"

#include "timpi/standard_type.h"
#include "timpi/packing.h"

#include <cstring> // for std::memcpy

namespace TIMPI
{
using MetaPhysicL::DynamicSparseNumberArray;

template <typename T, typename I>
class StandardType<DynamicSparseNumberArray<T, I>> : public NotADataType
{
public:
  StandardType(const DynamicSparseNumberArray<T, I> *) {}
};

} // namespace TIMPI

namespace libMesh
{
namespace Parallel
{
using MetaPhysicL::DynamicSparseNumberArray;

template <typename T, typename I>
class Packing<DynamicSparseNumberArray<T, I>>
{
public:
  static_assert(TIMPI::StandardType<T>::is_fixed_type &&
                TIMPI::StandardType<I>::is_fixed_type,
                "We cannot pack DynamicSparseNumberArray unless its "
                "template arguments are fixed size types");

  typedef std::size_t buffer_type;

  template <typename OutputIter, typename Context>
  static void pack(const DynamicSparseNumberArray<T, I> & dsna,
                   OutputIter data_out,
                   const Context * context);

  template <typename Context>
  static unsigned int packable_size(const DynamicSparseNumberArray<T, I> & dsna,
                                    const Context * context);

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter);

  template <typename BufferIter, typename Context>
  static DynamicSparseNumberArray<T, I> unpack(BufferIter in, Context * ctx);

  static const unsigned int sizetypes_per_T = (sizeof(T) + sizeof(std::size_t) - 1) /
    sizeof(std::size_t);
  static const unsigned int sizetypes_per_I = (sizeof(I) + sizeof(std::size_t) - 1) /
    sizeof(std::size_t);
};

template <typename T, typename I>
template <typename Context>
unsigned int
Packing<DynamicSparseNumberArray<T, I>>::
packable_size(const DynamicSparseNumberArray<T, I> & dsna,
              const Context *)
{
  // A spot to hold the size
  return 1 +
    // and then the data
    MetaPhysicL::cast_int<unsigned int>(dsna.size()) * (sizetypes_per_T + sizetypes_per_I);
}

template <typename T, typename I>
template <typename BufferIter>
unsigned int
Packing<DynamicSparseNumberArray<T, I>>::
packed_size(BufferIter iter)
{
  // We recorded the size in the first buffer entry
  return 1 + MetaPhysicL::cast_int<unsigned int>((*iter) * (sizetypes_per_T + sizetypes_per_I));
}

template <typename T, typename I>
template <typename OutputIter, typename Context>
void
Packing<DynamicSparseNumberArray<T, I>>::
pack(const DynamicSparseNumberArray<T, I> & dsna,
     OutputIter data_out,
     const Context *)
{
  *data_out++ = dsna.size();

  const auto & dsna_data = dsna.nude_data();
  const auto & indices = dsna.nude_indices();

  for (const T datum : dsna_data)
  {
    std::size_t T_as_sizetypes[sizetypes_per_T];
    std::memcpy(T_as_sizetypes, &datum, sizeof(T));
    for (unsigned int i = 0; i != sizetypes_per_T; ++i)
      *data_out++ = T_as_sizetypes[i];
  }

  for (const I index : indices)
  {
    std::size_t I_as_sizetypes[sizetypes_per_I];
    std::memcpy(I_as_sizetypes, &index, sizeof(I));
    for (unsigned int i = 0; i != sizetypes_per_I; ++i)
      *data_out++ = I_as_sizetypes[i];
  }
}

template <typename T, typename I>
template <typename BufferIter, typename Context>
DynamicSparseNumberArray<T, I>
Packing<DynamicSparseNumberArray<T, I>>::
unpack(BufferIter in, Context *)
{
  DynamicSparseNumberArray<T, I> dsna;

  dsna.resize(*in++);

  auto & dsna_data = dsna.nude_data();
  auto & indices = dsna.nude_indices();

  for (T & datum : dsna_data)
  {
    std::memcpy(&datum, &(*in), sizeof(T));
    in += sizetypes_per_T;
  }

  for (I & index : indices)
  {
    std::memcpy(&index, &(*in), sizeof(I));
    in += sizetypes_per_I;
  }

  return dsna;
}
} // namespace Parallel
} // namespace libMesh
#endif // METAPHYSICL_HAVE_TIMPI
#endif // METAPHYSICL_PARALLEL_DYNAMICSPARSENUMBERARRAY_H
