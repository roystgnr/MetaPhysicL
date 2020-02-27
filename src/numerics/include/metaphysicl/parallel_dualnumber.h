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

#ifndef METAPHYSICL_PARALLEL_DUALNUMBER_H
#define METAPHYSICL_PARALLEL_DUALNUMBER_H

#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_TIMPI

#include "metaphysicl/dualnumber.h"
#include "metaphysicl/metaphysicl_cast.h"

#include "timpi/standard_type.h"
#include "timpi/packing.h"

#include <cstring> // for std::memcpy

namespace TIMPI
{
using MetaPhysicL::DualNumber;
using MetaPhysicL::DynamicSparseNumberArray;

// Handle fixed size DualNumber first

template <typename T, typename D>
class StandardType<
    DualNumber<T, D>,
    typename std::enable_if<StandardType<T>::is_fixed_type && StandardType<D>::is_fixed_type>::type>
  : public DataType
{
public:
  explicit StandardType(const DualNumber<T, D> * example = nullptr)
  {
    // We need an example for MPI_Address to use
    static const DualNumber<T, D> p;
    if (!example)
      example = &p;

#ifdef TIMPI_HAVE_MPI

    // Get the sub-data-types, and make sure they live long enough
    // to construct the derived type
    StandardType<T> d1(&example->value());
    StandardType<D> d2(&example->derivatives());

    MPI_Datatype types[] = {(data_type)d1, (data_type)d2};
    int blocklengths[] = {1, 1};
    MPI_Aint displs[2], start;

    timpi_call_mpi(MPI_Get_address(const_cast<DualNumber<T, D> *>(example), &start));
    timpi_call_mpi(MPI_Get_address(const_cast<T *>(&example->value()), &displs[0]));
    timpi_call_mpi(MPI_Get_address(const_cast<D *>(&example->derivatives()), &displs[1]));
    displs[0] -= start;
    displs[1] -= start;

    // create a prototype structure
    MPI_Datatype tmptype;
    timpi_call_mpi(MPI_Type_create_struct(2, blocklengths, displs, types, &tmptype));
    timpi_call_mpi(MPI_Type_commit(&tmptype));

    // resize the structure type to account for padding, if any
    timpi_call_mpi(MPI_Type_create_resized(tmptype, 0, sizeof(DualNumber<T, D>), &_datatype));
    timpi_call_mpi(MPI_Type_free(&tmptype));

    this->commit();

#endif // TIMPI_HAVE_MPI
  }

  StandardType(const StandardType<DualNumber<T, D>> & timpi_mpi_var(t))
  {
    timpi_call_mpi(MPI_Type_dup(t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }

  static const bool is_fixed_type = true;
};

// Not fixed size DualNumber
template <typename T, typename D>
class StandardType<DualNumber<T, D>,
                   typename std::enable_if<!(StandardType<T>::is_fixed_type &&
                                             StandardType<D>::is_fixed_type)>::type>
{
public:
  static const bool is_fixed_type = false;

private:
  /**
   * we make the constructor private to catch mistakes at compile-time rather than link-time.
   */
  StandardType(const DualNumber<T, D> * example = nullptr);
};
} // namespace TIMPI

namespace libMesh
{
namespace Parallel
{
using MetaPhysicL::DualNumber;

template <typename T, typename D>
class Packing<DualNumber<T, D>,
              typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D>>::is_fixed_type>::type>
{
public:
  typedef std::size_t buffer_type;

  template <typename OutputIter, typename Context>
  static void pack(const DualNumber<T, D> & dn, OutputIter data_out, const Context * context);

  template <typename Context>
  static unsigned int packable_size(const DualNumber<T, D> & dn, const Context * context);

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter);

  template <typename BufferIter, typename Context>
  static DualNumber<T, D> unpack(BufferIter in, Context * ctx);

private:
  template <typename T2>
  struct IsFixed
  {
    static const bool value = TIMPI::StandardType<T2>::is_fixed_type;
  };
  template <typename T2>
  struct SizeTypesPer
  {
    static const unsigned int value = (sizeof(T2) + sizeof(std::size_t) - 1) / sizeof(std::size_t);
  };

  template <typename T2,
            typename Context,
            typename std::enable_if<IsFixed<T2>::value, int>::type = 0>
  static unsigned int packable_size_comp(const T2 &, const Context *)
  {
    return SizeTypesPer<T2>::value;
  }

  template <typename T2,
            typename Context,
            typename std::enable_if<!IsFixed<T2>::value, int>::type = 0>
  static unsigned int packable_size_comp(const T2 & comp, const Context * ctx)
  {
    return Packing<T2>::packable_size(comp, ctx);
  }

  template <typename T2,
            typename OutputIter,
            typename Context,
            typename std::enable_if<IsFixed<T2>::value, int>::type = 0>
  static void pack_comp(const T2 & comp, OutputIter data_out, const Context *)
  {
    std::size_t T2_as_sizetypes[SizeTypesPer<T2>::value];
    std::memcpy(T2_as_sizetypes, &comp, sizeof(T2));
    for (unsigned int i = 0; i != SizeTypesPer<T2>::value; ++i)
      *data_out++ = T2_as_sizetypes[i];
  }

  template <typename T2,
            typename OutputIter,
            typename Context,
            typename std::enable_if<!IsFixed<T2>::value, int>::type = 0>
  static void pack_comp(const T2 & comp, OutputIter data_out, const Context * ctx)
  {
    Packing<T2>::pack(comp, data_out, ctx);
  }

  template <typename T2,
            typename BufferIter,
            typename Context,
            typename std::enable_if<IsFixed<T2>::value, int>::type = 0>
  static void unpack_comp(T2 & comp, BufferIter in, Context *)
  {
    std::memcpy(&comp, &(*in), sizeof(T2));
  }

  template <typename T2,
            typename BufferIter,
            typename Context,
            typename std::enable_if<!IsFixed<T2>::value, int>::type = 0>
  static void unpack_comp(T2 & comp, BufferIter in, Context * ctx)
  {
    comp = Packing<T2>::unpack(in, ctx);
  }
};

template <typename T, typename D>
template <typename Context>
unsigned int
Packing<DualNumber<T, D>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D>>::is_fixed_type>::type>::
    packable_size(const DualNumber<T, D> & dn, const Context * ctx)
{
  return 1 + packable_size_comp(dn.value(), ctx) + packable_size_comp(dn.derivatives(), ctx);
}

template <typename T, typename D>
template <typename BufferIter>
unsigned int
Packing<DualNumber<T, D>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D>>::is_fixed_type>::type>::
    packed_size(BufferIter iter)
{
  // We recorded the size in the first buffer entry
  return *iter;
}

template <typename T, typename D>
template <typename OutputIter, typename Context>
void
Packing<DualNumber<T, D>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D>>::is_fixed_type>::type>::
    pack(const DualNumber<T, D> & dn, OutputIter data_out, const Context * ctx)
{
  unsigned int size = packable_size(dn, ctx);

  // First write out info about the buffer size
  *data_out++ = MetaPhysicL::cast_int<std::size_t>(size);

  // Now pack the data
  pack_comp(dn.value(), data_out, ctx);

  // TIMPI uses a back_inserter for `pack_range` so we don't (and can't)
  // actually increment the iterator with operator+=. operator++ is a no-op
  //
  // data_out += packable_size_comp(dn.value(), ctx);

  pack_comp(dn.derivatives(), data_out, ctx);
}

template <typename T, typename D>
template <typename BufferIter, typename Context>
DualNumber<T, D>
Packing<DualNumber<T, D>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D>>::is_fixed_type>::type>::
    unpack(BufferIter in, Context * ctx)
{
  DualNumber<T, D> dn;

  // We don't care about the size
  in++;

  // Unpack the data
  unpack_comp(dn.value(), in, ctx);

  // Make sure we increment the iterator
  in += packable_size_comp(dn.value(), ctx);

  unpack_comp(dn.derivatives(), in, ctx);

  return dn;
}

} // namespace Parallel
} // namespace libMesh

#endif // METAPHYSICL_HAVE_TIMPI
#endif // METAPHYSICL_PARALLEL_DUALNUMBER_H
