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

#include "timpi/op_function.h"
#include "timpi/standard_type.h"
#include "timpi/packing_decl.h"

#include <cstring> // for std::memcpy

namespace TIMPI
{
using MetaPhysicL::DualNumber;

// Handle fixed size DualNumber first

template <typename T, typename D, bool asd>
class StandardType<
    DualNumber<T, D, asd>,
    typename std::enable_if<StandardType<T>::is_fixed_type && StandardType<D>::is_fixed_type>::type>
  : public DataType
{
public:
  explicit StandardType(const DualNumber<T, D, asd> * example = nullptr)
  {
#ifdef TIMPI_HAVE_MPI
    static data_type static_type = MPI_DATATYPE_NULL;
    if (static_type == MPI_DATATYPE_NULL)
      {
        // We need an example for MPI_Address to use
        static const DualNumber<T, D, asd> p(0);
        if (!example)
          example = &p;

        // Get the sub-data-types, and make sure they live long enough
        // to construct the derived type
        StandardType<T> d1(&example->value());
        StandardType<D> d2(&example->derivatives());

        MPI_Datatype types[] = {(data_type)d1, (data_type)d2};
        int blocklengths[] = {1, 1};
        MPI_Aint displs[2], start;

        timpi_call_mpi(MPI_Get_address(const_cast<DualNumber<T, D, asd> *>(example), &start));
        timpi_call_mpi(MPI_Get_address(const_cast<T *>(&example->value()), &displs[0]));
        timpi_call_mpi(MPI_Get_address(const_cast<D *>(&example->derivatives()), &displs[1]));
        displs[0] -= start;
        displs[1] -= start;

        // create a prototype structure
        MPI_Datatype tmptype;
        timpi_call_mpi(MPI_Type_create_struct(2, blocklengths, displs, types, &tmptype));
        timpi_call_mpi(MPI_Type_commit(&tmptype));

        // resize the structure type to account for padding, if any
        timpi_call_mpi(MPI_Type_create_resized(tmptype, 0, sizeof(DualNumber<T, D, asd>), &static_type));
        timpi_call_mpi(MPI_Type_free(&tmptype));

        SemiPermanent::add
          (std::make_unique<ManageType>(static_type));
      }
    _datatype = static_type;
#else
    timpi_ignore(example);
#endif // TIMPI_HAVE_MPI
  }

  StandardType(const StandardType<DualNumber<T, D, asd>> & timpi_mpi_var(t)) :
    DataType(t._datatype) {}

  StandardType & operator=(StandardType & t)
  {
    _datatype = t._datatype;
    return *this;
  }

  static const bool is_fixed_type = true;
};

// Not fixed size DualNumber
template <typename T, typename D, bool asd>
class StandardType<DualNumber<T, D, asd>,
                   typename std::enable_if<!(StandardType<T>::is_fixed_type &&
                                             StandardType<D>::is_fixed_type)>::type>
  : public NotADataType
{
public:
  StandardType(const DualNumber<T, D, asd> *) {}
};


#ifdef TIMPI_HAVE_MPI

# define METAPHYSICL_DUALNUMBER_MPI_BINARY(funcname) \
static inline void \
timpi_mpi_metaphysicl_dualnumber_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  const MetaPhysicL::DualNumber<T,D,asd> * in = \
    static_cast<MetaPhysicL::DualNumber<T,D,asd> *>(a); \
  MetaPhysicL::DualNumber<T,D,asd> * inout = \
    static_cast<MetaPhysicL::DualNumber<T,D,asd> *>(b); \
  for (int i=0; i != size; ++i) \
    { \
      inout[i].value() = std::funcname(in[i].value(),inout[i].value()); \
      inout[i].derivatives() = \
        std::funcname(in[i].derivatives(),inout[i].derivatives()); \
    } \
}

# define METAPHYSICL_DUALNUMBER_MPI_PAIR_LOCATOR(funcname) \
static inline void \
timpi_mpi_metaphysicl_dualnumber_##funcname##_location(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  typedef std::pair<MetaPhysicL::DualNumber<T,D,asd>, int> dtype; \
 \
  dtype *in = static_cast<dtype*>(a); \
  dtype *inout = static_cast<dtype*>(b); \
  for (int i=0; i != size; ++i) \
    { \
      MetaPhysicL::DualNumber<T,D,asd> old_inout = inout[i].first; \
      inout[i].first.value() = \
        std::funcname(in[i].first.value(), inout[i].first.value()); \
      inout[i].first.derivatives() = \
        std::funcname(in[i].first.derivatives(),inout[i].first.derivatives()); \
      if (old_inout != inout[i].first) \
        inout[i].second = in[i].second; \
    } \
}



# define METAPHYSICL_DUALNUMBER_MPI_PAIR_BINARY_FUNCTOR(funcname) \
static inline void \
timpi_mpi_metaphysicl_dualnumber_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  typedef MetaPhysicL::DualNumber<T,D,asd> dtype; \
 \
  const dtype * in = \
    static_cast<dtype *>(a); \
  dtype * inout = \
    static_cast<dtype *>(b); \
  for (int i=0; i != size; ++i) \
    { \
      inout[i].value()  = std::funcname<T>()(in[i].value(), inout[i].value()); \
      inout[i].derivatives() = \
        std::funcname<D>()(in[i].derivatives(),inout[i].derivatives()); \
    } \
}


  template<typename T, typename D, bool asd>
  class OpFunction<MetaPhysicL::DualNumber<T,D,asd>>
  {
    METAPHYSICL_DUALNUMBER_MPI_BINARY(max)
    METAPHYSICL_DUALNUMBER_MPI_BINARY(min)
    METAPHYSICL_DUALNUMBER_MPI_PAIR_LOCATOR(max)
    METAPHYSICL_DUALNUMBER_MPI_PAIR_LOCATOR(min)
    METAPHYSICL_DUALNUMBER_MPI_PAIR_BINARY_FUNCTOR(plus)
    METAPHYSICL_DUALNUMBER_MPI_PAIR_BINARY_FUNCTOR(multiplies)

  public:
    TIMPI_MPI_OPFUNCTION(max, metaphysicl_dualnumber_max)
    TIMPI_MPI_OPFUNCTION(min, metaphysicl_dualnumber_min)
    TIMPI_MPI_OPFUNCTION(sum, metaphysicl_dualnumber_plus)
    TIMPI_MPI_OPFUNCTION(product, metaphysicl_dualnumber_multiplies)

    TIMPI_MPI_OPFUNCTION(max_location, metaphysicl_dualnumber_max_location)
    TIMPI_MPI_OPFUNCTION(min_location, metaphysicl_dualnumber_min_location)
  };
# else // TIMPI_HAVE_MPI
  template<typename T, typename U>
  class OpFunction<MetaPhysicL::DualNumber<T,D,asd>> {};
#endif




} // namespace TIMPI

namespace libMesh
{
namespace Parallel
{
using MetaPhysicL::DualNumber;

template <typename T, typename D, bool asd>
class Packing<
    DualNumber<T, D, asd>,
    typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D, asd>>::is_fixed_type>::type>
{
public:
  typedef std::size_t buffer_type;

  template <typename OutputIter, typename Context>
  static void pack(const DualNumber<T, D, asd> & dn, OutputIter data_out, const Context * context);

  template <typename Context>
  static unsigned int packable_size(const DualNumber<T, D, asd> & dn, const Context * context);

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter);

  template <typename BufferIter, typename Context>
  static DualNumber<T, D, asd> unpack(BufferIter in, Context * ctx);

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

template <typename T, typename D, bool asd>
template <typename Context>
inline unsigned int
Packing<DualNumber<T, D, asd>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D, asd>>::is_fixed_type>::type>::
    packable_size(const DualNumber<T, D, asd> & dn, const Context * ctx)
{
  return 1 + packable_size_comp(dn.value(), ctx) + packable_size_comp(dn.derivatives(), ctx);
}

template <typename T, typename D, bool asd>
template <typename BufferIter>
inline unsigned int
Packing<DualNumber<T, D, asd>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D, asd>>::is_fixed_type>::type>::
    packed_size(BufferIter iter)
{
  // We recorded the size in the first buffer entry
  return *iter;
}

template <typename T, typename D, bool asd>
template <typename OutputIter, typename Context>
inline void
Packing<DualNumber<T, D, asd>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D, asd>>::is_fixed_type>::type>::
    pack(const DualNumber<T, D, asd> & dn, OutputIter data_out, const Context * ctx)
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

template <typename T, typename D, bool asd>
template <typename BufferIter, typename Context>
inline DualNumber<T, D, asd>
Packing<DualNumber<T, D, asd>,
        typename std::enable_if<!TIMPI::StandardType<DualNumber<T, D, asd>>::is_fixed_type>::type>::
    unpack(BufferIter in, Context * ctx)
{
  DualNumber<T, D, asd> dn;

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
