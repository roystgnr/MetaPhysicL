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

#ifndef METAPHYSICL_PARALLEL_SEMIDYNAMICSPARSENUMBERARRAY_H
#define METAPHYSICL_PARALLEL_SEMIDYNAMICSPARSENUMBERARRAY_H

#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_TIMPI

#include "metaphysicl/semidynamicsparsenumberarray.h"
#include "metaphysicl/parallel_dynamic_std_array_wrapper.h"

#include "timpi/op_function.h"
#include "timpi/standard_type.h"

namespace TIMPI
{
template <typename T, typename I, typename N>
class StandardType<MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N>> : public DataType
{
public:
  explicit StandardType(
      const MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N> * example = nullptr)
  {
#ifdef TIMPI_HAVE_MPI
    static data_type static_type = MPI_DATATYPE_NULL;
    if (static_type == MPI_DATATYPE_NULL)
      {
        // We need an example for MPI_Address to use
        static const MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N> p;
        if (!example)
          example = &p;

        // Get the sub-data-types, and make sure they live long enough
        // to construct the derived type
        StandardType<MetaPhysicL::DynamicStdArrayWrapper<T, N>> d1(&example->nude_data());
        StandardType<MetaPhysicL::DynamicStdArrayWrapper<I, N>> d2(&example->nude_indices());

        MPI_Datatype types[] = {(data_type)d1, (data_type)d2};
        int blocklengths[] = {1, 1};
        MPI_Aint displs[2], start;

        timpi_call_mpi(MPI_Get_address(
            const_cast<MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N> *>(example), &start));
        timpi_call_mpi(MPI_Get_address(
            const_cast<MetaPhysicL::DynamicStdArrayWrapper<T, N> *>(&example->nude_data()),
            &displs[0]));
        timpi_call_mpi(MPI_Get_address(
            const_cast<MetaPhysicL::DynamicStdArrayWrapper<I, N> *>(&example->nude_indices()),
            &displs[1]));
        displs[0] -= start;
        displs[1] -= start;

        // create a prototype structure
        MPI_Datatype tmptype;
        timpi_call_mpi(MPI_Type_create_struct(2, blocklengths, displs, types, &tmptype));
        timpi_call_mpi(MPI_Type_commit(&tmptype));

        // resize the structure type to account for padding, if any
        timpi_call_mpi(MPI_Type_create_resized(
            tmptype, 0, sizeof(MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N>), &static_type));
        timpi_call_mpi(MPI_Type_free(&tmptype));

        SemiPermanent::add
          (std::make_unique<ManageType>(static_type));
      }
    _datatype = static_type;
#else
    timpi_ignore(example);
#endif // TIMPI_HAVE_MPI
  }

  StandardType(const StandardType & t) :
    DataType(t._datatype) {}

  StandardType & operator=(StandardType & t)
  {
    _datatype = t._datatype;
    return *this;
  }

  static const bool is_fixed_type = true;
};


#ifdef TIMPI_HAVE_MPI

# define METAPHYSICL_SDSNA_MPI_BINARY(funcname) \
static inline void \
timpi_mpi_metaphysicl_sdsna_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  const MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N> * in = \
    static_cast<MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N> *>(a); \
  MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N> * inout = \
    static_cast<MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N> *>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname(in[i],inout[i]); \
}


// MPI maxloc and minloc don't work in the array context - operator<
// returns an array of bools, but MPI doesn't let us return an array
// of locations
// # define METAPHYSICL_SDSNA_MPI_PAIR_LOCATOR(funcname)


# define METAPHYSICL_SDSNA_MPI_PAIR_BINARY_FUNCTOR(funcname) \
static inline void \
timpi_mpi_metaphysicl_sdsna_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  typedef MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N> dtype; \
 \
  const dtype * in = \
    static_cast<dtype *>(a); \
  dtype * inout = \
    static_cast<dtype *>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname<dtype>()(in[i], inout[i]); \
}


  template<typename T, typename I, typename N>
  class OpFunction<MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N>>
  {
    METAPHYSICL_SDSNA_MPI_BINARY(max)
    METAPHYSICL_SDSNA_MPI_BINARY(min)
    // METAPHYSICL_SDSNA_MPI_PAIR_LOCATOR(max)
    // METAPHYSICL_SDSNA_MPI_PAIR_LOCATOR(min)
    METAPHYSICL_SDSNA_MPI_PAIR_BINARY_FUNCTOR(plus)
    METAPHYSICL_SDSNA_MPI_PAIR_BINARY_FUNCTOR(multiplies)

  public:
    TIMPI_MPI_OPFUNCTION(max, metaphysicl_sdsna_max)
    TIMPI_MPI_OPFUNCTION(min, metaphysicl_sdsna_min)
    TIMPI_MPI_OPFUNCTION(sum, metaphysicl_sdsna_plus)
    TIMPI_MPI_OPFUNCTION(product, metaphysicl_sdsna_multiplies)

    // TIMPI_MPI_OPFUNCTION(max_location, metaphysicl_sdsna_max_location)
    // TIMPI_MPI_OPFUNCTION(min_location, metaphysicl_sdsna_min_location)
  };
# else // TIMPI_HAVE_MPI
  template<typename T, typename I, typename N>
  class OpFunction<MetaPhysicL::SemiDynamicSparseNumberArray<T,I,N>> {};
#endif




}

#endif // METAPHYSICL_HAVE_TIMPI
#endif // METAPHYSICL_PARALLEL_SEMIDYNAMICSPARSENUMBERARRAY_H
