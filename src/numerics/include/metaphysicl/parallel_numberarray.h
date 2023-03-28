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

#ifndef METAPHYSICL_PARALLEL_NUMBERARRAY_H
#define METAPHYSICL_PARALLEL_NUMBERARRAY_H

#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_TIMPI

#include "metaphysicl/numberarray.h"

#include "timpi/op_function.h"
#include "timpi/standard_type.h"

namespace TIMPI
{
template <std::size_t N, typename T>
class StandardType<MetaPhysicL::NumberArray<N, T>> : public DataType
{
public:
  explicit StandardType(
      const MetaPhysicL::NumberArray<N, T> * example = nullptr)
  {
#ifdef TIMPI_HAVE_MPI
    static data_type static_type = MPI_DATATYPE_NULL;
    if (static_type == MPI_DATATYPE_NULL)
      {
        // We need an example for MPI_Address to use
        static const MetaPhysicL::NumberArray<N, T> p {0};
        if (!example)
          example = &p;

        static_assert(N > 0, "Zero-length NumberArray not supported by TIMPI");
        StandardType<T> T_type(&((*example)[0]));

        int blocklength = N;
        MPI_Aint displs, start;
        MPI_Datatype tmptype, type = T_type;

        timpi_call_mpi
          (MPI_Get_address (example, &start));
        timpi_call_mpi
          (MPI_Get_address (&((*example)[0]), &displs));

        // subtract off offset to first value from the beginning of the structure
        displs -= start;

        // create a prototype structure
        timpi_call_mpi
          (MPI_Type_create_struct (1, &blocklength, &displs, &type,
                                   &tmptype));
        timpi_call_mpi
          (MPI_Type_commit (&tmptype));

        // resize the structure type to account for padding, if any
        timpi_call_mpi
          (MPI_Type_create_resized (tmptype, 0, sizeof(MetaPhysicL::NumberArray<N,T>),
                                    &static_type));

        timpi_call_mpi
          (MPI_Type_free (&tmptype));

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

# define METAPHYSICL_NA_MPI_BINARY(funcname) \
static inline void \
timpi_mpi_metaphysicl_na_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  const MetaPhysicL::NumberArray<N,T> * in = \
    static_cast<MetaPhysicL::NumberArray<N,T> *>(a); \
  MetaPhysicL::NumberArray<N,T> * inout = \
    static_cast<MetaPhysicL::NumberArray<N,T> *>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname(in[i],inout[i]); \
}


// MPI maxloc and minloc don't work in the array context - operator<
// returns an array of bools, but MPI doesn't let us return an array
// of locations
// # define METAPHYSICL_NA_MPI_PAIR_LOCATOR(funcname)


# define METAPHYSICL_NA_MPI_PAIR_BINARY_FUNCTOR(funcname) \
static inline void \
timpi_mpi_metaphysicl_na_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  typedef MetaPhysicL::NumberArray<N,T> dtype; \
 \
  const dtype * in = \
    static_cast<dtype *>(a); \
  dtype * inout = \
    static_cast<dtype *>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname<dtype>()(in[i], inout[i]); \
}


  template<std::size_t N, typename T>
  class OpFunction<MetaPhysicL::NumberArray<N,T>>
  {
    METAPHYSICL_NA_MPI_BINARY(max)
    METAPHYSICL_NA_MPI_BINARY(min)
    // METAPHYSICL_NA_MPI_PAIR_LOCATOR(max)
    // METAPHYSICL_NA_MPI_PAIR_LOCATOR(min)
    METAPHYSICL_NA_MPI_PAIR_BINARY_FUNCTOR(plus)
    METAPHYSICL_NA_MPI_PAIR_BINARY_FUNCTOR(multiplies)

  public:
    TIMPI_MPI_OPFUNCTION(max, metaphysicl_na_max)
    TIMPI_MPI_OPFUNCTION(min, metaphysicl_na_min)
    TIMPI_MPI_OPFUNCTION(sum, metaphysicl_na_plus)
    TIMPI_MPI_OPFUNCTION(product, metaphysicl_na_multiplies)

    // TIMPI_MPI_OPFUNCTION(max_location, metaphysicl_na_max_location)
    // TIMPI_MPI_OPFUNCTION(min_location, metaphysicl_na_min_location)
  };
# else // TIMPI_HAVE_MPI
  template<std::size_t N, typename T>
  class OpFunction<MetaPhysicL::NumberArray<N,T>> {};
#endif




}

#endif // METAPHYSICL_HAVE_TIMPI
#endif // METAPHYSICL_PARALLEL_NUMBERARRAY_H
