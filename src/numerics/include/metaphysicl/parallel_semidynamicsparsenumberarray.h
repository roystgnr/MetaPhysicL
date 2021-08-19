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
    // We need an example for MPI_Address to use
    static const MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N> p;
    if (!example)
      example = &p;

#ifdef TIMPI_HAVE_MPI

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
        tmptype, 0, sizeof(MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N>), &_datatype));
    timpi_call_mpi(MPI_Type_free(&tmptype));

    this->commit();

#endif // TIMPI_HAVE_MPI
  }

  StandardType(
      const StandardType<MetaPhysicL::SemiDynamicSparseNumberArray<T, I, N>> & timpi_mpi_var(t))
  {
    timpi_call_mpi(MPI_Type_dup(t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }

  static const bool is_fixed_type = true;
};
}

#endif // METAPHYSICL_HAVE_TIMPI
#endif // METAPHYSICL_PARALLEL_SEMIDYNAMICSPARSENUMBERARRAY_H
