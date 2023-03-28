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

#ifndef METAPHYSICL_PARALLEL_DYNAMIC_STD_ARRAY_WRAPPER_H
#define METAPHYSICL_PARALLEL_DYNAMIC_STD_ARRAY_WRAPPER_H

#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_TIMPI

#include "metaphysicl/dynamic_std_array_wrapper.h"
#include "timpi/standard_type.h"

namespace TIMPI
{
template <typename T, typename NType>
class StandardType<MetaPhysicL::DynamicStdArrayWrapper<T, NType>> : public DataType
{
public:
  explicit StandardType(const MetaPhysicL::DynamicStdArrayWrapper<T, NType> * example = nullptr)
  {
#ifdef TIMPI_HAVE_MPI
    static data_type static_type = MPI_DATATYPE_NULL;
    if (static_type == MPI_DATATYPE_NULL)
      {
        // We need an example for MPI_Address to use
        static const MetaPhysicL::DynamicStdArrayWrapper<T, NType> p{};
        if (!example)
          example = &p;

        // Get the sub-data-types, and make sure they live long enough
        // to construct the derived type
        StandardType<std::array<T, NType::size>> d1(&example->_data);
        StandardType<std::size_t> d2(&example->_dynamic_n);

        MPI_Datatype types[] = {(data_type)d1, (data_type)d2};
        int blocklengths[] = {1, 1};
        MPI_Aint displs[2], start;

        timpi_call_mpi(MPI_Get_address(
            const_cast<MetaPhysicL::DynamicStdArrayWrapper<T, NType> *>(example), &start));
        timpi_call_mpi(
            MPI_Get_address(const_cast<std::array<T, NType::size> *>(&example->_data), &displs[0]));
        timpi_call_mpi(MPI_Get_address(const_cast<std::size_t *>(&example->_dynamic_n), &displs[1]));
        displs[0] -= start;
        displs[1] -= start;

        // create a prototype structure
        MPI_Datatype tmptype;
        timpi_call_mpi(MPI_Type_create_struct(2, blocklengths, displs, types, &tmptype));
        timpi_call_mpi(MPI_Type_commit(&tmptype));

        // resize the structure type to account for padding, if any
        timpi_call_mpi(MPI_Type_create_resized(
            tmptype, 0, sizeof(MetaPhysicL::DynamicStdArrayWrapper<T, NType>), &static_type));
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
}

#endif // METAPHYSICL_HAVE_TIMPI
#endif // METAPHYSICL_PARALLEL_DYNAMIC_STD_ARRAY_WRAPPER_H
