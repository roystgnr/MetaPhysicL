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
// $Id: core.h 37197 2013-02-21 05:49:09Z roystgnr $
//
//--------------------------------------------------------------------------


#ifndef METAPHYSICL_SPARSENUMBERUTILS_H
#define METAPHYSICL_SPARSENUMBERUTILS_H

#include <algorithm>
#include <functional>

namespace MetaPhysicL {

// Utilities


  template <typename SubFunctor, typename IndexSet,
            typename IndexSetOut, typename T, typename Tout>
  struct UnaryVectorFunctor {
    UnaryVectorFunctor(SubFunctor f, const T* in, Tout* out) :
      _subfunctor(f), _datain(in), _dataout(out) {}

    template <typename ValueType>
    inline void operator()() const {
      const unsigned int
        indexin  = IndexSet::template IndexOf<ValueType>::index,
        indexout = IndexSetOut::template IndexOf<ValueType>::index;
      _dataout[indexout] = _subfunctor(_datain[indexin]);
    }

  private:
    SubFunctor _subfunctor;
    const T* _datain;
    Tout* _dataout;
  };


  template <typename SubFunctor,
            typename IndexSet1, typename IndexSet2, typename IndexSetOut,
            typename T1, typename T2, typename Tout>
  struct BinaryVectorFunctor {
    BinaryVectorFunctor(SubFunctor f, const T1* in1, const T2* in2,
		       Tout* out) :
      _subfunctor(f), _datain1(in1), _datain2(in2), _dataout(out) {}

    template <typename ValueType>
    inline void operator()() const {
      const unsigned int
        indexin1 = IndexSet1::template IndexOf<ValueType>::index,
        indexin2 = IndexSet2::template IndexOf<ValueType>::index,
        indexout = IndexSetOut::template IndexOf<ValueType>::index;
      _dataout[indexout] = _subfunctor(_datain1[indexin1], _datain2[indexin2]);
    }

  private:
    SubFunctor _subfunctor;
    const T1* _datain1;
    const T2* _datain2;
    Tout* _dataout;
  };

} // namespace std


#endif // METAPHYSICL_SPARSENUMBERUTILS_H
