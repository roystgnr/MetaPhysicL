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

#ifndef METAPHYSICL_DYNAMIC_STD_ARRAY_WRAPPER_H
#define METAPHYSICL_DYNAMIC_STD_ARRAY_WRAPPER_H

#include "metaphysicl/metaphysicl_asserts.h"
#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_TIMPI
#include "timpi/standard_type_forward.h"
#endif

#include <array>

namespace MetaPhysicL
{
template <typename T, typename NType>
class DynamicStdArrayWrapper
{
public:
  static const std::size_t N = NType::size;

  typedef std::size_t size_type;

  typedef typename std::array<T, N>::iterator iterator;
  typedef typename std::array<T, N>::const_iterator const_iterator;
  typedef typename std::array<T, N>::reverse_iterator reverse_iterator;
  typedef typename std::array<T, N>::const_reverse_iterator const_reverse_iterator;
  typedef typename std::array<T, N>::value_type value_type;
  typedef typename std::array<T, N>::reference reference;
  typedef typename std::array<T, N>::const_reference const_reference;

  DynamicStdArrayWrapper(const DynamicStdArrayWrapper & src)
  {
    _dynamic_n = src._dynamic_n;
    metaphysicl_assert(_dynamic_n <= N);
    std::copy(src.begin(), src.end(), _data.begin());
  }

  // A std::array isn't movable but it's contents might be
  DynamicStdArrayWrapper(DynamicStdArrayWrapper && src)
  {
    _dynamic_n = src._dynamic_n;
    metaphysicl_assert(_dynamic_n <= N);
    auto src_it = src.begin(), src_end = src.end(), this_it = _data.begin();
    for (; src_it != src_end; ++src_it, ++this_it)
      *this_it = std::move(*src_it);
  }

  DynamicStdArrayWrapper & operator=(const DynamicStdArrayWrapper & src)
  {
    _dynamic_n = src._dynamic_n;
    metaphysicl_assert(_dynamic_n <= N);
    std::copy(src.begin(), src.end(), _data.begin());
    return *this;
  }

  // A std::array isn't movable but it's contents might be
  DynamicStdArrayWrapper & operator=(DynamicStdArrayWrapper && src)
  {
    _dynamic_n = src._dynamic_n;
    metaphysicl_assert(_dynamic_n <= N);
    auto src_it = src.begin(), src_end = src.end(), this_it = _data.begin();
    for (; src_it != src_end; ++src_it, ++this_it)
      *this_it = std::move(*src_it);
    return *this;
  }

  DynamicStdArrayWrapper() = default;

  iterator begin() { return _data.begin(); }

  const_iterator begin() const { return _data.begin(); }

  iterator end()
  {
    metaphysicl_assert(_dynamic_n <= N);
    return _data.begin() + _dynamic_n;
  }

  const_iterator end() const
  {
    metaphysicl_assert(_dynamic_n <= N);
    return _data.begin() + _dynamic_n;
  }

  T & operator[](size_type i)
  {
    metaphysicl_assert(i < N);
    return _data[i];
  }

  const T & operator[](size_type i) const
  {
    metaphysicl_assert(i < N);
    return _data[i];
  }

  size_type size() const { return _dynamic_n; }

  void resize(size_type new_size)
  {
    if (new_size > N)
      metaphysicl_error();
    _dynamic_n = new_size;
  }

  reverse_iterator rbegin()
  {
    metaphysicl_assert(_dynamic_n <= N);
    return _data.rend() - _dynamic_n;
  }

  const_reverse_iterator rbegin() const
  {
    metaphysicl_assert(_dynamic_n <= N);
    return _data.rend() - _dynamic_n;
  }

  reverse_iterator rend() { return _data.rend(); }

  const_reverse_iterator rend() const { return _data.rend(); }

protected:
#ifdef METAPHYSICL_HAVE_TIMPI
  friend class TIMPI::StandardType<DynamicStdArrayWrapper<T, NType>>;
#endif
  std::array<T, N> _data;
  std::size_t _dynamic_n = 0;
};

template <std::size_t N>
struct NWrapper
{
  static const std::size_t size = N;
};
}

#endif // METAPHYSICL_DYNAMIC_STD_ARRAY_WRAPPER_H
