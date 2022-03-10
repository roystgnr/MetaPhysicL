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


#ifndef METAPHYSICL_DYNAMICSPARSENUMBERBASE_H
#define METAPHYSICL_DYNAMICSPARSENUMBERBASE_H

#include "metaphysicl/dynamicsparsenumberbase_decl.h"

namespace MetaPhysicL {

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
std::size_t
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::size() const
{ metaphysicl_assert_equal_to(_data.size(), _indices.size());
  return _data.size(); }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
void
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::resize(std::size_t s)
{ metaphysicl_assert_equal_to(_data.size(), _indices.size());
  _data.resize(s);
  _indices.resize(s); }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::DynamicSparseNumberBase() {}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <typename Data2, typename Indices2, class... SubTypeArgs2>
inline
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::
DynamicSparseNumberBase(const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & src)
{ this->resize(src.size());
  std::copy(src.nude_data().begin(), src.nude_data().end(), _data.begin());
  std::copy(src.nude_indices().begin(), src.nude_indices().end(), _indices.begin()); }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <typename Data2, typename Indices2, class... SubTypeArgs2>
inline DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::DynamicSparseNumberBase(
    DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> && src)
{
  _data = std::move(src.nude_data());
  _indices = std::move(src.nude_indices());
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
typename Data::value_type*
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::raw_data()
{ return size()?&_data[0]:NULL; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
const typename Data::value_type*
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::raw_data() const
{ return size()?&_data[0]:NULL; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
typename Data::reference
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::raw_at(unsigned int i)
{ return _data[i]; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
typename Data::const_reference
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::raw_at(unsigned int i) const
{ return _data[i]; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
typename Indices::value_type&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::raw_index(unsigned int i)
{ return _indices[i]; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
const typename Indices::value_type &
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::raw_index(unsigned int i) const
{ return _indices[i]; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
const Data&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::nude_data() const
{ return _data; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
Data&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::nude_data()
{ return _data; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
const Indices&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::nude_indices() const
{ return _indices; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
Indices&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::nude_indices()
{ return _indices; }

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
std::size_t
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::runtime_index_query(index_value_type i) const
{
  auto it =
    std::lower_bound(_indices.begin(), _indices.end(), i);
  if (it == _indices.end() || *it != i)
    return std::numeric_limits<std::size_t>::max();
  std::size_t offset = it - _indices.begin();
  metaphysicl_assert_equal_to(_indices[offset], i);
  return offset;
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
std::size_t
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::runtime_index_of(index_value_type i) const
{
  auto it =
    std::lower_bound(_indices.begin(), _indices.end(), i);
  metaphysicl_assert(it != _indices.end());
  std::size_t offset = it - _indices.begin();
  metaphysicl_assert_equal_to(_indices[offset], i);
  return offset;
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
typename Data::value_type &
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator[](index_value_type i)
{
  typedef typename Data::value_type T;
  static T zero = 0;

  // Bad user code could make this fail.  We'd prefer to catch OOB
  // writes at *write* time but at least we can catch at read time.
  metaphysicl_assert(zero == T(0));

  std::size_t rq = runtime_index_query(i);
  if (rq == std::numeric_limits<std::size_t>::max())
    return zero;
  return _data[rq];
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
const typename Data::value_type&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator[](index_value_type i) const
{
  typedef typename Data::value_type T;
  static const T zero = 0;
  std::size_t rq = runtime_index_query(i);
  if (rq == std::numeric_limits<std::size_t>::max())
    return zero;
  return _data[rq];
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <unsigned int i>
inline
typename DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::template entry_type<i>::type&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::get() {
  return _data[runtime_index_of(i)];
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <unsigned int i>
inline
const typename DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::template entry_type<i>::type&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::get() const {
  return _data[runtime_index_of(i)];
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
typename DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::value_type&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::insert(unsigned int i)
{
  auto upper_it =
    std::lower_bound(_indices.begin(), _indices.end(), i);
  std::size_t offset = upper_it - _indices.begin();

  // If we don't have entry i, insert it.  Yes this is O(N).
  if ((upper_it == _indices.end()) ||
      *upper_it != i)
    {
      std::size_t old_size = this->size();
      this->resize(old_size+1);
      std::copy_backward(_indices.begin()+offset, _indices.begin()+old_size, _indices.end());
      std::copy_backward(_data.begin()+offset, _data.begin()+old_size, _data.end());
      _indices[offset] = i;
      _data[offset] = 0;
    }

  // We have entry i now; return it
  return _data[offset];
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <unsigned int i>
inline
typename DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::template entry_type<i>::type&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::insert() {
  return this->insert(i);
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <unsigned int i, typename T2>
inline
void
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::set(const T2& val) {
  _data[runtime_index_of(i)] = val;
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
bool
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::boolean_test() const {
  std::size_t index_size = size();
  for (unsigned int i=0; i != index_size; ++i)
    if (_data[i])
      return true;
  return false;
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
SubType<SubTypeArgs...>
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator- () const {
  std::size_t index_size = size();
  SubType<SubTypeArgs...> returnval;
  returnval.resize(index_size);
  for (unsigned int i=0; i != index_size; ++i)
    {
      returnval.raw_index(i) = _indices[i];
      returnval.raw_at(i) = -_data[i];
    }
  return returnval;
}

  // Since this is a dynamically allocated sparsity pattern, we can
  // increase it as needed to support e.g. operator+=
template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <typename Indices2>
inline
void
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::sparsity_union (const Indices2 & new_indices)
{
  metaphysicl_assert
    (std::adjacent_find(_indices.begin(), _indices.end()) ==
     _indices.end());
  metaphysicl_assert
    (std::adjacent_find(new_indices.begin(), new_indices.end()) ==
     new_indices.end());
#ifdef METAPHYSICL_HAVE_CXX11
  metaphysicl_assert(std::is_sorted(_indices.begin(), _indices.end()));
  metaphysicl_assert(std::is_sorted(new_indices.begin(), new_indices.end()));
#endif

  auto index_it = _indices.begin();
  auto index2_it = new_indices.begin();

  typedef typename Indices2::value_type I2;
  typedef typename CompareTypes<I,I2>::supertype max_index_type;
  max_index_type unseen_indices = 0;

  const I maxI = std::numeric_limits<I>::max();

  while (index2_it != new_indices.end()) {
    I idx1 = (index_it == _indices.end()) ? maxI : *index_it;
    I2 idx2 = *index2_it;

    while (idx1 < idx2) {
      ++index_it;
      idx1 = (index_it == _indices.end()) ? maxI : *index_it;
    }

    while ((idx1 == idx2) &&
           (idx1 != maxI)) {
      ++index_it;
      idx1 = (index_it == _indices.end()) ? maxI : *index_it;
      ++index2_it;
      idx2 = (index2_it == new_indices.end()) ? maxI : *index2_it;
    }

    while (idx2 < idx1) {
      ++unseen_indices;
        ++index2_it;
      if (index2_it == new_indices.end())
        break;
      idx2 = *index2_it;
    }
  }

  // The common case is cheap
  if (!unseen_indices)
    return;

  std::size_t old_size = this->size();

  this->resize(old_size + unseen_indices);

  auto md_it = _data.rbegin();
  auto mi_it = _indices.rbegin();

  auto d_it =
    _data.rbegin() + unseen_indices;
  auto i_it =
    _indices.rbegin() + unseen_indices;
  auto i2_it = new_indices.rbegin();

  // Duplicate copies of rend() to work around
  // http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#179
  auto      mirend  = _indices.rend();
  auto  rend  = mirend;
  auto rend2 = new_indices.rend();
#ifndef NDEBUG
  auto      mdrend = _data.rend();
  auto drend = mdrend;
#endif

  for (; mi_it != mirend; ++md_it, ++mi_it) {
    if ((i_it == rend) ||
        ((i2_it != rend2) &&
         (*i2_it > *i_it))) {
      *md_it = 0;
      *mi_it = *i2_it;
      ++i2_it;
    } else {
      if ((i2_it != rend2) &&
          (*i2_it == *i_it))
        ++i2_it;
      metaphysicl_assert(d_it < drend);
      metaphysicl_assert(md_it < mdrend);
      *md_it = *d_it;
      *mi_it = *i_it;
      ++d_it;
      ++i_it;
    }
  }

  metaphysicl_assert(i_it  == rend);
  metaphysicl_assert(i2_it == rend2);
  metaphysicl_assert(d_it  == drend);
  metaphysicl_assert(md_it == mdrend);
}


  // Since this is a dynamically allocated sparsity pattern, we can
  // decrease it when possible for efficiency
template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <typename Indices2>
inline
void
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::
sparsity_intersection (const Indices2 & new_indices)
{
  metaphysicl_assert
    (std::adjacent_find(_indices.begin(), _indices.end()) ==
     _indices.end());
  metaphysicl_assert
    (std::adjacent_find(new_indices.begin(), new_indices.end()) ==
     new_indices.end());
#ifdef METAPHYSICL_HAVE_CXX11
  metaphysicl_assert(std::is_sorted(_indices.begin(), _indices.end()));
  metaphysicl_assert(std::is_sorted(new_indices.begin(), new_indices.end()));
#endif

#ifndef NDEBUG
  typedef typename Indices2::value_type I2;
  typedef typename CompareTypes<I,I2>::supertype max_index_type;
  auto index_it = _indices.begin();
  auto index2_it = new_indices.begin();

  max_index_type shared_indices = 0;

  const I maxI = std::numeric_limits<I>::max();

  while (index2_it != new_indices.end()) {
    I idx1 = (index_it == _indices.end()) ? maxI : *index_it;
    I2 idx2 = *index2_it;

    while (idx1 < idx2) {
      ++index_it;
      idx1 = (index_it == _indices.end()) ? maxI : *index_it;
    }

    while ((idx1 == idx2) &&
           (idx1 != maxI)) {
      ++index_it;
      idx1 = (index_it == _indices.end()) ? maxI : *index_it;
      ++index2_it;
      idx2 = (index2_it == new_indices.end()) ? maxI : *index2_it;
      ++shared_indices;
    }

    while (idx2 < idx1) {
      ++index2_it;
      if (index2_it == new_indices.end())
        break;
      idx2 = *index2_it;
    }
  }
#endif

  // We'll loop up through the arrays, copying indices (and
  // corresponding data) that should be there downward into place.

  // Merged values:
  auto md_it = _data.begin();
  auto mi_it = _indices.begin();

  // Our old values:
  auto d_it = _data.begin();
  auto i_it = _indices.begin();

  // Values to merge with:
  auto i2_it = new_indices.begin();

  for (; i_it != _indices.end() && i2_it != new_indices.end();
       ++md_it, ++mi_it, ++d_it, ++i_it, ++i2_it) {
    while (*i2_it < *i_it) {
      ++i2_it;
      if (i2_it == new_indices.end())
        break;
    }
    if (i2_it == new_indices.end())
      break;
    while (*i2_it > *i_it) {
        ++i_it;
      if (i_it == _indices.end())
        break;
    }
    if (i_it == _indices.end())
      break;

    *md_it = *d_it;
    *mi_it = *i_it;
  }

  metaphysicl_assert_equal_to(md_it - _data.begin(),
                              shared_indices);
  metaphysicl_assert_equal_to(mi_it - _indices.begin(),
                              shared_indices);

  const std::size_t n_indices = md_it - _data.begin();

  _indices.resize(n_indices);
  _data.resize(n_indices);
}



  // Since this is a dynamically allocated sparsity pattern, we can
  // decrease it when possible for efficiency
template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
inline
void
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::sparsity_trim (const value_type tolerance)
{
  metaphysicl_assert
    (std::adjacent_find(_indices.begin(), _indices.end()) ==
     _indices.end());
#ifdef METAPHYSICL_HAVE_CXX11
  metaphysicl_assert(std::is_sorted(_indices.begin(), _indices.end()));
#endif
  metaphysicl_assert(tolerance >= 0);

#ifndef NDEBUG
  I used_indices = 0;

  {
    auto index_it = _indices.begin();
    auto data_it = _data.begin();
    for (; index_it != _indices.end(); ++index_it, ++data_it)
      if (std::abs(*data_it) > tolerance)
        ++used_indices;
  }
#endif

  // We'll loop up through the arrays, copying indices (and
  // corresponding data) that should be there downward into place.

  // Downward-merged values:
  auto md_it = _data.begin();
  auto mi_it = _indices.begin();

  // Our old values:
  auto d_it = _data.begin();

  for (auto i_it = _indices.begin();
       i_it != _indices.end(); ++i_it, ++d_it)
    if (std::abs(*d_it) > tolerance)
      {
        *mi_it = *i_it;
        *md_it = *d_it;
        ++mi_it;
        ++md_it;
      }

  const std::size_t n_indices = md_it - _data.begin();

  metaphysicl_assert_equal_to(n_indices, used_indices);
  metaphysicl_assert_equal_to(mi_it - _indices.begin(),
                              used_indices);

  _indices.resize(n_indices);
  _data.resize(n_indices);
}

  // Not defineable since !0 != 0
  // SubType<SubTypeArgs...> operator! () const;

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <class... SubTypeArgs2>
inline
SubType<SubTypeArgs...>&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator+= (const SubType<SubTypeArgs2...>& a)
{
  // Resize if necessary
  this->sparsity_union(a.nude_indices());

  auto data_it  = _data.begin();
  auto index_it = _indices.begin();
  auto data2_it  =
    a.nude_data().begin();
  auto index2_it =
    a.nude_indices().begin();
  for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
    {
      auto idx1 = *index_it;
      auto idx2 = *index2_it;

      while (idx1 < idx2) {
        ++index_it;
        ++data_it;
        metaphysicl_assert(index_it != _indices.end());
        idx1 = *index_it;
      }
      metaphysicl_assert_equal_to(idx1, idx2);

      *data_it += *data2_it;
    }

  return static_cast<SubType<SubTypeArgs...>&>(*this);
}


template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <class... SubTypeArgs2>
inline
SubType<SubTypeArgs...>&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator-= (const SubType<SubTypeArgs2...>& a)
{
  // Resize if necessary
  this->sparsity_union(a.nude_indices());

  auto data_it  = _data.begin();
  auto index_it = _indices.begin();
  auto data2_it  =
    a.nude_data().begin();
  auto index2_it =
    a.nude_indices().begin();
  for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
    {
      auto idx1 = *index_it;
      auto idx2 = *index2_it;

      while (idx1 < idx2) {
        ++index_it;
        ++data_it;
        metaphysicl_assert(index_it != _indices.end());
        idx1 = *index_it;
      }
      metaphysicl_assert_equal_to(idx1, idx2);

      *data_it -= *data2_it;
    }

  return static_cast<SubType<SubTypeArgs...>&>(*this);
}


template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <class... SubTypeArgs2>
inline
SubType<SubTypeArgs...>&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator*= (const SubType<SubTypeArgs2...>& a)
{
  // Resize if possible
  this->sparsity_intersection(a.nude_indices());

  auto data_it  = _data.begin();
  auto index_it = _indices.begin();
  auto data2_it  =
    a.nude_data().begin();
  auto index2_it =
    a.nude_indices().begin();
  for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
    {
      auto idx1 = *index_it;
      auto idx2 = *index2_it;

      while (idx1 < idx2) {
        ++index_it;
        ++data_it;
        metaphysicl_assert(index_it != _indices.end());
        idx1 = *index_it;
      }

      if (idx1 == idx2)
        *data_it *= *data2_it;
    }

  return static_cast<SubType<SubTypeArgs...>&>(*this);
}


template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <class... SubTypeArgs2>
inline
SubType<SubTypeArgs...>&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator/= (const SubType<SubTypeArgs2...>& a)
{
  auto data_it  = _data.begin();
  auto index_it = _indices.begin();
  auto data2_it  =
    a.nude_data().begin();
  auto index2_it =
    a.nude_indices().begin();
  for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
    {
      auto idx1 = *index_it;
      auto idx2 = *index2_it;

      while (idx1 < idx2) {
        ++index_it;
        ++data_it;
        metaphysicl_assert(index_it != _indices.end());
        idx1 = *index_it;
      }

      if (idx1 == idx2)
        *data_it /= *data2_it;
    }

  return static_cast<SubType<SubTypeArgs...>&>(*this);
}


template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <typename T2>
inline
SubType<SubTypeArgs...>&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator*= (const T2& a)
{
  std::size_t index_size = size();
  for (unsigned int i=0; i != index_size; ++i)
    _data[i] *= a;
  return static_cast<SubType<SubTypeArgs...>&>(*this);
}

template <typename Data, typename Indices, template <class...> class SubType, class... SubTypeArgs>
template <typename T2>
inline
SubType<SubTypeArgs...>&
DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>::operator/= (const T2& a)
{
  std::size_t index_size = size();
  for (unsigned int i=0; i != index_size; ++i)
    _data[i] /= a;
  return static_cast<SubType<SubTypeArgs...>&>(*this);
}

//
// Non-member functions
//

template <template <class...> class SubType,
          typename BoolData,
          typename BoolIndices,
          class... BoolSubTypeArgs,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename Data2,
          typename Indices2,
          class... SubTypeArgs2>
inline typename SubType<SubTypeArgs...>::template rebind<
    typename CompareTypes<typename Data::value_type, typename Data2::value_type>::supertype,
    typename CompareTypes<typename Indices::value_type,
                          typename Indices2::value_type>::supertype>::other
if_else(
    const DynamicSparseNumberBase<BoolData, BoolIndices, SubType, BoolSubTypeArgs...> & condition,
    const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & if_true,
    const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & if_false)
{
  metaphysicl_assert
    (std::adjacent_find(condition.nude_indices().begin(), condition.nude_indices().end()) ==
     condition.nude_indices().end());
  metaphysicl_assert
    (std::adjacent_find(if_true.nude_indices().begin(), if_true.nude_indices().end()) ==
     if_true.nude_indices().end());
  metaphysicl_assert
    (std::adjacent_find(if_false.nude_indices().begin(), if_false.nude_indices().end()) ==
     if_false.nude_indices().end());
#ifdef METAPHYSICL_HAVE_CXX11
  metaphysicl_assert(std::is_sorted(condition.nude_indices().begin(), condition.nude_indices().end()));
  metaphysicl_assert(std::is_sorted(if_true.nude_indices().begin(), if_true.nude_indices().end()));
  metaphysicl_assert(std::is_sorted(if_false.nude_indices().begin(), if_false.nude_indices().end()));
#endif

  typedef
      typename CompareTypes<typename Data::value_type, typename Data2::value_type>::supertype TS;
  typedef
      typename CompareTypes<typename Indices::value_type, typename Indices2::value_type>::supertype
          IS;

  typename SubType<SubTypeArgs...>::template rebind<TS, IS>::other returnval;

  // First count returnval size
  IS required_size = 0;
  {
    auto indexcond_it      = condition.nude_indices().begin();
    auto datacond_it        = condition.nude_data().begin();
    auto indextrue_it       = if_true.nude_indices().begin();
    const auto endtrue_it   = if_true.nude_indices().end();
    auto datatrue_it        = if_true.nude_data().begin();
    auto indexfalse_it     = if_false.nude_indices().begin();
    const auto endfalse_it = if_false.nude_indices().end();
    auto datafalse_it      = if_false.nude_data().begin();

    for (; indexcond_it != condition.nude_indices().end(); ++indexcond_it, ++datacond_it)
     {
       while (indexfalse_it != endfalse_it &&
              *indexfalse_it < *indexcond_it)
         {
           if (*datafalse_it)
             ++required_size;

           ++indexfalse_it;
           ++datafalse_it;
         }

       if (*datacond_it)
         {
           while (indextrue_it != endtrue_it &&
                  *indextrue_it < *indexcond_it)
             {
               ++indextrue_it;
               ++datatrue_it;
             }
           if (indextrue_it != endtrue_it &&
               *indextrue_it == *indexcond_it &&
               *datatrue_it)
             {
               ++required_size;
               ++indextrue_it;
               ++datatrue_it;
             }
           if (*indexfalse_it == *indexcond_it)
             {
               ++indexfalse_it;
               ++datafalse_it;
             }
         }
       else
         {
           if (indexfalse_it != endfalse_it &&
               *indexfalse_it == *indexcond_it &&
               *datafalse_it)
             {
               ++required_size;
               ++indexfalse_it;
               ++datafalse_it;
             }
         }
     }
  }

  // Then fill returnval
  returnval.resize(required_size);
  {
    auto indexcond_it      = condition.nude_indices().begin();
    auto datacond_it        = condition.nude_data().begin();
    auto indextrue_it       = if_true.nude_indices().begin();
    const auto endtrue_it   = if_true.nude_indices().end();
    auto datatrue_it        = if_true.nude_data().begin();
    auto indexfalse_it     = if_false.nude_indices().begin();
    const auto endfalse_it = if_false.nude_indices().end();
    auto datafalse_it      = if_false.nude_data().begin();

    auto indexreturn_it          = returnval.nude_indices().begin();
    auto datareturn_it           = returnval.nude_data().begin();

    for (; indexcond_it != condition.nude_indices().end(); ++indexcond_it, ++datacond_it)
     {
       while (indexfalse_it != endfalse_it &&
              *indexfalse_it < *indexcond_it)
         {
           if (*datafalse_it)
             {
               *indexreturn_it = *indexfalse_it;
               *datareturn_it  = *datafalse_it;
               ++indexreturn_it;
               ++datareturn_it;
             }

           ++indexfalse_it;
           ++datafalse_it;
         }

       if (*datacond_it)
         {
           while (indextrue_it != endtrue_it &&
                  *indextrue_it < *indexcond_it)
             {
               ++indextrue_it;
               ++datatrue_it;
             }
           if (indextrue_it != endtrue_it &&
               *indextrue_it == *indexcond_it &&
               *datatrue_it)
             {
               *indexreturn_it = *indextrue_it;
               *datareturn_it  = *datatrue_it;
               ++indexreturn_it;
               ++datareturn_it;
               ++indextrue_it;
               ++datatrue_it;
             }
           if (*indexfalse_it == *indexcond_it)
             {
               ++indexfalse_it;
               ++datafalse_it;
             }
         }
       else
         {
           if (indexfalse_it != endfalse_it &&
               *indexfalse_it == *indexcond_it &&
               *datafalse_it)
             {
               *indexreturn_it = *indexfalse_it;
               *datareturn_it  = *datafalse_it;
               ++indexreturn_it;
               ++datareturn_it;
               ++indexfalse_it;
               ++datafalse_it;
             }
         }
     }
  }

  metaphysicl_assert
    (std::adjacent_find(returnval.nude_indices().begin(), returnval.nude_indices().end()) ==
     returnval.nude_indices().end());
#ifdef METAPHYSICL_HAVE_CXX11
  metaphysicl_assert(std::is_sorted(returnval.nude_indices().begin(), returnval.nude_indices().end()));
#endif

  return returnval;
}



#define DynamicSparseNumberBase_op_ab(opname, subtypename, functorname) \
  template <class... AArgs, class... BArgs> \
inline \
typename Symmetric##functorname##Type<subtypename<AArgs...>, \
                                      subtypename<BArgs...>>::supertype \
operator opname (const subtypename<AArgs...> & a, const subtypename<BArgs...> & b) \
{ \
  typedef typename Symmetric##functorname##Type<subtypename<AArgs...>,  \
                                                subtypename<BArgs...>>::supertype type; \
  type returnval = a; \
  returnval opname##= b; \
  return returnval; \
}


#if __cplusplus >= 201103L

#define DynamicSparseNumberBase_op(subtypename, opname, functorname) \
DynamicSparseNumberBase_op_ab(opname, subtypename, functorname) \
 \
template <class... AArgs, class... BArgs> \
inline typename Symmetric##functorname##Type<subtypename<AArgs...>, \
                                             subtypename<BArgs...>>::supertype \
operator opname(subtypename<AArgs...> && a, const subtypename<BArgs...> & b) \
{ \
  typedef typename Symmetric##functorname##Type<subtypename<AArgs...>, \
                                                subtypename<BArgs...>>::supertype type; \
  type returnval = std::move(a); \
  returnval opname##= b; \
  return returnval; \
}

#else

#define DynamicSparseNumberBase_op(subtypename, opname, functorname) \
DynamicSparseNumberBase_op_ab(opname, subtypename, functorname)

#endif

// Let's also allow scalar times vector.
// Scalar plus vector, etc. remain undefined in the sparse context.

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename MultipliesType<SubType<SubTypeArgs...>, T, true>::supertype
operator*(const T & a, const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & b)
{
  const unsigned int index_size = b.size();

  typename MultipliesType<SubType<SubTypeArgs...>,T,true>::supertype
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a * b.raw_at(i);
    returnval.raw_index(i) = b.raw_index(i);
  }

  return returnval;
}

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename MultipliesType<SubType<SubTypeArgs...>, T>::supertype
operator*(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, const T & b)
{
  const unsigned int index_size = a.size();

  typename MultipliesType<SubType<SubTypeArgs...>,T>::supertype
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i) * b;
    returnval.raw_index(i) = a.raw_index(i);
  }
  return returnval;
}

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename DividesType<SubType<SubTypeArgs...>, T>::supertype
operator/(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, const T & b)
{
  const unsigned int index_size = a.size();

  typename DividesType<SubType<SubTypeArgs...>,T>::supertype returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i) / b;
    returnval.raw_index(i) = a.raw_index(i);
  }

  return returnval;
}

#if __cplusplus >= 201103L
template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename MultipliesType<SubType<SubTypeArgs...>, T>::supertype
operator*(DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> && a, const T & b)
{
  typename MultipliesType<SubType<SubTypeArgs...>,T>::supertype
    returnval = std::move(static_cast<SubType<SubTypeArgs...>&&>(a));

  returnval *= b;

  return returnval;
}

template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T>
inline typename DividesType<SubType<SubTypeArgs...>, T>::supertype
operator/(DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> && a, const T & b)
{
  typename DividesType<SubType<SubTypeArgs...>,T>::supertype
    returnval = std::move(static_cast<SubType<SubTypeArgs...>&&>(a));

  returnval /= b;

  return returnval;
}
#endif


#define DynamicSparseNumberBase_operator_binary(opname, functorname) \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename Data2, \
          typename Indices2, \
          class... SubTypeArgs2> \
inline typename SubType<SubTypeArgs...>::template rebind< \
  bool, \
  typename CompareTypes<typename Indices::value_type, \
                        typename Indices2::value_type>::supertype>::other \
operator opname(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, \
                const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & b) \
{ \
  typedef typename Data::value_type T; \
  typedef typename Data2::value_type T2; \
  typedef typename Indices::value_type I; \
  typedef typename Indices2::value_type I2; \
  typedef typename SymmetricCompareTypes<T, T2>::supertype TS; \
  typedef typename CompareTypes<I, I2>::supertype IS; \
  typedef typename SubType<SubTypeArgs...>::template rebind< \
    bool, \
    typename CompareTypes<typename Indices::value_type, \
                          typename Indices2::value_type>::supertype>::other RetType; \
  RetType returnval; \
  returnval.nude_indices() = a.nude_indices(); \
  returnval.nude_data().resize(a.nude_indices().size()); \
  returnval.sparsity_union(b.nude_indices()); \
 \
  auto  index_a_it = a.nude_indices().begin(); \
  auto index_b_it = b.nude_indices().begin(); \
  auto     index_out_it = returnval.nude_indices().begin(); \
 \
  auto  data_a_it = a.nude_data().begin(); \
  auto data_b_it = b.nude_data().begin(); \
  auto     data_out_it = returnval.nude_data().begin(); \
 \
  const IS  maxIS  = std::numeric_limits<IS>::max(); \
 \
  for (; index_out_it != returnval.nude_indices().end(); ++index_out_it, ++data_out_it) { \
    const IS index_a = (index_a_it == a.nude_indices().end()) ? maxIS : *index_a_it; \
    const IS index_b = (index_b_it == b.nude_indices().end()) ? maxIS : *index_b_it; \
    const IS index_out = *index_out_it; \
    const TS data_a  = (index_a_it == a.nude_indices().end()) ? 0: *data_a_it; \
    const TS data_b  = (index_b_it == b.nude_indices().end()) ? 0: *data_b_it; \
    TS & data_out = *data_out_it; \
 \
    if (index_a == index_out) { \
      if (index_b == index_out) { \
        data_out = data_a opname data_b; \
        index_b_it++; \
        data_b_it++; \
      } else { \
        data_out = data_a opname 0; \
      } \
      index_a_it++; \
      data_a_it++; \
    } else { \
      metaphysicl_assert_equal_to(index_b, index_out); \
      data_out = 0 opname data_b; \
      index_b_it++; \
      data_b_it++; \
    } \
  } \
 \
  return returnval; \
} \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename T> \
inline typename SubType<SubTypeArgs...>::template rebind<bool>::other operator opname( \
  const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, const T & b) \
{ \
  typename SubType<SubTypeArgs...>::template rebind<bool>::other returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices() = a.nude_indices(); \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = (a.raw_at(i) opname b); \
 \
  return returnval; \
} \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename T> \
inline typename SubType<SubTypeArgs...>::template rebind<bool>::other operator opname( \
    const T & a, const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & b) \
{ \
  typename SubType<SubTypeArgs...>::template rebind<bool>::other returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.nude_indices() = a.nude_indices(); \
  returnval.nude_data().resize(index_size); \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = (a opname b.raw_at(i)); \
 \
  return returnval; \
}

// NOTE: unary functions for which 0-op-0 is true are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

DynamicSparseNumberBase_operator_binary(<, less)
// DynamicSparseNumberBase_operator_binary(<=)
DynamicSparseNumberBase_operator_binary(>, greater)
// DynamicSparseNumberBase_operator_binary(>=)
// DynamicSparseNumberBase_operator_binary(==)
DynamicSparseNumberBase_operator_binary(!=, not_equal_to)

// FIXME - make && an intersection rather than a union for efficiency
DynamicSparseNumberBase_operator_binary(&&, logical_and)
DynamicSparseNumberBase_operator_binary(||, logical_or)


template <template <class...> class SubType, typename Data, typename Indices, class... SubTypeArgs>
inline
std::ostream&
operator<< (std::ostream& output, const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...>& a)
{
  // Enclose the entire output in braces
  output << '{';

  std::size_t index_size = a.size();

  // Output the first value from a non-empty set
  // All values are given as ordered (index, value) pairs
  if (index_size)
    output << '(' << a.raw_index(0) << ',' <<
              a.raw_at(0) << ')';

  // Output the comma-separated subsequent values from a non-singleton
  // set
  for (unsigned int i = 1; i < index_size; ++i)
    {
      output << ", (" << a.raw_index(i) << ',' << a.raw_data()[i] << ')';
    }
  output << '}';
  return output;
}

} // namespace MetaPhysicL

namespace std {

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DynamicSparseNumberBase;
using MetaPhysicL::SymmetricCompareTypes;

#define DynamicSparseNumberBase_std_unary(funcname) \
template <template <class...> class SubType, \
          typename Data, typename Indices, class... SubTypeArgs> \
inline \
SubType<SubTypeArgs...> \
funcname (const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a) \
{ \
  std::size_t index_size = a.size(); \
  SubType<SubTypeArgs...> returnval; \
  returnval.nude_indices() = a.nude_indices(); \
  returnval.nude_data().resize(index_size); \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = std::funcname(a.raw_at(i)); \
 \
  return returnval; \
}


#define DynamicSparseNumberBase_fl_unary(funcname) \
DynamicSparseNumberBase_std_unary(funcname##f) \
DynamicSparseNumberBase_std_unary(funcname##l)


#define DynamicSparseNumberBase_stdfl_unary(funcname) \
DynamicSparseNumberBase_std_unary(funcname) \
DynamicSparseNumberBase_fl_unary(funcname)


#define DynamicSparseNumberBase_std_binary_union(funcname) \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename Data2, \
          typename Indices2, \
          class... SubTypeArgs2> \
inline typename SubType<SubTypeArgs...>::template rebind< \
    typename SymmetricCompareTypes<typename Data::value_type, \
                                   typename Data2::value_type>::supertype, \
    typename CompareTypes<typename Indices::value_type, \
                          typename Indices2::value_type>::supertype>::other \
funcname(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, \
         const DynamicSparseNumberBase<Data2, Indices2, SubType, SubTypeArgs2...> & b) \
{ \
  typedef typename Data::value_type T; \
  typedef typename Data2::value_type T2; \
  typedef typename Indices::value_type I; \
  typedef typename Indices2::value_type I2; \
  typedef typename SymmetricCompareTypes<T, T2>::supertype TS; \
  typedef typename CompareTypes<I, I2>::supertype IS; \
  typedef typename SubType<SubTypeArgs...>::template rebind< \
      typename SymmetricCompareTypes<typename Data::value_type, \
                                     typename Data2::value_type>::supertype, \
      typename CompareTypes<typename Indices::value_type, \
                            typename Indices2::value_type>::supertype>::other RetType; \
  RetType returnval; \
 \
  std::size_t index_size = a.nude_indices.size(); \
  returnval.nude_indices = a.nude_indices; \
  returnval.nude_data.resize(index_size); \
  returnval.sparsity_union(b.nude_indices); \
 \
  auto  index_a_it = a.nude_indices.begin(); \
  auto index_b_it = b.nude_indices.begin(); \
  auto     index_out_it = returnval.nude_indices.begin(); \
 \
  auto  data_a_it = a.nude_data.begin(); \
  auto data_b_it = b.nude_data.begin(); \
  auto     data_out_it = returnval.nude_data.begin(); \
 \
  const IS  maxIS  = std::numeric_limits<IS>::max(); \
 \
  for (; index_out_it != returnval.nude_indices.end(); ++index_out_it, ++data_out_it) { \
    const IS index_a = (index_a_it == a.nude_indices.end()) ? maxIS : *index_a_it; \
    const IS index_b = (index_b_it == b.nude_indices.end()) ? maxIS : *index_b_it; \
    const IS index_out = *index_out_it; \
    const TS data_a  = (index_a_it == a.nude_indices.end()) ? 0: *data_a_it; \
    const TS data_b  = (index_b_it == b.nude_indices.end()) ? 0: *data_b_it; \
    TS & data_out = *data_out_it; \
 \
    if (index_a == index_out) { \
      if (index_b == index_out) { \
        data_out = std::funcname(data_a, data_b); \
        index_b_it++; \
        data_b_it++; \
      } else { \
        data_out = std::funcname(data_a, 0); \
      } \
      index_a_it++; \
      data_a_it++; \
    } else { \
      metaphysicl_assert_equal_to(index_b, index_out); \
      data_out = std::funcname(0, data_b); \
      index_b_it++; \
      data_b_it++; \
    } \
  } \
 \
  return returnval; \
} \
 \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename T2> \
inline typename SubType<SubTypeArgs...>::template rebind< \
    typename SymmetricCompareTypes<typename Data::value_type, T2>::supertype, \
    typename Indices::value_type>::other \
funcname(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, \
         const T2 & b) \
{ \
  typedef typename Data::value_type T; \
  typedef typename SymmetricCompareTypes<T, T2>::supertype TS; \
  typedef \
      typename SubType<SubTypeArgs...>::template rebind<TS, typename Indices::value_type>::other \
          RetType; \
  RetType returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices = a.nude_indices; \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = std::funcname(a.raw_at(i), b); \
 \
  return returnval; \
} \
 \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs, \
          typename T> \
inline typename SubType<SubTypeArgs...>::template rebind< \
    typename SymmetricCompareTypes<T, typename Data::value_type>::supertype, \
    typename Indices::value_type>::other \
funcname(const T & a, const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & b) \
{ \
  typedef typename Data::value_type T2; \
  typedef typename SymmetricCompareTypes<T, T2>::supertype TS; \
  typedef \
      typename SubType<SubTypeArgs...>::template rebind<TS, typename Indices::value_type>::other \
          RetType; \
  RetType returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices = b.nude_indices; \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = std::funcname(a, b.raw_at(i)); \
 \
  return returnval; \
}


#define DynamicSparseNumberBase_fl_binary_union(funcname) \
DynamicSparseNumberBase_std_binary_union(funcname##f) \
DynamicSparseNumberBase_std_binary_union(funcname##l)


#define DynamicSparseNumberBase_stdfl_binary_union(funcname) \
DynamicSparseNumberBase_std_binary_union(funcname) \
DynamicSparseNumberBase_fl_binary_union(funcname)


// Pow needs its own specialization, both to avoid being confused by
// pow<T1,T2> and because pow(x,0) isn't 0.
template <template <class...> class SubType,
          typename Data,
          typename Indices,
          class... SubTypeArgs,
          typename T2>
inline typename SubType<SubTypeArgs...>::template rebind<
    typename SymmetricCompareTypes<typename Data::value_type, T2>::supertype,
    typename Indices::value_type>::other
pow(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & a, const T2 & b)
{
  typedef typename Data::value_type T;
  typedef typename SymmetricCompareTypes<T, T2>::supertype TS;
  typedef typename SubType<SubTypeArgs...>::template rebind<TS, typename Indices::value_type>::other
      RetType;
  RetType returnval;

  std::size_t index_size = a.size();
  returnval.nude_indices() = a.nude_indices();
  returnval.nude_data().resize(index_size);

  for (unsigned int i=0; i != index_size; ++i)
    returnval.raw_at(i) = std::pow(a.raw_at(i), b);

  return returnval;
}


// NOTE: unary functions for which f(0) != 0 are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

// DynamicSparseNumberBase_std_binary(pow) // separate definition
// DynamicSparseNumberBase_std_unary(exp)
// DynamicSparseNumberBase_std_unary(log)
// DynamicSparseNumberBase_std_unary(log10)
DynamicSparseNumberBase_std_unary(sin)
// DynamicSparseNumberBase_std_unary(cos)
DynamicSparseNumberBase_std_unary(tan)
DynamicSparseNumberBase_std_unary(asin)
// DynamicSparseNumberBase_std_unary(acos)
DynamicSparseNumberBase_std_unary(atan)
DynamicSparseNumberBase_std_binary_union(atan2)
DynamicSparseNumberBase_std_unary(sinh)
// DynamicSparseNumberBase_std_unary(cosh)
DynamicSparseNumberBase_std_unary(tanh)
DynamicSparseNumberBase_std_unary(sqrt)
DynamicSparseNumberBase_std_unary(abs)
DynamicSparseNumberBase_std_unary(fabs)
DynamicSparseNumberBase_std_binary_union(max)
DynamicSparseNumberBase_std_binary_union(min)
DynamicSparseNumberBase_std_unary(ceil)
DynamicSparseNumberBase_std_unary(floor)
DynamicSparseNumberBase_std_binary_union(fmod) // TODO: optimize this

#if __cplusplus >= 201103L
DynamicSparseNumberBase_std_unary(llabs)
DynamicSparseNumberBase_std_unary(imaxabs)
DynamicSparseNumberBase_fl_unary(fabs)
DynamicSparseNumberBase_stdfl_unary(expm1)
DynamicSparseNumberBase_fl_unary(sqrt)
DynamicSparseNumberBase_stdfl_unary(cbrt)
DynamicSparseNumberBase_fl_unary(sin)
DynamicSparseNumberBase_fl_unary(tan)
DynamicSparseNumberBase_fl_unary(asin)
DynamicSparseNumberBase_fl_unary(atan)
DynamicSparseNumberBase_stdfl_unary(asinh)
DynamicSparseNumberBase_stdfl_unary(atanh)
DynamicSparseNumberBase_stdfl_unary(erf)
DynamicSparseNumberBase_fl_unary(ceil)
DynamicSparseNumberBase_fl_unary(floor)
DynamicSparseNumberBase_stdfl_unary(trunc)
DynamicSparseNumberBase_stdfl_unary(round)
DynamicSparseNumberBase_stdfl_unary(nearbyint)
DynamicSparseNumberBase_stdfl_unary(rint)

DynamicSparseNumberBase_fl_binary_union(fmod)
DynamicSparseNumberBase_stdfl_binary_union(remainder) // TODO: optimize this
DynamicSparseNumberBase_stdfl_binary_union(fmax)
DynamicSparseNumberBase_stdfl_binary_union(fmin)
DynamicSparseNumberBase_stdfl_binary_union(fdim)
DynamicSparseNumberBase_stdfl_binary_union(hypot)
DynamicSparseNumberBase_fl_binary_union(atan2)
#endif // __cplusplus >= 201103L

#define DynamicSparseNumberBase_std_unary_complex(funcname) \
template <template <class...> class SubType, \
          typename Data, \
          typename Indices, \
          class... SubTypeArgs> \
inline auto funcname(const DynamicSparseNumberBase<Data, Indices, SubType, SubTypeArgs...> & in) \
    ->typename SubType<SubTypeArgs...>::template rebind<decltype(std::funcname( \
                                                            typename Data::value_type())), \
                                                        typename Indices::value_type>::other \
{ \
  typedef typename SubType<SubTypeArgs...>::template rebind< \
      decltype(std::funcname(typename Data::value_type())), \
      typename Indices::value_type>::other RetType; \
  RetType returnval; \
  auto size = in.size(); \
  returnval.nude_indices() = in.nude_indices(); \
  returnval.nude_data().resize(size); \
 \
  for (decltype(size) i = 0; i < size; ++i) \
    returnval.raw_at(i) = std::funcname(in.raw_at(i));  \
  return returnval; \
}

DynamicSparseNumberBase_std_unary_complex(real)
DynamicSparseNumberBase_std_unary_complex(imag)
DynamicSparseNumberBase_std_unary_complex(norm)
} // namespace std


#endif // METAPHYSICL_DYNAMICSPARSENUMBERARRAY_H
