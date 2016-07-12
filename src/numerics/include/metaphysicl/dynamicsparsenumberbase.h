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

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <ostream>
#include <vector>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/ct_set.h"
#include "metaphysicl/metaphysicl_asserts.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/sparsenumberutils.h"
#include "metaphysicl/testable.h"

namespace MetaPhysicL {

// Data type T, index type I
template <typename T, typename I, template <typename, typename> class SubType>
class DynamicSparseNumberBase
{
public:
  typedef T value_type;

  template <unsigned int i>
  struct entry_type {
    typedef value_type type;
  };

  typedef I index_value_type;

  std::size_t size() const
    { metaphysicl_assert_equal_to(_data.size(), _indices.size());
      return _data.size(); }

  void resize(std::size_t s)
    { metaphysicl_assert_equal_to(_data.size(), _indices.size());
      _data.resize(s);
      _indices.resize(s); }

  DynamicSparseNumberBase() {}

#if __cplusplus >= 201103L
  // Move constructors are useful when all your data is on the heap
  DynamicSparseNumberBase(DynamicSparseNumberBase<T, I, SubType> && src) = default;

  // Move assignment avoids heap operations too
  DynamicSparseNumberBase& operator= (DynamicSparseNumberBase<T, I, SubType> && src) = default;

  // Standard copy operations get implicitly deleted upon move
  // constructor definition, so we manually enable them.
  DynamicSparseNumberBase(const DynamicSparseNumberBase<T, I, SubType> & src) = default;

  DynamicSparseNumberBase& operator= (const DynamicSparseNumberBase<T, I, SubType> & src) = default;
#endif

  template <typename T2, typename I2>
  DynamicSparseNumberBase(const DynamicSparseNumberBase<T2, I2, SubType> & src)
    { this->resize(src.size());
      std::copy(src.nude_data().begin(), src.nude_data().end(), _data.begin());
      std::copy(src.nude_indices().begin(), src.nude_indices().end(), _indices.begin()); }

  T* raw_data()
    { return size()?&_data[0]:NULL; }

  const T* raw_data() const
    { return size()?&_data[0]:NULL; }

  typename std::vector<T>::reference raw_at(unsigned int i)
    { return _data[i]; }

  typename std::vector<T>::const_reference raw_at(unsigned int i) const
    { return _data[i]; }

  I& raw_index(unsigned int i)
    { return _indices[i]; }

  const I& raw_index(unsigned int i) const
    { return _indices[i]; }

  // FIXME: these encapsulation violations are necessary for std::pow
  // until I can figure out the right friend declaration.
  const std::vector<T>& nude_data() const
    { return _data; }

  std::vector<T>& nude_data()
    { return _data; }

  const std::vector<I>& nude_indices() const
    { return _indices; }

  std::vector<I>& nude_indices()
    { return _indices; }

  std::size_t runtime_index_query(index_value_type i) const
    {
      typename std::vector<I>::const_iterator it =
        std::lower_bound(_indices.begin(), _indices.end(), i);
      if (it == _indices.end())
        return std::numeric_limits<std::size_t>::max();
      std::size_t offset = it - _indices.begin();
      metaphysicl_assert_equal_to(_indices[offset], i);
      return offset;
    }

  std::size_t runtime_index_of(index_value_type i) const
    {
      typename std::vector<I>::const_iterator it =
        std::lower_bound(_indices.begin(), _indices.end(), i);
      metaphysicl_assert(it != _indices.end());
      std::size_t offset = it - _indices.begin();
      metaphysicl_assert_equal_to(_indices[offset], i);
      return offset;
    }

  T& operator[](index_value_type i)
    { return _data[runtime_index_of(i)]; }

  const T& operator[](index_value_type i) const
    { return _data[runtime_index_of(i)]; }

  T query(index_value_type i) const
    {
      std::size_t rq = runtime_index_query(i);
      if (rq == std::numeric_limits<std::size_t>::max())
        return 0;
      return _data[rq];
    }

  template <unsigned int i>
  typename entry_type<i>::type& get() {
    return _data[runtime_index_of(i)];
  }

  template <unsigned int i>
  const typename entry_type<i>::type& get() const {
    return _data[runtime_index_of(i)];
  }

  value_type& insert(unsigned int i)
  {
    typename std::vector<I>::const_iterator upper_it =
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

  template <unsigned int i>
  typename entry_type<i>::type& insert() {
    return this->insert(i);
  }

  template <unsigned int i, typename T2>
  void set(const T2& val) {
    _data[runtime_index_of(i)] = val;
  }

  bool boolean_test() const {
    std::size_t index_size = size();
    for (unsigned int i=0; i != index_size; ++i)
      if (_data[i])
        return true;
    return false;
  }

  SubType<T,I> operator- () const {
    std::size_t index_size = size();
    SubType<T,I> returnval;
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
  template <typename I2>
  void sparsity_union (const std::vector<I2>& new_indices) {

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

    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<I2>::const_iterator index2_it = new_indices.begin();

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

    typename std::vector<T>::reverse_iterator md_it = _data.rbegin();
    typename std::vector<I>::reverse_iterator mi_it = _indices.rbegin();

    typename std::vector<T>::const_reverse_iterator d_it =
      _data.rbegin() + unseen_indices;
    typename std::vector<I>::const_reverse_iterator i_it =
      _indices.rbegin() + unseen_indices;
    typename std::vector<I2>::const_reverse_iterator i2_it = new_indices.rbegin();

    // Duplicate copies of rend() to work around
    // http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#179
    typename std::vector<I>::reverse_iterator      mirend  = _indices.rend();
    typename std::vector<I>::const_reverse_iterator  rend  = mirend;
    typename std::vector<I2>::const_reverse_iterator rend2 = new_indices.rend();
#ifndef NDEBUG
    typename std::vector<T>::reverse_iterator      mdrend = _data.rend();
    typename std::vector<T>::const_reverse_iterator drend = mdrend;
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
  template <typename I2>
  void sparsity_intersection (const std::vector<I2>& new_indices) {

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

    typedef typename CompareTypes<I,I2>::supertype max_index_type;

#ifndef NDEBUG
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<I2>::const_iterator index2_it = new_indices.begin();

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
    typename std::vector<T>::iterator md_it = _data.begin();
    typename std::vector<I>::iterator mi_it = _indices.begin();

    // Our old values:
    typename std::vector<T>::const_iterator d_it = _data.begin();
    typename std::vector<I>::const_iterator i_it = _indices.begin();

    // Values to merge with:
    typename std::vector<I2>::const_iterator i2_it = new_indices.begin();

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
  void sparsity_trim () {

    metaphysicl_assert
      (std::adjacent_find(_indices.begin(), _indices.end()) ==
       _indices.end());
#ifdef METAPHYSICL_HAVE_CXX11
    metaphysicl_assert(std::is_sorted(_indices.begin(), _indices.end()));
#endif

#ifndef NDEBUG
    I used_indices = 0;

    {
      typename std::vector<I>::iterator index_it = _indices.begin();
      typename std::vector<T>::iterator data_it = _data.begin();
      for (; index_it != _indices.end(); ++index_it, ++data_it)
        if (*data_it)
          ++used_indices;
    }
#endif

    // We'll loop up through the arrays, copying indices (and
    // corresponding data) that should be there downward into place.

    // Downward-merged values:
    typename std::vector<T>::iterator md_it = _data.begin();
    typename std::vector<I>::iterator mi_it = _indices.begin();

    // Our old values:
    typename std::vector<T>::const_iterator d_it = _data.begin();

    for (typename std::vector<I>::const_iterator i_it = _indices.begin();
         i_it != _indices.end(); ++i_it, ++d_it)
      if (*d_it)
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
  // SubType<T,I> operator! () const;

  template <typename T2, typename I2>
  SubType<T,I>&
    operator+= (const SubType<T2,I2>& a) {
    // Resize if necessary
    this->sparsity_union(a.nude_indices());

    typename std::vector<T>::iterator data_it  = _data.begin();
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<T2>::const_iterator data2_it  =
      a.nude_data().begin();
    typename std::vector<I2>::const_iterator index2_it =
      a.nude_indices().begin();
    for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
      {
        I idx1 = *index_it;
        I2 idx2 = *index2_it;

        while (idx1 < idx2) {
          ++index_it;
          ++data_it;
          metaphysicl_assert(index_it != _indices.end());
          idx1 = *index_it;
        }
        metaphysicl_assert_equal_to(idx1, idx2);

        *data_it += *data2_it;
      }

    return static_cast<SubType<T,I>&>(*this);
  }


  template <typename T2, typename I2>
  SubType<T,I>&
    operator-= (const SubType<T2,I2>& a) {
    // Resize if necessary
    this->sparsity_union(a.nude_indices());

    typename std::vector<T>::iterator data_it  = _data.begin();
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<T2>::const_iterator data2_it  =
      a.nude_data().begin();
    typename std::vector<I2>::const_iterator index2_it =
      a.nude_indices().begin();
    for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
      {
        I idx1 = *index_it;
        I2 idx2 = *index2_it;

        while (idx1 < idx2) {
          ++index_it;
          ++data_it;
          metaphysicl_assert(index_it != _indices.end());
          idx1 = *index_it;
        }
        metaphysicl_assert_equal_to(idx1, idx2);

        *data_it -= *data2_it;
      }

    return static_cast<SubType<T,I>&>(*this);
  }


  template <typename T2, typename I2>
  SubType<T,I>&
    operator*= (const SubType<T2,I2>& a) {
    // Resize if possible
    this->sparsity_intersection(a.nude_indices());

    typename std::vector<T>::iterator data_it  = _data.begin();
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<T2>::const_iterator data2_it  =
      a.nude_data().begin();
    typename std::vector<I2>::const_iterator index2_it =
      a.nude_indices().begin();
    for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
      {
        I idx1 = *index_it;
        I2 idx2 = *index2_it;

        while (idx1 < idx2) {
          ++index_it;
          ++data_it;
          metaphysicl_assert(index_it != _indices.end());
          idx1 = *index_it;
        }

        if (idx1 == idx2)
          *data_it *= *data2_it;
      }

    return static_cast<SubType<T,I>&>(*this);
  }


  template <typename T2, typename I2>
  SubType<T,I>&
    operator/= (const SubType<T2,I2>& a) {
    typename std::vector<T>::iterator data_it  = _data.begin();
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<T2>::const_iterator data2_it  =
      a.nude_data().begin();
    typename std::vector<I2>::const_iterator index2_it =
      a.nude_indices().begin();
    for (; data2_it != a.nude_data().end(); ++data2_it, ++index2_it)
      {
        I idx1 = *index_it;
        I2 idx2 = *index2_it;

        while (idx1 < idx2) {
          ++index_it;
          ++data_it;
          metaphysicl_assert(index_it != _indices.end());
          idx1 = *index_it;
        }

        if (idx1 == idx2)
          *data_it /= *data2_it;
      }

    return static_cast<SubType<T,I>&>(*this);
  }


  template <typename T2>
  SubType<T,I>& operator*= (const T2& a) {
    std::size_t index_size = size();
    for (unsigned int i=0; i != index_size; ++i)
      _data[i] *= a;
    return static_cast<SubType<T,I>&>(*this);
  }

  template <typename T2>
  SubType<T,I>& operator/= (const T2& a) {
    std::size_t index_size = size();
    for (unsigned int i=0; i != index_size; ++i)
      _data[i] /= a;
    return static_cast<SubType<T,I>&>(*this);
  }

protected:

  std::vector<T> _data;
  std::vector<I> _indices;
};


//
// Non-member functions
//

template <template <typename, typename> class SubType,
          typename B, typename IB,
          typename T, typename I,
          typename T2, typename I2>
inline
SubType<typename CompareTypes<T,T2>::supertype,
        typename CompareTypes<IB,I2>::supertype>
if_else (const DynamicSparseNumberBase<B, IB,SubType> & condition,
         const DynamicSparseNumberBase<T, I, SubType> & if_true,
         const DynamicSparseNumberBase<T2,I2,SubType> & if_false)
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

  // <I,I2>::supertype would have worked too; every index from the
  // truth case will be in both I and IB.
  typedef typename CompareTypes<IB,I2>::supertype IS;
  typedef typename CompareTypes<T,T2>::supertype TS;

  SubType<TS, IS> returnval;

  // First count returnval size
  IS required_size = 0;
  {
    typename std::vector<IB>::const_iterator indexcond_it      = condition.nude_indices().begin();
    typename std::vector<B>::const_iterator datacond_it        = condition.nude_data().begin();
    typename std::vector<I>::const_iterator indextrue_it       = if_true.nude_indices().begin();
    const typename std::vector<I>::const_iterator endtrue_it   = if_true.nude_indices().end();
    typename std::vector<T>::const_iterator datatrue_it        = if_true.nude_data().begin();
    typename std::vector<I2>::const_iterator indexfalse_it     = if_false.nude_indices().begin();
    const typename std::vector<I2>::const_iterator endfalse_it = if_false.nude_indices().end();
    typename std::vector<T2>::const_iterator datafalse_it      = if_false.nude_data().begin();

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
    typename std::vector<IB>::const_iterator indexcond_it      = condition.nude_indices().begin();
    typename std::vector<B>::const_iterator datacond_it        = condition.nude_data().begin();
    typename std::vector<I>::const_iterator indextrue_it       = if_true.nude_indices().begin();
    const typename std::vector<I>::const_iterator endtrue_it   = if_true.nude_indices().end();
    typename std::vector<T>::const_iterator datatrue_it        = if_true.nude_data().begin();
    typename std::vector<I2>::const_iterator indexfalse_it     = if_false.nude_indices().begin();
    const typename std::vector<I2>::const_iterator endfalse_it = if_false.nude_indices().end();
    typename std::vector<T2>::const_iterator datafalse_it      = if_false.nude_data().begin();

    typename std::vector<IS>::iterator indexreturn_it          = returnval.nude_indices().begin();
    typename std::vector<TS>::iterator datareturn_it           = returnval.nude_data().begin();

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



#define DynamicSparseNumberBase_op_ab(opname, atype, btype, functorname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
typename Symmetric##functorname##Type<atype,btype>::supertype \
operator opname (const atype& a, const btype& b) \
{ \
  typedef typename Symmetric##functorname##Type<atype,btype>::supertype type; \
  type returnval = a; \
  returnval opname##= b; \
  return returnval; \
}


#if __cplusplus >= 201103L

#define DynamicSparseNumberBase_op(subtypename, opname, functorname) \
DynamicSparseNumberBase_op_ab(opname, subtypename<T MacroComma I>, subtypename<T2 MacroComma I2>, functorname) \
 \
template <typename T, typename T2, typename I, typename I2> \
inline \
typename Symmetric##functorname##Type<subtypename<T,I>,subtypename<T2,I2> >::supertype \
operator opname (subtypename<T,I>&& a, \
                 const subtypename<T2,I2>& b) \
{ \
  typedef typename \
    Symmetric##functorname##Type<subtypename<T,I>,subtypename<T2,I2> >::supertype \
    type; \
  type returnval = std::move(a); \
  returnval opname##= b; \
  return returnval; \
}

#else

#define DynamicSparseNumberBase_op(subtypename, opname, functorname) \
DynamicSparseNumberBase_op_ab(opname, subtypename<T MacroComma I>, subtypename<T2 MacroComma I2>, functorname)

#endif

// Let's also allow scalar times vector.
// Scalar plus vector, etc. remain undefined in the sparse context.

template <template <typename, typename> class SubType,
          typename T, typename T2, typename I>
inline
typename MultipliesType<SubType<T2,I>,T,true>::supertype
operator * (const T& a, const DynamicSparseNumberBase<T2,I,SubType>& b)
{
  const unsigned int index_size = b.size();

  typename MultipliesType<SubType<T2,I>,T,true>::supertype
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a * b.raw_at(i);
    returnval.raw_index(i) = b.raw_index(i);
  }

  return returnval;
}

template <template <typename, typename> class SubType,
          typename T, typename T2, typename I>
inline
typename MultipliesType<SubType<T,I>,T2>::supertype
operator * (const DynamicSparseNumberBase<T,I,SubType>& a, const T2& b)
{
  const unsigned int index_size = a.size();

  typename MultipliesType<SubType<T,I>,T2>::supertype
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i) * b;
    returnval.raw_index(i) = a.raw_index(i);
  }
  return returnval;
}

template <template <typename, typename> class SubType,
          typename T, typename T2, typename I>
inline
typename DividesType<SubType<T,I>,T2>::supertype
operator / (const DynamicSparseNumberBase<T,I,SubType>& a, const T2& b)
{
  const unsigned int index_size = a.size();

  typename DividesType<SubType<T,I>,T2>::supertype returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i) / b;
    returnval.raw_index(i) = a.raw_index(i);
  }

  return returnval;
}

#if __cplusplus >= 201103L
template <template <typename, typename> class SubType,
          typename T, typename T2, typename I>
inline
typename MultipliesType<SubType<T,I>,T2>::supertype
operator * (DynamicSparseNumberBase<T,I,SubType>&& a, const T2& b)
{
  typename MultipliesType<SubType<T,I>,T2>::supertype
    returnval = std::move(static_cast<SubType<T,I>&&>(a));

  returnval *= b;

  return returnval;
}

template <template <typename, typename> class SubType,
          typename T, typename T2, typename I>
inline
typename DividesType<SubType<T,I>,T2>::supertype
operator / (DynamicSparseNumberBase<T,I,SubType>&& a, const T2& b)
{
  typename DividesType<SubType<T,I>,T2>::supertype
    returnval = std::move(static_cast<SubType<T,I>&&>(a));

  returnval /= b;

  return returnval;
}
#endif


#define DynamicSparseNumberBase_operator_binary(opname, functorname) \
template <template <typename, typename> class SubType, \
          typename T, typename T2, typename I, typename I2> \
inline \
SubType<bool, typename CompareTypes<I,I2>::supertype> \
operator opname (const DynamicSparseNumberBase<T,I,SubType>& a, \
                 const DynamicSparseNumberBase<T2,I2,SubType>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<I,I2>::supertype IS; \
  SubType<bool, IS> returnval; \
  returnval.nude_indices() = a.nude_indices(); \
  returnval.nude_data().resize(a.nude_indices().size()); \
  returnval.sparsity_union(b.nude_indices()); \
 \
  typename std::vector<I>::const_iterator  index_a_it = a.nude_indices().begin(); \
  typename std::vector<I2>::const_iterator index_b_it = b.nude_indices().begin(); \
  typename std::vector<IS>::iterator     index_out_it = returnval.nude_indices().begin(); \
 \
  typename std::vector<T>::const_iterator  data_a_it = a.nude_data().begin(); \
  typename std::vector<T2>::const_iterator data_b_it = b.nude_data().begin(); \
  typename std::vector<TS>::iterator     data_out_it = returnval.nude_data().begin(); \
 \
  const IS  maxIS  = std::numeric_limits<IS>::max(); \
 \
  for (; index_out_it != returnval.nude_indices().end(); ++index_out_it, ++data_out_it) { \
    const IS index_a = (index_a_it == a.nude_indices().end()) ? maxIS : *index_a_it; \
    const IS index_b = (index_b_it == b.nude_indices().end()) ? maxIS : *index_b_it; \
    const IS index_out = *index_out_it; \
    const TS data_a  = (index_a_it == a.nude_indices().end()) ? 0: *data_a_it; \
    const TS data_b  = (index_b_it == b.nude_indices().end()) ? 0: *data_b_it; \
    typename std::vector<TS>::reference data_out = *data_out_it; \
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
template <template <typename, typename> class SubType, \
          typename T, typename T2, typename I> \
inline \
SubType<bool, I> \
operator opname (const DynamicSparseNumberBase<T,I,SubType>& a, const T2& b) \
{ \
  SubType<bool, I> returnval; \
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
template <template <typename, typename> class SubType, \
          typename T, typename T2, typename I> \
inline \
SubType<bool, I> \
operator opname (const T& a, const DynamicSparseNumberBase<T2,I,SubType>& b) \
{ \
  SubType<bool, I> returnval; \
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


template <template <typename, typename> class SubType,
          typename T, typename I>
inline
std::ostream&
operator<< (std::ostream& output, const DynamicSparseNumberBase<T,I,SubType>& a)
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


// CompareTypes, RawType, ValueType specializations

#define DynamicSparseNumberBase_comparisons(subtypename, templatename) \
template<typename T, typename I, bool reverseorder> \
struct templatename<subtypename<T,I>, subtypename<T,I>, reverseorder> { \
  typedef subtypename<T,I> supertype; \
}; \
 \
template<typename T, typename T2, typename I, typename I2, bool reverseorder> \
struct templatename<subtypename<T,I>, subtypename<T2,I2>, reverseorder> { \
  typedef subtypename<typename Symmetric##templatename<T, T2, reverseorder>::supertype, \
                      typename CompareTypes<I,I2>::supertype> supertype; \
}; \
 \
template<typename T, typename T2, typename I, bool reverseorder> \
struct templatename<subtypename<T, I>, T2, reverseorder, \
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> { \
  typedef subtypename<typename Symmetric##templatename<T, T2, reverseorder>::supertype, I> supertype; \
}

} // namespace MetaPhysicL

#if 0

namespace std {

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DynamicSparseNumberArray;
using MetaPhysicL::SymmetricCompareTypes;

#define DynamicSparseNumberArray_std_unary(funcname) \
template <typename T, typename I> \
inline \
DynamicSparseNumberArray<T, I> \
funcname (DynamicSparseNumberArray<T, I> a) \
{ \
  std::size_t index_size = a.size(); \
  for (unsigned int i=0; i != index_size; ++i) \
    a.raw_at(i) = std::funcname(a.raw_at(i)); \
 \
  return a; \
}


#define DynamicSparseNumberArray_std_binary_union(funcname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, \
                         typename CompareTypes<I,I2>::supertype> \
funcname (const DynamicSparseNumberArray<T, I>& a, \
          const DynamicSparseNumberArray<T2, I2>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<I,I2>::supertype IS; \
  DynamicSparseNumberArray<TS, IS> returnval; \
 \
  std::size_t index_size = a.nude_indices.size(); \
  returnval.nude_indices = a.nude_indices; \
  returnval.nude_data.resize(index_size); \
  returnval.sparsity_union(b.nude_indices); \
 \
  typename std::vector<I>::const_iterator  index_a_it = a.nude_indices.begin(); \
  typename std::vector<I2>::const_iterator index_b_it = b.nude_indices.begin(); \
  typename std::vector<IS>::iterator     index_out_it = returnval.nude_indices.begin(); \
 \
  typename std::vector<T>::const_iterator  data_a_it = a.nude_data.begin(); \
  typename std::vector<T2>::const_iterator data_b_it = b.nude_data.begin(); \
  typename std::vector<TS>::iterator     data_out_it = returnval.nude_data.begin(); \
 \
  const IS  maxIS  = std::numeric_limits<IS>::max(); \
 \
  for (; index_out_it != returnval.nude_indices.end(); ++index_out_it, ++data_out_it) { \
    const IS index_a = (index_a_it == a.nude_indices.end()) ? maxIS : *index_a_it; \
    const IS index_b = (index_b_it == b.nude_indices.end()) ? maxIS : *index_b_it; \
    const IS index_out = *index_out_it; \
    const TS data_a  = (index_a_it == a.nude_indices.end()) ? 0: *data_a_it; \
    const TS data_b  = (index_b_it == b.nude_indices.end()) ? 0: *data_b_it; \
    typename std::vector<TS>::reference data_out = *data_out_it; \
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
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const DynamicSparseNumberArray<T, I>& a, const T2& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  DynamicSparseNumberArray<TS, I> returnval; \
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
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const T& a, const DynamicSparseNumberArray<T2, I>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  DynamicSparseNumberArray<TS, I> returnval; \
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


// Pow needs its own specialization, both to avoid being confused by
// pow<T1,T2> and because pow(x,0) isn't 0.
template <typename T, typename T2, typename I>
inline
DynamicSparseNumberArray<typename SymmetricCompareTypes<T,T2>::supertype, I>
pow (const DynamicSparseNumberArray<T, I>& a, const T2& b)
{
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS;
  DynamicSparseNumberArray<TS, I> returnval;

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

// DynamicSparseNumberArray_std_binary(pow) // separate definition
// DynamicSparseNumberArray_std_unary(exp)
// DynamicSparseNumberArray_std_unary(log)
// DynamicSparseNumberArray_std_unary(log10)
DynamicSparseNumberArray_std_unary(sin)
// DynamicSparseNumberArray_std_unary(cos)
DynamicSparseNumberArray_std_unary(tan)
DynamicSparseNumberArray_std_unary(asin)
// DynamicSparseNumberArray_std_unary(acos)
DynamicSparseNumberArray_std_unary(atan)
DynamicSparseNumberArray_std_binary_union(atan2)
DynamicSparseNumberArray_std_unary(sinh)
// DynamicSparseNumberArray_std_unary(cosh)
DynamicSparseNumberArray_std_unary(tanh)
DynamicSparseNumberArray_std_unary(sqrt)
DynamicSparseNumberArray_std_unary(abs)
DynamicSparseNumberArray_std_unary(fabs)
DynamicSparseNumberArray_std_binary_union(max)
DynamicSparseNumberArray_std_binary_union(min)
DynamicSparseNumberArray_std_unary(ceil)
DynamicSparseNumberArray_std_unary(floor)
DynamicSparseNumberArray_std_binary_union(fmod) // TODO: optimize this


template <typename T, typename I>
class numeric_limits<DynamicSparseNumberArray<T, I> > :
  public MetaPhysicL::raw_numeric_limits<DynamicSparseNumberArray<T, I>, T> {};

} // namespace std

#endif // 0


#endif // METAPHYSICL_DYNAMICSPARSENUMBERARRAY_H
