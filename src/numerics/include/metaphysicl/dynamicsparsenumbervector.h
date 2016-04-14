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


#ifndef METAPHYSICL_DYNAMICSPARSENUMBERVECTOR_H
#define METAPHYSICL_DYNAMICSPARSENUMBERVECTOR_H

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <ostream>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/ct_set.h"
#include "metaphysicl/metaphysicl_asserts.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/sparsenumberutils.h"
#include "metaphysicl/testable.h"

namespace MetaPhysicL {

// Forward declarations

// Data type T, index type I
template <typename T, typename I>
class DynamicSparseNumberVector;

// Helper structs

template<typename I1, typename I2, typename S, typename T, bool reverseorder>
struct DotType<DynamicSparseNumberVector<S,I1>,
               DynamicSparseNumberVector<T,I2>, reverseorder> {
  typedef typename MultipliesType<S,T,reverseorder>::supertype supertype;
};

template<typename I1, typename I2, typename S, typename T, bool reverseorder>
struct OuterProductType<DynamicSparseNumberVector<S, I1>,
                        DynamicSparseNumberVector<T, I2>, reverseorder> {
  typedef
    DynamicSparseNumberVector<DynamicSparseNumberVector<
      typename MultipliesType<S,T,reverseorder>::supertype,
      I2>, I1> supertype;
};

template<typename S, typename I>
struct SumType<DynamicSparseNumberVector<S, I> > {
  typedef S supertype;
};


template <typename T, typename I>
class DynamicSparseNumberVector :
  public safe_bool<DynamicSparseNumberVector<T,I> >
{
public:
  typedef T value_type;

  template <unsigned int i>
  struct entry_type {
    typedef value_type type;
  };

  typedef I index_value_type;

  template <typename T2>
  struct rebind {
    typedef DynamicSparseNumberVector<T2, I> other;
  };

  std::size_t size() const
    { metaphysicl_assert_equal_to(_data.size(), _indices.size());
      return _data.size(); }

  void resize(std::size_t s)
    { metaphysicl_assert_equal_to(_data.size(), _indices.size());
      _data.resize(s);
      _indices.resize(s); }

  DynamicSparseNumberVector() {}

  DynamicSparseNumberVector(const T& val) {
    // This makes no sense unless val is 0!
#ifndef NDEBUG
    if (val)
      throw std::domain_error("Cannot initialize DynamicSparseNumberVector with non-zero scalar");
#endif
  }

  template <typename T2>
  DynamicSparseNumberVector(const T2& val) {
    // This makes no sense unless val is 0!
#ifndef NDEBUG
    if (val)
      throw std::domain_error("Cannot initialize DynamicSparseNumberVector with non-zero scalar");
#endif
  }

#if __cplusplus >= 201103L
  // Move constructors are useful when all your data is on the heap
  DynamicSparseNumberVector(DynamicSparseNumberVector<T, I> && src) = default;

  // Move assignment avoids heap operations too
  DynamicSparseNumberVector& operator= (DynamicSparseNumberVector<T, I> && src) = default;

  // Standard copy operations get implicitly deleted upon move
  // constructor definition, so we redefine them.
  DynamicSparseNumberVector(const DynamicSparseNumberVector<T, I> & src) = default;

  DynamicSparseNumberVector& operator= (const DynamicSparseNumberVector<T, I> & src) = default;
#endif

  template <typename T2, typename I2>
  DynamicSparseNumberVector(DynamicSparseNumberVector<T2, I2> src)
    { this->resize(src.size());
      std::copy(src.nude_data().begin(), src.nude_data().end(), _data.begin());
      std::copy(src.nude_indices().begin(), src.nude_indices().end(), _indices.begin()); }

  T* raw_data()
    { return size()?&_data[0]:NULL; }

  const T* raw_data() const
    { return size()?&_data[0]:NULL; }

  T& raw_at(unsigned int i)
    { return _data[i]; }

  const T& raw_at(unsigned int i) const
    { return _data[i]; }

  I& raw_index(unsigned int i)
    { return _indices[i]; }

  const I& raw_index(unsigned int i) const
    { return _indices[i]; }

  // FIXME: these encapsulation violations are necessary
  // until I can figure out the right friend declaration.
  const std::vector<T>& nude_data() const
    { return _data; }

  std::vector<T>& nude_data()
    { return _data; }

  const std::vector<I>& nude_indices() const
    { return _indices; }

  std::vector<I>& nude_indices()
    { return _indices; }

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

  DynamicSparseNumberVector<T,I> operator- () const {
    std::size_t index_size = size();
    DynamicSparseNumberVector<T,I> returnval;
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

#ifndef NDEBUG
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<I2>::const_iterator index2_it = new_indices.begin();

    typedef typename CompareTypes<I,I2>::supertype max_index_type;
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

    // We'll loop up through the array, copying indices (and
    // corresponding data) that should be there downward into place.

    // Merged values:
    typename std::vector<T>::iterator md_it = _data.begin();
    typename std::vector<I>::iterator mi_it = _indices.begin();

    // Our old values:
    typename std::vector<T>::iterator d_it = _data.begin();
    typename std::vector<I>::iterator i_it = _indices.begin();

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


  // Not defineable since !0 != 0
  // DynamicSparseNumberVector<T,I> operator! () const;


  template <typename T2, typename I2>
  DynamicSparseNumberVector<T,I>&
    operator+= (const DynamicSparseNumberVector<T2,I2>& a) {
    // Resize if necessary
    this->sparsity_union(a.nude_indices());

    typename std::vector<T>::iterator data_it  = _data.begin();
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<T2>::const_iterator data2_it  = a.nude_data().begin();
    typename std::vector<I2>::const_iterator index2_it = a.nude_indices().begin();
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

    return *this;
  }


  template <typename T2, typename I2>
  DynamicSparseNumberVector<T,I>&
    operator-= (const DynamicSparseNumberVector<T2,I2>& a) {
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

    return *this;
  }


  template <typename T2, typename I2>
  DynamicSparseNumberVector<T,I>&
    operator*= (const DynamicSparseNumberVector<T2,I2>& a) {
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

    return *this;
  }


  template <typename T2, typename I2>
  DynamicSparseNumberVector<T,I>&
    operator/= (const DynamicSparseNumberVector<T2,I2>& a) {
    typename std::vector<T>::iterator data_it  = _data.begin();
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<T2>::const_iterator data2_it  = a.nude_data().begin();
    typename std::vector<I2>::const_iterator index2_it = a.nude_indices().begin();
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

    return *this;
  }


  template <typename T2>
  DynamicSparseNumberVector<T,I>& operator*= (const T2& a) {
    std::size_t index_size = size();
    for (unsigned int i=0; i != index_size; ++i)
      _data[i] *= a;
    return *this;
  }

  template <typename T2>
  DynamicSparseNumberVector<T,I>& operator/= (const T2& a) {
    std::size_t index_size = size();
    for (unsigned int i=0; i != index_size; ++i)
      _data[i] /= a;
    return *this;
  }

  template <typename T2, typename I2>
  typename MultipliesType<T,T2>::supertype
  dot (const DynamicSparseNumberVector<T2,I2>& a) const
  {
    typename MultipliesType<T,T2>::supertype returnval = 0;

    for (I i1 = 0; i1 != _indices.size(); ++i1)
      {
        typename std::vector<I2>::const_iterator it2 =
          std::lower_bound(a.nude_indices().begin(),
                           a.nude_indices().end(),
                           _indices[i1]);

        if (it2 != a.nude_indices().end())
          {
            std::size_t i2 = it2 - a.nude_indices().begin();

            returnval += _data[i1] * a.raw_at(i2);
          }
      }

    return returnval;
  }

  template <typename T2, typename I2>
  DynamicSparseNumberVector<DynamicSparseNumberVector<
    typename MultipliesType<T,T2>::supertype,
    I2>, I>
  outerproduct (const DynamicSparseNumberVector<T2, I2>& a) const
  {
    DynamicSparseNumberVector<DynamicSparseNumberVector<
      typename MultipliesType<T,T2>::supertype,
      I2>, I> returnval;

    returnval.nude_indices() = this->_indices;

    std::size_t index_size = size();
    std::size_t index2_size = a.size();

    returnval.nude_data().resize(index_size);
    for (unsigned int i=0; i != index_size; ++i)
      {
        returnval.raw_at(i).nude_indices() = a.nude_indices();

        returnval.raw_at(i).nude_data().resize(index2_size);
        for (unsigned int j=0; j != index2_size; ++j)
          returnval.raw_at(i).raw_at(j) = _data[i] * a.raw_at(j);
      }

    return returnval;
  }

  static DynamicSparseNumberVector<DynamicSparseNumberVector<T, I>, I>
  identity(std::size_t n = 0)
  {
    DynamicSparseNumberVector<DynamicSparseNumberVector<T, I>, I>
      returnval;
    returnval.resize(n);
    for (unsigned int i=0; i != n; ++i)
      {
        returnval.raw_index(i) = i;
        returnval.raw_at(i).nude_indices().resize(1, i);
        returnval.raw_at(i).nude_data().resize(1, 1);
      }
    return returnval;
  }

private:

  std::vector<T> _data;
  std::vector<I> _indices;
};


//
// Non-member functions
//

template <unsigned int N,
          unsigned int index1=0, typename Data1=void,
          unsigned int index2=0, typename Data2=void,
          unsigned int index3=0, typename Data3=void,
          unsigned int index4=0, typename Data4=void,
          unsigned int index5=0, typename Data5=void,
          unsigned int index6=0, typename Data6=void,
          unsigned int index7=0, typename Data7=void,
          unsigned int index8=0, typename Data8=void>
struct DynamicSparseNumberVectorOf
{
  typedef
  typename SymmetricCompareTypes<Data1,
    typename SymmetricCompareTypes<Data2,
      typename SymmetricCompareTypes<Data3,
        typename SymmetricCompareTypes<Data4,
          typename SymmetricCompareTypes<Data5,
            typename SymmetricCompareTypes<Data6,
              typename SymmetricCompareTypes<Data7,Data8>::supertype
            >::supertype
          >::supertype
        >::supertype
      >::supertype
    >::supertype
  >::supertype supertype;

  typedef DynamicSparseNumberVector<supertype, unsigned int> type;
};



template <std::size_t N, unsigned int index, typename T>
struct DynamicSparseNumberVectorUnitVector
{
  typedef DynamicSparseNumberVector<T, unsigned int> type;

  static type value() {
    type returnval;
    returnval.resize(1);
    returnval.raw_at(0) = 1;
    returnval.raw_index(0) = index;
    return returnval;
  }
};


template <std::size_t N, typename T>
struct DynamicSparseNumberVectorFullVector
{
  typedef DynamicSparseNumberVector<T,unsigned int> type;

  static type value() {
    type returnval;
    returnval.resize(N);
    for (unsigned int i=0; i != N; ++i)
      {
        returnval.raw_at(i) = 1;
        returnval.raw_index(i) = i;
      }
    return returnval;
  }
};



template <typename T, typename I, typename I2>
inline
DynamicSparseNumberVector<DynamicSparseNumberVector<T, I>, I2>
transpose(const DynamicSparseNumberVector<DynamicSparseNumberVector<T, I2>, I>& a)
{
  DynamicSparseNumberVector<DynamicSparseNumberVector<T, I>, I2> returnval;

  metaphysicl_not_implemented();

  return returnval;
}


template <typename T, typename I>
T
sum (const DynamicSparseNumberVector<T, I> &a)
{
  std::size_t index_size = a.size();

  T returnval = 0;

  for (unsigned int i=0; i != index_size; ++i) {
    returnval += a.raw_at(i);
  }

  return returnval;
}



#define DynamicSparseNumberVector_op_ab(opname, atype, btype, functorname) \
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

#define DynamicSparseNumberVector_op(opname, functorname) \
DynamicSparseNumberVector_op_ab(opname, DynamicSparseNumberVector<T MacroComma I>, DynamicSparseNumberVector<T2 MacroComma I2>, functorname) \
 \
template <typename T, typename T2, typename I, typename I2> \
inline \
typename Symmetric##functorname##Type<DynamicSparseNumberVector<T,I>,DynamicSparseNumberVector<T2,I2> >::supertype \
operator opname (DynamicSparseNumberVector<T,I>&& a, \
                 const DynamicSparseNumberVector<T2,I2>& b) \
{ \
  typedef typename \
    Symmetric##functorname##Type<DynamicSparseNumberVector<T,I>,DynamicSparseNumberVector<T2,I2> >::supertype \
    type; \
  type returnval = std::move(a); \
  returnval opname##= b; \
  return returnval; \
}

#else

#define DynamicSparseNumberVector_op(opname, functorname) \
DynamicSparseNumberVector_op_ab(opname, DynamicSparseNumberVector<T MacroComma I>, DynamicSparseNumberVector<T2 MacroComma I2>, functorname)

#endif

DynamicSparseNumberVector_op(+, Plus)       // Union)
DynamicSparseNumberVector_op(-, Minus)      // Union)
DynamicSparseNumberVector_op(*, Multiplies) // Intersection)
DynamicSparseNumberVector_op(/, Divides)    // First)

// Let's also allow scalar times vector.
// Scalar plus vector, etc. remain undefined in the sparse context.

template <typename T, typename T2, typename I>
inline
typename MultipliesType<DynamicSparseNumberVector<T2,I>,T,true>::supertype
operator * (const T& a, const DynamicSparseNumberVector<T2,I>& b)
{
  const unsigned int index_size = b.size();

  typename MultipliesType<DynamicSparseNumberVector<T2,I>,T,true>::supertype
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a * b.raw_at(i);
    returnval.raw_index(i) = b.raw_index(i);
  }

  return returnval;
}

template <typename T, typename T2, typename I>
inline
typename MultipliesType<DynamicSparseNumberVector<T,I>,T2>::supertype
operator * (const DynamicSparseNumberVector<T,I>& a, const T2& b)
{
  const unsigned int index_size = a.size();

  typename MultipliesType<DynamicSparseNumberVector<T,I>,T2>::supertype
    returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i) * b;
    returnval.raw_index(i) = a.raw_index(i);
  }
  return returnval;
}

template <typename T, typename T2, typename I>
inline
typename DividesType<DynamicSparseNumberVector<T,I>,T2>::supertype
operator / (const DynamicSparseNumberVector<T,I>& a, const T2& b)
{
  const unsigned int index_size = a.size();

  typename DividesType<DynamicSparseNumberVector<T,I>,T2>::supertype returnval;
  returnval.resize(index_size);

  for (unsigned int i=0; i != index_size; ++i) {
    returnval.raw_at(i) = a.raw_at(i) / b;
    returnval.raw_index(i) = a.raw_index(i);
  }

  return returnval;
}

#if __cplusplus >= 201103L
template <typename T, typename T2, typename I>
inline
typename MultipliesType<DynamicSparseNumberVector<T,I>,T2>::supertype
operator * (DynamicSparseNumberVector<T,I>&& a, const T2& b)
{
  typename MultipliesType<DynamicSparseNumberVector<T,I>,T2>::supertype
    returnval = std::move(a);

  returnval *= b;

  return returnval;
}

template <typename T, typename T2, typename I>
inline
typename DividesType<DynamicSparseNumberVector<T,I>,T2>::supertype
operator / (DynamicSparseNumberVector<T,I>&& a, const T2& b)
{
  const unsigned int index_size = a.size();

  typename DividesType<DynamicSparseNumberVector<T,I>,T2>::supertype returnval;
  returnval = std::move(a);

  returnval /= b;

  return returnval;
}
#endif


#define DynamicSparseNumberVector_operator_binary(opname, functorname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
DynamicSparseNumberVector<bool, typename CompareTypes<I,I2>::supertype> \
operator opname (const DynamicSparseNumberVector<T,I>& a, \
                 const DynamicSparseNumberVector<T2,I2>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<I,I2>::supertype IS; \
  DynamicSparseNumberVector<bool, IS> returnval; \
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
    TS &   data_out  = *data_out_it; \
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
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberVector<bool, I> \
operator opname (const DynamicSparseNumberVector<T, I>& a, const T2& b) \
{ \
  DynamicSparseNumberVector<bool, I> returnval; \
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
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberVector<bool, I> \
operator opname (const T& a, const DynamicSparseNumberVector<T2,I>& b) \
{ \
  DynamicSparseNumberVector<bool, I> returnval; \
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

DynamicSparseNumberVector_operator_binary(<, less)
// DynamicSparseNumberVector_operator_binary(<=)
DynamicSparseNumberVector_operator_binary(>, greater)
// DynamicSparseNumberVector_operator_binary(>=)
// DynamicSparseNumberVector_operator_binary(==)
DynamicSparseNumberVector_operator_binary(!=, not_equal_to)

// FIXME - make && an intersection rather than a union for efficiency
DynamicSparseNumberVector_operator_binary(&&, logical_and)
DynamicSparseNumberVector_operator_binary(||, logical_or)

template <typename T, typename I>
inline
std::ostream&
operator<< (std::ostream& output, const DynamicSparseNumberVector<T, I>& a)
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

#define DynamicSparseNumberVector_comparisons(templatename, settype) \
template<typename T, typename I, bool reverseorder> \
struct templatename<DynamicSparseNumberVector<T,I>, DynamicSparseNumberVector<T,I>, reverseorder> { \
  typedef DynamicSparseNumberVector<T,I> supertype; \
}; \
 \
template<typename T, typename T2, typename I, typename I2, bool reverseorder> \
struct templatename<DynamicSparseNumberVector<T,I>, DynamicSparseNumberVector<T2,I2>, reverseorder> { \
  typedef DynamicSparseNumberVector<typename Symmetric##templatename<T, T2, reverseorder>::supertype, \
                            typename CompareTypes<I,I2>::supertype> supertype; \
}; \
 \
template<typename T, typename T2, typename I, bool reverseorder> \
struct templatename<DynamicSparseNumberVector<T, I>, T2, reverseorder, \
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> { \
  typedef DynamicSparseNumberVector<typename Symmetric##templatename<T, T2, reverseorder>::supertype, I> supertype; \
}

DynamicSparseNumberVector_comparisons(CompareTypes, Union);
DynamicSparseNumberVector_comparisons(PlusType, Union);
DynamicSparseNumberVector_comparisons(MinusType, Union);
DynamicSparseNumberVector_comparisons(MultipliesType, Intersection);
DynamicSparseNumberVector_comparisons(DividesType, First);
DynamicSparseNumberVector_comparisons(AndType, Intersection);
DynamicSparseNumberVector_comparisons(OrType, Union);


template <typename T, typename I>
struct RawType<DynamicSparseNumberVector<T, I> >
{
  typedef DynamicSparseNumberVector<typename RawType<T>::value_type, I> value_type;

  static value_type value(const DynamicSparseNumberVector<T, I>& a)
    {
      value_type returnval;
      returnval.nude_indices() = a.nude_indices();

      std::size_t index_size = a.size();
      returnval.nude_data().resize(index_size);

      for (unsigned int i=0; i != index_size; ++i)
        returnval.raw_at(i) = RawType<T>::value(a.raw_at(i));
      return returnval;
    }
};

template <typename T, typename I>
struct ValueType<DynamicSparseNumberVector<T, I> >
{
  typedef typename ValueType<T>::type type;
};

} // namespace MetaPhysicL


namespace std {

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DynamicSparseNumberVector;
using MetaPhysicL::SymmetricCompareTypes;

#define DynamicSparseNumberVector_std_unary(funcname) \
template <typename T, typename I> \
inline \
DynamicSparseNumberVector<T, I> \
funcname (DynamicSparseNumberVector<T, I> a) \
{ \
  std::size_t index_size = a.size(); \
  for (unsigned int i=0; i != index_size; ++i) \
    a.raw_at(i) = std::funcname(a.raw_at(i)); \
 \
  return a; \
}


#define DynamicSparseNumberVector_std_binary_union(funcname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
DynamicSparseNumberVector<typename SymmetricCompareTypes<T,T2>::supertype, \
                         typename CompareTypes<I,I2>::supertype> \
funcname (const DynamicSparseNumberVector<T, I>& a, \
          const DynamicSparseNumberVector<T2, I2>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<I,I2>::supertype IS; \
  DynamicSparseNumberVector<TS, IS> returnval; \
 \
  std::size_t index_size = a.nude_indices().size(); \
  returnval.nude_indices() = a.nude_indices(); \
  returnval.nude_data().resize(index_size); \
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
    TS &   data_out  = *data_out_it; \
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
DynamicSparseNumberVector<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const DynamicSparseNumberVector<T, I>& a, const T2& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  DynamicSparseNumberVector<TS, I> returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices() = a.nude_indices(); \
 \
  for (unsigned int i=0; i != index_size; ++i) \
    returnval.raw_at(i) = std::funcname(a.raw_at(i), b); \
 \
  return returnval; \
} \
 \
template <typename T, typename T2, typename I> \
inline \
DynamicSparseNumberVector<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const T& a, const DynamicSparseNumberVector<T2, I>& b) \
{ \
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS; \
  DynamicSparseNumberVector<TS, I> returnval; \
 \
  std::size_t index_size = a.size(); \
  returnval.resize(index_size); \
  returnval.nude_indices() = b.nude_indices(); \
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
DynamicSparseNumberVector<typename SymmetricCompareTypes<T,T2>::supertype, I>
pow (const DynamicSparseNumberVector<T, I>& a, const T2& b)
{
  typedef typename SymmetricCompareTypes<T,T2>::supertype TS;
  DynamicSparseNumberVector<TS, I> returnval;

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

// DynamicSparseNumberVector_std_binary(pow) // separate definition
// DynamicSparseNumberVector_std_unary(exp)
// DynamicSparseNumberVector_std_unary(log)
// DynamicSparseNumberVector_std_unary(log10)
DynamicSparseNumberVector_std_unary(sin)
// DynamicSparseNumberVector_std_unary(cos)
DynamicSparseNumberVector_std_unary(tan)
DynamicSparseNumberVector_std_unary(asin)
// DynamicSparseNumberVector_std_unary(acos)
DynamicSparseNumberVector_std_unary(atan)
DynamicSparseNumberVector_std_binary_union(atan2)
DynamicSparseNumberVector_std_unary(sinh)
// DynamicSparseNumberVector_std_unary(cosh)
DynamicSparseNumberVector_std_unary(tanh)
DynamicSparseNumberVector_std_unary(sqrt)
DynamicSparseNumberVector_std_unary(abs)
DynamicSparseNumberVector_std_unary(fabs)
DynamicSparseNumberVector_std_binary_union(max)
DynamicSparseNumberVector_std_binary_union(min)
DynamicSparseNumberVector_std_unary(ceil)
DynamicSparseNumberVector_std_unary(floor)
DynamicSparseNumberVector_std_binary_union(fmod) // TODO: optimize this


template <typename T, typename I>
class numeric_limits<DynamicSparseNumberVector<T, I> > :
  public MetaPhysicL::raw_numeric_limits<DynamicSparseNumberVector<T, I>, T> {};

} // namespace std


#endif // METAPHYSICL_DYNAMICSPARSENUMBERVECTOR_H
