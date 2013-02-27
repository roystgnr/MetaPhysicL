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


#ifndef METAPHYSICL_NUMBERARRAY_H
#define METAPHYSICL_NUMBERARRAY_H

#include <algorithm>
#include <ostream>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/raw_type.h"

namespace MetaPhysicL {

// Forward declarations
template<std::size_t size, typename T>
class NumberArray;

template<std::size_t size, typename S, typename T, bool reverseorder>
struct DotType<NumberArray<size,S>, NumberArray<size,T>, reverseorder> {
  typedef NumberArray<size, typename DotType<S,T,reverseorder>::supertype> supertype;
};

template<std::size_t size, typename S, typename T, bool reverseorder>
struct OuterProductType<NumberArray<size,S>, NumberArray<size,T>, reverseorder> {
  typedef 
  NumberArray<size, typename OuterProductType<S,T,reverseorder>::supertype> supertype;
};

template<std::size_t size, typename S>
struct SumType<NumberArray<size,S> > {
  typedef NumberArray<size, typename SumType<S>::supertype> supertype;
};

template <std::size_t size, typename T>
class NumberArray
{
public:
  typedef T value_type;

  template <std::size_t i>
  struct entry_type {
    typedef value_type type;
  };

  template <typename T2>
  struct rebind {
    typedef NumberArray<size, T2> other;
  };

  NumberArray() {}

  NumberArray(const T& val)
    { std::fill(_data, _data+size, val); }

  NumberArray(const T* vals)
    { std::copy(vals, vals+size, _data); }

  template <typename T2>
  NumberArray(NumberArray<size, T2> src)
    { if (size) std::copy(&src[0], &src[0]+size, _data); }

  template <typename T2>
  NumberArray(const T2& val)
    { std::fill(_data, _data+size, T(val)); }

  T& operator[](std::size_t i)
    { return _data[i]; }

  const T& operator[](std::size_t i) const
    { return _data[i]; }

  template <std::size_t i>
  typename entry_type<i>::type& get()
    { return _data[i]; }

  template <std::size_t i>
  const typename entry_type<i>::type& get() const
    { return _data[i]; }

  NumberArray<size,T> operator- () const {
    NumberArray<size,T> returnval;
    for (std::size_t i=0; i != size; ++i) returnval[i] = -_data[i];
    return returnval;
  }

  template <typename T2>
  NumberArray<size,T>& operator+= (const NumberArray<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] += a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator+= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] += a; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator-= (const NumberArray<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] -= a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator-= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] -= a; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator*= (const NumberArray<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] *= a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator*= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] *= a; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator/= (const NumberArray<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] /= a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator/= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] /= a; return *this; }

  template <typename T2>
  NumberArray<size, typename DotType<T,T2>::supertype>
  dot (const NumberArray<size,T2>& a) const
  {
    NumberArray<size, typename DotType<T,T2>::supertype> returnval;
    for (std::size_t i=0; i != size; ++i)
      returnval[i] = _data[i].dot(a[i]);
    return returnval;
  }

  template <typename T2>
  typename OuterProductType<NumberArray<size,T>,NumberArray<size,T2> >::supertype
  outerproduct (const NumberArray<size,T2>& a) const
  {
    typename OuterProductType<NumberArray<size,T>,NumberArray<size,T2> >::supertype
      returnval;

    for (std::size_t i=0; i != size; ++i)
      returnval[i] = _data[i].outerproduct(a[i]);

    return returnval;
  }

  T
  sum () const
  {
    T returnval = 0;
    
    for (std::size_t i=0; i != size; ++i)
      returnval[i] = _data[i].sum();

    return returnval;
  }


private:
  T _data[size];
};



//
// Non-member functions
//

template <std::size_t size,
          unsigned int index1=0, typename Data1=void,
          unsigned int index2=0, typename Data2=void,
          unsigned int index3=0, typename Data3=void,
          unsigned int index4=0, typename Data4=void,
          unsigned int index5=0, typename Data5=void,
          unsigned int index6=0, typename Data6=void,
          unsigned int index7=0, typename Data7=void,
          unsigned int index8=0, typename Data8=void>
struct NumberArrayOf
{
  typedef
  typename CompareTypes<Data1,
    typename CompareTypes<Data2,
      typename CompareTypes<Data3,
        typename CompareTypes<Data4,
          typename CompareTypes<Data5,
            typename CompareTypes<Data6,
              typename CompareTypes<Data7,Data8>::supertype
            >::supertype
          >::supertype
        >::supertype
      >::supertype
    >::supertype
  >::supertype supertype;

  typedef NumberArray<size, supertype> type;
};



template <std::size_t size, std::size_t index, typename T>
struct NumberArrayUnitVector
{
  typedef NumberArray<size, T> type;

  static const type value() {
    type returnval = 0;
    returnval[index] = 1;
    return returnval;
  }
};


template <std::size_t size, typename T>
struct NumberArrayFullVector
{
  typedef NumberArray<size,T> type;

  static const type value() {
    type returnval;
    for (std::size_t i=0; i != size; ++i)
      returnval[i] = 1;
    return returnval;
  }
};



template <std::size_t size, typename T>
inline
NumberArray<size, NumberArray<size, T> >
transpose(const NumberArray<size, NumberArray<size, T> > a)
{
  for (std::size_t i=0; i != size; ++i)
    for (std::size_t j=i+1; j != size; ++j)
      std::swap(a[i][j], a[j][i]);

  return a;
}



#define NumberArray_op_ab(opname, atype, btype, newtype) \
template <std::size_t size, typename T, typename T2> \
inline \
typename newtype::supertype \
operator opname (const atype& a, const btype& b) \
{ \
  typedef typename newtype::supertype TS; \
  TS returnval(a); \
  returnval opname##= b; \
  return returnval; \
}

#define NumberArray_op(opname, typecomparison) \
NumberArray_op_ab(opname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, \
                  typecomparison##Type<NumberArray<size MacroComma T> MacroComma NumberArray<size MacroComma T2> >) \
NumberArray_op_ab(opname,                             T , NumberArray<size MacroComma T2>, \
                  typecomparison##Type<NumberArray<size MacroComma T2> MacroComma T MacroComma true>) \
NumberArray_op_ab(opname, NumberArray<size MacroComma T>,                             T2 , \
                  typecomparison##Type<NumberArray<size MacroComma T> MacroComma T2>)

NumberArray_op(+,Plus)
NumberArray_op(-,Minus)
NumberArray_op(*,Multiplies)
NumberArray_op(/,Divides)


#define NumberArray_operator_binary_abab(opname, atype, btype, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
inline \
NumberArray<size, bool> \
operator opname (const atype& a, const btype& b) \
{ \
  NumberArray<size, bool> returnval; \
 \
  for (std::size_t i=0; i != size; ++i) \
    returnval[i] = (aarg opname barg); \
 \
  return returnval; \
}

#define NumberArray_operator_binary(opname) \
NumberArray_operator_binary_abab(opname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, a[i], b[i]) \
NumberArray_operator_binary_abab(opname,                             T , NumberArray<size MacroComma T2>, a,    b[i]) \
NumberArray_operator_binary_abab(opname, NumberArray<size MacroComma T>,                             T2 , a[i], b)

NumberArray_operator_binary(<)
NumberArray_operator_binary(<=)
NumberArray_operator_binary(>)
NumberArray_operator_binary(>=)
NumberArray_operator_binary(==)
NumberArray_operator_binary(!=)

template <std::size_t size, typename T>
inline
std::ostream&      
operator<< (std::ostream& output, const NumberArray<size,T>& a)
{
  output << '{';
  if (size)
    output << a[0];
  for (std::size_t i=1; i<size; ++i)
    output << ',' << a[i];
  output << '}';
  return output;
}


// CompareTypes, RawType specializations

#define NumberArray_comparisons(templatename) \
template<std::size_t size, typename T, bool reverseorder> \
struct templatename<NumberArray<size,T>, NumberArray<size,T>, reverseorder> { \
  typedef NumberArray<size, T> supertype; \
}; \
 \
template<std::size_t size, typename T, typename T2, bool reverseorder> \
struct templatename<NumberArray<size,T>, NumberArray<size,T2>, reverseorder> { \
  typedef NumberArray<size, typename Symmetric##templatename<T, T2, reverseorder>::supertype> supertype; \
}; \
 \
template<std::size_t size, typename T, typename T2, bool reverseorder> \
struct templatename<NumberArray<size, T>, T2, reverseorder> { \
  typedef NumberArray<size, typename Symmetric##templatename<T, T2, reverseorder>::supertype> supertype; \
}

NumberArray_comparisons(CompareTypes);
NumberArray_comparisons(PlusType);
NumberArray_comparisons(MinusType);
NumberArray_comparisons(MultipliesType);
NumberArray_comparisons(DividesType);

template <std::size_t size, typename T>
struct RawType<NumberArray<size, T> >
{
  typedef NumberArray<size, typename RawType<T>::value_type> value_type;

  static value_type value(const NumberArray<size, T>& a)
    {
      value_type returnval;
      for (std::size_t i=0; i != size; ++i)
        returnval[i] = RawType<T>::value(a[i]);
      return returnval;
    }
};

} // namespace MetaPhysicL



namespace std {

using MetaPhysicL::NumberArray;
using MetaPhysicL::CompareTypes;

#define NumberArray_std_unary(funcname) \
template <std::size_t size, typename T> \
inline \
NumberArray<size, T> \
funcname (NumberArray<size, T> a) \
{ \
  for (std::size_t i=0; i != size; ++i) \
    a[i] = std::funcname(a[i]); \
 \
  return a; \
}


#define NumberArray_std_binary_abab(funcname, atype, btype, abtypes, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
inline \
typename CompareTypes<abtypes>::supertype \
funcname (const atype& a, const btype& b) \
{ \
  typedef typename CompareTypes<abtypes>::supertype TS; \
  TS returnval; \
 \
  for (std::size_t i=0; i != size; ++i) \
    returnval[i] = std::funcname(aarg, barg); \
 \
  return returnval; \
}

#define NumberArray_std_binary(funcname) \
NumberArray_std_binary_abab(funcname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, \
                            NumberArray<size MacroComma T> MacroComma NumberArray<size MacroComma T2>, a[i], b[i]) \
NumberArray_std_binary_abab(funcname,                             T , NumberArray<size MacroComma T2>, \
                            NumberArray<size MacroComma T2> MacroComma T,                              a,    b[i]) \
NumberArray_std_binary_abab(funcname, NumberArray<size MacroComma T>,                             T2 , \
                            NumberArray<size MacroComma T> MacroComma T2,                              a[i],    b)

NumberArray_std_binary(pow)
NumberArray_std_unary(exp)
NumberArray_std_unary(log)
NumberArray_std_unary(log10)
NumberArray_std_unary(sin)
NumberArray_std_unary(cos)
NumberArray_std_unary(tan)
NumberArray_std_unary(asin)
NumberArray_std_unary(acos)
NumberArray_std_unary(atan)
NumberArray_std_binary(atan2)
NumberArray_std_unary(sinh)
NumberArray_std_unary(cosh)
NumberArray_std_unary(tanh)
NumberArray_std_unary(sqrt)
NumberArray_std_unary(abs)
NumberArray_std_binary(max)
NumberArray_std_binary(min)
NumberArray_std_unary(ceil)
NumberArray_std_unary(floor)
NumberArray_std_binary(fmod)


template <std::size_t size, typename T>
class numeric_limits<NumberArray<size, T> > : 
  public MetaPhysicL::raw_numeric_limits<NumberArray<size, T>, T> {};

} // namespace std


#endif // METAPHYSICL_NUMBERARRAY_H
