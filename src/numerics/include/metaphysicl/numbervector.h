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


#ifndef METAPHYSICL_NUMBERVECTOR_H
#define METAPHYSICL_NUMBERVECTOR_H

#include <algorithm>
#include <ostream>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/raw_type.h"

namespace MetaPhysicL {

template <std::size_t size, typename T>
class NumberVector;

template<std::size_t size, typename S, typename T, bool reverseorder>
struct DotType<NumberVector<size,S>, NumberVector<size,T>, reverseorder> {
  typedef typename MultipliesType<S,T,reverseorder>::supertype supertype;
};

template<std::size_t size1, std::size_t size2, typename S, typename T, bool reverseorder>
struct OuterProductType<NumberVector<size1,S>, NumberVector<size2,T>, reverseorder> {
  typedef 
  NumberVector<size1, NumberVector<size2,
    typename MultipliesType<S,T,reverseorder>::supertype> > supertype;
};

template<std::size_t size, typename S>
struct SumType<NumberVector<size,S> > {
  typedef S supertype;
};

template <std::size_t size, typename T>
class NumberVector
{
public:
  typedef T value_type;

  template <std::size_t i>
  struct entry_type {
    typedef value_type type;
  };

  template <typename T2>
  struct rebind {
    typedef NumberVector<size, T2> other;
  };

  NumberVector() {}

  NumberVector(const T& val)
    { std::fill(_data, _data+size, val); }

  NumberVector(const T* vals)
    { std::copy(vals, vals+size, _data); }

  template <typename T2>
  NumberVector(NumberVector<size, T2> src)
    { if (size) std::copy(&src[0], &src[0]+size, _data); }

  template <typename T2>
  NumberVector(const T2& val)
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

  NumberVector<size,T> operator- () const {
    NumberVector<size,T> returnval;
    for (std::size_t i=0; i != size; ++i) returnval[i] = -_data[i];
    return returnval;
  }

  template <typename T2>
  NumberVector<size,T>& operator+= (const NumberVector<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] += a[i]; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator+= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] += a; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator-= (const NumberVector<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] -= a[i]; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator-= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] -= a; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator*= (const NumberVector<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] *= a[i]; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator*= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] *= a; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator/= (const NumberVector<size,T2>& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] /= a[i]; return *this; }

  template <typename T2>
  NumberVector<size,T>& operator/= (const T2& a)
    { for (std::size_t i=0; i != size; ++i) _data[i] /= a; return *this; }

  template <typename T2>
  typename SymmetricMultipliesType<T,T2>::supertype
  dot (const NumberVector<size,T2>& a) const
  {
    typename SymmetricMultipliesType<T,T2>::supertype returnval = 0;
    for (std::size_t i=0; i != size; ++i)
      returnval += _data[i] * a[i];
    return returnval;
  }

  template <typename T2>
  NumberVector<size, NumberVector<size, typename SymmetricMultipliesType<T,T2>::supertype> >
  outerproduct (const NumberVector<size,T2>& a) const
  {
    NumberVector<size, NumberVector<size, typename SymmetricMultipliesType<T,T2>::supertype> > returnval;

    for (std::size_t i=0; i != size; ++i)
      for (std::size_t j=0; j != size; ++j)
        returnval[i][j] = _data[i] * a[j];

    return returnval;
  }

  static NumberVector<size, NumberVector<size, T> > identity()
  {
    NumberVector<size, NumberVector<size, T> > returnval(0);
  
    for (std::size_t i=0; i != size; ++i)
      returnval[i][i] = 1;

    return returnval;
  }

  T sum () const
  {
    T returnval = 0;
    
    for (std::size_t i=0; i != size; ++i)
      returnval += _data[i];

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
struct NumberVectorOf
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

  typedef NumberVector<size, supertype> type;
};



template <std::size_t size, std::size_t index, typename T>
struct NumberVectorUnitVector
{
  typedef NumberVector<size, T> type;

  static const type value() {
    type returnval = 0;
    returnval[index] = 1;
    return returnval;
  }
};


template <std::size_t size, typename T>
struct NumberVectorFullVector
{
  typedef NumberVector<size,T> type;

  static const type value() {
    type returnval;
    for (std::size_t i=0; i != size; ++i)
      returnval[i] = 1;
    return returnval;
  }
};



template <std::size_t size, typename T>
inline
NumberVector<size, NumberVector<size, T> >
transpose(NumberVector<size, NumberVector<size, T> > a)
{
  for (std::size_t i=0; i != size; ++i)
    for (std::size_t j=i+1; j != size; ++j)
      std::swap(a[i][j], a[j][i]);

  return a;
}



#define NumberVector_op_ab(opname, atype, btype, newtype) \
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

#define NumberVector_op(opname, typecomparison) \
NumberVector_op_ab(opname, NumberVector<size MacroComma T>, NumberVector<size MacroComma T2>, \
                  typecomparison##Type<NumberVector<size MacroComma T> MacroComma NumberVector<size MacroComma T2> >) \
NumberVector_op_ab(opname,                             T , NumberVector<size MacroComma T2>, \
                  typecomparison##Type<NumberVector<size MacroComma T2> MacroComma T MacroComma true>) \
NumberVector_op_ab(opname, NumberVector<size MacroComma T>,                             T2 , \
                  typecomparison##Type<NumberVector<size MacroComma T> MacroComma T2>)

NumberVector_op(+,Plus)
NumberVector_op(-,Minus)
NumberVector_op(*,Multiplies)
NumberVector_op(/,Divides)

namespace std {

#define NumberVector_std_unary(funcname) \
template <std::size_t size, typename T> \
inline \
NumberVector<size, T> \
funcname (NumberVector<size, T> a) \
{ \
  for (std::size_t i=0; i != size; ++i) \
    a[i] = std::funcname(a[i]); \
 \
  return a; \
}


#define NumberVector_std_binary_abab(funcname, atype, btype, abtypes, aarg, barg) \
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

#define NumberVector_std_binary(funcname) \
NumberVector_std_binary_abab(funcname, NumberVector<size MacroComma T>, NumberVector<size MacroComma T2>, \
                            NumberVector<size MacroComma T> MacroComma NumberVector<size MacroComma T2>, a[i], b[i]) \
NumberVector_std_binary_abab(funcname,                             T , NumberVector<size MacroComma T2>, \
                            NumberVector<size MacroComma T2> MacroComma T,                              a,    b[i]) \
NumberVector_std_binary_abab(funcname, NumberVector<size MacroComma T>,                             T2 , \
                            NumberVector<size MacroComma T> MacroComma T2,                              a[i],    b)

NumberVector_std_binary(pow)
NumberVector_std_unary(exp)
NumberVector_std_unary(log)
NumberVector_std_unary(log10)
NumberVector_std_unary(sin)
NumberVector_std_unary(cos)
NumberVector_std_unary(tan)
NumberVector_std_unary(asin)
NumberVector_std_unary(acos)
NumberVector_std_unary(atan)
NumberVector_std_binary(atan2)
NumberVector_std_unary(sinh)
NumberVector_std_unary(cosh)
NumberVector_std_unary(tanh)
NumberVector_std_unary(sqrt)
NumberVector_std_unary(abs)
NumberVector_std_binary(max)
NumberVector_std_binary(min)
NumberVector_std_unary(ceil)
NumberVector_std_unary(floor)
NumberVector_std_binary(fmod)


template <std::size_t size, typename T>
class numeric_limits<NumberVector<size, T> > : 
  public raw_numeric_limits<NumberVector<size, T>, T> {};

} // namespace std

#define NumberVector_operator_binary_abab(opname, atype, btype, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
inline \
NumberVector<size, bool> \
operator opname (const atype& a, const btype& b) \
{ \
  NumberVector<size, bool> returnval; \
 \
  for (std::size_t i=0; i != size; ++i) \
    returnval[i] = (aarg opname barg); \
 \
  return returnval; \
}

#define NumberVector_operator_binary(opname) \
NumberVector_operator_binary_abab(opname, NumberVector<size MacroComma T>, NumberVector<size MacroComma T2>, a[i], b[i]) \
NumberVector_operator_binary_abab(opname,                             T , NumberVector<size MacroComma T2>, a,    b[i]) \
NumberVector_operator_binary_abab(opname, NumberVector<size MacroComma T>,                             T2 , a[i], b)

NumberVector_operator_binary(<)
NumberVector_operator_binary(<=)
NumberVector_operator_binary(>)
NumberVector_operator_binary(>=)
NumberVector_operator_binary(==)
NumberVector_operator_binary(!=)

template <std::size_t size, typename T>
inline
std::ostream&      
operator<< (std::ostream& output, const NumberVector<size,T>& a)
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

#define NumberVector_comparisons(templatename) \
template<std::size_t size, typename T, bool reverseorder> \
struct templatename<NumberVector<size,T>, NumberVector<size,T>, reverseorder> { \
  typedef NumberVector<size, T> supertype; \
}; \
 \
template<std::size_t size, typename T, typename T2, bool reverseorder> \
struct templatename<NumberVector<size,T>, NumberVector<size,T2>, reverseorder> { \
  typedef NumberVector<size, typename Symmetric##templatename<T, T2, reverseorder>::supertype> supertype; \
}; \
 \
template<std::size_t size, typename T, typename T2, bool reverseorder> \
struct templatename<NumberVector<size, T>, T2, reverseorder, \
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> { \
  typedef NumberVector<size, typename Symmetric##templatename<T, T2, reverseorder>::supertype> supertype; \
}

NumberVector_comparisons(CompareTypes);
NumberVector_comparisons(PlusType);
NumberVector_comparisons(MinusType);
NumberVector_comparisons(MultipliesType);
NumberVector_comparisons(DividesType);

template <std::size_t size, typename T>
struct RawType<NumberVector<size, T> >
{
  typedef NumberVector<size, typename RawType<T>::value_type> value_type;

  static value_type value(const NumberVector<size, T>& a)
    {
      value_type returnval;
      for (std::size_t i=0; i != size; ++i)
        returnval[i] = RawType<T>::value(a[i]);
      return returnval;
    }
};

} // namespace MetaPhysicL

#endif // METAPHYSICL_NUMBERVECTOR_H
