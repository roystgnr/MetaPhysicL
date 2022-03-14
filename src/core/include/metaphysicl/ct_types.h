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


#ifndef METAPHYSICL_CT_TYPES_H
#define METAPHYSICL_CT_TYPES_H

#include <array>
#include <vector>

// Compile-time type functions

namespace MetaPhysicL {

// IfElse takes a boolean condition as the first parameter.
// If the condition is true, then IfElse::type is the same type as parameter 2.
// If the condition is false, then IfElse::type is the same type as parameter 3.
template <bool Condition, typename TrueResult, typename FalseResult>
struct IfElse;

template <typename TrueResult, typename FalseResult>
struct IfElse<true, TrueResult, FalseResult> {
  typedef TrueResult type;
};

template <typename TrueResult, typename FalseResult>
struct IfElse<false, TrueResult, FalseResult> {
  typedef FalseResult type;
};

// TypesEqual takes two types as parameters.
// If they are the exact same type, then TypesEqual::value is the boolean true,
// otherwise TypesEqual::value is the boolean false.
template <typename T1, typename T2>
struct TypesEqual {
  static const bool value = false;
};

template <typename T>
struct TypesEqual<T,T> {
  static const bool value = true;
};

// NullType is used as the equivalent of "NULL" in some compile-time algorithms
struct NullType {
  template <typename T>
  struct rebind { typedef NullType other; };
};

// TrueType is used as the equivalent of boolean true in some compile-time algorithms
struct TrueType
{
  typedef bool value_type;
  static const bool value = true;
};

// FalseType is used as the equivalent of boolean true in some compile-time algorithms
struct FalseType
{
  typedef bool value_type;
  static const bool value = false;
};

/**
 * Structure that can be used for conveniently replacing types in containers. Calling this structure
 * where the first template argument (T) is a container will result in replacement of the inner
 * container type with the second structure template argument (U) resulting in a container of U's.
 * If the first template argument is a nuclear/atomic algebraic type, such as a tensor, vector, or
 * scalar then the resulting type will simply be U
 */
template <typename T, typename U>
struct ReplaceAlgebraicType;

#define METAPHYSICL_BUILTIN_REPLACE_TYPE(Type)  \
  template <typename U>                         \
  struct ReplaceAlgebraicType<Type, U>          \
  {                                             \
    typedef U type;                             \
  }

METAPHYSICL_BUILTIN_REPLACE_TYPE(char);
METAPHYSICL_BUILTIN_REPLACE_TYPE(signed char);
METAPHYSICL_BUILTIN_REPLACE_TYPE(unsigned char);
METAPHYSICL_BUILTIN_REPLACE_TYPE(short int);
METAPHYSICL_BUILTIN_REPLACE_TYPE(unsigned short int);
METAPHYSICL_BUILTIN_REPLACE_TYPE(int);
METAPHYSICL_BUILTIN_REPLACE_TYPE(unsigned int);
METAPHYSICL_BUILTIN_REPLACE_TYPE(long);
METAPHYSICL_BUILTIN_REPLACE_TYPE(long long);
METAPHYSICL_BUILTIN_REPLACE_TYPE(unsigned long);
METAPHYSICL_BUILTIN_REPLACE_TYPE(unsigned long long);
METAPHYSICL_BUILTIN_REPLACE_TYPE(float);
METAPHYSICL_BUILTIN_REPLACE_TYPE(double);
METAPHYSICL_BUILTIN_REPLACE_TYPE(long double);

template <typename T, typename A, typename U>
struct ReplaceAlgebraicType<std::vector<T, A>, U>
{
  typedef std::vector<U, A> type;
};

template <typename T, std::size_t N, typename U>
struct ReplaceAlgebraicType<std::array<T, N>, U>
{
  typedef std::array<U, N> type;
};

} // namespace MetaPhysicL

#endif // METAPHYSICL_CT_TYPES_H
