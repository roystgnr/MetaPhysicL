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
// $Id: metaprogramming.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef METAPHYSICL_METAPROGRAMMING_H
#define METAPHYSICL_METAPROGRAMMING_H

namespace MetaPhysicL
{
  // Helper metafunctions
  template <bool B, class T = void>
  struct enable_if_c {
    typedef T type;
  };

  template <class T>
  struct enable_if_c<false, T> {};

  template <typename T>
  class has_size
  {
    typedef char no;
    typedef char yes[2];
    template <class C> static yes& test(char (*)[sizeof(&C::size)]);
    template <class C> static no& test(...);
  public:
    const static bool value = (sizeof(test<T>(0)) == sizeof(yes&));
  };

  template <typename T>
  class has_supertype
  {
    typedef char no;
    typedef char yes[2];

    struct Fallback { struct supertype { }; };
    struct Derived : T, Fallback { };

    template < class C >
    static no& test ( typename C::supertype* );
    template < typename C >
    static yes& test ( C* );

  public:
    const static bool value = (sizeof(test<Derived>(0)) == sizeof(yes));
  };

} // end namespace MetaPhysicL

#endif // METAPHYSICL_METAPROGRAMMING_H
