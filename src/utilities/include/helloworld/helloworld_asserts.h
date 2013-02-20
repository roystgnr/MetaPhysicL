//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// HelloWorld - An Autotools library template
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
//--------------------------------------------------------------------------

#ifndef HELLOWORLD_ASSERTS_H
#define HELLOWORLD_ASSERTS_H

// C++
#include <iostream>
#include <iomanip>

// HelloWorld
#include "helloworld/helloworld_exceptions.h"

#define helloworld_here()     do { std::cerr << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; } while (0)

// The helloworld_assert() macro acts like C's assert(), but throws a
// helloworld_error() (including stack trace, etc) instead of just exiting
#ifdef NDEBUG
#define helloworld_assert(asserted)  ((void) 0)
#define helloworld_assert_msg(asserted, msg)  ((void) 0)
#define helloworld_assert_equal_to(expr1,expr2)  ((void) 0)
#define helloworld_assert_not_equal_to(expr1,expr2)  ((void) 0)
#define helloworld_assert_less(expr1,expr2)  ((void) 0)
#define helloworld_assert_greater(expr1,expr2)  ((void) 0)
#define helloworld_assert_less_equal(expr1,expr2)  ((void) 0)
#define helloworld_assert_greater_equal(expr1,expr2)  ((void) 0)
#else
#define helloworld_assert(asserted)  do { if (!(asserted)) { std::cerr << "Assertion `" #asserted "' failed." << std::endl; helloworld_error(); } } while(0)
#define helloworld_assert_equal_to(expr1,expr2)  do { if (!(expr1 == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; helloworld_error(); } } while(0)
#define helloworld_assert_not_equal_to(expr1,expr2)  do { if (!(expr1 != expr2)) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; helloworld_error(); } } while(0)
#define helloworld_assert_less(expr1,expr2)  do { if (!(expr1 < expr2)) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; helloworld_error(); } } while(0)
#define helloworld_assert_greater(expr1,expr2)  do { if (!(expr1 > expr2)) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; helloworld_error(); } } while(0)
#define helloworld_assert_less_equal(expr1,expr2)  do { if (!(expr1 <= expr2)) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; helloworld_error(); } } while(0)
#define helloworld_assert_greater_equal(expr1,expr2)  do { if (!(expr1 >= expr2)) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; helloworld_error(); } } while(0)
#endif


#define helloworld_error()    do { helloworld_here(); HELLOWORLD_THROW(HelloWorld::LogicError()); } while(0)
#define helloworld_not_implemented()    do { helloworld_here(); HELLOWORLD_THROW(HelloWorld::NotImplemented()); } while(0)
#define helloworld_file_error(filename)    do { helloworld_here(); HELLOWORLD_THROW(HelloWorld::FileError(filename)); } while(0)

#endif // HELLOWORLD_ASSERTS_H
