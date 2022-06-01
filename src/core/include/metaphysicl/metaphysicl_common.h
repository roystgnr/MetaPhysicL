#ifndef METAPHYSICL_COMMON_H
#define METAPHYSICL_COMMON_H

// The libmesh_dbg_var() macro indicates that an argument to a function
// is used only in debug mode (i.e., when NDEBUG is not defined).
#ifndef NDEBUG
#define metaphysicl_dbg_var(var) var
#else
#define metaphysicl_dbg_var(var)
#endif

// If we have C++11, we'd like to define as many operators and
// functions as possible with rvalue reference parameters, so as to
// optimize heap access.
#if __cplusplus >= 201103L

// But if we have e.g. clang 9, then we've seen SIGFPE triggered when
// using optimizations on rvalue references, in code that raises no
// warnings and no valgrind errors from -O0, or from any configuration
// of five other compilers ... so we're just going to make things as
// easy on clang as possible, and if we want more optimized code we're
// going to want a better compiler.
//
// I don't have clang++-10 handy to check with, so I'm not going to
// trust it either.
#ifndef __clang__
#define METAPHYSICL_USE_STD_MOVE 1
#else
#  if __clang_major__ > 10
#  define METAPHYSICL_USE_STD_MOVE 1
#  else
#  undef METAPHYSICL_USE_STD_MOVE
#  endif

#endif

#endif

#endif // METAPHYSICL_COMMON_H
