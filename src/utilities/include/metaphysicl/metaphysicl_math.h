// Code created with assistance of ChatGPT-5

// When we can, call C math in the global namespace (device-safe on
// CUDA/HIP). Otherwise make unqualified call to enable ADL (nested
// duals, complexes, etc.). We avoid self-recursion by routing ADL
// through a helper in the nested detail namespace.

#ifndef METAPHYSICL_MATH_H
#define METAPHYSICL_MATH_H

#include "metaphysicl/metaphysicl_device.h"
#ifdef METAPHYSICL_KOKKOS_COMPILATION
#include "metaphysicl/ignore_warnings.h"
#include <Kokkos_Core.hpp>
#include "metaphysicl/restore_warnings.h"
#endif

#include <cmath>
#include <type_traits>

namespace MetaPhysicL {
namespace math {

#define METAPHYSICL_MATH_DETAIL_UNQUALIFIED_UNARY_STD(NAME)                    \
  namespace detail {                                                           \
  template <class X> METAPHYSICL_INLINE auto unqualified_##NAME(const X &x) {  \
    if constexpr (std::is_arithmetic_v<X>)                                     \
      return std::NAME(x);                                                     \
    else {                                                                     \
      using std::NAME;                                                         \
      return NAME(x);                                                          \
    }                                                                          \
  }                                                                            \
  }

#ifdef METAPHYSICL_KOKKOS_COMPILATION
#define METAPHYSICL_MATH_DETAIL_UNQUALIFIED_UNARY(NAME)                        \
  namespace detail {                                                           \
  template <class X> METAPHYSICL_INLINE auto unqualified_##NAME(const X &x) {  \
    if constexpr (std::is_arithmetic_v<X>)                                     \
      return Kokkos::NAME(x);                                                  \
    else {                                                                     \
      using Kokkos::NAME;                                                      \
      return NAME(x);                                                          \
    }                                                                          \
  }                                                                            \
  }
#else
#define METAPHYSICL_MATH_DETAIL_UNQUALIFIED_UNARY(NAME)                        \
  METAPHYSICL_MATH_DETAIL_UNQUALIFIED_UNARY_STD(NAME)
#endif

#ifdef METAPHYSICL_KOKKOS_COMPILATION
#define METAPHYSICL_MATH_DETAIL_UNQUALIFIED_BINARY(NAME)                       \
  namespace detail {                                                           \
  template <class X, class Y>                                                  \
  METAPHYSICL_INLINE auto unqualified_##NAME(const X &x, const Y &y) {         \
    if constexpr (std::is_arithmetic_v<X> && std::is_arithmetic_v<Y>)          \
      return Kokkos::NAME(x, y);                                               \
    else {                                                                     \
      using Kokkos::NAME;                                                      \
      return NAME(x, y);                                                       \
    }                                                                          \
  }                                                                            \
  }
#else
#define METAPHYSICL_MATH_DETAIL_UNQUALIFIED_BINARY(NAME)                       \
  namespace detail {                                                           \
  template <class X, class Y>                                                  \
  METAPHYSICL_INLINE auto unqualified_##NAME(const X &x, const Y &y) {         \
    if constexpr (std::is_arithmetic_v<X> && std::is_arithmetic_v<Y>)          \
      return std::NAME(x, y);                                                  \
    else {                                                                     \
      using std::NAME;                                                         \
      return NAME(x, y);                                                       \
    }                                                                          \
  }                                                                            \
  }
#endif

#define METAPHYSICL_MATH_UNARY(NAME)                                           \
  METAPHYSICL_MATH_DETAIL_UNQUALIFIED_UNARY(NAME)                              \
  template <class X> METAPHYSICL_INLINE auto NAME(const X &x) {                \
    return detail::unqualified_##NAME(x);                                      \
  }

#define METAPHYSICL_MATH_UNARY_STD(NAME)                                       \
  METAPHYSICL_MATH_DETAIL_UNQUALIFIED_UNARY_STD(NAME)                          \
  template <class X> METAPHYSICL_INLINE auto NAME(const X &x) {                \
    return detail::unqualified_##NAME(x);                                      \
  }

// Macro to call when there is a backing math function in c
#define METAPHYSICL_MATH_BINARY(NAME)                                          \
  METAPHYSICL_MATH_DETAIL_UNQUALIFIED_BINARY(NAME)                             \
  template <class X, class Y>                                                  \
  METAPHYSICL_INLINE auto NAME(const X &x, const Y &y) {                       \
    return detail::unqualified_##NAME(x, y);                                   \
  }

METAPHYSICL_MATH_UNARY(sin)
METAPHYSICL_MATH_UNARY(cos)
METAPHYSICL_MATH_UNARY(tanh)
METAPHYSICL_MATH_UNARY(exp)
METAPHYSICL_MATH_UNARY(log)
METAPHYSICL_MATH_UNARY(log2)
METAPHYSICL_MATH_UNARY(sqrt)
METAPHYSICL_MATH_UNARY(asin)
METAPHYSICL_MATH_UNARY(acos)
METAPHYSICL_MATH_UNARY(atan)
METAPHYSICL_MATH_UNARY(sinh)
METAPHYSICL_MATH_UNARY(cosh)
METAPHYSICL_MATH_UNARY(abs)
METAPHYSICL_MATH_UNARY(isnan)
METAPHYSICL_MATH_UNARY(isinf)
METAPHYSICL_MATH_UNARY(log10)
METAPHYSICL_MATH_UNARY(tan)
METAPHYSICL_MATH_UNARY(ceil)
METAPHYSICL_MATH_UNARY(floor)
METAPHYSICL_MATH_UNARY(fabs)
METAPHYSICL_MATH_UNARY(exp2)
METAPHYSICL_MATH_UNARY(expm1)
METAPHYSICL_MATH_UNARY(log1p)
METAPHYSICL_MATH_UNARY(cbrt)
METAPHYSICL_MATH_UNARY(asinh)
METAPHYSICL_MATH_UNARY(acosh)
METAPHYSICL_MATH_UNARY(atanh)
METAPHYSICL_MATH_UNARY(erf)
METAPHYSICL_MATH_UNARY(erfc)
METAPHYSICL_MATH_UNARY(trunc)
METAPHYSICL_MATH_UNARY(round)
#ifndef METAPHYSICL_KOKKOS_COMPILATION
METAPHYSICL_MATH_UNARY(nearbyint)
#endif
METAPHYSICL_MATH_UNARY_STD(rint)
METAPHYSICL_MATH_UNARY(real)
METAPHYSICL_MATH_UNARY(imag)
METAPHYSICL_MATH_UNARY_STD(norm)
METAPHYSICL_MATH_UNARY(conj)
METAPHYSICL_MATH_UNARY(tgamma)
METAPHYSICL_MATH_UNARY(lgamma)

METAPHYSICL_MATH_BINARY(pow)
METAPHYSICL_MATH_BINARY(fmod)
METAPHYSICL_MATH_BINARY(atan2)
METAPHYSICL_MATH_BINARY(hypot)
METAPHYSICL_MATH_BINARY(remainder)
METAPHYSICL_MATH_BINARY(fmax)
METAPHYSICL_MATH_BINARY(fmin)
METAPHYSICL_MATH_BINARY(fdim)
METAPHYSICL_MATH_BINARY(max)
METAPHYSICL_MATH_BINARY(min)

} // namespace math
} // namespace MetaPhysicL

#endif // METAPHYSICL_MATH_H
