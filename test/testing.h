
#ifndef __testing_h__
#define __testing_h__

// Order of header inclusion matters here: ShadowNumbers need to
// precede DualNumbers which need to precede *Vector

#define USE_SHADOW
#define USE_SPARSE
#define USE_STRUCT

const unsigned int NDIM = 2;

#define TESTPRINT(a) do { std::cout << #a " = " << a << std::endl; } while (0)

#ifdef USE_SPARSE
  #ifdef USE_SHADOW
    #ifdef USE_STRUCT
      #include "dualshadowsparsestruct.h"
    #else // USE_STRUCT
      #include "dualshadowsparsevector.h"
    #endif // USE_STRUCT
    typedef ShadowNumber<double, long double> RawScalar;
  #else // USE_SHADOW
    #ifdef USE_STRUCT
      #include "dualsparsenumberstruct.h"
    #else // USE_STRUCT
      #include "dualsparsenumbervector.h"
    #endif // USE_STRUCT
    typedef double RawScalar;
  #endif // USE_SHADOW
  #ifdef USE_STRUCT
    typedef CompileTimeContainer::SetConstructor<CompileTimeContainer::UnsignedIntType<0,RawScalar>, CompileTimeContainer::UnsignedIntType<1,RawScalar> >::type IndexSet;
    typedef SparseNumberStruct<IndexSet> RawVector;
    #define VectorUnitVector SparseNumberStructUnitVector
    #define VectorOf SparseNumberStructOf
  #else // USE_STRUCT
    typedef CompileTimeContainer::SetConstructor<CompileTimeContainer::UnsignedIntType<0>, CompileTimeContainer::UnsignedIntType<1> >::type IndexSet;
    typedef SparseNumberVector<RawScalar, IndexSet> RawVector;
    #define VectorUnitVector SparseNumberVectorUnitVector
    #define VectorOf SparseNumberVectorOf
  #endif // USE_STRUCT
#else // USE_SPARSE
  #ifdef USE_SHADOW
    #include "dualshadowvector.h"
    typedef ShadowNumber<double, long double> RawScalar;
  #else // USE_SHADOW
    #include "dualnumbervector.h"
    typedef double RawScalar;
  #endif // USE_SHADOW
  typedef NumberVector<NDIM, RawScalar> RawVector;
  #define VectorUnitVector NumberVectorUnitVector
  #define VectorOf NumberVectorOf
#endif // USE_SPARSE

#endif // __testing_h__
