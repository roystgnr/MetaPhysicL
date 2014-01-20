
#ifndef __testing_h__
#define __testing_h__

// Order of header inclusion matters here: ShadowNumbers need to
// precede DualNumbers which need to precede *Vector

const unsigned int NDIM = 2;

#define TESTPRINT(a) do { std::cout << #a " = " << a << std::endl; } while (0)

#ifdef USE_SPARSE
  #ifdef USE_SHADOW
    #ifdef USE_STRUCT
      #include "metaphysicl/dualshadowsparsestruct.h"
    #else // USE_STRUCT
      #include "metaphysicl/dualshadowsparsevector.h"
    #endif // USE_STRUCT
    typedef MetaPhysicL::ShadowNumber<double, long double> RawScalar;
  #else // USE_SHADOW
    #ifdef USE_STRUCT
      #include "metaphysicl/dualsparsenumberstruct.h"
    #else // USE_STRUCT
      #include "metaphysicl/dualsparsenumbervector.h"
    #endif // USE_STRUCT
    typedef double RawScalar;
  #endif // USE_SHADOW
  #ifdef USE_STRUCT
    typedef
    MetaPhysicL::SetConstructor<MetaPhysicL::UnsignedIntType<0,RawScalar>, MetaPhysicL::UnsignedIntType<1,RawScalar> >::type IndexSet;
    typedef MetaPhysicL::SparseNumberStruct<IndexSet> RawVector;
    #define VectorUnitVector SparseNumberStructUnitVector
    #define VectorOf SparseNumberStructOf
  #else // USE_STRUCT
    typedef MetaPhysicL::SetConstructor<MetaPhysicL::UnsignedIntType<0>, MetaPhysicL::UnsignedIntType<1> >::type IndexSet;
    typedef MetaPhysicL::SparseNumberVector<RawScalar, IndexSet> RawVector;
    #define VectorUnitVector SparseNumberVectorUnitVector
    #define VectorOf SparseNumberVectorOf
  #endif // USE_STRUCT
#else // USE_SPARSE
  #ifdef USE_SHADOW
    #include "metaphysicl/dualshadowvector.h"
    typedef MetaPhysicL::ShadowNumber<double, long double> RawScalar;
  #else // USE_SHADOW
    #include "metaphysicl/dualnumbervector.h"
    typedef double RawScalar;
  #endif // USE_SHADOW
  typedef MetaPhysicL::NumberVector<NDIM, RawScalar> RawVector;
  #define VectorUnitVector NumberVectorUnitVector
  #define VectorOf NumberVectorOf
#endif // USE_SPARSE

#endif // __testing_h__
