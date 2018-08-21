#ifndef METAPHYSICL_NDDUALNUMBER_H
#define METAPHYSICL_NDDUALNUMBER_H

#include "metaphysicl/dualnumber_decl.h"

namespace MetaPhysicL {

template <typename T, typename D, typename Enable>
class NotADuckDualNumber : public DualNumber<T, D>
{
public:
  using DualNumber<T, D>::DualNumber;

  NotADuckDualNumber<T, D> operator-() const { return NotADuckDualNumber<T, D>(-this->value(), -this->derivatives()); }
  NotADuckDualNumber<T, D> operator!() const { return NotADuckDualNumber<T, D>(!this->value(), !this->derivatives()); }
};

template <typename T, typename D, typename Enable>
using NDDualNumber = NotADuckDualNumber<T,D,Enable>;

}

#endif
