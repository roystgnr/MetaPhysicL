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
//--------------------------------------------------------------------------

#include "metaphysicl/metaphysicl_exceptions.h"

#include "metaphysicl/metaphysicl_config.h"

#ifdef METAPHYSICL_HAVE_FEENABLEEXCEPT
#include <fenv.h>
#endif

namespace MetaPhysicL
{

  void enableFPE(bool on)
  {
#ifdef METAPHYSICL_HAVE_FEENABLEEXCEPT
    if (on)
      feenableexcept(FE_DIVBYZERO | FE_INVALID);
    else
      fedisableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
    (void)on; // avoid unused var warnings
  }

} // end namespace MetaPhysicL
