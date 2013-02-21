//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// HelloWorld - A template for autotools applications
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

#include "helloworld/helloworld_version.h"

namespace HelloWorld
{

  void helloworld_version_stdout()
  {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "helloworld Package: Version = " << HELLOWORLD_LIB_VERSION;
    std::cout << " (" << get_helloworld_version() << ")" << std::endl << std::endl;
  
    std::cout << HELLOWORLD_LIB_RELEASE << std::endl << std::endl;
  
    std::cout << "Build Date   = " << HELLOWORLD_BUILD_DATE     << std::endl;
    std::cout << "Build Host   = " << HELLOWORLD_BUILD_HOST     << std::endl;
    std::cout << "Build User   = " << HELLOWORLD_BUILD_USER     << std::endl;
    std::cout << "Build Arch   = " << HELLOWORLD_BUILD_ARCH     << std::endl;
    std::cout << "Build Rev    = " << HELLOWORLD_BUILD_VERSION  << std::endl << std::endl;
  
    std::cout << "C++ Config   = " << HELLOWORLD_CXX << " " << HELLOWORLD_CXXFLAGS << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  
    return;
  }

  int get_helloworld_version()
  {
    /* Note: return format follows the versioning convention xx.yy.zz where
   
       xx = major version number
       yy = minor version number
       zz = micro version number
     
       For example:
       v.   0.23  -> 002300 = 2300
       v   0.23.1 -> 002301 = 2301
       v. 10.23.2 -> 102302         */

    int major_version = 0;
    int minor_version = 0;
    int micro_version = 0;

#ifdef HELLOWORLD_MAJOR_VERSION
    major_version = HELLOWORLD_MAJOR_VERSION;
#endif

#ifdef HELLOWORLD_MINOR_VERSION
    minor_version = HELLOWORLD_MINOR_VERSION;
#endif

#ifdef HELLOWORLD_MICRO_VERSION
    micro_version = HELLOWORLD_MICRO_VERSION;
#endif
      
    return major_version*10000 + minor_version*100 + micro_version;
  }

} // end namespace HelloWorld
