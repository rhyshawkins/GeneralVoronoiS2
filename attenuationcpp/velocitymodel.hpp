#pragma once
#ifndef velocitymodel_hpp
#define velocitymodel_hpp

#include "ak135.hpp"

template
<typename value>
value pwave_velocity(value r)
{
  //
  // Example polynomial (put in your one and uncomment)
  //
  return ak135_Vp(6371.0 - r);

}


#endif // velocitymodel_hpp
