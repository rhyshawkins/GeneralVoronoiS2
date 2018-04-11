#pragma once
#ifndef tomographyutil_hpp
#define tomographyutil_hpp

#include <vector>
#include <math.h>

#include "coordinate.hpp"

double greatcircledistkm(double lon1, double lat1,
			 double lon2, double lat2,
			 double r)
{
  sphericalcoordinate<double> c1((90.0 - lat1) * M_PI/180.0, lon1 * M_PI/180.0);
  sphericalcoordinate<double> c2((90.0 - lat2) * M_PI/180.0, lon2 * M_PI/180.0);

  return r * c1.distance(c2);
}

double cartesiandistkm(double lon1, double lat1, double r1,
		       double lon2, double lat2, double r2)
{
  sphericalcoordinate<double> c1((90.0 - lat1) * M_PI/180.0, lon1 * M_PI/180.0);
  sphericalcoordinate<double> c2((90.0 - lat2) * M_PI/180.0, lon2 * M_PI/180.0);

  vector3<double> v1, v2;
  
  sphericalcoordinate<double>::sphericaltocartesian(c1, v1);
  v1 *= r1;

  sphericalcoordinate<double>::sphericaltocartesian(c2, v2);
  v2 *= r2;

  return (v1 - v2).length();
}


void compute_weights(const std::vector<double> &lons,
		     const std::vector<double> &lats,
		     const std::vector<double> &rs,
		     std::vector<double> &weights)
{
  weights.resize(lons.size());
  for (auto &w : weights) {
    w = 0.0;
  }

  for (int i = 1; i < (int)lons.size(); i ++) {

    if (rs[i - 1] == rs[i]) {
      //
      // Great circle path
      //
      double d = greatcircledistkm(lons[i - 1], lats[i - 1],
				   lons[i], lats[i],
				   rs[i]);

      weights[i - 1] += d/2.0;
      weights[i] += d/2.0;

    } else {
      //
      // Cartesian distance
      //
      double d = cartesiandistkm(lons[i - 1], lats[i - 1], rs[i - 1],
				 lons[i], lats[i], rs[i]);
      
      weights[i - 1] += d/2.0;
      weights[i] += d/2.0;

    }
  }
}

#endif // tomographyutil_hpp
