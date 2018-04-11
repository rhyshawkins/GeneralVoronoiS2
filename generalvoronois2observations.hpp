//
//    GeneralVoronoiS2 : A general Trans-dimensional Voronoic Cell program
//    for surface spherical problems.
//
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#pragma once
#ifndef generalvoronois2observations_hpp
#define generalvoronois2observations_hpp

#include <vector>
#include "coordinate.hpp"

typedef double (*synthetic_model_f)(double phi, double theta);

class GeneralVoronoiS2Observations {
public:

  struct point {

    sphericalcoordinate<double> p;

    double pred;
    double weight;
    int idx;
  };
  
  struct observation {
    double pred;
    double res;

    std::vector<sphericalcoordinate<double>> points;

    std::vector<double> values;
    std::vector<double> weights;
    std::vector<int> idx;
  };

  
  GeneralVoronoiS2Observations()
  {
  }

  void add(int npoints, const double *lon, const double *lat)
  {
    int n = obs.size();
    obs.push_back(observation());

    for (int i = 0; i < npoints; i ++) {
    
      double phi = (90.0 - lat[i]) * M_PI/180.0;
      double theta = lon[i] * M_PI/180.0;

      obs[n].points.push_back(sphericalcoordinate<double>(phi,theta));
    }

    obs[n].pred = 0.0;
    obs[n].res = 0.0;

    obs[n].values.resize(npoints);
    obs[n].weights.resize(npoints);
    obs[n].idx.resize(npoints);
  }

  std::vector<observation> obs;
};

#endif // generalvoronois2observations_hpp
