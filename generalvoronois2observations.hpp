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
