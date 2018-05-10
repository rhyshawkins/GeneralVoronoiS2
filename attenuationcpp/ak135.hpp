#pragma once
#ifndef ak135_hpp
#define ak135_hpp

#include <vector>

struct ak135entry {
  double depth;   // km
  double density; // Mg/km3
  double Vp;      // km/s
  double Vs;      // km/s
};

extern std::vector<ak135entry> ak135_table;

double ak135_density(double depth);
double ak135_Vp(double depth);
double ak135_Vs(double depth);


#endif // ak135_hpp
