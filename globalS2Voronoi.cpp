
#include "globalS2Voronoi.hpp"

globalS2Voronoi *
globalS2Voronoi::current_state = nullptr;

int
globalS2Voronoi::addobservation(int *npoints,
				double *lons,
				double *lats)
{
  if (current_state != nullptr) {

    current_state->data.add(*npoints, lons, lats);
    return 0;
    
  }

  return -1;
}
