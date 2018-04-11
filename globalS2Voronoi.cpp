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
