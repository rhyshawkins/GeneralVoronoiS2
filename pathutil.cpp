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

#include <stdio.h>

#include "pathutil.hpp"

void mkpath(const char *prefix, const char *filename, char *path)
{
  if (prefix != nullptr) {
    sprintf(path, "%s%s", prefix, filename);
  } else {
    sprintf(path, "%s", filename);
  }
}

void mkrankpath(int rank, const char *prefix, const char *filename, char *path)
{
  if (prefix != nullptr) {
    sprintf(path, "%s%s-%03d", prefix, filename, rank);
  } else {
    sprintf(path, "%s-%03d", filename, rank);
  }
}
