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
#ifndef sphericalvoronoimodel_hpp
#define sphericalvoronoimodel_hpp

#include <vector>

#include "coordinate.hpp"

extern "C" {
  #include "slog.h"
};

template
<
  typename value
>
class sphericalvoronoimodel {
public:
  typedef sphericalcoordinate<value> coord_t;

  typedef struct cell {

    cell() :
      v(0.0)
    {
    }

    cell(const coord_t &_c,
	 const value &_v) :
      c(_c),
      v(_v)
    {
    }

    value distance(const coord_t &p) const
    {
      return c.distance(p);
    }
    
    coord_t c;
    value v;
  } cell_t;

  sphericalvoronoimodel(bool _logspace) :
    logspace(_logspace)
  {
  }
  
  ~sphericalvoronoimodel()
  {
  }

  void reset()
  {
    cells.clear();
  }

  int ncells() const
  {
    return cells.size();
  }

  void dump() const
  {
    int i = 0;
    for (auto &c : cells) {
      printf("  %3d: %10.6f ", i, c.v);
      c.c.writetext(stdout);
      printf("\n");
      i ++;
    }
  }
  
  void add_cell(const coord_t &p, const value &v)
  {
    cells.push_back(cell_t(p, v));
  }

  void pop()
  {
    if (cells.size() <= 1) {
      throw GENERALVORONOIS2EXCEPTION("Not enough cells");
    }

    cells.pop_back();
  }

  void delete_cell(int index)
  {
    if (cells.size() <= 1 || index < 0 || index >= (int)cells.size()) {
      throw GENERALVORONOIS2EXCEPTION("Not enough cells or index out of range");
    }

    cells.erase(cells.begin() + index);
  }

  void insert_cell(int index, const coord_t &p, const value &v)
  {
    if (index < 0 || index > (int)cells.size()) {
      throw GENERALVORONOIS2EXCEPTION("Index out of range");
    }

    cells.insert(cells.begin() + index, cell_t(p, v));
  }

  void nearest(const coord_t &p, coord_t &cell_centre, value &cell_value) const
  {
    if (cells.size() == 0) {
      throw GENERALVORONOIS2EXCEPTION("No nodes\n");
    }

    value mindist = cells[0].distance(p);
    cell_value = cells[0].v;
    cell_centre = cells[0].c;
    
    for (int i = 1; i < (int)cells.size(); i ++) {

      value d = cells[i].distance(p);
      if (d < mindist) {
	cell_value = cells[i].v;
	cell_centre = cells[i].c;
	mindist = d;
      }
    }

    if (logspace) {
      cell_value = exp(cell_value);
    }
  }

  value value_at_point(const coord_t &p) const
  {
    coord_t centre;
    value v;

    nearest(p, centre, v);

    return v;
  }

  cell_t *get_cell_by_index(size_t i)
  {
    if (i >= cells.size()) {
      throw GENERALVORONOIS2EXCEPTION("Index out of range %d (%d)", (int)i, (int)cells.size());
    }

    return &(cells[i]);
  }
  
  cell_t &operator[](size_t i) {
    if (i >= cells.size()) {
      throw GENERALVORONOIS2EXCEPTION("Index out of range");
    }

    return cells[i];
  }

  bool save(const char *filename)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }

    fprintf(fp, "%d %d\n", (int)cells.size(), (int)logspace);
    
    for (auto &c : cells) {

      fprintf(fp, "%15.9f %15.9f %.9g\n", c.c.phi, c.c.theta, c.v);

    }

    fclose(fp);

    return true;
  }

  bool load(const char *filename)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      return false;
    }

    int n, l;
    if (fscanf(fp, "%d %d\n", &n, &l) != 2) {
      ERROR("Failed to read ncells and logspace flag");
      return false;
    }

    logspace = (bool)l;

    cells.clear();
    
    for (int i = 0; i < n; i ++) {

      double phi, theta, v;

      if (fscanf(fp, "%lf %lf %lf\n", &phi, &theta, &v) != 3) {
	ERROR("Failed to read cell %d", i);
	return false;
      }

      cells.push_back(cell_t(coord_t(phi, theta), v));
    }

    fclose(fp);
    return true;
  }

  std::vector<cell_t> cells;

  bool logspace;
};

#endif //sphericalvoronoimodel_hpp
  
  
