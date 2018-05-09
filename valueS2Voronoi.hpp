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
#ifndef valueS2Voronoi_hpp
#define valueS2Voronoi_hpp

#include <mpi.h>

#include <string>

#include "globalS2Voronoi.hpp"
#include "perturbationS2Voronoi.hpp"

class ValueS2Voronoi : public PerturbationS2Voronoi {
public:

  typedef sphericalcoordinate<double> coord_t;
  typedef typename sphericalvoronoimodel<double>::cell_t cell_t;
  typedef deltaVoronoi<coord_t, double> delta_t;
  typedef model_deltaVoronoi<coord_t, double> model_delta_t;
  
  ValueS2Voronoi() :
    undo_cell(nullptr),
    undo_v(0.0),
    last_log_proposal_ratio(0.0),
    p(0),
    a(0)
  {
  }
  
  ~ValueS2Voronoi()
  {
  }

  virtual bool propose(int maxcells,
		       int nobs,
		       Rng &random,
		       PriorProposal &prior,
		       SphericalPriorProposal &position_prior,
		       sphericalvoronoimodel<double> &model,
		       PriorProposal &hierarchical_prior,
		       hierarchical_model &hierarchical,
		       double temperature,
		       double &log_prior_ratio,
		       delta_t *&perturbation,
		       bool &relocate)
  {
    bool validproposal = false;
    int cell = -1;
    double oldv = 0.0;
    double newv;

    relocate = false;
    
    if (this->primary()) {

      //
      // Primary does generation of new value
      //
      
      p ++;
      cell = random.uniform(model.ncells());
      
      //
      // Next get the active node
      //
      cell_t *c = model.get_cell_by_index(cell);
      
      oldv = c->v;
    
      if (prior.propose(random, temperature, oldv, newv, log_prior_ratio)) {
	
	perturbation = model_delta_t::mkvalue(cell, oldv, newv);
	
	validproposal = true;
	
      } else {
	perturbation = model_delta_t::mkvalue(cell, oldv, 0.0);
      }
	
    }

    this->communicate(validproposal);

    if (validproposal) {
      this->communicate(cell);
      this->communicate(newv);

      cell_t *c = model.get_cell_by_index(cell);
      undo_v = c->v;
      undo_cell = c;

      c->v = newv;
      
      last_log_proposal_ratio = prior.log_proposal_ratio(random, temperature, undo_v, newv);
    }
    
    return validproposal;
  }

  virtual double log_proposal_ratio(Rng &random,
				    PriorProposal &prior,
				    SphericalPriorProposal &position_prior,
				    sphericalvoronoimodel<double> &proposed_model,
				    PriorProposal &hierarchical_prior,
				    hierarchical_model &proposed_hierarchical,
				    double temperature)
  {
    return last_log_proposal_ratio;
  }

  void accept()
  {
    a ++;
    
    if (undo_cell == nullptr) {
      throw GENERALVORONOIS2EXCEPTION("No undo information\n");
    }
    
    undo_cell = nullptr;
    undo_v = 0.0;
  }
  
  void reject(sphericalvoronoimodel<double> &model)
  {
    if (undo_cell == nullptr) {
      throw GENERALVORONOIS2EXCEPTION("No undo information\n");
    }
    
    undo_cell->v = undo_v;
    
    undo_cell = nullptr;
    undo_v = 0.0;
  }
  
  virtual int proposal_count() const
  {
    return p;
  }
  
  virtual int acceptance_count() const
  {
    return a;
  }

  virtual const char *displayname() const
  {
    return "Value";
  }
  
private:

  cell_t *undo_cell;
  double undo_v;

  double last_log_proposal_ratio;

  int p;
  int a;
  
};

#endif // valueS2Voronoi_hpp
