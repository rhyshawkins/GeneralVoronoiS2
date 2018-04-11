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
#ifndef perturbationcollectionS2Voronoi_hpp
#define perturbationcollectionS2_hpp

#include "globalS2Voronoi.hpp"
#include "chainhistoryVoronoi.hpp"
#include "perturbationS2Voronoi.hpp"

class PerturbationCollectionS2Voronoi {
public:

  typedef sphericalcoordinate<double> coord_t;
  typedef deltaVoronoi<coord_t, double> delta_t;
  
  
  PerturbationCollectionS2Voronoi() :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1),
    weight_sum(0.0),
    active(-1)
  {
  }
  
  ~PerturbationCollectionS2Voronoi()
  {
    for (auto &i : perturbations) {
      delete i.p;
    }
  }

  void initialize_mpi(MPI_Comm _communicator)
  {
    MPI_Comm_dup(_communicator, &communicator);
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
  }
  
  void add(PerturbationS2Voronoi *p, double weight)
  {
    if (weight <= 0.0) {
      throw GENERALVORONOIS2EXCEPTION("Invalid weight: %f\n", weight);
    }
    
    weight_sum += weight;
    perturbations.push_back(WeightedPerturbation(p, weight));
  }

  bool propose(globalS2Voronoi &g, double &log_prior_ratio, delta_t *&perturbation)
  {
    if (active >= 0) {
      throw GENERALVORONOIS2EXCEPTION("Already have active perturbation\n");
    }
    
    if (primary()) {
      double u = g.random.uniform();
      active = 0;
      for (auto &wp: perturbations) {
	double w = wp.w/weight_sum;
	
	if (u < w) {
	  break;
	} else {
	  active ++;
	  u -= w;
	}
      }
    
      if (active >= (int)perturbations.size()) {
	throw GENERALVORONOIS2EXCEPTION("Failed to determine active perturbation\n");
      }
    }

    communicate(active);
    
    WeightedPerturbation &wp = perturbations[active];
    
    int nobs = g.data.obs.size();

    bool r = wp.p->propose(g.maxcells,
			   nobs,
			   g.random,
			   *g.prior,
			   *g.positionprior,
			   *g.model,
			   *g.hierarchicalprior,
			   *g.hierarchical,
			   g.temperature,
			   log_prior_ratio,
			   perturbation);
    if (!r) {
      active = -1;
    }
    
    return r;
  }

  double log_proposal_ratio(globalS2Voronoi &g) 
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw GENERALVORONOIS2EXCEPTION("Invalid active perturbation %d (%d)\n", active, (int)perturbations.size());
    }
    
    WeightedPerturbation &wp = perturbations[active];

    return wp.p->log_proposal_ratio(g.random,
				    *g.prior,
				    *g.positionprior,
				    *g.model,
				    *g.hierarchicalprior,
				    *g.hierarchical,
				    1.0);
  }

  void accept(globalS2Voronoi &g)
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw GENERALVORONOIS2EXCEPTION("Invalid active perturbation %d (%d)\n", active, (int)perturbations.size());
    }
    
    WeightedPerturbation &wp = perturbations[active];
    wp.p->accept();
    
    active = -1;
  }
  

  void reject(globalS2Voronoi &g)
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw GENERALVORONOIS2EXCEPTION("Invalid active perturbation %d\n", active);
    }
    
    WeightedPerturbation &wp = perturbations[active];
    wp.p->reject(*g.model);
    
    active = -1;
  }
  

  void writeacceptancereport(FILE *fp)
  {
    fprintf(fp, "%s", generateacceptancereport().c_str());
  }

  std::string generateacceptancereport()
  {
    std::string s;
    char linebuffer[1024];
    
    for (auto &wp: perturbations) {
      int p = wp.p->proposal_count();
      int a = wp.p->acceptance_count();
      
      double f = 0.0;
      if (p > 0) {
	f = (double)a/(double)p * 100.0;
      }
      
      sprintf(linebuffer, "  %12s: %6d %6d : %6.2f\n", wp.p->displayname(), p, a, f);

      s += linebuffer;
    }

    return s;
  }
	       
private:

  bool primary()
  {
    return (communicator == MPI_COMM_NULL) || (rank == 0);
  }

  void communicate(int &i)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&i, 1, MPI_INT, 0, communicator);
    }
  }

  void communicate(double &d)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&d, 1, MPI_DOUBLE, 0, communicator);
    }
  }
  

  struct WeightedPerturbation {
    WeightedPerturbation(PerturbationS2Voronoi *_p, double _w) :
      p(_p),
      w(_w)
    {
    }
    
    PerturbationS2Voronoi *p;
    double w;
  };

  MPI_Comm communicator;
  int rank;
  int size;
  
  double weight_sum;
  std::vector<WeightedPerturbation> perturbations;
  int active;

};

#endif // perturbationcollectionS2Voronoi_hpp
