#pragma once
#ifndef hierarchicalS2Voronoi_hpp
#define hierarchicalS2Voronoi_hpp

#include <string>

#include "globalS2Voronoi.hpp"
#include "perturbationS2Voronoi.hpp"

class HierarchicalS2Voronoi : public PerturbationS2Voronoi {
public:

  typedef sphericalcoordinate<double> coord_t;
  typedef deltaVoronoi<coord_t, double> delta_t;
  typedef model_deltaVoronoi<coord_t, double> model_delta_t;

  HierarchicalS2Voronoi() :
    undo_h(nullptr),
    undo_index(-1),
    undo_v(0.0),
    p(0),
    a(0)
  {
  }
  
  ~HierarchicalS2Voronoi()
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
		       delta_t *&perturbation)
  {
    bool validproposal = false;
    int hindex;
    double newv;

    if (undo_h != nullptr) {
      throw GENERALVORONOIS2EXCEPTION("Undo information existing for new proposal\n");
    }

    if (this->primary()) {
      
    
      p ++;
      
      hindex = 0;
      if (hierarchical.get_nhierarchical() > 1) {
	hindex = random.uniform(hierarchical.get_nhierarchical());
      }
    
      double oldv = hierarchical.get(hindex);
      if (hierarchical_prior.propose(random, temperature, oldv, newv, log_prior_ratio)) {
	
	perturbation = new hierarchical_deltaVoronoi<coord_t, double>(1, &hindex, &oldv, &newv);
	validproposal = true;
	
      } else {
	perturbation = new hierarchical_deltaVoronoi<coord_t, double>(1, &hindex, &oldv, &newv);
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      this->communicate(hindex);
      this->communicate(newv);
			
      undo_h = &hierarchical;
      undo_index = hindex;
      undo_v = hierarchical.get(hindex);
	
      log_prior_ratio += 
	(double)nobs*(log(undo_v) - log(newv));
	
      hierarchical.set(hindex, newv);
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
    return 0.0;
  }

  virtual void accept()
  {
    a ++;
    
    if (undo_h == nullptr) {
      throw GENERALVORONOIS2EXCEPTION("No undo information\n");
    }
    
    undo_h = nullptr;
    undo_index = -1;
    undo_v = 0.0;
  }    

  virtual void reject(sphericalvoronoimodel<double> &model)
  {
    if (undo_h == nullptr) {
      throw GENERALVORONOIS2EXCEPTION("No undo information\n");
    }
    
    undo_h->set(undo_index, undo_v);
    
    undo_h = nullptr;
    undo_index = -1;
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
    return "Hierarchical";
  }
  
private:

  hierarchical_model *undo_h;
  int undo_index;
  double undo_v;

  int p;
  int a;
  
  
};

#endif // hierarchicalS2_hpp
