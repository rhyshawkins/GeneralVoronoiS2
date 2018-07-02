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
#ifndef globalspherical_hpp
#define globalspherical_hpp

#include <mpi.h>

#include "sphericalvoronoimodel.hpp"

#include "generalvoronois2observations.hpp"
#include "prior.hpp"
#include "sphericalprior.hpp"

#include "hierarchical_model.hpp"

#include "genericinterface.hpp"

#include "simulated_annealing_scales.hpp"

extern "C" {
  #include "slog.h"
};

class globalS2Voronoi {
public:

  typedef sphericalcoordinate<double> coord_t;

  globalS2Voronoi(const char *input,
		  const char *prior_file,
		  const char *hierarchical_prior_file,
		  const char *position_prior_file,
		  const char *birthdeath_proposal_file,
		  int _maxcells,
		  double _lambda,
		  double _temperature,
		  int seed,
		  bool posterior,
		  bool logspace) :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1),
    mpi_counts(nullptr),
    mpi_offsets(nullptr),
    model(nullptr),
    prior(nullptr),
    positionprior(nullptr),
    hierarchicalprior(nullptr),
    birthdeathvalueproposal(nullptr),
    birthdeathpositionproposal(nullptr),
    hierarchical(new singlescaling_hierarchical_model(_lambda)),
    temperature(_temperature),
    residual_size(0),
    mean_residual_n(0),
    maxcells(_maxcells),
    random(seed)
  {

    if (prior_file == nullptr) {
      Prior *p = new UniformPrior(-1.5, 1.5);
      Proposal *pp = new GaussianProposal(*p, 0.05);
      prior = new PriorProposal(p, pp);
    } else {
      prior = PriorProposal::load(prior_file);
      if (prior == nullptr) {
        throw GENERALVORONOIS2EXCEPTION("Failed to open prior/proposal file\n");
      }
    }

    if (hierarchical_prior_file == nullptr) {
      Prior *p = new UniformPrior(0.1, 5.0);
      Proposal *pp = new GaussianProposal(*p, 0.05);
      hierarchicalprior = new PriorProposal(p, pp);
    } else {
      hierarchicalprior = PriorProposal::load(hierarchical_prior_file);
      if (hierarchicalprior == nullptr) {
	throw GENERALVORONOIS2EXCEPTION("Failed to open hierarchical prior/proposal file\n");
      }
    }

    if (position_prior_file == nullptr) {
      SphericalPrior *p = new UniformSphericalPrior();
      SphericalProposal *pp = new VonMisesSphericalProposal(*p, 1.0);
      positionprior = new SphericalPriorProposal(p, pp);
    } else {
      positionprior = SphericalPriorProposal::load(position_prior_file);
      if (positionprior == nullptr) {
	throw GENERALVORONOIS2EXCEPTION("Failed to open position prior/proposal file\n");
      }
    }

    if (birthdeath_proposal_file == nullptr) {
      birthdeathvalueproposal = new PriorSampleProposal(*(prior->get_prior()));
      birthdeathpositionproposal = new PriorSampleSphericalProposal(*(positionprior->get_prior()));
    } else {
      FILE *fp = fopen(birthdeath_proposal_file, "r");
      if (fp == NULL) {
        throw GENERALVORONOIS2EXCEPTION("Failed to open birth/death proposal file\n");
      }
      
      birthdeathvalueproposal = Proposal::load(fp, *(prior->get_prior()));
      birthdeathpositionproposal = SphericalProposal::load(fp, *(positionprior->get_prior()));
      
      fclose(fp);

      if (birthdeathvalueproposal == nullptr ||
          birthdeathpositionproposal == nullptr) {
        throw GENERALVORONOIS2EXCEPTION("Failed to load birth/death proposal file\n");
      }
    }

    if (!posterior) {

      current_state = this;
      int n = strlen(input);
      if (gvs2_loaddata_(&n, input, addobservation) < 0) {
	throw GENERALVORONOIS2EXCEPTION("Failed to load observations\n");
      }
      current_state = nullptr;
      
      residual_size = data.obs.size();
      predictions.resize(residual_size);
      residuals.resize(residual_size);
      mean_residuals.resize(residual_size);
      last_valid_residuals.resize(residual_size);
      mean_residual_n = 0;
      for (int i = 0; i < residual_size; i ++) {
	predictions[i] = 0.0;
	residuals[i] = 0.0;
	mean_residuals[i] = 0.0;
	last_valid_residuals[i] = 0.0;
      }

      modelweights.resize(_maxcells);
      
    } else {
      
    }
    
    model = new sphericalvoronoimodel<double>(logspace);

  }

  
  void initialize_mpi(MPI_Comm _communicator, double _temperature)
  {
    MPI_Comm_dup(_communicator, &communicator);

    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    temperature = _temperature;

    int observations = (int)data.obs.size();
    int processes = size;

    mpi_offsets = new int[size];
    mpi_counts = new int[size];
      
    for (int i = 0; i < size; i ++) {
      mpi_counts[i] = observations/processes;
      observations -= mpi_counts[i];
      processes --;
    }

    mpi_offsets[0] = 0;
    for (int i = 1; i < size; i ++) {
      mpi_offsets[i] = mpi_offsets[i - 1] + mpi_counts[i - 1];
    }

    if (mpi_offsets[size - 1] + mpi_counts[size - 1] != (int)data.obs.size()) {
      throw GENERALVORONOIS2EXCEPTION("Failed to distribute data points properly");
    }
  }

  void initialize(const char *initial_model, int initial_cells)
  {
    if (initial_model == nullptr) {

      if (communicator == MPI_COMM_NULL) {

	for (int j = 0; j < initial_cells; j ++) {

	
	  double phi, theta;
	  positionprior->sample(random, phi, theta);
	  model->add_cell(globalS2Voronoi::coord_t(phi, theta),
			  prior->sample(random));

	}
	
      } else {

	for (int j = 0; j < initial_cells; j ++) {

	  double mp[3];

	  if (rank == 0) {
	    double phi, theta;
	    positionprior->sample(random, phi, theta);

	    double v = prior->sample(random);

	    mp[0] = phi;
	    mp[1] = theta;
	    mp[2] = v;
	  }

	  MPI_Bcast(mp, 3, MPI_DOUBLE, 0, communicator);
	    
	  model->add_cell(globalS2Voronoi::coord_t(mp[0], mp[1]),
			  mp[2]);

	}
      }
    } else {
    
      if (initial_cells > 1) {
	INFO("Warning: initial model file provided and initial cells specified: initial cells ignored");
      }

      if (!model->load(initial_model)) {
	throw GENERALVORONOIS2EXCEPTION("Failed to load initial model from %s", initial_model);
      }

      INFO("Loaded model with %d cells", model->ncells());
    }
  }

  double optimize_sa(int iterations, double Tmax, double acceptance, int verbose)
  {
    if (communicator == MPI_COMM_NULL) {

      //
      // Auto scaling for efficient jump sizes
      //
      double scale_factor = 1.0e-3;
      
      //
      // Initial likelihood
      //
      double current_N = 0.0;
      double current_P = likelihood(current_N, true);

      //
      // Prune unused model parameters
      //
      std::vector<bool> used(model->cells.size(), false);
      for (auto &o : data.obs) {
	for (auto &mi : o.idx) {
	  used[mi] = true;
	}
      }
      
      int pruned_count = 0;
      for (int j = used.size() - 1; j >= 0; j --) {
	if (!used[j]) {
	  model->cells.erase(model->cells.begin() + j);
	  pruned_count ++;
	}
      }
      
      if (pruned_count > 0) {
	printf("Pre  prune: %12.6f %6d\n", current_P, (int)used.size());
	current_P = likelihood(current_N, true);
	printf("Post prune: %12.6f %6d %d\n", current_P, (int)model->cells.size(), pruned_count);
      }
      
      std::vector<double> scale(model->cells.size(), 1.0);
      std::vector<int> local_proposed(model->cells.size(), 0);
      std::vector<int> local_accepted(model->cells.size(), 0);
      
      int proposed = 0;
      int accepted = 0;
      
      for (int i = 0; i < iterations; i ++) {
	
	double T = temperature_power(i, iterations, Tmax, 2.0);
	
	//
	// Perturb model parameters
	//
	// If the perturbation violates prior bounds, we just keep going hoping
	// that the next perturbation won't. The number of prior violations is
	// counted which can be used as a diagnostic. (Proposal width is too large)
	//
	int j = 0;
	for (auto &c : model->cells) {
	  
	  double oldv = c.v;
	  double newv;
	  double logpriorratio;
	  
	  if (prior->propose(random,
			     T * scale[j],
			     oldv,
			     newv,
			     logpriorratio)) {
	    local_proposed[j] ++;
	    proposed ++;
	    
	    c.v = newv;
	    
	    //
	    // Compute perturbed likelihood
	    //
	    double proposed_N = 0.0;
	    double proposed_P = likelihood(proposed_N, false);
	    
	    //
	    // Accept/reject
	    //
	    double u = log(random.uniform());
	    if (u < (current_P - proposed_P)/T) {
	      //
	      // Accept: set likelihood and keep model as is
	      //
	      current_P = proposed_P;
	      accepted ++;
	      local_accepted[j] ++;
	      
	    } else {
	      //
	      // Reject: restore old model
	      //
	      c.v = oldv;
	    }
	    
	    //
	    // Update adaptive step size
	    //
	    if (local_proposed[j] >= 10 && local_proposed[j] % 10 == 0) {
	      
	      double local_acceptance = (double)local_accepted[j]/(double)local_proposed[j];
	      if (local_acceptance < 0.5) {
		scale[j] *= (1.0 - scale_factor);
	      } else {
		scale[j] *= (1.0 + scale_factor);
	      }
	    }
	    
	  }
	  
	}
	
	
	if (proposed > 0) {
	  acceptance = (double)accepted/(double)proposed;
	}
	
	if (verbose > 0 && (i + 1) % verbose == 0) {

	  INFO("%6d %12.6f %10.6f (%8d/%8d)\n",
	       i + 1,
	       current_P,
	       acceptance,
	       accepted, proposed);
	  
	}
	
      }
      
      return current_P;
    } else {

      //
      // Auto scaling for efficient jump sizes
      //
      double scale_factor = 1.0e-3;
      
      //
      // Initial likelihood
      //
      double current_N = 0.0;
      double current_P = likelihood(current_N, true);      

      //
      // Prune unused model parameters
      //
      std::vector<int> used(model->cells.size(), 0);
      for (auto &o : data.obs) {
	for (auto &mi : o.idx) {
	  used[mi] = 1;
	}
      }

      if (rank == 0) {
	MPI_Reduce(MPI_IN_PLACE, used.data(), used.size(), MPI_INT, MPI_LOR, 0, communicator);
      } else {
	MPI_Reduce(used.data(), NULL, used.size(), MPI_INT, MPI_LOR, 0, communicator);
      }
      MPI_Bcast(used.data(), used.size(), MPI_INT, 0, communicator);

      int pruned_count = 0;
      for (int j = (int)(used.size() - 1); j >= 0; j --) {
	if (used[j] == 0) {
	  model->cells.erase(model->cells.begin() + j);
	  pruned_count ++;
	}
      }

      if (model->cells.size() == 0) {
	throw GENERALVORONOIS2EXCEPTION("Pruned all cells!");
      }

      if (rank == 0 && pruned_count > 0) {
	INFO("Pre  prune: %12.6f %6d\n", current_P, (int)used.size());
	current_P = likelihood(current_N, true);
	INFO("Post prune: %12.6f %6d %d\n", current_P, (int)model->cells.size(), pruned_count);
      }

      std::vector<double> scale(model->cells.size(), 1.0);
      std::vector<int> local_proposed(model->cells.size(), 0);
      std::vector<int> local_accepted(model->cells.size(), 0);
      
      int proposed = 0;
      int accepted = 0;
     
	  
      for (int i = 0; i < iterations; i ++) {

	double T = temperature_power(i, iterations, Tmax, 2.0);

	int j = 0;

	for (auto &c : model->cells) {

	  double oldv = c.v;
	  double newv;
	  double logpriorratio;
	  
	  if (rank == 0) {
	    while (!prior->propose(random,
				   T * scale[j],
				   oldv,
				   newv,
				   logpriorratio)) {
	      // Keep proposing until we get a valid proposal
	    }
					   
	  }

	  MPI_Bcast(&newv, 1, MPI_DOUBLE, 0, communicator);

	  proposed ++;
	  local_proposed[j] ++;
	  c.v = newv;

	  //
	  // Compute perturbed likelihood
	  //
	  double proposed_N = 0.0;
	  double proposed_P = likelihood(proposed_N, false);

	  int accept;

	  if (rank == 0) {
	    double u = log(random.uniform());
	    if (u < (current_P - proposed_P)/T) {
	      accept = 1;
	    } else {
	      accept = 0;
	    }
	  }

	  MPI_Bcast(&accept, 1, MPI_INT, 0, communicator);

	  if (accept == 1) {
	    current_P = proposed_P;
	    accepted ++;
	    local_accepted[j] ++;
	    
	  } else {
	    //
	    // Reject: restore old model
	    //
	    c.v = oldv;
	  }
	  
	  if (rank == 0) {

	    //
	    // Update adaptive step size
	    //
	    if (local_proposed[j] >= 10 && local_proposed[j] % 10 == 0) {
	      
	      double local_acceptance = (double)local_accepted[j]/(double)local_proposed[j];
	      if (local_acceptance < 0.5) {
		scale[j] *= (1.0 - scale_factor);
	      } else {
		scale[j] *= (1.0 + scale_factor);
	      }
	    }

	  }

	  j ++;
	}

	if (proposed > 0) {
	  acceptance = (double)accepted/(double)proposed;
	}
	
	if (rank == 0 && verbose > 0 && (i + 1) % verbose == 0) {
	  
	  INFO("%6d %12.6f %10.6f (%8d/%8d)\n",
	       i + 1,
	       current_P,
	       acceptance,
	       accepted, proposed);
	  
	}
      }	

      return current_P;
    }
  }

  double likelihood(double &norm, bool relocate)
  {
    if (data.obs.size() > 0) {
      if (communicator == MPI_COMM_NULL) {

	for (int i = 0; i < (int)data.obs.size(); i ++) {
	  int npoints = data.obs[i].points.size();

	  //
	  // Look up points
	  //
	  if (relocate) {
	    // Slow method.
	    for (int j = 0; j < npoints; j ++) {
	      data.obs[i].idx[j] = -1;
	      data.obs[i].values[j] = model->value_at_point(data.obs[i].points[j],
							    data.obs[i].idx[j]);
	    }
	  } else {
	    // If cache is ok, look up value by index
	    for (int j = 0; j < npoints; j ++) {
	      data.obs[i].values[j] = model->value_at_index(data.obs[i].idx[j]);

#if 0	      
	      int tidx = -1;
	      double tval;

	      tval = model->value_at_point(data.obs[i].points[j], tidx);
	      if (data.obs[i].idx[j] != tidx) {
		throw GENERALVORONOIS2EXCEPTION("Mismatch in index: %d %d",
						data.obs[i].idx[j], tidx);
	      }

	      if (data.obs[i].values[j] != tval) {
		throw GENERALVORONOIS2EXCEPTION("Mismatch in value: %f %f",
						data.obs[i].values[j], tval);
	      }
#endif
	      
	    }
	  }

	  //
	  // Custom forwardmodel
	  //
	  if (gvs2_compute_prediction_(&i,
				       &npoints,
				       data.obs[i].values.data(),
				       data.obs[i].weights.data(),
				       &predictions[i]) < 0) {
	    throw GENERALVORONOIS2EXCEPTION("Failed to compute predictions");
	  }
	}

	//
	// Likelihood
	//
	int nobs = data.obs.size();
	double hvalue = hierarchical->get(0);
	double like;
	
	if (gvs2_compute_likelihood_(&nobs,
				     &hvalue,
				     predictions.data(),
				     residuals.data(),
				     modelweights.data(),
				     &like,
				     &norm) < 0) {
	  throw GENERALVORONOIS2EXCEPTION("Failed to compute likelihood");
	}

	return like;
      } else {

	for (int i = 0; i < mpi_counts[rank]; i ++) {
	  int o = mpi_offsets[rank] + i;
	  int npoints = data.obs[o].points.size();

	  //
	  // Look up points
	  //
	  if (relocate) {
	    
	    for (int j = 0; j < npoints; j ++) {
	      data.obs[o].values[j] = model->value_at_point(data.obs[o].points[j],
							    data.obs[o].idx[j]);
	    }

	  } else {
	    // If cache is ok, look up value by index
	    for (int j = 0; j < npoints; j ++) {
	      data.obs[o].values[j] = model->value_at_index(data.obs[o].idx[j]);
	    }
	  }

	  //
	  // Custom forwardmodel
	  //
	  
	  if (gvs2_compute_prediction_(&o,
				       &npoints,
				       data.obs[o].values.data(),
				       data.obs[o].weights.data(),
				       &predictions[o]) < 0) {
	    throw GENERALVORONOIS2EXCEPTION("Failed to compute predictions");
	  }
	}
	  
	//
	// Collect predictions to root
	//
	MPI_Gatherv(predictions.data() + mpi_offsets[rank],
		    mpi_counts[rank],
		    MPI_DOUBLE,
		    predictions.data(),
		    mpi_counts,
		    mpi_offsets,
		    MPI_DOUBLE,
		    0,
		    communicator);

	double like;
	
	if (rank == 0) {
	  int nobs = data.obs.size();
	  double hvalue = hierarchical->get(0);
	  if (gvs2_compute_likelihood_(&nobs,
				       &hvalue,
				       predictions.data(),
				       residuals.data(),
				       modelweights.data(),
				       &like,
				       &norm) < 0) {
	    throw GENERALVORONOIS2EXCEPTION("Failed to compute likelihood");
	  }
	}

	MPI_Bcast(&like, 1, MPI_DOUBLE, 0, communicator);
	MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, communicator);

	return like;
	
      }
    } else {
      norm = 0.0;
      return 1.0;
    }
  }

  void accept()
  {
    for (int i = 0; i < residual_size; i ++) {
      last_valid_residuals[i] = residuals[i];
    }
    update_mean_residual();
  }

  void reject()
  {
    update_mean_residual();
  }

  void update_mean_residual()
  {
    mean_residual_n ++;

    for (int i = 0; i < residual_size; i ++) {
      double delta = last_valid_residuals[i] - mean_residuals[i];
      mean_residuals[i] += delta/(double)mean_residual_n;
    }
  }

  static globalS2Voronoi *current_state;
  static int addobservation(int *npoints, double *lons, double *lats);

  MPI_Comm communicator;
  int rank;
  int size;
  int *mpi_counts;
  int *mpi_offsets;

  GeneralVoronoiS2Observations data;
  sphericalvoronoimodel<double> *model;

  PriorProposal *prior;
  SphericalPriorProposal *positionprior;
  PriorProposal *hierarchicalprior;

  Proposal *birthdeathvalueproposal;
  SphericalProposal *birthdeathpositionproposal;

  hierarchical_model *hierarchical;
  double temperature;

  int residual_size;
  int mean_residual_n;

  std::vector<double> mean_residuals;
  std::vector<double> predictions;
  std::vector<double> residuals;
  std::vector<double> last_valid_residuals;

  std::vector<double> modelweights;
  
  int maxcells;
  
  Rng random;
  
};

#endif // globalspherical_hpp
