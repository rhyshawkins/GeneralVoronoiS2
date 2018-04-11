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

extern "C" {
  #include "slog.h"
};

class globalS2Voronoi {
public:

  typedef sphericalcoordinate<double> coord_t;

  globalS2Voronoi(const char *input,
		  const char *initial_model,
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
    }

    if (hierarchical_prior_file == nullptr) {
      Prior *p = new UniformPrior(0.1, 5.0);
      Proposal *pp = new GaussianProposal(*p, 0.05);
      hierarchicalprior = new PriorProposal(p, pp);
    } else {
      hierarchicalprior = PriorProposal::load(hierarchical_prior_file);
    }

    if (position_prior_file == nullptr) {
      SphericalPrior *p = new UniformSphericalPrior();
      SphericalProposal *pp = new VonMisesSphericalProposal(*p, 1.0);
      positionprior = new SphericalPriorProposal(p, pp);
    } else {
      positionprior = SphericalPriorProposal::load(position_prior_file);
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

    if (initial_model == nullptr) {

      //
      // Initialize to single cell sampled from prior
      //
      model->add_cell(coord_t(0.0, 0.0), prior->sample(random));

    } else {

      if (!model->load(initial_model)) {
	throw GENERALVORONOIS2EXCEPTION("Failed to load initial model from %s", initial_model);
      }

      INFO("Loaded model with %d cells", model->ncells());

    }
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

  double likelihood(double &norm)
  {
    if (data.obs.size() > 0) {
      if (communicator == MPI_COMM_NULL) {

	for (int i = 0; i < (int)data.obs.size(); i ++) {
	  int npoints = data.obs[i].points.size();

	  //
	  // Look up points
	  //
	  for (int j = 0; j < npoints; j ++) {
	    data.obs[i].values[j] = model->value_at_point(data.obs[i].points[j]);
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
	  for (int j = 0; j < npoints; j ++) {
	    data.obs[o].values[j] = model->value_at_point(data.obs[o].points[j]);
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
	// Share predictions
	//
	MPI_Allgatherv(predictions.data() + mpi_offsets[rank],
		       mpi_counts[rank],
		       MPI_DOUBLE,
		       predictions.data(),
		       mpi_counts,
		       mpi_offsets,
		       MPI_DOUBLE,
		       communicator);

	
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
