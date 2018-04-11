#pragma once
#ifndef genericinterface_hpp
#define genericinterface_hpp

extern "C" {

  
  typedef int (gvs2_addobservation_t)(int *npoints,
				      double *lons,
				      double *lats);
  //
  // Load data and for each observation, call a callback to tell
  // the inversion about the points involved
  //
  int gvs2_loaddata_(int *n,
		     const char *filename,
		     gvs2_addobservation_t addobs);

  //
  // For a single observation, compute the prediction given
  // model values for each point
  //
  int gvs2_compute_prediction_(int *observation,
			       int *npoints,
			       const double *value,
			       double *unused,
			       double *prediction);
  //
  // For the observations, given predictions, compute residuals, likelihood
  // and norm
  //
  int gvs2_compute_likelihood_(int *nobservation,
			       double *hierarchical,
			       double *predictions,
			       double *residuals,
			       double *unused,
			       double *like,
			       double *norm);

  //
  // Used for making synthetic datasets, save data in correct format
  // using predictions to overwrite observations
  //
  int gvs2_savedata_(int *n,
		     const char *filename,
		     double *hierarchical,
		     int *nobservations,
		     double *predictions);

  //
  // Finish (save state for restarting)
  //
  int gvs2_finish_(const char *path);

  //
  // Starting (pre load)
  //
  int gvs2_start_(const char *path);

};

#endif // genericinterface_hpp

		 
