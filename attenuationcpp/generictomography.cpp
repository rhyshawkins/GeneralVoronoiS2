#include <vector>
#include <stdio.h>
#include <math.h>

#include "genericinterface.hpp"

#include "tomographyutil.hpp"

struct observation {
  double ttime;
  double sigma;

  std::vector<double> lons;
  std::vector<double> lats;
  std::vector<double> rs;
  std::vector<double> weights;
};

static std::vector<struct observation> obs;


extern "C" {
  
  int gvs2_loaddata_(int *n,
		     const char *filename,
		     gvs2_addobservation_t addobs)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      return -1;
    }

    int nobs;
    if (fscanf(fp, "%d\n", &nobs) != 1) {
      return -1;
    }

    obs.resize(nobs);
    for (int i = 0; i < nobs; i ++) {

      int npoints;
      if (fscanf(fp, "%lf %lf %d\n",
		 &obs[i].ttime,
		 &obs[i].sigma,
		 &npoints) != 3) {
	return -1;
      }

      obs[i].lons.resize(npoints);
      obs[i].lats.resize(npoints);
      obs[i].rs.resize(npoints);
      for (int j = 0; j < npoints; j ++) {

	if (fscanf(fp, "%lf %lf %lf\n",
		   &obs[i].lons[j],
		   &obs[i].lats[j],
		   &obs[i].rs[j]) != 3) {
	  return -1;
	}
      }

      if (addobs(&npoints, obs[i].lons.data(), obs[i].lats.data()) < 0) {
	return -1;
      }

      //
      // Precompute distance between points for computing predictions
      // from models.
      //
      compute_weights(obs[i].lons, obs[i].lats, obs[i].rs, obs[i].weights);
    }

    fclose(fp);
    return 0;
  }

  int gvs2_compute_prediction_(int *observation,
			       int *npoints,
			       const double *values,
			       double *unused,
			       double *prediction)
  {
    if (*observation < 0 || *observation >= (int)obs.size()) {
      fprintf(stderr, "gvs2_compute_prediction_: observation out of range\n");
      return -1;
    }
    
    if (*npoints != (int)obs[*observation].weights.size()) {
      fprintf(stderr, "gvs2_compute_prediction_: npoints mismatch\n");
      return -1;
    }

    //
    // Here the values are log(Q) and the weights are the distances/Vp so
    // t* is a sum over weights/exp(values).
    //
    double ttsum = 0.0;
    
    for (int j = 0; j < (*npoints); j ++) {

      ttsum += obs[*observation].weights[j]/exp(values[j]);

    }

    // printf("pred: %2d %10.6f %10.6f\n", *observation, ttsum, obs[*observation].ttime);
    *prediction = ttsum;
    return 0;
  }

  int gvs2_compute_likelihood_(int *nobservation,
                               double *hierarchical,
                               double *predictions,
                               double *residuals,
                               double *unused,
                               double *_like,
                               double *_norm)

  {
    double sum = 0.0;
    double norm = -0.5*log(2.0*M_PI);
    for (int i = 0; i < (*nobservation); i ++) {
      double res = predictions[i] - obs[i].ttime;
      residuals[i] = res;
      double n = hierarchical[0] * obs[i].sigma;

      sum += (res*res)/(2.0*n*n);
      norm += log(n);
    }

    *_like = sum;
    *_norm = norm;

    return 0;
  }

  int gvs2_savedata_(int *n,
		     const char *filename,
                     double *noiselevel,
                     int *nobservations,
                     double *predictions)
  {
    FILE *fp;

    fp = fopen(filename, "w");
    if (fp == NULL) {
      return -1;
    }

    fprintf(fp, "%d\n", (int)obs.size());
    int i = 0;
    for (auto &o : obs) {
      fprintf(fp, "%16.9f %16.9f %d\n",
              predictions[i], noiselevel[0], (int)o.lons.size());

      for (int j = 0; j < (int)o.lons.size(); j ++) {

	fprintf(fp, "%16.9f %16.9f %16.9f\n", o.lons[j], o.lats[j], o.rs[j]);

      }
      i ++;
    }

    fclose(fp);
    return 0;
  }

  
}

