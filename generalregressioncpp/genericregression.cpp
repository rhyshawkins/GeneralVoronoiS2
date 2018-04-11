#include <vector>
#include <stdio.h>
#include <math.h>

#include "genericinterface.hpp"

struct observation {
  double lon;
  double lat;
  double value;
  double sigma;
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
      if (fscanf(fp, "%lf %lf %lf %lf\n",
		 &obs[i].lon,
		 &obs[i].lat,
		 &obs[i].value,
		 &obs[i].sigma) != 4) {
	return -1;
      }

      int n;
      n = 1;
      if (addobs(&n, &obs[i].lon, &obs[i].lat) < 0) {
	return -1;
      }
    }

    fclose(fp);
    return 0;
  }

  int gvs2_compute_prediction_(int *observation,
			       int *npoints,
			       const double *value,
			       double *unused,
			       double *prediction)
  {
    prediction[0] = value[0];
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
      double res = predictions[i] - obs[i].value;
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
      fprintf(fp, "%16.9f %16.9f %16.9f %16.9f\n",
	      o.lon, o.lat, predictions[i], noiselevel[0]);
      i ++;
    }

    fclose(fp);
    return 0;
  }

  
}

