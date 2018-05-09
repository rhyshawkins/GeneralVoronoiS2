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

#include <math.h>

#include <map>
#include <string>

#include <getopt.h>
#include <string.h>

#include "generalvoronois2observations.hpp"
#include "coordinate.hpp"
#include "rng.hpp"

#include "genericinterface.hpp"

double synthetic_constant(double phi, double theta)
{
  return 0.5;
}

double synthetic_eastwest(double phi, double theta)
{
  if (theta < 0.0) {
    return 0.25;
  } else {
    return 0.75;
  }
}

double synthetic_northsouth(double phi, double theta)
{
  if (phi < M_PI/2.0) {
    return 0.33;
  } else {
    return 0.66;
  }
}

double synthetic_cubedsphere(double phi, double theta)
{
  vector3<double> v;

  sphericalcoordinate<double>::sphericaltocartesian(phi, theta, v);

  if (fabs(v.x) > fabs(v.y)) {
    if (fabs(v.x) > fabs(v.z)) {
      if (v.x < 0.0) {
	return 0.25;
      } else {
	return 0.625;
      }
    } else {
      if (v.z < 0.0) {
	return 0.5;
      } else {
	return 0.375;
      }
    }
  } else if (fabs(v.y) > fabs(v.z)) {
    if (v.y < 0.0) {
      return 0.125;
    } else {
      return 0.750;
    }
  } else {

    if (v.z < 0.0) {
      return 0.5;
    } else {
      return 0.375;
    }
  }
}

//
// Extra synthetic tests added by Sima Mousavi, and modified by Rhys Hawkins
//
double synthetic_checkerboard(double phi, double theta)
{
  static constexpr double RADIUS = M_PI/8.0;
  static constexpr double CLAT = M_PI/3.0;
  static constexpr double CLON = M_PI/4.0;
  
  static constexpr double HIGH = 0.75;
  static constexpr double LOW = 0.25;
  sphericalcoordinate<double> p(phi, theta);
  
  if (theta < 0.0) {

    //
    // Low background hemisphere
    //

    if (theta > -M_PI/2.0) {

      if (phi < M_PI/2.0) {
	if (p.distance(sphericalcoordinate<double>(CLAT, -CLON)) < RADIUS) {
	  return HIGH;
	}
      } else {
	if (p.distance(sphericalcoordinate<double>(M_PI - CLAT, -CLON)) < RADIUS) {
	  return HIGH;
	}
      }
	
    } else {

      if (phi < M_PI/2.0) {
	if (p.distance(sphericalcoordinate<double>(CLAT, -M_PI + CLON)) < RADIUS) {
	  return HIGH;
	}
      } else {
	if (p.distance(sphericalcoordinate<double>(M_PI - CLAT, -M_PI + CLON)) < RADIUS) {
	  return HIGH;
	}
      }

    }

    return LOW;

  } else {

    //
    // Low background hemisphere
    //

    if (theta < M_PI/2.0) {

      if (phi < M_PI/2.0) {
	if (p.distance(sphericalcoordinate<double>(CLAT, CLON)) < RADIUS) {
	  return LOW;
	}
      } else {
	if (p.distance(sphericalcoordinate<double>(M_PI - CLAT, CLON)) < RADIUS) {
	  return LOW;
	}
      }
	
    } else {

      if (phi < M_PI/2.0) {
	if (p.distance(sphericalcoordinate<double>(CLAT, M_PI - CLON)) < RADIUS) {
	  return LOW;
	}
      } else {
	if (p.distance(sphericalcoordinate<double>(M_PI - CLAT, M_PI - CLON)) < RADIUS) {
	  return LOW;
	}
      }

    }

    return HIGH;

  }
}

std::map<std::string, synthetic_model_f> models = { {"Constant", synthetic_constant},
						    {"EastWest", synthetic_eastwest},
						    {"NorthSouth", synthetic_northsouth},
						    {"CubedSphere", synthetic_cubedsphere},
						    {"CheckerBoard", synthetic_checkerboard}};

static char short_options[] = "i:o:O:I:m:ln:S:W:H:s:f:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"output-true", required_argument, 0, 'O'},

  {"model", required_argument, 0, 'm'},
  {"list-models", no_argument, 0, 'l'},
  
  {"noise", required_argument, 0, 'n'},
  {"seed", required_argument, 0, 'S'},

  {"image-output", required_argument, 0, 'I'},
  {"image-width", required_argument, 0, 'W'},
  {"image-height", required_argument, 0, 'H'},

  {"scale", required_argument, 0, 's'},
  {"offset", required_argument, 0, 'f'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static GeneralVoronoiS2Observations observations;

static int addobservation(int *npoints,
			  double *lons,
			  double *lats)
{
  observations.add(*npoints, lons, lats);

  return 0;
}


static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  char *source_data;
  char *output_file;
  char *output_true;
  int seed;
  double noise_sigma;
  std::string model_name;

  char *image_output;
  int image_width;
  int image_height;

  double scale;
  double offset;
    
  //
  // Defaults
  //
  model_name = "Constant";
  source_data = nullptr;
  output_file = nullptr;
  output_true = nullptr;
  noise_sigma = 0.1;
  seed = 983;

  image_output = nullptr;
  image_width = 128;
  image_height = 64;

  scale = 1.0;
  offset = 0.0;

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      source_data = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'O':
      output_true = optarg;
      break;
      
    case 'I':
      image_output = optarg;
      break;

    case 'm':
      model_name = optarg;
      break;

    case 'l':
      {
	fprintf(stderr, "Models:\n");
	for (auto &mn : models) {
	  fprintf(stderr, "  `%s'\n", mn.first.c_str());
	}
	return -1;
      }
      break;

    case 'n':
      noise_sigma = atof(optarg);
      if (noise_sigma <= 0.0) {
	fprintf(stderr, "error: noise must be greater than 0\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;

    case 'W':
      image_width = atoi(optarg);
      if (image_width < 16) {
	fprintf(stderr, "error: image width must be 16 or greater\n");
	return -1;
      }
      break;

    case 'H':
      image_height = atoi(optarg);
      if (image_height < 16) {
	fprintf(stderr, "error: image_height must be 16 or greater\n");
	return -1;
      }
      break;

    case 's':
      scale = atof(optarg);
      if (scale <= 0.0) {
	fprintf(stderr, "error: scale must be greater than 0\n");
	return -1;
      }
      break;

    case 'f':
      offset = atof(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;

    }
  }
  
  if (source_data == nullptr) {
    fprintf(stderr, "error: required parameter input not set\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required parameter output not set\n");
    return -1;
  }

  
  std::map<std::string, synthetic_model_f>::const_iterator i = models.find(model_name);

  if (i == models.end()) {
    fprintf(stderr, "error: invalid model name: %s\n", model_name.c_str());
    return -1;
  }
  synthetic_model_f model = i->second;

  //
  // Load data
  //
  int n = strlen(source_data);
  if (gvs2_loaddata_(&n, source_data, addobservation) < 0) {
    fprintf(stderr, "error: failed to load data from %s\n", source_data);
    return -1;
  }

  printf("%d observations\n", (int)observations.obs.size());


  //
  // Compute predictions
  //
  int oi = 0;
  std::vector<double> predictions;
  for (auto &o : observations.obs) {

    //
    // Look up synthetic model values
    //
    int pi = 0;
    int npoints = o.points.size();
    for (auto &p : o.points) {

      o.values[pi] = model(p.phi, p.theta) * scale + offset;
      pi ++;
    }
	   
    //
    // Compute prediction
    //
    if (gvs2_compute_prediction_(&oi,
				 &npoints,
				 o.values.data(),
				 o.weights.data(),
				 &o.pred)) {
      fprintf(stderr, "error: received error from compute prediction\n");
      return -1;
    }

    predictions.push_back(o.pred);
    oi ++;
  }				

  //
  // Save true observations
  //
  if (output_true) {

    double noise = 0.0;
    int nobs = observations.obs.size();
    int n = strlen(output_true);
    
    if (gvs2_savedata_(&n,
		       output_true,
		       &noise,
		       &nobs,
		       predictions.data()) < 0) {
      fprintf(stderr, "error: failed to save synthetic true observations\n");
      return -1;
    }
  }

  //
  // Add noise
  //
  Rng random(seed);
  
  for (auto &p : predictions) {
    p += random.normal(noise_sigma);
  }

  //
  // Save noisy observations
  //
  int nobs = observations.obs.size();

  n = strlen(output_file);
  if (gvs2_savedata_(&n,
		     output_file,
		     &noise_sigma,
		     &nobs,
		     predictions.data()) < 0) {
    fprintf(stderr, "error: failed to save synthetic observations\n");
    return -1;
  }

  if (image_output != nullptr) {

    FILE *fp_image = fopen(image_output, "w");
    if (fp_image == NULL) {
      fprintf(stderr, "error: failed to create image output file\n");
      return -1;
    }

    for (int j = 0; j < image_height; j ++) {
      
      // North Pole to South Pole
      double imagephi = ((double)j + 0.5)/(double)image_height * M_PI;
      
      for (int i = 0; i < image_width; i ++) {
	
	// -180 .. 180
	double imagetheta = ((double)i + 0.5)/(double)image_width * 2.0 * M_PI - M_PI;

	fprintf(fp_image, "%15.9f ", model(imagephi, imagetheta) * scale + offset);

      }

      fprintf(fp_image, "\n");
    }

    fclose(fp_image);
  }
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of\n"
	  "\n"
	  " -i | --input <filename>           Input data to base recompute with synthetic model\n"
	  " -o | --output <filename>          Output file to write synthetic noisy observations to\n"
	  " -O | --output-true <filename>     Output file to write synthetic true observations to\n"
	  "\n"
	  " -m | --model <name>               Synthetic model to use\n"
	  " -l | --list-models                List available synthetic models and exit\n"
	  "\n"
	  " -n | --noise <float>              Std dev of gaussian noise to add to observations\n"
	  " -S | --seed <int>                 Random seed\n"
	  "\n"
	  " -I | --image-output <filename>    Write synthetic model image\n"
	  " -W | --image-width <int>          Image width\n"
	  " -H | --image-height <int>         Image height\n"
	  "\n"
	  " -h | --help                       Usage\n"
	  "\n",
	  pname);
}
