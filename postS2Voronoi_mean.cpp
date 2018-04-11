
#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "coordinate.hpp"
#include "sphericalvoronoimodel.hpp"
#include "chainhistoryVoronoi.hpp"

typedef sphericalcoordinate<double> coord_t;
typedef chainhistoryreaderVoronoi<coord_t, double> chainhistoryreader_t;

static char short_options[] = "i:fo:m:M:T:V:e:E:g:b:I:z:Z:W:H:t:s:Lh";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"fake", required_argument, 0, 'f'},
  {"output", required_argument, 0, 'o'},
  {"median", required_argument, 0, 'm'},
  {"mode", required_argument, 0, 'M'},
  {"stddev", required_argument, 0, 'T'},
  {"variance", required_argument, 0, 'V'},
  {"credmin", required_argument, 0, 'e'},
  {"credmax", required_argument, 0, 'E'},
  
  {"histogram", required_argument, 0, 'g'},
  {"histogram-bins", required_argument, 0, 'b'},

  {"credinterval", required_argument, 0, 'I'},

  {"zmin", required_argument, 0, 'z'},
  {"zmax", required_argument, 0, 'Z'},

  {"lonsamples", required_argument, 0, 'W'},
  {"latsamples", required_argument, 0, 'H'},
  
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"logspace", no_argument, 0, 'L'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
 
};

static void usage(const char *pname);

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins);
static double median_from_histogram(int *hist, double vmin, double vmax, int bins);
static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);
static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);

static int saveimage(const char *filename, double *image, int width, int height);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  int lonsamples;
  int latsamples;
  
  //
  // Chain processing
  //
  int skip;
  int thin;

  //
  // Input files
  //
  char *input;

  //
  // Output Files
  //
  char *output; // mean
  char *histogram_file;
  char *median_file;
  char *mode_file;
  char *stddev_file;
  char *variance_file;
  char *credmin_file;
  char *credmax_file;

  int *histogram;
  double histmin;
  double histmax;
  int histbins;
  int histsize;
  int histrows;
  int histcols;

  double credinterval;

  int fake;

  bool logspace;

  //
  // State
  //
  int imagesize;
  double *image;

  int meann;
  double delta;
  double *mean;
  double *variance;
  
  //
  // Defaults
  //

  lonsamples = 16;
  latsamples = 16;
  
  input = nullptr;
  fake = 0;
  output = nullptr;
  histogram_file = nullptr;
  median_file = nullptr;
  mode_file = nullptr;
  stddev_file = nullptr;
  variance_file = nullptr;
  credmin_file = nullptr;
  credmax_file = nullptr;

  histmin = 0.0;
  histmax = 1000.0;
  histbins = 500;
  histsize = -1;
  histrows = 0;
  histcols = 0;

  credinterval = 0.90;
  
  skip = 0;
  thin = 0;

  logspace = false;
  
  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input = optarg;
      break;

    case 'f':
      fake = 1;
      break;

    case 'o':
      output = optarg;
      break;

    case 'g':
      histogram_file = optarg;
      break;

    case 'm':
      median_file = optarg;
      break;

    case 'M':
      mode_file = optarg;
      break;

    case 'T':
      stddev_file = optarg;
      break;

    case 'V':
      variance_file = optarg;
      break;

    case 'e':
      credmin_file = optarg;
      break;

    case 'E':
      credmax_file = optarg;
      break;

    case 'b':
      histbins = atoi(optarg);
      if (histbins < 10) {
	fprintf(stderr, "error: bins must be larger than 10\n");
	return -1;
      }
      break;
      
    case 'I':
      credinterval = atof(optarg);
      if (credinterval <= 0.0 ||
	  credinterval >= 1.0) {
	fprintf(stderr, "error: credible interval must be greater than 0 and less than 1\n");
	return -1;
      }
      break;

    case 'W':
      lonsamples = atoi(optarg);
      if (lonsamples < 1) {
	fprintf(stderr, "error: xsamples must be 1 or greater\n");
	return -1;
      }
      break;

    case 'H':
      latsamples = atoi(optarg);
      if (latsamples < 1) {
	fprintf(stderr, "error: latsamples must be 1 or greater\n");
	return -1;
      }
      break;

    case 'z':
      histmin = atof(optarg);
      break;

    case 'Z':
      histmax = atof(optarg);
      break;
      
    case 's':
      skip = atoi(optarg);
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 'L':
      logspace = true;
      break;

    default:
      fprintf(stderr, "error: invalid option '%c'\n", c);
      
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (input == nullptr && !fake) {
    fprintf(stderr, "error: required input file parameter missing\n");
    return -1;
  }

  if (output == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }

  //
  // Initialize state
  //
  histogram = nullptr;
  histrows = lonsamples;
  histcols = latsamples;
  
  histsize = histrows * histcols * histbins;
  histogram = new int[histsize];
  for (int i = 0; i < histsize; i ++) {
    histogram[i] = 0;
  }
			  

  meann = 0;
  imagesize = lonsamples * latsamples;
  mean = new double[imagesize];
  image = new double[imagesize];
  variance = new double[imagesize];

  for (int i = 0; i < imagesize; i ++) {
    image[i] = 0.0;
    mean[i] = 0.0;
    variance[i] = 0.0;
  }

  if (fake) {
    //
    // Fake the sub-sampled image
    //
    for (int j = 0; j < latsamples; j ++) {
      
      // North Pole to South Pole
      double imagephi = ((double)j + 0.5)/(double)latsamples * M_PI;
      double lat = 90.0 - imagephi * 180.0/M_PI;
      
      for (int i = 0; i < lonsamples; i ++) {
	
	// -180 .. 180
	double imagetheta = ((double)i + 0.5)/(double)lonsamples * 2.0 * M_PI - M_PI;
	double lon = imagetheta * 180.0/M_PI;
	      
	mean[j * lonsamples + i] = 500.0 * exp(-((lat + 35.0)*(lat + 35.0) + (lon - 140.0)*(lon - 140.0))/(2.0 * 10.0 * 10.0));
	      
      }
    }

  } else {
    chainhistoryreader_t reader(input);
    
    sphericalvoronoimodel<double> model(logspace);
    singlescaling_hierarchical_model hierarchical;
    double likelihood;
    double norm;
    
    int status = reader.step(model, hierarchical, likelihood, norm);
    int step = 0;
    while (status > 0) {
      
      if (step >= skip) {
	
	if (thin <= 1 || (step - skip) % thin == 0) {
	  //
	  // Compute the sub-sampled image
	  //
	  for (int j = 0; j < latsamples; j ++) {
	    
	    // North Pole to South Pole
	    double imagephi = ((double)j + 0.5)/(double)latsamples * M_PI;
	    for (int i = 0; i < lonsamples; i ++) {
	      
	      // -180 .. 180
	      double imagetheta = ((double)i + 0.5)/(double)lonsamples * 2.0 * M_PI - M_PI;
	      
	      image[j * lonsamples + i] = model.value_at_point(coord_t(imagephi, imagetheta));
	      
	    }
	  }
	  
	  //
	  // Update mean/variance and hist counts
	  //
	  meann ++;
	  for (int i = 0; i < imagesize; i ++) {
	    
	    delta = image[i] - mean[i];
	    mean[i] += delta/(double)meann;
	    variance[i] += delta * (image[i] - mean[i]);
	    
	    int hi = (image[i] - histmin)/(histmax - histmin) * (double)(histbins);
	    if (hi >= 0 && hi < histbins) {
	      histogram[i * histbins + hi] ++;
	    }
	  }
	}
      }

      status = reader.step(model, hierarchical, likelihood, norm);
      step ++;
      
      if ((step % 100000) == 0) {
	printf("%d\n", step);
      }
    }
    
    if (status < 0) {
      fprintf(stderr, "error: failed to step through chain history\n");
      return -1;
    }
  }
  
  //
  // Always save the mean
  //
  if (saveimage(output, mean, lonsamples, latsamples) < 0) {
    fprintf(stderr, "error: failed to save mean\n");
    return -1;
  }

  delete [] mean;
  mean = nullptr;

  //
  // Finalize variance
  //
  if (meann > 1) {
    for (int i = 0; i < imagesize; i ++) {
      variance[i] /= (double)(meann - 1);
    }
  }

  if (variance_file != nullptr) {
    if (saveimage(variance_file, variance, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save variance\n");
      return -1;
    }
  }

  if (stddev_file != nullptr) {
    for (int i = 0; i < imagesize; i ++) {
      variance[i] = sqrt(variance[i]);
    }

    if (saveimage(stddev_file, variance, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save std dev\n");
      return -1;
    }
  }

  delete [] variance;
  variance = nullptr;

  //
  // Save the histogram if required
  //
  if (histogram_file != nullptr) {

    printf("Saving histogram: %s\n", histogram_file);
    
    FILE *fp = fopen(histogram_file, "w");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to create histogram\n");
      return -1;
    }

    fprintf(fp, "%d %d %d\n", histrows, histcols, histbins);
    fprintf(fp, "%10.6f %10.6f\n", histmin, histmax);

    for (int i = 0; i < imagesize; i ++) {
      for (int j = 0; j < histbins; j ++) {
	fprintf(fp, "%d ", histogram[i * histbins + j]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);

  }

  if (median_file != nullptr) {

    for (int i = 0; i < imagesize; i ++) {
      image[i] = median_from_histogram(histogram + i * histbins, histmin, histmax, histbins);
    }

    if (saveimage(median_file, image, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save median\n");
      return -1;
    }

  }

  if (mode_file != nullptr) {

    for (int i = 0; i < imagesize; i ++) {
      image[i] = mode_from_histogram(histogram + i * histbins, histmin, histmax, histbins);
    }

    if (saveimage(mode_file, image, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save median\n");
      return -1;
    }

  }

  if (credmin_file != nullptr) {
    int credible_drop = (int)(((double)meann * (1.0 - credinterval))/2.0);
    for (int i = 0; i < imagesize; i ++) {
      image[i] = head_from_histogram(histogram + i * histbins, histmin, histmax, histbins, credible_drop);
    }

    if (saveimage(credmin_file, image, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save credible min\n");
      return -1;
    }

  }
    
  if (credmax_file != nullptr) {
    int credible_drop = (int)(((double)meann * (1.0 - credinterval))/2.0);
    for (int i = 0; i < imagesize; i ++) {
      image[i] = tail_from_histogram(histogram + i * histbins, histmin, histmax, histbins, credible_drop);
    }

    if (saveimage(credmax_file, image, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save credible max\n");
      return -1;
    }

  }

  delete [] histogram;

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>       Input observations file (required)\n"
	  " -o|--output <filename>      Output likelihood file (required)\n"
	  "\n"
	  " -m|--median <filename>      Median output\n"
	  " -M|--mode <filename>        Modal output\n"
	  " -T|--stddev <filename>      Std dev. output\n"
	  " -V|--variance <filename>    Variance output\n"
	  " -e|--credmin <filename>     Credible min\n"
	  " -E|--credmax <filename>     Credible max\n"
	  "\n"
	  " -g|--histogram <filename>   Histogram output\n"
	  " -b|--histogram-bins <int>   No. bins in histogram\n"
	  "\n"
	  " -I|--credinterval <float>   Credible interval (default = 0.95)\n"
	  "\n"
	  " -z|--zmin <float>           Min value of histogram\n"
	  " -Z|--zmax <float>           Max value of histogram\n"
	  "\n"
	  " -W|--lonsamples <int>       No. samples in longitude direction\n"
	  " -H|--latsamples <int>       No. samples in latitude direction\n"
	  "\n"
	  " -t|--thin <int>             Only use every nth model\n"
	  " -s|--skip <int>             Skip first n models\n"
	  "\n"
	  " -L|--logspace               Chain models are in logspace\n"
	  "\n"
	  " -h|--help                   Usage\n"
	  "\n",
	  pname);
}

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins)
{
  int i;
  int m;
  int mi;

  m = 0;
  mi = -1;

  for (i = 0; i < bins; i ++) {
    if (hist[i] > m) {
      m = hist[i];
      mi = i;
    }
  }
  
  if (mi < 0) {
    return 0.0;
  }

  return ((double)mi + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double median_from_histogram(int *hist, double vmin, double vmax, int bins)
{
  int i;
  int j;
  int ci;
  int cj;

  i = 0;
  j = bins - 1;
  ci = 0;
  cj = 0;

  while (i != j) {
    if (ci < cj) {
      ci += hist[i];
      i ++;
    } else {
      cj += hist[j];
      j --;
    }
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
  int i;
  int ci;

  i = 0; 
  ci = 0;
  while(i < bins && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i ++;
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
  int i;
  int ci;

  i = bins - 1; 
  ci = 0;
  while(i > 0 && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i --;
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static int saveimage(const char *filename, double *image, int width, int height)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return -1;
  }

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {

      fprintf(fp, "%10.6f ", image[width * j + i]);

    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}
