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

#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "coordinate.hpp"
#include "sphericalvoronoimodel.hpp"
#include "chainhistoryVoronoi.hpp"

typedef sphericalcoordinate<double> coord_t;
typedef chainhistoryreaderVoronoi<coord_t, double> chainhistoryreader_t;

static char short_options[] = "i:o:H:K:t:s:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"hierarchical", required_argument, 0, 'H'},
  {"khistory", required_argument, 0, 'K'},
  
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
 
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
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
  char *output; // Likelihood
  char *hierarchical;
  char *khistory;
  
  //
  // Defaults
  //

  input = nullptr;
  output = nullptr;
  hierarchical = nullptr;
  khistory = nullptr;
  
  skip = 0;
  thin = 0;
  
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

    case 'o':
      output = optarg;
      break;

    case 'H':
      hierarchical = optarg;
      break;

    case 'K':
      khistory = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      if (thin < 0) {
	fprintf(stderr, "error: thin must be 0 or positive\n");
	return -1;
      }
      break;

    case 's':
      skip = atoi(optarg);
      if (skip < 0) {
	fprintf(stderr, "error: skip must be 0 or greater\n");
	return -1;
      }
      break;

    default:
      fprintf(stderr, "error: invalid option '%c'\n", c);
      
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (input == nullptr) {
    fprintf(stderr, "error: required input file parameter missing\n");
    return -1;
  }

  if (output == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }

  chainhistoryreader_t reader(input);
    
  sphericalvoronoimodel<double> model(false);
  singlescaling_hierarchical_model hierarchical_model;
  double likelihood;
  double norm;
    
  int status = reader.step(model, hierarchical_model, likelihood, norm);
  int step = 0;

  FILE *fp_like = fopen(output, "w");
  if (fp_like == NULL) {
    fprintf(stderr, "error: failed to create likelihood file\n");
    return -1;
  }

  FILE *fp_hierarchical = NULL;
  if (hierarchical) {
    fp_hierarchical = fopen(hierarchical, "w");
    if (fp_hierarchical == NULL) {
      fprintf(stderr, "error: failed to create hierarchical file\n");
      return -1;
    }
  }

  FILE *fp_khistory = NULL;
  if (khistory) {
    fp_khistory = fopen(khistory, "w");
    if (fp_khistory == NULL) {
      fprintf(stderr, "error: failed to create khistory file\n");
      return -1;

    }
  }
  
  while (status > 0) {
      
    if (step >= skip) {
	
      if (thin <= 1 || (step - skip) % thin == 0) {
	
	fprintf(fp_like, "%15.9f %15.9f\n", likelihood, norm);
	
	if (fp_hierarchical != NULL) {
	  for (int i = 0; i < hierarchical_model.get_nhierarchical(); i ++) {
	    fprintf(fp_hierarchical, "%15.9f ", hierarchical_model.get(i));
	  }
	  fprintf(fp_hierarchical, "\n");
	}

	if (fp_khistory != NULL) {
	  fprintf(fp_khistory, "%d\n", (int)model.cells.size());
	}
      }
    }

    status = reader.step(model, hierarchical_model, likelihood, norm);
    step ++;
    
    if ((step % 100000) == 0) {
      printf("%d\n", step);
    }
  }
  
  if (status < 0) {
    fprintf(stderr, "error: failed to step through chain history\n");
    return -1;
  }

  fclose(fp_like);
  if (fp_hierarchical != NULL) {
    fclose(fp_hierarchical);
  }
  if (fp_khistory != NULL) {
    fclose(fp_khistory);
  }

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>        Input observations file (required)\n"
	  " -o|--output <filename>       Output likelihood file (required)\n"
	  " -H|--hierarchical <filename> Output hierarchical file (optional)\n"
	  "\n"
	  " -t|--thin <int>             Only use every nth model\n"
	  " -s|--skip <int>             Skip first n models\n"
	  "\n"
	  " -h|--help                   Usage\n"
	  "\n",
	  pname);
}
