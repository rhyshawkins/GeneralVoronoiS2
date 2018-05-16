
TDTBASE=../TDTbase

INCLUDES = \
	$(shell mpicc -showme:compile) \
	-I$(TDTBASE)/log \
	-I$(TDTBASE)/tracking

EXTRA_LIBS = \
	-L$(TDTBASE)/tracking -ltracking \
	-L$(TDTBASE)/log -llog


CXX = g++
CXXFLAGS = -c -g -Wall --std=c++11 $(INCLUDES)

CXXFLAGS += -O3

INSTALL = install
INSTALLFLAGS = -D

LIBS = $(EXTRA_LIBS) \
	-lm \
	$(shell gsl-config --libs) \
	$(shell mpicxx -showme:link)

OBJS = rng.o \
	prior.o \
	sphericalprior.o \
	hierarchical_model.o \
	pathutil.o \
	generalvoronois2exception.o \
	globalS2Voronoi.o \
	simulated_annealing_scales.o

MPI_OBJS = 


SRCS = Makefile \
	birthgenericS2Voronoi.hpp \
	chainhistoryVoronoi.hpp \
	coordinate.hpp \
	deathgenericS2Voronoi.hpp \
	generalvoronois2exception.hpp \
	generalvoronois2observations.hpp \
	generalvoronois2util.hpp \
	genericinterface.hpp \
	globalS2Voronoi.hpp \
	hierarchicalS2Voronoi.hpp \
	hierarchical_model.hpp \
	moveS2Voronoi.hpp \
	pathutil.hpp \
	perturbationS2Voronoi.hpp \
	perturbationcollectionS2Voronoi.hpp \
	prior.hpp \
	rng.hpp \
	sphericalprior.hpp \
	sphericalvoronoimodel.hpp \
	util.hpp \
	valueS2Voronoi.hpp \
	generalvoronois2.cpp \
	generalvoronois2exception.cpp \
	generalvoronois2pt.cpp \
	globalS2Voronoi.cpp \
	hierarchical_model.cpp \
	mksynthetic.cpp \
	pathutil.cpp \
	postS2Voronoi_likelihood.cpp \
	postS2Voronoi_mean.cpp \
	postS2Voronoi_mean_mpi.cpp \
	postS2Voronoi_text.cpp \
	prior.cpp \
	rng.cpp \
	sphericalprior.cpp \
	simulated_annealing.hpp \
	simulated_annealing_scales.cpp \
	simulated_annealing_scales.hpp

EXAMPLES = generalregressioncpp/Makefile \
	generalregressioncpp/genericregression.cpp \
	generalregressioncpp/example/Makefile \
	generaltomographycpp/Makefile \
	generaltomographycpp/generictomography.cpp \
	generaltomographycpp/tomographyutil.hpp \
	generaltomographycpp/example/Makefile \
	generalregressionf/Makefile \
	generalregressionf/genericregression.f90 \
	generalregressionf/example/Makefile \
	generaltomographyf/Makefile \
	generaltomographyf/generictomography.f90 \
	generaltomographyf/example/Makefile

EXTRADIST = \
	documentation/manual.tex \
	documentation/manual.pdf \
	scripts/mksyntheticpathtemplate.py \
	scripts/mksyntheticpointstemplate.py \
	scripts/plot_likelihood_converge.py \
	scripts/plot_image_ortho.py \
	scripts/plot_map.py \
	terrawulf/pbs_singlechain_constant_init.sh \
	terrawulf/pbs_singlechain_constant_cont.sh

TARGETS = postS2Voronoi_mean \
	postS2Voronoi_mean_mpi \
	postS2Voronoi_likelihood \
	postS2Voronoi_text

documentation/manual.pdf : documentation/manual.tex
	cd documentation && pdflatex manual.tex

all : $(TARGETS)

postS2Voronoi_mean : postS2Voronoi_mean.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_mean.o $(OBJS) $(LIBS)

postS2Voronoi_mean_mpi : postS2Voronoi_mean_mpi.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_mean_mpi.o $(OBJS) $(LIBS)

postS2Voronoi_likelihood : postS2Voronoi_likelihood.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_likelihood.o $(OBJS) $(LIBS)

postS2Voronoi_text : postS2Voronoi_text.o $(OBJS)
	$(CXX) -o $@ postS2Voronoi_text.o $(OBJS) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

DATE = $(shell date +"%Y%m%d%H%M")
DIR = GeneralVoronoiS2
TGZ = $(DIR).tar.gz

dist : all documentation/manual.pdf
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS) $(EXAMPLES) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f $(TARGETS) *.o




