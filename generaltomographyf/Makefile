
FCC=mpif90
FCFLAGS=-c -g -Wall -std=f2008

CXX=mpicxx
CXXFLAGS=-c -g -Wall -I..

CXXFLAGS+=-O3

BASEOBJS = \
	../generalvoronois2exception.o \
	../hierarchical_model.o \
	../pathutil.o \
	../prior.o \
	../rng.o \
	../sphericalprior.o \
	../globalS2Voronoi.o \
	../simulated_annealing_scales.o

TDTBASE=../../TDTbase
LIBS = -L$(TDTBASE)/log -llog \
	$(shell gsl-config --libs) \
	-lgfortran

TARGETOBJS = generictomography.o

TARGETS = mksynthetictomography \
	tomographyvoronois2 \
	tomographyvoronois2pt

all : $(TARGETS)

generictomography.o : generictomography.f90
	$(FCC) $(FCFLAGS) -o $@ $<

../%.o : ../%.cpp
	make -C ..  $*.o

mksynthetictomography : $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o $(LIBS)

tomographyvoronois2 : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o $(LIBS)

tomographyvoronois2pt : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o $(LIBS)

clean :
	rm -f $(TARGETS)
