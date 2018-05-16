
CXX=mpicxx
CXXFLAGS=-c -g -Wall --std=c++11 -I..

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
	$(shell gsl-config --libs)

TARGETOBJS = genericregression.o

TARGETS = mksyntheticregression \
	regressionvoronois2 \
	regressionvoronois2pt

all : $(TARGETS)

genericregression.o : genericregression.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

../%.o : ../%.cpp
	make -C ..  $*.o

mksyntheticregression : $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o $(LIBS)

regressionvoronois2 : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o $(LIBS)

regressionvoronois2pt : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o $(LIBS)

clean :
	rm -f $(TARGETS) *.o
