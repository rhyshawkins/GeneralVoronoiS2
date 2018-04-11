
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
	../globalS2Voronoi.o

TDTBASE=../../TDTbase
LIBS = -L$(TDTBASE)/log -llog \
	$(shell gsl-config --libs) \
	-lgfortran

TARGETOBJS = genericregression.o

TARGETS = mksyntheticregression \
	regressionvoronois2 \
	regressionvoronois2pt

all : $(TARGETS)

genericregression.o : genericregression.f90
	$(FCC) $(FCFLAGS) -o $@ $<

../%.o : ../%.cpp
	make -C ..  $*.o

mksyntheticregression : $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o $(LIBS)

regressionvoronois2 : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o $(LIBS)

regressionvoronois2pt : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o $(LIBS)

clean :
	rm -f $(TARGETS)
