
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
	$(shell gsl-config --libs)

TARGETOBJS = attenuationtomography.o \
	ak135.o

TARGETS = mksynthetictomography \
	attenuationvoronois2 \
	attenuationvoronois2pt

all : $(TARGETS)

attenuationtomography.o : attenuationtomography.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

../%.o : ../%.cpp
	make -C ..  $*.o

mksynthetictomography : $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../mksynthetic.o $(LIBS)

attenuationvoronois2 : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2.o $(LIBS)

attenuationvoronois2pt : $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../generalvoronois2pt.o $(LIBS)

clean :
	rm -f $(TARGETS) *.o
