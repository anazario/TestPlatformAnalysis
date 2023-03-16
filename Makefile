ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)
#ROOTLIBS      = $(shell root-config --libs)

#ROOTBINDIR    = $(shell root-config --bindir)

CXX           = g++ -g

CXXFLAGS      = $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))
CXXFLAGS      += $(filter-out -stdlib=libc++ -pthread , $(RFCFLAGS))

GLIBS         = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
GLIBS         += $(filter-out -stdlib=libc++ -pthread , $(RFGLIBS))
#GLIBS         += $(filter-out -stdlib=libc++ -pthread , $(ROOTLIBS))
GLIBS         += -lRooFit -lMinuit -lRooFitCore -lhdf5 -lhdf5_cpp -larmadillo -lfftw3

#HDF5          = -lhdf5 -lhdf5_cpp
#ARMA          = -larmadillo

INCLUDEDIR       = ./include/
SRCDIR           = ./src/
CXX	         += -I$(INCLUDEDIR) -I.
OUTOBJ	         = ./obj/

CC_FILES := $(wildcard src/*.cc)
H_FILES := $(wildcard include/*.h)
OBJDIR := $(OUTOBJ)
OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))

all: test.x h5test.x PulseFit.x ProcessScan.x ScanPulseVariation.x

test.x: $(SRCDIR)test.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o test.x $(OUTOBJ)*.o $(GLIBS) $ $<
	touch test.x

h5test.x: $(SRCDIR)h5test.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o h5test.x $(OUTOBJ)*.o $(GLIBS) $ $<
	touch h5test.x

PulseFit.x: $(SRCDIR)PulseFit.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o PulseFit.x $(OUTOBJ)*.o $(GLIBS) $ $<
	touch PulseFit.x

ProcessScan.x: $(SRCDIR)ProcessScan.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o ProcessScan.x $(OUTOBJ)*.o $(GLIBS) $ $<
	touch ProcessScan.x

ScanPulseVariation.x: $(SRCDIR)ScanPulseVariation.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o ScanPulseVariation.x $(OUTOBJ)*.o $(GLIBS) $ $<
	touch ScanPulseVariation.x

$(OUTOBJ)%.o: src/%.cc include/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_FILES): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(OUTOBJ)*.o 
	rm -f *.x
