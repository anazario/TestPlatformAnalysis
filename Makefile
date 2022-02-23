ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++

CXXFLAGS      = $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))
CXXFLAGS      += $(filter-out -stdlib=libc++ -pthread , $(RFCFLAGS))

GLIBS         = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
GLIBS         += $(filter-out -stdlib=libc++ -pthread , $(RFGLIBS))
GLIBS         += -lRooFit -lMinuit -lRooFitCore

HDF5          = -lhdf5 -lhdf5_cpp

INCLUDEDIR       = ./include/
SRCDIR           = ./src/
CXX	         += -I$(INCLUDEDIR) -I.
OUTOBJ	         = ./obj/

CC_FILES := $(wildcard src/*.cc)
H_FILES := $(wildcard include/*.h)
OBJDIR := $(OUTOBJ)
OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))

all: h5test.x ProcessScan.x ScanPulseVariation.x

h5test.x: $(SRCDIR)h5test.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o h5test.x $(OUTOBJ)*.o $(GLIBS) $(HDF5) $ $<
	touch h5test.x

ProcessScan.x: $(SRCDIR)ProcessScan.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o ProcessScan.x $(OUTOBJ)*.o $(GLIBS) $(HDF5) $ $<
	touch ProcessScan.x

ScanPulseVariation.x: $(SRCDIR)ScanPulseVariation.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) -o ScanPulseVariation.x $(OUTOBJ)*.o $(GLIBS) $(HDF5) $ $<
	touch ScanPulseVariation.x

$(OUTOBJ)%.o: src/%.cc include/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_FILES): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(OUTOBJ)*.o 
	rm -f *.x
