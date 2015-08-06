PROGRAMS	 = analyze_events
#PROGRAMS        = jpt_structure
#PROGRAMS        = diffcross_v2
#PROGRAMS        = kinematics

SHELL		= /bin/sh
CPP	 	= /lib/cpp -P

CXX		= g++
LD		= g++
LDFLAGS		= -g
SOFLAGS		= -shared
#PROF_FLAGS	= -pg2
ROOTCFLAGS	= $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS	= $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTGLIBS	= $(shell $(ROOTSYS)/bin/root-config --glibs)
USERINCLUDE	= /Users/nsonmez/Delphes-3.1.2/external/fastjet/plugins/ATLASCone
WARNING		= -ansi -Wall \
#		  -Wconversion \
# 		  -Wstrict-prototypes -Wmissing-prototypes \
# 		  -Wmissing-declarations

NGLIBB         = $(ROOTGLIBS) 
NGLIBB        += -lMinuit
GLIBB          = $(filter-out -lNew, $(NGLIBB))

CXXFLAGS	= -g -O3 -fPIC $(WARNING) $(ROOTCFLAGS) -I$(USERINCLUDE)
LDLIBS		= $(PROF_FLAGS) -g $(ROOTLIBS) -lExRootAnalysis
# -lGui -lRooFitCore -lRooFit -lFoam -lMinuit

all: $(PROGRAMS)

clean:
	rm -f *.o 
	rm -f *Dict* 
	rm -f *LinkDef.h 
	rm -f $(PROGRAMS)
