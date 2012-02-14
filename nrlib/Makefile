
include ../Makeheader

# NRLib directories where we want all files are given here. 
FULLMODULES := grid iotools segy stormgrid  volume tinyxml boost random math trend

SRC := 

# Files with sources for each of the modules.
include $(patsubst %, %/module.mk, $(FULLMODULES))

# This file contains the files we want from NRLib where we do not want entire directory.
include submodulescrava.mk

OBJDIR = ../obj/nrlib/

OBJECTS := $(patsubst %.cpp, $(OBJDIR)%.o, $(filter %.cpp, $(SRC)))

DIRS := $(sort $(dir $(OBJECTS)))

INCLUDE := -I./
CPPFLAGS += $(INCLUDE)

$(OBJDIR)%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(DIRS) $(OBJECTS)

$(DIRS):
	install -d $@
