
include ../Makeheader

# Directories where we want all files are given here. 
FULLMODULES := flens

SRC := 

# Files with sources for each of the modules.
include $(patsubst %, %/module.mk, $(FULLMODULES))

# This file contains the files we want from NRLib where we do not want entire directory.
#include submodulescrava.mk

OBJDIR = ../obj/nrlib/

OBJECTS := $(patsubst %.cc, $(OBJDIR)%.o, $(filter %.cc, $(SRC)))

DIRS := $(sort $(dir $(OBJECTS)))

INCLUDE := -I./
CPPFLAGS += $(INCLUDE)

$(OBJDIR)%.o : %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(DIRS) $(OBJECTS)

$(DIRS):
	install -d $@
