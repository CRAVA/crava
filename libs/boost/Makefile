
include ../../Makeheader

# Boost directories where we want all files are given here.
FULLMODULES := system filesystem

SRC :=

# Files with sources for each of the modules.
include $(patsubst %, %/module.mk, $(FULLMODULES))

# This file contains the files we want from Boost where we do not want entire directory.
#include submodulescrava.mk

OBJDIR = ../../obj/libs/boost/

OBJECTS := $(patsubst %.cpp, $(OBJDIR)%.o, $(filter %.cpp, $(SRC)))

DIRS := $(sort $(dir $(OBJECTS)))

INCLUDE := -I../
CPPFLAGS += $(INCLUDE)

$(OBJDIR)%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(DIRS) $(OBJECTS)

$(DIRS):
	install -d $@
