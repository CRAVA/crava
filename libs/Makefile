
include ../Makeheader

# NRLib directories where we want all files are given here.
NRLIB_MODULES := geometry \
                 grid \
                 iotools \
                 math \
                 random \
                 segy \
                 stormgrid \
                 tinyxml \
                 trend \
                 volume

#MODULES := boost flens $(addprefix nrlib/, $(NRLIB_MODULES))

MODULES := $(addprefix nrlib/, $(NRLIB_MODULES))

NRLIB_BASE_DIR := nrlib/

SRC :=

# Files with sources for each of the modules.
include $(patsubst %, %/module.mk, $(MODULES))

# This file contains the files we want from NRLib where we do not want entire directory.
include submodulescrava.mk

OBJDIR = ../obj/libs/

OBJECTS := $(patsubst %.cpp, $(OBJDIR)%.o, $(filter %.cpp, $(SRC)))

DIRS := $(sort $(dir $(OBJECTS)))

#INCLUDE := -I../libs/
CPPFLAGS += $(INCLUDE)

$(OBJDIR)%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(DIRS) $(OBJECTS)

$(DIRS):
	install -d $@
