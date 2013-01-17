
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

MODULES := $(addprefix nrlib/, $(NRLIB_MODULES)) flens

NRLIB_BASE_DIR := nrlib/

SRC :=
# Files with sources for each of the modules.
include $(patsubst %, %/module.mk, $(MODULES))
# This file contains the files we want from NRLib where we do not want entire directory.
include submodulescrava.mk

OBJDIR = ../obj/libs/

OBJ_NRLIB := $(patsubst %.cpp, $(OBJDIR)%.o, $(filter %.cpp, $(SRC)))
OBJ_FLENS := $(patsubst %.cc, $(OBJDIR)%.o, $(filter %.cc,  $(SRC)))

DIRS := $(sort $(dir $(OBJ_NRLIB)))
DIRS += $(sort $(dir $(OBJ_FLENS)))

INCLUDE := -I../libs

CPPFLAGS += $(INCLUDE)

$(OBJDIR)%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(OBJDIR)%.o : %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLENS) $< -o $@

all: $(DIRS) $(OBJ_NRLIB) $(OBJ_FLENS)

help:
	@echo 'MODULES:'
	@echo $(MODULES)
	@echo 'A:'
	@echo $(A)

$(DIRS):
	install -d $@
