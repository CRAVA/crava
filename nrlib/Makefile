
include ../Makeheader

MODULES := grid iotools segy stormgrid surface volume 

SRC := 

# Files with sources for each of the modules.
include $(patsubst %, %/module.mk, $(MODULES))

OBJDIR = ../obj/nrlib/

OBJECTS := $(patsubst %.cpp, $(OBJDIR)%.o, $(filter %.cpp, $(SRC)))

DIRS := $(addprefix $(OBJDIR), $(MODULES)) 

INCLUDE := 
CPPFLAGS += $(INCLUDE)

$(OBJDIR)%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(DIRS) $(OBJECTS)

$(DIRS):
	install -d $@
