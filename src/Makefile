
include ../Makeheader

INCLUDE = -I.. -I../libs -I../libs/fft/include
CPPFLAGS += $(INCLUDE)

OBJDIR = ../obj
SRCS = $(wildcard *.cpp)
OBJECTS = $(SRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(OBJECTS)