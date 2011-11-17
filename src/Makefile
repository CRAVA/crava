
include ../Makeheader

INCLUDE = -I.. -I../fft/include -I../flens -I $(MKLINCLUDE)
CPPFLAGS += $(INCLUDE)

OBJDIR = ../obj
SRCS = $(wildcard *.cpp)
OBJECTS = $(SRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: $(OBJECTS)