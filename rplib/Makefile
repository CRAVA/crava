include ../Makeheader

INCLUDE =  -I.. -I ../fft/include
CPPFLAGS += $(INCLUDE) -I../flens

OBJDIR = ../obj
CXXSRCS = $(wildcard *.cpp)
OBJECTS = $(CXXSRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@


all: $(OBJECTS) 