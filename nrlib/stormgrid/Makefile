include ../../Makeheader

INCLUDE =  
CPPFLAGS += $(INCLUDE)

OBJDIR = ../../obj
CSRCS = $(wildcard *.c)
CXXSRCS = $(wildcard *.cpp)
OBJECTS = $(CSRCS:%.c=$(OBJDIR)/%.o) $(CXXSRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

all: $(OBJECTS) 
