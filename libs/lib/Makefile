include ../../Makeheader

INCLUDE = -I../.. -I../../libs -I../../libs/fft/include
CPPFLAGS += $(INCLUDE)

OBJDIR = ../../obj/libs/lib

CSRCS = $(wildcard *.c)
CXXSRCS = $(wildcard *.cpp)
OBJECTS = $(CSRCS:%.c=$(OBJDIR)/%.o) $(CXXSRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

all: $(OBJDIR) $(OBJECTS)

$(OBJDIR):
	install -d $@
