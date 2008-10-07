#
# Makefile for CRAVA
#
include Makeheader

DIRS        = src lib nrlib fft/fftw fft/rfftw
OBJDIR      = obj
OBJSUBDIR   = obj/fft
OBJNRLIBDIR = obj/nrlib
INCLUDE     = -I. -I./fft/include
CPPFLAGS   += $(INCLUDE) 

all:	$(PROGRAM)

$(PROGRAM): $(DIRS) main.o
	$(PURIFY) $(CXX) $(LFLAGS) -o $@ $(OBJDIR)/*.o $(OBJSUBDIR)/*.o $(OBJNRLIBDIR)/*/*.o main.o

$(OBJDIR):
	install -d $(OBJDIR)

$(OBJSUBDIR):
	install -d $(OBJSUBDIR)

.PHONY: clean $(DIRS)

$(DIRS): $(OBJDIR) $(OBJSUBDIR)
	cd $@ && $(MAKE)

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(PROGRAM) main.o

veryclean:
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJSUBDIR)/*.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(PROGRAM) main.o

test:	$(PROGRAM) 
	cd test_suite; chmod +x TestScript.pl; perl -s ./TestScript.pl ../$(PROGRAM) $(case); cd ..
