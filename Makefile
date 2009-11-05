#
# Makefile for CRAVA
#
include Makeheader

DIRS        = src lib nrlib fft/fftw fft/rfftw
OBJDIR      = obj
OBJSUBDIR   = obj/fft
OBJNRLIBDIR = obj/nrlib
INCLUDE     = -I. -I./fft/include -I./nrlib
CPPFLAGS   += $(INCLUDE) 

all:	$(PROGRAM)

$(PROGRAM): $(DIRS) main.o
	$(PURIFY) $(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $(OBJDIR)/*.o $(OBJSUBDIR)/*.o $(OBJNRLIBDIR)/*/*.o $(OBJNRLIBDIR)/*/*/*.o main.o

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

cleanlib:
	rm -f $(OBJDIR)/*.o
	rm -f $(PROGRAM) main.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(OBJNRLIBDIR)/*/*/*.o

cleanall:
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJSUBDIR)/*.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(OBJNRLIBDIR)/*/*/*.o
	rm -f $(PROGRAM) main.o

test:	$(PROGRAM) 
	cd test_suite; chmod +x TestScript.pl; perl -s ./TestScript.pl ../$(PROGRAM) $(case); cd ..

help:
	@echo ''
	@echo 'Usage:  make type [mode=...] [case=...]'
	@echo ''
	@echo 'types'
	@echo '  clean     : Remove object files generated from  src'
	@echo '  cleanlib  : Remove object files generated from  src + NRLib'
	@echo '  cleanall  : Remove object files generated from  src + NRLib + fft'
	@echo '  test      : Run CRAVA in test suite'
	@echo '  all       : Make CRAVA'
	@echo ''
	@echo ''
	@echo 'modes'
	@echo '  debug     : Compile and link with -g -O0'
	@echo '  profile   : Compile and link with -g -pg'
	@echo '  purify    : Compile and link with -g -p0 and link with purify. Executable becomes cravarun.purify'
	@echo ''
	@echo 'cases'
	@echo '  n         : Test case number (number given first in the test case directory name)'
	@echo ''
