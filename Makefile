#
# Makefile for CRAVA
#
include Makeheader

DIRS        = src libs libs/lib libs/boost libs/fft/fftw libs/fft/rfftw
OBJDIR      = obj
OBJLIBDIR   = obj/libs/lib
OBJFFTDIR   = obj/libs/fft
OBJBOOSTDIR = obj/libs/boost
OBJNRLIBDIR = obj/libs/nrlib
OBJFINDGRAM = findgrammar/findgrammar.o
OBJGRAMMAR  = $(OBJNRLIBDIR)/iotools/fileio.o         \
              $(OBJNRLIBDIR)/tinyxml/tinyxml.o        \
              $(OBJNRLIBDIR)/tinyxml/tinyxmlerror.o   \
              $(OBJNRLIBDIR)/tinyxml/tinyxmlparser.o  \
              $(OBJBOOSTDIR)/system/error_code.o      \
              $(OBJBOOSTDIR)/filesystem/path.o        \
              $(OBJBOOSTDIR)/filesystem/operations.o  \
              $(OBJBOOSTDIR)/filesystem/portability.o
INCLUDE     = -I. -I./libs -I./libs/nrlib -I./libs/fft/include
CPPFLAGS   += $(INCLUDE)

all:	$(PROGRAM)

$(PROGRAM): $(DIRS) main.o
	$(PURIFY) $(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $(OBJDIR)/*.o $(OBJLIBDIR)/*.o $(OBJFFTDIR)/*.o $(OBJBOOSTDIR)/*/*.o $(OBJNRLIBDIR)/*/*.o  main.o

$(GRAMMAR): findgrammar/findgrammar.o
	$(PURIFY) $(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $(OBJGRAMMAR) $(OBJFINDGRAM)

$(OBJDIR):
	install -d $(OBJDIR)

$(OBJFFTDIR):
	install -d $(OBJFFTDIR)

.PHONY: clean $(DIRS)

$(DIRS): $(OBJDIR) $(OBJFFTDIR)
	cd $@ && $(MAKE)

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(PROGRAM) main.o

cleanlib:
	rm -f $(OBJDIR)/*.o
	rm -f $(PROGRAM) main.o
	rm -f $(OBJLIBDIR)/*.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(OBJBOOSTDIR)/*/*.o

cleanall:
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJLIBDIR)/*.o
	rm -f $(OBJFFTDIR)/*.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(OBJBOOSTDIR)/*/*.o
	rm -f $(PROGRAM) $(GRAMMAR) $(OBJFINDGRAM) main.o

test:	$(PROGRAM) $(GRAMMAR)
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
	@echo 'modes'
	@echo '  debug     : Compile and link with -g -O0'
	@echo '  profile   : Compile and link with -g -pg'
	@echo '  purify    : Compile and link with -g -p0 and link with purify. Executable becomes cravarun.purify'
	@echo ''
	@echo 'cases'
	@echo '  n         : Comma-separated list of test case number (number given first in the test case directory name)'
	@echo ''
