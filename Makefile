#
# Makefile for CRAVA
#
include Makeheader

DIRS        = src libs libs/lib libs/boost libs/fft/fftw libs/fft/rfftw rplib
OBJDIR      = obj
OBJLIBDIR   = obj/libs/lib
OBJFFTDIR   = obj/libs/fft
OBJFLENSDIR = obj/libs/flens
OBJBOOSTDIR = obj/libs/boost
OBJNRLIBDIR = obj/libs/nrlib

OBJGRAMMAR  = findgrammar/findgrammar.o               \
	      $(OBJNRLIBDIR)/iotools/fileio.o         \
              $(OBJNRLIBDIR)/tinyxml/tinyxml.o        \
              $(OBJNRLIBDIR)/tinyxml/tinyxmlerror.o   \
              $(OBJNRLIBDIR)/tinyxml/tinyxmlparser.o  \
              $(OBJBOOSTDIR)/system/error_code.o      \
              $(OBJBOOSTDIR)/filesystem/path.o        \
              $(OBJBOOSTDIR)/filesystem/operations.o  \
              $(OBJBOOSTDIR)/filesystem/portability.o

OBJCOMPARE  = compare_storm_binary_volumes/compare.o   \
	      $(OBJNRLIBDIR)/iotools/fileio.o          \
              $(OBJNRLIBDIR)/iotools/stringtools.o     \
              $(OBJNRLIBDIR)/stormgrid/stormcontgrid.o \
              $(OBJNRLIBDIR)/volume/volume.o           \
              $(OBJNRLIBDIR)/surface/surfaceio.o       \
              $(OBJNRLIBDIR)/math/constants.o          \
              $(OBJBOOSTDIR)/system/error_code.o       \
              $(OBJBOOSTDIR)/filesystem/path.o         \
              $(OBJBOOSTDIR)/filesystem/operations.o   \
              $(OBJBOOSTDIR)/filesystem/portability.o

INCLUDE     = -I. -I./libs -I./libs/nrlib -I./libs/flens -I./libs/fft/include
CPPFLAGS   += $(INCLUDE)

all:	$(PROGRAM)

comp:	$(COMPARE)

$(PROGRAM): $(DIRS) main.o
	$(PURIFY) $(CXX) $(OBJDIR)/*.o $(OBJLIBDIR)/*.o $(OBJNRLIBDIR)/*/*.o $(OBJFFTDIR)/*.o $(OBJBOOSTDIR)/*/*.o $(OBJFLENSDIR)/*.o main.o $(LFLAGS) -o $@

$(GRAMMAR): findgrammar/findgrammar.o
	$(PURIFY) $(CXX) $(OBJGRAMMAR) $(LFLAGS) -o $@

$(COMPARE): compare_storm_binary_volumes/compare.o
	$(PURIFY) $(CXX) $(OBJCOMPARE) $(LFLAGS) -o $@

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
	rm -f $(OBJLIBDIR)/*.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(OBJBOOSTDIR)/*/*.o
	rm -f $(OBJFLENSDIR)/*.o
	rm -f $(PROGRAM) main.o

cleanall:
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJLIBDIR)/*.o
	rm -f $(OBJFFTDIR)/*.o
	rm -f $(OBJNRLIBDIR)/*/*.o
	rm -f $(OBJBOOSTDIR)/*/*.o
	rm -f $(OBJFLENSDIR)/*.o
	rm -f $(OBJGRAMMAR)/*.o
	rm -f $(OBJCOMPARE)/*.o
	rm -f $(GRAMMAR) findgrammar/findgrammar.o
	rm -f $(COMPARE) compare_storm_binary_volumes/compare.o
	rm -f $(PROGRAM) main.o

test:	$(PROGRAM) $(GRAMMAR) $(COMPARE)
	cd test_suite; chmod +x TestScript.pl; perl -s ./TestScript.pl ../$(PROGRAM) $(passive) $(case); cd ..

help:
	@echo ''
	@echo 'Usage:  make type [mode=...] [case=...] [passive=...] [at=...]'
	@echo ''
	@echo 'types'
	@echo '  clean     : Remove object files generated from  src'
	@echo '  cleanlib  : Remove object files generated from  src + boost + flens + NRLib'
	@echo '  cleanall  : Remove object files generated from  src + boost + flens + NRLib + fft'
	@echo '  test      : Run CRAVA in test suite'
	@echo '  all       : Make CRAVA'
	@echo ''
	@echo 'modes'
	@echo '  debug     : Compile and link with -g -O0'
	@echo '  profile   : Compile and link with -g -pg'
	@echo '  purify    : Compile and link with -g -p0 and link with purify. Executable becomes cravarun.purify'
	@echo ''
	@echo 'passive'
	@echo '  true      : Run test suite in passive mode (do not run CRAVA)'
	@echo ''
	@echo 'case'
	@echo '  n         : Comma-separated list of test case numbers (number given first in the test case'
	@echo '              directory name) or a range give as 1-5'
	@echo ''
	@echo 'at'
	@echo '  nr        : Needed to compile under Ubuntu at NR'
	@echo ''
