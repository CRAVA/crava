These are the release notes for CRAVA

CRAVA is a software package for seismic inversion and conditioning of geological reservoir models. CRAVA is copyrighted by Norwegian Computing Centre and Statoil and licensed under GPLv3+. See COPYING for details.

ON WHAT PLATFORMS DOES IT RUN?
CRAVA is developep for GNU/Linux.

DOCUMENTATION
No efforts to provide documentation has been made yet.

DOWNLOADING CRAVA
git clone git://github.com/CRAVA/nrlib.git
git clone git://github.com/CRAVA/crava.git

BUILDING CRAVA
Make a softlink to nrlib, e.g.,
 ln -s PATH/nrlib/nrlib PATH/crava/lib
Compile CRAVA:
 cd PATH/crava
 make
If you have the test-suite available, unpack it to PATH/crava and run
 cd PATH/crava
 make test
