#ifndef FFTFILEGRID_H
#define FFTFILEGRID_H

#include <string>
#include "fft/include/fftw.h"

#include "fftgrid.h"

class Corr;
class Wavelet;
class Simbox;
class GridMapping;

class FFTFileGrid : public FFTGrid
{ 
public:
  FFTFileGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp);
  FFTFileGrid(FFTFileGrid * FFTFileGrid);
  ~FFTFileGrid();

  fftw_complex getNextComplex() ;
  float        getNextReal() ;
  float        getRealValue(int i, int j, int k);
  float        getRealValueInterpolated(int i, int j, float kindex);
  int          setRealValue(int i, int j, int k, float value);
  int          setNextComplex(fftw_complex);
  int          setNextReal(float);
  int          square();
  int          expTransf();
  int          logTransf();
  void         multiplyByScalar(float scalar);
  int          collapseAndAdd(float*);
  void         add(FFTGrid* fftGrid);
  void         multiply(FFTGrid* fftGrid);              // pointwise multiplication! 
  void         fillInComplexNoise(RandomGen * ranGen);
  void         fftInPlace();	
  void         invFFTInPlace();
  void         createRealGrid();
  void         createComplexGrid();

  void         setAccessMode(int mode);
  void         endAccess();
  void         writeFile(const std::string & fileName, const Simbox * simbox, 
                         const std::string sgriLabel = "NO_LABEL", float z0 = 0.0, 
                         GridMapping * depthMap = NULL, GridMapping * timeMap = NULL);  //Use this instead of the ones below.
  void         writeStormFile(const std::string & fileName, const Simbox * simbox, bool expTrans = false,
                              bool ascii = false, bool padding = false, bool flat = false);
  int          writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0);
  int          writeSgriFile(const std::string & fileName, const Simbox *simbox, const std::string label);
  void         writeResampledStormCube(GridMapping *gridmapping, const std::string & fileName, 
                                       const Simbox *simbox, const int format, bool expTrans);
  void         writeDirectFile(const std::string & fileName);
  std::string  readDirectFile(const std::string & fileName);

  bool         isFile() {return(1);}

private:
  int accMode_;
  int modified_;   //Tells if grid is modified during RANDOMACCESS.
  char * fNameIn_; //Temporary names, switches whenever a write has occured.
  char * fNameOut_;
  FILE * inFile_;
  FILE * outFile_;

  static int gNum; //Number used for generating temporary files.

  void genFileName();
  void load();
  void unload();
  void save();
};
#endif
