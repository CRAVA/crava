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
  float        getRealValue(int i, int j, int k, bool extSimbox = false);
  float        getRealValueInterpolated(int i, int j, float kindex);
  int          setRealValue(int i, int j, int k, float value, bool extSimbox = false);
  int          setNextComplex(fftw_complex);
  int          setNextReal(float);
  int          square();
  int          expTransf();
  int          logTransf();
  void         multiplyByScalar(float scalar);
  int          collapseAndAdd(float*);
  void         add(FFTGrid* fftGrid);
  void         subtract(FFTGrid* fftGrid);
  void         changeSign();
  void         multiply(FFTGrid* fftGrid);              // pointwise multiplication! 
  void         fillInComplexNoise(RandomGen * ranGen);
  void         fftInPlace();
  void         invFFTInPlace();
  void         createRealGrid(bool add = true);
  void         createComplexGrid();
  void         setAccessMode(int mode);
  void         endAccess();
  void         writeFile(const std::string & fileName, 
                         const std::string & subDir, 
                         const Simbox      * simbox, 
                         const std::string   sgriLabel = "NO_LABEL", 
                         const float         z0        = 0.0, 
                         GridMapping       * depthMap  = NULL, 
                         GridMapping       * timeMap   = NULL,
                         const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));  //Use this instead of the ones below.
  void         writeStormFile(const std::string & fileName, const Simbox * simbox, bool expTrans = false,
                              bool ascii = false, bool padding = false, bool flat = false);
  int          writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0, 
                             const TraceHeaderFormat &thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));
  int          writeSgriFile(const std::string & fileName, const Simbox *simbox, const std::string label);
  void         writeResampledStormCube(GridMapping *gridmapping, const std::string & fileName, 
                                       const Simbox *simbox, const int format, bool expTrans);
  void         writeCravaFile(const std::string & fileName, const Simbox * simbox);
  void         readCravaFile(const std::string & fileName, std::string & error, bool nopadding = false);

  bool         isFile() {return(1);}
  float        *getRealTrace(int i, int j);
  int          setRealTrace(int i, int j, float *value);
  void         fillInFromRealFFTGrid(FFTGrid& fftGrid);
private:
  void         genFileName();
  void         load();
  void         unload();
  void         save();

  int          accMode_;
  int          modified_;   //Tells if grid is modified during RANDOMACCESS.
  char       * fNameIn_; //Temporary names, switches whenever a write has occured.
  char       * fNameOut_;
  FILE       * inFile_;
  FILE       * outFile_;

  static int   gNum; //Number used for generating temporary files.
};
#endif
