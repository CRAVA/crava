#ifndef FFTGRID_H
#define FFTGRID_H

#include <assert.h>
#include <string>

#include "fft/include/fftw.h"
#include "definitions.h"

class Corr;
class Wavelet;
class Simbox;
class RandomGen;
class GridMapping;

class FFTGrid{
public:

  FFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp);
  FFTGrid(FFTGrid * fftGrid);
  FFTGrid() {} //Dummy constructor needed for FFTFileGrid
  virtual ~FFTGrid();

  void setType(int cubeType) {cubetype_ = cubeType;}
  void setAngle(float angle) {theta_ = angle;}


  void                 fillInFromSegY(SegY * segy, Simbox *simbox );            // No mode

  void                 fillInFromStorm(Simbox     * actSimBox,
                                       StormContGrid      * grid,
                                       const char * parName);    // No mode
  void                 fillInFromRealFFTGrid(FFTGrid& fftGrid);  // No mode
  void                 fillInConstant(float value);              // No mode
  fftw_real*           fillInParamCorr(Corr* corr,int minIntFq,
                                       float gradI, float gradJ);// No mode
  void                 fillInErrCorr(Corr* parCorr,              // No mode
                                     float gradI, float gradJ);  // No mode
  virtual void         fillInComplexNoise(RandomGen * ranGen);   // No mode/randomaccess

  void                 fillInTest(float value1, float value2);   // No mode /DEBUG
  void                 fillInFromArray(float *value);

  virtual fftw_complex getNextComplex() ;                       // Accessmode read/readandwrite
  virtual float        getNextReal() ;                          // Accessmode read/readandwrite
  virtual int          setNextComplex(fftw_complex);            // Accessmode write/readandwrite
  virtual int          setNextReal(float);                      // Accessmode write/readandwrite
  float                getRealValue(int i, int j, int k, bool extSimbox = false);  // Accessmode randomaccess
  float                getRealValueCyclic(int i, int j, int k);
  float                getRealValueInterpolated(int i, int j, float kindex, bool extSimbox = false);
  fftw_complex         getComplexValue(int i, int j, int k, bool extSimbox = false) const;
  int                  setRealValue(int i, int j, int k, float value, bool extSimbox = false);  // Accessmode randomaccess
  int                  setComplexValue(int i, int j ,int k, fftw_complex value, bool extSimbox = false);
  fftw_complex         getFirstComplexValue();    
  virtual int          square();                                // No mode/randomaccess
  virtual int          expTransf();                             // No mode/randomaccess
  virtual int          logTransf();                             // No mode/randomaccess
  virtual void         realAbs();
  virtual int          collapseAndAdd(float* grid);             // No mode/randomaccess
  virtual void         fftInPlace();	                        // No mode/randomaccess
  virtual void         invFFTInPlace();                         // No mode/randomaccess


  virtual void         add(FFTGrid* fftGrid);                   // No mode/randomaccess
  virtual void         subtract(FFTGrid* fftGrid);                   // No mode/randomaccess
  virtual void         changeSign();                   // No mode/randomaccess
  virtual void         multiply(FFTGrid* fftGrid);              // pointwise multiplication!    
  bool                 consistentSize(int nx,int ny, int nz, int nxp, int nyp, int nzp);
  int                  getCounterForGet() const {return(counterForGet_);}
  int                  getCounterForSet() const {return(counterForSet_);}
  int                  getNx()    const {return(nx_);}
  int                  getNy()    const {return(ny_);}
  int                  getNz()    const {return(nz_);}
  int                  getNxp()   const {return(nxp_);}
  int                  getNyp()   const {return(nyp_);}
  int                  getNzp()   const {return(nzp_);}
  int                  getRNxp()  const {return(rnxp_);}
  int                  getcsize() const {return(csize_);}
  int                  getrsize() const {return(rsize_);}
  float                getTheta() const {return(theta_);}  
  float                getScale() const {return(scale_);}
  bool                 getIsTransformed() const {return(istransformed_);}
  enum                 gridTypes{CTMISSING,DATA,PARAMETER,COVARIANCE,VELOCITY};
  enum                 accessMode{NONE,READ,WRITE,READANDWRITE,RANDOMACCESS};
  virtual void         multiplyByScalar(float scalar);      //No mode/randomaccess
  int                  getType() const {return(cubetype_);}
  virtual void         setAccessMode(int mode){assert(mode>=0);}
  virtual void         endAccess(){}
  virtual void         writeFile(const std::string & fileName, const Simbox * simbox, 
                                 const std::string sgriLabel = "NO_LABEL", float z0 = 0.0, 
                                 GridMapping * depthMap = NULL, GridMapping * timeMap = NULL); 
  //Use this instead of the ones below.
  virtual void         writeStormFile(const std::string & fileName, const Simbox * simbox, bool expTrans = false, 
                                      bool ascii = false, bool padding = false, bool flat = false);//No mode/randomaccess
  virtual int          writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0);   //No mode/randomaccess
  virtual int          writeSgriFile(const std::string & fileName, const Simbox * simbox, const std::string label);
  virtual void         writeAsciiFile(const std::string & fileName);
  virtual void         writeAsciiRaw(const std::string & fileName);
  virtual void         writeResampledStormCube(GridMapping *gridmapping, const std::string & fileName, 
                                               const Simbox *simbox, const int format, bool expTrans);
  virtual void         writeDirectFile(const std::string & fileName, const Simbox * simbox);
  virtual std::string  readDirectFile(const std::string & fileName);

  virtual bool         isFile() {return(0);}    // indicates wether the grid is in memory or on disk  

  static void          setOutputFlags(int format, int domain) {formatFlag_ = format;domainFlag_=domain;};
  static void          setOutputFormat(int format) {formatFlag_ = format;} 
  int                  getOutputFormat() {return(formatFlag_);} 
  static void          setOutputDomain(int domain) {domainFlag_ = domain;} 
  int                  getOutputDomain() {return(domainFlag_);} 

  virtual void         createRealGrid();
  virtual void         createComplexGrid();

  //This function interpolates seismic in all missing traces inside area, and mirrors to padding.
  //Also interpolates in traces where energy is lower than treshold.
  virtual void         interpolateSeismic(float energyTreshold = 0);

  void                 checkNaN(); //NBNB Ragnar: For debug purpose. Negative number = OK.
  float                getDistToBoundary(int i, int n , int np); 
  float               *getRealTrace(int i, int j);
  int                  setRealTrace(int i, int j, float *value);
protected:
  int    cubetype_;        // see enum gridtypes above
  float  theta_;           // angle in angle gather (case of data)
  float  scale_;           // To keep track of the scalings after fourier transforms
  int    nx_;              // size of original grid in depth (time)
  int    ny_;              // size of original grid in lateral x direction 
  int    nz_;              // size of original grid in lateral y direction
  int    nxp_;             // size of padded FFT grid in depth (time) 
  int    nyp_;             // size of padded FFT grid in lateral x direction 
  int    nzp_;             // size of padded FFT grid in lateral y direction

  int    cnxp_;	           // size in x direction for storage inplace algorithm (complex grid) nxp_/2+1
  int    rnxp_;            // expansion in x direction for storage inplace algorithm (real grid) 2*(nxp_/2+1)

  int    csize_;           // size of complex grid, cnxp_*nyp_*nzp_
  int    rsize_;           // size of real grid rnxp_*nyp_*nzp_
  int    counterForGet_;   // active cell in grid
  int    counterForSet_;   // active cell in grid

  bool   istransformed_;   // true if the grid contain Fourier values (i.e complex variables)

  fftw_complex* cvalue_;   // values of complex parameter in grid points
  fftw_real*    rvalue_;   // values of real parameter in grid points

  static int  formatFlag_; // Decides format of output (see ModelSettings).
  static int  domainFlag_; // Decides domain of output (see ModelSettings).

  //int                 setPaddingSize(int n, float p); 
  int                   getFillNumber(int i, int n, int np );

  int                   getXSimboxIndex(int i){return(getFillNumber(i, nx_, nxp_ ));}
  int                   getYSimboxIndex(int j){return(getFillNumber(j, ny_, nyp_ ));}
  int                   getZSimboxIndex(int k);
  void                  computeCircCorrT(Corr* corr,fftw_real* CircCorrT);
  void                  makeCircCorrTPosDef(fftw_real* CircularCorrT,int minIntFq);
  fftw_complex*         fft1DzInPlace(fftw_real*  in);
  fftw_real*            invFFT1DzInPlace(fftw_complex* in);

  //Interpolation into SegY and sgri
  float                 getRegularZInterpolatedRealValue(int i, int j, double z0Reg,
                                                         double dzReg, int kReg,
                                                         double z0Grid, double dzGrid);

  //Supporting functions for interpolateSeismic
  int                   interpolateTrace(int index, short int * flags, int i, int j);
  void                  extrapolateSeismic(int imin, int imax, int jmin, int jmax);

  /// Called from writeResampledStormCube
  void writeSegyFromStorm(StormContGrid *data, std::string fileName);
  void makeDepthCubeForSegy(Simbox *simbox,const std::string & fileName);

};
#endif
