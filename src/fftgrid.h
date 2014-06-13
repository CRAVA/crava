/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef FFTGRID_H
#define FFTGRID_H

#include <assert.h>
#include <complex>
#include <string>

#include "fftw.h"
#include "rfftw.h"
#include "definitions.h"

class Wavelet;
class Simbox;
class RandomGen;
class GridMapping;
class SeismicParametersHolder;

class FFTGrid
{
public:

  FFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp);
  FFTGrid(FFTGrid * fftGrid, bool expTrans = false);
  FFTGrid() {} //Dummy constructor needed for FFTFileGrid
  virtual ~FFTGrid();

  void setType(int cubeType) {cubetype_ = cubeType;}
  void setAngle(float angle) {theta_ = angle;}

  void                 fillInData(const Simbox  * timeSimbox,
                                  StormContGrid * grid,
                                  const SegY   *  segy,
                                  float           smooth_length,
                                  int           & missingTracesSimbox,
                                  int           & missingTracesPadding,
                                  int           & deadTracesSimbox,
                                  std::string   & errTxt,
                                  bool            scale = false,
                                  bool            is_segy = true);
  void                 smoothTraceInGuardZone(std::vector<float> & data_trace,
                                              //float                z0_data,
                                              //float                zn_data,
                                              float                dz_data,
                                              float                smooth_length);
                                              //std::string        & errTxt);
  void                 resampleTrace(const std::vector<float> & data_trace,
                                     const rfftwnd_plan       & fftplan1,
                                     const rfftwnd_plan       & fftplan2,
                                     fftw_real                * rAmpData,
                                     fftw_real                * rAmpFine,
                                     int                        cnt,
                                     int                        rnt,
                                     int                        cmt,
                                     int                        rmt);

  double               InterpolateTrilinear(double  x_min,
                                   double           x_max,
                                   double           y_min,
                                   double           y_max,
                                   double           z_min,
                                   double           z_max,
                                   double           x,
                                   double           y,
                                   double           z);

  double              InterpolateBilinearXY(double x_min,
                                            double x_max,
                                            double y_min,
                                            double y_max,
                                            double x,
                                            double y);

  void                 interpolateGridValues(std::vector<float> & grid_trace,
                                             float                z0_grid,
                                             float                dz_grid,
                                             fftw_real          * rAmpFine,
                                             float                z0_data,
                                             float                dz_fine,
                                             int                  n_fine);
  void                 interpolateAndShiftTrend(std::vector<float>       & interpolated_trend,
                                                float                      z0_grid,
                                                float                      dz_grid,
                                                const std::vector<float> & trend_long,
                                                float                      z0_data,
                                                float                      dz_fine,
                                                int                        n_fine);
  void                 setTrace(const std::vector<float> & trace, size_t i, size_t j);
  void                 setTrace(float value, size_t i, size_t j);

  void                 fillInConstant(float value, bool add = true);              // No mode

  void                 fillInErrCorr(const Surface * priorCorrXY,
                                     float           gradI,
                                     float           gradJ);

  void                 fillInParamCorr(const Surface   * priorCorrXY,
                                       const fftw_real * circCorrT,
                                       float             gradI,
                                       float             gradJ);// No mode
  void                fillInTimeCov(const fftw_real * circCorrT);


  void                 fillInGenExpCorr(double Rx,double Ry,double Rz,
                          float             gradI,
                          float             gradJ);// No mode



  virtual void         fillInComplexNoise(RandomGen * ranGen);   // No mode/randomaccess

  void                 fillInFromArray(float *value);
  void                 calculateStatistics();                    // min,max, avg
  void                 setUndefinedCellsToGlobalAverage();      // For BG model

  virtual fftw_complex getNextComplex() ;                       // Accessmode read/readandwrite
  virtual float        getNextReal() ;                          // Accessmode read/readandwrite
  virtual int          setNextComplex(fftw_complex);            // Accessmode write/readandwrite
  virtual int          SetNextComplex(std::complex<double> & v);// Accessmode write/readandwrite
  virtual int          setNextReal(float);                      // Accessmode write/readandwrite
  float                getRealValue(int i, int j, int k, bool extSimbox = false) const;  // Accessmode randomaccess
  float                getRealValueCyclic(int i, int j, int k);
  float                getRealValueInterpolated(int i, int j, float kindex, bool extSimbox = false);
  fftw_complex         getComplexValue(int i, int j, int k, bool extSimbox = false) const;
  virtual int          setRealValue(int i, int j, int k, float value, bool extSimbox = false);  // Accessmode randomaccess
  int                  setComplexValue(int i, int j ,int k, fftw_complex value, bool extSimbox = false);
  fftw_complex         getFirstComplexValue();
  float                getFirstRealValue();                     // No mode/randomaccess
  virtual int          square();                                // No mode/randomaccess
  virtual int          expTransf();                             // No mode/randomaccess
  virtual int          logTransf();                             // No mode/randomaccess
  virtual void         realAbs();
  virtual int          collapseAndAdd(float* grid);             // No mode/randomaccess
  virtual void         fftInPlace();                            // No mode/randomaccess
  virtual void         invFFTInPlace();                         // No mode/randomaccess


  virtual void         add(FFTGrid* fftGrid);                   // No mode/randomaccess
  virtual void         addScalar(float scalar);                 // No mode/randomaccess, only for real grids
  virtual void         subtract(FFTGrid* fftGrid);              // No mode/randomaccess
  virtual void         changeSign();                   // No mode/randomaccess
  virtual void         multiply(FFTGrid* fftGrid);              // pointwise multiplication!
  virtual void         conjugate();                             // No mode/randomaccess
  bool                 consistentSize(int nx,int ny, int nz, int nxp, int nyp, int nzp);
  int                  getCounterForGet() const {return(counterForGet_);}
  int                  getCounterForSet() const {return(counterForSet_);}
  int                  getNx()      const {return(nx_);}
  int                  getNy()      const {return(ny_);}
  int                  getNz()      const {return(nz_);}
  int                  getNxp()     const {return(nxp_);}
  int                  getNyp()     const {return(nyp_);}
  int                  getNzp()     const {return(nzp_);}
  int                  getRNxp()    const {return(rnxp_);}
  int                  getCNxp()    const {return(cnxp_);}
  int                  getcsize()   const {return(csize_);}
  int                  getrsize()   const {return(rsize_);}
  float                getTheta()   const {return(theta_);}
  float                getScale()   const {return(scale_);}
  float                getMinReal() const {return rValMin_;}
  float                getMaxReal() const {return rValMax_;}
  float                getAvgReal() const {return rValAvg_;}

  bool                 getIsTransformed() const {return(istransformed_);}
  //For use when writing to a grid that may be in the wrong state.
  void                 setTransformedStatus(bool status) {istransformed_ = status;}

  enum                 gridTypes{CTMISSING, DATA, PARAMETER, COVARIANCE, VELOCITY};
  enum                 accessMode{NONE, READ, WRITE, READANDWRITE, RANDOMACCESS};

  virtual void         multiplyByScalar(float scalar);      //No mode/randomaccess
  int                  getType() const {return(cubetype_);}
  virtual void         setAccessMode(int mode){assert(mode>=0);}
  virtual void         endAccess(){counterForGet_ = 0; counterForSet_ = 0;}
  virtual void         writeFile(const std::string              & fileName,
                                 const std::string              & subDir,
                                 const Simbox                   * simbox,
                                 const std::string                sgriLabel = "NO_LABEL",
                                 const float                      z0        = 0.0,
                                 const GridMapping              * depthMap  = NULL,
                                 const GridMapping              * timeMap   = NULL,
                                 const TraceHeaderFormat        & thf       = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS),
                                 bool                             padding   = false,
                                 bool                             scientific_format = false,
                                 const std::vector<std::string> & headerText = std::vector<std::string>());
  //Use this instead of the ones below.
  virtual void         writeStormFile(const std::string & fileName, const Simbox * simbox, bool ascii = false,
                                      bool padding = false, bool flat = false, bool scientific_format = false);//No mode/randomaccess
  virtual int          writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0,
                                     const TraceHeaderFormat &thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS),
                                     const std::vector<std::string> & headerText = std::vector<std::string>());   //No mode/randomaccess
  virtual int          writeSgriFile(const std::string & fileName, const Simbox * simbox, const std::string label);
  virtual void         writeAsciiFile(const std::string & fileName);
  virtual void         writeAsciiRaw(const std::string & fileName);
  virtual void         writeResampledStormCube(const GridMapping *gridmapping, const std::string & fileName,
                                               const Simbox *simbox, const int format);
  virtual void         writeCravaFile(const std::string & fileName, const Simbox * simbox);
  virtual void         readCravaFile(const std::string & fileName, std::string & errText, bool nopadding = false);

  virtual bool         isFile() {return(0);}    // indicates wether the grid is in memory or on disk

  static void          setOutputFlags(int format, int domain) {formatFlag_ = format;domainFlag_=domain;};
  static void          setOutputFormat(int format) {formatFlag_ = format;}
  int                  getOutputFormat() {return(formatFlag_);}
  static void          setOutputDomain(int domain) {domainFlag_ = domain;}
  int                  getOutputDomain() {return(domainFlag_);}
  static void          setMaxAllowedGrids(int maxAllowedGrids) {maxAllowedGrids_ = maxAllowedGrids ;}
  static int           getMaxAllowedGrids()   { return maxAllowedGrids_   ;}
  static int           getMaxAllocatedGrids() { return maxAllocatedGrids_ ;}
  static void          setTerminateOnMaxGrid(bool terminate) {terminateOnMaxGrid_ = terminate ;}
  static int           findClosestFactorableNumber(int leastint);

  static fftw_complex* fft1DzInPlace(fftw_real*  in, int nzp);
  static fftw_real*    invFFT1DzInPlace(fftw_complex* in, int nzp);

  virtual void         createRealGrid(bool add = true);
  virtual void         createComplexGrid();

  //This function interpolates seismic in all missing traces inside area, and mirrors to padding.
  //Also interpolates in traces where energy is lower than treshold.
  virtual void         interpolateSeismic(float energyTreshold = 0);

  void                 checkNaN(); //NBNB Ragnar: For debug purpose. Negative number = OK.
  float                getDistToBoundary(int i, int n , int np);
  virtual void         getRealTrace(float * value, int i, int j);
  virtual int          setRealTrace(int i, int j, float *value);
  std::vector<float>   getRealTrace2(int i, int j) const;


  static void          reportFFTMemoryAndWait(const std::string & msg) {
                         LogKit::LogFormatted(LogKit::High, "%s: %2d grids, %10.2f MB\n", msg.c_str(), nGrids_, FFTMemUse_/(1024.0f*1024.0f));
                         float tmp;
                         std::cin >> tmp;
                         LogKit::LogFormatted(LogKit::High, "Memory used %4.0f MB, used outside grid %4.0f MB\n", tmp, tmp-FFTMemUse_/(1024.0f*1024.0f));
                       }

  void                 createGrid();
protected:
  //int                setPaddingSize(int n, float p);
  int                  getFillNumber(int i, int n, int np );

  int                  getXSimboxIndex(int i) { return (getFillNumber(i, nx_, nxp_ )) ;}
  int                  getYSimboxIndex(int j) { return (getFillNumber(j, ny_, nyp_ )) ;}
  int                  getZSimboxIndex(int k);

  //Interpolation into SegY and sgri
  float                getRegularZInterpolatedRealValue(int i, int j, double z0Reg,
                                                         double dzReg, int kReg,
                                                         double z0Grid, double dzGrid);

  //Supporting functions for interpolateSeismic
  int                  interpolateTrace(int index, short int * flags, int i, int j);
  void                 extrapolateSeismic(int imin, int imax, int jmin, int jmax);

  /// Called from writeResampledStormCube
  void                 writeSegyFromStorm(Simbox * simbox, StormContGrid *data, std::string fileName);
  void                 makeDepthCubeForSegy(Simbox *simbox,const std::string & fileName);

  int                  cubetype_;          // see enum gridtypes above
  float                theta_;             // angle in angle gather (case of data)
  float                scale_;             // To keep track of the scalings after fourier transforms
  int                  nx_;                // size of original grid in lateral x direction
  int                  ny_;                // size of original grid in lateral y direction
  int                  nz_;                // size of original grid in depth (time)
  int                  nxp_;               // size of padded FFT grid in lateral x direction
  int                  nyp_;               // size of padded FFT grid in lateral y direction
  int                  nzp_;               // size of padded FFT grid in depth (time)

  int                  cnxp_;              // size in x direction for storage inplace algorithm (complex grid) nxp_/2+1
  int                  rnxp_;              // expansion in x direction for storage inplace algorithm (real grid) 2*(nxp_/2+1)

  int                  csize_;             // size of complex grid, cnxp_*nyp_*nzp_
  int                  rsize_;             // size of real grid rnxp_*nyp_*nzp_
  int                  counterForGet_;     // active cell in grid
  int                  counterForSet_;     // active cell in grid

  bool                 istransformed_;     // true if the grid contain Fourier values (i.e complex variables)

  fftw_complex       * cvalue_;            // values of complex parameter in grid points
  fftw_real          * rvalue_;            // values of real parameter in grid points

  float                rValMin_;           // minimum real value
  float                rValMax_;           // maximum real value
  float                rValAvg_;           // average real value

  static int           formatFlag_;        // Decides format of output (see ModelSettings).
  static int           domainFlag_;        // Decides domain of output (see ModelSettings).

  static int           maxAllowedGrids_;   // The maximum number of grids we are allowed to allocate.
  static int           maxAllocatedGrids_; // The maximum number of grids that has actually been allocated.
  static int           nGrids_;            // The actually number of grids allocated (varies as crava runs).
  static bool          terminateOnMaxGrid_; // If true, terminate when we try to allocate more than maxAllowedGrids.
  bool                 add_;                // Tells whether we should change nGrids_ or not

  static float         maxFFTMemUse_;
  static float         FFTMemUse_;

};
#endif
