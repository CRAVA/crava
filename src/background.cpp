#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#include "lib/lib_misc.h"
#include "lib/global_def.h"
#include "lib/kriging1d.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/model.h"
#include "src/modelsettings.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/welldata.h"
#include "src/background.h"
#include "src/covgrid2d.h"
#include "src/krigingdata2d.h"
#include "src/kriging2d.h"
#include "src/krigingdata3d.h"
#include "src/covgridseparated.h"
#include "src/krigingadmin.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/io.h"

Background::Background(FFTGrid       ** grids,
                       WellData      ** wells,
                       FFTGrid       *& velocity,
                       Simbox         * timeSimbox,
                       Simbox         * timeBGSimbox,
                       ModelSettings  * modelSettings)
  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
    vsvp_(RMISSING)
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];

  FFTGrid * bgAlpha;
  FFTGrid * bgBeta;
  FFTGrid * bgRho;
  
  if (timeBGSimbox == NULL)
  {
    generateBackgroundModel(bgAlpha,bgBeta,bgRho,
                            velocity,wells,
                            timeSimbox,
                            modelSettings);
  }
  else
  {
    generateBackgroundModel(bgAlpha,bgBeta,bgRho,
                            velocity,wells,
                            timeBGSimbox,
                            modelSettings);
    resampleBackgroundModel(bgAlpha,bgBeta,bgRho,
                            timeBGSimbox,
                            timeSimbox,
                            modelSettings);
  }  
  padAndSetBackgroundModel(bgAlpha,bgBeta,bgRho);    

  delete bgAlpha;
  delete bgBeta;
  delete bgRho;

  findMeanVsVp(backModel_[0],
               backModel_[1]);
}

//-------------------------------------------------------------------------------
Background::Background(FFTGrid ** grids) 
  : DataTarget_(IMISSING),
    vsvp_(RMISSING)
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];
  findMeanVsVp(backModel_[0],
               backModel_[1]);
}

//-------------------------------------------------------------------------------
Background::~Background(void)
{
  //
  // Hmmm, denne gir "*** glibc detected *** ./crava.exe: double free or corruption ***"
  //
  //for (int i=0 ; i<3 ; i++)
  //  delete backModel_[i];
  //delete backModel_;
}

//-------------------------------------------------------------------------------
void
Background::generateBackgroundModel(FFTGrid      *& bgAlpha,
                                    FFTGrid      *& bgBeta,
                                    FFTGrid      *& bgRho,
                                    FFTGrid      *& velocity,
                                    WellData     ** wells,
                                    Simbox        * simbox,
                                    ModelSettings * modelSettings)
{
  const int nz = simbox->getnz();
  float * trendAlpha = new float[nz];
  float * trendBeta  = new float[nz];
  float * trendRho   = new float[nz];
  float * trendVel   = new float[nz]; // Allocate (for simplicity) although not always needed

  const int nWells = modelSettings->getNumberOfWells();
  float * avgDevAlpha = new float[nWells];
  float * avgDevBeta  = new float[nWells];
  float * avgDevRho   = new float[nWells];
  float * avgDevVel   = new float[nWells]; // Allocate (for simplicity) although not always needed

  bool hasVelocityTrend = velocity != NULL;
  bool write1D = ((modelSettings->getOtherOutputFlag()& ModelSettings::BACKGROUND_TREND_1D) > 0);
  bool write3D = ((modelSettings->getGridOutputFlag() & ModelSettings::BACKGROUND_TREND) > 0);
  calculateBackgroundTrend(trendAlpha, avgDevAlpha,
                           wells, simbox, 
                           modelSettings->getAlphaMin(), 
                           modelSettings->getAlphaMax(),
                           modelSettings->getMaxHzBackground(), 
                           write1D, write3D,
                           nWells, hasVelocityTrend,
                           std::string("Vp"));
  calculateBackgroundTrend(trendBeta, avgDevBeta, 
                           wells, simbox, 
                           modelSettings->getBetaMin(), 
                           modelSettings->getBetaMax(),
                           modelSettings->getMaxHzBackground(), 
                           write1D, write3D,
                           nWells, hasVelocityTrend,
                           std::string("Vs"));
  calculateBackgroundTrend(trendRho, avgDevRho,
                           wells, simbox, 
                           modelSettings->getRhoMin(), 
                           modelSettings->getRhoMax(),
                           modelSettings->getMaxHzBackground(), 
                           write1D, write3D,
                           nWells, hasVelocityTrend,
                           std::string("Rho"));

  if (velocity != NULL) {
    //
    // We still want calculateBackgroundTrend() for alpha above. By calculating
    // avgDevAlpha we can check that the bgAlpha calculated from velocity is as 
    // good as or better than that calculated by crava.
    //
    calculateVelocityDeviations(velocity, wells, simbox, 
                                trendVel, avgDevVel, avgDevAlpha,
                                modelSettings->getGridOutputFlag(),
                                nWells);
    velocity->logTransf();
    delete bgAlpha;
    bgAlpha = velocity;
    velocity = NULL;
    writeDeviationsFromVerticalTrend(avgDevVel, avgDevBeta, avgDevRho, 
                                     trendVel, trendBeta, trendRho, 
                                     wells, nWells, nz);
  }
  else {
    writeDeviationsFromVerticalTrend(avgDevAlpha, avgDevBeta, avgDevRho,
                                     trendAlpha, trendBeta, trendRho,
                                     wells, nWells, nz);
  }

  bool useNewAlgorithm = true;
  if (useNewAlgorithm) 
  {
    std::vector<KrigingData2D> krigingDataAlpha(nz);
    std::vector<KrigingData2D> krigingDataBeta(nz);
    std::vector<KrigingData2D> krigingDataRho(nz);
    
    setupKrigingData2D(krigingDataAlpha,krigingDataBeta,krigingDataRho,
                       trendAlpha,trendBeta,trendRho,
                       modelSettings->getGridOutputFlag(),
                       simbox,wells,nWells);
    
    const CovGrid2D & covGrid2D = makeCovGrid2D(simbox,
                                                modelSettings->getBackgroundVario(),
                                                modelSettings->getDebugFlag());
    
    makeKrigedBackground(krigingDataAlpha, bgAlpha, trendAlpha, simbox, covGrid2D, "Vp");
    makeKrigedBackground(krigingDataBeta, bgBeta, trendBeta, simbox, covGrid2D, "Vs");
    makeKrigedBackground(krigingDataRho, bgRho, trendRho, simbox, covGrid2D, "Rho");
    delete &covGrid2D;
  }
  else
  {
    //
    // NBNB-PAL: Delete this code soon (2009.02.06) ...
    //
    KrigingData3D * krigingData = NULL;
    setupKrigingData3D(krigingData,wells,simbox,trendAlpha,trendBeta,trendRho,nWells);
    
    const int nx = simbox->getnx();
    const int ny = simbox->getny();
    bgAlpha = new FFTGrid(nx, ny, nz, nx, ny, nz);
    bgBeta  = new FFTGrid(nx, ny, nz, nx, ny, nz);
    bgRho   = new FFTGrid(nx, ny, nz, nx, ny, nz);
    fillInVerticalTrend(bgAlpha, trendAlpha);
    fillInVerticalTrend(bgBeta, trendBeta);
    fillInVerticalTrend(bgRho, trendRho);
    
    interpolateBackgroundTrend(krigingData,
                               bgAlpha, 
                               bgBeta, 
                               bgRho,
                               simbox,
                               modelSettings->getBackgroundVario(),
                               modelSettings->getDebugFlag());
  delete krigingData;
  }

  delete [] avgDevAlpha;
  delete [] avgDevBeta;
  delete [] avgDevRho;
  delete [] avgDevVel;

  delete [] trendAlpha;
  delete [] trendBeta;
  delete [] trendRho;
  delete [] trendVel;
}
 
//---------------------------------------------------------------------------
void
Background::calculateVelocityDeviations(FFTGrid   * velocity,                           
                                        WellData ** wells,
                                        Simbox    * simbox,
                                        float    *& trendVel,
                                        float    *& avgDevVel,
                                        float     * avgDevAlpha,
                                        int         outputFlag,
                                        int         nWells)
{
  if((outputFlag & ModelSettings::BACKGROUND_TREND) > 0) {
    std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + "VpFromFile";
    velocity->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
  }
    
  //
  // Calculate deviation between well data and trend
  //
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * velocityLog = new float[maxBlocks];

  const int nz = simbox->getnz();
  float * vtAlpha    = new float[nz];
  float * vtVelocity = new float[nz];

  for (int k=0 ; k<nz ; k++)
    trendVel[k]=0.0;

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    const float * alphaLog = bl->getAlphaHighCutBackground();
    bl->getVerticalTrend(alphaLog, vtAlpha);
    bl->getBlockedGrid(velocity, velocityLog);
    bl->getVerticalTrend(velocityLog, vtVelocity);
    float sumDev = 0.0f;
    int count = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (vtAlpha[k] != RMISSING) {
        trendVel[k] += vtVelocity[k]; 
        float diff = exp(vtAlpha[k]) - vtVelocity[k]; // Velocity trend is in exp-domain
        sumDev += diff*diff;
        count++;
      }
    }
    if (count > 0) 
      sumDev /= count;
    avgDevVel[w] = sqrt(sumDev);
  }
  delete [] vtVelocity;
  delete [] vtAlpha;
  delete [] velocityLog;

  for (int k=0 ; k<nz ; k++)
    trendVel[k] /= nWells;

  LogKit::LogFormatted(LogKit::LOW,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
  LogKit::LogFormatted(LogKit::LOW,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
  LogKit::LogFormatted(LogKit::LOW,"Well             TrendFromFile  TrendFromData\n");
  LogKit::LogFormatted(LogKit::LOW,"---------------------------------------------\n");
  for (int i=0 ; i<nWells ; i++)
    LogKit::LogFormatted(LogKit::LOW,"%-24s %5.1f          %5.1f\n",
                         wells[i]->getWellname(),avgDevVel[i],avgDevAlpha[i]);
}

//---------------------------------------------------------------------------
void
Background::calculateBackgroundTrend(float             * trend,
                                     float             * avgDev,
                                     WellData         ** wells,
                                     Simbox            * simbox,
                                     float               logMin, 
                                     float               logMax,
                                     float               maxHz, 
                                     bool                write1D,
                                     bool                write3D,
                                     int                 nWells, 
                                     bool                hasVelocityTrend,
                                     const std::string & name)
{
  const int   nz = simbox->getnz();
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  calculateVerticalTrend(wells, trend, 
                         logMin, logMax, 
                         maxHz, nWells, 
                         nz, dz, name);
  
  if(write1D == true) {
    writeVerticalTrend(trend, dz, nz, name);
  }
  
  calculateDeviationFromVerticalTrend(wells, trend, avgDev,
                                      nWells, nz, name);
  
  if(write3D == true && !(name=="Vp" && hasVelocityTrend))
  {
    const int nx = simbox->getnx();
    const int ny = simbox->getny();
    FFTGrid * trendGrid = new FFTGrid(nx, ny, nz, nx, ny, nz);
    fillInVerticalTrend(trendGrid, trend);

    std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + name;
    trendGrid->writeFile(fileName, IO::PathToBackground(), simbox, "exptrans");
    //trendGrid->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
    delete trendGrid;
  }
}

//-------------------------------------------------------------------------------
void
Background::setupKrigingData2D(std::vector<KrigingData2D> & krigingDataAlpha,
                               std::vector<KrigingData2D> & krigingDataBeta,
                               std::vector<KrigingData2D> & krigingDataRho,
                               float                      * trendAlpha,
                               float                      * trendBeta, 
                               float                      * trendRho , 
                               int                          outputFlag,
                               Simbox                     * simbox,
                               WellData                  ** wells,
                               const int                    nWells)
{
  //
  // Although unnecessary, we have chosen to set up kriging data fro 
  // Vp, Vs and Rho simultaneously. This gives code easier to read.
  //
  int totBlocks = 0;
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsExtendedBG()->getNumberOfBlocks();
    totBlocks += nBlocks;
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }  
  
  KrigingData3D forLogging(totBlocks); 
  
  const int   nz = simbox->getnz();
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  float * blAlpha = new float[maxBlocks];   // bl = blocked logs
  float * blBeta  = new float[maxBlocks];
  float * blRho   = new float[maxBlocks];

  float * vtAlpha = new float[nz];          // vt = vertical trend
  float * vtBeta  = new float[nz];
  float * vtRho   = new float[nz];

  for (int w = 0 ; w < nWells ; w++)
  {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    const int nBlocks = bl->getNumberOfBlocks();

    Utils::copyVector(bl->getAlphaHighCutBackground(), blAlpha, nBlocks);
    Utils::copyVector(bl->getBetaHighCutBackground(), blBeta, nBlocks);
    Utils::copyVector(bl->getRhoHighCutBackground(), blRho, nBlocks);
    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    bl->getVerticalTrend(blAlpha, vtAlpha);
    bl->getVerticalTrend(blBeta, vtBeta);
    bl->getVerticalTrend(blRho, vtRho);
    //
    // Kriging vertical trend (vt....) against global vertical trend (trend...)
    //
    Kriging1D::krigVector(vtAlpha, trendAlpha, nz, dz);
    Kriging1D::krigVector(vtBeta, trendBeta, nz, dz);
    Kriging1D::krigVector(vtRho, trendRho, nz, dz);
    //
    // Use kriged vertical trend where original log is not defined.
    //
    const int * ipos = bl->getIpos();
    const int * jpos = bl->getJpos();
    const int * kpos = bl->getKpos();

    for (int m = 0 ; m < nBlocks ; m++) 
    {
      int i = ipos[m];
      int j = jpos[m];
      int k = kpos[m];

      if (blAlpha[m] == RMISSING) 
      {
        blAlpha[m] = vtAlpha[k];
      }
      if (blBeta[m] == RMISSING) 
      {
        blBeta[m] = vtBeta[k];
      }
      if (blRho[m] == RMISSING) 
      {
        blRho[m] = vtRho[k];
      }
      krigingDataAlpha[k].addData(i, j, blAlpha[m]);
      krigingDataBeta[k].addData(i, j, blBeta[m]);
      krigingDataRho[k].addData(i, j, blRho[m]);
    }

    forLogging.addData(blAlpha, blBeta, blRho,
                       ipos, jpos, kpos,
                       nBlocks);
  }

  for (int k=0 ; k<nz ; k++) 
  {
    krigingDataAlpha[k].findMeanValues();
    krigingDataBeta[k].findMeanValues();
    krigingDataRho[k].findMeanValues();
  }

  if((outputFlag & ModelSettings::BACKGROUND) > 0) {
    forLogging.divide();
    std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + IO::SuffixGeneralData();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    forLogging.writeToFile(fileName);
  }
  
  delete [] vtAlpha;
  delete [] vtBeta;
  delete [] vtRho;

  delete [] blAlpha;
  delete [] blBeta;
  delete [] blRho;
}

//---------------------------------------------------------------------------
const CovGrid2D &
Background::makeCovGrid2D(Simbox * simbox,
                          Vario  * vario, 
                          int      debugFlag)
{
  //
  // Pretabulate all needed covariances
  //
  const int    nx = simbox->getnx();
  const int    ny = simbox->getny();
  
  const float  dx = static_cast<float>(simbox->getdx());
  const float  dy = static_cast<float>(simbox->getdy());
  
  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);
  
  if(debugFlag == 1) {
    std::string baseName = IO::PrefixBackground() + "covGrid2D" + IO::SuffixAsciiIrapClassic();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    cov->writeToFile(fileName);
  }
  return (*cov); 
}

//---------------------------------------------------------------------------
void
Background::makeKrigedBackground(const std::vector<KrigingData2D> & krigingData,
                                 FFTGrid                         *& bgGrid,
                                 float                            * trend,
                                 Simbox                           * simbox,
                                 const CovGrid2D                  & covGrid2D,
                                 const std::string                & type)
{
  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::LOW,text);

  const int     nx  = simbox->getnx();
  const int     ny  = simbox->getny();
  const int     nz  = simbox->getnz();

  const int    nxp  = nx;
  const int    nyp  = ny;
  const int    nzp  = nz;
  const int    rnxp = 2*(nxp/2 + 1);

  const double  x0  = simbox->getx0();
  const double  y0  = simbox->gety0();
  const double  lx  = simbox->getlx();
  const double  ly  = simbox->getly();

  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);

  float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float nextMonitor = monitorSize;
  std::cout 
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^"; 

  bgGrid = new FFTGrid(nx, ny, nz, nxp, nyp, nzp);  
  bgGrid->createRealGrid();
  bgGrid->setType(FFTGrid::PARAMETER);
  bgGrid->setAccessMode(FFTGrid::WRITE);

  for (int k=0 ; k<nzp ; k++)
  {
    // Set trend for layer
    surface.Assign(trend[k]);                                    

    // Kriging of layer
    Kriging2D::krigSurface(surface, krigingData[k], covGrid2D);  

    // Set layer in background model from surface
    for(int j=0 ; j<nyp ; j++) {        
      for(int i=0 ; i<rnxp ; i++) {
        if(i<nxp)
          bgGrid->setNextReal(float(surface(i,j)));
        else
          bgGrid->setNextReal(0);  //dummy in padding
      }
    }
    
    // Log progress
    if (static_cast<float>(k+1) >= nextMonitor) 
    { 
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  bgGrid->endAccess();
}

//-------------------------------------------------------------------------------
void
Background::setupKrigingData3D(KrigingData3D *& krigingData,
                               WellData      ** wells,
                               Simbox         * simbox,
                               float          * trendAlpha,
                               float          * trendBeta, 
                               float          * trendRho , 
                               const int        nWells)
{
  //
  // NBNB-PAL: Foreløpig setter vi opp krigingData for Vp, Vs og Rho samtidig.
  // Dersom vi lager en 2D kriging, kan disse splittes opp. Pr. idag må vi
  // imidlertid bruke KrigingData3D klassen.
  //
  int totBlocks = 0;
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsExtendedBG()->getNumberOfBlocks();
    totBlocks += nBlocks;
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }  
  
  krigingData = new KrigingData3D(totBlocks); 

  const int   nz = simbox->getnz();
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  float * blAlpha = new float[maxBlocks];   // bl = blocked logs
  float * blBeta  = new float[maxBlocks];
  float * blRho   = new float[maxBlocks];

  float * vtAlpha = new float[nz];          // vt = vertical trend
  float * vtBeta  = new float[nz];
  float * vtRho   = new float[nz];

  for (int w = 0 ; w < nWells ; w++)
  {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    const int nBlocks = bl->getNumberOfBlocks();

    Utils::copyVector(bl->getAlphaHighCutBackground(), blAlpha, nBlocks);
    Utils::copyVector(bl->getBetaHighCutBackground(), blBeta, nBlocks);
    Utils::copyVector(bl->getRhoHighCutBackground(), blRho, nBlocks);
    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    bl->getVerticalTrend(blAlpha, vtAlpha);
    bl->getVerticalTrend(blBeta, vtBeta);
    bl->getVerticalTrend(blRho, vtRho);
    //
    // Kriging vertical trend (vt....) against global vertical trend (trend...)
    //
    Kriging1D::krigVector(vtAlpha, trendAlpha, nz, dz);
    Kriging1D::krigVector(vtBeta, trendBeta, nz, dz);
    Kriging1D::krigVector(vtRho, trendRho, nz, dz);
    //
    // Use kriged vertical trend where original log is not defined.
    //
    const int * ipos = bl->getIpos();
    const int * jpos = bl->getJpos();
    const int * kpos = bl->getKpos();

    for (int m = 0 ; m < nBlocks ; m++) 
    {
      if (blAlpha[m] == RMISSING) 
      {
        blAlpha[m] = vtAlpha[kpos[m]];
      }
      if (blBeta[m] == RMISSING) 
      {
        blBeta[m] = vtBeta[kpos[m]];
      }
      if (blRho[m] == RMISSING) 
      {
        blRho[m] = vtRho[kpos[m]];
      }
    }

    //
    // Add to kriging data object
    //
    krigingData->addData(blAlpha, blBeta, blRho,
                         ipos, jpos, kpos,
                         nBlocks);
  }
  krigingData->divide();

  std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + IO::SuffixGeneralData();
  std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
  krigingData->writeToFile(fileName);

  delete [] vtAlpha;
  delete [] vtBeta;
  delete [] vtRho;

  delete [] blAlpha;
  delete [] blBeta;
  delete [] blRho;
}

//---------------------------------------------------------------------------
void
Background::interpolateBackgroundTrend(KrigingData3D * krigingData,
                                       FFTGrid       * bgAlpha,
                                       FFTGrid       * bgBeta,
                                       FFTGrid       * bgRho,
                                       Simbox        * simbox,
                                       Vario         * vario,
                                       int             debugFlag)
{
  GenExpVario* pVario = dynamic_cast<GenExpVario*>(vario);
  float power;
  if (pVario)
    power = pVario->getPower();
  else {
    LogKit::LogFormatted(LogKit::LOW,"ERROR: Only the General Exponential Variogram type is supported when creating background model.\n");
    exit(1);
  }
  float rangeX, rangeY;
  float rotAngle;
  vario->getParams(rangeX, rangeY, rotAngle);
  float rangeZ = 0.00001f;  // Simulate 2D kriging using 3D algorithm

  const int nx = simbox->getnx();
  const int ny = simbox->getny();
  const int nz = simbox->getnz();
  
  const float dx = static_cast<float>(simbox->getdx());
  const float dy = static_cast<float>(simbox->getdy());
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  CovGridSeparated covAlpha(nx, ny, nz, dx, dy, dz, rangeX, rangeY, rangeZ, power, rotAngle, false);
  CovGridSeparated covBeta (nx, ny, nz, dx, dy, dz, rangeX, rangeY, rangeZ, power, rotAngle, false);
  CovGridSeparated covRho  (nx, ny, nz, dx, dy, dz, rangeX, rangeY, rangeZ, power, rotAngle, false);
  CovGridSeparated covCrAlphaBeta(nx, ny, nz);
  CovGridSeparated covCrAlphaRho(nx, ny, nz);
  CovGridSeparated covCrBetaRho(nx, ny, nz);

  if(debugFlag == 1)
  {
    covAlpha.writeXYGrid("covAlpha");
    covBeta.writeXYGrid("covBeta");
    covRho.writeXYGrid("covRho");
  }

  CKrigingAdmin kriging(*simbox, 
                        krigingData->getData(),
                        krigingData->getNumberOfData(),
                        covAlpha, covBeta, covRho, 
                        covCrAlphaBeta, covCrAlphaRho, covCrBetaRho, 
                       DataTarget_, true);

  LogKit::LogFormatted(LogKit::LOW,"\nFill volumes using kriging:");
  kriging.KrigAll(*bgAlpha, *bgBeta, *bgRho, false);
}

//-------------------------------------------------------------------------------
void 
Background::calculateVerticalTrend(WellData         ** wells,
                                   float             * trend, 
                                   float               logMin,
                                   float               logMax,
                                   float               maxHz,
                                   int                 nWells,
                                   int                 nz,
                                   float               dz,
                                   const std::string & name)
{  
  float * filtered_log = new float[nz]; 
  float * wellTrend    = filtered_log;   // Use same memory twice
  int   * count = new int[nz];
  //
  // Calculate the average values of well log
  // ----------------------------------------
  // When calculating the vertical trend, we do not want each well
  // to contribute with more than one value to each layer 'k',
  // and therefore we calculate average values for each well 
  // first. This way we avoid strange behaviour caused by
  // deviated/horizontal wells.
  //
  for (int k = 0 ; k < nz ; k++) {
    trend[k] = 0.0f;
    count[k] = 0;
  }
  int iWells = 0;
  for (int w = 0 ; w < nWells ; w++) { 
    if (wells[w]->getUseForBackgroundTrend()) {
      BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
      if (name == "Vp")
        bl->getVerticalTrend(bl->getAlpha(), wellTrend);
      else if (name == "Vs")
        bl->getVerticalTrend(bl->getBeta(), wellTrend);
      else if (name == "Rho")
        bl->getVerticalTrend(bl->getRho(), wellTrend);
      else {
        LogKit::LogFormatted(LogKit::LOW,"ERROR in Background::calculateVerticalTrend(): ");
        LogKit::LogFormatted(LogKit::LOW,"Log \'"+name+"\' requested, but no such log exists.\n");
        exit(1);
      }
      for (int k = 0 ; k < nz ; k++) {
        if (wellTrend[k] != RMISSING) {
          trend[k] += exp(wellTrend[k]);
          count[k]++;
        }
      }
      iWells++;
    }
  }
  if (iWells > 0) {
    for (int k = 0 ; k < nz ; k++) {
      if (count[k] > 0) {
        trend[k] = trend[k]/count[k]; 
      }
    }
  }
  else {
    LogKit::LogFormatted(LogKit::LOW,"\nERROR in Background::calculateVerticalTrend(): There are no wells\n");
    LogKit::LogFormatted(LogKit::LOW,"available for the estimation of background trend.\n");
    exit(1);
  }

  bool use_local_linear_regression = true;
  if (use_local_linear_regression) {
    smoothTrendWithLocalLinearRegression(trend, count,
                                         iWells, nz, dz, 
                                         logMin,
                                         logMax,
                                         name);
    WellData::applyFilter(filtered_log, 
                          trend, 
                          nz, 
                          dz, 
                          maxHz);
  }
  else {
    float * interpolated_log = new float[nz]; 
    smoothTrendWithMovingAverage(trend, 
                                 count, 
                                 iWells, 
                                 nz);
    WellData::interpolateLog(interpolated_log, trend, nz); 
    WellData::applyFilter(filtered_log, interpolated_log, nz, dz, maxHz);
    extrapolateTrend(name,  filtered_log, nz); 
    delete [] interpolated_log;
  }
  
  for (int i=0 ; i<nz ; i++) {
    trend[i] = filtered_log[i];
  }
  
  delete [] filtered_log; 
  delete [] count;
}


//-------------------------------------------------------------------------------
void
Background::smoothTrendWithLocalLinearRegression(float      * trend, 
                                                 int        * count,
                                                 int          iWells,
                                                 int          nz,
                                                 float        dz,
                                                 float        min_value, 
                                                 float        max_value,
                                                 std::string  parName) 
{
  bool debug = false;
  //
  // 1. Center-parts of scatter plots
  //
  // In the center parts of the scatter plots the average value should be
  // accepted as a trend value if the number of data points behind each 
  // average is fraction * iWells, where fraction is the acceptance
  // fraction, possibly larger than one.
  //
  // Sometimes we have only one, two, or three wells available. In the
  // case of two and three wells, the logs for these may differ considerably,
  // in which case the trend need to be stabilised by requiring a minimum
  // number of data points behind each trend value. The denser the log is 
  // sampled, the more prone it is to numerical instabilities. This 
  // minimum value is therefore linked to the sampling density 'dz'.
  //
  // Finally, we must require a definite minimum number of data points 
  // that should be behind each trend value.
  // 
  // nDataMin = max(min_req_points, max(min_time_sample, fraction * iWells))
  //
  // 2. End-parts of scatter plot
  //
  // If the blocked wells do not contain any values for the lower layers of 
  // the simbox, the end part of the scatter plot needs special attention.
  //
  // We must possibly require a larger min_req_points to avoid an 
  // "arbitrary" end behaviour.
  //

  float fraction   = 3.0f;                      // Require minimum 3*iWells
  int   nTimeLimit = static_cast<int>(50.0/dz); // The smaller sampling density, the more values are needed.
  int   nLowLimit  = 10;                        // Require minimum 10
  int   nDataMin   = std::max(nLowLimit, std::max(nTimeLimit, int(fraction * iWells)));

  bool  use_weights = true;
  bool  error = false;

  //
  // Copy the average values (stored in array 'trend') to the array 'mean'.
  //
  float * mean = new float[nz];
  for (int k = 0 ; k < nz ; k++) {
    mean[k] = trend[k];
  }

  float * x = new float[nz];  // Time indices
  float * y = new float[nz];  // Log values
  float * w = new float[nz];  // Weights (number of data behind each avg.)

  for (int k = 0 ; k < nz ; k++) {
    int n = 0;
    int nData = 0;
    if (debug) 
      LogKit::LogFormatted(LogKit::LOW,"k=%d\n",k);
    //
    // 1. Add current data point to local data set if present.
    //
    if (count[k] > 0) {
      w[0]   = static_cast<float>(count[k]);
      x[0]   = static_cast<float>(k);
      y[0]   = trend[k];
      nData += count[k];
      if (debug) 
        LogKit::LogFormatted(LogKit::LOW,"   A:t=%.2f   x[0] y[0]  %d   %.2f\n",dz*(x[0] + 0.5f),int(x[0]),y[0]);
      n++;
    }

    //
    // 2. Add local data points to get 'nDataMin' points behind each trend 
    //    value. Note that the bandwidth varies
    //
    int i = 0;
    while (nData < nDataMin) {
      i++;
      if (k - i >= 0 && count[k - i] > 0) {
        w[n]   = static_cast<float>(count[k - i]);
        x[n]   = static_cast<float>(k - i);
        y[n]   = mean [k - i];
        nData += count[k - i];
        if (debug) 
          LogKit::LogFormatted(LogKit::LOW,"   B:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k + i < nz  && count[k + i] > 0) {
        w[n]   = static_cast<float>(count[k + i]);
        x[n]   = static_cast<float>(k + i);
        y[n]   = mean [k + i];
        nData += count[k + i];
        if (debug) 
          LogKit::LogFormatted(LogKit::LOW,"   C:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k-i < 0 && k+i >= nz) { // We will never find enough data
        break;
      }
    }

    //
    // Calculate normalised weights
    //
    if (use_weights) 
      for (i = 0 ; i < n ; i++)
        w[i] /= nData;
    else
      for (i = 0 ; i < n ; i++) 
        w[i] = 1.0f/n;
    
    //
    // We need at least two points to make a line.
    //
    if (n > 1) {
      //
      // Estimate local regression line: y = bx + a
      //
      float Sx  = x[0]*w[0];
      float Sy  = y[0]*w[0];
      float Sxx = x[0]*w[0]*x[0];          
      float Sxy = x[0]*w[0]*y[0];
      for (i = 1 ; i < n ; i++) {
        Sx  += x[i]*w[i];
        Sy  += y[i]*w[i];
        Sxx += x[i]*w[i]*x[i];
        Sxy += x[i]*w[i]*y[i];
      }
      float b = (Sxy - Sx*Sy)/(Sxx - Sx*Sx);
      float a = (Sy - b*Sx);
      if (debug) 
        LogKit::LogFormatted(LogKit::LOW,"Sx, Sy, Sxx, Sxy  : %.4f, %.4f, %.4f, %.4f     a, b : %.4f %.4f\n",Sx,Sy,Sxx,Sxy,a,b);
      
      //
      // Estimate value of regression line at requested point.
      //
      float value = a + b*static_cast<float>(k);
      if (value < min_value || value > max_value) {
        if (debug)
          LogKit::LogFormatted(LogKit::LOW,"   TREND: trend[k] = %.2f\n",value);
        error = true;
        break;
      }
      trend[k] = log(value);
    }
    else {
      trend[k] = log(y[0]);
    }
    if (debug) 
      LogKit::LogFormatted(LogKit::LOW,"   TREND: trend[k] = %.2f        (minLog/maxLog = %.2f / %.2f)\n",exp(trend[k]),min_value,max_value);
  }

  if (error) {
    //
    // NBNB-PAL: Here we should possibly first try a global linear regression...
    //
    LogKit::LogFormatted(LogKit::LOW,"\nWARNING : The calculation of the vertical trend for parameter "+parName+" using local linear\n");
    LogKit::LogFormatted(LogKit::LOW,"          regression failed - trying global mean instead. Possible causes: \n");
    LogKit::LogFormatted(LogKit::LOW,"          1) Available logs cover too small a part of inversion grid giving extrapolation problems.\n");
    LogKit::LogFormatted(LogKit::LOW,"          2) There are too many layers in grid compared to well logs available.\n");
    float sum = 0.0f;
    int nData = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (count[k] > 0) {
        if (use_weights) {
          sum   += mean[k]*count[k];
          nData += count[k];
          if (debug)           
            LogKit::LogFormatted(LogKit::LOW,"k=%d  count[k], mean[k]  nData, sum  %d  %8.3f     %d  %8.3f\n",
                             k,count[k],mean[k],nData,sum);
        }
        else {
          sum   += mean[k];
          nData += 1;
        }
      }
    }
    float global_mean = sum/nData;
    for (int k = 0 ; k < nz ; k++) {
      trend[k] = log(global_mean);
      if (debug) 
        LogKit::LogFormatted(LogKit::LOW,"   TREND: k = %d   trend[k] = %.2f\n",k,exp(trend[k]));
    }
    LogKit::LogFormatted(LogKit::LOW,"\nGlobal mean for parameter %s = %.2f\n\n",parName.c_str(),global_mean);
  }
  delete [] x;
  delete [] y;
  delete [] w;
  delete [] mean;
}
 
//-------------------------------------------------------------------------------
void
Background::smoothTrendWithMovingAverage(float * trend, 
                                         int   * count,
                                         int     iWells,
                                         int     nz) 
{
  //
  // It is unreasonable to require nDataMin = iWells, and we use
  // an acceptance fraction to adjust the minimum required number 
  // of data in each trend value.
  //
  // Useful range :  0.5 < fraction < 1.0
  //
  float fraction = 0.51f;                          
  int   nDataMin = std::max(1, int(fraction * iWells + 0.5f)); 

  float * mean = new float[nz];
  for (int k = 0 ; k < nz ; k++) {
    mean[k] = trend[k];
  }

  for (int k = 0 ; k < nz ; k++) {
    if (count[k] == 0) {
      //
      // Do NOT estimate a value here. Fix later with interpolation 
      // or extrapolation. If the lower part of the simbox contains no 
      // BWs, making these values from layers above will affect the  
      // Fourier filtering. This should be avoided since these data 
      // are made from other data points.
      //
      trend[k] = RMISSING;
    }
    else {
      //
      // We require nDataMin data points behind each trend value. 
      //
      int   nData = count[k];
      float value = mean[k];
      LogKit::LogFormatted(LogKit::LOW,"k=%d      mean[k] count[k]  %.2f %d\n",k,mean[k],count[k]);
      
      int nelms = 1;
      if (count[k] == 0)
        nelms = 0;
      
      int i = 0;
      while (nData < nDataMin) {
        i++;
        if (k - i >= 0 && count[k - i] > 0) {
          nData += count[k - i];
          value += mean [k - i];
          nelms++;
          LogKit::LogFormatted(LogKit::LOW,"   k-i = %d  mean[k-i] count[k-i]  %.2f %d\n",k-i,mean[k-i],count[k-i]);
        }
        if (k + i < nz  && count[k + i] > 0) {
          nData += count[k + i];
          value += mean [k + i];
          nelms++;
          LogKit::LogFormatted(LogKit::LOW,"   k+i = %d  mean[k+i] count[k+i]  %.2f %d\n",k+i,mean[k+i],count[k+i]);
        }
      }
      trend[k] = log(value/nelms);
      LogKit::LogFormatted(LogKit::LOW,"   trend[k]  %.2f\n",exp(trend[k])); 
    } // end if
  }  
  delete [] mean;
}

//-------------------------------------------------------------------------------
void
Background::extrapolateTrend(std::string  pName, 
                             float      * log,
                             int          nz) 
{
  //
  // Extrapolate log[] in both ends if needed
  //  
  int i=0;
  while (i<nz && log[i]==RMISSING)
    i++;
  if (i < nz - 1) 
  { 
    int first_nonmissing = i;
    i = nz - 1;
    while (i>0 && log[i]==RMISSING)
      i--;
    int last_nonmissing = i;
    
    for(int i=0 ; i < first_nonmissing ; i++) { 
      log[i] = log[first_nonmissing]; 
    }
    for(int i=last_nonmissing + 1 ; i < nz ; i++) { 
      log[i] = log[last_nonmissing]; 
    }
    if (first_nonmissing > 0)
      LogKit::LogFormatted(LogKit::LOW,"Vertical trend for %s extrapolated first %d cells.\n",
                           pName.c_str(), first_nonmissing + 1);
    if (nz - 1 - last_nonmissing > 0)
      LogKit::LogFormatted(LogKit::LOW,"Vertical trend for %s extrapolated last %d cells.\n",
                           pName.c_str(), nz - 1 - last_nonmissing);
    //for (int i=0 ; i<nz ; i++) {
    //  printf("i log[i]  %d %.3f\n",i,log[i]);
    //}
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"WARNING: All ... trend values are missing.\n");
  }
}

//-------------------------------------------------------------------------------
void 
Background::writeVerticalTrend(float      * trend, 
                               float        dz,
                               int          nz,
                               std::string  name) 
{  
  float z0 = dz/2.0f;
  std::string baseName = IO::PrefixBackground() + IO::PrefixTrend() + name + IO::SuffixAsciiIrapClassic();
  std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  for (int i=0 ; i<nz ; i++) {
    file << std::fixed 
         << std::setprecision(2)
         << std::setw(8) << (z0 + i*dz)     << " "
         << std::setprecision(3)
         << std::setw(8) << exp( trend[i] ) << " "
         << "0.00\n";
  }
  file << "999.00 999.00 999.00" << std::endl;
  file.close();
}

//-------------------------------------------------------------------------------
void
Background::calculateDeviationFromVerticalTrend(WellData    ** wells,
                                                const float  * globalTrend, 
                                                float        * avg_dev,
                                                int            nWells,
                                                int            nz,
                                                std::string    name)
{
  float * wellTrend = new float[nz];

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    if (name == "Vp")
      bl->getVerticalTrend(bl->getAlphaHighCutBackground(), wellTrend);
    else if (name == "Vs")
      bl->getVerticalTrend(bl->getBetaHighCutBackground(), wellTrend);
    else if (name == "Rho")
      bl->getVerticalTrend(bl->getRhoHighCutBackground(), wellTrend);
    else {
      LogKit::LogFormatted(LogKit::LOW,"ERROR in Background::calculateVerticalTrend(): ");
      LogKit::LogFormatted(LogKit::LOW,"Log \'"+name+"\' requested, but no such log exists.\n");
      exit(1);
    }
    float sum_dev = 0.0f;
    int count = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (wellTrend[k] != RMISSING) {
        float diff = exp(wellTrend[k]) - exp(globalTrend[k]);
        sum_dev += diff*diff;
        count++;
      }
    }
    if (count > 0) 
      sum_dev /= count;
    avg_dev[w] = sqrt(sum_dev);
  }
  delete [] wellTrend;
}

//-------------------------------------------------------------------------------
void
Background::writeDeviationsFromVerticalTrend(const float *  avg_dev_alpha,
                                             const float *  avg_dev_beta,
                                             const float *  avg_dev_rho,
                                             const float *  trend_alpha,
                                             const float *  trend_beta,
                                             const float *  trend_rho,
                                             WellData    ** wells,
                                             const int      nWells,
                                             const int      nz)
{
  float global_mean_alpha = 0.0f;
  float global_mean_beta  = 0.0f;
  float global_mean_rho   = 0.0f;

  for (int k=0 ; k<nz ; k++)
  {
    global_mean_alpha += exp(trend_alpha[k]);
    global_mean_beta  += exp(trend_beta[k]);
    global_mean_rho   += exp(trend_rho[k]);
  }
  global_mean_alpha /= nz;
  global_mean_beta  /= nz;
  global_mean_rho   /= nz;

  //
  // Find the relative average deviations (mean of Vp,Vs and Rho deviations).
  //
  float * rel_avg_dev = new float[nWells];
  for (int i=0 ; i<nWells ; i++)
  {
    float rel_dev_alpha = avg_dev_alpha[i]/global_mean_alpha;
    float rel_dev_beta  = avg_dev_beta[i]/global_mean_beta;
    float rel_dev_rho   = avg_dev_rho[i]/global_mean_rho;
    rel_avg_dev[i] = (rel_dev_alpha + rel_dev_beta + rel_dev_rho)/3;
  }
  //
  // Sort deviations to find worst well.
  //
  int * index = new int[nWells];
  float cur_min = 99999.0f;
  int   cur_pos = 0;
  for (int i=0 ; i<nWells ; i++)
  {
    float cur_max = 0.0f;
    for (int j=0 ; j<nWells ; j++)
    {
      if (rel_avg_dev[j] > cur_max && rel_avg_dev[j] < cur_min) 
      {   
        cur_max = rel_avg_dev[j];
        cur_pos = j;
      }
    }
    index[i] = cur_pos;
    cur_min = cur_max;
  }
  //
  // Print results
  //
  if (nWells > 0) 
  {
    LogKit::LogFormatted(LogKit::LOW,"\nSummary of average deviation from vertical trend (well with largest misfit listed first):\n\n");
    LogKit::LogFormatted(LogKit::LOW,"Well                        Vp       Vs      Rho\n");
    LogKit::LogFormatted(LogKit::LOW,"------------------------------------------------\n");
  }  
  for (int i=0 ; i<nWells ; i++)
  {
    int ii = index[i];
    LogKit::LogFormatted(LogKit::LOW,"%-24s %5.1f    %5.1f    %5.3f\n", wells[ii]->getWellname(),
                     avg_dev_alpha[ii], avg_dev_beta[ii], avg_dev_rho[ii]);
  }

  if (nWells == 1)
  {
    LogKit::LogFormatted(LogKit::HIGH,"\nNOTE: A deviation may be observed even with one well since the global trend is");
    LogKit::LogFormatted(LogKit::HIGH,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
  }
  delete [] rel_avg_dev;
  delete [] index;
}

//-------------------------------------------------------------------------------
void           
Background::fillInVerticalTrend(FFTGrid     * grid, 
                                const float * trend) 
{
  const int nzp = grid->getNzp();  // equals nx
  const int nyp = grid->getNyp();  // equals ny
  const int nxp = grid->getNxp();  // equals nz

  grid->createRealGrid();
  grid->setType(FFTGrid::PARAMETER);
  grid->setAccessMode(FFTGrid::WRITE);

  int rnxp = 2*(nxp/2 + 1);
  for (int k = 0; k < nzp; k++) 
    for (int j = 0; j < nyp; j++)
      for (int i = 0; i < rnxp; i++)
        grid->setNextReal(trend[k]);

  grid->endAccess();
}

//-------------------------------------------------------------------------------
void
Background::findMeanVsVp(FFTGrid * Vp,
                         FFTGrid * Vs)
{
  Vp->setAccessMode(FFTGrid::READ);
  Vs->setAccessMode(FFTGrid::READ);
  double mean = 0;
  int nxp = 2*(Vp->getNxp()/2+1);
  int nyp = Vp->getNyp();
  int nzp = Vp->getNzp();
  int nx  = Vp->getNx();
  int ny  = Vp->getNy();
  int nz  = Vp->getNz();
  for(int k=0;k<nzp;k++) 
    for(int j=0;j<nyp;j++)
      for(int i=0;i<nxp;i++) {
        float v1 = Vp->getNextReal();
        float v2 = Vs->getNextReal();
        if(i < nx && j < ny && k < nz)
          mean += exp(v2-v1);
      }
  mean = mean/double(nx*ny*nz);
  
  Vp->endAccess();
  Vs->endAccess();
  
  vsvp_ = mean;
}
        
//-------------------------------------------------------------------------------
void
Background::setClassicVsVp()
{
  int nxp = backModel_[0]->getNxp();
  int nyp = backModel_[0]->getNyp();
  int nzp = backModel_[0]->getNzp();
  float vp = backModel_[0]->getFirstComplexValue().re;
  vp = float(exp(vp/sqrt(float(nxp*nyp*nzp))));
  float vs = backModel_[1]->getFirstComplexValue().re;
  vs = float(exp(vs/sqrt(float(nxp*nyp*nzp))));
  vsvp_ = vs/vp;
}

//-------------------------------------------------------------------------------
void
Background::resampleBackgroundModel(FFTGrid      *& bgAlpha,
                                    FFTGrid      *& bgBeta,
                                    FFTGrid      *& bgRho,
                                    Simbox        * timeBGSimbox,
                                    Simbox        * timeSimbox,
                                    ModelSettings * modelSettings)
{
  if((modelSettings->getGridOutputFlag() & ModelSettings::EXTRA_GRIDS) > 0) {
    std::string fileName1 = IO::PrefixBackground() + "Vp_BackgroundGrid";
    std::string fileName2 = IO::PrefixBackground() + "Vs_BackgroundGrid";
    std::string fileName3 = IO::PrefixBackground() + "Rho_BackgroundGrid";
    bgAlpha->writeFile(fileName1, IO::PathToBackground(), timeBGSimbox, "exptrans");
    bgBeta->writeFile(fileName2, IO::PathToBackground(), timeBGSimbox, "exptrans");
    bgRho->writeFile(fileName3, IO::PathToBackground(), timeBGSimbox, "exptrans");
  }

  FFTGrid * resBgAlpha = NULL;
  FFTGrid * resBgBeta = NULL;
  FFTGrid * resBgRho = NULL;
  
  LogKit::LogFormatted(LogKit::LOW,"\nResampling background model...\n");
  resampleParameter(resBgAlpha,bgAlpha,timeSimbox, timeBGSimbox);
  resampleParameter(resBgBeta ,bgBeta ,timeSimbox, timeBGSimbox);
  resampleParameter(resBgRho  ,bgRho  ,timeSimbox, timeBGSimbox);
  
  if((modelSettings->getGridOutputFlag() & ModelSettings::EXTRA_GRIDS) > 0) {
    std::string fileName1 = IO::PrefixBackground() + "Vp_InversionGrid";
    std::string fileName2 = IO::PrefixBackground() + "Vs_InversionGrid";
    std::string fileName3 = IO::PrefixBackground() + "Rho_InversionGrid";
    resBgAlpha->writeFile(fileName1, IO::PathToBackground(), timeSimbox, "exptrans");
    resBgBeta->writeFile(fileName2, IO::PathToBackground(), timeSimbox, "exptrans");
    resBgRho->writeFile(fileName3, IO::PathToBackground(), timeSimbox, "exptrans");
  }
  
  delete bgAlpha;
  delete bgBeta;
  delete bgRho;
  
  bgAlpha = resBgAlpha;
  bgBeta  = resBgBeta;
  bgRho   = resBgRho;
 }

//-------------------------------------------------------------------------------
void
Background::resampleParameter(FFTGrid *& pNew,        // Resample to 
                              FFTGrid  * pOld,        // Resample from 
                              Simbox   * simboxNew,
                              Simbox   * simboxOld)
{
  int nx  = simboxNew->getnx();
  int ny  = simboxNew->getny();
  int nz  = simboxNew->getnz();
  //
  // Use same padding as for nonresampled cubes
  //
  int nxp = nx + (pOld->getNxp() - pOld->getNxp()); 
  int nyp = ny + (pOld->getNyp() - pOld->getNyp());
  int nzp = nz + (pOld->getNzp() - pOld->getNzp());

  //
  // Set up relation between old layer index and new layer index using
  //
  // k2 = dz1/dz2 * k1 + (z02 - z01)/dz2    (from dz2*k2 + z02 = dz1*k1 + z01) 
  //
  double * a = new double[nx*ny];
  double * b = new double[nx*ny];

  int ij = 0;
  for(int j=0;j<ny;j++) {
    for(int i=0;i<nx;i++) {
      double dzNew = simboxNew->getdz(i,j); 
      double dzOld = simboxOld->getdz(i,j); 
      double z0New = simboxNew->getTop(i,j); 
      double z0Old = simboxOld->getTop(i,j); 
        a[ij] = dzNew/dzOld;
      b[ij] = (z0New - z0Old)/dzOld; 
      ij++;
    }
  }

  //
  // Resample parameter
  //
  pNew = new FFTGrid(nx,ny,nz,nxp,nyp,nzp);
  pNew->createRealGrid();
  pNew->setType(FFTGrid::PARAMETER);
  pNew->setAccessMode(FFTGrid::WRITE);

  pOld->setAccessMode(FFTGrid::RANDOMACCESS);

  int rnxp = 2*(nxp/2 + 1);

  for(int k=0 ; k<nzp ; k++) {
    int ij=0;
    for(int j=0 ; j<nyp ; j++) {
      for(int i=0 ; i<rnxp ; i++) {
        float value;
        if(i < nx && j < ny && k < nz) {
          int kOld = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]); 
          ij++;
          value = pOld->getRealValue(i, j, kOld);
        }
        else {
          value = RMISSING;
        }
        pNew->setNextReal(value);
      }
    }
  }
  pOld->endAccess();
  pNew->endAccess();

  delete [] a;
  delete [] b;
}

//-------------------------------------------------------------------------------
void
Background::padAndSetBackgroundModel(FFTGrid * bgAlpha,
                                     FFTGrid * bgBeta,
                                     FFTGrid * bgRho)
{
  LogKit::LogFormatted(LogKit::LOW,"\nPadding background model...\n");
  backModel_[0]->fillInFromRealFFTGrid(*bgAlpha);
  backModel_[1]->fillInFromRealFFTGrid(*bgBeta);
  backModel_[2]->fillInFromRealFFTGrid(*bgRho);
}

//-------------------------------------------------------------------------------
void
Background::writeBackgrounds(Simbox      * simbox, 
                             GridMapping * depthMapping, 
                             GridMapping * timeMapping) const 
{
  if(depthMapping != NULL && depthMapping->getSimbox() == NULL) {
    const Simbox * timeSimbox = simbox;
    if(timeMapping != NULL)
      timeSimbox = timeMapping->getSimbox();
    backModel_[0]->setAccessMode(FFTGrid::RANDOMACCESS);
    depthMapping->setMappingFromVelocity(backModel_[0], timeSimbox);
    backModel_[0]->endAccess();
  }

  std::string fileName1 = IO::PrefixBackground() + "Vp" ;
  std::string fileName2 = IO::PrefixBackground() + "Vs" ;
  std::string fileName3 = IO::PrefixBackground() + "Rho";
  backModel_[0]->writeFile(fileName1, IO::PathToBackground(), simbox, "exptrans", 0, depthMapping, timeMapping);
  backModel_[1]->writeFile(fileName2, IO::PathToBackground(), simbox, "exptrans", 0, depthMapping, timeMapping);
  backModel_[2]->writeFile(fileName3, IO::PathToBackground(), simbox, "exptrans", 0, depthMapping, timeMapping);
  //
  // For debugging: write cubes not in ASCII, with padding, and with flat top.
  //
  //backModel_[0]->writeStormFile(fileName1, IO::PathToBackground(), simbox, true, false, true, true);
  //backModel_[1]->writeStormFile(fileName2, IO::PathToBackground(), simbox, true, false, true, true);
  //backModel_[2]->writeStormFile(fileName3, IO::PathToBackground(), simbox, true, false, true, true);
}
