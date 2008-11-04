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
#include "src/krigingdata3d.h"
#include "src/covgridseparated.h"
#include "src/krigingadmin.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"

Background::Background(FFTGrid       ** grids,
                       WellData      ** wells,
                       Simbox         * simbox,
                       ModelSettings  * modelSettings)
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];

  DataTarget_ = 250;

  generateBackgroundModel(wells,
                          simbox, 
                          modelSettings);
  findMeanVsVp();
}

//-------------------------------------------------------------------------------
Background::Background(FFTGrid ** grids) 
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];
  findMeanVsVp();
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
Background::writeBackgrounds(Simbox * simbox) const 
{
  printf("\nExp-transforming cubes...\n");
  backModel_[0]->expTransf();
  backModel_[1]->expTransf();
  backModel_[2]->expTransf();
  backModel_[0]->writeFile("BG_Vp",  simbox);
  backModel_[1]->writeFile("BG_Vs",  simbox);
  backModel_[2]->writeFile("BG_Rho", simbox);
  printf("\nBacktransforming cubes...\n");
  backModel_[0]->logTransf();
  backModel_[1]->logTransf();
  backModel_[2]->logTransf();
}

//-------------------------------------------------------------------------------
void
Background::generateBackgroundModel(WellData      ** wells,
                                    Simbox         * simbox,
                                    ModelSettings  * modelSettings)
{
  const int   nz = simbox->getnz();
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());
  //const double dz    = simbox->getdz();
  const int nWells = modelSettings->getNumberOfWells();

  float * trendAlpha   = new float[nz];
  float * trendBeta    = new float[nz];
  float * trendRho     = new float[nz];

  calculateVerticalTrend(wells, trendAlpha, 
                         modelSettings->getAlphaMin(), 
                         modelSettings->getAlphaMax(),
                         modelSettings->getMaxHzBackground(), 
                         nWells, nz, dz, "Vp");
  calculateVerticalTrend(wells, trendBeta, 
                         modelSettings->getBetaMin(), 
                         modelSettings->getBetaMax(),
                         modelSettings->getMaxHzBackground(), 
                         nWells, nz, dz, "Vs");
  calculateVerticalTrend(wells, trendRho, 
                         modelSettings->getRhoMin(), 
                         modelSettings->getRhoMax(),
                         modelSettings->getMaxHzBackground(), 
                         nWells, nz, dz, "Rho");

  if((modelSettings->getOutputFlag() & ModelSettings::BACKGROUND) > 0) {
    writeVerticalTrend(trendAlpha, dz, nz, "Vp");
    writeVerticalTrend(trendBeta, dz, nz, "Vs");
    writeVerticalTrend(trendRho, dz, nz, "Rho");
  }

  float * avgDevAlpha = new float[nWells];
  float * avgDevBeta  = new float[nWells];
  float * avgDevRho   = new float[nWells];

  calculateDeviationFromVerticalTrend(wells, 
                                      trendAlpha,
                                      avgDevAlpha,
                                      nWells,
                                      nz,
                                      "Vp");
  calculateDeviationFromVerticalTrend(wells, 
                                      trendBeta,
                                      avgDevBeta,
                                      nWells,
                                      nz,
                                      "Vs");
  calculateDeviationFromVerticalTrend(wells, 
                                      trendRho,
                                      avgDevRho,
                                      nWells,
                                      nz,
                                      "Rho");
  
  writeDeviationsFromVerticalTrend(wells, 
                                   avgDevAlpha, avgDevBeta, avgDevRho,
                                   trendAlpha, trendBeta, trendRho,
                                   nWells, nz);
  delete [] avgDevAlpha;
  delete [] avgDevBeta;
  delete [] avgDevRho;

  int totBlocks = 0;
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsPropThick()->getNumberOfBlocks();
    totBlocks += nBlocks;
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  KrigingData3D kd(totBlocks); 

  float * blAlpha = new float[maxBlocks];   // bl = blocked logs
  float * blBeta  = new float[maxBlocks];
  float * blRho   = new float[maxBlocks];

  float * vtAlpha = new float[nz];          // vt = vertical trend
  float * vtBeta  = new float[nz];
  float * vtRho   = new float[nz];

  for (int w = 0 ; w < nWells ; w++)
  {
    BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
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
    kd.addData(blAlpha, blBeta, blRho,
               ipos, jpos, kpos,
               nBlocks);
  }
  delete [] vtAlpha;
  delete [] vtBeta;
  delete [] vtRho;

  delete [] blAlpha;
  delete [] blBeta;
  delete [] blRho;

  kd.divide();
  kd.writeToFile("BG");


  FFTGrid *trendCubeAlpha, *trendCubeBeta, *trendCubeRho;
  createTrendCube(simbox, trendCubeAlpha, trendAlpha);
  createTrendCube(simbox, trendCubeBeta, trendBeta);
  createTrendCube(simbox, trendCubeRho, trendRho);

  delete [] trendAlpha;
  delete [] trendBeta;
  delete [] trendRho;

  if(modelSettings->getDebugFlag() == 1)
  {
    printf("\nExp-transforming cubes...\n");
    trendCubeAlpha->expTransf();
    trendCubeBeta->expTransf();
    trendCubeRho->expTransf();
    trendCubeAlpha->writeStormFile("BG_trendVp", simbox);
    trendCubeBeta->writeStormFile("BG_trendVs", simbox);
    trendCubeRho->writeStormFile("BG_trendRho", simbox);
    printf("\nBacktransforming cubes...\n");
    trendCubeAlpha->logTransf();
    trendCubeBeta->logTransf();
    trendCubeRho->logTransf();
  }

  float power;
  Vario      * vario  = modelSettings->getBackgroundVario();
  GenExpVario* pVario = dynamic_cast<GenExpVario*>(vario);
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

  const int   nx = simbox->getnx();
  const int   ny = simbox->getny();
  const float dx = static_cast<float>(simbox->getdx());
  const float dy = static_cast<float>(simbox->getdy());

  CovGridSeparated covAlpha(nx, ny, nz, dx, dy, dz, rangeX, rangeY, rangeZ, power, rotAngle,false);
  CovGridSeparated covBeta (nx, ny, nz, dx, dy, dz, rangeX, rangeY, rangeZ, power, rotAngle, false);
  CovGridSeparated covRho  (nx, ny, nz, dx, dy, dz, rangeX, rangeY, rangeZ, power, rotAngle, false);

  if(modelSettings->getDebugFlag() == 1)
  {
    covAlpha.writeXYGrid("covAlpha");
    covBeta.writeXYGrid("covBeta");
    covRho.writeXYGrid("covRho");
  }

  CovGridSeparated covCrAlphaBeta(nx, ny, nz);
  CovGridSeparated covCrAlphaRho(nx, ny, nz);
  CovGridSeparated covCrBetaRho(nx, ny, nz);

  CKrigingAdmin kriging(*simbox, 
                        kd.getData(),
                        kd.getNumberOfData(),
                        covAlpha, covBeta, covRho, 
                        covCrAlphaBeta, covCrAlphaRho, covCrBetaRho, 
                        DataTarget_, true);

  LogKit::LogFormatted(LogKit::LOW,"\nFill volumes using kriging:");
  kriging.KrigAll(*trendCubeAlpha, *trendCubeBeta, *trendCubeRho,false);

  backModel_[0]->fillInFromRealFFTGrid(*trendCubeAlpha);
  backModel_[1]->fillInFromRealFFTGrid(*trendCubeBeta);
  backModel_[2]->fillInFromRealFFTGrid(*trendCubeRho);

  //backModel_[0] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);
  //backModel_[0]->setType(FFTGrid::PARAMETER);
  //backModel_[1] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);
  //backModel_[1]->setType(FFTGrid::PARAMETER);
  //backModel_[2] = new FFTGrid(nx, ny, nz, nxPad, nyPad, nzPad);
  //backModel_[2]->setType(FFTGrid::PARAMETER);

  delete trendCubeAlpha; 
  delete trendCubeBeta; 
  delete trendCubeRho;
}

//-------------------------------------------------------------------------------
void 
Background::calculateVerticalTrend(WellData  ** wells,
                                   float      * trend, 
                                   float        logMin,
                                   float        logMax,
                                   float        maxHz,
                                   int          nWells,
                                   int          nz,
                                   float        dz,
                                   const char * name)
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
      BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
      if (strcmp(name,"Vp") == 0)
        bl->getVerticalTrend(bl->getAlpha(), wellTrend);
      else if (strcmp(name,"Vs") == 0)
        bl->getVerticalTrend(bl->getBeta(), wellTrend);
      else if (strcmp(name,"Rho") == 0)
        bl->getVerticalTrend(bl->getRho(), wellTrend);
      else {
        LogKit::LogFormatted(LogKit::LOW,"ERROR in Background::calculateVerticalTrend(): ");
        LogKit::LogFormatted(LogKit::LOW,"Log \'%s\' requested, but no such log exists.\n",name);
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
                                                 const char * parName) 
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
  // We need must possibly require a larger min_req_points to avoid an 
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
    LogKit::LogFormatted(LogKit::LOW,"\nWARNING : The calculation of the vertical trend for parameter %s using local linear\n",parName);
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
    LogKit::LogFormatted(LogKit::LOW,"\nGlobal mean for parameter %s = %.2f\n\n",parName,global_mean);
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
Background::extrapolateTrend(const char * pName, 
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
                       pName, first_nonmissing + 1);
    if (nz - 1 - last_nonmissing > 0)
      LogKit::LogFormatted(LogKit::LOW,"Vertical trend for %s extrapolated last %d cells.\n",
                       pName, nz - 1 - last_nonmissing);
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
  std::string tmp = "BG_verticalTrend_"+name;
  std::string filename = ModelSettings::makeFullFileName(tmp+".irap");
  FILE * file = fopen(filename.c_str(), "w");
  for (int i=0 ; i<nz ; i++) {
    fprintf(file,"%8.2f %8.3f 0.00\n",(z0+i*dz),exp( trend[i] ));
  }
  fprintf(file,"999.00 999.00 999.00\n");
  fclose(file);
}

//-------------------------------------------------------------------------------
void
Background::calculateDeviationFromVerticalTrend(WellData    ** wells,
                                                const float  * globalTrend, 
                                                float        * avg_dev,
                                                int            nWells,
                                                int            nz,
                                                const char   * name)
{
  float * wellTrend = new float[nz];

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
    if (strcmp(name,"Vp") == 0)
      bl->getVerticalTrend(bl->getAlphaHighCutBackground(), wellTrend);
    else if (strcmp(name,"Vs") == 0)
      bl->getVerticalTrend(bl->getBetaHighCutBackground(), wellTrend);
    else if (strcmp(name,"Rho") == 0)
      bl->getVerticalTrend(bl->getRhoHighCutBackground(), wellTrend);
    else {
      LogKit::LogFormatted(LogKit::LOW,"ERROR in Background::calculateVerticalTrend(): ");
      LogKit::LogFormatted(LogKit::LOW,"Log \'%s\' requested, but no such log exists.\n",name);
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
Background::writeDeviationsFromVerticalTrend(WellData    ** wells,
                                             const float *  avg_dev_alpha,
                                             const float *  avg_dev_beta,
                                             const float *  avg_dev_rho,
                                             const float *  trend_alpha,
                                             const float *  trend_beta,
                                             const float *  trend_rho,
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
    LogKit::LogFormatted(LogKit::MEDIUM,"\nNOTE: A deviation may be observed even with one well since the global trend is");
    LogKit::LogFormatted(LogKit::MEDIUM,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
  }
  delete [] rel_avg_dev;
  delete [] index;
}

//-------------------------------------------------------------------------------
void           
Background::createTrendCube(Simbox      *  simbox,
                            FFTGrid     *& pFFTGrid, 
                            const float *  trend) 
{
  const int nx = simbox->getnx();
  const int ny = simbox->getny();
  const int nz = simbox->getnz();
  pFFTGrid = new FFTGrid(nx, ny, nz, nx, ny, nz);
  
  const int nzp = pFFTGrid->getNzp();
  const int nyp = pFFTGrid->getNyp();
  const int nxp = pFFTGrid->getNxp();
  pFFTGrid->createRealGrid();
  pFFTGrid->setAccessMode(FFTGrid::RANDOMACCESS);

  int i,j,k;
  for (k = 0; k < nzp; k++) 
    for (j = 0; j < nyp; j++)
      for (i = 0; i < nxp; i++)
        pFFTGrid->setRealValue(i, j, k, trend[k]);
        
  pFFTGrid->endAccess();
}

//-------------------------------------------------------------------------------
void
Background::findMeanVsVp()
{
  backModel_[0]->setAccessMode(FFTGrid::READ);
  backModel_[1]->setAccessMode(FFTGrid::READ);
  int i, j, k;
  float v1,v2;
  double mean = 0;
  int nxp = 2*(backModel_[0]->getNxp()/2+1);
  int nyp = backModel_[0]->getNyp();
  int nzp = backModel_[0]->getNzp();
  int nx = backModel_[0]->getNx();
  int ny = backModel_[0]->getNy();
  int nz = backModel_[0]->getNz();
  for(k=0;k<nzp;k++) 
    for(j=0;j<nyp;j++)
      for(i=0;i<nxp;i++) {
        v1 = backModel_[0]->getNextReal();
        v2 = backModel_[1]->getNextReal();
        if(i < nx && j < ny && k < nz)
          mean += exp(v2-v1);
      }
  mean = mean/double(nx*ny*nz);

  backModel_[0]->endAccess();
  backModel_[1]->endAccess();

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
