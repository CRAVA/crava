#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#include "lib/kriging1d.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/random/beta.hpp"
#include "nrlib/random/distribution.hpp"

#include "src/blockedlogsforzone.h"
#include "src/definitions.h"
#include "src/modelgeneral.h"
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

Background::Background(FFTGrid                       ** grids,
                       WellData                      ** wells,
                       Simbox                         * timeSimbox,
                       ModelSettings                  * modelSettings,
                       const std::vector<std::string> & surface_files)
  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
    vsvp_(RMISSING)
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];

  FFTGrid * bgAlpha;
  FFTGrid * bgBeta;
  FFTGrid * bgRho;

  generateMultizoneBackgroundModel(bgAlpha,bgBeta,bgRho,
                                   wells,
                                   timeSimbox,
                                   modelSettings,
                                   surface_files);

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
  for (int i=0 ; i<3 ; i++)
    delete backModel_[i];
}

//-------------------------------------------------------------------------------
void
Background::releaseGrids()
{
  for (int i=0 ; i<3 ; i++)
    backModel_[i] = NULL;
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
  const int   nz     = simbox->getnz();
  const int   nWells = modelSettings->getNumberOfWells();
  const float dz     = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  std::string name_vp  = "Vp";
  std::string name_vs  = "Vs";
  std::string name_rho = "Rho";

  std::vector<float *> wellTrendAlpha(nWells);
  std::vector<float *> wellTrendBeta(nWells);
  std::vector<float *> wellTrendRho(nWells);
  std::vector<float *> highCutWellTrendAlpha(nWells);
  std::vector<float *> highCutWellTrendBeta(nWells);
  std::vector<float *> highCutWellTrendRho(nWells);

  getWellTrends(wellTrendAlpha,highCutWellTrendAlpha,wells,nz,name_vp);
  getWellTrends(wellTrendBeta, highCutWellTrendBeta, wells,nz,name_vs);
  getWellTrends(wellTrendRho,  highCutWellTrendRho,  wells,nz,name_rho);

  float * trendAlpha = new float[nz];
  float * trendBeta  = new float[nz];
  float * trendRho   = new float[nz];
  float * trendVel   = new float[nz]; // Allocate (for simplicity) although not always needed

  float * avgDevAlpha = new float[nWells];
  float * avgDevBeta  = new float[nWells];
  float * avgDevRho   = new float[nWells];
  float * avgDevVel   = new float[nWells]; // Allocate (for simplicity) although not always needed

  calculateBackgroundTrend(trendAlpha,
                           avgDevAlpha,
                           nz,
                           dz,
                           modelSettings->getAlphaMin(),
                           modelSettings->getAlphaMax(),
                           modelSettings->getMaxHzBackground(),
                           wellTrendAlpha,
                           highCutWellTrendAlpha,
                           name_vp);
  calculateBackgroundTrend(trendBeta,
                           avgDevBeta,
                           nz,
                           dz,
                           modelSettings->getBetaMin(),
                           modelSettings->getBetaMax(),
                           modelSettings->getMaxHzBackground(),
                           wellTrendBeta,
                           highCutWellTrendBeta,
                           name_vs);
  calculateBackgroundTrend(trendRho,
                           avgDevRho,
                           nz,
                           dz,
                           modelSettings->getRhoMin(),
                           modelSettings->getRhoMax(),
                           modelSettings->getMaxHzBackground(),
                           wellTrendRho,
                           highCutWellTrendRho,
                           name_rho);

  bool hasVelocityTrend = velocity != NULL;
  bool write1D          = ((modelSettings->getOtherOutputFlag()& IO::BACKGROUND_TREND_1D) > 0);
  bool write3D          = ((modelSettings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);

  writeTrendsToFile(trendAlpha,simbox, write1D, write3D, hasVelocityTrend, name_vp, modelSettings->getFileGrid());
  writeTrendsToFile(trendBeta, simbox, write1D, write3D, hasVelocityTrend, name_vs, modelSettings->getFileGrid());
  writeTrendsToFile(trendRho,  simbox, write1D, write3D, hasVelocityTrend, name_rho,modelSettings->getFileGrid());

  if (velocity != NULL) {
    //
    // We still want calculateBackgroundTrend() for alpha above. By calculating
    // avgDevAlpha we can check that the bgAlpha calculated from velocity is as
    // good as or better than that calculated by crava.
    //
    calculateVelocityDeviations(velocity, wells, simbox,
                                trendVel, avgDevVel, avgDevAlpha,
                                modelSettings->getOutputGridsElastic(),
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

  std::vector<KrigingData2D> krigingDataAlpha(nz);
  std::vector<KrigingData2D> krigingDataBeta(nz);
  std::vector<KrigingData2D> krigingDataRho(nz);

  std::vector<float *>     blAlpha(nWells);   // bl = blocked logs
  std::vector<float *>     blBeta(nWells);
  std::vector<float *>     blRho(nWells);
  std::vector<float *>     vtAlpha(nWells);   // vt = vertical trend
  std::vector<float *>     vtBeta(nWells);
  std::vector<float *>     vtRho(nWells);
  std::vector<const int *> ipos(nWells);
  std::vector<const int *> jpos(nWells);
  std::vector<const int *> kpos(nWells);

  for(int i=0; i<nWells; i++) {
    vtAlpha[i] = new float[nz];
    vtBeta[i]  = new float[nz];
    vtRho[i]   = new float[nz];
  }

  std::vector<int> nBlocks(nWells);
  int              totBlocks;

  getKrigingWellTrends(blAlpha,blBeta,blRho,
                       vtAlpha,vtBeta,vtRho,
                       ipos,jpos,kpos,
                       nBlocks,totBlocks,
                       wells,nWells);

  setupKrigingData2D(krigingDataAlpha,krigingDataBeta,krigingDataRho,
                     trendAlpha,trendBeta,trendRho,
                     modelSettings->getOutputGridsElastic(),
                     nz,dz,totBlocks,nBlocks,
                     blAlpha,blBeta,blRho,
                     vtAlpha,vtBeta,vtRho,
                     ipos,jpos,kpos);

  const CovGrid2D & covGrid2D = makeCovGrid2D(simbox,
                                              modelSettings->getBackgroundVario(),
                                              modelSettings->getDebugFlag());

  makeKrigedBackground(krigingDataAlpha, bgAlpha, trendAlpha, simbox, covGrid2D, "Vp" , modelSettings->getFileGrid());
  makeKrigedBackground(krigingDataBeta , bgBeta , trendBeta , simbox, covGrid2D, "Vs" , modelSettings->getFileGrid());
  makeKrigedBackground(krigingDataRho  , bgRho  , trendRho  , simbox, covGrid2D, "Rho", modelSettings->getFileGrid());

  delete &covGrid2D;

  delete [] avgDevAlpha;
  delete [] avgDevBeta;
  delete [] avgDevRho;
  delete [] avgDevVel;

  delete [] trendAlpha;
  delete [] trendBeta;
  delete [] trendRho;
  delete [] trendVel;

  for(int i=0; i<nWells; i++) {

    delete [] wellTrendAlpha[i];
    delete [] wellTrendBeta[i];
    delete [] wellTrendRho[i];

    delete [] highCutWellTrendAlpha[i];
    delete [] highCutWellTrendBeta[i];
    delete [] highCutWellTrendRho[i];

    delete [] blAlpha[i];
    delete [] blBeta[i];
    delete [] blRho[i];

    delete [] vtAlpha[i];
    delete [] vtBeta[i];
    delete [] vtRho[i];

  }
}
//-------------------------------------------------------------------------------
void
Background::generateMultizoneBackgroundModel(FFTGrid                       *& bgAlpha,
                                             FFTGrid                       *& bgBeta,
                                             FFTGrid                       *& bgRho,
                                             WellData                      ** wells,
                                             Simbox                         * simbox,
                                             ModelSettings                  * modelSettings,
                                             const std::vector<std::string> & surface_files)
{
  LogKit::LogFormatted(LogKit::Low,"Multizone background model:\n");

  std::vector<int> correlation_structure = modelSettings->getCorrelationStructure();
  std::vector<int> erosion_priority      = modelSettings->getErosionPriority();

  int    nWells    = modelSettings->getNumberOfWells();
  int    nZones    = static_cast<int>(correlation_structure.size()) - 1;
  float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save mamory

  std::vector<Surface> surface(nZones+1);
  for(int i=0; i<nZones+1; i++)
    surface[i] = Surface(surface_files[i]);

  //Check that surfaces cover inversion area
  checkSurfaces(surface, simbox);

  std::vector<StormContGrid> alpha_zones(nZones);
  std::vector<StormContGrid> beta_zones(nZones);
  std::vector<StormContGrid> rho_zones(nZones);

  BuildSeismicPropertyZones(alpha_zones,
                            beta_zones,
                            rho_zones,
                            surface,
                            correlation_structure,
                            simbox,
                            dz);

  std::vector<Surface *> eroded_surfaces(nZones+1);
  for(int i=0; i<nZones+1; i++)
    eroded_surfaces[i] = NULL;

  ErodeAllSurfaces(eroded_surfaces,
                   erosion_priority,
                   surface,
                   simbox);

  const CovGrid2D & covGrid2D = makeCovGrid2D(simbox, modelSettings->getBackgroundVario(), modelSettings->getDebugFlag());

  std::string name_vp  = "Vp";
  std::string name_vs  = "Vs";
  std::string name_rho = "Rho";

  std::vector<float *> trendAlphaZone(nZones);
  std::vector<float *> trendBetaZone(nZones);
  std::vector<float *> trendRhoZone(nZones);

  for(int i=0; i<nZones; i++) {
    LogKit::LogFormatted(LogKit::Low,"\nZone%2d:",i+1);

    int nz = static_cast<int>(alpha_zones[i].GetNK());

    std::vector<float *> wellTrendAlpha(nWells);
    std::vector<float *> wellTrendBeta(nWells);
    std::vector<float *> wellTrendRho(nWells);

    std::vector<float *> highCutWellTrendAlpha(nWells);
    std::vector<float *> highCutWellTrendBeta(nWells);
    std::vector<float *> highCutWellTrendRho(nWells);

    StormContGrid eroded_zone;

    BuildErodedZones(eroded_zone,
                     eroded_surfaces,
                     nz,
                     simbox,
                     i);

    std::vector<bool> hitZone(nWells);
    checkWellHitsZone(hitZone, wells, eroded_zone, nWells);

    std::vector<BlockedLogsForZone *> blocked_logs(nWells);

    getWellTrendsZone(blocked_logs, wellTrendAlpha, highCutWellTrendAlpha, wells, eroded_zone, hitZone, nz, name_vp,  i);
    getWellTrendsZone(blocked_logs, wellTrendBeta,  highCutWellTrendBeta,  wells, eroded_zone, hitZone, nz, name_vs,  i);
    getWellTrendsZone(blocked_logs, wellTrendRho,   highCutWellTrendRho,   wells, eroded_zone, hitZone, nz, name_rho, i);

    trendAlphaZone[i] = new float[nz];
    trendBetaZone[i]  = new float[nz];
    trendRhoZone[i]   = new float[nz];

    float * avgDevAlphaZone = new float[nWells];
    float * avgDevBetaZone  = new float[nWells];
    float * avgDevRhoZone   = new float[nWells];

    calculateBackgroundTrend(trendAlphaZone[i],
                             avgDevAlphaZone,
                             nz,
                             dz,
                             modelSettings->getAlphaMin(),
                             modelSettings->getAlphaMax(),
                             modelSettings->getMaxHzBackground(),
                             wellTrendAlpha,
                             highCutWellTrendAlpha,
                             name_vp);
    calculateBackgroundTrend(trendBetaZone[i],
                             avgDevBetaZone,
                             nz,
                             dz,
                             modelSettings->getBetaMin(),
                             modelSettings->getBetaMax(),
                             modelSettings->getMaxHzBackground(),
                             wellTrendBeta,
                             highCutWellTrendBeta,
                             name_vs);
    calculateBackgroundTrend(trendRhoZone[i],
                             avgDevRhoZone,
                             nz,
                             dz,
                             modelSettings->getRhoMin(),
                             modelSettings->getRhoMax(),
                             modelSettings->getMaxHzBackground(),
                             wellTrendRho,
                             highCutWellTrendRho,
                             name_rho);

    writeDeviationsFromVerticalTrend(avgDevAlphaZone,
                                     avgDevBetaZone,
                                     avgDevRhoZone,
                                     trendAlphaZone[i],
                                     trendBetaZone[i],
                                     trendRhoZone[i],
                                     wells,
                                     nWells,
                                     nz);

    std::vector<float *>     blAlpha(nWells);   // bl = blocked logs
    std::vector<float *>     blBeta(nWells);
    std::vector<float *>     blRho(nWells);

    std::vector<float *>     vtAlpha(nWells);   // vt = vertical trend
    std::vector<float *>     vtBeta(nWells);
    std::vector<float *>     vtRho(nWells);

    std::vector<const int *> ipos(nWells);
    std::vector<const int *> jpos(nWells);
    std::vector<const int *> kpos(nWells);

    int              totBlocks;
    std::vector<int> nBlocks(nWells);

    getKrigingWellTrendsZone(blocked_logs,
                             blAlpha,blBeta,blRho,
                             vtAlpha,vtBeta,vtRho,
                             ipos,jpos,kpos,
                             nBlocks,totBlocks,nz);

    std::vector<KrigingData2D> krigingDataAlpha(nz);
    std::vector<KrigingData2D> krigingDataBeta(nz);
    std::vector<KrigingData2D> krigingDataRho(nz);

    setupKrigingData2D(krigingDataAlpha,krigingDataBeta,krigingDataRho,
                       trendAlphaZone[i],trendBetaZone[i],trendRhoZone[i],
                       modelSettings->getOutputGridsElastic(),
                       nz,dz,totBlocks,nBlocks,
                       blAlpha,blBeta,blRho,
                       vtAlpha,vtBeta,vtRho,
                       ipos,jpos,kpos);

    makeKrigedZone(krigingDataAlpha, trendAlphaZone[i], alpha_zones[i], covGrid2D);
    makeKrigedZone(krigingDataBeta , trendBetaZone[i] , beta_zones[i],  covGrid2D);
    makeKrigedZone(krigingDataRho  , trendRhoZone[i]  , rho_zones[i],   covGrid2D);

    delete [] avgDevAlphaZone;
    delete [] avgDevBetaZone;
    delete [] avgDevRhoZone;

    for(int j=0; j<nWells; j++) {
      delete [] wellTrendAlpha[j];
      delete [] wellTrendBeta[j];
      delete [] wellTrendRho[j];

      delete [] highCutWellTrendAlpha[j];
      delete [] highCutWellTrendBeta[j];
      delete [] highCutWellTrendRho[j];

      delete [] blAlpha[j];
      delete [] blBeta[j];
      delete [] blRho[j];

      delete [] vtAlpha[j];
      delete [] vtBeta[j];
      delete [] vtRho[j];

      delete blocked_logs[j];
    }
  }

  MakeMultizoneBackground(bgAlpha,bgBeta,bgRho,
                          alpha_zones,beta_zones,rho_zones,
                          simbox,
                          erosion_priority,
                          surface,
                          modelSettings->getSurfaceUncertainty(),
                          modelSettings->getFileGrid(),
                          "multizone");


  bool write3D = ((modelSettings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);

  if(write3D) {

    std::vector<int> nz_zone(nZones);
    for(int i=0; i<nZones; i++)
      nz_zone[i] = static_cast<int>(alpha_zones[i].GetNK());

    alpha_zones.clear();
    beta_zones.clear();
    rho_zones.clear();

    writeMultizoneTrendsToFile(trendAlphaZone,trendBetaZone,trendRhoZone,
                               nz_zone,
                               simbox,
                               erosion_priority,
                               surface,
                               modelSettings->getSurfaceUncertainty(),
                               modelSettings->getFileGrid());

  }

  delete &covGrid2D;

  for(int i=0; i<nZones; i++) {
    delete [] trendAlphaZone[i];
    delete [] trendBetaZone[i];
    delete [] trendRhoZone[i];

    delete eroded_surfaces[i];
  }

}

void
Background::MakeMultizoneBackground(FFTGrid                         *& bgAlpha,
                                    FFTGrid                         *& bgBeta,
                                    FFTGrid                         *& bgRho,
                                    const std::vector<StormContGrid> & alpha_zones,
                                    const std::vector<StormContGrid> & beta_zones,
                                    const std::vector<StormContGrid> & rho_zones,
                                    const Simbox                     * simbox,
                                    const std::vector<int>           & erosion_priority,
                                    const std::vector<Surface>       & surface,
                                    const std::vector<double>        & surface_uncertainty,
                                    const bool                         isFile,
                                    const std::string                & type) const
{

  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  int nZones = static_cast<int>(alpha_zones.size());

  int nx = simbox->getnx();
  int ny = simbox->getny();
  int nz = simbox->getnz();

  const int nxp  = nx;
  const int nyp  = ny;
  const int nzp  = nz;
  const int rnxp = 2*(nxp/2 + 1);

  float monitorSize = std::max(1.0f, static_cast<float>(rnxp)*0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";




  bgAlpha = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  bgBeta  = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  bgRho   = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);

  bgAlpha->createRealGrid();
  bgBeta ->createRealGrid();
  bgRho  ->createRealGrid();

  bgAlpha->setType(FFTGrid::PARAMETER);
  bgBeta ->setType(FFTGrid::PARAMETER);
  bgRho  ->setType(FFTGrid::PARAMETER);

  bgAlpha->setAccessMode(FFTGrid::WRITE);
  bgBeta ->setAccessMode(FFTGrid::WRITE);
  bgRho  ->setAccessMode(FFTGrid::WRITE);

  // Beta distributed uncertainty on each surface
  std::vector<NRLib::Beta> horizon_distributions(nZones+1);
  for(int zone=0; zone<nZones+1; zone++) {
    horizon_distributions[zone] = NRLib::Beta(-surface_uncertainty[zone], surface_uncertainty[zone], 2, 2);
  }

  for(int i=0; i<rnxp; i++) {

    for(int j=0; j<nyp; j++) {

      for(int k=0; k<nzp; k++) {

        if(i<nx) {

          double x;
          double y;
          double z;
          simbox->getCoord(i,j,k,x,y,z);

          std::vector<double> z_relative(nZones+1);
          for(int zone=0; zone<nZones+1; zone++)
            z_relative[zone] = z - surface[zone].GetZ(x,y);

          std::vector<double> zone_probability(nZones);

          ComputeZoneProbability(z_relative, horizon_distributions, erosion_priority, zone_probability);

          double vp  = 0;
          double vs  = 0;
          double rho = 0;

          for(int zone=0; zone<nZones; zone++) {
            // if-test to speed up computations
            if(zone_probability[zone] > 0) {
              double vp_zone  = alpha_zones[zone].GetValueZInterpolated(x,y,z);
              double vs_zone  = beta_zones[zone].GetValueZInterpolated(x,y,z);
              double rho_zone = rho_zones[zone].GetValueZInterpolated(x,y,z);

              vp  +=  vp_zone  * zone_probability[zone];
              vs  +=  vs_zone  * zone_probability[zone];
              rho +=  rho_zone * zone_probability[zone];
            }
          }

          bgAlpha->setRealValue(i, j, k, float(vp));
          bgBeta ->setRealValue(i, j, k, float(vs));
          bgRho  ->setRealValue(i, j, k, float(rho));
        }

        else {
          bgAlpha->setRealValue(i, j, k, 0);
          bgBeta ->setRealValue(i, j, k, 0);
          bgRho  ->setRealValue(i, j, k, 0);
        }
      }
    }

    // Log progress
    if (i+1 >= static_cast<int>(nextMonitor)) {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  bgAlpha->endAccess();
  bgBeta ->endAccess();
  bgRho  ->endAccess();
}
//---------------------------------------------------------------------------
void
Background::ComputeZoneProbability(const std::vector<double>      & z,
                                   const std::vector<NRLib::Beta> & horizon_distributions,
                                   const std::vector<int>         & erosion_priority,
                                   std::vector<double>            & zone_probability) const
{

  int nZones = static_cast<int>(zone_probability.size());

  for(int zone=0; zone<nZones; zone++) {
    //Initialize with probability that we are below top surface for zone
    double prob = horizon_distributions[zone].Cdf(z[zone]);

    //Multiply with probability that we are above base surface for zone
    prob *= (1-horizon_distributions[zone+1].Cdf(z[zone+1]));

    //We may be eroded from above. Must consider the surfaces that
    //1. Are above top in the standard sequence.
    //2. Have lower erosion priority number than the top.
    //3. Have no horizons with lower erosion priority number between it and top.
    int min_erosion = erosion_priority[zone];
    for(int prev_hor = zone-1; prev_hor >=0; prev_hor--) {
      if(erosion_priority[prev_hor] < min_erosion) {
        prob *= horizon_distributions[prev_hor].Cdf(z[prev_hor]);
        min_erosion = erosion_priority[prev_hor]; //Those with higher number stop in this
      }
    }

    //We may be eroded from below. Must consider the surfaces that
    //1. Are below base in the standard sequence.
    //2. Have lower erosion priority number than the base.
    //3. Have no horizons with lower erosion priority number between it and base.
    min_erosion = erosion_priority[zone+1];
    for(int late_hor = zone+2; late_hor < nZones+1; late_hor++) {
      if(erosion_priority[late_hor] < min_erosion) {
        prob *= (1-horizon_distributions[late_hor].Cdf(z[late_hor]));
        min_erosion = erosion_priority[late_hor]; //Those with higher number stop in this
      }
    }

    zone_probability[zone] = prob;
  }
}

//---------------------------------------------------------------------------
void
Background::BuildErodedZones(StormContGrid                & eroded_zone,
                             const std::vector<Surface *> & eroded_surfaces,
                             const int                    & nz,
                             const Simbox                 * simbox,
                             const int                    & i) const
{
  int    nx        = simbox->getnx();
  int    ny        = simbox->getny();
  double x_min     = simbox->GetXMin();
  double y_min     = simbox->GetYMin();
  double lx        = simbox->GetLX();
  double ly        = simbox->GetLY();
  double angle     = simbox->getAngle();


  NRLib::Volume volume(x_min, y_min, lx, ly, *eroded_surfaces[i], *eroded_surfaces[i+1], angle);

  eroded_zone = StormContGrid(volume, nx, ny, nz);

}
//---------------------------------------------------------------------------
void
Background::ErodeAllSurfaces(std::vector<Surface *>     & eroded_surfaces,
                             const std::vector<int>     & erosion_priority,
                             const std::vector<Surface> & surface,
                             const Simbox               * simbox) const
{
  int    nSurf     = static_cast<int>(eroded_surfaces.size());

  for(int i=0; i<nSurf; i++) {
    int l=0;
    while(i+1 != erosion_priority[l])
      l++;

    Surface * temp_surface = new Surface(surface[l]);

    //Find closest eroded surface downward
    for(int k=l+1; k<nSurf; k++) {
      if(eroded_surfaces[k] != NULL) {
        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, false);
        break;
      }
    }
    //Find closest eroded surface upward
    for(int k=l-1; k>=0; k--) {
      if(eroded_surfaces[k] != NULL) {
        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, true);
        break;
      }
    }
    eroded_surfaces[l] = temp_surface;
  }
}

//---------------------------------------------------------------------------

void
Background::BuildSeismicPropertyZones(std::vector<StormContGrid> & alpha_zones,
                                      std::vector<StormContGrid> & beta_zones,
                                      std::vector<StormContGrid> & rho_zones,
                                      const std::vector<Surface> & surface,
                                      const std::vector<int>     & correlation_structure,
                                      const Simbox               * simbox,
                                      const float                & dz) const
{
  int    nZones    = static_cast<int>(alpha_zones.size());
  int    nx        = simbox->getnx();
  int    ny        = simbox->getny();
  double x_min     = simbox->GetXMin();
  double y_min     = simbox->GetYMin();
  double lx        = simbox->GetLX();
  double ly        = simbox->GetLY();
  double angle     = simbox->getAngle();

  for(int i=1; i<nZones+1; i++) {
    Surface temp_top;
    Surface temp_base;
    double  x;
    double  y;
    double  z_top;
    double  z_base;

    Surface top  = surface[i-1];
    Surface base = surface[i];

    double top_missing  = top.GetMissingValue();
    double base_missing = base.GetMissingValue();

    //Find maximum distance between the surfaces
    double max_distance = 0;

    for(int j=0; j<nx; j++) {
      for(int k=0; k<ny; k++) {
        simbox->getXYCoord(j,k,x,y);

        z_top  = top.GetZ(x,y);
        z_base = base.GetZ(x,y);

        if(z_base-z_top > max_distance) {
          if(z_top != top_missing && z_base != base_missing)
            max_distance = z_base-z_top;
        }
      }
    }

    //Make new top and base surfaces
    if(correlation_structure[i] == ModelSettings::TOP) {
      temp_top  = top;
      temp_base = top;
      temp_base.Add(max_distance);
    }
    else if(correlation_structure[i] == ModelSettings::BASE) {
      temp_top  = base;
      temp_top.Subtract(max_distance);
      temp_base = base;
    }
    else {
      temp_top  = top;
      temp_base = base;
    }

    NRLib::Volume volume(x_min, y_min, lx, ly, temp_top, temp_base, angle);

    int nz_zone = static_cast<int>(std::ceil(max_distance/dz));

    alpha_zones[i-1] = StormContGrid(volume, nx, ny, nz_zone);
    beta_zones[i-1]  = StormContGrid(volume, nx, ny, nz_zone);
    rho_zones[i-1]   = StormContGrid(volume, nx, ny, nz_zone);
  }
}

//---------------------------------------------------------------------------
void
Background::checkSurfaces(const std::vector<Surface> & surface,
                          const Simbox               * simbox) const
{
  int    nSurf     = static_cast<int>(surface.size());
  int    nx        = simbox->getnx();
  int    ny        = simbox->getny();

  for(int i=0; i<nSurf; i++) {

    double  x;
    double  y;
    double  z;

    Surface surf  = surface[i];

    double surf_missing  = surf.GetMissingValue();

    for(int j=0; j<nx; j++) {
      for(int k=0; k<ny; k++) {

        simbox->getXYCoord(j,k,x,y);

        z  = surf.GetZInside(x,y);

        if(z == surf_missing) {
          const std::string name = surf.GetName();

          LogKit::LogFormatted(LogKit::Low,"ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }
      }
    }
  }
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
  if((outputFlag & IO::BACKGROUND_TREND) > 0) {
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

  LogKit::LogFormatted(LogKit::Low,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
  LogKit::LogFormatted(LogKit::Low,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
  LogKit::LogFormatted(LogKit::Low,"Well             TrendFromFile  TrendFromData\n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------\n");
  for (int i=0 ; i<nWells ; i++)
    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f          %5.1f\n",
                         wells[i]->getWellname().c_str(),avgDevVel[i],avgDevAlpha[i]);
}

//---------------------------------------------------------------------------
void
Background::calculateBackgroundTrend(float              * trend,
                                     float              * avgDev,
                                     const int            nz,
                                     const float          dz,
                                     float                logMin,
                                     float                logMax,
                                     float                maxHz,
                                     std::vector<float *> wellTrend,
                                     std::vector<float *> highCutWellTrend,
                                     const std::string  & name)
{

  calculateVerticalTrend(wellTrend,
                         trend,
                         logMin,
                         logMax,
                         maxHz,
                         nz,
                         dz,
                         name);


  calculateDeviationFromVerticalTrend(highCutWellTrend, trend, avgDev, nz);


}
//---------------------------------------------------------------------------
void
Background::getKrigingWellTrends(std::vector<float *>     & blAlpha,
                                 std::vector<float *>     & blBeta,
                                 std::vector<float *>     & blRho,
                                 std::vector<float *>     & vtAlpha,
                                 std::vector<float *>     & vtBeta,
                                 std::vector<float *>     & vtRho,
                                 std::vector<const int *> & ipos,
                                 std::vector<const int *> & jpos,
                                 std::vector<const int *> & kpos,
                                 std::vector<int>         & nBlocks,
                                 int                      & totBlocks,
                                 WellData                ** wells,
                                 const int                & nWells) const
{
  int maxBlocks = 0;
  totBlocks     = 0;

  for (int w = 0 ; w < nWells ; w++) {
    nBlocks[w] = wells[w]->getBlockedLogsExtendedBG()->getNumberOfBlocks();
    totBlocks += nBlocks[w];
    if (nBlocks[w] > maxBlocks)
      maxBlocks = nBlocks[w];
  }

  for(int i=0; i<nWells; i++) {
    blAlpha[i] = new float[maxBlocks];
    blBeta[i]  = new float[maxBlocks];
    blRho[i]   = new float[maxBlocks];
  }

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();

    Utils::copyVector(bl->getAlphaHighCutBackground(), blAlpha[w], nBlocks[w]);
    Utils::copyVector(bl->getBetaHighCutBackground(),  blBeta[w],  nBlocks[w]);
    Utils::copyVector(bl->getRhoHighCutBackground(),   blRho[w],   nBlocks[w]);
    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    bl->getVerticalTrend(blAlpha[w], vtAlpha[w]);
    bl->getVerticalTrend(blBeta[w],  vtBeta[w]);
    bl->getVerticalTrend(blRho[w],   vtRho[w]);

    ipos[w] = bl->getIpos();
    jpos[w] = bl->getJpos();
    kpos[w] = bl->getKpos();

  }
}
//---------------------------------------------------------------------------
void
Background::getKrigingWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
                                     std::vector<float *>              & blAlpha,
                                     std::vector<float *>              & blBeta,
                                     std::vector<float *>              & blRho,
                                     std::vector<float *>              & vtAlpha,
                                     std::vector<float *>              & vtBeta,
                                     std::vector<float *>              & vtRho,
                                     std::vector<const int *>          & ipos,
                                     std::vector<const int *>          & jpos,
                                     std::vector<const int *>          & kpos,
                                     std::vector<int>                  & nBlocks,
                                     int                               & totBlocks,
                                     const int                         & nz) const
{
  int nWells    = static_cast<int>(bl.size());
  int maxBlocks = 0;
  totBlocks     = 0;

  for (int w = 0 ; w < nWells ; w++) {
    if(bl[w] != NULL) {
      nBlocks[w] = bl[w]->getNumberOfBlocks();
      totBlocks += nBlocks[w];
      if (nBlocks[w] > maxBlocks)
        maxBlocks = nBlocks[w];
    }
    else
      nBlocks[w] = 0;
  }

  for (int w = 0; w < nWells; w++) {
    if(bl[w] != NULL) {
      std::vector<float> blAlphaHighCut = bl[w]->getAlphaHighCutBackground();
      std::vector<float> blBetaHighCut  = bl[w]->getBetaHighCutBackground();
      std::vector<float> blRhoHighCut   = bl[w]->getRhoHighCutBackground();

      float * blAlphaCopy = new float[maxBlocks];
      float * blBetaCopy  = new float[maxBlocks];
      float * blRhoCopy   = new float[maxBlocks];

      for(int i=0; i<nBlocks[w]; i++) {
        blAlphaCopy[i] = blAlphaHighCut[i];
        blBetaCopy[i]  = blBetaHighCut[i];
        blRhoCopy[i]   = blRhoHighCut[i];
      }
      blAlpha[w] = blAlphaCopy;
      blBeta[w]  = blBetaCopy;
      blRho[w]   = blRhoCopy;
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      vtAlpha[w] = new float[nz];
      vtBeta[w]  = new float[nz];
      vtRho[w]   = new float[nz];

      bl[w]->getVerticalTrend(blAlphaHighCut, vtAlpha[w]);
      bl[w]->getVerticalTrend(blBetaHighCut,  vtBeta[w]);
      bl[w]->getVerticalTrend(blRhoHighCut,   vtRho[w]);

      ipos[w] = bl[w]->getIpos();
      jpos[w] = bl[w]->getJpos();
      kpos[w] = bl[w]->getKpos();
    }
    else {
      vtAlpha[w] = NULL;
      vtBeta[w]  = NULL;
      vtRho[w]   = NULL;

      blAlpha[w] = NULL;
      blBeta[w]  = NULL;
      blRho[w]   = NULL;
    }
  }
}
//---------------------------------------------------------------------------
void
Background::getWellTrends(std::vector<float *> & wellTrend,
                          std::vector<float *> & highCutWellTrend,
                          WellData            ** wells,
                          const int            & nz,
                          const std::string    & name) const
{
  int nWells = static_cast<int>(wellTrend.size());
  int iWells = 0;

  for (int w = 0 ; w < nWells ; w++) {
    if (wells[w]->getUseForBackgroundTrend()) {
      BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
      if(bl != NULL) {
        wellTrend[w] = new float[nz];
        if (name == "Vp")
          bl->getVerticalTrend(bl->getAlpha(), wellTrend[w]);
        else if (name == "Vs")
          bl->getVerticalTrend(bl->getBeta(), wellTrend[w]);
        else if (name == "Rho")
          bl->getVerticalTrend(bl->getRho(), wellTrend[w]);
        else {
          LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrends(): ");
          LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
          exit(1);
        }
        iWells++;
      }
      else
        wellTrend[w] = NULL;
    }
    else
      wellTrend[w] = NULL;
  }
  if(iWells == 0) {
    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::getWellTrends(): There are no wells\n");
    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
    exit(1);
  }

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    if(bl != NULL) {
      highCutWellTrend[w] = new float[nz];
      if (name == "Vp")
        bl->getVerticalTrend(bl->getAlphaHighCutBackground(), highCutWellTrend[w]);
      else if (name == "Vs")
        bl->getVerticalTrend(bl->getBetaHighCutBackground(), highCutWellTrend[w]);
      else if (name == "Rho")
        bl->getVerticalTrend(bl->getRhoHighCutBackground(), highCutWellTrend[w]);
      else {
        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrends(): ");
        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
        exit(1);
      }
    }
    else
      highCutWellTrend[w] = NULL;
  }
}
//---------------------------------------------------------------------------
void
Background::getWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
                              std::vector<float *>              & wellTrend,
                              std::vector<float *>              & highCutWellTrend,
                              WellData                         ** wells,
                              StormContGrid                     & eroded_zone,
                              const std::vector<bool>           & hitZone,
                              const int                         & nz,
                              const std::string                 & name,
                              const int                         & i) const
{

  int nValidWellsInZone = 0;
  int nWells            = static_cast<int>(bl.size());

  std::vector<bool> use_for_background(nWells);

  for(int w=0; w<nWells; w++) {
    if(hitZone[w] == true) {
      bl[w] = new BlockedLogsForZone(wells[w], eroded_zone);
      nValidWellsInZone++;
    }
    else
      bl[w] = NULL;

    if (wells[w]->getUseForBackgroundTrend())
      use_for_background[w] = true;
    else
      use_for_background[w] = false;
  }

  if(nValidWellsInZone == 0) {
    LogKit::LogFormatted(LogKit::Low, "Invalid multizone background estimation: No well hits zone number "+NRLib::ToString(i+1)+"\n");
    exit(1);
  }

  int iWells = 0;

  for (int w = 0 ;w < nWells; w++) {
    if (use_for_background[w] == true) {
      if(bl[w] != NULL) {
        wellTrend[w] = new float[nz];
        if (name == "Vp")
          bl[w]->getVerticalTrend(bl[w]->getAlpha(), wellTrend[w]);
        else if (name == "Vs")
          bl[w]->getVerticalTrend(bl[w]->getBeta(), wellTrend[w]);
        else if (name == "Rho")
          bl[w]->getVerticalTrend(bl[w]->getRho(), wellTrend[w]);
        else {
          LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrendsZone(): ");
          LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
          exit(1);
        }
        iWells++;
      }
      else wellTrend[w] = NULL;
    }
    else
      wellTrend[w] = NULL;
  }
  if(iWells == 0) {
    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::getWellTrendsZone(): There are no wells\n");
    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
    exit(1);
  }


  for (int w = 0; w < nWells; w++) {
    if(bl[w] != NULL) {
      highCutWellTrend[w] = new float[nz];
      if (name == "Vp")
        bl[w]->getVerticalTrend(bl[w]->getAlphaHighCutBackground(), highCutWellTrend[w]);
      else if (name == "Vs")
        bl[w]->getVerticalTrend(bl[w]->getBetaHighCutBackground(), highCutWellTrend[w]);
      else if (name == "Rho")
        bl[w]->getVerticalTrend(bl[w]->getRhoHighCutBackground(), highCutWellTrend[w]);
      else {
        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrendsZone(): ");
        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
        exit(1);
      }
    }
    else
      highCutWellTrend[w] = NULL;
  }
}
//---------------------------------------------------------------------------
void
Background::checkWellHitsZone(std::vector<bool>   & hitZone,
                              WellData           ** wells,
                              StormContGrid       & eroded_zone,
                              const int           & nWells) const
{
  for(int w=0; w<nWells; w++) {
    if(wells[w]->checkStormgrid(eroded_zone) == 0) {
      hitZone[w] = true;
    }
    else
      hitZone[w] = false;
  }
}
//---------------------------------------------------------------------------
void
Background::writeTrendsToFile(float             * trend,
                              Simbox            * simbox,
                              bool                write1D,
                              bool                write3D,
                              bool                hasVelocityTrend,
                              const std::string & name,
                              bool                isFile)
{
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());
  const int   nz = simbox->getnz();

  if(write1D == true) {
    writeVerticalTrend(trend, dz, nz, name);
  }

  if(write3D == true && !(name=="Vp" && hasVelocityTrend))
  {
    const int nx = simbox->getnx();
    const int ny = simbox->getny();
    FFTGrid * trendGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nx, ny, nz, isFile);
    fillInVerticalTrend(trendGrid, trend);
    FFTGrid * expTrend = copyFFTGrid(trendGrid, true, isFile);
    delete trendGrid;
    std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + name;
    expTrend->writeFile(fileName, IO::PathToBackground(), simbox);
    delete expTrend;
  }
}
//---------------------------------------------------------------------------
void
Background::writeMultizoneTrendsToFile(const std::vector<float *>   alpha_zones,
                                       const std::vector<float *>   beta_zones,
                                       const std::vector<float *>   rho_zones,
                                       const std::vector<int>     & nz,
                                       Simbox                     * simbox,
                                       const std::vector<int>     & erosion_priority,
                                       const std::vector<Surface> & surface,
                                       const std::vector<double>  & surface_uncertainty,
                                       const bool                   isFile) const
{
  int nZones = static_cast<int>(alpha_zones.size());

  std::vector<StormContGrid> alpha_trend_zone(nZones);
  std::vector<StormContGrid> beta_trend_zone(nZones);
  std::vector<StormContGrid> rho_trend_zone(nZones);

  for(int i=0; i<nZones; i++) {
    makeTrendZoneGrid(alpha_trend_zone[i], surface[i], surface[i+1], nz[i], simbox);
    makeTrendZoneGrid(beta_trend_zone[i],  surface[i], surface[i+1], nz[i], simbox);
    makeTrendZoneGrid(rho_trend_zone[i],   surface[i], surface[i+1], nz[i], simbox);
  }

  for(int i=0; i<nZones; i++) {
    makeTrendZone(alpha_zones[i], alpha_trend_zone[i]);
    makeTrendZone(beta_zones[i],  beta_trend_zone[i]);
    makeTrendZone(rho_zones[i],   rho_trend_zone[i]);
  }

  FFTGrid * trendAlpha;
  FFTGrid * trendBeta;
  FFTGrid * trendRho;

  MakeMultizoneBackground(trendAlpha,trendBeta,trendRho,
                          alpha_trend_zone,beta_trend_zone,rho_trend_zone,
                          simbox,
                          erosion_priority,
                          surface,
                          surface_uncertainty,
                          isFile,
                          "trend in multizone");



  FFTGrid * expTrendAlpha = copyFFTGrid(trendAlpha, true, isFile);
  FFTGrid * expTrendBeta  = copyFFTGrid(trendBeta,  true, isFile);
  FFTGrid * expTrendRho   = copyFFTGrid(trendRho,   true, isFile);

  std::string fileNameAlpha = IO::PrefixBackground() + IO::PrefixTrend() + "Vp";
  std::string fileNameBeta  = IO::PrefixBackground() + IO::PrefixTrend() + "Vs";
  std::string fileNameRho   = IO::PrefixBackground() + IO::PrefixTrend() + "Rho";

  expTrendAlpha->writeFile(fileNameAlpha, IO::PathToBackground(), simbox);
  expTrendBeta ->writeFile(fileNameBeta,  IO::PathToBackground(), simbox);
  expTrendRho  ->writeFile(fileNameRho,   IO::PathToBackground(), simbox);

  delete expTrendAlpha;
  delete expTrendBeta;
  delete expTrendRho;

  delete trendAlpha;
  delete trendBeta;
  delete trendRho;
}

//-------------------------------------------------------------------------------
void
Background::setupKrigingData2D(std::vector<KrigingData2D>     & krigingDataAlpha,
                               std::vector<KrigingData2D>     & krigingDataBeta,
                               std::vector<KrigingData2D>     & krigingDataRho,
                               float                          * trendAlpha,
                               float                          * trendBeta,
                               float                          * trendRho,
                               const int                        outputFlag,
                               const int                      & nz,
                               const float                    & dz,
                               const int                      & totBlocks,
                               const std::vector<int>         & nBlocks,
                               const std::vector<float *>     & blAlpha,
                               const std::vector<float *>     & blBeta,
                               const std::vector<float *>     & blRho,
                               const std::vector<float *>     & vtAlpha,
                               const std::vector<float *>     & vtBeta,
                               const std::vector<float *>     & vtRho,
                               const std::vector<const int *>   ipos,
                               const std::vector<const int *>   jpos,
                               const std::vector<const int *>   kpos) const
{
  //
  // Although unnecessary, we have chosen to set up kriging data from
  // Vp, Vs and Rho simultaneously. This gives code easier to read.
  //
  const int nWells = static_cast<int>(blAlpha.size());

  KrigingData3D forLogging(totBlocks);

  for (int w = 0 ; w < nWells ; w++) {

    if(vtAlpha[w] != NULL) {
      float * vtAlphaWell = vtAlpha[w];
      float * vtBetaWell  = vtBeta[w];
      float * vtRhoWell   = vtRho[w];

      float * blAlphaWell = blAlpha[w];
      float * blBetaWell  = blBeta[w];
      float * blRhoWell   = blRho[w];
      //
      // Kriging vertical trend (vt....) against global vertical trend (trend...)
      //
      Kriging1D::krigVector(vtAlphaWell, trendAlpha, nz, dz);
      Kriging1D::krigVector(vtBetaWell,  trendBeta,  nz, dz);
      Kriging1D::krigVector(vtRhoWell,   trendRho,   nz, dz);
      //
      // Use kriged vertical trend where original log is not defined.
      //
      const int * iposWell = ipos[w];
      const int * jposWell = jpos[w];
      const int * kposWell = kpos[w];

      for (int m = 0 ; m < nBlocks[w] ; m++) {
        int i = iposWell[m];
        int j = jposWell[m];
        int k = kposWell[m];

        if (blAlphaWell[m] == RMISSING)
          blAlphaWell[m] = vtAlphaWell[k];
        if (blBetaWell[m] == RMISSING)
          blBetaWell[m] = vtBetaWell[k];
        if (blRhoWell[m] == RMISSING)
          blRhoWell[m] = vtRhoWell[k];

        krigingDataAlpha[k].addData(i, j, blAlphaWell[m]);
        krigingDataBeta[k].addData(i, j, blBetaWell[m]);
        krigingDataRho[k].addData(i, j, blRhoWell[m]);
      }

      forLogging.addData(blAlphaWell, blBetaWell, blRhoWell,
        iposWell,jposWell,kposWell,
        nBlocks[w]);
    }

    for (int k=0 ; k<nz ; k++) {
      krigingDataAlpha[k].findMeanValues();
      krigingDataBeta[k].findMeanValues();
      krigingDataRho[k].findMeanValues();
    }
  }

  if((outputFlag & IO::BACKGROUND) > 0) {
    forLogging.divide();
    std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + IO::SuffixGeneralData();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    forLogging.writeToFile(fileName);
  }
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
                                 const std::string                & type,
                                 bool                               isFile) const
{
  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  const int    nx   = simbox->getnx();
  const int    ny   = simbox->getny();
  const int    nz   = simbox->getnz();

  const int    nxp  = nx;
  const int    nyp  = ny;
  const int    nzp  = nz;
  const int    rnxp = 2*(nxp/2 + 1);

  const double x0   = simbox->getx0();
  const double y0   = simbox->gety0();
  const double lx   = simbox->getlx();
  const double ly   = simbox->getly();

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

  bgGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
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
          bgGrid->setNextReal(0);  //dummy in padding (but there is no padding)
      }
    }

    // Log progress
    if (k+1 >= static_cast<int>(nextMonitor))
    {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  bgGrid->endAccess();
}
//---------------------------------------------------------------------------
void
Background::makeTrendZoneGrid(StormContGrid & grid,
                              const Surface & top_surface,
                              const Surface & base_surface,
                              const int     & nz,
                              const Simbox  * simbox) const
{
  double x_min     = simbox->GetXMin();
  double y_min     = simbox->GetYMin();
  double lx        = simbox->GetLX();
  double ly        = simbox->GetLY();
  double angle     = simbox->getAngle();
  int    nx        = simbox->getnx();
  int    ny        = simbox->getny();

  NRLib::Volume volume(x_min, y_min, lx, ly, top_surface, base_surface, angle);

  grid = StormContGrid(volume, nx, ny, nz);

  }
//---------------------------------------------------------------------------
void
Background::makeTrendZone(const float   * trend,
                          StormContGrid & trend_zone) const
{
  const size_t nx   = trend_zone.GetNI();
  const size_t ny   = trend_zone.GetNJ();
  const size_t nz   = trend_zone.GetNK();

  const double x0   = trend_zone.GetXMin();
  const double y0   = trend_zone.GetYMin();
  const double lx   = trend_zone.GetLX();
  const double ly   = trend_zone.GetLY();
  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);
  for (size_t k=0; k<nz; k++) {

    // Set trend for layer
    surface.Assign(trend[k]);

    // Set layer in background model from surface
    for(size_t j=0 ; j<ny; j++) {
      for(size_t i=0 ; i<nx; i++)
        trend_zone(i,j,k) = float(surface(i,j));
    }
  }
}

//---------------------------------------------------------------------------
void
Background::makeKrigedZone(const std::vector<KrigingData2D> & krigingData,
                           const float                      * trend,
                           StormContGrid                    & kriged_zone,
                           const CovGrid2D                  & covGrid2D) const
{
  const size_t nx   = kriged_zone.GetNI();
  const size_t ny   = kriged_zone.GetNJ();
  const size_t nz   = kriged_zone.GetNK();

  const double x0   = kriged_zone.GetXMin();
  const double y0   = kriged_zone.GetYMin();
  const double lx   = kriged_zone.GetLX();
  const double ly   = kriged_zone.GetLY();
  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);
  for (size_t k=0; k<nz; k++) {

    // Set trend for layer
    surface.Assign(trend[k]);

    // Kriging of layer
    Kriging2D::krigSurface(surface, krigingData[k], covGrid2D);

    // Set layer in background model from surface
    for(size_t j=0 ; j<ny; j++) {
      for(size_t i=0 ; i<nx; i++)
        kriged_zone(i,j,k) = float(surface(i,j));
    }
  }
}

//-------------------------------------------------------------------------------
void
Background::calculateVerticalTrend(std::vector<float *> wellTrend,
                                   float              * trend,
                                   float                logMin,
                                   float                logMax,
                                   float                maxHz,
                                   int                  nz,
                                   float                dz,
                                   const std::string  & name)
{
  int     nWells       = static_cast<int>(wellTrend.size());
  float * filtered_log = new float[nz];
  int   * count        = new int[nz];
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
    if (wellTrend[w] != NULL) {
      float * well_trend = wellTrend[w];

      for (int k = 0 ; k < nz ; k++) {
        if (well_trend[k] != RMISSING) {
          trend[k] += exp(well_trend[k]);
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

  //Utils::writeVectorToFile(std::string("trend_mean_values_") + name, trend, nz);

  smoothTrendWithLocalLinearRegression(trend, count,
                                       iWells, nz, dz,
                                       logMin,
                                       logMax,
                                       name);

  //Utils::writeVectorToFile(std::string("trend_after_linreg_") + name, trend, nz);

  WellData::applyFilter(filtered_log,
                        trend,
                        nz,
                        dz,
                        maxHz);

  for (int i=0 ; i<nz ; i++) {
    trend[i] = filtered_log[i];
  }

  //Utils::writeVectorToFile(std::string("trend_after_filter_") + name, trend, nz);


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
  // fraction, typically larger than one.
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

  float fraction   = 5.0f;                      // Require minimum 5*iWells
  int   nTimeLimit = static_cast<int>(50.0/dz); // The smaller sampling density, the more values are needed.
  int   nLowLimit  = 10;                        // Require minimum 10
  int   nDataMin   = std::max(nLowLimit, std::max(nTimeLimit, int(fraction * iWells)));

  bool  use_weights = true;
  bool  errorMid    = false;
  bool  errorHead   = false;
  bool  errorTrail  = false;

  //
  // Copy the average values (stored in array 'trend') to the array 'mean'.
  //
  float * mean = new float[nz];
  for (int k = 0 ; k < nz ; k++) {
    mean[k] = trend[k];
  }

  //
  // Find first non-missing value
  //
  int firstNonmissing = 0;
  for (int k = 0 ; k < nz ; k++) {
    if (trend[k] > 0.0f) {
      firstNonmissing = k;
      break;
    }
  }

  //
  // Find last non-missing value
  //
  int lastNonmissing = nz - 1;
  for (int k = nz - 1 ; k > 0 ; k--) {
    if (trend[k] > 0.0f) {
      lastNonmissing = k;
      break;
    }
  }

  float * x = new float[nz];  // Time indices
  float * y = new float[nz];  // Log values
  float * w = new float[nz];  // Weights (number of data behind each avg.)

  for (int k = 0 ; k < nz ; k++) {
    int nCurDataMin = nDataMin;
    if (k < firstNonmissing || k > lastNonmissing) {
      nCurDataMin *= 2;
    }

    int n = 0;
    int nData = 0;
    if (debug)
      LogKit::LogFormatted(LogKit::Low,"k=%d\n",k);
    //
    // 1. Add current data point to local data set if present.
    //
    if (count[k] > 0) {
      w[0]   = static_cast<float>(count[k]);
      x[0]   = static_cast<float>(k);
      y[0]   = trend[k];
      nData += count[k];
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   A:t=%.2f   x[0] y[0]  %d   %.2f\n",dz*(x[0] + 0.5f),int(x[0]),y[0]);
      n++;
    }

    //
    // 2. Add local data points to get 'nCurDataMin' points behind each trend
    //    value. Note that the bandwidth varies
    //
    int i = 0;
    while (nData < nCurDataMin) {
      i++;
      if (k - i >= 0 && count[k - i] > 0) {
        w[n]   = static_cast<float>(count[k - i]);
        x[n]   = static_cast<float>(k - i);
        y[n]   = mean [k - i];
        nData += count[k - i];
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   B:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k + i < nz  && count[k + i] > 0) {
        w[n]   = static_cast<float>(count[k + i]);
        x[n]   = static_cast<float>(k + i);
        y[n]   = mean [k + i];
        nData += count[k + i];
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   C:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
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
        LogKit::LogFormatted(LogKit::Low,"Sx, Sy, Sxx, Sxy  : %.4f, %.4f, %.4f, %.4f     a, b : %.4f %.4f\n",Sx,Sy,Sxx,Sxy,a,b);

      //
      // Estimate value of regression line at requested point.
      //
      float value = a + b*static_cast<float>(k);
      if (value < min_value || value > max_value) {
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f\n",value);
        if (k < firstNonmissing)
          errorHead = true;
        else if (k > lastNonmissing)
          errorTrail = true;
        else {
          errorMid   = true;
          break;
        }
      }
      trend[k] = log(value);
    }
    else {
      trend[k] = log(y[0]);
    }
    if (debug)
      LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f        (minLog/maxLog = %.2f / %.2f)\n",exp(trend[k]),min_value,max_value);
  }

  if (errorMid) {
    // Big problem ...
    LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+parName+" using local linear\n");
    LogKit::LogFormatted(LogKit::Low,"          regression failed - trying global mean instead. Possible causes: \n");
    LogKit::LogFormatted(LogKit::Low,"          1) Available logs cover too small a part of inversion grid giving extrapolation problems.\n");
    LogKit::LogFormatted(LogKit::Low,"          2) There are too many layers in grid compared to well logs available.\n");
    float sum = 0.0f;
    int nData = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (count[k] > 0) {
        if (use_weights) {
          sum   += mean[k]*count[k];
          nData += count[k];
          if (debug)
            LogKit::LogFormatted(LogKit::Low,"k=%d  count[k], mean[k]  nData, sum  %d  %8.3f     %d  %8.3f\n",
                                 k,count[k],mean[k],nData,sum);
        }
        else {
          sum += mean[k];
          nData += 1;
        }
      }
    }
    float global_mean = sum/nData;
    for (int k = 0 ; k < nz ; k++) {
      trend[k] = log(global_mean);
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   TREND: k = %d   trend[k] = %.2f\n",k,exp(trend[k]));
    }
    LogKit::LogFormatted(LogKit::Low,"\nGlobal mean for parameter %s = %.2f\n\n",parName.c_str(),global_mean);
  }
  else {
    if (errorHead) {
      // Fix first part of trend containing missing-values.
      float firstValue = trend[firstNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter \'"+parName+"\' using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [0-%d] where the log is undefined. The first\n",firstNonmissing-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d we be used throughout this region.\n",exp(firstValue),firstNonmissing);
      for (int k = 0 ; k < firstNonmissing ; k++) {
        trend[k] = firstValue;
      }
    }
    if (errorTrail) {
      // Fix last part of trend containing missing-values.
      float lastValue = trend[lastNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+parName+" using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [%d,%d] where the log is undefined. The last\n",lastNonmissing+1,nz-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d we be used throughout this region.\n",exp(lastValue),lastNonmissing);
      for (int k = lastNonmissing + 1 ; k < nz ; k++) {
        trend[k] = lastValue;
      }
    }
  }
  delete [] x;
  delete [] y;
  delete [] w;
  delete [] mean;
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
  file.close();
}

//-------------------------------------------------------------------------------
void
Background::calculateDeviationFromVerticalTrend(std::vector<float *>   wellTrend,
                                                const float          * globalTrend,
                                                float                * avg_dev,
                                                const int              nz)
{
  int nWells = static_cast<int>(wellTrend.size());

  for (int w = 0 ; w < nWells ; w++) {
    if(wellTrend[w] != NULL) {
      float * well_trend = wellTrend[w];
      float sum_dev = 0.0f;
      int count = 0;
      for (int k = 0 ; k < nz ; k++) {
        if (well_trend[k] != RMISSING) {
          float diff = exp(well_trend[k]) - exp(globalTrend[k]);
          sum_dev += diff*diff;
          count++;
        }
      }
      if (count > 0)
        sum_dev /= count;
      avg_dev[w] = sqrt(sum_dev);
    }
    else
      avg_dev[w] = RMISSING;
  }
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
  for(int i=0;i<nWells;i++)
    index[i] = i;

  for (int i=0 ; i<nWells ; i++)
  {
    for (int j=i; j<nWells ; j++)
    {
      if (rel_avg_dev[index[j]] > rel_avg_dev[index[i]])
      {
        int tmp = index[i];
        index[i] = index[j];
        index[j] = tmp;
      }
    }
  }
  //
  // Print results
  //
  if (nWells > 0)
  {
    LogKit::LogFormatted(LogKit::Low,"\nSummary of average deviation from vertical trend (well with largest misfit listed first):\n\n");
    LogKit::LogFormatted(LogKit::Low,"Well                        Vp       Vs      Rho\n");
    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------\n");
  }
  for (int i=0 ; i<nWells ; i++)
  {
    int ii = index[i];
    if(avg_dev_alpha[ii] != RMISSING) {
      LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f    %5.1f    %5.3f\n", wells[ii]->getWellname().c_str(),
                           avg_dev_alpha[ii], avg_dev_beta[ii], avg_dev_rho[ii]);
    }
  }

  if (nWells == 1)
  {
    LogKit::LogFormatted(LogKit::High,"\nNOTE: A deviation may be observed even with one well since the global trend is");
    LogKit::LogFormatted(LogKit::High,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
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
  bool isFile = modelSettings->getFileGrid();

  if((modelSettings->getOutputGridsOther() & IO::EXTRA_GRIDS) > 0) {
    std::string fileName1 = IO::PrefixBackground() + "Vp_BackgroundGrid";
    std::string fileName2 = IO::PrefixBackground() + "Vs_BackgroundGrid";
    std::string fileName3 = IO::PrefixBackground() + "Rho_BackgroundGrid";

    FFTGrid * expAlpha = copyFFTGrid(bgAlpha, true, isFile);
    expAlpha->writeFile(fileName1, IO::PathToBackground(), timeBGSimbox);
    delete expAlpha;

    FFTGrid * expBeta = copyFFTGrid(bgBeta, true, isFile);
    expBeta->writeFile(fileName2, IO::PathToBackground(), timeBGSimbox);
    delete expBeta;

    FFTGrid * expRho = copyFFTGrid(bgRho, true, isFile);
    expRho->writeFile(fileName3, IO::PathToBackground(), timeBGSimbox);
    delete expRho;
  }

  FFTGrid * resBgAlpha = NULL;
  FFTGrid * resBgBeta = NULL;
  FFTGrid * resBgRho = NULL;

  LogKit::LogFormatted(LogKit::Low,"\nResampling background model...\n");
  resampleParameter(resBgAlpha,bgAlpha,timeSimbox, timeBGSimbox, modelSettings->getFileGrid());
  resampleParameter(resBgBeta ,bgBeta ,timeSimbox, timeBGSimbox, modelSettings->getFileGrid());
  resampleParameter(resBgRho  ,bgRho  ,timeSimbox, timeBGSimbox, modelSettings->getFileGrid());

  if((modelSettings->getOutputGridsOther() & IO::EXTRA_GRIDS) > 0) {
    std::string fileName1 = IO::PrefixBackground() + "Vp_InversionGrid";
    std::string fileName2 = IO::PrefixBackground() + "Vs_InversionGrid";
    std::string fileName3 = IO::PrefixBackground() + "Rho_InversionGrid";

    FFTGrid * expResAlpha = copyFFTGrid(resBgAlpha, true, isFile);
    expResAlpha->writeFile(fileName1, IO::PathToBackground(), timeSimbox);
    delete expResAlpha;

    FFTGrid * expResBeta = copyFFTGrid(resBgBeta, true, isFile);
    expResBeta->writeFile(fileName2, IO::PathToBackground(), timeSimbox);
    delete expResBeta;

    FFTGrid * expResRho = copyFFTGrid(resBgRho, true, isFile);
    expResRho->writeFile(fileName3, IO::PathToBackground(), timeSimbox);
    delete expResRho;
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
                              Simbox   * simboxOld,
                              bool       isFile)
{
  int nx  = simboxNew->getnx();
  int ny  = simboxNew->getny();
  int nz  = simboxNew->getnz();
  //
  // Use same padding as for nonresampled cubes
  //
  // NBNB-PAL: These grids are unpadded, so all nxp, nyp, ... would probably
  //           better be replaced by nx, ny, ... to avoid confusion...
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
  pNew = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  pNew->createRealGrid();
  pNew->setType(FFTGrid::PARAMETER);
  pNew->setAccessMode(FFTGrid::WRITE);

  pOld->setAccessMode(FFTGrid::RANDOMACCESS);

  int rnxp = 2*(nxp/2 + 1);

  double * layer = new double[nx*ny];

  for(int k=0 ; k<nzp ; k++) {
    //
    // Map a layer
    //
    int ij=0;
    for (int j=0 ; j<nyp ; j++) {
      for (int i=0 ; i<rnxp ; i++) {
        if (i < nx && j < ny && k < nz) {
          int kOld = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]);
          layer[ij] = pOld->getRealValue(i, j, kOld);
          ij++;
        }
      }
    }
    //
    // Smooth the layer (equal weighting of all neighbouring cells)
    //
    float value;
    for (int j=0 ; j<nyp ; j++) {
      for (int i=0 ; i<rnxp ; i++) {
        if (i < nx && j < ny && k < nz) {
          int n = 1;
          double sum = layer[j*nx + i];
          if (i>1) {
            sum += layer[j*nx + i - 1];
            n++;
          }
          if (j>1) {
            sum += layer[(j - 1)*nx + i];
            n++;
          }
          if (i>1 && j>1) {
            sum += layer[(j - 1)*nx + i - 1];
            n++;
          }
          if (i<nx-1) {
            sum += layer[j*nx + i + 1];
            n++;
          }
          if (j<ny-1) {
            sum += layer[(j + 1)*nx + i];
            n++;
          }
          if (i<nx-1 && j<ny-1) {
            sum += layer[(j + 1)*nx + i + 1];
            n++;
          }
          value = static_cast<float>(sum)/static_cast<float>(n);
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

  delete [] layer;
  delete [] a;
  delete [] b;
}

//-------------------------------------------------------------------------------
void
Background::padAndSetBackgroundModel(FFTGrid * bgAlpha,
                                     FFTGrid * bgBeta,
                                     FFTGrid * bgRho)
{
  //LogKit::LogFormatted(LogKit::Low,"\nPadding background model...\n");
  createPaddedParameter(backModel_[0], bgAlpha);
  createPaddedParameter(backModel_[1], bgBeta);
  createPaddedParameter(backModel_[2], bgRho);
}

//-------------------------------------------------------------------------------
void
Background::createPaddedParameter(FFTGrid *& pNew,     // Padded
                                  FFTGrid  * pOld)     // Non-padded
{
  //
  // Fill padding using linear interpolation between edges.
  //
  // When we fill the z-padding, we assume that x- and y-padding is
  // already filled. The loop structure ensures this. Likewise, it
  // is assumed that the x-padding is filled when we fill the
  // y-padding.
  //
  // The linear algortihm is not "perfect", but should be more
  // than good enough for padding the smooth background model.
  //
  int nx   = pNew->getNx();
  int ny   = pNew->getNy();
  int nz   = pNew->getNz();
  int nxp  = pNew->getNxp();
  int nyp  = pNew->getNyp();
  int nzp  = pNew->getNzp();
  int rnxp = pNew->getRNxp();

  pNew->createRealGrid();
  pNew->setType(FFTGrid::PARAMETER);

  pNew->setAccessMode(FFTGrid::RANDOMACCESS);
  pOld->setAccessMode(FFTGrid::RANDOMACCESS);

  float sum_c = 1.0f/static_cast<float>(nzp - nz + 1);
  float sum_b = 1.0f/static_cast<float>(nyp - ny + 1);
  float sum_a = 1.0f/static_cast<float>(nxp - nx + 1);

  for(int k = 0 ; k < nzp ; k++) {
    for(int j = 0 ; j < nyp ; j++) {
      for(int i = 0 ; i < rnxp ; i++) { // Must fill entire grid to avoid UMR.

        float value = RMISSING;
        if(i < nx && j < ny && k < nz) { // Not in padding
          value = pOld->getRealValue(i, j, k);
        }
        else {
          if(i >= nxp)       //In dummy area for real grid, but fill to avoid UMR.
            value = 0;
          else if(k >= nz) { // In z-padding (x- and y- padding is filled in pNew)
            float c1 = pNew->getRealValue(i, j, 0     , true);
            float c2 = pNew->getRealValue(i, j, nz - 1, true);
            float w1 = sum_c*static_cast<float>(k - nz + 1);
            float w2 = sum_c*static_cast<float>(nzp - k);
            value = c1*w1 + c2*w2;
          }
          else if(j >= ny) { // In y-padding (x-padding is filled in pNew)
            float b1 = pNew->getRealValue(i, 0     , k, true);
            float b2 = pNew->getRealValue(i, ny - 1, k, true);
            float w1 = sum_b*static_cast<float>(j - ny + 1);
            float w2 = sum_b*static_cast<float>(nyp - j);
            value = b1*w1 + b2*w2;
          }
          else if(i >= nx) { // In x-padding
            float a1 = pNew->getRealValue(     0, j, k, true);
            float a2 = pNew->getRealValue(nx - 1, j, k, true);
            float w1 = sum_a*static_cast<float>(i - nx + 1);
            float w2 = sum_a*static_cast<float>(nxp - i);
            value = a1*w1 + a2*w2;
          }
        }
        pNew->setRealValue(i,j,k,value,true);
      }
    }
  }
  pNew->endAccess();
  pOld->endAccess();
}

//-------------------------------------------------------------------------------
void
Background::writeBackgrounds(Simbox                  * simbox,
                             GridMapping             * depthMapping,
                             GridMapping             * timeMapping,
                             const bool                isFile,
                             const TraceHeaderFormat & thf) const
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

  FFTGrid * expAlpha = copyFFTGrid(backModel_[0], true, isFile);
  expAlpha->writeFile(fileName1, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expAlpha;

  FFTGrid * expBeta = copyFFTGrid(backModel_[1], true, isFile);
  expBeta->writeFile(fileName2, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expBeta;

  FFTGrid * expRho = copyFFTGrid(backModel_[2], true, isFile);
  expRho->writeFile(fileName3, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expRho;

  //
  // For debugging: write cubes not in ASCII, with padding, and with flat top.
  //
  //backModel_[0]->writeStormFile(fileName1, simbox, true, false, true, true);
  //backModel_[1]->writeStormFile(fileName2, simbox, true, false, true, true);
  //backModel_[2]->writeStormFile(fileName3, simbox, true, false, true, true);
}

FFTGrid *
Background::copyFFTGrid(FFTGrid   * origGrid,
                        const bool  expTrans,
                        const bool  fileGrid) const
{
  FFTGrid * newGrid;

  if (fileGrid)
    newGrid = new FFTFileGrid(static_cast<FFTFileGrid *>(origGrid), expTrans);
  else
    newGrid = new FFTGrid(origGrid, expTrans);

  return (newGrid);
}

//---------------------------------------------------------------------------

void
Background::ErodeSurface(Surface       *& surface,
                         const Surface *  priority_surface,
                         const Simbox  *  simbox,
                         const bool    &  compare_upward) const
{
  int nx    = simbox->getnx();
  int ny    = simbox->getny();
  double x0 = simbox->GetXMin();
  double y0 = simbox->GetYMin();
  double lx = simbox->GetLX();
  double ly = simbox->GetLY();

  NRLib::Grid2D<double> eroded_surface(nx,ny,0);
  double x;
  double y;
  double z;
  double z_priority;

  double missing = surface->GetMissingValue();
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      simbox->getXYCoord(i,j,x,y);

      z_priority = priority_surface->GetZ(x,y);
      z          = surface->GetZ(x,y);

      if(compare_upward) {
        if(z < z_priority && z != missing)
          eroded_surface(i,j) = z_priority;
        else
          eroded_surface(i,j) = z;
      }

      else {
        if(z > z_priority && z_priority != missing)
          eroded_surface(i,j) = z_priority;
        else
          eroded_surface(i,j) = z;
      }
    }
  }
  delete surface;

  surface = new Surface(x0, y0, lx, ly, eroded_surface);
}
