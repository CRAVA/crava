/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

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
//#include "nrlib/well/well.hpp"
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

//Background::Background(FFTGrid                       ** grids,
//                       const std::vector<WellData *>  & wells,
//                       FFTGrid                       *& velocity,
//                       const Simbox                   * timeSimbox,
//                       const Simbox                   * timeBGSimbox,
//                       const ModelSettings            * modelSettings)
//  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
//    vsvp_(RMISSING)
//{
//  for(int i=0 ; i<3 ; i++)
//    backModel_[i] = grids[i];
//
//  FFTGrid * bgAlpha;
//  FFTGrid * bgBeta;
//  FFTGrid * bgRho;
//
//  if (timeBGSimbox == NULL)
//  {
//    generateBackgroundModel(bgAlpha,bgBeta,bgRho,
//                            velocity,wells,
//                            timeSimbox,
//                            modelSettings);
//  }
//  else
//  {
//    generateBackgroundModel(bgAlpha,bgBeta,bgRho,
//                            velocity,wells,
//                            timeBGSimbox,
//                            modelSettings);
//    resampleBackgroundModel(bgAlpha,bgBeta,bgRho,
//                            timeBGSimbox,
//                            timeSimbox,
//                            modelSettings);
//  }
//  padAndSetBackgroundModel(bgAlpha,bgBeta,bgRho);
//
//  delete bgAlpha;
//  delete bgBeta;
//  delete bgRho;
//
//  findMeanVsVp(backModel_[0],
//               backModel_[1]);
//}

//-------------------------------------------------------------------------------

Background::Background(std::vector<NRLib::Grid<double> >          & parameters,
                       const std::vector<NRLib::Well>             & wells,
                       NRLib::Grid<double>                        & velocity,
                       const Simbox                               * time_simbox,
                       const Simbox                               * time_bg_simbox,
                       std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                       std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                       const ModelSettings                        * model_settings)
  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
    vsvp_(RMISSING)
{
  const int nx    = time_simbox->getnx();
  const int ny    = time_simbox->getny();
  const int nz    = time_simbox->getnz();

  for(int i=0 ; i<3 ; i++)
    parameters[i].Resize(nx, ny, nz);

  NRLib::Grid<double> & bg_vp = parameters[0];
  NRLib::Grid<double> & bg_vs = parameters[1];
  NRLib::Grid<double> & bg_rho = parameters[2];

  if (time_bg_simbox == NULL)
  {
    generateBackgroundModel(bg_vp, bg_vs, bg_rho,
                            velocity, wells,
                            time_simbox,
                            blocked_logs,
                            bg_blocked_logs,
                            model_settings);
  }
  else
  {
    generateBackgroundModel(bg_vp, bg_vs, bg_rho,
                            velocity, wells,
                            time_bg_simbox,
                            blocked_logs,
                            bg_blocked_logs,
                            model_settings);

    resampleBackgroundModel(bg_vp, bg_vs, bg_rho,
                            time_bg_simbox,
                            time_simbox,
                            model_settings);
  }

  //padAndSetBackgroundModel(bg_vp, bg_vs, bg_rho);

  //findMeanVsVp(back_model_[0],  //Vp_vs from background not used in CommonData
  //             back_model_[1]);
}

//-------------------------------------------------------------------------------

//Background::Background(FFTGrid                       ** grids,
//                       const std::vector<WellData *>  & wells,
//                       const Simbox                   * timeSimbox,
//                       const ModelSettings            * modelSettings,
//                       const std::vector<std::string> & surface_files)
//  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
//    vsvp_(RMISSING)
//{
//  for(int i=0 ; i<3 ; i++)
//    backModel_[i] = grids[i];
//
//  FFTGrid * bgAlpha;
//  FFTGrid * bgBeta;
//  FFTGrid * bgRho;
//
//  generateMultizoneBackgroundModel(bgAlpha,bgBeta,bgRho,
//                                   wells,
//                                   timeSimbox,
//                                   modelSettings,
//                                   surface_files);
//
//  padAndSetBackgroundModel(bgAlpha,bgBeta,bgRho);
//
//  delete bgAlpha;
//  delete bgBeta;
//  delete bgRho;
//
//  findMeanVsVp(backModel_[0],
//               backModel_[1]);
//}


//-------------------------------------------------------------------------------
Background::Background(std::vector<NRLib::Grid<double> > & parameters,
                       const std::vector<NRLib::Well>    & wells,
                       const Simbox                      * interval_simbox,
                       const ModelSettings               * model_settings,
                       const std::vector<std::string>    & surface_files,
                       std::string                       & err_text)
  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
    vsvp_(RMISSING)
{

  const int nx    = interval_simbox->getnx();
  const int ny    = interval_simbox->getny();
  const int nz    = interval_simbox->getnz();

  for(int i=0 ; i<3 ; i++)
    parameters[i].Resize(nx, ny, nz);

  NRLib::Grid<double> & bg_vp = parameters[0];
  NRLib::Grid<double> & bg_vs = parameters[1];
  NRLib::Grid<double> & bg_rho = parameters[2];

  GenerateMultizoneBackgroundModel(bg_vp, bg_vs, bg_rho,
                                   wells,
                                   interval_simbox,
                                   model_settings,
                                   surface_files,
                                   err_text);

  //padAndSetBackgroundModel(bg_vp, bg_vs, bg_rho);

  //findMeanVsVp(back_model_[0],
  //             back_model_[1]);
}

//-------------------------------------------------------------------------------
Background::Background(std::vector<std::vector<NRLib::Grid<double> > > & parameters,
                       const std::vector<NRLib::Well>                  & wells,
                       MultiIntervalGrid                               * multiple_interval_grid,
                       const ModelSettings                             * model_settings,
                       std::string                                     & err_text)
  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
    vsvp_(RMISSING)
{

  //Background for multiple intervals
  int n_intervals = multiple_interval_grid->GetNIntervals();

  std::vector<NRLib::Grid<double> > bg_vp(n_intervals);
  std::vector<NRLib::Grid<double> > bg_vs(n_intervals);
  std::vector<NRLib::Grid<double> > bg_rho(n_intervals);

  GenerateMultiIntervalBackgroundModel(bg_vp, bg_vs, bg_rho,
                                       wells,
                                       multiple_interval_grid,
                                       model_settings,
                                       err_text);

  parameters.push_back(bg_vp);
  parameters.push_back(bg_vs);
  parameters.push_back(bg_rho);

  //padAndSetBackgroundModelInterval(bg_vp, bg_vs, bg_rho);

  //findMeanVsVp(back_model_[0], //Vp_vs from background not used in CommonData
               //back_model_[1]);
}

//-------------------------------------------------------------------------------
Background::Background(FFTGrid ** grids)
  : DataTarget_(IMISSING),
    vsvp_(RMISSING)
{
  for(int i=0 ; i<3 ; i++)
    back_model_[i] = grids[i];
  findMeanVsVp(back_model_[0],
               back_model_[1]);
}

//-------------------------------------------------------------------------------
Background::~Background(void)
{
  for (int i=0 ; i<3 ; i++)
    delete back_model_[i];
}

//-------------------------------------------------------------------------------
void
Background::releaseGrids()
{
  for (int i=0 ; i<3 ; i++)
    back_model_[i] = NULL;
}

//-------------------------------------------------------------------------------
//void
//Background::generateBackgroundModel(FFTGrid                      *& bgAlpha,
//                                    FFTGrid                      *& bgBeta,
//                                    FFTGrid                      *& bgRho,
//                                    FFTGrid                      *& velocity,
//                                    const std::vector<WellData *> & wells,
//                                    const Simbox                  * simbox,
//                                    const ModelSettings           * modelSettings)
//{
//  const int   nz     = simbox->getnz();
//  const int   nWells = modelSettings->getNumberOfWells();
//  const float dz     = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());
//
//  std::string name_vp  = "Vp";
//  std::string name_vs  = "Vs";
//  std::string name_rho = "Rho";
//
//  std::vector<float *> wellTrendAlpha(nWells);
//  std::vector<float *> wellTrendBeta(nWells);
//  std::vector<float *> wellTrendRho(nWells);
//  std::vector<float *> highCutWellTrendAlpha(nWells);
//  std::vector<float *> highCutWellTrendBeta(nWells);
//  std::vector<float *> highCutWellTrendRho(nWells);
//
//  getWellTrends(wellTrendAlpha,highCutWellTrendAlpha,wells,nz,name_vp);
//  getWellTrends(wellTrendBeta, highCutWellTrendBeta, wells,nz,name_vs);
//  getWellTrends(wellTrendRho,  highCutWellTrendRho,  wells,nz,name_rho);
//
//  float * trendAlpha = new float[nz];
//  float * trendBeta  = new float[nz];
//  float * trendRho   = new float[nz];
//  float * trendVel   = new float[nz]; // Allocate (for simplicity) although not always needed
//
//  float * avgDevAlpha = new float[nWells];
//  float * avgDevBeta  = new float[nWells];
//  float * avgDevRho   = new float[nWells];
//  float * avgDevVel   = new float[nWells]; // Allocate (for simplicity) although not always needed
//
//  calculateBackgroundTrend(trendAlpha,
//                           avgDevAlpha,
//                           nz,
//                           dz,
//                           modelSettings->getAlphaMin(),
//                           modelSettings->getAlphaMax(),
//                           modelSettings->getMaxHzBackground(),
//                           wellTrendAlpha,
//                           highCutWellTrendAlpha,
//                           name_vp);
//  calculateBackgroundTrend(trendBeta,
//                           avgDevBeta,
//                           nz,
//                           dz,
//                           modelSettings->getBetaMin(),
//                           modelSettings->getBetaMax(),
//                           modelSettings->getMaxHzBackground(),
//                           wellTrendBeta,
//                           highCutWellTrendBeta,
//                           name_vs);
//  calculateBackgroundTrend(trendRho,
//                           avgDevRho,
//                           nz,
//                           dz,
//                           modelSettings->getRhoMin(),
//                           modelSettings->getRhoMax(),
//                           modelSettings->getMaxHzBackground(),
//                           wellTrendRho,
//                           highCutWellTrendRho,
//                           name_rho);
//
//  bool hasVelocityTrend = velocity != NULL;
//  bool write1D          = ((modelSettings->getOtherOutputFlag()& IO::BACKGROUND_TREND_1D) > 0);
//  bool write3D          = ((modelSettings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);
//
//  writeTrendsToFile(trendAlpha,simbox, write1D, write3D, hasVelocityTrend, name_vp, modelSettings->getFileGrid());
//  writeTrendsToFile(trendBeta, simbox, write1D, write3D, hasVelocityTrend, name_vs, modelSettings->getFileGrid());
//  writeTrendsToFile(trendRho,  simbox, write1D, write3D, hasVelocityTrend, name_rho,modelSettings->getFileGrid());
//
//  if (velocity != NULL) {
//    //
//    // We still want calculateBackgroundTrend() for alpha above. By calculating
//    // avgDevAlpha we can check that the bgAlpha calculated from velocity is as
//    // good as or better than that calculated by crava.
//    //
//    calculateVelocityDeviations(velocity, wells, simbox,
//                                trendVel, avgDevVel, avgDevAlpha,
//                                modelSettings->getOutputGridsElastic(),
//                                nWells);
//    velocity->logTransf();
//    delete bgAlpha;
//    bgAlpha = velocity;
//    velocity = NULL;
//    writeDeviationsFromVerticalTrend(avgDevVel, avgDevBeta, avgDevRho,
//                                     trendVel, trendBeta, trendRho,
//                                     wells, nWells, nz);
//  }
//  else {
//    writeDeviationsFromVerticalTrend(avgDevAlpha, avgDevBeta, avgDevRho,
//                                     trendAlpha, trendBeta, trendRho,
//                                     wells, nWells, nz);
//  }
//
//  std::vector<KrigingData2D> krigingDataAlpha(nz);
//  std::vector<KrigingData2D> krigingDataBeta(nz);
//  std::vector<KrigingData2D> krigingDataRho(nz);
//
//  std::vector<float *>     blAlpha(nWells);   // bl = blocked logs
//  std::vector<float *>     blBeta(nWells);
//  std::vector<float *>     blRho(nWells);
//  std::vector<float *>     vtAlpha(nWells);   // vt = vertical trend
//  std::vector<float *>     vtBeta(nWells);
//  std::vector<float *>     vtRho(nWells);
//  std::vector<const int *> ipos(nWells);
//  std::vector<const int *> jpos(nWells);
//  std::vector<const int *> kpos(nWells);
//
//  for(int i=0; i<nWells; i++) {
//    vtAlpha[i] = new float[nz];
//    vtBeta[i]  = new float[nz];
//    vtRho[i]   = new float[nz];
//  }
//
//  std::vector<int> nBlocks(nWells);
//  int              totBlocks;
//
//  getKrigingWellTrends(blAlpha,blBeta,blRho,
//                       vtAlpha,vtBeta,vtRho,
//                       ipos,jpos,kpos,
//                       nBlocks,totBlocks,
//                       wells,nWells);
//
//  setupKrigingData2D(krigingDataAlpha,krigingDataBeta,krigingDataRho,
//                     trendAlpha,trendBeta,trendRho,
//                     modelSettings->getOutputGridsElastic(),
//                     nz,dz,totBlocks,nBlocks,
//                     blAlpha,blBeta,blRho,
//                     vtAlpha,vtBeta,vtRho,
//                     ipos,jpos,kpos);
//
//  const CovGrid2D & covGrid2D = makeCovGrid2D(simbox,
//                                              modelSettings->getBackgroundVario(),
//                                              modelSettings->getDebugFlag());
//
//  makeKrigedBackground(krigingDataAlpha, bgAlpha, trendAlpha, simbox, covGrid2D, "Vp" , modelSettings->getFileGrid());
//  makeKrigedBackground(krigingDataBeta , bgBeta , trendBeta , simbox, covGrid2D, "Vs" , modelSettings->getFileGrid());
//  makeKrigedBackground(krigingDataRho  , bgRho  , trendRho  , simbox, covGrid2D, "Rho", modelSettings->getFileGrid());
//
//  delete &covGrid2D;
//
//  delete [] avgDevAlpha;
//  delete [] avgDevBeta;
//  delete [] avgDevRho;
//  delete [] avgDevVel;
//
//  delete [] trendAlpha;
//  delete [] trendBeta;
//  delete [] trendRho;
//  delete [] trendVel;
//
//  for(int i=0; i<nWells; i++) {
//
//    delete [] wellTrendAlpha[i];
//    delete [] wellTrendBeta[i];
//    delete [] wellTrendRho[i];
//
//    delete [] highCutWellTrendAlpha[i];
//    delete [] highCutWellTrendBeta[i];
//    delete [] highCutWellTrendRho[i];
//
//    delete [] blAlpha[i];
//    delete [] blBeta[i];
//    delete [] blRho[i];
//
//    delete [] vtAlpha[i];
//    delete [] vtBeta[i];
//    delete [] vtRho[i];
//
//  }
//}

//-------------------------------------------------------------------------------
void
Background::generateBackgroundModel(NRLib::Grid<double>                        & bg_vp,
                                    NRLib::Grid<double>                        & bg_vs,
                                    NRLib::Grid<double>                        & bg_rho,
                                    NRLib::Grid<double>                        & velocity,
                                    const std::vector<NRLib::Well>             & wells,
                                    const Simbox                               * simbox,
                                    std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                    std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                                    const ModelSettings                        * model_settings)
{
  const int   nz     = simbox->getnz();
  const int   n_wells = model_settings->getNumberOfWells();
  const float dz     = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  std::string name_vp  = "Vp";
  std::string name_vs  = "Vs";
  std::string name_rho = "Rho";

  std::vector<float *> well_trend_vp(n_wells);
  std::vector<float *> well_trend_vs(n_wells);
  std::vector<float *> well_trend_rho(n_wells);
  std::vector<float *> high_cut_well_trend_vp(n_wells);
  std::vector<float *> high_cut_well_trend_vs(n_wells);
  std::vector<float *> high_cut_well_trend_rho(n_wells);

  getWellTrends(well_trend_vp, high_cut_well_trend_vp, wells, bg_blocked_logs, nz, name_vp);
  getWellTrends(well_trend_vs, high_cut_well_trend_vs, wells, bg_blocked_logs, nz, name_vs);
  getWellTrends(well_trend_rho, high_cut_well_trend_rho, wells, bg_blocked_logs, nz, name_rho);

  float * trend_vp = new float[nz];
  float * trend_vs = new float[nz];
  float * trend_rho  = new float[nz];
  float * trend_vel  = new float[nz]; // Allocate (for simplicity) although not always needed

  float * avg_dev_vp = new float[n_wells];
  float * avg_dev_vs = new float[n_wells];
  float * avg_dev_rho = new float[n_wells];
  float * avg_dev_vel = new float[n_wells]; // Allocate (for simplicity) although not always needed

  calculateBackgroundTrend(trend_vp,
                           avg_dev_vp,
                           nz,
                           dz,
                           model_settings->getAlphaMin(),
                           model_settings->getAlphaMax(),
                           model_settings->getMaxHzBackground(),
                           well_trend_vp,
                           high_cut_well_trend_vp,
                           name_vp);
  calculateBackgroundTrend(trend_vs,
                           avg_dev_vs,
                           nz,
                           dz,
                           model_settings->getBetaMin(),
                           model_settings->getBetaMax(),
                           model_settings->getMaxHzBackground(),
                           well_trend_vs,
                           high_cut_well_trend_vs,
                           name_vs);
  calculateBackgroundTrend(trend_rho,
                           avg_dev_rho,
                           nz,
                           dz,
                           model_settings->getRhoMin(),
                           model_settings->getRhoMax(),
                           model_settings->getMaxHzBackground(),
                           well_trend_rho,
                           high_cut_well_trend_rho,
                           name_rho);

  bool has_velocity_trend = velocity.GetN() != 0; // != NULL;
  bool write1D          = ((model_settings->getOtherOutputFlag()& IO::BACKGROUND_TREND_1D) > 0);
  bool write3D          = ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);

  writeTrendsToFile(trend_vp, simbox, write1D, write3D, has_velocity_trend, name_vp, model_settings->getFileGrid());
  writeTrendsToFile(trend_vs, simbox, write1D, write3D, has_velocity_trend, name_vs, model_settings->getFileGrid());
  writeTrendsToFile(trend_rho, simbox, write1D, write3D, has_velocity_trend, name_rho, model_settings->getFileGrid());

  if (velocity.GetN() != 0) { //!= NULL) {
    //
    // We still want calculateBackgroundTrend() for alpha above. By calculating
    // avgDevAlpha we can check that the bgAlpha calculated from velocity is as
    // good as or better than that calculated by crava.
    //
    calculateVelocityDeviations(velocity, wells, simbox, blocked_logs, bg_blocked_logs,
                                trend_vel, avg_dev_vel, avg_dev_vp,
                                model_settings->getOutputGridsElastic(),
                                n_wells);

    velocity.LogTransform(RMISSING);
    //velocity->logTransf();

    //delete bg_vp;
    bg_vp = velocity;
    velocity.Resize(0, 0, 0, 0.0); // = NULL;
    writeDeviationsFromVerticalTrend(avg_dev_vel, avg_dev_vs, avg_dev_rho,
                                     trend_vel, trend_vs, trend_rho,
                                     wells, n_wells, nz);
  }
  else {
    writeDeviationsFromVerticalTrend(avg_dev_vp, avg_dev_vs, avg_dev_rho,
                                     trend_vp, trend_vs, trend_rho,
                                     wells, n_wells, nz);
  }

  std::vector<KrigingData2D> kriging_data_vp(nz);
  std::vector<KrigingData2D> kriging_data_vs(nz);
  std::vector<KrigingData2D> kriging_data_rho(nz);

  std::vector<std::vector<double > > bl_vp(n_wells); // bl = blocked logs
  std::vector<std::vector<double > > bl_vs(n_wells);
  std::vector<std::vector<double > > bl_rho(n_wells);
  std::vector<std::vector<double > > vt_vp(n_wells);   // vt = vertical trend
  std::vector<std::vector<double > > vt_vs(n_wells);
  std::vector<std::vector<double > > vt_rho(n_wells);

  std::vector<const std::vector<int> > ipos(n_wells); //std::vector<const int *> ipos(n_wells);
  std::vector<const std::vector<int> > jpos(n_wells);
  std::vector<const std::vector<int> > kpos(n_wells);

  for(int i=0; i<n_wells; i++) {
    vt_vp[i] = std::vector<double>(nz);
    vt_vs[i] = std::vector<double>(nz);
    vt_rho[i] = std::vector<double>(nz);
  }

  std::vector<int> n_blocks(n_wells);
  int              tot_blocks;

  getKrigingWellTrends(bl_vp,bl_vs,bl_rho,
                       vt_vp,vt_vs,vt_rho,
                       ipos,jpos,kpos,
                       n_blocks,tot_blocks,
                       wells, bg_blocked_logs,
                       n_wells);

  setupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
                     trend_vp,trend_vs,trend_rho,
                     model_settings->getOutputGridsElastic(),
                     nz,dz,tot_blocks,n_blocks,
                     bl_vp,bl_vs,bl_rho,
                     vt_vp,vt_vs,vt_rho,
                     ipos,jpos,kpos);

  const CovGrid2D & covGrid2D = makeCovGrid2D(simbox,
                                              model_settings->getBackgroundVario(),
                                              model_settings->getDebugFlag());

  makeKrigedBackground(kriging_data_vp, bg_vp, trend_vp, simbox, covGrid2D, "Vp");
  makeKrigedBackground(kriging_data_vs , bg_vs , trend_vs , simbox, covGrid2D, "Vs");
  makeKrigedBackground(kriging_data_rho  , bg_rho  , trend_rho  , simbox, covGrid2D, "Rho");

  delete &covGrid2D;

  delete [] avg_dev_vp;
  delete [] avg_dev_vs;
  delete [] avg_dev_rho;
  delete [] avg_dev_vel;

  delete [] trend_vp;
  delete [] trend_vs;
  delete [] trend_rho;
  delete [] trend_vel;

  for(int i=0; i<n_wells; i++) {

    delete [] well_trend_vp[i];
    delete [] well_trend_vs[i];
    delete [] well_trend_rho[i];

    delete [] high_cut_well_trend_vp[i];
    delete [] high_cut_well_trend_vs[i];
    delete [] high_cut_well_trend_rho[i];

  }
}


//-------------------------------------------------------------------------------
//void
//Background::generateMultizoneBackgroundModel(FFTGrid                       *& bgAlpha,
//                                             FFTGrid                       *& bgBeta,
//                                             FFTGrid                       *& bgRho,
//                                             const std::vector<WellData *>  & wells,
//                                             const Simbox                   * simbox,
//                                             const ModelSettings            * modelSettings,
//                                             const std::vector<std::string> & surface_files)
//{
//  LogKit::LogFormatted(LogKit::Low,"Multizone background model:\n");
//
//  std::vector<int> correlation_structure = modelSettings->getCorrelationStructure();
//  std::vector<int> erosion_priority      = modelSettings->getErosionPriority();
//
//  int    nWells    = modelSettings->getNumberOfWells();
//  int    nZones    = static_cast<int>(correlation_structure.size()) - 1;
//  float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory
//
//  std::vector<Surface> surface(nZones+1);
//  for(int i=0; i<nZones+1; i++)
//    surface[i] = Surface(surface_files[i]);
//
//  std::vector<StormContGrid> alpha_zones(nZones);
//  std::vector<StormContGrid> beta_zones(nZones);
//  std::vector<StormContGrid> rho_zones(nZones);
//
//  BuildSeismicPropertyZones(alpha_zones,
//                            beta_zones,
//                            rho_zones,
//                            surface,
//                            correlation_structure,
//                            simbox,
//                            dz);
//
//  std::vector<Surface *> eroded_surfaces(nZones+1);
//  for(int i=0; i<nZones+1; i++)
//    eroded_surfaces[i] = NULL;
//
//  ErodeAllSurfaces(eroded_surfaces,
//                   erosion_priority,
//                   surface,
//                   simbox);
//
//  const CovGrid2D & covGrid2D = makeCovGrid2D(simbox, modelSettings->getBackgroundVario(), modelSettings->getDebugFlag());
//
//  std::string name_vp  = "Vp";
//  std::string name_vs  = "Vs";
//  std::string name_rho = "Rho";
//
//  std::vector<float *> trendAlphaZone(nZones);
//  std::vector<float *> trendBetaZone(nZones);
//  std::vector<float *> trendRhoZone(nZones);
//
//  for(int i=0; i<nZones; i++) {
//    LogKit::LogFormatted(LogKit::Low,"\nZone%2d:",i+1);
//
//    int nz = static_cast<int>(alpha_zones[i].GetNK());
//
//    std::vector<float *> wellTrendAlpha(nWells);
//    std::vector<float *> wellTrendBeta(nWells);
//    std::vector<float *> wellTrendRho(nWells);
//
//    std::vector<float *> highCutWellTrendAlpha(nWells);
//    std::vector<float *> highCutWellTrendBeta(nWells);
//    std::vector<float *> highCutWellTrendRho(nWells);
//
//    StormContGrid eroded_zone;
//
//    BuildErodedZones(eroded_zone,
//                     eroded_surfaces,
//                     nz,
//                     simbox,
//                     i);
//
//    std::vector<bool> hitZone(nWells);
//    checkWellHitsZone(hitZone, wells, eroded_zone, nWells);
//
//    std::vector<BlockedLogsForZone *> blocked_logs(nWells);
//
//    getWellTrendsZone(blocked_logs, wellTrendAlpha, highCutWellTrendAlpha, wells, eroded_zone, hitZone, nz, name_vp,  i);
//    getWellTrendsZone(blocked_logs, wellTrendBeta,  highCutWellTrendBeta,  wells, eroded_zone, hitZone, nz, name_vs,  i);
//    getWellTrendsZone(blocked_logs, wellTrendRho,   highCutWellTrendRho,   wells, eroded_zone, hitZone, nz, name_rho, i);
//
//    trendAlphaZone[i] = new float[nz];
//    trendBetaZone[i]  = new float[nz];
//    trendRhoZone[i]   = new float[nz];
//
//    float * avgDevAlphaZone = new float[nWells];
//    float * avgDevBetaZone  = new float[nWells];
//    float * avgDevRhoZone   = new float[nWells];
//
//    calculateBackgroundTrend(trendAlphaZone[i],
//                             avgDevAlphaZone,
//                             nz,
//                             dz,
//                             modelSettings->getAlphaMin(),
//                             modelSettings->getAlphaMax(),
//                             modelSettings->getMaxHzBackground(),
//                             wellTrendAlpha,
//                             highCutWellTrendAlpha,
//                             name_vp);
//    calculateBackgroundTrend(trendBetaZone[i],
//                             avgDevBetaZone,
//                             nz,
//                             dz,
//                             modelSettings->getBetaMin(),
//                             modelSettings->getBetaMax(),
//                             modelSettings->getMaxHzBackground(),
//                             wellTrendBeta,
//                             highCutWellTrendBeta,
//                             name_vs);
//    calculateBackgroundTrend(trendRhoZone[i],
//                             avgDevRhoZone,
//                             nz,
//                             dz,
//                             modelSettings->getRhoMin(),
//                             modelSettings->getRhoMax(),
//                             modelSettings->getMaxHzBackground(),
//                             wellTrendRho,
//                             highCutWellTrendRho,
//                             name_rho);
//
//    writeDeviationsFromVerticalTrend(avgDevAlphaZone,
//                                     avgDevBetaZone,
//                                     avgDevRhoZone,
//                                     trendAlphaZone[i],
//                                     trendBetaZone[i],
//                                     trendRhoZone[i],
//                                     wells,
//                                     nWells,
//                                     nz);
//
//    std::vector<float *>     blAlpha(nWells);   // bl = blocked logs
//    std::vector<float *>     blBeta(nWells);
//    std::vector<float *>     blRho(nWells);
//
//    std::vector<float *>     vtAlpha(nWells);   // vt = vertical trend
//    std::vector<float *>     vtBeta(nWells);
//    std::vector<float *>     vtRho(nWells);
//
//    std::vector<const int *> ipos(nWells);
//    std::vector<const int *> jpos(nWells);
//    std::vector<const int *> kpos(nWells);
//
//    int              totBlocks;
//    std::vector<int> nBlocks(nWells);
//
//    getKrigingWellTrendsZone(blocked_logs,
//                             blAlpha,blBeta,blRho,
//                             vtAlpha,vtBeta,vtRho,
//                             ipos,jpos,kpos,
//                             nBlocks,totBlocks,nz);
//
//    std::vector<KrigingData2D> krigingDataAlpha(nz);
//    std::vector<KrigingData2D> krigingDataBeta(nz);
//    std::vector<KrigingData2D> krigingDataRho(nz);
//
//    setupKrigingData2D(krigingDataAlpha,krigingDataBeta,krigingDataRho,
//                       trendAlphaZone[i],trendBetaZone[i],trendRhoZone[i],
//                       modelSettings->getOutputGridsElastic(),
//                       nz,dz,totBlocks,nBlocks,
//                       blAlpha,blBeta,blRho,
//                       vtAlpha,vtBeta,vtRho,
//                       ipos,jpos,kpos);
//
//    makeKrigedZone(krigingDataAlpha, trendAlphaZone[i], alpha_zones[i], covGrid2D);
//    makeKrigedZone(krigingDataBeta , trendBetaZone[i] , beta_zones[i],  covGrid2D);
//    makeKrigedZone(krigingDataRho  , trendRhoZone[i]  , rho_zones[i],   covGrid2D);
//
//    delete [] avgDevAlphaZone;
//    delete [] avgDevBetaZone;
//    delete [] avgDevRhoZone;
//
//    for(int j=0; j<nWells; j++) {
//      delete [] wellTrendAlpha[j];
//      delete [] wellTrendBeta[j];
//      delete [] wellTrendRho[j];
//
//      delete [] highCutWellTrendAlpha[j];
//      delete [] highCutWellTrendBeta[j];
//      delete [] highCutWellTrendRho[j];
//
//      delete [] blAlpha[j];
//      delete [] blBeta[j];
//      delete [] blRho[j];
//
//      delete [] vtAlpha[j];
//      delete [] vtBeta[j];
//      delete [] vtRho[j];
//
//      delete blocked_logs[j];
//    }
//  }
//
//  MakeMultizoneBackground(bgAlpha,bgBeta,bgRho,
//                          alpha_zones, beta_zones, rho_zones,
//                          simbox,
//                          erosion_priority,
//                          surface,
//                          modelSettings->getSurfaceUncertainty(),
//                          modelSettings->getFileGrid(),
//                          "multizone");
//
//
//  bool write3D = ((modelSettings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);
//
//  if(write3D) {
//    writeMultizoneTrendsToFile(trendAlphaZone, trendBetaZone, trendRhoZone,
//                               alpha_zones, beta_zones, rho_zones,
//                               simbox,
//                               erosion_priority,
//                               surface,
//                               modelSettings->getSurfaceUncertainty(),
//                               modelSettings->getFileGrid());
//
//  }
//
//  delete &covGrid2D;
//
//  for(int i=0; i<nZones; i++) {
//    delete [] trendAlphaZone[i];
//    delete [] trendBetaZone[i];
//    delete [] trendRhoZone[i];
//
//    delete eroded_surfaces[i];
//  }
//
//}

//-------------------------------------------------------------------------------
void
Background::GenerateMultizoneBackgroundModel(NRLib::Grid<double>            & bg_vp, //FFTGrid                       *& bg_vp,
                                             NRLib::Grid<double>            & bg_vs,
                                             NRLib::Grid<double>            & bg_rho,
                                             const std::vector<NRLib::Well> & wells,
                                             const Simbox                   * simbox,
                                             const ModelSettings            * model_settings,
                                             const std::vector<std::string> & surface_files,
                                             std::string                    & err_text)
{
  LogKit::LogFormatted(LogKit::Low,"Multizone background model:\n");

  std::vector<int> correlation_structure = model_settings->getCorrelationStructure();
  std::vector<int> erosion_priority      = model_settings->getErosionPriority();

  int    n_wells    = model_settings->getNumberOfWells();
  int    n_zones    = static_cast<int>(correlation_structure.size()) - 1;
  float  dz         = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

  std::vector<Surface> surface(n_zones+1);
  for(int i=0; i<n_zones+1; i++)
    surface[i] = Surface(surface_files[i]);

  std::vector<StormContGrid> vp_zones(n_zones);
  std::vector<StormContGrid> vs_zones(n_zones);
  std::vector<StormContGrid> rho_zones(n_zones);

  std::string err_text_tmp = "";
  BuildSeismicPropertyZones(vp_zones,
                            vs_zones,
                            rho_zones,
                            surface,
                            correlation_structure,
                            simbox,
                            dz,
                            err_text_tmp);

  if(err_text_tmp != "") {
    err_text += err_text_tmp;
  }
  else {

    std::vector<Surface *> eroded_surfaces(n_zones+1);
    for(int i=0; i<n_zones+1; i++)
      eroded_surfaces[i] = NULL;

    ErodeAllSurfaces(eroded_surfaces,
                     erosion_priority,
                     surface,
                     simbox);

    const CovGrid2D & covGrid2D = makeCovGrid2D(simbox, model_settings->getBackgroundVario(), model_settings->getDebugFlag());

    std::string name_vp  = "Vp";
    std::string name_vs  = "Vs";
    std::string name_rho = "Rho";

    std::vector<float *> trend_vp_zone(n_zones);
    std::vector<float *> trend_vs_zone(n_zones);
    std::vector<float *> trend_rho_zone(n_zones);

    for(int i=0; i<n_zones; i++) {
      LogKit::LogFormatted(LogKit::Low,"\nZone%2d:",i+1);

      int nz = static_cast<int>(vp_zones[i].GetNK());

      std::vector<float *> well_trend_vp(n_wells);
      std::vector<float *> well_trend_vs(n_wells);
      std::vector<float *> well_trend_rho(n_wells);

      std::vector<float *> high_cut_well_trend_vp(n_wells);
      std::vector<float *> high_cut_well_trend_vs(n_wells);
      std::vector<float *> high_cut_well_trend_rho(n_wells);

      StormContGrid eroded_zone;

      BuildErodedZones(eroded_zone,
                       eroded_surfaces,
                       nz,
                       simbox,
                       i);

      std::vector<bool> hit_zone(n_wells);
      checkWellHitsZone(hit_zone, wells, eroded_zone, n_wells);

      std::vector<BlockedLogsCommon *> blocked_logs(n_wells);

      getWellTrendsZone(model_settings, blocked_logs, well_trend_vp, high_cut_well_trend_vp, wells, eroded_zone, hit_zone, nz, name_vp,  i);
      getWellTrendsZone(model_settings, blocked_logs, well_trend_vs, high_cut_well_trend_vs, wells, eroded_zone, hit_zone, nz, name_vs,  i);
      getWellTrendsZone(model_settings, blocked_logs, well_trend_rho, high_cut_well_trend_rho, wells, eroded_zone, hit_zone, nz, name_rho, i);

      trend_vp_zone[i] = new float[nz];
      trend_vs_zone[i]  = new float[nz];
      trend_rho_zone[i]   = new float[nz];

      float * avg_dev_vp_zone = new float[n_wells];
      float * avg_dev_vs_zone = new float[n_wells];
      float * avg_dev_rho_zone = new float[n_wells];

      calculateBackgroundTrend(trend_vp_zone[i],
                               avg_dev_vp_zone,
                               nz,
                               dz,
                               model_settings->getAlphaMin(),
                               model_settings->getAlphaMax(),
                               model_settings->getMaxHzBackground(),
                               well_trend_vp,
                               high_cut_well_trend_vp,
                               name_vp);
      calculateBackgroundTrend(trend_vs_zone[i],
                               avg_dev_vs_zone,
                               nz,
                               dz,
                               model_settings->getBetaMin(),
                               model_settings->getBetaMax(),
                               model_settings->getMaxHzBackground(),
                               well_trend_vs,
                               high_cut_well_trend_vs,
                               name_vs);
      calculateBackgroundTrend(trend_rho_zone[i],
                               avg_dev_rho_zone,
                               nz,
                               dz,
                               model_settings->getRhoMin(),
                               model_settings->getRhoMax(),
                               model_settings->getMaxHzBackground(),
                               well_trend_rho,
                               high_cut_well_trend_rho,
                               name_rho);

      writeDeviationsFromVerticalTrend(avg_dev_vp_zone,
                                       avg_dev_vs_zone,
                                       avg_dev_rho_zone,
                                       trend_vp_zone[i],
                                       trend_vs_zone[i],
                                       trend_rho_zone[i],
                                       wells,
                                       n_wells,
                                       nz);

      std::vector<std::vector<double > > bl_vp(n_wells); // bl = blocked logs
      std::vector<std::vector<double > > bl_vs(n_wells);
      std::vector<std::vector<double > > bl_rho(n_wells);
      std::vector<std::vector<double > > vt_vp(n_wells);   // vt = vertical trend
      std::vector<std::vector<double > > vt_vs(n_wells);
      std::vector<std::vector<double > > vt_rho(n_wells);

      std::vector<const std::vector<int> > ipos(n_wells); //std::vector<const int *> ipos(n_wells);
      std::vector<const std::vector<int> > jpos(n_wells);
      std::vector<const std::vector<int> > kpos(n_wells);

      int              totBlocks;
      std::vector<int> nBlocks(n_wells);

      getKrigingWellTrendsZone(blocked_logs,
                               bl_vp,bl_vs,bl_rho,
                               vt_vp,vt_vs,vt_rho,
                               ipos,jpos,kpos,
                               nBlocks,totBlocks,nz);

      std::vector<KrigingData2D> kriging_data_vp(nz);
      std::vector<KrigingData2D> kriging_data_vs(nz);
      std::vector<KrigingData2D> kriging_data_rho(nz);

      setupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
                         trend_vp_zone[i],trend_vs_zone[i],trend_rho_zone[i],
                         model_settings->getOutputGridsElastic(),
                         nz,dz,totBlocks,nBlocks,
                         bl_vp,bl_vs,bl_rho,
                         vt_vp,vt_vs,vt_rho,
                         ipos,jpos,kpos);

      makeKrigedZone(kriging_data_vp, trend_vp_zone[i], vp_zones[i], covGrid2D);
      makeKrigedZone(kriging_data_vs, trend_vs_zone[i], vs_zones[i], covGrid2D);
      makeKrigedZone(kriging_data_rho, trend_rho_zone[i], rho_zones[i], covGrid2D);

      delete [] avg_dev_vp_zone;
      delete [] avg_dev_vs_zone;
      delete [] avg_dev_rho_zone;

      for(int j=0; j<n_wells; j++) {
        delete [] well_trend_vp[j];
        delete [] well_trend_vs[j];
        delete [] well_trend_rho[j];

        delete [] high_cut_well_trend_vp[j];
        delete [] high_cut_well_trend_vs[j];
        delete [] high_cut_well_trend_rho[j];

        delete blocked_logs[j];
      }
    }

    MakeMultizoneBackground(bg_vp, bg_vs, bg_rho,
                            vp_zones, vs_zones, rho_zones,
                            simbox,
                            erosion_priority,
                            surface,
                            model_settings->getSurfaceUncertainty(),
                            model_settings->getFileGrid(),
                            "multizone");


    bool write3D = ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);

    if(write3D) {
      writeMultizoneTrendsToFile(trend_vp_zone, trend_vs_zone, trend_rho_zone, //H Not working
                                 vp_zones, vs_zones, rho_zones,
                                 simbox,
                                 erosion_priority,
                                 surface,
                                 model_settings->getSurfaceUncertainty(),
                                 model_settings->getFileGrid());

    }

    delete &covGrid2D;

    for(int i=0; i<n_zones; i++) {
      delete [] trend_vp_zone[i];
      delete [] trend_vs_zone[i];
      delete [] trend_rho_zone[i];

      delete eroded_surfaces[i];
    }
  }

}


//-------------------------------------------------------------------------------
void
Background::GenerateMultiIntervalBackgroundModel(std::vector<NRLib::Grid<double> > & bg_vp,
                                                 std::vector<NRLib::Grid<double> > & bg_vs,
                                                 std::vector<NRLib::Grid<double> > & bg_rho,
                                                 const std::vector<NRLib::Well>    & wells,
                                                 MultiIntervalGrid                 * multiple_interval_grid,
                                                 const ModelSettings               * model_settings,
                                                 std::string                       & err_text)
{
  LogKit::LogFormatted(LogKit::Low,"MultiInterval background models:\n");

  //std::vector<int> correlation_structure = model_settings->getCorrelationStructure();
  //std::vector<int> erosion_priority      = model_settings->getErosionPriority();

  int    n_wells    = wells.size(); //model_settings->getNumberOfWells();
  int    n_intervals = multiple_interval_grid->GetNIntervals();
  //float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

  //std::vector<Surface> surface(n_zones+1);
  //for(int i=0; i<n_zones+1; i++)
  //  surface[i] = Surface(surface_files[i]);

  std::vector<StormContGrid> vp_zones(n_intervals); //NRLib::Grid? Get from multiintervalgrid?
  std::vector<StormContGrid> vs_zones(n_intervals);
  std::vector<StormContGrid> rho_zones(n_intervals);

  std::string err_text_tmp = "";
  BuildSeismicPropertyIntervals(vp_zones,
                                vs_zones,
                                rho_zones,
                                multiple_interval_grid,
                                err_text_tmp);

  if(err_text_tmp != "") {
    err_text += err_text_tmp;
  }
  else {

    std::string name_vp  = "Vp";
    std::string name_vs  = "Vs";
    std::string name_rho = "Rho";

    std::vector<float *> trend_vp_zone(n_intervals);
    std::vector<float *> trend_vs_zone(n_intervals);
    std::vector<float *> trend_rho_zone(n_intervals);

    std::vector<const NRLib::Surface<double> *> surfaces; //(n_intervals+1);

    const NRLib::Surface<double> & top_surface = multiple_interval_grid->GetIntervalSimbox(0)->GetTopSurface();
    surfaces.push_back(&top_surface);

    //std::vector<const NRLib::Surface<double>& > surfaces;
    //surfaces.push_back(tmp);

    //surfaces[0] = tmp; //multiple_interval_grid->GetSimbox(0)->GetTopSurface();
    //surfaces.push_back(&multiple_interval_grid->GetSimbox(0)->GetTopSurface());

    for(int i=0; i<n_intervals; i++) {

      const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i);
      surfaces.push_back(&simbox->GetBotSurface());

      const CovGrid2D & cov_grid_2D = makeCovGrid2D(simbox, model_settings->getBackgroundVario(), model_settings->getDebugFlag());
      float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

      LogKit::LogFormatted(LogKit::Low,"\nInterval%2d:",i+1);

      int nz = static_cast<int>(vp_zones[i].GetNK());

      std::vector<float *> well_trend_vp(n_wells);
      std::vector<float *> well_trend_vs(n_wells);
      std::vector<float *> well_trend_rho(n_wells);

      std::vector<float *> high_cut_well_trend_vp(n_wells);
      std::vector<float *> high_cut_well_trend_vs(n_wells);
      std::vector<float *> high_cut_well_trend_rho(n_wells);

      StormContGrid eroded_zone;

      BuildErodedIntervals(eroded_zone,
                           nz,
                           simbox);

      std::vector<bool> hit_zone(n_wells);
      checkWellHitsZone(hit_zone, wells, eroded_zone, n_wells);

      std::vector<BlockedLogsCommon *> blocked_logs(n_wells);

      getWellTrendsZone(model_settings, blocked_logs, well_trend_vp, high_cut_well_trend_vp, wells, eroded_zone, hit_zone, nz, name_vp,  i);
      getWellTrendsZone(model_settings, blocked_logs, well_trend_vs, high_cut_well_trend_vs, wells, eroded_zone, hit_zone, nz, name_vs,  i);
      getWellTrendsZone(model_settings, blocked_logs, well_trend_rho, high_cut_well_trend_rho, wells, eroded_zone, hit_zone, nz, name_rho, i);

      trend_vp_zone[i] = new float[nz];
      trend_vs_zone[i] = new float[nz];
      trend_rho_zone[i] = new float[nz];

      float * avg_dev_vp_zone = new float[n_wells];
      float * avg_dev_vs_zone = new float[n_wells];
      float * avg_dev_rho_zone = new float[n_wells];

      calculateBackgroundTrend(trend_vp_zone[i],
                               avg_dev_vp_zone,
                               nz,
                               dz,
                               model_settings->getAlphaMin(),
                               model_settings->getAlphaMax(),
                               model_settings->getMaxHzBackground(),
                               well_trend_vp,
                               high_cut_well_trend_vp,
                               name_vp);
      calculateBackgroundTrend(trend_vs_zone[i],
                               avg_dev_vs_zone,
                               nz,
                               dz,
                               model_settings->getBetaMin(),
                               model_settings->getBetaMax(),
                               model_settings->getMaxHzBackground(),
                               well_trend_vs,
                               high_cut_well_trend_vs,
                               name_vs);
      calculateBackgroundTrend(trend_rho_zone[i],
                               avg_dev_rho_zone,
                               nz,
                               dz,
                               model_settings->getRhoMin(),
                               model_settings->getRhoMax(),
                               model_settings->getMaxHzBackground(),
                               well_trend_rho,
                               high_cut_well_trend_rho,
                               name_rho);

      writeDeviationsFromVerticalTrend(avg_dev_vp_zone,
                                       avg_dev_vs_zone,
                                       avg_dev_rho_zone,
                                       trend_vp_zone[i],
                                       trend_vs_zone[i],
                                       trend_rho_zone[i],
                                       wells,
                                       n_wells,
                                       nz);

      std::vector<std::vector<double > > bl_vp(n_wells); // bl = blocked logs
      std::vector<std::vector<double > > bl_vs(n_wells);
      std::vector<std::vector<double > > bl_rho(n_wells);
      std::vector<std::vector<double > > vt_vp(n_wells);   // vt = vertical trend
      std::vector<std::vector<double > > vt_vs(n_wells);
      std::vector<std::vector<double > > vt_rho(n_wells);

      std::vector<const std::vector<int> > ipos(n_wells); //std::vector<const int *> ipos(n_wells);
      std::vector<const std::vector<int> > jpos(n_wells);
      std::vector<const std::vector<int> > kpos(n_wells);

      int              tot_blocks;
      std::vector<int> n_blocks(n_wells);

      getKrigingWellTrendsZone(blocked_logs,
                               bl_vp,bl_vs,bl_rho,
                               vt_vp,vt_vs,vt_rho,
                               ipos,jpos,kpos,
                               n_blocks,tot_blocks,nz);

      std::vector<KrigingData2D> kriging_data_vp(nz);
      std::vector<KrigingData2D> kriging_data_vs(nz);
      std::vector<KrigingData2D> kriging_data_rho(nz);

      setupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
                         trend_vp_zone[i],trend_vs_zone[i],trend_rho_zone[i],
                         model_settings->getOutputGridsElastic(),
                         nz,dz,tot_blocks,n_blocks,
                         bl_vp,bl_vs,bl_rho,
                         vt_vp,vt_vs,vt_rho,
                         ipos,jpos,kpos);

      makeKrigedZone(kriging_data_vp, trend_vp_zone[i], vp_zones[i], cov_grid_2D);
      makeKrigedZone(kriging_data_vs, trend_vs_zone[i], vs_zones[i], cov_grid_2D);
      makeKrigedZone(kriging_data_rho, trend_rho_zone[i], rho_zones[i], cov_grid_2D);

      delete [] avg_dev_vp_zone;
      delete [] avg_dev_vs_zone;
      delete [] avg_dev_rho_zone;

      for(int j=0; j<n_wells; j++) {
        delete [] well_trend_vp[j];
        delete [] well_trend_vs[j];
        delete [] well_trend_rho[j];

        delete [] high_cut_well_trend_vp[j];
        delete [] high_cut_well_trend_vs[j];
        delete [] high_cut_well_trend_rho[j];

        delete blocked_logs[j];
      }

    delete &cov_grid_2D;

    } //Intervals

    MakeMultiIntervalBackground(bg_vp, bg_vs, bg_rho,
                                vp_zones, vs_zones, rho_zones,
                                multiple_interval_grid,
                                surfaces,
                                model_settings->getSurfaceUncertainty(),
                                model_settings->getFileGrid(),
                                "multiinterval");


    bool write3D = ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);

    if(write3D) {
      writeMultiIntervalTrendsToFile(trend_vp_zone, trend_vs_zone, trend_rho_zone,
                                     vp_zones, vs_zones, rho_zones,
                                     multiple_interval_grid,
                                     surfaces,
                                     model_settings->getSurfaceUncertainty(),
                                     model_settings->getFileGrid());

    }

    //delete &cov_grid_2D;

    for(int i=0; i<n_intervals; i++) {
      delete [] trend_vp_zone[i];
      delete [] trend_vs_zone[i];
      delete [] trend_rho_zone[i];

      //delete eroded_surfaces[i];
    }
  }

}

void
Background::MakeMultizoneBackground(NRLib::Grid<double>              & bg_vp,
                                    NRLib::Grid<double>              & bg_vs,
                                    NRLib::Grid<double>              & bg_rho,
                                    const std::vector<StormContGrid> & vp_zones,
                                    const std::vector<StormContGrid> & vs_zones,
                                    const std::vector<StormContGrid> & rho_zones,
                                    const Simbox                     * simbox,
                                    const std::vector<int>           & erosion_priority,
                                    const std::vector<Surface>       & surface,
                                    const std::vector<double>        & surface_uncertainty,
                                    const bool                         is_file,
                                    const std::string                & type) const
{

  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  int n_zones = static_cast<int>(vp_zones.size());

  int nx = simbox->getnx();
  int ny = simbox->getny();
  int nz = simbox->getnz();

  const int nxp  = nx;
  const int nyp  = ny;
  const int nzp  = nz;
  const int rnxp = 2*(nxp/2 + 1);

  float monitorSize = std::max(1.0f, static_cast<float>(nzp)*0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  bg_vp.Resize(nx, ny, nz);
  bg_vs.Resize(nx, ny, nz);
  bg_rho.Resize(nx, ny, nz);

  //bg_vp   = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, is_file);
  //bg_vs   = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, is_file);
  //bg_rho  = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, is_file);

  //bg_vp->createRealGrid();
  //bg_vs->createRealGrid();
  //bg_rho->createRealGrid();

  //bg_vp->setType(FFTGrid::PARAMETER);
  //bg_vs->setType(FFTGrid::PARAMETER);
  //bg_rho->setType(FFTGrid::PARAMETER);

  //bg_vp->setAccessMode(FFTGrid::RANDOMACCESS);
  //bg_vs->setAccessMode(FFTGrid::RANDOMACCESS);
  //bg_vs->setAccessMode(FFTGrid::RANDOMACCESS);

  // Beta distributed uncertainty on each surface
  // Note that the upper and lower surfaces not are assigned Beta distributions as these have zero uncertainty
  std::vector<NRLib::Beta> horizon_distributions(n_zones+1);
  for(int zone=1; zone<n_zones; zone++) {
    horizon_distributions[zone] = NRLib::Beta(-surface_uncertainty[zone], surface_uncertainty[zone], 2, 2);
  }

  NRLib::Grid<double> z_surface(rnxp, nyp, n_zones+1);
  for(int i=0; i<rnxp; i++) {
    for(int j=0; j<nyp; j++) {
      double x;
      double y;
      simbox->getXYCoord(i, j, x, y);

      for(int k=0; k<n_zones+1; k++)
        z_surface(i, j, k) = surface[k].GetZ(x,y);
    }
  }

 for(int k=0; k<nzp; k++) {

    for(int j=0; j<nyp; j++) {

      for(int i=0; i<rnxp; i++) {

        if(i<nx) {

          double x;
          double y;

          simbox->getXYCoord(i, j, x, y);

          // Calculate z directly to decrease computation time
          double z = z_surface(i,j,0)+(z_surface(i,j,n_zones)-z_surface(i,j,0))*static_cast<double>(k+0.5)/static_cast<double>(nzp);

          std::vector<double> z_relative(n_zones+1);
          for(int zone=0; zone<n_zones+1; zone++)
            z_relative[zone] = z - z_surface(i,j,zone);

          std::vector<double> zone_probability(n_zones);

          ComputeZoneProbability(z_relative, horizon_distributions, erosion_priority, zone_probability);

          double vp  = 0;
          double vs  = 0;
          double rho = 0;

          for(int zone=0; zone<n_zones; zone++) {

            if(zone_probability[zone] > 0) {
              size_t ind1;
              size_t ind2;
              double t;

              vp_zones[zone].FindZInterpolatedIndex(x, y, z, ind1, ind2, t);

              double vp_zone_new  = vp_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
              double vs_zone_new  = vs_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
              double rho_zone_new = rho_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);

              vp  +=  vp_zone_new  * zone_probability[zone];
              vs  +=  vs_zone_new  * zone_probability[zone];
              rho +=  rho_zone_new * zone_probability[zone];
            }
          }

          bg_vp(i,j,k) = vp;
          bg_vs(i,j,k) = vs;
          bg_rho(i,j,k) = rho;
        }

        else {
          bg_vp(i,j,k) = 0;
          bg_vs(i,j,k) = 0;
          bg_rho(i,j,k) = 0;
        }
      }
    }

    // Log progress
    if (k+1 >= static_cast<int>(nextMonitor)) {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }

  //bg_vp->endAccess();
  //bg_vs->endAccess();
  //bg_rho->endAccess();
}

void
Background::MakeMultiIntervalBackground(std::vector<NRLib::Grid<double> >          & bg_vp,
                                        std::vector<NRLib::Grid<double> >          & bg_vs,
                                        std::vector<NRLib::Grid<double> >          & bg_rho,
                                        const std::vector<StormContGrid>           & vp_zones,
                                        const std::vector<StormContGrid>           & vs_zones,
                                        const std::vector<StormContGrid>           & rho_zones,
                                        MultiIntervalGrid                          * multiple_interval_grid,
                                        std::vector<const NRLib::Surface<double> *>  surfaces,
                                        const std::vector<double>                  & surface_uncertainty,
                                        const bool                                   is_file,
                                        const std::string                          & type) const
{

  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  int n_intervals = multiple_interval_grid->GetNIntervals();
  const std::vector<int> & erosion_priority = multiple_interval_grid->GetErosionPriorities();

  for(int i_interval = 0; i_interval < n_intervals; i_interval++) {

    const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i_interval);

    int nx = simbox->getnx();
    int ny = simbox->getny();
    int nz = simbox->getnz();

    const int nxp  = nx;
    const int nyp  = ny;
    const int nzp  = nz;
    const int rnxp = 2*(nxp/2 + 1);

    float monitorSize = std::max(1.0f, static_cast<float>(nzp)*0.02f);
    float nextMonitor = monitorSize;
    std::cout
      << "\n  0%       20%       40%       60%       80%      100%"
      << "\n  |    |    |    |    |    |    |    |    |    |    |  "
      << "\n  ^";

    bg_vp[i_interval].Resize(nx, ny, nz);
    bg_vs[i_interval].Resize(nx, ny, nz);
    bg_rho[i_interval].Resize(nx, ny, nz);

    // Beta distributed uncertainty on each surface
    // Note that the upper and lower surfaces not are assigned Beta distributions as these have zero uncertainty
    std::vector<NRLib::Beta> horizon_distributions(n_intervals+1);
    for(int zone=1; zone<n_intervals; zone++) {
      horizon_distributions[zone] = NRLib::Beta(-surface_uncertainty[zone], surface_uncertainty[zone], 2, 2);
    }

    NRLib::Grid<double> z_surface(rnxp, nyp, n_intervals+1);
    for(int i=0; i<rnxp; i++) {
      for(int j=0; j<nyp; j++) {
        double x;
        double y;
        simbox->getXYCoord(i, j, x, y);

        for(int k=0; k<n_intervals+1; k++)
          z_surface(i, j, k) = surfaces[k]->GetZ(x,y); //H?
      }
    }

   for(int k=0; k<nzp; k++) {

      for(int j=0; j<nyp; j++) {

        for(int i=0; i<rnxp; i++) {

          if(i<nx) {

            double x;
            double y;

            simbox->getXYCoord(i, j, x, y);

            // Calculate z directly to decrease computation time
            double z = z_surface(i,j,0)+(z_surface(i,j,n_intervals)-z_surface(i,j,0))*static_cast<double>(k+0.5)/static_cast<double>(nzp);

            std::vector<double> z_relative(n_intervals+1);
            for(int zone=0; zone<n_intervals+1; zone++)
              z_relative[zone] = z - z_surface(i,j,zone);

            std::vector<double> zone_probability(n_intervals);

            ComputeZoneProbability(z_relative, horizon_distributions, erosion_priority, zone_probability);

            double vp  = 0;
            double vs  = 0;
            double rho = 0;

            //for(int zone=0; zone<n_intervals; zone++) { //H For Multizone background this is used to create a a weighted vp from all zones (?). In Multiinterval we store vp per interval, remove this weightings?

            //  if(zone_probability[zone] > 0) {
            //    size_t ind1;
            //    size_t ind2;
            //    double t;

            //    vp_zones[zone].FindZInterpolatedIndex(x, y, z, ind1, ind2, t);

            //    double vp_zone_new  = vp_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
            //    double vs_zone_new  = vs_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
            //    double rho_zone_new = rho_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);

            //    vp  +=  vp_zone_new  * zone_probability[zone];
            //    vs  +=  vs_zone_new  * zone_probability[zone];
            //    rho +=  rho_zone_new * zone_probability[zone];
            //  }
            //}

            //H For multinterval. No weightings between intervals, no zone_probability. Correct?
            size_t ind1;
            size_t ind2;
            double t;

            vp_zones[i_interval].FindZInterpolatedIndex(x, y, z, ind1, ind2, t);
            vp = vp_zones[i_interval].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
            vs = vs_zones[i_interval].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
            rho = rho_zones[i_interval].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);

            bg_vp[i_interval](i, j, k) = vp;
            bg_vs[i_interval](i, j, k) = vs;
            bg_rho[i_interval](i, j, k) = rho;
          }
          else {
            bg_vp[i_interval](i, j, k) = 0;
            bg_vs[i_interval](i, j, k) = 0;
            bg_rho[i_interval](i, j, k) = 0;
          }
        }
      }

      // Log progress
      if (k+1 >= static_cast<int>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }
  }
}

//---------------------------------------------------------------------------
void
Background::ComputeZoneProbability(const std::vector<double>      & z,
                                   const std::vector<NRLib::Beta> & horizon_distributions,
                                   const std::vector<int>         & erosion_priority,
                                   std::vector<double>            & zone_probability) const
{

  int nZones = static_cast<int>(zone_probability.size());

  std::vector<double> horizon_cdf(nZones+1, 0);
  horizon_cdf[0] = 1; //The lower surface has cdf 1, whereas the upper surface has cdf 0
  for(int i=1; i<nZones; i++)
    horizon_cdf[i] = horizon_distributions[i].Cdf(z[i]);

  for(int zone=0; zone<nZones; zone++) {
    //Initialize with probability that we are below top surface for zone
    double prob = horizon_cdf[zone];

    //Multiply with probability that we are above base surface for zone
    prob *= (1-horizon_cdf[zone+1]);

    //We may be eroded from above. Must consider the surfaces that
    //1. Are above top in the standard sequence.
    //2. Have lower erosion priority number than the top.
    //3. Have no horizons with lower erosion priority number between it and top.
    int min_erosion = erosion_priority[zone];
    for(int prev_hor = zone-1; prev_hor >=0; prev_hor--) {
      if(erosion_priority[prev_hor] < min_erosion) {
        prob        *= horizon_cdf[prev_hor];
        min_erosion  = erosion_priority[prev_hor]; //Those with higher number stop in this
      }
    }

    //We may be eroded from below. Must consider the surfaces that
    //1. Are below base in the standard sequence.
    //2. Have lower erosion priority number than the base.
    //3. Have no horizons with lower erosion priority number between it and base.
    min_erosion = erosion_priority[zone+1];
    for(int late_hor = zone+2; late_hor < nZones+1; late_hor++) {
      if(erosion_priority[late_hor] < min_erosion) {
        prob        *= (1-horizon_cdf[late_hor]);
        min_erosion  = erosion_priority[late_hor]; //Those with higher number stop in this
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
Background::BuildErodedIntervals(StormContGrid                & eroded_interval,
                                 //const std::vector<Surface>   & eroded_surfaces,
                                 const int                    & nz,
                                 const Simbox                 * simbox) const
{
  int    nx        = simbox->getnx();
  int    ny        = simbox->getny();
  double x_min     = simbox->GetXMin();
  double y_min     = simbox->GetYMin();
  double lx        = simbox->GetLX();
  double ly        = simbox->GetLY();
  double angle     = simbox->getAngle();

  const Surface * top_eroded_surface = simbox->GetTopErodedSurface();
  const Surface * base_eroded_surface = simbox->GetBaseErodedSurface();

  NRLib::Volume volume(x_min, y_min, lx, ly, *top_eroded_surface, *base_eroded_surface, angle);

  eroded_interval = StormContGrid(volume, nx, ny, nz);

}
//---------------------------------------------------------------------------
void
Background::ErodeAllSurfaces(std::vector<Surface *>     & eroded_surfaces,
                             const std::vector<int>     & erosion_priority,
                             const std::vector<Surface> & surface,
                             const Simbox               * simbox) const
{
  int    n_surf     = static_cast<int>(eroded_surfaces.size());

  for(int i=0; i<n_surf; i++) {
    int l=0;
    while(i+1 != erosion_priority[l])
      l++;

    Surface * temp_surface = new Surface(surface[l]);

    //Find closest eroded surface downward
    for(int k=l+1; k<n_surf; k++) {
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
Background::BuildSeismicPropertyZones(std::vector<StormContGrid> & vp_zones,
                                      std::vector<StormContGrid> & vs_zones,
                                      std::vector<StormContGrid> & rho_zones,
                                      const std::vector<Surface> & surface,
                                      const std::vector<int>     & correlation_structure,
                                      const Simbox               * simbox,
                                      const float                & dz,
                                      std::string                & err_text) const
{
  int    n_zones    = static_cast<int>(vp_zones.size());
  int    nx        = simbox->getnx();
  int    ny        = simbox->getny();
  double x_min     = simbox->GetXMin();
  double y_min     = simbox->GetYMin();
  double lx        = simbox->GetLX();
  double ly        = simbox->GetLY();
  double angle     = simbox->getAngle();

  for(int i=1; i<n_zones+1; i++) {
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

        if(z_top == top_missing) {
          const std::string name = top.GetName();
          err_text += "ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n";
          //LogKit::LogFormatted(LogKit::Low,"ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n");
          //exit(1);
        }
        else if(z_base == base_missing) {
          const std::string name = base.GetName();
          err_text += "ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n";
          //LogKit::LogFormatted(LogKit::Low,"ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n");
          //exit(1);
        }

        if(z_base-z_top > max_distance) {
          if(z_top != top_missing && z_base != base_missing)
            max_distance = z_base-z_top;
        }
      }
    }

    if(max_distance == 0) {
      err_text += "ERROR: Zone number "+NRLib::ToString(i)+" has size zero. Check the that surface "+NRLib::ToString(i)+" is above surface "+NRLib::ToString(i+1)+".\n";
      //LogKit::LogFormatted(LogKit::Low,"ERROR: Zone number "+NRLib::ToString(i)+" has size zero. Check the that surface "+NRLib::ToString(i)+" is above surface "+NRLib::ToString(i+1)+".\n");
      //exit(1);
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

    vp_zones[i-1] = StormContGrid(volume, nx, ny, nz_zone);
    vs_zones[i-1] = StormContGrid(volume, nx, ny, nz_zone);
    rho_zones[i-1] = StormContGrid(volume, nx, ny, nz_zone);

  }
}

//---------------------------------------------------------------------------

void
Background::BuildSeismicPropertyIntervals(std::vector<StormContGrid> & vp_zones,
                                          std::vector<StormContGrid> & vs_zones,
                                          std::vector<StormContGrid> & rho_zones,
                                          MultiIntervalGrid          * multiple_interval_grid,
                                          std::string                & err_text) const
{

  int    n_intervals = multiple_interval_grid->GetNIntervals();

  for(int i=0; i<n_intervals; i++) {

    const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i);

    int    nx        = simbox->getnx();
    int    ny        = simbox->getny();
    double x_min     = simbox->GetXMin();
    double y_min     = simbox->GetYMin();
    double lx        = simbox->GetLX();
    double ly        = simbox->GetLY();
    double angle     = simbox->getAngle();

    float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

    //Surface temp_top;
    //Surface temp_base;
    double  x;
    double  y;
    double  z_top;
    double  z_base;

    //Surface top  = surface[i-1];
    //Surface base = surface[i];

    const NRLib::Surface<double> & top = simbox->GetTopSurface();
    const NRLib::Surface<double> & base = simbox->GetBotSurface();

    double top_missing  = RMISSING; //top.IsMissing GetMissingValue();
    double base_missing = RMISSING; //base.GetMissingValue();

    //Find maximum distance between the surfaces
    double max_distance = 0;

    for(int j=0; j<nx; j++) {
      for(int k=0; k<ny; k++) {
        simbox->getXYCoord(j,k,x,y);

        z_top  = top.GetZ(x,y);
        z_base = base.GetZ(x,y);

        if(z_top == top_missing) {
          //const std::string name = top.GetName();

          LogKit::LogFormatted(LogKit::Low,"ERROR: The top surface for interval \'"+multiple_interval_grid->GetIntervalName(i)+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }
        else if(z_base == base_missing) {
          //const std::string name = base.GetName();

          LogKit::LogFormatted(LogKit::Low,"ERROR: The base surface for interval \'"+multiple_interval_grid->GetIntervalName(i)+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }

        if(z_base-z_top > max_distance) {
          if(z_top != top_missing && z_base != base_missing)
            max_distance = z_base-z_top;
        }
      }
    }

    if(max_distance == 0) {
      LogKit::LogFormatted(LogKit::Low,"ERROR: Zone number "+NRLib::ToString(i)+" has size zero. Check the that surface "+NRLib::ToString(i)+" is above surface "+NRLib::ToString(i+1)+".\n");
      exit(1);
    }

    ////Make new top and base surfaces
    //if(correlation_structure[i] == ModelSettings::TOP) {
    //  temp_top  = top;
    //  temp_base = top;
    //  temp_base.Add(max_distance);
    //}
    //else if(correlation_structure[i] == ModelSettings::BASE) {
    //  temp_top  = base;
    //  temp_top.Subtract(max_distance);
    //  temp_base = base;
    //}
    //else {
    //  temp_top  = top;
    //  temp_base = base;
    //}

    NRLib::Volume volume(x_min, y_min, lx, ly, top, base, angle);

    int nz_zone = static_cast<int>(std::ceil(max_distance/dz));

    vp_zones[i] = StormContGrid(volume, nx, ny, nz_zone);
    vs_zones[i] = StormContGrid(volume, nx, ny, nz_zone);
    rho_zones[i] = StormContGrid(volume, nx, ny, nz_zone);

  }
}

//---------------------------------------------------------------------------
//void
//Background::calculateVelocityDeviations(FFTGrid               * velocity,
//                                        const std::vector<WellData *> & wells,
//                                        const Simbox          * simbox,
//                                        float                *& trendVel,
//                                        float                *& avgDevVel,
//                                        float                 * avgDevAlpha,
//                                        int                     outputFlag,
//                                        int                     nWells)
//{
//  if((outputFlag & IO::BACKGROUND_TREND) > 0) {
//    std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + "VpFromFile";
//    velocity->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
//  }
//
//  //
//  // Calculate deviation between well data and trend
//  //
//  int maxBlocks = 0;
//  for (int w = 0 ; w < nWells ; w++) {
//    int nBlocks = wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks();
//    if (nBlocks > maxBlocks)
//      maxBlocks = nBlocks;
//  }
//  float * velocityLog = new float[maxBlocks];
//
//  const int nz = simbox->getnz();
//  float * vtAlpha    = new float[nz];
//  float * vtVelocity = new float[nz];
//
//  for (int k=0 ; k<nz ; k++)
//    trendVel[k]=0.0;
//
//  for (int w = 0 ; w < nWells ; w++) {
//    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
//    const float * alphaLog = bl->getAlphaHighCutBackground();
//    bl->getVerticalTrend(alphaLog, vtAlpha);
//    bl->getBlockedGrid(velocity, velocityLog);
//    bl->getVerticalTrend(velocityLog, vtVelocity);
//    float sumDev = 0.0f;
//    int count = 0;
//    for (int k = 0 ; k < nz ; k++) {
//      if (vtAlpha[k] != RMISSING) {
//        trendVel[k] += vtVelocity[k];
//        float diff = exp(vtAlpha[k]) - vtVelocity[k]; // Velocity trend is in exp-domain
//        sumDev += diff*diff;
//        count++;
//      }
//    }
//    if (count > 0)
//      sumDev /= count;
//    avgDevVel[w] = sqrt(sumDev);
//  }
//  delete [] vtVelocity;
//  delete [] vtAlpha;
//  delete [] velocityLog;
//
//  for (int k=0 ; k<nz ; k++)
//    trendVel[k] /= nWells;
//
//  LogKit::LogFormatted(LogKit::Low,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
//  LogKit::LogFormatted(LogKit::Low,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
//  LogKit::LogFormatted(LogKit::Low,"Well             TrendFromFile  TrendFromData\n");
//  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------\n");
//  for (int i=0 ; i<nWells ; i++)
//    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f          %5.1f\n",
//                         wells[i]->getWellname().c_str(),avgDevVel[i],avgDevAlpha[i]);
//}

void
Background::calculateVelocityDeviations(NRLib::Grid<double>                        & velocity, //FFTGrid                                    * velocity,
                                        const std::vector<NRLib::Well>             & wells,
                                        const Simbox                               * simbox,
                                        std::map<std::string, BlockedLogsCommon *> & bl,
                                        std::map<std::string, BlockedLogsCommon *> & bg_bl,
                                        float                                     *& trend_vel,
                                        float                                     *& avg_dev_vel,
                                        float                                      * avg_dev_vp,
                                        int                                          output_flag,
                                        int                                          n_wells)
{
  //H Writing of NRLib::Grid missing.
  //if((output_flag & IO::BACKGROUND_TREND) > 0) {
  //  std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + "VpFromFile";
  //  velocity->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
  //}

  //
  // Calculate deviation between well data and trend
  //
  int max_blocks = 0;
  for (int w = 0 ; w < n_wells ; w++) {
    int n_blocks = bl.find(wells[w].GetWellName())->second->GetNumberOfBlocks(); //wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (n_blocks > max_blocks)
      max_blocks = n_blocks;
  }
  //float * velocity_log = new float[max_blocks];
  std::vector<double> velocity_log(max_blocks);

  const int nz = simbox->getnz();
  float * vt_vp    = new float[nz];
  float * vt_velocity = new float[nz];

  for (int k=0 ; k<nz ; k++)
    trend_vel[k]=0.0;

  for (int w = 0 ; w < n_wells ; w++) {
    //BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    BlockedLogsCommon * blocked_log = bg_bl.find(wells[w].GetWellName())->second;

    //const float * alphaLog = bl->getAlphaHighCutBackground();
    const std::vector<double> vp_log = blocked_log->GetVpHighCutBackground(); //H make sure HighCutBackground is created.
    blocked_log->GetVerticalTrend(vp_log, vt_vp);
    blocked_log->GetBlockedGrid(velocity, velocity_log);
    blocked_log->GetVerticalTrend(velocity_log, vt_velocity);
    float sumDev = 0.0f;
    int count = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (vt_vp[k] != RMISSING) {
        trend_vel[k] += vt_velocity[k];
        float diff = exp(vt_vp[k]) - vt_velocity[k]; // Velocity trend is in exp-domain
        sumDev += diff*diff;
        count++;
      }
    }
    if (count > 0)
      sumDev /= count;
    avg_dev_vel[w] = sqrt(sumDev);
  }
  delete [] vt_velocity;
  delete [] vt_vp;
  //delete [] velocity_log;

  for (int k=0 ; k<nz ; k++)
    trend_vel[k] /= n_wells;

  LogKit::LogFormatted(LogKit::Low,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
  LogKit::LogFormatted(LogKit::Low,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
  LogKit::LogFormatted(LogKit::Low,"Well             TrendFromFile  TrendFromData\n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------\n");
  for (int i=0 ; i<n_wells ; i++)
    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f          %5.1f\n",
                         wells[i].GetWellName().c_str(),avg_dev_vel[i],avg_dev_vp[i]);
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
//void
//Background::getKrigingWellTrends(std::vector<float *>          & blAlpha,
//                                 std::vector<float *>          & blBeta,
//                                 std::vector<float *>          & blRho,
//                                 std::vector<float *>          & vtAlpha,
//                                 std::vector<float *>          & vtBeta,
//                                 std::vector<float *>          & vtRho,
//                                 std::vector<const int *>      & ipos,
//                                 std::vector<const int *>      & jpos,
//                                 std::vector<const int *>      & kpos,
//                                 std::vector<int>              & nBlocks,
//                                 int                           & totBlocks,
//                                 const std::vector<WellData *> & wells,
//                                 const int                     & nWells) const
//{
//  int maxBlocks = 0;
//  totBlocks     = 0;
//
//  for (int w = 0 ; w < nWells ; w++) {
//    nBlocks[w] = wells[w]->getBlockedLogsExtendedBG()->getNumberOfBlocks();
//    totBlocks += nBlocks[w];
//    if (nBlocks[w] > maxBlocks)
//      maxBlocks = nBlocks[w];
//  }
//
//  for(int i=0; i<nWells; i++) {
//    blAlpha[i] = new float[maxBlocks];
//    blBeta[i]  = new float[maxBlocks];
//    blRho[i]   = new float[maxBlocks];
//  }
//
//  for (int w = 0 ; w < nWells ; w++) {
//    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
//
//    Utils::copyVector(bl->getAlphaHighCutBackground(), blAlpha[w], nBlocks[w]);
//    Utils::copyVector(bl->getBetaHighCutBackground(),  blBeta[w],  nBlocks[w]);
//    Utils::copyVector(bl->getRhoHighCutBackground(),   blRho[w],   nBlocks[w]);
//    //
//    // Extract a one-value-for-each-layer array of blocked logs
//    //
//    bl->getVerticalTrend(blAlpha[w], vtAlpha[w]);
//    bl->getVerticalTrend(blBeta[w],  vtBeta[w]);
//    bl->getVerticalTrend(blRho[w],   vtRho[w]);
//
//    ipos[w] = bl->getIpos();
//    jpos[w] = bl->getJpos();
//    kpos[w] = bl->getKpos();
//
//  }
//}
//---------------------------------------------------------------------------
void
Background::getKrigingWellTrends(std::vector<std::vector<double> >          & bl_vp, // std::vector<float *>           & bl_vp,
                                 std::vector<std::vector<double> >          & bl_vs,
                                 std::vector<std::vector<double> >          & bl_rho,
                                 std::vector<std::vector<double> >          & vt_vp,
                                 std::vector<std::vector<double> >          & vt_vs,
                                 std::vector<std::vector<double> >          & vt_rho,
                                 std::vector<const std::vector<int> >       & ipos,
                                 //std::vector<const int *>       & ipos,
                                 std::vector<const std::vector<int> >       & jpos,
                                 std::vector<const std::vector<int> >       & kpos,
                                 std::vector<int>                           & n_blocks,
                                 int                                        & tot_blocks,
                                 const std::vector<NRLib::Well>             & wells,
                                 std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                                 const int                                  & n_wells) const
{
  int max_blocks = 0;
  tot_blocks     = 0;

  for (int w = 0 ; w < n_wells ; w++) {
    n_blocks[w] = bg_blocked_logs.find(wells[w].GetWellName())->second->GetNumberOfBlocks(); //wells[w]->getBlockedLogsExtendedBG()->getNumberOfBlocks();
    tot_blocks += n_blocks[w];
    if (n_blocks[w] > max_blocks)
      max_blocks = n_blocks[w];
  }

  for(int i=0; i<n_wells; i++) {
    bl_vp[i] = std::vector<double>(max_blocks); //new float[maxBlocks];
    bl_vs[i] = std::vector<double>(max_blocks); //new float[maxBlocks];
    bl_rho[i] = std::vector<double>(max_blocks);  //new float[maxBlocks];
  }

  for (int w = 0 ; w < n_wells ; w++) {
    //BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    BlockedLogsCommon * bg_blocked_log  = bg_blocked_logs.find(wells[w].GetWellName())->second;

    bl_vp[w] = bg_blocked_log->GetVpHighCutBackground();
    bl_vs[w] = bg_blocked_log->GetVsHighCutBackground();
    bl_rho[w] = bg_blocked_log->GetRhoHighCutBackground();

    //Utils::copyVector(bg_blocked_log->GetVpHighCutBackground(), bl_vp[w], n_blocks[w]);
    //Utils::copyVector(bg_blocked_log->GetVsHighCutBackground(), bl_vs[w],  n_blocks[w]);
    //Utils::copyVector(bg_blocked_log->GetRhoHighCutBackground(), bl_rho[w],   n_blocks[w]);
    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    bg_blocked_log->GetVerticalTrend(bl_vp[w], vt_vp[w]);
    bg_blocked_log->GetVerticalTrend(bl_vs[w], vt_vs[w]);
    bg_blocked_log->GetVerticalTrend(bl_rho[w], vt_rho[w]);

    ipos[w] = bg_blocked_log->GetIposVector(); //bl->getIpos();
    jpos[w] = bg_blocked_log->GetJposVector(); //bl->getJpos();
    kpos[w] = bg_blocked_log->GetKposVector(); //bl->getKpos();

  }
}
//---------------------------------------------------------------------------
//void
//Background::getKrigingWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
//                                     std::vector<float *>              & blAlpha,
//                                     std::vector<float *>              & blBeta,
//                                     std::vector<float *>              & blRho,
//                                     std::vector<float *>              & vtAlpha,
//                                     std::vector<float *>              & vtBeta,
//                                     std::vector<float *>              & vtRho,
//                                     std::vector<const int *>          & ipos,
//                                     std::vector<const int *>          & jpos,
//                                     std::vector<const int *>          & kpos,
//                                     std::vector<int>                  & nBlocks,
//                                     int                               & totBlocks,
//                                     const int                         & nz) const
//{
//  int nWells    = static_cast<int>(bl.size());
//  int maxBlocks = 0;
//  totBlocks     = 0;
//
//  for (int w = 0 ; w < nWells ; w++) {
//    if(bl[w] != NULL) {
//      nBlocks[w] = bl[w]->getNumberOfBlocks();
//      totBlocks += nBlocks[w];
//      if (nBlocks[w] > maxBlocks)
//        maxBlocks = nBlocks[w];
//    }
//    else
//      nBlocks[w] = 0;
//  }
//
//  for (int w = 0; w < nWells; w++) {
//    if(bl[w] != NULL) {
//      std::vector<float> blAlphaHighCut = bl[w]->getAlphaHighCutBackground();
//      std::vector<float> blBetaHighCut  = bl[w]->getBetaHighCutBackground();
//      std::vector<float> blRhoHighCut   = bl[w]->getRhoHighCutBackground();
//
//      float * blAlphaCopy = new float[maxBlocks];
//      float * blBetaCopy  = new float[maxBlocks];
//      float * blRhoCopy   = new float[maxBlocks];
//
//      for(int i=0; i<nBlocks[w]; i++) {
//        blAlphaCopy[i] = blAlphaHighCut[i];
//        blBetaCopy[i]  = blBetaHighCut[i];
//        blRhoCopy[i]   = blRhoHighCut[i];
//      }
//      blAlpha[w] = blAlphaCopy;
//      blBeta[w]  = blBetaCopy;
//      blRho[w]   = blRhoCopy;
//      //
//      // Extract a one-value-for-each-layer array of blocked logs
//      //
//      vtAlpha[w] = new float[nz];
//      vtBeta[w]  = new float[nz];
//      vtRho[w]   = new float[nz];
//
//      bl[w]->getVerticalTrend(blAlphaHighCut, vtAlpha[w]);
//      bl[w]->getVerticalTrend(blBetaHighCut,  vtBeta[w]);
//      bl[w]->getVerticalTrend(blRhoHighCut,   vtRho[w]);
//
//      ipos[w] = bl[w]->getIpos();
//      jpos[w] = bl[w]->getJpos();
//      kpos[w] = bl[w]->getKpos();
//    }
//    else {
//      vtAlpha[w] = NULL;
//      vtBeta[w]  = NULL;
//      vtRho[w]   = NULL;
//
//      blAlpha[w] = NULL;
//      blBeta[w]  = NULL;
//      blRho[w]   = NULL;
//    }
//  }
//}

//---------------------------------------------------------------------------
void
Background::getKrigingWellTrendsZone(std::vector<BlockedLogsCommon *>  & blocked_logs,
                                     std::vector<std::vector<double> >           & bl_vp,
                                     std::vector<std::vector<double> >           & bl_vs,
                                     std::vector<std::vector<double> >           & bl_rho,
                                     std::vector<std::vector<double> >           & vt_vp,
                                     std::vector<std::vector<double> >           & vt_vs,
                                     std::vector<std::vector<double> >           & vt_rho,
                                     //std::vector<float *>              & blAlpha,
                                     //std::vector<float *>              & blBeta,
                                     //std::vector<float *>              & blRho,
                                     //std::vector<float *>              & vtAlpha,
                                     //std::vector<float *>              & vtBeta,
                                     //std::vector<float *>              & vtRho,
                                     std::vector<const std::vector<int> >       & ipos,
                                     std::vector<const std::vector<int> >       & jpos,
                                     std::vector<const std::vector<int> >       & kpos,
                                     //std::vector<const int *>          & ipos,
                                     //std::vector<const int *>          & jpos,
                                     //std::vector<const int *>          & kpos,
                                     std::vector<int>                  & n_blocks,
                                     int                               & tot_blocks,
                                     const int                         & nz) const
{
  int n_wells    = static_cast<int>(blocked_logs.size());
  int max_blocks = 0;
  tot_blocks     = 0;

  for (int w = 0 ; w < n_wells ; w++) {
    if(blocked_logs[w] != NULL) {
      n_blocks[w] = blocked_logs[w]->GetNumberOfBlocks();
      tot_blocks += n_blocks[w];
      if (n_blocks[w] > max_blocks)
        max_blocks = n_blocks[w];
    }
    else
      n_blocks[w] = 0;
  }

  for (int w = 0; w < n_wells; w++) {
    if(blocked_logs[w] != NULL) {
      std::vector<double> bl_vp_high_cut = blocked_logs[w]->GetVpHighCutBackground();
      std::vector<double> bl_vs_high_cut = blocked_logs[w]->GetVsHighCutBackground();
      std::vector<double> bl_rho_high_cut = blocked_logs[w]->GetRhoHighCutBackground();

      //float * blAlphaCopy = new float[max_blocks];
      //float * blBetaCopy  = new float[max_blocks];
      //float * blRhoCopy   = new float[max_blocks];

      std::vector<double> bl_vp_copy(max_blocks);
      std::vector<double> bl_vs_copy(max_blocks);
      std::vector<double> bl_rho_copy(max_blocks);

      for(int i=0; i<n_blocks[w]; i++) {
        bl_vp_copy[i] = bl_vp_high_cut[i];
        bl_vs_copy[i] = bl_vs_high_cut[i];
        bl_rho_copy[i] = bl_rho_high_cut[i];
      }
      bl_vp[w] = bl_vp_copy;
      bl_vs[w]  = bl_vs_copy;
      bl_rho[w]   = bl_rho_copy;
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      vt_vp[w] = std::vector<double>(nz); //new float[nz];
      vt_vs[w]  = std::vector<double>(nz); //new float[nz];
      vt_rho[w]   = std::vector<double>(nz); //new float[nz];

      blocked_logs[w]->GetVerticalTrend(bl_vp_high_cut, vt_vp[w]);
      blocked_logs[w]->GetVerticalTrend(bl_vs_high_cut,  vt_vs[w]);
      blocked_logs[w]->GetVerticalTrend(bl_rho_high_cut,   vt_rho[w]);

      const std::vector<int> & ipos_tmp = blocked_logs[w]->GetIposVector();
      const std::vector<int> & jpos_tmp = blocked_logs[w]->GetJposVector();
      const std::vector<int> & kpos_tmp = blocked_logs[w]->GetKposVector();

      //int * i_tmp;
      //int * j_tmp;
      //int * k_tmp;

      std::vector<int> i_tmp;
      std::vector<int> j_tmp;
      std::vector<int> k_tmp;

      for(size_t i = 0; i < ipos_tmp.size(); i++)
        i_tmp[i] = ipos_tmp[i];
      for(size_t j = 0; j < jpos_tmp.size(); j++)
        i_tmp[j] = ipos_tmp[j];
      for(size_t k = 0; k < kpos_tmp.size(); k++)
        i_tmp[k] = ipos_tmp[k];

      ipos[w] = i_tmp;
      jpos[w] = j_tmp;
      kpos[w] = k_tmp;

      //ipos[w] = bl[w]->GetIposVector();
      //jpos[w] = bl[w]->GetJposVector();
      //kpos[w] = bl[w]->GetKposVector();
    }
    else {
      vt_vp[w] = std::vector<double>(0); //NULL;
      vt_vs[w] = std::vector<double>(0); //NULL;
      vt_rho[w] = std::vector<double>(0); //NULL;

      bl_vp[w] = std::vector<double>(0); //NULL;
      bl_vs[w] = std::vector<double>(0); //NULL;
      bl_rho[w] = std::vector<double>(0); //NULL;
    }
  }
}

//---------------------------------------------------------------------------
//void
//Background::getWellTrends(std::vector<float *>          & wellTrend,
//                          std::vector<float *>          & highCutWellTrend,
//                          const std::vector<WellData *> & wells,
//                          const int                     & nz,
//                          const std::string             & name) const
//{
//  int nWells = static_cast<int>(wellTrend.size());
//  int iWells = 0;
//
//  for (int w = 0 ; w < nWells ; w++) {
//    if (wells[w]->getUseForBackgroundTrend()) {
//      BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
//      if(bl != NULL) {
//        wellTrend[w] = new float[nz];
//        if (name == "Vp")
//          bl->getVerticalTrend(bl->getAlpha(), wellTrend[w]);
//        else if (name == "Vs")
//          bl->getVerticalTrend(bl->getBeta(), wellTrend[w]);
//        else if (name == "Rho")
//          bl->getVerticalTrend(bl->getRho(), wellTrend[w]);
//        else {
//          LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrends(): ");
//          LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
//          exit(1);
//        }
//        iWells++;
//      }
//      else
//        wellTrend[w] = NULL;
//    }
//    else
//      wellTrend[w] = NULL;
//  }
//  if(iWells == 0) {
//    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::getWellTrends(): There are no wells\n");
//    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
//    exit(1);
//  }
//
//  for (int w = 0 ; w < nWells ; w++) {
//    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
//    if(bl != NULL) {
//      highCutWellTrend[w] = new float[nz];
//      if (name == "Vp")
//        bl->getVerticalTrend(bl->getAlphaHighCutBackground(), highCutWellTrend[w]);
//      else if (name == "Vs")
//        bl->getVerticalTrend(bl->getBetaHighCutBackground(), highCutWellTrend[w]);
//      else if (name == "Rho")
//        bl->getVerticalTrend(bl->getRhoHighCutBackground(), highCutWellTrend[w]);
//      else {
//        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrends(): ");
//        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
//        exit(1);
//      }
//    }
//    else
//      highCutWellTrend[w] = NULL;
//  }
//}
//---------------------------------------------------------------------------
void
Background::getWellTrends(std::vector<float *>                       & well_trend,
                          std::vector<float *>                       & high_cut_well_trend,
                          const std::vector<NRLib::Well>             & wells,
                          std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                          const int                                  & nz,
                          const std::string                          & name) const
{
  int n_wells = static_cast<int>(well_trend.size());
  int i_wells = 0;

  for (int w = 0 ; w < n_wells ; w++) {
    if (wells[w].getUseForBackgroundTrend()) {
      //BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
      BlockedLogsCommon * blocked_log = bg_blocked_logs.find(wells[w].GetWellName())->second;
      if(blocked_log != NULL) {
        well_trend[w] = new float[nz];
        if (name == "Vp")
          blocked_log->GetVerticalTrend(blocked_log->GetVpBlocked(), well_trend[w]); //bl->getAlpha()
        else if (name == "Vs")
          blocked_log->GetVerticalTrend(blocked_log->GetVsBlocked(), well_trend[w]);
        else if (name == "Rho")
          blocked_log->GetVerticalTrend(blocked_log->GetRhoBlocked(), well_trend[w]);
        else {
          LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrends(): ");
          LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
          exit(1);
        }
        i_wells++;
      }
      else
        well_trend[w] = NULL;
    }
    else
      well_trend[w] = NULL;
  }
  if(i_wells == 0) {
    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::getWellTrends(): There are no wells\n");
    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
    exit(1);
  }

  for (int w = 0 ; w < n_wells ; w++) {
    //BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    BlockedLogsCommon * blocked_log = bg_blocked_logs.find(wells[w].GetWellName())->second;
    if(blocked_log != NULL) {
      high_cut_well_trend[w] = new float[nz];
      if (name == "Vp")
        blocked_log->GetVerticalTrend(blocked_log->GetVpHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Vs")
        blocked_log->GetVerticalTrend(blocked_log->GetVsHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Rho")
        blocked_log->GetVerticalTrend(blocked_log->GetRhoHighCutBackground(), high_cut_well_trend[w]);
      else {
        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrends(): ");
        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
        exit(1);
      }
    }
    else
      high_cut_well_trend[w] = NULL;
  }
}
//---------------------------------------------------------------------------
//void
//Background::getWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
//                              std::vector<float *>              & wellTrend,
//                              std::vector<float *>              & highCutWellTrend,
//                              const std::vector<WellData *>     & wells,
//                              StormContGrid                     & eroded_zone,
//                              const std::vector<bool>           & hitZone,
//                              const int                         & nz,
//                              const std::string                 & name,
//                              const int                         & i) const
//{
//
//  int nValidWellsInZone = 0;
//  int nWells            = static_cast<int>(bl.size());
//
//  std::vector<bool> use_for_background(nWells);
//
//  for(int w=0; w<nWells; w++) {
//    if(hitZone[w] == true) {
//      bl[w] = new BlockedLogsForZone(wells[w], eroded_zone);
//      nValidWellsInZone++;
//    }
//    else
//      bl[w] = NULL;
//
//    if (wells[w]->getUseForBackgroundTrend())
//      use_for_background[w] = true;
//    else
//      use_for_background[w] = false;
//  }
//
//  if(nValidWellsInZone == 0) {
//    LogKit::LogFormatted(LogKit::Low, "Invalid multizone background estimation: No well hits zone number "+NRLib::ToString(i+1)+"\n");
//    exit(1);
//  }
//
//  int iWells = 0;
//
//  for (int w = 0 ;w < nWells; w++) {
//    if (use_for_background[w] == true) {
//      if(bl[w] != NULL) {
//        wellTrend[w] = new float[nz];
//        if (name == "Vp")
//          bl[w]->getVerticalTrend(bl[w]->getAlpha(), wellTrend[w]);
//        else if (name == "Vs")
//          bl[w]->getVerticalTrend(bl[w]->getBeta(), wellTrend[w]);
//        else if (name == "Rho")
//          bl[w]->getVerticalTrend(bl[w]->getRho(), wellTrend[w]);
//        else {
//          LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrendsZone(): ");
//          LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
//          exit(1);
//        }
//        iWells++;
//      }
//      else wellTrend[w] = NULL;
//    }
//    else
//      wellTrend[w] = NULL;
//  }
//  if(iWells == 0) {
//    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::getWellTrendsZone(): There are no wells\n");
//    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
//    exit(1);
//  }
//
//
//  for (int w = 0; w < nWells; w++) {
//    if(bl[w] != NULL) {
//      highCutWellTrend[w] = new float[nz];
//      if (name == "Vp")
//        bl[w]->getVerticalTrend(bl[w]->getAlphaHighCutBackground(), highCutWellTrend[w]);
//      else if (name == "Vs")
//        bl[w]->getVerticalTrend(bl[w]->getBetaHighCutBackground(), highCutWellTrend[w]);
//      else if (name == "Rho")
//        bl[w]->getVerticalTrend(bl[w]->getRhoHighCutBackground(), highCutWellTrend[w]);
//      else {
//        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrendsZone(): ");
//        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
//        exit(1);
//      }
//    }
//    else
//      highCutWellTrend[w] = NULL;
//  }
//}
//---------------------------------------------------------------------------
void
Background::getWellTrendsZone(const ModelSettings                              * model_settings,
                              std::vector<BlockedLogsCommon *>                 & bl,
                              std::vector<float *>                             & well_trend,
                              std::vector<float *>                             & high_cut_well_trend,
                              const std::vector<NRLib::Well>                   & wells,
                              StormContGrid                                    & eroded_zone,
                              //const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                              const std::vector<bool>                          & hit_zone,
                              const int                                        & nz,
                              const std::string                                & name,
                              const int                                        & i) const
{

  int n_valid_wells_in_zone = 0;
  int n_wells            = static_cast<int>(bl.size());

  std::vector<bool> use_for_background(n_wells);

  //NRLib::Well tmp_well = wells[0];
  //BlockedLogsCommon * tmp = mapped_blocked_logs.find(wells[0].GetWellName())->second;

  for(int w=0; w<n_wells; w++) {
    if(hit_zone[w] == true) {
      //bl[w] = new BlockedLogsForZone(wells[w], eroded_zone, tmp);
      //bl[w] = new BlockedLogsForZone(wells[w], eroded_zone, mapped_blocked_logs.find(wells[w].GetWellName())->second); ///H Use BlockedLogsCommon instead of this modified BlockedLogsForZone?

      bl[w] = new BlockedLogsCommon(&wells[w], eroded_zone, model_settings->getMaxHzBackground(), model_settings->getMaxHzSeismic());

      n_valid_wells_in_zone++;
    }
    else
      bl[w] = NULL;

    if (wells[w].getUseForBackgroundTrend())
      use_for_background[w] = true;
    else
      use_for_background[w] = false;
  }

  if(n_valid_wells_in_zone == 0) {
    LogKit::LogFormatted(LogKit::Low, "Invalid multizone background estimation: No well hits zone number "+NRLib::ToString(i+1)+"\n");
    exit(1);
  }

  int iWells = 0;

  for (int w = 0 ;w < n_wells; w++) {
    if (use_for_background[w] == true) {
      if(bl[w] != NULL) {
        well_trend[w] = new float[nz];
        //wellTrend[w].resize(nz);

        if (name == "Vp")
          bl[w]->GetVerticalTrend(bl[w]->GetVpBlocked(), well_trend[w]);
        else if (name == "Vs")
          bl[w]->GetVerticalTrend(bl[w]->GetVsBlocked(), well_trend[w]);
        else if (name == "Rho")
          bl[w]->GetVerticalTrend(bl[w]->GetRhoBlocked(), well_trend[w]);
        else {
          LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrendsZone(): ");
          LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
          exit(1);
        }
        iWells++;
      }
      else well_trend[w] = NULL;
    }
    else
      well_trend[w] = NULL;
  }
  if(iWells == 0) {
    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::getWellTrendsZone(): There are no wells\n");
    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
    exit(1);
  }

  for (int w = 0; w < n_wells; w++) {
    if(bl[w] != NULL) {
      high_cut_well_trend[w] = new float[nz];
      //highCutWellTrend[w].resize(nz);
      if (name == "Vp")
        bl[w]->GetVerticalTrend(bl[w]->GetVpHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Vs")
        bl[w]->GetVerticalTrend(bl[w]->GetVsHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Rho")
        bl[w]->GetVerticalTrend(bl[w]->GetRhoHighCutBackground(), high_cut_well_trend[w]);
      else {
        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::getWellTrendsZone(): ");
        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
        exit(1);
      }
    }
    else
      high_cut_well_trend[w] = NULL;
  }
}
//---------------------------------------------------------------------------
//void
//Background::checkWellHitsZone(std::vector<bool>             & hitZone,
//                              const std::vector<WellData *> & wells,
//                              StormContGrid                 & eroded_zone,
//                              const int                     & nWells) const
//{
//  for(int w=0; w<nWells; w++) {
//    if(wells[w]->checkStormgrid(eroded_zone) == 0) {
//      hitZone[w] = true;
//    }
//    else
//      hitZone[w] = false;
//  }
//}
//---------------------------------------------------------------------------
void
Background::checkWellHitsZone(std::vector<bool>              & hitZone,
                              const std::vector<NRLib::Well> & wells,
                              StormContGrid                  & eroded_zone,
                              const int                      & nWells) const
{
  for(int w=0; w<nWells; w++) {
    if(wells[w].CheckStormgrid(eroded_zone) == 0) {
      hitZone[w] = true;
    }
    else
      hitZone[w] = false;
  }
}
//---------------------------------------------------------------------------
void
Background::writeTrendsToFile(float             * trend,
                              const Simbox      * simbox,
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
Background::writeMultizoneTrendsToFile(const std::vector<float *>   vp_zones,
                                       const std::vector<float *>   vs_zones,
                                       const std::vector<float *>   rho_zones,
                                       std::vector<StormContGrid> & vp_trend_zone,
                                       std::vector<StormContGrid> & vs_trend_zone,
                                       std::vector<StormContGrid> & rho_trend_zone,
                                       const Simbox               * simbox,
                                       const std::vector<int>     & erosion_priority,
                                       const std::vector<Surface> & surface,
                                       const std::vector<double>  & surface_uncertainty,
                                       const bool                   isFile) const
{
  int n_zones = static_cast<int>(vp_zones.size());

  for(int i=0; i<n_zones; i++) {
    makeTrendZone(vp_zones[i], vp_trend_zone[i]);
    makeTrendZone(vs_zones[i], vs_trend_zone[i]);
    makeTrendZone(rho_zones[i], rho_trend_zone[i]);
  }

  //FFTGrid * trend_vp;
  //FFTGrid * trend_vs;
  //FFTGrid * trend_rho;
  NRLib::Grid<double> trend_vp;
  NRLib::Grid<double> trend_vs;
  NRLib::Grid<double> trend_rho;

  MakeMultizoneBackground(trend_vp,
                          trend_vs,
                          trend_rho,
                          vp_trend_zone,
                          vs_trend_zone,
                          rho_trend_zone,
                          simbox,
                          erosion_priority,
                          surface,
                          surface_uncertainty,
                          isFile,
                          "trend in multizone");

  //FFTGrid * exp_trend_vp = copyFFTGrid(trend_vp, true, isFile);
  //FFTGrid * exp_trend_vs = copyFFTGrid(trend_vs, true, isFile);
  //FFTGrid * exp_trend_rho = copyFFTGrid(trend_rho, true, isFile);

  //std::string file_name_vp = IO::PrefixBackground() + IO::PrefixTrend() + "Vp";
  //std::string file_name_vs = IO::PrefixBackground() + IO::PrefixTrend() + "Vs";
  //std::string file_name_rho = IO::PrefixBackground() + IO::PrefixTrend() + "Rho";

  //exp_trend_vp->writeFile(file_name_vp, IO::PathToBackground(), simbox);
  //exp_trend_vs->writeFile(file_name_vs, IO::PathToBackground(), simbox);
  //exp_trend_rho->writeFile(file_name_rho, IO::PathToBackground(), simbox);

  //delete exp_trend_vp;
  //delete exp_trend_vs;
  //delete exp_trend_rho;

  //delete trend_vp;
  //delete trend_vs;
  //delete trend_rho;
}

//---------------------------------------------------------------------------
void
Background::writeMultiIntervalTrendsToFile(const std::vector<float *>   vp_zones,
                                           const std::vector<float *>   vs_zones,
                                           const std::vector<float *>   rho_zones,
                                           std::vector<StormContGrid> & vp_trend_zone,
                                           std::vector<StormContGrid> & vs_trend_zone,
                                           std::vector<StormContGrid> & rho_trend_zone,
                                           //const Simbox                     * simbox,
                                           //const std::vector<int>           & erosion_priority,
                                           MultiIntervalGrid          * multiple_interval_grid,
                                           //const std::vector<int>     & erosion_priority,
                                           //const std::vector<NRLib::Surface<double> > & surfaces,
                                           //std::vector<const NRLib::Surface<double>& > surfaces,
                                           std::vector<const NRLib::Surface<double> *> surfaces,
                                           const std::vector<double>  & surface_uncertainty,
                                           const bool                   is_file) const
{
  int n_intervals = multiple_interval_grid->GetNIntervals();

  for(int i=0; i<n_intervals; i++) {
    makeTrendZone(vp_zones[i], vp_trend_zone[i]);
    makeTrendZone(vs_zones[i], vs_trend_zone[i]);
    makeTrendZone(rho_zones[i], rho_trend_zone[i]);
  }

  //std::vector<FFTGrid *> trend_vp(n_intervals);
  //std::vector<FFTGrid *> trend_vs(n_intervals);
  //std::vector<FFTGrid *> trend_rho(n_intervals);

  std::vector<NRLib::Grid<double> > trend_vp(n_intervals);
  std::vector<NRLib::Grid<double> > trend_vs(n_intervals);
  std::vector<NRLib::Grid<double> > trend_rho(n_intervals);

  MakeMultiIntervalBackground(trend_vp,
                              trend_vs,
                              trend_rho,
                              vp_trend_zone,
                              vs_trend_zone,
                              rho_trend_zone,
                              //simbox,
                              multiple_interval_grid,
                              //erosion_priority,
                              surfaces,
                              surface_uncertainty,
                              is_file,
                              "trend in multiinterval");

  for(int i=0; i < n_intervals; i++) {

    std::string interval_name = multiple_interval_grid->GetIntervalName(i);
    const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i);

    //FFTGrid * exp_trend_vp = copyFFTGrid(trend_vp[i], true, is_file);
    //FFTGrid * exp_trend_vs = copyFFTGrid(trend_vs[i], true, is_file);
    //FFTGrid * exp_trend_rho = copyFFTGrid(trend_rho[i], true, is_file);

    //std::string file_name_vp = "Interval_" + interval_name + "_" + IO::PrefixBackground() + IO::PrefixTrend() + "Vp";
    //std::string file_name_vs = "Interval_" + interval_name + "_" + IO::PrefixBackground() + IO::PrefixTrend() + "Vs";
    //std::string file_name_rho = "Interval_" + interval_name + "_" + IO::PrefixBackground() + IO::PrefixTrend() + "Rho";

    //exp_trend_vp->writeFile(file_name_vp, IO::PathToBackground(), simbox);
    //exp_trend_vs->writeFile(file_name_vs, IO::PathToBackground(), simbox);
    //exp_trend_rho->writeFile(file_name_rho, IO::PathToBackground(), simbox);

    //delete exp_trend_vp;
    //delete exp_trend_vs;
    //delete exp_trend_rho;

    //delete trend_vp[i];
    //delete trend_vs[i];
    //delete trend_rho[i];

  }
}

//-------------------------------------------------------------------------------
void
Background::setupKrigingData2D(std::vector<KrigingData2D>     & kriging_data_vp,
                               std::vector<KrigingData2D>     & kriging_data_vs,
                               std::vector<KrigingData2D>     & kriging_data_rho,
                               float                          * trend_vp,
                               float                          * trend_vs,
                               float                          * trend_rho,
                               const int                        output_flag,
                               const int                      & nz,
                               const float                    & dz,
                               const int                      & tot_blocks,
                               const std::vector<int>         & n_blocks,
                               const std::vector<std::vector<double> > & bl_vp,
                               //const std::vector<float *>     & blAlpha,
                               const std::vector<std::vector<double> > & bl_vs,
                               const std::vector<std::vector<double> > & bl_rho,
                               const std::vector<std::vector<double> > & vt_vp,
                               const std::vector<std::vector<double> > & vt_vs,
                               const std::vector<std::vector<double> > & vt_rho,
                               const std::vector<const std::vector<int> >   ipos,
                               const std::vector<const std::vector<int> >   jpos,
                               const std::vector<const std::vector<int> >   kpos) const
{
  //
  // Although unnecessary, we have chosen to set up kriging data from
  // Vp, Vs and Rho simultaneously. This gives code easier to read.
  //
  const int n_wells = static_cast<int>(bl_vp.size());

  KrigingData3D forLogging(tot_blocks);

  for (int w = 0 ; w < n_wells ; w++) {

    //if(vt_vp[w] != NULL) {
    if(vt_vp[w].size() > 0) {
      std::vector<double> vt_vp_well = vt_vp[w];
      std::vector<double> vt_vs_well = vt_vs[w];
      std::vector<double> vt_rho_well = vt_rho[w];

      //float * vtAlphaWell = vtAlpha[w];
      //float * vtBetaWell  = vtBeta[w];
      //float * vtRhoWell   = vtRho[w];

      std::vector<double> bl_vp_well = bl_vp[w];
      std::vector<double> bl_vs_well = bl_vs[w];
      std::vector<double> bl_rho_well = bl_rho[w];

      //float * blAlphaWell = blAlpha[w];
      //float * blBetaWell  = blBeta[w];
      //float * blRhoWell   = blRho[w];
      //
      // Kriging vertical trend (vt....) against global vertical trend (trend...)
      //
      Kriging1D::krigVector(&vt_vp_well[0], &trend_vp[0], nz, dz);
      Kriging1D::krigVector(&vt_vs_well[0], &trend_vs[0], nz, dz);
      Kriging1D::krigVector(&vt_rho_well[0], &trend_rho[0], nz, dz);
      //
      // Use kriged vertical trend where original log is not defined.
      //
      const std::vector<int> ipos_well = ipos[w];
      const std::vector<int> jpos_well = jpos[w];
      const std::vector<int> kpos_well = kpos[w];

      //const int * iposWell = ipos[w];
      //const int * jposWell = jpos[w];
      //const int * kposWell = kpos[w];

      for (int m = 0 ; m < n_blocks[w] ; m++) {
        int i = ipos_well[m];
        int j = jpos_well[m];
        int k = kpos_well[m];

        if (bl_vp_well[m] == RMISSING)
          bl_vp_well[m] = vt_vp_well[k];
        if (bl_vs_well[m] == RMISSING)
          bl_vs_well[m] = vt_vs_well[k];
        if (bl_rho_well[m] == RMISSING)
          bl_rho_well[m] = vt_rho_well[k];

        kriging_data_vp[k].addData(i, j, static_cast<float>(bl_vp_well[m]));
        kriging_data_vs[k].addData(i, j, static_cast<float>(bl_vs_well[m]));
        kriging_data_rho[k].addData(i, j, static_cast<float>(bl_rho_well[m]));
      }

      forLogging.addData(&bl_vp_well[0], &bl_vs_well[0], &bl_rho_well[0],
        &ipos_well[0],&jpos_well[0],&kpos_well[0],
        n_blocks[w]);
    }

    for (int k=0 ; k<nz ; k++) {
      kriging_data_vp[k].findMeanValues();
      kriging_data_vs[k].findMeanValues();
      kriging_data_rho[k].findMeanValues();
    }
  }

  if((output_flag & IO::BACKGROUND) > 0) {
    forLogging.divide();
    std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + IO::SuffixGeneralData();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    forLogging.writeToFile(fileName);
  }
}

//---------------------------------------------------------------------------
const CovGrid2D &
Background::makeCovGrid2D(const Simbox * simbox,
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
//void
//Background::makeKrigedBackground(const std::vector<KrigingData2D> & krigingData,
//                                 FFTGrid                         *& bgGrid,
//                                 const float                      * trend,
//                                 const Simbox                     * simbox,
//                                 const CovGrid2D                  & covGrid2D,
//                                 const std::string                & type,
//                                 bool                               isFile) const
//{
//  std::string text = "\nBuilding "+type+" background:";
//  LogKit::LogFormatted(LogKit::Low,text);
//
//  const int    nx   = simbox->getnx();
//  const int    ny   = simbox->getny();
//  const int    nz   = simbox->getnz();
//
//  const int    nxp  = nx;
//  const int    nyp  = ny;
//  const int    nzp  = nz;
//  const int    rnxp = 2*(nxp/2 + 1);
//
//  const double x0   = simbox->getx0();
//  const double y0   = simbox->gety0();
//  const double lx   = simbox->getlx();
//  const double ly   = simbox->getly();
//
//  //
//  // Template surface to be kriged
//  //
//  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);
//
//  float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
//  float nextMonitor = monitorSize;
//  std::cout
//    << "\n  0%       20%       40%       60%       80%      100%"
//    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
//    << "\n  ^";
//
//  bgGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
//  bgGrid->createRealGrid();
//  bgGrid->setType(FFTGrid::PARAMETER);
//  bgGrid->setAccessMode(FFTGrid::WRITE);
//
//  for (int k=0 ; k<nzp ; k++)
//  {
//    // Set trend for layer
//    surface.Assign(trend[k]);
//
//    // Kriging of layer
//    Kriging2D::krigSurface(surface, krigingData[k], covGrid2D);
//
//    // Set layer in background model from surface
//    for(int j=0 ; j<nyp ; j++) {
//      for(int i=0 ; i<rnxp ; i++) {
//        if(i<nxp)
//          bgGrid->setNextReal(float(surface(i,j)));
//        else
//          bgGrid->setNextReal(0);  //dummy in padding (but there is no padding)
//      }
//    }
//
//    // Log progress
//    if (k+1 >= static_cast<int>(nextMonitor))
//    {
//      nextMonitor += monitorSize;
//      std::cout << "^";
//      fflush(stdout);
//    }
//  }
//  bgGrid->endAccess();
//}

void
Background::makeKrigedBackground(const std::vector<KrigingData2D> & kriging_data,
                                 NRLib::Grid<double>              & bg_grid, //FFTGrid                         *& bgGrid,
                                 const float                      * trend,
                                 const Simbox                     * simbox,
                                 const CovGrid2D                  & cov_grid_2D,
                                 const std::string                & type) const
                                 //bool                               is_file) const
{
  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  const int    nx   = simbox->getnx();
  const int    ny   = simbox->getny();
  const int    nz   = simbox->getnz();

  //const int    nxp  = nx;
  //const int    nyp  = ny;
  //const int    nzp  = nz;
  //const int    rnxp = 2*(nxp/2 + 1);

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

  //bgGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  //bgGrid->createRealGrid();
  //bgGrid->setType(FFTGrid::PARAMETER);
  //bgGrid->setAccessMode(FFTGrid::WRITE);

  for (int k=0; k < nz; k++)
  {
    // Set trend for layer
    surface.Assign(trend[k]);

    // Kriging of layer
    Kriging2D::krigSurface(surface, kriging_data[k], cov_grid_2D);

    // Set layer in background model from surface
    for(int j=0; j < ny; j++) {
      for(int i=0 ; i < nx ; i++) {
        bg_grid(i, j, k) = surface(i, j);

        //if(i<nxp)
        //  bgGrid->setNextReal(float(surface(i,j)));
        //else
        //  bgGrid->setNextReal(0);  //dummy in padding (but there is no padding)
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
  //bgGrid->endAccess();
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
Background::calculateVerticalTrend(std::vector<float *>   wellTrend,
                                   float                * trend,
                                   float                  logMin,
                                   float                  logMax,
                                   float                  maxHz,
                                   int                    nz,
                                   float                  dz,
                                   const std::string    & name)
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
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(firstValue),firstNonmissing);
      for (int k = 0 ; k < firstNonmissing ; k++) {
        trend[k] = firstValue;
      }
    }
    if (errorTrail) {
      // Fix last part of trend containing missing-values.
      float lastValue = trend[lastNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+parName+" using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [%d,%d] where the log is undefined. The last\n",lastNonmissing+1,nz-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(lastValue),lastNonmissing);
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
//void
//Background::writeDeviationsFromVerticalTrend(const float                   * avg_dev_alpha,
//                                             const float                   * avg_dev_beta,
//                                             const float                   * avg_dev_rho,
//                                             const float                   * trend_alpha,
//                                             const float                   * trend_beta,
//                                             const float                   * trend_rho,
//                                             const std::vector<WellData *> & wells,
//                                             const int                       nWells,
//                                             const int                       nz)
//{
//  float global_mean_alpha = 0.0f;
//  float global_mean_beta  = 0.0f;
//  float global_mean_rho   = 0.0f;
//
//  for (int k=0 ; k<nz ; k++)
//  {
//    global_mean_alpha += exp(trend_alpha[k]);
//    global_mean_beta  += exp(trend_beta[k]);
//    global_mean_rho   += exp(trend_rho[k]);
//  }
//  global_mean_alpha /= nz;
//  global_mean_beta  /= nz;
//  global_mean_rho   /= nz;
//
//  //
//  // Find the relative average deviations (mean of Vp,Vs and Rho deviations).
//  //
//  float * rel_avg_dev = new float[nWells];
//  for (int i=0 ; i<nWells ; i++)
//  {
//    float rel_dev_alpha = avg_dev_alpha[i]/global_mean_alpha;
//    float rel_dev_beta  = avg_dev_beta[i]/global_mean_beta;
//    float rel_dev_rho   = avg_dev_rho[i]/global_mean_rho;
//    rel_avg_dev[i] = (rel_dev_alpha + rel_dev_beta + rel_dev_rho)/3;
//  }
//  //
//  // Sort deviations to find worst well.
//  //
//
//  int * index = new int[nWells];
//  for(int i=0;i<nWells;i++)
//    index[i] = i;
//
//  for (int i=0 ; i<nWells ; i++)
//  {
//    for (int j=i; j<nWells ; j++)
//    {
//      if (rel_avg_dev[index[j]] > rel_avg_dev[index[i]])
//      {
//        int tmp = index[i];
//        index[i] = index[j];
//        index[j] = tmp;
//      }
//    }
//  }
//  //
//  // Print results
//  //
//  if (nWells > 0)
//  {
//    LogKit::LogFormatted(LogKit::Low,"\nSummary of average deviation from vertical trend (well with largest misfit listed first):\n\n");
//    LogKit::LogFormatted(LogKit::Low,"Well                        Vp       Vs      Rho\n");
//    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------\n");
//  }
//  for (int i=0 ; i<nWells ; i++)
//  {
//    int ii = index[i];
//    if(avg_dev_alpha[ii] != RMISSING) {
//      LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f    %5.1f    %5.3f\n", wells[ii]->getWellname().c_str(),
//                           avg_dev_alpha[ii], avg_dev_beta[ii], avg_dev_rho[ii]);
//    }
//  }
//
//  if (nWells == 1)
//  {
//    LogKit::LogFormatted(LogKit::High,"\nNOTE: A deviation may be observed even with one well since the global trend is");
//    LogKit::LogFormatted(LogKit::High,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
//  }
//  delete [] rel_avg_dev;
//  delete [] index;
//}

//-------------------------------------------------------------------------------
void
Background::writeDeviationsFromVerticalTrend(const float                    * avg_dev_alpha,
                                             const float                    * avg_dev_beta,
                                             const float                    * avg_dev_rho,
                                             const float                    * trend_alpha,
                                             const float                    * trend_beta,
                                             const float                    * trend_rho,
                                             const std::vector<NRLib::Well> & wells,
                                             const int                        nWells,
                                             const int                        nz)
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
      LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f    %5.1f    %5.3f\n", wells[ii].GetWellName().c_str(),
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
  int nxp = back_model_[0]->getNxp();
  int nyp = back_model_[0]->getNyp();
  int nzp = back_model_[0]->getNzp();
  float vp = back_model_[0]->getFirstComplexValue().re;
  vp = float(exp(vp/sqrt(float(nxp*nyp*nzp))));
  float vs = back_model_[1]->getFirstComplexValue().re;
  vs = float(exp(vs/sqrt(float(nxp*nyp*nzp))));
  vsvp_ = vs/vp;
}

//-------------------------------------------------------------------------------
void
Background::resampleBackgroundModel(NRLib::Grid<double> & bg_vp, //FFTGrid            *& bg_vp,
                                    NRLib::Grid<double> & bg_vs,
                                    NRLib::Grid<double> & bg_rho,
                                    const Simbox        * time_bg_simbox,
                                    const Simbox        * time_simbox,
                                    const ModelSettings * model_settings)
{
  bool is_file = model_settings->getFileGrid();

  //H Writing of grids missing
  //if((model_settings->getOutputGridsOther() & IO::EXTRA_GRIDS) > 0) {
  //  std::string file_name_1 = IO::PrefixBackground() + "Vp_BackgroundGrid";
  //  std::string file_name_2 = IO::PrefixBackground() + "Vs_BackgroundGrid";
  //  std::string file_name_3 = IO::PrefixBackground() + "Rho_BackgroundGrid";

  //  FFTGrid * exp_vp = copyFFTGrid(bg_vp, true, is_file);
  //  exp_vp->writeFile(file_name_1, IO::PathToBackground(), time_bg_simbox);
  //  delete exp_vp;

  //  FFTGrid * exp_vs = copyFFTGrid(bg_vs, true, is_file);
  //  exp_vs->writeFile(file_name_2, IO::PathToBackground(), time_bg_simbox);
  //  delete exp_vs;

  //  FFTGrid * exp_rho = copyFFTGrid(bg_rho, true, is_file);
  //  exp_rho->writeFile(file_name_3, IO::PathToBackground(), time_bg_simbox);
  //  delete exp_rho;
  //}

  //FFTGrid * res_bg_vp = NULL;
  //FFTGrid * res_bg_vs = NULL;
  //FFTGrid * res_bg_rho = NULL;
  NRLib::Grid<double> res_bg_vp;
  NRLib::Grid<double> res_bg_vs;
  NRLib::Grid<double> res_bg_rho;

  LogKit::LogFormatted(LogKit::Low,"\nResampling background model...\n");
  resampleParameter(res_bg_vp, bg_vp, time_simbox, time_bg_simbox);
  resampleParameter(res_bg_vs, bg_vs, time_simbox, time_bg_simbox);
  resampleParameter(res_bg_rho, bg_rho, time_simbox, time_bg_simbox);

  //if((model_settings->getOutputGridsOther() & IO::EXTRA_GRIDS) > 0) {
  //  std::string file_name_1 = IO::PrefixBackground() + "Vp_InversionGrid";
  //  std::string file_name_2 = IO::PrefixBackground() + "Vs_InversionGrid";
  //  std::string file_name_3 = IO::PrefixBackground() + "Rho_InversionGrid";

  //  FFTGrid * exp_res_vp = copyFFTGrid(res_bg_vp, true, is_file);
  //  exp_res_vp->writeFile(file_name_1, IO::PathToBackground(), time_simbox);
  //  delete exp_res_vp;

  //  FFTGrid * exp_res_vs = copyFFTGrid(res_bg_vs, true, is_file);
  //  exp_res_vs->writeFile(file_name_2, IO::PathToBackground(), time_simbox);
  //  delete exp_res_vs;

  //  FFTGrid * exp_res_rho = copyFFTGrid(res_bg_rho, true, is_file);
  //  exp_res_rho->writeFile(file_name_3, IO::PathToBackground(), time_simbox);
  //  delete exp_res_rho;
  //}

  //delete bg_vp;
  //delete bg_vs;
  //delete bg_rho;

  bg_vp = res_bg_vp;
  bg_vs = res_bg_vs;
  bg_rho = res_bg_rho;
 }

//-------------------------------------------------------------------------------
//void
//Background::resampleParameter(FFTGrid *& pNew,        // Resample to
//                              FFTGrid  * pOld,        // Resample from
//                              const Simbox   * simboxNew,
//                              const Simbox   * simboxOld,
//                              bool       isFile)
//{
//  int nx  = simboxNew->getnx();
//  int ny  = simboxNew->getny();
//  int nz  = simboxNew->getnz();
//  //
//  // Use same padding as for nonresampled cubes
//  //
//  // NBNB-PAL: These grids are unpadded, so all nxp, nyp, ... would probably
//  //           better be replaced by nx, ny, ... to avoid confusion...
//  //
//  int nxp = nx + (pOld->getNxp() - pOld->getNxp());
//  int nyp = ny + (pOld->getNyp() - pOld->getNyp());
//  int nzp = nz + (pOld->getNzp() - pOld->getNzp());
//
//  //
//  // Set up relation between old layer index and new layer index using
//  //
//  // k2 = dz1/dz2 * k1 + (z02 - z01)/dz2    (from dz2*k2 + z02 = dz1*k1 + z01)
//  //
//  double * a = new double[nx*ny];
//  double * b = new double[nx*ny];
//
//  int ij = 0;
//  for(int j=0;j<ny;j++) {
//    for(int i=0;i<nx;i++) {
//      double dzNew = simboxNew->getdz(i,j);
//      double dzOld = simboxOld->getdz(i,j);
//      double z0New = simboxNew->getTop(i,j);
//      double z0Old = simboxOld->getTop(i,j);
//        a[ij] = dzNew/dzOld;
//      b[ij] = (z0New - z0Old)/dzOld;
//      ij++;
//    }
//  }
//
//  //
//  // Resample parameter
//  //
//  pNew = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
//  pNew->createRealGrid();
//  pNew->setType(FFTGrid::PARAMETER);
//  pNew->setAccessMode(FFTGrid::WRITE);
//
//  pOld->setAccessMode(FFTGrid::RANDOMACCESS);
//
//  int rnxp = 2*(nxp/2 + 1);
//
//  double * layer = new double[nx*ny];
//
//  for(int k=0 ; k<nzp ; k++) {
//    //
//    // Map a layer
//    //
//    int ij=0;
//    for (int j=0 ; j<nyp ; j++) {
//      for (int i=0 ; i<rnxp ; i++) {
//        if (i < nx && j < ny && k < nz) {
//          int kOld = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]);
//          layer[ij] = pOld->getRealValue(i, j, kOld);
//          ij++;
//        }
//      }
//    }
//    //
//    // Smooth the layer (equal weighting of all neighbouring cells)
//    //
//    float value;
//    for (int j=0 ; j<nyp ; j++) {
//      for (int i=0 ; i<rnxp ; i++) {
//        if (i < nx && j < ny && k < nz) {
//          int n = 1;
//          double sum = layer[j*nx + i];
//          if (i>1) {
//            sum += layer[j*nx + i - 1];
//            n++;
//          }
//          if (j>1) {
//            sum += layer[(j - 1)*nx + i];
//            n++;
//          }
//          if (i>1 && j>1) {
//            sum += layer[(j - 1)*nx + i - 1];
//            n++;
//          }
//          if (i<nx-1) {
//            sum += layer[j*nx + i + 1];
//            n++;
//          }
//          if (j<ny-1) {
//            sum += layer[(j + 1)*nx + i];
//            n++;
//          }
//          if (i<nx-1 && j<ny-1) {
//            sum += layer[(j + 1)*nx + i + 1];
//            n++;
//          }
//          value = static_cast<float>(sum)/static_cast<float>(n);
//        }
//        else {
//          value = RMISSING;
//        }
//        pNew->setNextReal(value);
//      }
//    }
//  }
//  pOld->endAccess();
//  pNew->endAccess();
//
//  delete [] layer;
//  delete [] a;
//  delete [] b;
//}

void
Background::resampleParameter(NRLib::Grid<double> & p_new, //FFTGrid *& pNew,        // Resample to
                              NRLib::Grid<double> & p_old, //FFTGrid  * pOld,        // Resample from
                              const Simbox        * simbox_new,
                              const Simbox        * simbox_old)
                              //bool          isFile)
{
  int nx  = simbox_new->getnx();
  int ny  = simbox_new->getny();
  int nz  = simbox_new->getnz();
  //
  // Use same padding as for nonresampled cubes
  //
  // NBNB-PAL: These grids are unpadded, so all nxp, nyp, ... would probably
  //           better be replaced by nx, ny, ... to avoid confusion...
  //
  //int nxp = nx + (pOld->getNxp() - pOld->getNxp());
  //int nyp = ny + (pOld->getNyp() - pOld->getNyp());
  //int nzp = nz + (pOld->getNzp() - pOld->getNzp());

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
      double dzNew = simbox_new->getdz(i,j);
      double dzOld = simbox_old->getdz(i,j);
      double z0New = simbox_new->getTop(i,j);
      double z0Old = simbox_old->getTop(i,j);
        a[ij] = dzNew/dzOld;
      b[ij] = (z0New - z0Old)/dzOld;
      ij++;
    }
  }

  //
  // Resample parameter
  //
  p_new.Resize(nx, ny, nz);
  //pNew = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  //pNew->createRealGrid();
  //pNew->setType(FFTGrid::PARAMETER);
  //pNew->setAccessMode(FFTGrid::WRITE);

  //pOld->setAccessMode(FFTGrid::RANDOMACCESS);

  //int rnxp = 2*(nxp/2 + 1);

  double * layer = new double[nx*ny];

  for(int k=0 ; k<nz; k++) {
    //
    // Map a layer
    //
    int ij=0;
    for (int j=0 ; j<ny; j++) {
      for (int i=0 ; i<nx; i++) {
        if (i < nx && j < ny && k < nz) {
          int k_old = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]);
          layer[ij] = p_old(i, j, k_old);//pOld->getRealValue(i, j, kOld);
          ij++;
        }
      }
    }
    //
    // Smooth the layer (equal weighting of all neighbouring cells)
    //
    float value;
    for (int j=0 ; j<ny; j++) {
      for (int i=0 ; i<nx; i++) {
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
        p_new(i, j, k) = value;
        //pNew->setNextReal(value);
      }
    }
  }
  //pOld->endAccess();
  //pNew->endAccess();

  delete [] layer;
  delete [] a;
  delete [] b;
}

//-------------------------------------------------------------------------------
void
Background::padAndSetBackgroundModel(FFTGrid * bg_vp,
                                     FFTGrid * bg_vs,
                                     FFTGrid * bg_rho)
{
  //LogKit::LogFormatted(LogKit::Low,"\nPadding background model...\n");
  createPaddedParameter(back_model_[0], bg_vp);
  createPaddedParameter(back_model_[1], bg_vs);
  createPaddedParameter(back_model_[2], bg_rho);
}

//-------------------------------------------------------------------------------
//void
//Background::padAndSetBackgroundModelInterval(std::vector<FFTGrid *> bg_vp,
//                                             std::vector<FFTGrid *> bg_vs,
//                                             std::vector<FFTGrid *> bg_rho)
//{
//  //LogKit::LogFormatted(LogKit::Low,"\nPadding background model...\n");
//  int n_intervals = bg_vp.size();
//
//  for(int i = 0; i < n_intervals; i++) {
//    createPaddedParameter(back_model_interval_[i][0], bg_vp[i]);
//    createPaddedParameter(back_model_interval_[i][1], bg_vs[i]);
//    createPaddedParameter(back_model_interval_[i][2], bg_rho[i]);
//  }
//}

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
Background::writeBackgrounds(const Simbox            * simbox,
                             GridMapping             * depthMapping,
                             const GridMapping       * timeMapping,
                             const bool                isFile,
                             const TraceHeaderFormat & thf) const
{
  if(depthMapping != NULL && depthMapping->getSimbox() == NULL) {
    const Simbox * timeSimbox = simbox;
    if(timeMapping != NULL)
      timeSimbox = timeMapping->getSimbox();
    back_model_[0]->setAccessMode(FFTGrid::RANDOMACCESS);
    depthMapping->setMappingFromVelocity(back_model_[0], timeSimbox);
    back_model_[0]->endAccess();
  }

  std::string fileName1 = IO::PrefixBackground() + "Vp" ;
  std::string fileName2 = IO::PrefixBackground() + "Vs" ;
  std::string fileName3 = IO::PrefixBackground() + "Rho";

  FFTGrid * expAlpha = copyFFTGrid(back_model_[0], true, isFile);
  expAlpha->writeFile(fileName1, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expAlpha;

  FFTGrid * expBeta = copyFFTGrid(back_model_[1], true, isFile);
  expBeta->writeFile(fileName2, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expBeta;

  FFTGrid * expRho = copyFFTGrid(back_model_[2], true, isFile);
  expRho->writeFile(fileName3, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expRho;

  //
  // For debugging: write cubes not in ASCII, with padding, and with flat top.
  //
  //back_model_[0]->writeStormFile(fileName1, simbox, true, false, true, true);
  //back_model_[1]->writeStormFile(fileName2, simbox, true, false, true, true);
  //back_model_[2]->writeStormFile(fileName3, simbox, true, false, true, true);
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
