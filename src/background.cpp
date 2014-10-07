/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#ifdef PARALLEL
#include <omp.h>
#endif

#include "lib/kriging1d.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/random/beta.hpp"
#include "nrlib/random/distribution.hpp"

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelsettings.h"
#include "src/vario.h"
#include "src/simbox.h"
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
#include "src/multiintervalgrid.h"

Background::Background(std::vector<NRLib::Grid<float> *>                & parameters,
                       NRLib::Grid<float>                               * velocity,
                       const Simbox                                     * simbox,
                       const Simbox                                     * bg_simbox,
                       const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                       const std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                       const ModelSettings                              * model_settings,
                       std::string                                      & err_text)
  : DataTarget_(250) // For kriging: Increase surrounding until 250 data points is aquired
{

  int nx_pad, ny_pad, nz_pad;

  if (bg_simbox == NULL) {
    nx_pad = simbox->getnx();
    ny_pad = simbox->getny();
    nz_pad = simbox->getnz();
  }
  else {
    nx_pad = bg_simbox->getnx();
    ny_pad = bg_simbox->getny();
    nz_pad = bg_simbox->getnz();
  }

  for (int i=0 ; i<3 ; i++)
    parameters[i]->Resize(nx_pad, ny_pad, nz_pad);

  if (bg_simbox == NULL) {
    GenerateBackgroundModel(parameters[0], parameters[1], parameters[2],
                            velocity,
                            simbox,
                            blocked_logs,
                            model_settings,
                            err_text);
  }
  else {
    GenerateBackgroundModel(parameters[0], parameters[1], parameters[2],
                            velocity,
                            bg_simbox,
                            bg_blocked_logs,
                            model_settings,
                            err_text);

    ResampleBackgroundModel(parameters[0], parameters[1], parameters[2],
                            bg_simbox,
                            simbox);

  }

  //padAndSetBackgroundModel(bg_vp, bg_vs, bg_rho);

  //vsvp_ = findMeanVsVp(bg_vp, bg_vs); //Moved to Commondata
}

//-------------------------------------------------------------------------------
//Background::Background(std::vector<NRLib::Grid<float> *> & parameters,
//                       const std::vector<NRLib::Well>    & wells,
//                       const Simbox                      * interval_simbox,
//                       const ModelSettings               * model_settings,
//                       const std::vector<std::string>    & surface_files,
//                       std::string                       & err_text)
//  : DataTarget_(250) // For kriging: Increase surrounding until 250 data points is aquired
//{
//
//  const int nx    = interval_simbox->getnx();
//  const int ny    = interval_simbox->getny();
//  const int nz    = interval_simbox->getnz();
//
//  for (int i=0 ; i<3 ; i++)
//    parameters[i]->Resize(nx, ny, nz);
//
//  GenerateMultizoneBackgroundModel(parameters[0], //vp
//                                   parameters[1], //vs
//                                   parameters[2], //rho
//                                   wells,
//                                   interval_simbox,
//                                   model_settings,
//                                   surface_files,
//                                   err_text);
//
//  //padAndSetBackgroundModel(bg_vp, bg_vs, bg_rho);
//
//  //vsvp_ = findMeanVsVp(bg_vp, bg_vs); //Moved to Commondata
//}

//-------------------------------------------------------------------------------
//Background::Background(std::vector<std::vector<NRLib::Grid<float> *> > & parameters, //vector(intervals) vector(parameter)
//                       const std::vector<NRLib::Well>                  & wells,
//                       MultiIntervalGrid                               * multiple_interval_grid,
//                       const ModelSettings                             * model_settings,
//                       std::string                                     & err_text)
//  : DataTarget_(250) // For kriging: Increase surrounding until 250 data points is aquired
//{
//
//  //Background for multiple intervals
//  int n_intervals = multiple_interval_grid->GetNIntervals();
//
//  for (int i = 0; i < n_intervals; i++) {
//    parameters[i].resize(3);
//  }
//
//  GenerateMultiIntervalBackgroundModel(parameters,
//                                       wells,
//                                       multiple_interval_grid,
//                                       model_settings,
//                                       err_text);
//
//  //for (int i = 0; i < n_intervals; i++) { //Moved to Commondata
//  //  vs_vp_ratios[i] = findMeanVsVp(bg_vp[i], bg_vs[i]);
//  //}
//
//}

//-------------------------------------------------------------------------------
//Background::Background(FFTGrid ** grids)
//  : DataTarget_(IMISSING),
//    vsvp_(RMISSING)
//{
//  for (int i=0 ; i<3 ; i++)
//    back_model_[i] = grids[i];
//  findMeanVsVp(back_model_[0],
//               back_model_[1]);
//}

//-------------------------------------------------------------------------------
Background::~Background(void)
{
  //for (int i=0 ; i<3 ; i++)
  //  if (back_model_[i] != NULL)
  //    delete back_model_[i];
}

//-------------------------------------------------------------------------------
//void
//Background::releaseGrids()
//{
//  //for (int i=0 ; i<3 ; i++)
//  //  back_model_[i] = NULL;
//}

//-------------------------------------------------------------------------------
void
Background::GenerateBackgroundModel(NRLib::Grid<float>                               * bg_vp,
                                    NRLib::Grid<float>                               * bg_vs,
                                    NRLib::Grid<float>                               * bg_rho,
                                    NRLib::Grid<float>                               * velocity,
                                    const Simbox                                     * simbox,
                                    const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                    const ModelSettings                              * model_settings,
                                    std::string                                      & err_text)
{
  const int   nz      = simbox->getnz();
  const int   n_wells = model_settings->getNumberOfWells();
  const float dz      = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  std::string name_vp  = "Vp";
  std::string name_vs  = "Vs";
  std::string name_rho = "Rho";

  std::vector<std::vector<double> > well_trend_vp(n_wells);
  std::vector<std::vector<double> > well_trend_vs(n_wells);
  std::vector<std::vector<double> > well_trend_rho(n_wells);
  std::vector<std::vector<double> > high_cut_well_trend_vp(n_wells);
  std::vector<std::vector<double> > high_cut_well_trend_vs(n_wells);
  std::vector<std::vector<double> > high_cut_well_trend_rho(n_wells);

  std::string err_text_tmp = "";

  GetWellTrends(well_trend_vp,  high_cut_well_trend_vp,  blocked_logs, nz, name_vp,  err_text_tmp);
  GetWellTrends(well_trend_vs,  high_cut_well_trend_vs,  blocked_logs, nz, name_vs,  err_text_tmp);
  GetWellTrends(well_trend_rho, high_cut_well_trend_rho, blocked_logs, nz, name_rho, err_text_tmp);

  if (err_text_tmp != "") {
    err_text += err_text_tmp;
  }
  else {
    std::vector<double> trend_vp(nz);
    std::vector<double> trend_vs(nz);
    std::vector<double> trend_rho(nz);
    std::vector<double> trend_vel(nz);

    //float * trend_vel  = new float[nz]; // Allocate (for simplicity) although not always needed

    std::vector<double> avg_dev_vp(n_wells);
    std::vector<double> avg_dev_vs(n_wells);
    std::vector<double> avg_dev_rho(n_wells);
    std::vector<double> avg_dev_vel(n_wells);

    //float * avg_dev_vel = new float[n_wells]; // Allocate (for simplicity) although not always needed

    CalculateBackgroundTrend(trend_vp,
                             avg_dev_vp,
                             nz,
                             dz,
                             model_settings->getVpMin(),
                             model_settings->getVpMax(),
                             model_settings->getMaxHzBackground(),
                             well_trend_vp,
                             high_cut_well_trend_vp,
                             name_vp);
    CalculateBackgroundTrend(trend_vs,
                             avg_dev_vs,
                             nz,
                             dz,
                             model_settings->getVsMin(),
                             model_settings->getVsMax(),
                             model_settings->getMaxHzBackground(),
                             well_trend_vs,
                             high_cut_well_trend_vs,
                             name_vs);
    CalculateBackgroundTrend(trend_rho,
                             avg_dev_rho,
                             nz,
                             dz,
                             model_settings->getRhoMin(),
                             model_settings->getRhoMax(),
                             model_settings->getMaxHzBackground(),
                             well_trend_rho,
                             high_cut_well_trend_rho,
                             name_rho);

    bool has_velocity_trend = velocity->GetN() != 0;
    bool write1D            = ((model_settings->getOtherOutputFlag()& IO::BACKGROUND_TREND_1D) > 0);
    bool write3D            = ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);

    WriteTrendsToFile(trend_vp,  simbox, write1D, write3D, has_velocity_trend, name_vp,  model_settings->getFileGrid());
    WriteTrendsToFile(trend_vs,  simbox, write1D, write3D, has_velocity_trend, name_vs,  model_settings->getFileGrid());
    WriteTrendsToFile(trend_rho, simbox, write1D, write3D, has_velocity_trend, name_rho, model_settings->getFileGrid());

    if (velocity->GetN() != 0) {
      //
      // We still want CalculateBackgroundTrend() for alpha above. By calculating
      // avgDevAlpha we can check that the bgAlpha calculated from velocity is as
      // good as or better than that calculated by crava.
      //
      CalculateVelocityDeviations(velocity,
                                  simbox, blocked_logs,
                                  trend_vel, avg_dev_vel, avg_dev_vp,
                                  //model_settings->getOutputGridsElastic(),
                                  n_wells);

      velocity->LogTransform(RMISSING);

      delete bg_vp;
      bg_vp = velocity;
      velocity->Resize(0, 0, 0, 0.0);
      WriteDeviationsFromVerticalTrend(avg_dev_vel, avg_dev_vs, avg_dev_rho,
                                       trend_vel, trend_vs, trend_rho,
                                       blocked_logs, n_wells, nz);
    }
    else {
      WriteDeviationsFromVerticalTrend(avg_dev_vp, avg_dev_vs, avg_dev_rho,
                                       trend_vp, trend_vs, trend_rho,
                                       blocked_logs, n_wells, nz);
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

    std::vector<const std::vector<int> *> ipos(n_wells);
    std::vector<const std::vector<int> *> jpos(n_wells);
    std::vector<const std::vector<int> *> kpos(n_wells);

    for (int i=0; i < n_wells; i++) {
      vt_vp[i]  = std::vector<double>(nz);
      vt_vs[i]  = std::vector<double>(nz);
      vt_rho[i] = std::vector<double>(nz);
    }

    std::vector<int> n_blocks(n_wells);
    int              tot_blocks;

    GetKrigingWellTrends(bl_vp,bl_vs,bl_rho,
                         vt_vp,vt_vs,vt_rho,
                         ipos,jpos,kpos,
                         n_blocks,tot_blocks,
                         blocked_logs, n_wells);

    SetupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
                       trend_vp,trend_vs,trend_rho,
                       model_settings->getOutputGridsElastic(),
                       nz,dz,tot_blocks,n_blocks,
                       bl_vp,bl_vs,bl_rho,
                       vt_vp,vt_vs,vt_rho,
                       ipos,jpos,kpos);

    const CovGrid2D & covGrid2D = MakeCovGrid2D(simbox,
                                                model_settings->getBackgroundVario(),
                                                model_settings->getDebugFlag());

    MakeKrigedBackground(kriging_data_vp,  bg_vp,  trend_vp,  simbox, covGrid2D, "Vp",model_settings->getNumberOfThreads());
    MakeKrigedBackground(kriging_data_vs,  bg_vs,  trend_vs,  simbox, covGrid2D, "Vs",model_settings->getNumberOfThreads());
    MakeKrigedBackground(kriging_data_rho, bg_rho, trend_rho, simbox, covGrid2D, "Rho",model_settings->getNumberOfThreads());

    delete &covGrid2D;
  }
}

//-------------------------------------------------------------------------------
//void
//Background::GenerateMultizoneBackgroundModel(NRLib::Grid<float>             * bg_vp,
//                                             NRLib::Grid<float>             * bg_vs,
//                                             NRLib::Grid<float>             * bg_rho,
//                                             const std::vector<NRLib::Well> & wells,
//                                             const Simbox                   * simbox,
//                                             const ModelSettings            * model_settings,
//                                             const std::vector<std::string> & surface_files,
//                                             std::string                    & err_text)
//{
//  LogKit::LogFormatted(LogKit::Low,"Multizone background model:\n");
//
//  std::vector<int> correlation_structure = model_settings->getCorrelationStructure();
//  std::vector<int> erosion_priority      = model_settings->getErosionPriority();
//
//  int    n_wells    = model_settings->getNumberOfWells();
//  int    n_zones    = static_cast<int>(correlation_structure.size()) - 1;
//  float  dz         = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory
//
//  std::vector<Surface> surface(n_zones+1);
//  for (int i=0; i<n_zones+1; i++)
//    surface[i] = Surface(surface_files[i]);
//
//  std::vector<StormContGrid> vp_zones(n_zones);
//  std::vector<StormContGrid> vs_zones(n_zones);
//  std::vector<StormContGrid> rho_zones(n_zones);
//
//  std::string err_text_tmp = "";
//  BuildSeismicPropertyZones(vp_zones,
//                            vs_zones,
//                            rho_zones,
//                            surface,
//                            correlation_structure,
//                            simbox,
//                            dz,
//                            err_text_tmp);
//
//  if (err_text_tmp != "") {
//    err_text += err_text_tmp;
//  }
//  else {
//
//    std::vector<Surface *> eroded_surfaces(n_zones+1);
//    for (int i=0; i<n_zones+1; i++)
//      eroded_surfaces[i] = NULL;
//
//    ErodeAllSurfaces(eroded_surfaces,
//                     erosion_priority,
//                     surface,
//                     simbox);
//
//    const CovGrid2D & covGrid2D = MakeCovGrid2D(simbox, model_settings->getBackgroundVario(), model_settings->getDebugFlag());
//
//    std::string name_vp  = "Vp";
//    std::string name_vs  = "Vs";
//    std::string name_rho = "Rho";
//
//    std::vector<std::vector<double> > trend_vp_zone(n_zones);
//    std::vector<std::vector<double> > trend_vs_zone(n_zones);
//    std::vector<std::vector<double> > trend_rho_zone(n_zones);
//
//    for (int i=0; i < n_zones; i++) {
//      LogKit::LogFormatted(LogKit::Low,"\nZone%2d:",i+1);
//
//      int nz = static_cast<int>(vp_zones[i].GetNK());
//
//      std::vector<std::vector<double> > well_trend_vp(n_wells);
//      std::vector<std::vector<double> > well_trend_vs(n_wells);
//      std::vector<std::vector<double> > well_trend_rho(n_wells);
//
//      std::vector<std::vector<double> > high_cut_well_trend_vp(n_wells);
//      std::vector<std::vector<double> > high_cut_well_trend_vs(n_wells);
//      std::vector<std::vector<double> > high_cut_well_trend_rho(n_wells);
//
//      StormContGrid eroded_zone;
//
//      BuildErodedZones(eroded_zone,
//                       eroded_surfaces,
//                       nz,
//                       simbox,
//                       i);
//
//      std::vector<bool> hit_zone(n_wells);
//      checkWellHitsZone(hit_zone, wells, eroded_zone, n_wells); //Check blocked logs?
//
//      std::vector<BlockedLogsCommon *> blocked_logs(n_wells);
//
//      std::string err_text_tmp = "";
//
//      getWellTrendsZone(blocked_logs, well_trend_vp,  high_cut_well_trend_vp,  wells, eroded_zone, hit_zone, nz, name_vp,  i, err_text_tmp);
//      getWellTrendsZone(blocked_logs, well_trend_vs,  high_cut_well_trend_vs,  wells, eroded_zone, hit_zone, nz, name_vs,  i, err_text_tmp);
//      getWellTrendsZone(blocked_logs, well_trend_rho, high_cut_well_trend_rho, wells, eroded_zone, hit_zone, nz, name_rho, i, err_text_tmp);
//
//      if (err_text_tmp != "") {
//        err_text += err_text_tmp;
//        break;
//      }
//      else {
//
//        trend_vp_zone[i].resize(nz);
//        trend_vs_zone[i].resize(nz);
//        trend_rho_zone[i].resize(nz);
//
//        std::vector<double> avg_dev_vp_zone(n_wells);
//        std::vector<double> avg_dev_vs_zone(n_wells);
//        std::vector<double> avg_dev_rho_zone(n_wells);
//
//        CalculateBackgroundTrend(trend_vp_zone[i],
//                                 avg_dev_vp_zone,
//                                 nz,
//                                 dz,
//                                 model_settings->getAlphaMin(),
//                                 model_settings->getAlphaMax(),
//                                 model_settings->getMaxHzBackground(),
//                                 well_trend_vp,
//                                 high_cut_well_trend_vp,
//                                 name_vp);
//        CalculateBackgroundTrend(trend_vs_zone[i],
//                                 avg_dev_vs_zone,
//                                 nz,
//                                 dz,
//                                 model_settings->getBetaMin(),
//                                 model_settings->getBetaMax(),
//                                 model_settings->getMaxHzBackground(),
//                                 well_trend_vs,
//                                 high_cut_well_trend_vs,
//                                 name_vs);
//        CalculateBackgroundTrend(trend_rho_zone[i],
//                                 avg_dev_rho_zone,
//                                 nz,
//                                 dz,
//                                 model_settings->getRhoMin(),
//                                 model_settings->getRhoMax(),
//                                 model_settings->getMaxHzBackground(),
//                                 well_trend_rho,
//                                 high_cut_well_trend_rho,
//                                 name_rho);
//
//        WriteDeviationsFromVerticalTrend(avg_dev_vp_zone,
//                                         avg_dev_vs_zone,
//                                         avg_dev_rho_zone,
//                                         trend_vp_zone[i],
//                                         trend_vs_zone[i],
//                                         trend_rho_zone[i],
//                                         wells,
//                                         n_wells,
//                                         nz);
//
//        std::vector<std::vector<double > > bl_vp(n_wells); // bl = blocked logs
//        std::vector<std::vector<double > > bl_vs(n_wells);
//        std::vector<std::vector<double > > bl_rho(n_wells);
//        std::vector<std::vector<double > > vt_vp(n_wells);   // vt = vertical trend
//        std::vector<std::vector<double > > vt_vs(n_wells);
//        std::vector<std::vector<double > > vt_rho(n_wells);
//
//        std::vector<const std::vector<int> > ipos(n_wells);
//        std::vector<const std::vector<int> > jpos(n_wells);
//        std::vector<const std::vector<int> > kpos(n_wells);
//
//        int              totBlocks;
//        std::vector<int> nBlocks(n_wells);
//
//        getKrigingWellTrendsZone(blocked_logs,
//                                 bl_vp,bl_vs,bl_rho,
//                                 vt_vp,vt_vs,vt_rho,
//                                 ipos,jpos,kpos,
//                                 nBlocks,totBlocks,nz);
//
//        std::vector<KrigingData2D> kriging_data_vp(nz);
//        std::vector<KrigingData2D> kriging_data_vs(nz);
//        std::vector<KrigingData2D> kriging_data_rho(nz);
//
//        SetupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
//                           trend_vp_zone[i],trend_vs_zone[i],trend_rho_zone[i],
//                           model_settings->getOutputGridsElastic(),
//                           nz,dz,totBlocks,nBlocks,
//                           bl_vp,bl_vs,bl_rho,
//                           vt_vp,vt_vs,vt_rho,
//                           ipos,jpos,kpos);
//
//        makeKrigedZone(kriging_data_vp,  trend_vp_zone[i],  vp_zones[i],  covGrid2D);
//        makeKrigedZone(kriging_data_vs,  trend_vs_zone[i],  vs_zones[i],  covGrid2D);
//        makeKrigedZone(kriging_data_rho, trend_rho_zone[i], rho_zones[i], covGrid2D);
//
//      }
//    }
//
//    MakeMultizoneBackground(bg_vp, bg_vs, bg_rho,
//                            vp_zones, vs_zones, rho_zones,
//                            simbox,
//                            erosion_priority,
//                            surface,
//                            model_settings->getSurfaceUncertainty(),
//                            "multizone");
//
//
//    bool write3D = ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);
//
//    if (write3D) {
//      writeMultizoneTrendsToFile(trend_vp_zone, trend_vs_zone, trend_rho_zone,
//                                 vp_zones, vs_zones, rho_zones,
//                                 simbox,
//                                 erosion_priority,
//                                 surface,
//                                 model_settings->getSurfaceUncertainty(),
//                                 model_settings->getFileGrid());
//
//    }
//
//    delete &covGrid2D;
//
//    for (int i=0; i<n_zones; i++) {
//      delete eroded_surfaces[i];
//    }
//  }
//
//}


//-------------------------------------------------------------------------------
//void
//Background::GenerateMultiIntervalBackgroundModel(std::vector<std::vector<NRLib::Grid<float> *> > & parameters, //vector(intervals) vector(parameters)
//                                                 const std::vector<NRLib::Well>                  & wells,
//                                                 MultiIntervalGrid                               * multiple_interval_grid,
//                                                 const ModelSettings                             * model_settings,
//                                                 std::string                                     & err_text)
//{
//  LogKit::LogFormatted(LogKit::Low,"MultiInterval background models:\n");
//
//  //std::vector<int> correlation_structure = model_settings->getCorrelationStructure();
//  //std::vector<int> erosion_priority      = model_settings->getErosionPriority();
//
//  int    n_wells     = wells.size();
//  int    n_intervals = multiple_interval_grid->GetNIntervals();
//  //float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory
//
//  //std::vector<Surface> surface(n_zones+1);
//  //for (int i=0; i<n_zones+1; i++)
//  //  surface[i] = Surface(surface_files[i]);
//
//  std::vector<StormContGrid> vp_zones(n_intervals);
//  std::vector<StormContGrid> vs_zones(n_intervals);
//  std::vector<StormContGrid> rho_zones(n_intervals);
//
//  std::string err_text_tmp = "";
//  BuildSeismicPropertyIntervals(vp_zones,
//                                vs_zones,
//                                rho_zones,
//                                multiple_interval_grid,
//                                err_text_tmp);
//
//  if (err_text_tmp != "") {
//    err_text += err_text_tmp;
//  }
//  else {
//
//    std::string name_vp  = "Vp";
//    std::string name_vs  = "Vs";
//    std::string name_rho = "Rho";
//
//    std::vector<std::vector<double> > trend_vp_zone(n_intervals);
//    std::vector<std::vector<double> > trend_vs_zone(n_intervals);
//    std::vector<std::vector<double> > trend_rho_zone(n_intervals);
//
//    std::vector<const NRLib::Surface<double> *> surfaces; //(n_intervals+1);
//
//    const NRLib::Surface<double> & top_surface = multiple_interval_grid->GetIntervalSimbox(0)->GetTopSurface();
//    surfaces.push_back(&top_surface);
//
//    for (int i=0; i < n_intervals; i++) {
//
//      const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i);
//      surfaces.push_back(&simbox->GetBotSurface());
//
//      const CovGrid2D & cov_grid_2D = MakeCovGrid2D(simbox, model_settings->getBackgroundVario(), model_settings->getDebugFlag());
//      float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory
//
//      LogKit::LogFormatted(LogKit::Low,"\nInterval%2d:",i+1);
//
//      int nz = static_cast<int>(vp_zones[i].GetNK());
//
//      std::vector<std::vector<double> > well_trend_vp(n_wells);
//      std::vector<std::vector<double> > well_trend_vs(n_wells);
//      std::vector<std::vector<double> > well_trend_rho(n_wells);
//
//      std::vector<std::vector<double> > high_cut_well_trend_vp(n_wells);
//      std::vector<std::vector<double> > high_cut_well_trend_vs(n_wells);
//      std::vector<std::vector<double> > high_cut_well_trend_rho(n_wells);
//
//      StormContGrid eroded_zone;
//
//      BuildErodedIntervals(eroded_zone,
//                           nz,
//                           simbox);
//
//      std::vector<bool> hit_zone(n_wells);
//      checkWellHitsZone(hit_zone, wells, eroded_zone, n_wells);
//
//      std::vector<BlockedLogsCommon *> blocked_logs(n_wells);
//
//      std::string err_text_tmp = "";
//
//      getWellTrendsZone(blocked_logs, well_trend_vp,  high_cut_well_trend_vp,  wells, eroded_zone, hit_zone, nz, name_vp,  i, err_text_tmp);
//      getWellTrendsZone(blocked_logs, well_trend_vs,  high_cut_well_trend_vs,  wells, eroded_zone, hit_zone, nz, name_vs,  i, err_text_tmp);
//      getWellTrendsZone(blocked_logs, well_trend_rho, high_cut_well_trend_rho, wells, eroded_zone, hit_zone, nz, name_rho, i, err_text_tmp);
//
//      if (err_text_tmp != "") {
//        err_text += err_text_tmp;
//        break;
//      }
//      else {
//
//        trend_vp_zone[i].resize(nz);
//        trend_vs_zone[i].resize(nz);
//        trend_rho_zone[i].resize(nz);
//
//        std::vector<double> avg_dev_vp_zone(n_wells);
//        std::vector<double> avg_dev_vs_zone(n_wells);
//        std::vector<double> avg_dev_rho_zone(n_wells);
//
//        CalculateBackgroundTrend(trend_vp_zone[i],
//                                 avg_dev_vp_zone,
//                                 nz,
//                                 dz,
//                                 model_settings->getAlphaMin(),
//                                 model_settings->getAlphaMax(),
//                                 model_settings->getMaxHzBackground(),
//                                 well_trend_vp,
//                                 high_cut_well_trend_vp,
//                                 name_vp);
//        CalculateBackgroundTrend(trend_vs_zone[i],
//                                 avg_dev_vs_zone,
//                                 nz,
//                                 dz,
//                                 model_settings->getBetaMin(),
//                                 model_settings->getBetaMax(),
//                                 model_settings->getMaxHzBackground(),
//                                 well_trend_vs,
//                                 high_cut_well_trend_vs,
//                                 name_vs);
//        CalculateBackgroundTrend(trend_rho_zone[i],
//                                 avg_dev_rho_zone,
//                                 nz,
//                                 dz,
//                                 model_settings->getRhoMin(),
//                                 model_settings->getRhoMax(),
//                                 model_settings->getMaxHzBackground(),
//                                 well_trend_rho,
//                                 high_cut_well_trend_rho,
//                                 name_rho);
//
//        WriteDeviationsFromVerticalTrend(avg_dev_vp_zone,
//                                         avg_dev_vs_zone,
//                                         avg_dev_rho_zone,
//                                         trend_vp_zone[i],
//                                         trend_vs_zone[i],
//                                         trend_rho_zone[i],
//                                         wells,
//                                         n_wells,
//                                         nz);
//
//        std::vector<std::vector<double > > bl_vp(n_wells); // bl = blocked logs
//        std::vector<std::vector<double > > bl_vs(n_wells);
//        std::vector<std::vector<double > > bl_rho(n_wells);
//        std::vector<std::vector<double > > vt_vp(n_wells);   // vt = vertical trend
//        std::vector<std::vector<double > > vt_vs(n_wells);
//        std::vector<std::vector<double > > vt_rho(n_wells);
//
//        std::vector<const std::vector<int> > ipos(n_wells);
//        std::vector<const std::vector<int> > jpos(n_wells);
//        std::vector<const std::vector<int> > kpos(n_wells);
//
//        int              tot_blocks;
//        std::vector<int> n_blocks(n_wells);
//
//        getKrigingWellTrendsZone(blocked_logs,
//                                 bl_vp,bl_vs,bl_rho,
//                                 vt_vp,vt_vs,vt_rho,
//                                 ipos,jpos,kpos,
//                                 n_blocks,tot_blocks,nz);
//
//        std::vector<KrigingData2D> kriging_data_vp(nz);
//        std::vector<KrigingData2D> kriging_data_vs(nz);
//        std::vector<KrigingData2D> kriging_data_rho(nz);
//
//        SetupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
//                           trend_vp_zone[i],trend_vs_zone[i],trend_rho_zone[i],
//                           model_settings->getOutputGridsElastic(),
//                           nz,dz,tot_blocks,n_blocks,
//                           bl_vp,bl_vs,bl_rho,
//                           vt_vp,vt_vs,vt_rho,
//                           ipos,jpos,kpos);
//
//        makeKrigedZone(kriging_data_vp, trend_vp_zone[i], vp_zones[i], cov_grid_2D);
//        makeKrigedZone(kriging_data_vs, trend_vs_zone[i], vs_zones[i], cov_grid_2D);
//        makeKrigedZone(kriging_data_rho, trend_rho_zone[i], rho_zones[i], cov_grid_2D);
//
//
//        for (int j=0; j<n_wells; j++) {
//          delete blocked_logs[j];
//        }
//
//        delete &cov_grid_2D;
//      }
//    } //Intervals
//
//
//    //H model_settings->getSurfaceUncertainty() only given with multizone background not background created in a multi interval settings
//    // The part with zone_probability is removed.
//    MakeMultiIntervalBackground(parameters,
//                                vp_zones, vs_zones, rho_zones,
//                                multiple_interval_grid,
//                                surfaces,
//                                //model_settings->getSurfaceUncertainty(),
//                                "multiinterval");
//
//
//    bool write3D = ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);
//
//    if (write3D) {
//      writeMultiIntervalTrendsToFile(trend_vp_zone, trend_vs_zone, trend_rho_zone,
//                                     vp_zones, vs_zones, rho_zones,
//                                     multiple_interval_grid,
//                                     surfaces,
//                                     //model_settings->getSurfaceUncertainty(),
//                                     model_settings->getFileGrid());
//
//    }
//
//    //delete &cov_grid_2D;
//
//  }
//
//}

//void
//Background::MakeMultizoneBackground(NRLib::Grid<float>               * bg_vp,
//                                    NRLib::Grid<float>               * bg_vs,
//                                    NRLib::Grid<float>               * bg_rho,
//                                    const std::vector<StormContGrid> & vp_zones,
//                                    const std::vector<StormContGrid> & vs_zones,
//                                    const std::vector<StormContGrid> & rho_zones,
//                                    const Simbox                     * simbox,
//                                    const std::vector<int>           & erosion_priority,
//                                    const std::vector<Surface>       & surface,
//                                    const std::vector<double>        & surface_uncertainty,
//                                    const std::string                & type) const
//{
//
//  std::string text = "\nBuilding "+type+" background:";
//  LogKit::LogFormatted(LogKit::Low,text);
//
//  int n_zones = static_cast<int>(vp_zones.size());
//
//  int nx = simbox->getnx();
//  int ny = simbox->getny();
//  int nz = simbox->getnz();
//
//  float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
//  float nextMonitor = monitorSize;
//  std::cout
//    << "\n  0%       20%       40%       60%       80%      100%"
//    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
//    << "\n  ^";
//
//  bg_vp->Resize(nx, ny, nz);
//  bg_vs->Resize(nx, ny, nz);
//  bg_rho->Resize(nx, ny, nz);
//
//  // Beta distributed uncertainty on each surface
//  // Note that the upper and lower surfaces not are assigned Beta distributions as these have zero uncertainty
//  std::vector<NRLib::Beta> horizon_distributions(n_zones+1);
//  for (int zone=1; zone<n_zones; zone++) {
//    horizon_distributions[zone] = NRLib::Beta(-surface_uncertainty[zone], surface_uncertainty[zone], 2, 2);
//  }
//
//  NRLib::Grid<double> z_surface(nx, ny, n_zones+1);
//  for (int i = 0; i < nx; i++) {
//    for (int j = 0; j < ny; j++) {
//      double x = 0.0;
//      double y = 0.0;
//      simbox->getXYCoord(i, j, x, y);
//
//      for (int k=0; k<n_zones+1; k++)
//        z_surface(i, j, k) = surface[k].GetZ(x,y);
//    }
//  }
//
// for (int k = 0; k < nz; k++) {
//
//    for (int j = 0; j < ny; j++) {
//
//      for (int i = 0; i < nx; i++) {
//
//        double x = 0.0;
//        double y = 0.0;
//
//        simbox->getXYCoord(i, j, x, y);
//
//        // Calculate z directly to decrease computation time
//        double z = z_surface(i,j,0)+(z_surface(i,j,n_zones)-z_surface(i,j,0))*static_cast<double>(k+0.5)/static_cast<double>(nz);
//
//        std::vector<double> z_relative(n_zones+1);
//        for (int zone=0; zone<n_zones+1; zone++)
//          z_relative[zone] = z - z_surface(i,j,zone);
//
//        std::vector<double> zone_probability(n_zones);
//
//        ComputeZoneProbability(z_relative, horizon_distributions, erosion_priority, zone_probability);
//
//        double vp  = 0.0;
//        double vs  = 0.0;
//        double rho = 0.0;
//
//        for (int zone=0; zone<n_zones; zone++) {
//
//          if (zone_probability[zone] > 0) {
//            size_t ind1;
//            size_t ind2;
//            double t;
//
//            vp_zones[zone].FindZInterpolatedIndex(x, y, z, ind1, ind2, t);
//
//            double vp_zone_new  = vp_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//            double vs_zone_new  = vs_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//            double rho_zone_new = rho_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//
//            vp  +=  vp_zone_new  * zone_probability[zone];
//            vs  +=  vs_zone_new  * zone_probability[zone];
//            rho +=  rho_zone_new * zone_probability[zone];
//          }
//        }
//
//        bg_vp->SetValue(i, j, k, static_cast<float>(vp));
//        bg_vs->SetValue(i, j, k, static_cast<float>(vs));
//        bg_rho->SetValue(i, j, k, static_cast<float>(rho));
//      }
//    }
//
//    // Log progress
//    if (k+1 >= static_cast<int>(nextMonitor)) {
//      nextMonitor += monitorSize;
//      std::cout << "^";
//      fflush(stdout);
//    }
//  }
//}

//void
//Background::MakeMultiIntervalBackground(std::vector<std::vector<NRLib::Grid<float> *> > & parameters, //vector(intervals) vector(parameters)
//                                        const std::vector<StormContGrid>                & vp_zones,
//                                        const std::vector<StormContGrid>                & vs_zones,
//                                        const std::vector<StormContGrid>                & rho_zones,
//                                        MultiIntervalGrid                               * multiple_interval_grid,
//                                        std::vector<const NRLib::Surface<double> *>       surfaces,
//                                        //const std::vector<double>                       & surface_uncertainty,
//                                        const std::string                               & type) const
//{
//
//  std::string text = "\nBuilding "+type+" background:";
//  LogKit::LogFormatted(LogKit::Low,text);
//
//  int n_intervals = multiple_interval_grid->GetNIntervals();
//  //const std::vector<int> & erosion_priority = multiple_interval_grid->GetErosionPriorities();
//
//  for (int i_interval = 0; i_interval < n_intervals; i_interval++) {
//
//    const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i_interval);
//
//    int nx = simbox->getnx();
//    int ny = simbox->getny();
//    int nz = simbox->getnz();
//
//    //const int nxp  = nx;
//    const int nyp  = ny;
//    const int nzp  = nz;
//    const int rnxp = nx; //2*(nxp/2 + 1);
//
//    float monitorSize = std::max(1.0f, static_cast<float>(nzp)*0.02f);
//    float nextMonitor = monitorSize;
//    std::cout
//      << "\n  0%       20%       40%       60%       80%      100%"
//      << "\n  |    |    |    |    |    |    |    |    |    |    |  "
//      << "\n  ^";
//
//    parameters[i_interval][0]->Resize(nx, ny, nz);
//    parameters[i_interval][1]->Resize(nx, ny, nz);
//    parameters[i_interval][2]->Resize(nx, ny, nz);
//
//    // Beta distributed uncertainty on each surface
//    // Note that the upper and lower surfaces not are assigned Beta distributions as these have zero uncertainty
//    //std::vector<NRLib::Beta> horizon_distributions(n_intervals+1);
//    //for (int zone=1; zone<n_intervals; zone++) {
//    //  horizon_distributions[zone] = NRLib::Beta(-surface_uncertainty[zone], surface_uncertainty[zone], 2, 2);
//    //}
//
//    NRLib::Grid<double> z_surface(rnxp, nyp, n_intervals+1);
//    for (int i=0; i<rnxp; i++) {
//      for (int j=0; j<nyp; j++) {
//        double x;
//        double y;
//        simbox->getXYCoord(i, j, x, y);
//
//        for (int k=0; k<n_intervals+1; k++)
//          z_surface(i, j, k) = surfaces[k]->GetZ(x,y);
//      }
//    }
//
//   for (int k=0; k<nzp; k++) {
//
//      for (int j=0; j<nyp; j++) {
//
//        for (int i=0; i<rnxp; i++) {
//
//          if (i<nx) {
//
//            double x = 0.0;
//            double y = 0.0;
//
//            simbox->getXYCoord(i, j, x, y);
//
//            // Calculate z directly to decrease computation time
//            double z = z_surface(i,j,0)+(z_surface(i,j,n_intervals)-z_surface(i,j,0))*static_cast<double>(k+0.5)/static_cast<double>(nzp);
//
//            //std::vector<double> z_relative(n_intervals+1);
//            //for (int zone=0; zone<n_intervals+1; zone++)
//            //  z_relative[zone] = z - z_surface(i,j,zone);
//
//            //std::vector<double> zone_probability(n_intervals);
//
//            //ComputeZoneProbability(z_relative, horizon_distributions, erosion_priority, zone_probability);
//
//            double vp  = 0.0;
//            double vs  = 0.0;
//            double rho = 0.0;
//
//            //for (int zone=0; zone<n_intervals; zone++) {
//
//            //  if (zone_probability[zone] > 0) {
//            //    size_t ind1;
//            //    size_t ind2;
//            //    double t;
//
//            //    vp_zones[zone].FindZInterpolatedIndex(x, y, z, ind1, ind2, t);
//
//            //    double vp_zone_new  = vp_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//            //    double vs_zone_new  = vs_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//            //    double rho_zone_new = rho_zones[zone].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//
//            //    vp  +=  vp_zone_new  * zone_probability[zone];
//            //    vs  +=  vs_zone_new  * zone_probability[zone];
//            //    rho +=  rho_zone_new * zone_probability[zone];
//            //  }
//            //}
//
//            //H For multinterval. No weightings between intervals, no zone_probability. Store each interval grid seperatly.
//            size_t ind1;
//            size_t ind2;
//            double t = 0.0;
//
//            vp_zones[i_interval].FindZInterpolatedIndex(x, y, z, ind1, ind2, t);
//            vp  = vp_zones[i_interval].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//            vs  = vs_zones[i_interval].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//            rho = rho_zones[i_interval].GetValueZInterpolatedFromIndexNoMissing(ind1, ind2, t);
//
//            parameters[i_interval][0]->SetValue(i, j, k, static_cast<float>(vp));
//            parameters[i_interval][1]->SetValue(i, j, k, static_cast<float>(vs));
//            parameters[i_interval][2]->SetValue(i, j, k, static_cast<float>(rho));
//
//          }
//          else {
//            parameters[i_interval][0]->SetValue(i, j, k, 0.0f);
//            parameters[i_interval][1]->SetValue(i, j, k, 0.0f);
//            parameters[i_interval][2]->SetValue(i, j, k, 0.0f);
//          }
//        }
//      }
//
//      // Log progress
//      if (k+1 >= static_cast<int>(nextMonitor)) {
//        nextMonitor += monitorSize;
//        std::cout << "^";
//        fflush(stdout);
//      }
//    }
//  }
//}

//---------------------------------------------------------------------------
//void
//Background::ComputeZoneProbability(const std::vector<double>      & z,
//                                   const std::vector<NRLib::Beta> & horizon_distributions,
//                                   const std::vector<int>         & erosion_priority,
//                                   std::vector<double>            & zone_probability) const
//{
//
//  int nZones = static_cast<int>(zone_probability.size());
//
//  std::vector<double> horizon_cdf(nZones+1, 0);
//  horizon_cdf[0] = 1; //The lower surface has cdf 1, whereas the upper surface has cdf 0
//  for (int i=1; i<nZones; i++)
//    horizon_cdf[i] = horizon_distributions[i].Cdf(z[i]);
//
//  for (int zone=0; zone<nZones; zone++) {
//    //Initialize with probability that we are below top surface for zone
//    double prob = horizon_cdf[zone];
//
//    //Multiply with probability that we are above base surface for zone
//    prob *= (1-horizon_cdf[zone+1]);
//
//    //We may be eroded from above. Must consider the surfaces that
//    //1. Are above top in the standard sequence.
//    //2. Have lower erosion priority number than the top.
//    //3. Have no horizons with lower erosion priority number between it and top.
//    int min_erosion = erosion_priority[zone];
//    for (int prev_hor = zone-1; prev_hor >=0; prev_hor--) {
//      if (erosion_priority[prev_hor] < min_erosion) {
//        prob        *= horizon_cdf[prev_hor];
//        min_erosion  = erosion_priority[prev_hor]; //Those with higher number stop in this
//      }
//    }
//
//    //We may be eroded from below. Must consider the surfaces that
//    //1. Are below base in the standard sequence.
//    //2. Have lower erosion priority number than the base.
//    //3. Have no horizons with lower erosion priority number between it and base.
//    min_erosion = erosion_priority[zone+1];
//    for (int late_hor = zone+2; late_hor < nZones+1; late_hor++) {
//      if (erosion_priority[late_hor] < min_erosion) {
//        prob        *= (1-horizon_cdf[late_hor]);
//        min_erosion  = erosion_priority[late_hor]; //Those with higher number stop in this
//      }
//    }
//
//    zone_probability[zone] = prob;
//  }
//}

//---------------------------------------------------------------------------
//void
//Background::BuildErodedZones(StormContGrid                & eroded_zone,
//                             const std::vector<Surface *> & eroded_surfaces,
//                             const int                    & nz,
//                             const Simbox                 * simbox,
//                             const int                    & i) const
//{
//  int    nx        = simbox->getnx();
//  int    ny        = simbox->getny();
//  double x_min     = simbox->GetXMin();
//  double y_min     = simbox->GetYMin();
//  double lx        = simbox->GetLX();
//  double ly        = simbox->GetLY();
//  double angle     = simbox->getAngle();
//
//
//  NRLib::Volume volume(x_min, y_min, lx, ly, *eroded_surfaces[i], *eroded_surfaces[i+1], angle);
//
//  eroded_zone = StormContGrid(volume, nx, ny, nz);
//
//}

//---------------------------------------------------------------------------
//void
//Background::BuildErodedIntervals(StormContGrid & eroded_interval,
//                                 const int     & nz,
//                                 const Simbox  * simbox) const
//{
//  int    nx        = simbox->getnx();
//  int    ny        = simbox->getny();
//  double x_min     = simbox->GetXMin();
//  double y_min     = simbox->GetYMin();
//  double lx        = simbox->GetLX();
//  double ly        = simbox->GetLY();
//  double angle     = simbox->getAngle();
//
//  const NRLib::Surface<double> * top_eroded_surface  = simbox->GetTopErodedSurface();
//  const NRLib::Surface<double> * base_eroded_surface = simbox->GetBaseErodedSurface();
//
//  NRLib::Volume volume(x_min, y_min, lx, ly, *top_eroded_surface, *base_eroded_surface, angle);
//
//  eroded_interval = StormContGrid(volume, nx, ny, nz);
//
//}
//---------------------------------------------------------------------------
//void
//Background::ErodeAllSurfaces(std::vector<Surface *>     & eroded_surfaces,
//                             const std::vector<int>     & erosion_priority,
//                             const std::vector<Surface> & surface,
//                             const Simbox               * simbox) const
//{
//  int    n_surf     = static_cast<int>(eroded_surfaces.size());
//
//  for (int i=0; i<n_surf; i++) {
//    int l=0;
//    while(i+1 != erosion_priority[l])
//      l++;
//
//    Surface * temp_surface = new Surface(surface[l]);
//
//    //Find closest eroded surface downward
//    for (int k=l+1; k<n_surf; k++) {
//      if (eroded_surfaces[k] != NULL) {
//        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, false);
//        break;
//      }
//    }
//    //Find closest eroded surface upward
//    for (int k=l-1; k>=0; k--) {
//      if (eroded_surfaces[k] != NULL) {
//        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, true);
//        break;
//      }
//    }
//    eroded_surfaces[l] = temp_surface;
//  }
//}


//---------------------------------------------------------------------------

//void
//Background::BuildSeismicPropertyZones(std::vector<StormContGrid> & vp_zones,
//                                      std::vector<StormContGrid> & vs_zones,
//                                      std::vector<StormContGrid> & rho_zones,
//                                      const std::vector<Surface> & surface,
//                                      const std::vector<int>     & correlation_structure,
//                                      const Simbox               * simbox,
//                                      const float                & dz,
//                                      std::string                & err_text) const
//{
//  int    n_zones    = static_cast<int>(vp_zones.size());
//  int    nx        = simbox->getnx();
//  int    ny        = simbox->getny();
//  double x_min     = simbox->GetXMin();
//  double y_min     = simbox->GetYMin();
//  double lx        = simbox->GetLX();
//  double ly        = simbox->GetLY();
//  double angle     = simbox->getAngle();
//
//  for (int i=1; i<n_zones+1; i++) {
//    Surface temp_top;
//    Surface temp_base;
//    double  x;
//    double  y;
//    double  z_top;
//    double  z_base;
//
//    Surface top  = surface[i-1];
//    Surface base = surface[i];
//
//    double top_missing  = top.GetMissingValue();
//    double base_missing = base.GetMissingValue();
//
//    //Find maximum distance between the surfaces
//    double max_distance = 0;
//
//    for (int j=0; j<nx; j++) {
//      for (int k=0; k<ny; k++) {
//        simbox->getXYCoord(j,k,x,y);
//
//        z_top  = top.GetZ(x,y);
//        z_base = base.GetZ(x,y);
//
//        if (z_top == top_missing) {
//          const std::string name = top.GetName();
//          err_text += "ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n";
//        }
//        else if (z_base == base_missing) {
//          const std::string name = base.GetName();
//          err_text += "ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n";
//        }
//
//        if (z_base-z_top > max_distance) {
//          if (z_top != top_missing && z_base != base_missing)
//            max_distance = z_base-z_top;
//        }
//      }
//    }
//
//    if (max_distance == 0) {
//      err_text += "ERROR: Zone number "+NRLib::ToString(i)+" has size zero. Check the that surface "+NRLib::ToString(i)+" is above surface "+NRLib::ToString(i+1)+".\n";
//    }
//
//    //Make new top and base surfaces
//    if (correlation_structure[i] == ModelSettings::TOP) {
//      temp_top  = top;
//      temp_base = top;
//      temp_base.Add(max_distance);
//    }
//    else if (correlation_structure[i] == ModelSettings::BASE) {
//      temp_top  = base;
//      temp_top.Subtract(max_distance);
//      temp_base = base;
//    }
//    else {
//      temp_top  = top;
//      temp_base = base;
//    }
//
//    NRLib::Volume volume(x_min, y_min, lx, ly, temp_top, temp_base, angle);
//
//    int nz_zone = static_cast<int>(std::ceil(max_distance/dz));
//
//    vp_zones[i-1]  = StormContGrid(volume, nx, ny, nz_zone);
//    vs_zones[i-1]  = StormContGrid(volume, nx, ny, nz_zone);
//    rho_zones[i-1] = StormContGrid(volume, nx, ny, nz_zone);
//
//  }
//}

//---------------------------------------------------------------------------

//void
//Background::BuildSeismicPropertyIntervals(std::vector<StormContGrid> & vp_zones,
//                                          std::vector<StormContGrid> & vs_zones,
//                                          std::vector<StormContGrid> & rho_zones,
//                                          MultiIntervalGrid          * multiple_interval_grid,
//                                          std::string                & err_text) const
//{
//
//  int    n_intervals = multiple_interval_grid->GetNIntervals();
//
//  for (int i=0; i<n_intervals; i++) {
//
//    const Simbox * simbox = multiple_interval_grid->GetIntervalSimbox(i);
//
//    int    nx        = simbox->getnx();
//    int    ny        = simbox->getny();
//    double x_min     = simbox->GetXMin();
//    double y_min     = simbox->GetYMin();
//    double lx        = simbox->GetLX();
//    double ly        = simbox->GetLY();
//    double angle     = simbox->getAngle();
//
//    float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory
//
//    double  x;
//    double  y;
//    double  z_top;
//    double  z_base;
//
//    const NRLib::Surface<double> & top  = simbox->GetTopSurface();
//    const NRLib::Surface<double> & base = simbox->GetBotSurface();
//
//    double top_missing  = RMISSING; //top.IsMissing GetMissingValue();
//    double base_missing = RMISSING; //base.GetMissingValue();
//
//    //Find maximum distance between the surfaces
//    double max_distance = 0;
//
//    for (int j=0; j<nx; j++) {
//      for (int k=0; k<ny; k++) {
//        simbox->getXYCoord(j,k,x,y);
//
//        z_top  = top.GetZ(x,y);
//        z_base = base.GetZ(x,y);
//
//        if (z_top == top_missing) {
//          err_text += "ERROR: The top surface for interval \'"+multiple_interval_grid->GetIntervalName(i)+"\' does not cover the inversion grid, or it contains missing values.\n";
//        }
//        else if (z_base == base_missing) {
//          err_text += "ERROR: The base surface for interval \'"+multiple_interval_grid->GetIntervalName(i)+"\' does not cover the inversion grid, or it contains missing values.\n";
//        }
//
//        if (z_base-z_top > max_distance) {
//          if (z_top != top_missing && z_base != base_missing)
//            max_distance = z_base-z_top;
//        }
//      }
//    }
//
//    if (max_distance == 0) {
//      err_text += "ERROR: Zone number "+NRLib::ToString(i)+" has size zero. Check the that surface "+NRLib::ToString(i)+" is above surface "+NRLib::ToString(i+1)+".\n";
//    }
//
//    ////Make new top and base surfaces
//    //if (correlation_structure[i] == ModelSettings::TOP) {
//    //  temp_top  = top;
//    //  temp_base = top;
//    //  temp_base.Add(max_distance);
//    //}
//    //else if (correlation_structure[i] == ModelSettings::BASE) {
//    //  temp_top  = base;
//    //  temp_top.Subtract(max_distance);
//    //  temp_base = base;
//    //}
//    //else {
//    //  temp_top  = top;
//    //  temp_base = base;
//    //}
//
//    NRLib::Volume volume(x_min, y_min, lx, ly, top, base, angle);
//
//    int nz_zone = static_cast<int>(std::ceil(max_distance/dz));
//
//    vp_zones[i]  = StormContGrid(volume, nx, ny, nz_zone);
//    vs_zones[i]  = StormContGrid(volume, nx, ny, nz_zone);
//    rho_zones[i] = StormContGrid(volume, nx, ny, nz_zone);
//
//  }
//}

void
Background::CalculateVelocityDeviations(NRLib::Grid<float>                               * velocity,
                                        const Simbox                                     * simbox,
                                        const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                        std::vector<double>                              & trend_vel,
                                        std::vector<double>                              & avg_dev_vel,
                                        std::vector<double>                              & avg_dev_vp,
                                        //int                                                output_flag,
                                        int                                                n_wells)
{
  //H Writing of NRLib::Grid missing.
  //if ((output_flag & IO::BACKGROUND_TREND) > 0) {
  //  std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + "VpFromFile";
  //  velocity->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
  //}

  //
  // Calculate deviation between well data and trend
  //
  int max_blocks = 0;
  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    int n_blocks = iter->second->GetNumberOfBlocks();

    if (n_blocks > max_blocks)
      max_blocks = n_blocks;

    w++;
  }

  std::vector<double> velocity_log(max_blocks);

  const int nz        = simbox->getnz();
  std::vector<double> vt_vp(nz);
  std::vector<double> vt_velocity(nz);

  for (int k=0 ; k<nz ; k++)
    trend_vel[k]=0.0;

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    const std::vector<double> & vp_log = blocked_log->GetVpHighCutBackground();
    blocked_log->GetVerticalTrend(vp_log, vt_vp);
    blocked_log->GetBlockedGrid(velocity, velocity_log);
    blocked_log->GetVerticalTrend(velocity_log, vt_velocity);
    double sumDev = 0.0;
    int count = 0;
    for (int k = 0; k < nz; k++) {
      if (vt_vp[k] != RMISSING) {
        trend_vel[k] += vt_velocity[k];
        double diff = exp(vt_vp[k]) - vt_velocity[k]; // Velocity trend is in exp-domain
        sumDev += diff*diff;
        count++;
      }
    }
    if (count > 0)
      sumDev /= count;
    avg_dev_vel[w] = sqrt(sumDev);

    w++;
  }

  for (int k = 0; k < nz; k++)
    trend_vel[k] /= n_wells;

  LogKit::LogFormatted(LogKit::Low,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
  LogKit::LogFormatted(LogKit::Low,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
  LogKit::LogFormatted(LogKit::Low,"Well             TrendFromFile  TrendFromData\n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------\n");

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f          %5.1f\n",
                         iter->second->GetWellName().c_str(),avg_dev_vel[w],avg_dev_vp[w]);

    w++;
  }
}

//---------------------------------------------------------------------------
void
Background::CalculateBackgroundTrend(std::vector<double>               & trend,
                                     std::vector<double>               & avg_dev,
                                     const int                           nz,
                                     const float                         dz,
                                     float                               log_min,
                                     float                               log_max,
                                     float                               max_hz,
                                     std::vector<std::vector<double> > & well_trend,
                                     std::vector<std::vector<double> > & high_cut_well_trend,
                                     const std::string                 & name)
{

  CalculateVerticalTrend(well_trend,
                         trend,
                         log_min,
                         log_max,
                         max_hz,
                         nz,
                         dz,
                         name);


  CalculateDeviationFromVerticalTrend(high_cut_well_trend, trend, avg_dev, nz);


}

//---------------------------------------------------------------------------
void
Background::GetKrigingWellTrends(std::vector<std::vector<double> >                & bl_vp,
                                 std::vector<std::vector<double> >                & bl_vs,
                                 std::vector<std::vector<double> >                & bl_rho,
                                 std::vector<std::vector<double> >                & vt_vp,
                                 std::vector<std::vector<double> >                & vt_vs,
                                 std::vector<std::vector<double> >                & vt_rho,
                                 std::vector<const std::vector<int> *>            & ipos,
                                 std::vector<const std::vector<int> *>            & jpos,
                                 std::vector<const std::vector<int> *>            & kpos,
                                 std::vector<int>                                 & n_blocks,
                                 int                                              & tot_blocks,
                                 const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                 const int                                        & n_wells) const
{
  int max_blocks = 0;
  tot_blocks     = 0;

  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);

    n_blocks[w] = iter->second->GetNumberOfBlocks();
    tot_blocks += n_blocks[w];
    if (n_blocks[w] > max_blocks)
      max_blocks = n_blocks[w];

    w++;
  }

  for (int i = 0; i < n_wells; i++) {
    bl_vp[i]  = std::vector<double>(max_blocks);
    bl_vs[i]  = std::vector<double>(max_blocks);
    bl_rho[i] = std::vector<double>(max_blocks);
  }

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    bl_vp[w]  = blocked_log->GetVpHighCutBackground();
    bl_vs[w]  = blocked_log->GetVsHighCutBackground();
    bl_rho[w] = blocked_log->GetRhoHighCutBackground();

    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    blocked_log->GetVerticalTrend(bl_vp[w],  vt_vp[w]);
    blocked_log->GetVerticalTrend(bl_vs[w],  vt_vs[w]);
    blocked_log->GetVerticalTrend(bl_rho[w], vt_rho[w]);

    ipos[w] = &(blocked_log->GetIposVector());
    jpos[w] = &(blocked_log->GetJposVector());
    kpos[w] = &(blocked_log->GetKposVector());

    w++;
  }
}

//---------------------------------------------------------------------------
//void
//Background::getKrigingWellTrendsZone(std::vector<BlockedLogsCommon *>     & blocked_logs,
//                                     std::vector<std::vector<double> >    & bl_vp,
//                                     std::vector<std::vector<double> >    & bl_vs,
//                                     std::vector<std::vector<double> >    & bl_rho,
//                                     std::vector<std::vector<double> >    & vt_vp,
//                                     std::vector<std::vector<double> >    & vt_vs,
//                                     std::vector<std::vector<double> >    & vt_rho,
//                                     std::vector<const std::vector<int> > & ipos,
//                                     std::vector<const std::vector<int> > & jpos,
//                                     std::vector<const std::vector<int> > & kpos,
//                                     std::vector<int>                     & n_blocks,
//                                     int                                  & tot_blocks,
//                                     const int                            & nz) const
//{
//  int n_wells    = static_cast<int>(blocked_logs.size());
//  int max_blocks = 0;
//  tot_blocks     = 0;
//
//  for (int w = 0; w < n_wells; w++) {
//    if (blocked_logs[w] != NULL) {
//      n_blocks[w] = blocked_logs[w]->GetNumberOfBlocks();
//      tot_blocks += n_blocks[w];
//      if (n_blocks[w] > max_blocks)
//        max_blocks = n_blocks[w];
//    }
//    else {
//      n_blocks[w] = 0;
//    }
//  }
//
//  for (int w = 0; w < n_wells; w++) {
//    if (blocked_logs[w] != NULL) {
//
//      const std::vector<double> & bl_vp_high_cut  = blocked_logs[w]->GetVpHighCutBackground();
//      const std::vector<double> & bl_vs_high_cut  = blocked_logs[w]->GetVsHighCutBackground();
//      const std::vector<double> & bl_rho_high_cut = blocked_logs[w]->GetRhoHighCutBackground();
//
//      std::vector<double> bl_vp_copy(max_blocks);
//      std::vector<double> bl_vs_copy(max_blocks);
//      std::vector<double> bl_rho_copy(max_blocks);
//
//      for (int i = 0; i < n_blocks[w]; i++) {
//        bl_vp_copy[i]  = bl_vp_high_cut[i];
//        bl_vs_copy[i]  = bl_vs_high_cut[i];
//        bl_rho_copy[i] = bl_rho_high_cut[i];
//      }
//
//      bl_vp[w]  = bl_vp_copy;
//      bl_vs[w]  = bl_vs_copy;
//      bl_rho[w] = bl_rho_copy;
//      //
//      // Extract a one-value-for-each-layer array of blocked logs
//      //
//      vt_vp[w]  = std::vector<double>(nz);
//      vt_vs[w]  = std::vector<double>(nz);
//      vt_rho[w] = std::vector<double>(nz);
//
//      blocked_logs[w]->GetVerticalTrend(bl_vp_high_cut,  vt_vp[w]);
//      blocked_logs[w]->GetVerticalTrend(bl_vs_high_cut,  vt_vs[w]);
//      blocked_logs[w]->GetVerticalTrend(bl_rho_high_cut, vt_rho[w]);
//
//      ipos[w] = blocked_logs[w]->GetIposVector();
//      jpos[w] = blocked_logs[w]->GetJposVector();
//      kpos[w] = blocked_logs[w]->GetKposVector();
//
//    }
//    else {
//      vt_vp[w]  = std::vector<double>(0);
//      vt_vs[w]  = std::vector<double>(0);
//      vt_rho[w] = std::vector<double>(0);
//
//      bl_vp[w]  = std::vector<double>(0);
//      bl_vs[w]  = std::vector<double>(0);
//      bl_rho[w] = std::vector<double>(0);
//    }
//  }
//}

//---------------------------------------------------------------------------
void
Background::GetWellTrends(std::vector<std::vector<double> >                & well_trend,
                          std::vector<std::vector<double> >                & high_cut_well_trend,
                          const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                          const int                                        & nz,
                          const std::string                                & name,
                          std::string                                      & err_text) const
{
  //int n_wells = static_cast<int>(well_trend.size());
  int i_wells = 0;

  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log->GetUseForBackgroundTrend()) {

      if (blocked_log != NULL) {
        well_trend[w].resize(nz);

        if (name == "Vp")
          blocked_log->GetVerticalTrend(blocked_log->GetVpBlocked(), well_trend[w]);
        else if (name == "Vs")
          blocked_log->GetVerticalTrend(blocked_log->GetVsBlocked(), well_trend[w]);
        else if (name == "Rho")
          blocked_log->GetVerticalTrend(blocked_log->GetRhoBlocked(), well_trend[w]);
        else {
          err_text += "ERROR in Background::GetWellTrends(): ";
          err_text += "Log \'"+name+"\' requested, but no such log exists.\n";
        }
        i_wells++;
      }
      else
        well_trend[w].resize(0);
    }
    else
      well_trend[w].resize(0);

    w++;
  }
  if (i_wells == 0) {
    err_text += "\nERROR in Background::GetWellTrends(): There are no wells\n";
    err_text += "available for the estimation of background trend.\n";
  }

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log != NULL) {
      high_cut_well_trend[w].resize(nz);
      if (name == "Vp")
        blocked_log->GetVerticalTrend(blocked_log->GetVpHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Vs")
        blocked_log->GetVerticalTrend(blocked_log->GetVsHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Rho")
        blocked_log->GetVerticalTrend(blocked_log->GetRhoHighCutBackground(), high_cut_well_trend[w]);
      else {
        err_text += "ERROR in Background::GetWellTrends(): ";
        err_text += "Log \'"+name+"\' requested, but no such log exists.\n";
      }
    }
    else
      high_cut_well_trend[w].resize(0);

    w++;
  }
}

//---------------------------------------------------------------------------
//void
//Background::getWellTrendsZone(std::vector<BlockedLogsCommon *>  & bl,
//                              std::vector<std::vector<double> > & well_trend,
//                              std::vector<std::vector<double> > & high_cut_well_trend,
//                              const std::vector<NRLib::Well>    & wells,
//                              StormContGrid                     & eroded_zone,
//                              const std::vector<bool>           & hit_zone,
//                              const int                         & nz,
//                              const std::string                 & name,
//                              const int                         & i,
//                              std::string                       & err_text) const
//{
//
//  int n_valid_wells_in_zone = 0;
//  int n_wells               = static_cast<int>(bl.size());
//
//  std::vector<bool> use_for_background(n_wells);
//
//  for (int w = 0; w < n_wells; w++) {
//    if (hit_zone[w] == true) {
//
//      bl[w] = new BlockedLogsCommon(&wells[w], eroded_zone);
//      n_valid_wells_in_zone++;
//    }
//    else
//      bl[w] = NULL;
//
//    if (bl[w]->GetUseForBackgroundTrend())
//      use_for_background[w] = true;
//    else
//      use_for_background[w] = false;
//  }
//
//  if (n_valid_wells_in_zone == 0) {
//    err_text += "Invalid multizone background estimation: No well hits zone number "+NRLib::ToString(i+1)+"\n";
//  }
//
//  int i_wells = 0;
//
//  for (int w = 0; w < n_wells; w++) {
//    if (use_for_background[w] == true) {
//      if (bl[w] != NULL) {
//        well_trend[w].resize(nz);
//
//        if (name == "Vp")
//          bl[w]->GetVerticalTrend(bl[w]->GetVpBlocked(), well_trend[w]);
//        else if (name == "Vs")
//          bl[w]->GetVerticalTrend(bl[w]->GetVsBlocked(), well_trend[w]);
//        else if (name == "Rho")
//          bl[w]->GetVerticalTrend(bl[w]->GetRhoBlocked(), well_trend[w]);
//        else {
//          err_text += "ERROR in Background::getWellTrendsZone(): ";
//          err_text += "Log \'"+name+"\' requested, but no such log exists.\n";
//        }
//        i_wells++;
//      }
//      else well_trend[w].resize(0);
//    }
//    else
//      well_trend[w].resize(0);
//  }
//  if (i_wells == 0) {
//    err_text += "\nERROR in Background::getWellTrendsZone(): There are no wells\n";
//    err_text += "available for the estimation of background trend.\n";
//  }
//
//  for (int w = 0; w < n_wells; w++) {
//    if (bl[w] != NULL) {
//      high_cut_well_trend[w].resize(nz);
//      if (name == "Vp")
//        bl[w]->GetVerticalTrend(bl[w]->GetVpHighCutBackground(), high_cut_well_trend[w]);
//      else if (name == "Vs")
//        bl[w]->GetVerticalTrend(bl[w]->GetVsHighCutBackground(), high_cut_well_trend[w]);
//      else if (name == "Rho")
//        bl[w]->GetVerticalTrend(bl[w]->GetRhoHighCutBackground(), high_cut_well_trend[w]);
//      else {
//        err_text += "ERROR in Background::getWellTrendsZone(): ";
//        err_text += "Log \'"+name+"\' requested, but no such log exists.\n";
//      }
//    }
//    else
//      high_cut_well_trend[w].resize(0);
//  }
//}

//---------------------------------------------------------------------------
//void
//Background::checkWellHitsZone(std::vector<bool>              & hitZone,
//                              const std::vector<NRLib::Well> & wells,
//                              StormContGrid                  & eroded_zone,
//                              const int                      & nWells) const
//{
//  for (int w=0; w<nWells; w++) {
//    if (wells[w].CheckStormgrid(eroded_zone) == 0) {
//      hitZone[w] = true;
//    }
//    else
//      hitZone[w] = false;
//  }
//}
//---------------------------------------------------------------------------
void
Background::WriteTrendsToFile(std::vector<double> & trend,
                              const Simbox        * simbox,
                              bool                  write1D,
                              bool                  write3D,
                              bool                  has_velocity_trend,
                              const std::string   & name,
                              bool                  is_file)
{
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());
  const int   nz = simbox->getnz();

  if (write1D == true) {
    WriteVerticalTrend(trend, dz, nz, name);
  }

  if (write3D == true && !(name=="Vp" && has_velocity_trend)) {
    const int nx = simbox->getnx();
    const int ny = simbox->getny();
    FFTGrid * trend_grid = ModelGeneral::CreateFFTGrid(nx, ny, nz, nx, ny, nz, is_file);
    FillInVerticalTrend(trend_grid, trend);
    FFTGrid * exp_trend = CopyFFTGrid(trend_grid, true, is_file);
    delete trend_grid;
    std::string file_name = IO::PrefixBackground() + IO::PrefixTrend() + name;
    exp_trend->writeFile(file_name, IO::PathToBackground(), simbox);
    delete exp_trend;
  }
}
//---------------------------------------------------------------------------
//void
//Background::writeMultizoneTrendsToFile(const std::vector<std::vector<double> > & vp_zones,
//                                       const std::vector<std::vector<double> > & vs_zones,
//                                       const std::vector<std::vector<double> > & rho_zones,
//                                       std::vector<StormContGrid>              & vp_trend_zone,
//                                       std::vector<StormContGrid>              & vs_trend_zone,
//                                       std::vector<StormContGrid>              & rho_trend_zone,
//                                       const Simbox                            * simbox,
//                                       const std::vector<int>                  & erosion_priority,
//                                       const std::vector<Surface>              & surface,
//                                       const std::vector<double>               & surface_uncertainty,
//                                       const bool                                isFile) const
//{
//  int n_zones = static_cast<int>(vp_zones.size());
//
//  for (int i=0; i<n_zones; i++) {
//    makeTrendZone(vp_zones[i],  vp_trend_zone[i]);
//    makeTrendZone(vs_zones[i],  vs_trend_zone[i]);
//    makeTrendZone(rho_zones[i], rho_trend_zone[i]);
//  }
//
//  NRLib::Grid<float> * trend_vp  = new NRLib::Grid<float>();
//  NRLib::Grid<float> * trend_vs  = new NRLib::Grid<float>();
//  NRLib::Grid<float> * trend_rho = new NRLib::Grid<float>();
//
//  MakeMultizoneBackground(trend_vp,
//                          trend_vs,
//                          trend_rho,
//                          vp_trend_zone,
//                          vs_trend_zone,
//                          rho_trend_zone,
//                          simbox,
//                          erosion_priority,
//                          surface,
//                          surface_uncertainty,
//                          "trend in multizone");
//
//  //FFTGrid * exp_trend_vp = CopyFFTGrid(trend_vp, true, isFile);
//  //FFTGrid * exp_trend_vs = CopyFFTGrid(trend_vs, true, isFile);
//  //FFTGrid * exp_trend_rho = CopyFFTGrid(trend_rho, true, isFile);
//
//  //std::string file_name_vp = IO::PrefixBackground() + IO::PrefixTrend() + "Vp";
//  //std::string file_name_vs = IO::PrefixBackground() + IO::PrefixTrend() + "Vs";
//  //std::string file_name_rho = IO::PrefixBackground() + IO::PrefixTrend() + "Rho";
//
//  //exp_trend_vp->writeFile(file_name_vp, IO::PathToBackground(), simbox);
//  //exp_trend_vs->writeFile(file_name_vs, IO::PathToBackground(), simbox);
//  //exp_trend_rho->writeFile(file_name_rho, IO::PathToBackground(), simbox);
//
//  //delete exp_trend_vp;
//  //delete exp_trend_vs;
//  //delete exp_trend_rho;
//
//  //delete trend_vp;
//  //delete trend_vs;
//  //delete trend_rho;
//}

//---------------------------------------------------------------------------
//void
//Background::writeMultiIntervalTrendsToFile(const std::vector<std::vector<double> >   & vp_zones,
//                                           const std::vector<std::vector<double> >   & vs_zones,
//                                           const std::vector<std::vector<double> >   & rho_zones,
//                                           std::vector<StormContGrid>                & vp_trend_zone,
//                                           std::vector<StormContGrid>                & vs_trend_zone,
//                                           std::vector<StormContGrid>                & rho_trend_zone,
//                                           MultiIntervalGrid                         * multiple_interval_grid,
//                                           std::vector<const NRLib::Surface<double> *> surfaces,
//                                           //const std::vector<double>                 & surface_uncertainty,
//                                           const bool                                  is_file) const
//{
//  int n_intervals = multiple_interval_grid->GetNIntervals();
//
//  for (int i=0; i<n_intervals; i++) {
//    makeTrendZone(vp_zones[i],  vp_trend_zone[i]);
//    makeTrendZone(vs_zones[i],  vs_trend_zone[i]);
//    makeTrendZone(rho_zones[i], rho_trend_zone[i]);
//  }
//
//  std::vector<std::vector<NRLib::Grid<float> *> > trend_parameters(n_intervals);
//  for (int i = 0; i < n_intervals; i++)
//    trend_parameters[i].resize(3);
//
//  MakeMultiIntervalBackground(trend_parameters,
//                              vp_trend_zone,
//                              vs_trend_zone,
//                              rho_trend_zone,
//                              multiple_interval_grid,
//                              surfaces,
//                              //surface_uncertainty,
//                              //is_file,
//                              "trend in multiinterval");
//
//  //for (int i=0; i < n_intervals; i++) {
//
//  //  std::string interval_name = multiple_interval_grid->GetIntervalName(i);
//  //  const Simbox * simbox     = multiple_interval_grid->GetIntervalSimbox(i);
//
//    //FFTGrid * exp_trend_vp = CopyFFTGrid(trend_vp[i], true, is_file);
//    //FFTGrid * exp_trend_vs = CopyFFTGrid(trend_vs[i], true, is_file);
//    //FFTGrid * exp_trend_rho = CopyFFTGrid(trend_rho[i], true, is_file);
//
//    //std::string file_name_vp = "Interval_" + interval_name + "_" + IO::PrefixBackground() + IO::PrefixTrend() + "Vp";
//    //std::string file_name_vs = "Interval_" + interval_name + "_" + IO::PrefixBackground() + IO::PrefixTrend() + "Vs";
//    //std::string file_name_rho = "Interval_" + interval_name + "_" + IO::PrefixBackground() + IO::PrefixTrend() + "Rho";
//
//    //exp_trend_vp->writeFile(file_name_vp, IO::PathToBackground(), simbox);
//    //exp_trend_vs->writeFile(file_name_vs, IO::PathToBackground(), simbox);
//    //exp_trend_rho->writeFile(file_name_rho, IO::PathToBackground(), simbox);
//
//    //delete exp_trend_vp;
//    //delete exp_trend_vs;
//    //delete exp_trend_rho;
//
//    //delete trend_vp[i];
//    //delete trend_vs[i];
//    //delete trend_rho[i];
//
//  //}
//}

//-------------------------------------------------------------------------------
void
Background::SetupKrigingData2D(std::vector<KrigingData2D>                & kriging_data_vp,
                               std::vector<KrigingData2D>                & kriging_data_vs,
                               std::vector<KrigingData2D>                & kriging_data_rho,
                               std::vector<double>                       & trend_vp,
                               std::vector<double>                       & trend_vs,
                               std::vector<double>                       & trend_rho,
                               const int                                   output_flag,
                               const int                                 & nz,
                               const float                               & dz,
                               const int                                 & tot_blocks,
                               const std::vector<int>                    & n_blocks,
                               const std::vector<std::vector<double> >   & bl_vp,
                               const std::vector<std::vector<double> >   & bl_vs,
                               const std::vector<std::vector<double> >   & bl_rho,
                               const std::vector<std::vector<double> >   & vt_vp,
                               const std::vector<std::vector<double> >   & vt_vs,
                               const std::vector<std::vector<double> >   & vt_rho,
                               const std::vector<const std::vector<int> *> ipos,
                               const std::vector<const std::vector<int> *> jpos,
                               const std::vector<const std::vector<int> *> kpos) const
{
  //
  // Although unnecessary, we have chosen to set up kriging data from
  // Vp, Vs and Rho simultaneously. This gives code easier to read.
  //
  const int n_wells = static_cast<int>(bl_vp.size());

  KrigingData3D forLogging(tot_blocks);

  for (int w = 0; w < n_wells; w++) {

    if (vt_vp[w].size() > 0) {
      std::vector<double> vt_vp_well  = vt_vp[w];
      std::vector<double> vt_vs_well  = vt_vs[w];
      std::vector<double> vt_rho_well = vt_rho[w];

      std::vector<double> bl_vp_well  = bl_vp[w];
      std::vector<double> bl_vs_well  = bl_vs[w];
      std::vector<double> bl_rho_well = bl_rho[w];

      //
      // Kriging vertical trend (vt....) against global vertical trend (trend...)
      //
      Kriging1D::krigVector(&vt_vp_well[0],  &trend_vp[0],  nz, dz);
      Kriging1D::krigVector(&vt_vs_well[0],  &trend_vs[0],  nz, dz);
      Kriging1D::krigVector(&vt_rho_well[0], &trend_rho[0], nz, dz);
      //
      // Use kriged vertical trend where original log is not defined.
      //
      const std::vector<int> ipos_well = *(ipos[w]);
      const std::vector<int> jpos_well = *(jpos[w]);
      const std::vector<int> kpos_well = *(kpos[w]);

      for (int m = 0; m < n_blocks[w]; m++) {
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

      forLogging.addData(bl_vp_well, bl_vs_well, bl_rho_well,
                         ipos_well, jpos_well, kpos_well,
                         n_blocks[w]);
    }

    for (int k = 0; k < nz; k++) {
      kriging_data_vp[k].findMeanValues();
      kriging_data_vs[k].findMeanValues();
      kriging_data_rho[k].findMeanValues();
    }
  }

  if ((output_flag & IO::BACKGROUND) > 0) {
    forLogging.divide();
    std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + IO::SuffixGeneralData();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    forLogging.writeToFile(fileName);
  }
}

//---------------------------------------------------------------------------
const CovGrid2D &
Background::MakeCovGrid2D(const Simbox * simbox,
                          Vario        * vario,
                          int            debug_flag)
{
  //
  // Pretabulate all needed covariances
  //
  const int    nx = simbox->getnx();
  const int    ny = simbox->getny();

  const float  dx = static_cast<float>(simbox->getdx());
  const float  dy = static_cast<float>(simbox->getdy());

  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);

  if (debug_flag == 1) {
    std::string base_name = IO::PrefixBackground() + "covGrid2D" + IO::SuffixAsciiIrapClassic();
    std::string file_name = IO::makeFullFileName(IO::PathToBackground(), base_name);
    cov->writeToFile(file_name);
  }
  return (*cov);
}

void
Background::MakeKrigedBackground(const std::vector<KrigingData2D> & kriging_data,
                                 NRLib::Grid<float>               * bg_grid,
                                 std::vector<double>              & trend,
                                 const Simbox                     * simbox,
                                 const CovGrid2D                  & cov_grid_2D,
                                 const std::string                & type,
                                 int                                n_threads) const
{
  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  const int    nx   = simbox->getnx();
  const int    ny   = simbox->getny();
  const int    nz   = simbox->getnz();

  const double x0   = simbox->getx0();
  const double y0   = simbox->gety0();
  const double lx   = simbox->getlx();
  const double ly   = simbox->getly();

  //
  // Store a surface for each layer (needed for parallelization)
  //
  Surface tmp(x0, y0, lx, ly, nx, ny, RMISSING);
  std::vector<Surface> surfaces(0);
  surfaces.reserve(nz);
  for (int k=0 ; k<nz ; k++)
    surfaces.push_back(tmp);

  float monitor_size = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float next_monitor = monitor_size;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  bg_grid->Resize(nx, ny, nz);
#ifdef PARALLEL
  int  chunk_size = 1;
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif

  for (int k=0 ; k<nz ; k++) {
    // Set trend for layer
    surfaces[k].Assign(trend[k]);

    // Kriging of layer
    Kriging2D::krigSurface(surfaces[k], kriging_data[k], cov_grid_2D);


    // Log progress
    if (k+1 >= static_cast<int>(next_monitor)) {
      next_monitor += monitor_size;
      std::cout << "^";
      fflush(stdout);
    }
  }

  // Set layers in background model from surface
  for (int k=0 ; k<nz ; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0 ; i < nx ; i++) {
        bg_grid->SetValue(i, j, k, static_cast<float>(surfaces[k](i,j)));
      }
    }
  }
}

//---------------------------------------------------------------------------
//void
//Background::makeTrendZone(const std::vector<double> & trend,
//                          StormContGrid             & trend_zone) const
//{
//  const size_t nx   = trend_zone.GetNI();
//  const size_t ny   = trend_zone.GetNJ();
//  const size_t nz   = trend_zone.GetNK();
//
//  const double x0   = trend_zone.GetXMin();
//  const double y0   = trend_zone.GetYMin();
//  const double lx   = trend_zone.GetLX();
//  const double ly   = trend_zone.GetLY();
//  //
//  // Template surface to be kriged
//  //
//  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);
//  for (size_t k = 0; k < nz; k++) {
//
//    // Set trend for layer
//    surface.Assign(trend[k]);
//
//    // Set layer in background model from surface
//    for (size_t j = 0; j < ny; j++) {
//      for (size_t i = 0; i < nx; i++)
//        trend_zone(i,j,k) = float(surface(i,j));
//    }
//  }
//}

//---------------------------------------------------------------------------
//void
//Background::makeKrigedZone(const std::vector<KrigingData2D> & kriging_data,
//                           const std::vector<double>        & trend,
//                           StormContGrid                    & kriged_zone,
//                           const CovGrid2D                  & covGrid2D) const
//{
//  const size_t nx   = kriged_zone.GetNI();
//  const size_t ny   = kriged_zone.GetNJ();
//  const size_t nz   = kriged_zone.GetNK();
//
//  const double x0   = kriged_zone.GetXMin();
//  const double y0   = kriged_zone.GetYMin();
//  const double lx   = kriged_zone.GetLX();
//  const double ly   = kriged_zone.GetLY();
//  //
//  // Template surface to be kriged
//  //
//  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);
//  for (size_t k = 0; k < nz; k++) {
//
//    // Set trend for layer
//    surface.Assign(trend[k]);
//
//    // Kriging of layer
//    Kriging2D::krigSurface(surface, kriging_data[k], covGrid2D);
//
//    // Set layer in background model from surface
//    for (size_t j = 0; j < ny; j++) {
//      for (size_t i = 0; i < nx; i++)
//        kriged_zone(i,j,k) = float(surface(i,j));
//    }
//  }
//}

//-------------------------------------------------------------------------------
void
Background::CalculateVerticalTrend(std::vector<std::vector<double> > & well_trend,
                                   std::vector<double>               & trend,
                                   float                               log_min,
                                   float                               log_max,
                                   float                               max_hz,
                                   int                                 nz,
                                   float                               dz,
                                   const std::string                 & name)
{
  int     n_wells      = static_cast<int>(well_trend.size());
  std::vector<double> filtered_log(nz);
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
  for (int k = 0; k < nz; k++) {
    trend[k] = 0.0f;
    count[k] = 0;
  }
  int i_wells = 0;
  for (int w = 0; w < n_wells; w++) {
    if (well_trend[w].size() > 0) {
      std::vector<double> & w_trend = well_trend[w];

      for (int k = 0; k < nz; k++) {
        if (w_trend[k] != RMISSING) {
          trend[k] += exp(w_trend[k]);
          count[k]++;
        }
      }
      i_wells++;
    }
  }
  if (i_wells > 0) {
    for (int k = 0; k < nz; k++) {
      if (count[k] > 0) {
        trend[k] = trend[k]/static_cast<float>(count[k]);
      }
    }
  }

  //Utils::writeVectorToFile(std::string("trend_mean_values_") + name, trend, nz);

  SmoothTrendWithLocalLinearRegression(trend, count,
                                       i_wells, nz, dz,
                                       log_min,
                                       log_max,
                                       name);

  //Utils::writeVectorToFile(std::string("trend_after_linreg_") + name, trend, nz);

  CommonData::ApplyFilter(filtered_log,
                          trend,
                          nz,
                          dz,
                          max_hz);

  for (int i = 0; i < nz; i++) {
    trend[i] = filtered_log[i];
  }

  //Utils::writeVectorToFile(std::string("trend_after_filter_") + name, trend, nz);

  delete [] count;
}

//-------------------------------------------------------------------------------
void
Background::SmoothTrendWithLocalLinearRegression(std::vector<double> & trend,
                                                 int                 * count,
                                                 int                   i_wells,
                                                 int                   nz,
                                                 float                 dz,
                                                 float                 min_value,
                                                 float                 max_value,
                                                 std::string           par_name)
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

  float fraction     = 5.0f;                      // Require minimum 5*iWells
  int   n_time_limit = static_cast<int>(50.0/dz); // The smaller sampling density, the more values are needed.
  int   n_low_limit  = 10;                        // Require minimum 10
  int   n_data_min   = std::max(n_low_limit, std::max(n_time_limit, int(fraction * i_wells)));

  bool  use_weights  = true;
  bool  error_mid    = false;
  bool  error_head   = false;
  bool  error_trail  = false;

  //
  // Copy the average values (stored in array 'trend') to the array 'mean'.
  //
  float * mean = new float[nz];
  for (int k = 0; k < nz; k++) {
    mean[k] = static_cast<float>(trend[k]);
  }

  //
  // Find first non-missing value
  //
  int first_nonmissing = 0;
  for (int k = 0; k < nz; k++) {
    if (trend[k] > 0.0f) {
      first_nonmissing = k;
      break;
    }
  }

  //
  // Find last non-missing value
  //
  int last_nonmissing = nz - 1;
  for (int k = nz - 1; k > 0; k--) {
    if (trend[k] > 0.0f) {
      last_nonmissing = k;
      break;
    }
  }

  float * x = new float[nz];  // Time indices
  float * y = new float[nz];  // Log values
  float * w = new float[nz];  // Weights (number of data behind each avg.)

  for (int k = 0; k < nz; k++) {
    int n_cur_data_min = n_data_min;
    if (k < first_nonmissing || k > last_nonmissing) {
      n_cur_data_min *= 2;
    }

    int n      = 0;
    int n_data = 0;
    if (debug)
      LogKit::LogFormatted(LogKit::Low,"k=%d\n",k);
    //
    // 1. Add current data point to local data set if present.
    //
    if (count[k] > 0) {
      w[0]   = static_cast<float>(count[k]);
      x[0]   = static_cast<float>(k);
      y[0]   = static_cast<float>(trend[k]);
      n_data += count[k];
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   A:t=%.2f   x[0] y[0]  %d   %.2f\n",dz*(x[0] + 0.5f),int(x[0]),y[0]);
      n++;
    }

    //
    // 2. Add local data points to get 'nCurDataMin' points behind each trend
    //    value. Note that the bandwidth varies
    //
    int i = 0;
    while (n_data < n_cur_data_min) {
      i++;
      if (k - i >= 0 && count[k - i] > 0) {
        w[n]    = static_cast<float>(count[k - i]);
        x[n]    = static_cast<float>(k - i);
        y[n]    = mean [k - i];
        n_data += count[k - i];
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   B:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k + i < nz  && count[k + i] > 0) {
        w[n]    = static_cast<float>(count[k + i]);
        x[n]    = static_cast<float>(k + i);
        y[n]    = mean [k + i];
        n_data += count[k + i];
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
      for (i = 0; i < n; i++)
        w[i] /= n_data;
    else
      for (i = 0; i < n; i++)
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
        if (k < first_nonmissing)
          error_head = true;
        else if (k > last_nonmissing)
          error_trail = true;
        else {
          error_mid   = true;
          break;
        }
      }
      trend[k] = log(static_cast<double>(value));
    }
    else {
      trend[k] = log(static_cast<double>(y[0]));
    }
    if (debug)
      LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f        (minLog/maxLog = %.2f / %.2f)\n",exp(trend[k]),min_value,max_value);
  }

  if (error_mid) {
    // Big problem ...
    LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+par_name+" using local linear\n");
    LogKit::LogFormatted(LogKit::Low,"          regression failed - trying global mean instead. Possible causes: \n");
    LogKit::LogFormatted(LogKit::Low,"          1) Available logs cover too small a part of inversion grid giving extrapolation problems.\n");
    LogKit::LogFormatted(LogKit::Low,"          2) There are too many layers in grid compared to well logs available.\n");
    float sum  = 0.0f;
    int n_data = 0;
    for (int k = 0; k < nz; k++) {
      if (count[k] > 0) {
        if (use_weights) {
          sum    += mean[k]*count[k];
          n_data += count[k];
          if (debug)
            LogKit::LogFormatted(LogKit::Low,"k=%d  count[k], mean[k]  nData, sum  %d  %8.3f     %d  %8.3f\n",
                                 k,count[k],mean[k],n_data,sum);
        }
        else {
          sum    += mean[k];
          n_data += 1;
        }
      }
    }
    float global_mean = sum/n_data;
    for (int k = 0; k < nz; k++) {
      trend[k] = log(global_mean);
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   TREND: k = %d   trend[k] = %.2f\n",k,exp(trend[k]));
    }
    LogKit::LogFormatted(LogKit::Low,"\nGlobal mean for parameter %s = %.2f\n\n",par_name.c_str(),global_mean);
  }
  else {
    if (error_head) {
      // Fix first part of trend containing missing-values.
      double first_value = trend[first_nonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter \'"+par_name+"\' using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [0-%d] where the log is undefined. The first\n",first_nonmissing-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(first_value),first_nonmissing);
      for (int k = 0 ; k < first_nonmissing ; k++) {
        trend[k] = first_value;
      }
    }
    if (error_trail) {
      // Fix last part of trend containing missing-values.
      double last_value = trend[last_nonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+par_name+" using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [%d,%d] where the log is undefined. The last\n",last_nonmissing+1,nz-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(last_value),last_nonmissing);
      for (int k = last_nonmissing + 1; k < nz; k++) {
        trend[k] = last_value;
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
Background::WriteVerticalTrend(std::vector<double> & trend,
                               float                 dz,
                               int                   nz,
                               std::string           name)
{
  float z0 = dz/2.0f;
  std::string base_name = IO::PrefixBackground() + IO::PrefixTrend() + name + IO::SuffixAsciiIrapClassic();
  std::string file_name = IO::makeFullFileName(IO::PathToBackground(), base_name);
  std::ofstream file;
  NRLib::OpenWrite(file, file_name);
  for (int i = 0; i < nz; i++) {
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
Background::CalculateDeviationFromVerticalTrend(std::vector<std::vector<double> > & well_trend,
                                                const std::vector<double>         & global_trend,
                                                std::vector<double>               & avg_dev,
                                                const int                           nz)
{
  int n_wells = static_cast<int>(well_trend.size());

  for (int w = 0; w < n_wells; w++) {
    if (well_trend[w].size() > 0) {
      std::vector<double> & w_trend = well_trend[w];
      double sum_dev = 0.0f;
      int count = 0;
      for (int k = 0 ; k < nz ; k++) {
        if (w_trend[k] != RMISSING) {
          double diff = exp(w_trend[k]) - exp(global_trend[k]);
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
Background::WriteDeviationsFromVerticalTrend(const std::vector<double>                        & avg_dev_vp,
                                             const std::vector<double>                        & avg_dev_vs,
                                             const std::vector<double>                        & avg_dev_rho,
                                             const std::vector<double>                        & trend_vp,
                                             const std::vector<double>                        & trend_vs,
                                             const std::vector<double>                        & trend_rho,
                                             const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                             const int                                          n_wells,
                                             const int                                          nz)
{
  double global_mean_vp  = 0.0;
  double global_mean_vs  = 0.0;
  double global_mean_rho = 0.0;

  for (int k = 0; k < nz; k++) {
    global_mean_vp  += exp(trend_vp[k]);
    global_mean_vs  += exp(trend_vs[k]);
    global_mean_rho += exp(trend_rho[k]);
  }
  global_mean_vp  /= nz;
  global_mean_vs  /= nz;
  global_mean_rho /= nz;

  //
  // Find the relative average deviations (mean of Vp,Vs and Rho deviations).
  //
  double * rel_avg_dev = new double[n_wells];
  for (int i = 0; i < n_wells; i++) {
    double rel_dev_vp  = avg_dev_vp[i]/global_mean_vp;
    double rel_dev_vs  = avg_dev_vs[i]/global_mean_vs;
    double rel_dev_rho = avg_dev_rho[i]/global_mean_rho;
    rel_avg_dev[i] = (rel_dev_vp + rel_dev_vs + rel_dev_rho)/3;
  }
  //
  // Sort deviations to find worst well.
  //

  std::vector<int> index(n_wells);
  std::vector<std::string> well_names(n_wells);

  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);

    index[w]      = w;
    well_names[w] = iter->second->GetWellName();
    w++;
  }

  for (int i = 0; i < n_wells; i++) {
    for (int j = i; j < n_wells; j++) {
      if (rel_avg_dev[index[j]] > rel_avg_dev[index[i]]) {
        int tmp = index[i];
        index[i] = index[j];
        index[j] = tmp;
      }
    }
  }
  //
  // Print results
  //
  if (n_wells > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nSummary of average deviation from vertical trend (well with largest misfit listed first):\n\n");
    LogKit::LogFormatted(LogKit::Low,"Well                        Vp       Vs      Rho\n");
    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------\n");
  }
  for (int i = 0; i < n_wells; i++) {

    int ii = index[i];
    if (avg_dev_vp[ii] != RMISSING) {
      LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f    %5.1f    %5.3f\n", well_names[ii].c_str(),
                           avg_dev_vp[ii], avg_dev_vs[ii], avg_dev_rho[ii]);
    }
  }

  if (n_wells == 1) {
    LogKit::LogFormatted(LogKit::High,"\nNOTE: A deviation may be observed even with one well since the global trend is");
    LogKit::LogFormatted(LogKit::High,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
  }
  delete [] rel_avg_dev;
}

//-------------------------------------------------------------------------------
void
Background::FillInVerticalTrend(FFTGrid                   * grid,
                                const std::vector<double> & trend)
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
        grid->setNextReal(static_cast<float>(trend[k]));

  grid->endAccess();
}

//-------------------------------------------------------------------------------
//void
//Background::findMeanVsVp(FFTGrid * Vp,
//                         FFTGrid * Vs)
//{
//  Vp->setAccessMode(FFTGrid::READ);
//  Vs->setAccessMode(FFTGrid::READ);
//  double mean = 0;
//  int nxp = 2*(Vp->getNxp()/2+1);
//  int nyp = Vp->getNyp();
//  int nzp = Vp->getNzp();
//  int nx  = Vp->getNx();
//  int ny  = Vp->getNy();
//  int nz  = Vp->getNz();
//  for (int k=0;k<nzp;k++)
//    for (int j=0;j<nyp;j++)
//      for (int i=0;i<nxp;i++) {
//        float v1 = Vp->getNextReal();
//        float v2 = Vs->getNextReal();
//        if (i < nx && j < ny && k < nz)
//          mean += exp(v2-v1);
//      }
//  mean = mean/double(nx*ny*nz);
//
//  Vp->endAccess();
//  Vs->endAccess();
//
//  vsvp_ = mean;
//}


//-------------------------------------------------------------------------------
//void
//Background::setClassicVsVp()
//{
//  int nxp = back_model_[0]->getNxp();
//  int nyp = back_model_[0]->getNyp();
//  int nzp = back_model_[0]->getNzp();
//  float vp = back_model_[0]->getFirstComplexValue().re;
//  vp = float(exp(vp/sqrt(float(nxp*nyp*nzp))));
//  float vs = back_model_[1]->getFirstComplexValue().re;
//  vs = float(exp(vs/sqrt(float(nxp*nyp*nzp))));
//  vsvp_ = vs/vp;
//}

//-------------------------------------------------------------------------------
void
Background::ResampleBackgroundModel(NRLib::Grid<float>  * & bg_vp,
                                    NRLib::Grid<float>  * & bg_vs,
                                    NRLib::Grid<float>  * & bg_rho,
                                    const Simbox        *   bg_simbox,
                                    const Simbox        *   simbox)
{
  NRLib::Grid<float> * res_bg_vp  = new NRLib::Grid<float>();
  NRLib::Grid<float> * res_bg_vs  = new NRLib::Grid<float>();
  NRLib::Grid<float> * res_bg_rho = new NRLib::Grid<float>();

  LogKit::LogFormatted(LogKit::Low,"\nResampling background model...\n");
  ResampleParameter(res_bg_vp,  bg_vp,  simbox, bg_simbox);
  ResampleParameter(res_bg_vs,  bg_vs,  simbox, bg_simbox);
  ResampleParameter(res_bg_rho, bg_rho, simbox, bg_simbox);

  delete bg_vp;
  delete bg_vs;
  delete bg_rho;

  bg_vp  = res_bg_vp;
  bg_vs  = res_bg_vs;
  bg_rho = res_bg_rho;
}

void
Background::ResampleParameter(NRLib::Grid<float> *& p_new, // Resample to
                              NRLib::Grid<float> *  p_old, // Resample from
                              const Simbox       *  simbox_new,
                              const Simbox       *  simbox_old)
{
  int nx = simbox_new->getnx();
  int ny = simbox_new->getny();
  int nz = simbox_new->getnz();
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
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double dz_new = simbox_new->getdz(i,j);
      double dz_old = simbox_old->getdz(i,j);
      double z0_new = simbox_new->getTop(i,j);
      double z0_old = simbox_old->getTop(i,j);
        a[ij] = dz_new/dz_old;
      b[ij] = (z0_new - z0_old)/dz_old;
      ij++;
    }
  }

  //
  // Resample parameter
  //
  p_new->Resize(nx, ny, nz);

  double * layer = new double[nx*ny];

  for (int k = 0; k < nz; k++) {
    //
    // Map a layer
    //
    int ij=0;
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int k_old = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]);
        layer[ij] = p_old->GetValue(i, j, k_old);
        ij++;
      }
    }
    //
    // Smooth the layer (equal weighting of all neighbouring cells)
    //
    double value = 0.0;
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int n = 1;
        double sum = layer[j*nx + i];
        if (i > 1) {
          sum += layer[j*nx + i - 1];
          n++;
        }
        if (j > 1) {
          sum += layer[(j - 1)*nx + i];
          n++;
        }
        if (i > 1 && j > 1) {
          sum += layer[(j - 1)*nx + i - 1];
          n++;
        }
        if (i < nx-1) {
          sum += layer[j*nx + i + 1];
          n++;
        }
        if (j < ny-1) {
          sum += layer[(j + 1)*nx + i];
          n++;
        }
        if (i < nx-1 && j < ny-1) {
          sum += layer[(j + 1)*nx + i + 1];
          n++;
        }
        value = sum/static_cast<double>(n);

        p_new->SetValue(i, j, k, static_cast<float>(value));
      }
    }
  }

  delete [] layer;
  delete [] a;
  delete [] b;
}

//-------------------------------------------------------------------------------
//void
//Background::padAndSetBackgroundModel(FFTGrid * bg_vp,
//                                     FFTGrid * bg_vs,
//                                     FFTGrid * bg_rho)
//{
//  //LogKit::LogFormatted(LogKit::Low,"\nPadding background model...\n");
//  createPaddedParameter(back_model_[0], bg_vp);
//  createPaddedParameter(back_model_[1], bg_vs);
//  createPaddedParameter(back_model_[2], bg_rho);
//}

//-------------------------------------------------------------------------------
//void
//Background::createPaddedParameter(FFTGrid *& pNew,     // Padded
//                                  FFTGrid  * pOld)     // Non-padded
//{
//  //
//  // Fill padding using linear interpolation between edges.
//  //
//  // When we fill the z-padding, we assume that x- and y-padding is
//  // already filled. The loop structure ensures this. Likewise, it
//  // is assumed that the x-padding is filled when we fill the
//  // y-padding.
//  //
//  // The linear algortihm is not "perfect", but should be more
//  // than good enough for padding the smooth background model.
//  //
//  int nx   = pNew->getNx();
//  int ny   = pNew->getNy();
//  int nz   = pNew->getNz();
//  int nxp  = pNew->getNxp();
//  int nyp  = pNew->getNyp();
//  int nzp  = pNew->getNzp();
//  int rnxp = pNew->getRNxp();
//
//  pNew->createRealGrid();
//  pNew->setType(FFTGrid::PARAMETER);
//
//  pNew->setAccessMode(FFTGrid::RANDOMACCESS);
//  pOld->setAccessMode(FFTGrid::RANDOMACCESS);
//
//  float sum_c = 1.0f/static_cast<float>(nzp - nz + 1);
//  float sum_b = 1.0f/static_cast<float>(nyp - ny + 1);
//  float sum_a = 1.0f/static_cast<float>(nxp - nx + 1);
//
//  for (int k = 0 ; k < nzp ; k++) {
//    for (int j = 0 ; j < nyp ; j++) {
//      for (int i = 0 ; i < rnxp ; i++) { // Must fill entire grid to avoid UMR.
//
//        float value = RMISSING;
//        if (i < nx && j < ny && k < nz) { // Not in padding
//          value = pOld->getRealValue(i, j, k);
//        }
//        else {
//          if (i >= nxp)       //In dummy area for real grid, but fill to avoid UMR.
//            value = 0;
//          else if (k >= nz) { // In z-padding (x- and y- padding is filled in pNew)
//            float c1 = pNew->getRealValue(i, j, 0     , true);
//            float c2 = pNew->getRealValue(i, j, nz - 1, true);
//            float w1 = sum_c*static_cast<float>(k - nz + 1);
//            float w2 = sum_c*static_cast<float>(nzp - k);
//            value = c1*w1 + c2*w2;
//          }
//          else if (j >= ny) { // In y-padding (x-padding is filled in pNew)
//            float b1 = pNew->getRealValue(i, 0     , k, true);
//            float b2 = pNew->getRealValue(i, ny - 1, k, true);
//            float w1 = sum_b*static_cast<float>(j - ny + 1);
//            float w2 = sum_b*static_cast<float>(nyp - j);
//            value = b1*w1 + b2*w2;
//          }
//          else if (i >= nx) { // In x-padding
//            float a1 = pNew->getRealValue(     0, j, k, true);
//            float a2 = pNew->getRealValue(nx - 1, j, k, true);
//            float w1 = sum_a*static_cast<float>(i - nx + 1);
//            float w2 = sum_a*static_cast<float>(nxp - i);
//            value = a1*w1 + a2*w2;
//          }
//        }
//        pNew->setRealValue(i,j,k,value,true);
//      }
//    }
//  }
//  pNew->endAccess();
//  pOld->endAccess();
//}

//-------------------------------------------------------------------------------
//void
//Background::writeBackgrounds(const Simbox            * simbox,
//                             GridMapping             * depthMapping,
//                             const GridMapping       * timeMapping,
//                             const bool                isFile,
//                             const TraceHeaderFormat & thf) const
//{
//  if (depthMapping != NULL && depthMapping->getSimbox() == NULL) {
//    const Simbox * timeSimbox = simbox;
//    if (timeMapping != NULL)
//      timeSimbox = timeMapping->getSimbox();
//    back_model_[0]->setAccessMode(FFTGrid::RANDOMACCESS);
//    depthMapping->setMappingFromVelocity(back_model_[0], timeSimbox);
//    back_model_[0]->endAccess();
//  }
//
//  std::string fileName1 = IO::PrefixBackground() + "Vp" ;
//  std::string fileName2 = IO::PrefixBackground() + "Vs" ;
//  std::string fileName3 = IO::PrefixBackground() + "Rho";
//
//  FFTGrid * expAlpha = CopyFFTGrid(back_model_[0], true, isFile);
//  expAlpha->writeFile(fileName1, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
//  delete expAlpha;
//
//  FFTGrid * expBeta = CopyFFTGrid(back_model_[1], true, isFile);
//  expBeta->writeFile(fileName2, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
//  delete expBeta;
//
//  FFTGrid * expRho = CopyFFTGrid(back_model_[2], true, isFile);
//  expRho->writeFile(fileName3, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
//  delete expRho;
//
//  //
//  // For debugging: write cubes not in ASCII, with padding, and with flat top.
//  //
//  //back_model_[0]->writeStormFile(fileName1, simbox, true, false, true, true);
//  //back_model_[1]->writeStormFile(fileName2, simbox, true, false, true, true);
//  //back_model_[2]->writeStormFile(fileName3, simbox, true, false, true, true);
//}

FFTGrid *
Background::CopyFFTGrid(FFTGrid   * orig_grid,
                        const bool  exp_trans,
                        const bool  file_grid) const
{
  FFTGrid * new_grid;

  if (file_grid)
    new_grid = new FFTFileGrid(static_cast<FFTFileGrid *>(orig_grid), exp_trans);
  else
    new_grid = new FFTGrid(orig_grid, exp_trans);

  return (new_grid);
}

//---------------------------------------------------------------------------

//void
//Background::ErodeSurface(Surface       *& surface,
//                         const Surface *  priority_surface,
//                         const Simbox  *  simbox,
//                         const bool    &  compare_upward) const
//{
//  int nx    = simbox->getnx();
//  int ny    = simbox->getny();
//  double x0 = simbox->GetXMin();
//  double y0 = simbox->GetYMin();
//  double lx = simbox->GetLX();
//  double ly = simbox->GetLY();
//
//  NRLib::Grid2D<double> eroded_surface(nx,ny,0);
//  double x;
//  double y;
//  double z;
//  double z_priority;
//
//  double missing = surface->GetMissingValue();
//  for (int i = 0; i < nx; i++) {
//    for (int j = 0; j < ny; j++) {
//      simbox->getXYCoord(i,j,x,y);
//
//      z_priority = priority_surface->GetZ(x,y);
//      z          = surface->GetZ(x,y);
//
//      if (compare_upward) {
//        if (z < z_priority && z != missing)
//          eroded_surface(i,j) = z_priority;
//        else
//          eroded_surface(i,j) = z;
//      }
//
//      else {
//        if (z > z_priority && z_priority != missing)
//          eroded_surface(i,j) = z_priority;
//        else
//          eroded_surface(i,j) = z;
//      }
//    }
//  }
//  delete surface;
//
//  surface = new Surface(x0, y0, lx, ly, eroded_surface);
//}
