/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <typeinfo>
#include <algorithm>

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/analyzelog.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/io.h"
#include "src/waveletfilter.h"
#include "src/tasklist.h"
#include "src/timeline.h"
#include "src/state4d.h"
#include "src/modelavodynamic.h"
#include "src/cravatrend.h"
#include "src/seismicparametersholder.h"
#include "src/parameteroutput.h"

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/surface/regularsurfacerotated.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"

#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsrock.h"


//ModelGeneral::ModelGeneral(ModelSettings           *& model_settings,
//                           const InputFiles         * input_files,
//                           SeismicParametersHolder  & seismicParameters,
//                           Simbox                  *& timeBGSimbox)
//                           :do_4D_inversion_(model_settings->getDo4DInversion()),
//                          do_4D_rock_physics_vnversion_(model_settings->getDo4DRockPhysicsInversion())
//{
//  simbox_             = new Simbox();
//  timeSimboxConstThick_   = NULL;
//
//  correlationDirection_   = NULL;
//  random_gen_              = NULL;
//  failed_                 = false;
//  gradX_                  = 0.0;
//  gradY_                  = 0.0;
//
//  timeDepthMapping_       = NULL;
//  timeCutMapping_         = NULL;
//  velocityFromInversion_  = false;
//
//  bool failedSimbox       = false;
//  bool failedDepthConv    = false;
//  bool failedRockPhysics  = false;
//
//  bool failedLoadingModel = false;
//
//  bool failedWells        = false;
//
//  bool failedBackground   = false;
//
//  Simbox * timeCutSimbox  = NULL;
//  time_line_               = NULL;
//
//  forwardModeling_        = model_settings->getForwardModeling();
//  numberOfWells_          = model_settings->getNumberOfWells();
//
//
//  {
//    int debugLevel = model_settings->getLogLevel();
//    if (model_settings->getDebugLevel() == 1)
//      debugLevel = LogKit::L_DebugLow;
//    else if (model_settings->getDebugLevel() == 2)
//      debugLevel = LogKit::L_DebugHigh;
//
//    LogKit::SetScreenLog(debugLevel);
//
//    std::string logFileName = IO::makeFullFileName("",IO::FileLog()+IO::SuffixTextFiles());
//    LogKit::SetFileLog(logFileName,model_settings->getLogLevel());
//
//    if (model_settings->getDebugFlag() > 0)
//    {
//      std::string fName = IO::makeFullFileName("",IO::FileDebug()+IO::SuffixTextFiles());
//      LogKit::SetFileLog(fName, debugLevel);
//    }
//
//    if (model_settings->getErrorFileFlag() == true)
//    {
//      std::string fName = IO::makeFullFileName("",IO::FileError()+IO::SuffixTextFiles());
//      LogKit::SetFileLog(fName, LogKit::Error);
//    }
//    LogKit::EndBuffering();
//
//    if (input_files->getSeedFile() == "")
//      random_gen_ = new RandomGen(model_settings->getSeed());
//    else
//      random_gen_ = new RandomGen(input_files->getSeedFile().c_str());
//
//    if (model_settings->getNumberOfSimulations() == 0)
//      model_settings->setWritePrediction(true); //write predicted grids.
//
//    printSettings(model_settings, input_files);
//
//    //Set output for all FFTGrids.
//    FFTGrid::setOutputFlags(model_settings->getOutputGridFormat(),
//                            model_settings->getOutputGridDomain());
//
//    std::string errText("");
//
//    LogKit::WriteHeader("Defining modelling grid");
//    makeTimeSimboxes(simbox_, timeCutSimbox, timeBGSimbox, timeSimboxConstThick_,  //Handles correlation direction too.
//                     correlationDirection_, model_settings, input_files,
//                     errText, failedSimbox);
//
//    if (!failedSimbox)
//    {
//      //
//      // FORWARD MODELLING
//      //
//      if (model_settings->getForwardModeling() == true)
//      {
//  //      checkAvailableMemory(simbox_, model_settings, input_files);
//      }
//      else {
//        //
//        // INVERSION/ESTIMATION
//        //
//        if (timeCutSimbox!=NULL)  {
//          timeCutMapping_ = new GridMapping();
//          timeCutMapping_->makeTimeTimeMapping(timeCutSimbox);
//        }
//
//        //checkAvailableMemory(simbox_, model_settings, input_files);
//
//        bool estimationMode = model_settings->getEstimationMode();
//
//        if (estimationMode == false && model_settings->getDoDepthConversion() == true)
//        {
//          processDepthConversion(timeCutSimbox, simbox_, model_settings,
//                                 input_files, errText, failedDepthConv);
//        }
//
//        processWells(wells_, simbox_, model_settings, input_files, errText, failedWells);
//
//        if (failedDepthConv == false)
//          processRockPhysics(simbox_, timeCutSimbox, model_settings, failedRockPhysics, errText, input_files);
//
//        //Set up timeline.
//        time_line_ = new TimeLine();
//        //Activate below when gravity data are ready.
//        //Do gravity first.
//        //for (int i=0;i<model_settings->getNumberOfGravityData();i++) {
//        //  int time = computeTime(model_settings->getGravityYear[i],
//        //                         model_settings->getGravityMonth[i],
//        //                         model_settings->getGravityDay[i]);
//        //  time_line_->AddEvent(time, TimeLine::GRAVITY, i);
//
//        bool firstGravimetricEvent = true;
//        for (int i=0;i<model_settings->getNumberOfVintages();i++) {
//          //Vintages may have both travel time and AVO
//          int time = computeTime(model_settings->getVintageYear(i),
//                                 model_settings->getVintageMonth(i),
//                                 model_settings->getVintageDay(i));
//           // Do gravity first
//           if (model_settings->getGravityTimeLapse(i)){
//             if (firstGravimetricEvent){
//               // Do not save first gravity event in timeline
//               firstGravimetricEvent = false;
//             }
//             else{
//               time_line_->AddEvent(time, TimeLine::GRAVITY, i);
//             }
//          }
//          //Activate below when travel time is ready.
//          //Travel time ebefore AVO for same vintage.
//          //if (travel time for this vintage)
//          //time_line_->AddEvent(time, TimeLine::TRAVEL_TIME, i);
//          if (model_settings->getNumberOfAngles(i) > 0) //Check for AVO data, could be pure travel time.
//            time_line_->AddEvent(time, TimeLine::AVO, i);
//        }
//
//        if (model_settings->getDo4DInversion() && failedRockPhysics == false){
//
//          setFaciesNamesFromRockPhysics();
//
//          NRLib::Vector initialMean(6);
//          NRLib::Matrix initialCov(6,6);
//          process4DBackground(model_settings, input_files, seismicParameters, errText, failedBackground,initialMean,initialCov);
//
//          time_evolution_ = TimeEvolution(10000, *time_line_, rock_distributions_.begin()->second); //NBNB OK 10000->1000 for speed during testing
//          time_evolution_.SetInitialMean(initialMean);
//          time_evolution_.SetInitialCov(initialCov);
//        }
//      }
//    }
//    failedLoadingModel = failedSimbox  || failedDepthConv || failedWells || failedBackground || failedRockPhysics;
//
//    if (failedLoadingModel) {
//      LogKit::WriteHeader("Error(s) while loading data");
//      LogKit::LogFormatted(LogKit::Error,"\n"+errText);
//      LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
//    }
//  }
//
//  failed_ = failedLoadingModel;
//  failed_details_.push_back(failedSimbox);
//  failed_details_.push_back(failedDepthConv);
//  failed_details_.push_back(failedWells);
//  failed_details_.push_back(failedBackground);
//  failed_details_.push_back(failedRockPhysics);
//
//  if (timeCutSimbox != NULL)
//    delete timeCutSimbox;
//}

ModelGeneral::ModelGeneral(ModelSettings           *& model_settings, //Multiple intervals
                           const InputFiles         * input_files,
                           SeismicParametersHolder  & seismic_parameters,
                           CommonData               * common_data,
                           int                        i_interval)
                           :do_4D_inversion_(model_settings->getDo4DInversion()),
                            do_4D_rock_physics_vnversion_(model_settings->getDo4DRockPhysicsInversion())
{
  simbox_ = common_data->GetMultipleIntervalGrid()->GetIntervalSimbox(i_interval),
  random_gen_    = NULL;
  time_line_     = NULL;
  multi_interval_ = false;

  {
    int debug_level = model_settings->getLogLevel();
    if (model_settings->getDebugLevel() == 1)
      debug_level = LogKit::L_DebugLow;
    else if (model_settings->getDebugLevel() == 2)
      debug_level = LogKit::L_DebugHigh;

    LogKit::SetScreenLog(debug_level);

    std::string log_file_name = IO::makeFullFileName("",IO::FileLog()+IO::SuffixTextFiles());
    LogKit::SetFileLog(log_file_name, model_settings->getLogLevel());

    if (model_settings->getDebugFlag() > 0) {
      std::string f_name = IO::makeFullFileName("",IO::FileDebug()+IO::SuffixTextFiles());
      LogKit::SetFileLog(f_name, debug_level);
    }

    if (model_settings->getErrorFileFlag() == true) {
      std::string f_name = IO::makeFullFileName("",IO::FileError()+IO::SuffixTextFiles());
      LogKit::SetFileLog(f_name, LogKit::Error);
    }
    LogKit::EndBuffering();

    if (input_files->getSeedFile() == "")
      random_gen_ = new RandomGen(model_settings->getSeed());
    else
      random_gen_ = new RandomGen(input_files->getSeedFile().c_str());

    if (model_settings->getNumberOfSimulations() == 0)
      model_settings->setWritePrediction(true); //write predicted grids.

    printSettings(model_settings, input_files);

    //Set output for all FFTGrids.
    FFTGrid::setOutputFlags(model_settings->getOutputGridFormat(),
                            model_settings->getOutputGridDomain());

    std::string errText("");

    LogKit::WriteHeader("Defining modelling grid");

    MultiIntervalGrid * multiple_interval_grid = common_data->GetMultipleIntervalGrid();
    //simbox_                                    = multiple_interval_grid->GetIntervalSimboxE(i_interval);

    //
    // FORWARD MODELLING
    //
    if (model_settings->getForwardModeling() == true) {
//      checkAvailableMemory(simbox_, model_settings, input_files);
    }
    else {
      //
      // INVERSION/ESTIMATION
      //
      if (model_settings->getIntervalNames().size() > 0)
        multi_interval_ = true;

      //checkAvailableMemory(simbox_, model_settings, input_files);

      //Facies-names
      facies_names_  = common_data->GetFaciesNames();
      facies_labels_ = common_data->GetFaciesLabels();

      //Priorfacies
      if (common_data->GetPriorFacies().size() > 0) {
        prior_facies_ = common_data->GetPriorFaciesInterval(i_interval);

        prior_facies_prob_cubes_.resize(3);
        prior_facies_prob_cubes_[0] = new FFTGrid(multiple_interval_grid->GetPriorFaciesProbCubesInterval(i_interval)[0], simbox_->GetNXpad(), simbox_->GetNYpad(), simbox_->GetNZpad());
        prior_facies_prob_cubes_[1] = new FFTGrid(multiple_interval_grid->GetPriorFaciesProbCubesInterval(i_interval)[1], simbox_->GetNXpad(), simbox_->GetNYpad(), simbox_->GetNZpad());
        prior_facies_prob_cubes_[2] = new FFTGrid(multiple_interval_grid->GetPriorFaciesProbCubesInterval(i_interval)[2], simbox_->GetNXpad(), simbox_->GetNYpad(), simbox_->GetNZpad());
      }

      bool estimation_mode = model_settings->getEstimationMode();

      //TimeDepthMapping if intervals isn't used.
      if (common_data->GetMultipleIntervalGrid()->GetNIntervals() == 1) {
        time_depth_mapping_      = common_data->GetTimeDepthMapping();
        velocity_from_inversion_ = common_data->GetVelocityFromInversion(); //Needed?
      }

      //Replace wells with blocked_logs
      blocked_logs_ = common_data->GetBlockedLogs();

      if(multiple_interval_grid->GetTrendCubes().size() > 0)
        trend_cubes_ = multiple_interval_grid->GetTrendCube(i_interval);

      rock_distributions_  = common_data->GetDistributionsRock();
      reservoir_variables_ = common_data->GetReservoirVariables();

      //Set up timeline
      time_line_ = common_data->GetTimeLine();

      if (model_settings->getDo4DInversion()) {

        //setFaciesNamesFromRockPhysics();
        NRLib::Vector initial_mean(6);
        NRLib::Matrix initial_cov(6,6);

        SetupState4D(model_settings, seismic_parameters, initial_mean, initial_cov);

        time_evolution_ = TimeEvolution(10000, *time_line_, rock_distributions_.begin()->second); //NBNB OK 10000->1000 for speed during testing
        time_evolution_.SetInitialMean(initial_mean);
        time_evolution_.SetInitialCov(initial_cov);
      }
    }
  }

}


ModelGeneral::~ModelGeneral(void)
{
  if (time_depth_mapping_!=NULL)
    delete time_depth_mapping_;

  //if (timeCutMapping_!=NULL)
  //  delete timeCutMapping_;

  //if (correlationDirection_ !=NULL)
  //  delete correlationDirection_;

  for (std::map<std::string, std::vector<DistributionsRock *> >::iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
    std::vector<DistributionsRock *> rock = it->second;
    for (size_t i=0; i<rock.size(); i++)
      delete rock[i];
  }

  for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator it = reservoir_variables_.begin(); it != reservoir_variables_.end(); it++) {
    std::vector<DistributionWithTrend *> variable = it->second;
    for (size_t i=0; i<variable.size(); i++)
      delete variable[i];
  }

  if (time_line_ != NULL)
    delete time_line_;

  delete random_gen_;
  delete simbox_;
  //delete timeSimboxConstThick_;

  // if (!forwardModeling_)
  //{
  //  for (int i=0 ; i<numberOfWells_ ; i++)
  //    if (wells_[i] != NULL)
  //      delete wells_[i];
  //}

  // delete priorCorrXY_;

}

//void
//ModelGeneral::readSegyFile(const std::string       & fileName,
//                           FFTGrid                *& target,
//                           const Simbox            * timeSimbox,
//                           const Simbox            * timeCutSimbox,
//                           const ModelSettings     * model_settings,
//                           const SegyGeometry     *& geometry,
//                           int                       gridType,
//                           const std::string       & parName,
//                           float                     offset,
//                           const TraceHeaderFormat * format,
//                           std::string             & errText,
//                           bool                      nopadding)
//{
//  SegY * segy = NULL;
//  bool failed = false;
//  target = NULL;
//
//  try
//  {
//    //
//    // Currently we have only one optional TraceHeaderFormat, but this can
//    // be augmented to a list with several formats ...
//    //
//    if (format == NULL) { //Unknown format
//      std::vector<TraceHeaderFormat*> traceHeaderFormats(0);
//      if (model_settings->getTraceHeaderFormat() != NULL)
//      {
//        traceHeaderFormats.push_back(model_settings->getTraceHeaderFormat());
//      }
//      segy = new SegY(fileName,
//                      offset,
//                      traceHeaderFormats,
//                      true); // Add standard formats to format search
//    }
//    else //Known format, read directly.
//      segy = new SegY(fileName, offset, *format);
//
//    float guard_zone = model_settings->getGuardZone();
//
//    std::string errTxt = "";
//    checkThatDataCoverGrid(segy,
//                           offset,
//                           timeCutSimbox,
//                           guard_zone,
//                           errTxt);
//
//    if (errTxt == "") {
//      bool  onlyVolume      = true;
//      // This is *not* the same as FFT-grid padding. If the padding
//      // size is changed from 2*guard_zone, the smoothing done in
//      // FFTGrid::smoothTraceInGuardZone() will become incorrect.
//      float padding         = 2*guard_zone;
//      bool  relativePadding = false;
//
//      segy->ReadAllTraces(timeCutSimbox,
//                          padding,
//                          onlyVolume,
//                          relativePadding);
//      segy->CreateRegularGrid();
//    }
//    else {
//      errText += errTxt;
//      failed = true;
//    }
//  }
//  catch (NRLib::Exception & e)
//  {
//    errText += e.what();
//    failed = true;
//  }
//
//  if (!failed)
//  {
//    int missingTracesSimbox  = 0;
//    int missingTracesPadding = 0;
//    int deadTracesSimbox     = 0;
//
//    const SegyGeometry * geo;
//    geo = segy->GetGeometry();
//    geo->WriteGeometry();
//    if (gridType == FFTGrid::DATA)
//      geometry = new SegyGeometry(geo);
//
//    int xpad, ypad, zpad;
//    if (nopadding)
//    {
//      xpad = timeSimbox->getnx();
//      ypad = timeSimbox->getny();
//      zpad = timeSimbox->getnz();
//    }
//    else
//    {
//      xpad = model_settings->getNXpad();
//      ypad = model_settings->getNYpad();
//      zpad = model_settings->getNZpad();
//    }
//    target = createFFTGrid(timeSimbox->getnx(),
//                           timeSimbox->getny(),
//                           timeSimbox->getnz(),
//                           xpad,
//                           ypad,
//                           zpad,
//                           model_settings->getFileGrid());
//    target->setType(gridType);
//
//    if (gridType == FFTGrid::DATA) {
//      target->fillInSeismicDataFromSegY(segy,
//                                        timeSimbox,
//                                        model_settings->getSmoothLength(),
//                                        missingTracesSimbox,
//                                        missingTracesPadding,
//                                        deadTracesSimbox,
//                                        errText);
//    }
//    else {
//      missingTracesSimbox = target->fillInFromSegY(segy,
//                                                   timeSimbox,
//                                                   parName,
//                                                   nopadding);
//    }
//
//    if (missingTracesSimbox > 0) {
//      if (missingTracesSimbox == timeSimbox->getnx()*timeSimbox->getny()) {
//        errText += "Error: Data in file "+fileName+" was completely outside the inversion area.\n";
//        failed = true;
//      }
//      else {
//        if (gridType == FFTGrid::PARAMETER) {
//          errText += "Grid in file "+fileName+" does not cover the inversion area.\n";
//        }
//        else {
//          LogKit::LogMessage(LogKit::Warning, "WARNING: "+NRLib::ToString(missingTracesSimbox)
//                             +" grid columns are outside the area defined by the seismic data.\n");
//          std::string text;
//          text += "Check seismic volumes and inversion area: A part of the inversion area is outside\n";
//          text += "   the seismic data specified in file \'"+fileName+"\'.";
//          TaskList::addTask(text);
//        }
//      }
//    }
//    if (missingTracesPadding > 0) {
//      int nx     = timeSimbox->getnx();
//      int ny     = timeSimbox->getny();
//      int nxpad  = xpad - nx;
//      int nypad  = ypad - ny;
//      int nxypad = nxpad*ny + nx*nypad - nxpad*nypad;
//      LogKit::LogMessage(LogKit::High, "Number of grid columns in padding that are outside area defined by seismic data : "
//                         +NRLib::ToString(missingTracesPadding)+" of "+NRLib::ToString(nxypad)+"\n");
//    }
//    if (deadTracesSimbox > 0) {
//      LogKit::LogMessage(LogKit::High, "Number of grid columns with no seismic data (nearest trace is dead) : "
//                         +NRLib::ToString(deadTracesSimbox)+" of "+NRLib::ToString(timeSimbox->getnx()*timeSimbox->getny())+"\n");
//    }
//  }
//  if (segy != NULL)
//    delete segy;
//}


//void
//ModelGeneral::checkThatDataCoverGrid(const SegY   * segy,
//                                     float         offset,
//                                     const Simbox * timeCutSimbox,
//                                     float         guard_zone,
//                                     std::string & errText)
//{
//  // Seismic data coverage (translate to CRAVA grid by adding half a grid cell)
//  float dz = segy->GetDz();
//  float z0 = offset + 0.5f*dz;
//  float zn = z0 + (segy->GetNz() - 1)*dz;
//
//  // Top and base of interval of interest
//  float top_grid = static_cast<float>(timeCutSimbox->getTopZMin());
//  float bot_grid = static_cast<float>(timeCutSimbox->getBotZMax());
//
//  // Find guard zone
//  float top_guard = top_grid - guard_zone;
//  float bot_guard = bot_grid + guard_zone;
//
//  if (top_guard < z0) {
//    float z0_new = z0 - ceil((z0 - top_guard)/dz)*dz;
//    errText += "\nThere is not enough seismic data above the interval of interest. The seismic data\n";
//    errText += "must start at "+NRLib::ToString(z0_new)+"ms (in CRAVA grid) to allow for a ";
//    errText += NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n";
//    errText += "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n";
//    errText += "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n";
//    errText += "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n";
//    errText += "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n";
//    errText += "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n";
//    errText += "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n";
//  }
//  if (bot_guard > zn) {
//    float zn_new = zn + ceil((bot_guard - zn)/dz)*dz;
//    errText += "\nThere is not enough seismic data below the interval of interest. The seismic data\n";
//    errText += "must end at "+NRLib::ToString(zn_new)+"ms (in CRAVA grid) to allow for a ";
//    errText += NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n";
//    errText += "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n";
//    errText += "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n";
//    errText += "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n";
//    errText += "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n";
//    errText += "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n";
//    errText += "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n";
//  }
//}


//void
//ModelGeneral::readStormFile(const std::string   & fName,
//                            FFTGrid            *& target,
//                            const int             gridType,
//                            const std::string   & parName,
//                            const Simbox        * timeSimbox,
//                            const ModelSettings * model_settings,
//                            std::string         & errText,
//                            bool                  scale,
//                            bool                  nopadding)
//{
//  StormContGrid * stormgrid = NULL;
//  bool failed = false;
//
//  try
//  {
//    stormgrid = new StormContGrid(0,0,0);
//    stormgrid->ReadFromFile(fName);
//  }
//  catch (NRLib::Exception & e)
//  {
//    errText += e.what();
//    failed = true;
//  }
//  int xpad, ypad, zpad;
//  if (nopadding==false)
//  {
//    xpad = model_settings->getNXpad();
//    ypad = model_settings->getNYpad();
//    zpad = model_settings->getNZpad();
//  }
//  else
//  {
//    xpad = timeSimbox->getnx();
//    ypad = timeSimbox->getny();
//    zpad = timeSimbox->getnz();
//  }
//
//  int outsideTraces = 0;
//  if (failed == false)
//  {
//    target = createFFTGrid(timeSimbox->getnx(),
//                           timeSimbox->getny(),
//                           timeSimbox->getnz(),
//                           xpad,
//                           ypad,
//                           zpad,
//                           model_settings->getFileGrid());
//    target->setType(gridType);
//
//    //try {
//    //  outsideTraces = target->fillInFromStorm(timeSimbox,stormgrid, parName, scale, nopadding);
//    //}
//    //catch (NRLib::Exception & e) {
//    //  errText += std::string(e.what());
//    //}
//  }
//
//  if (stormgrid != NULL)
//    delete stormgrid;
//
//  if (outsideTraces > 0) {
//    if (outsideTraces == timeSimbox->getnx()*timeSimbox->getny()) {
//      errText += "Error: Data in file \'"+fName+"\' was completely outside the inversion area.\n";
//      failed = true;
//    }
//    else {
//      if (gridType == FFTGrid::PARAMETER) {
//        errText += "Error: Data read from file \'"+fName+"\' does not cover the inversion area.\n";
//      }
//      else {
//        LogKit::LogMessage(LogKit::Warning, "WARNING: "+NRLib::ToString(outsideTraces)
//                           + " grid columns were outside the seismic data in file \'"+fName+"\'.\n");
//        TaskList::addTask("Check seismic data and inversion area: One or volumes did not have data enough to cover entire grid.\n");
//     }
//    }
//  }
//}

int
ModelGeneral::setPaddingSize(int nx, double px)
{
  int leastint    = static_cast<int>(ceil(nx*(1.0f+px)));
  int closestprod = FFTGrid::findClosestFactorableNumber(leastint);
  return(closestprod);
}

/*
void
ModelGeneral::makeTimeSimboxes(Simbox          *& timeSimbox,
                               Simbox          *& timeCutSimbox,
                               Simbox          *& timeBGSimbox,
                               Simbox          *& timeSimboxConstThick,
                               Surface         *& correlationDirection,
                               ModelSettings   *& modelSettings,
                               const InputFiles * inputFiles,
                               std::string      & errText,
                               bool             & failed)
{
  std::string gridFile("");


//void
//ModelGeneral::logIntervalInformation(const Simbox      * simbox,
//                                     const std::string & header_text1,
//                                     const std::string & header_text2)
//{
//  LogKit::LogFormatted(LogKit::Low,"\n"+header_text1+"\n");
//  double zmin, zmax;
//  simbox->getMinMaxZ(zmin,zmax);
//  LogKit::LogFormatted(LogKit::Low," %13s          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
//                       header_text2.c_str(),
//                       zmin+simbox->getlz()*simbox->getAvgRelThick()*0.5,
//                       zmin,zmax);
//  LogKit::LogFormatted(LogKit::Low,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n",
//                       simbox->getlz()*simbox->getAvgRelThick(),
//                       simbox->getlz()*simbox->getMinRelThick(),
//                       simbox->getlz());
//  LogKit::LogFormatted(LogKit::Low,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n",
//                       simbox->getdz()*simbox->getAvgRelThick(),
//                       simbox->getdz(),
//                       simbox->getdz()*simbox->getMinRelThick());
//}

  if(modelSettings->getForwardModeling())
    gridFile = inputFiles->getBackFile(0);    // Get geometry from earth model (Vp)
  else {
    if (modelSettings->getEstimationMode() == false || estimationModeNeedILXL)
      gridFile = inputFiles->getSeismicFile(0,0); // Get area from first seismic data volume
  }
  SegyGeometry * ILXLGeometry = NULL; //Geometry with correct XL and IL settings.

  //
  // Set area geometry information
  // -----------------------------
  //
  std::string areaType;
  if (areaSpecification == ModelSettings::AREA_FROM_UTM)
  {
    // The geometry is already present in modelSettings (read from model file).
    LogKit::LogFormatted(LogKit::High,"\nArea information has been taken from model file\n");
    areaType = "Model file";
  }
  else if(areaSpecification == ModelSettings::AREA_FROM_SURFACE)
  {
    LogKit::LogFormatted(LogKit::High,"\nFinding area information from surface \'"+inputFiles->getAreaSurfaceFile()+"\'\n");
    areaType = "Surface";
    RotatedSurface surf(inputFiles->getAreaSurfaceFile());
    SegyGeometry geometry(surf);
    modelSettings->setAreaParameters(&geometry);
  }
  else if(areaSpecification == ModelSettings::AREA_FROM_GRID_DATA         ||
          areaSpecification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM ||
          areaSpecification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE)
  {
    LogKit::LogFormatted(LogKit::High,"\nFinding inversion area from grid data in file \'"+gridFile+"\'\n");
    areaType = "Grid data";
    std::string tmpErrText;
    SegyGeometry * geometry;
    getGeometryFromGridOnFile(gridFile,
                              modelSettings->getTraceHeaderFormat(0,0), //Trace header format is the same for all time lapses
                              geometry,
                              tmpErrText);

    modelSettings->setSeismicDataAreaParameters(geometry);
    if(geometry != NULL) {
      geometry->WriteGeometry();

      if (modelSettings->getAreaILXL().size() > 0 || modelSettings->getSnapGridToSeismicData()) {
        SegyGeometry * fullGeometry = geometry;

        std::vector<int> areaILXL;
        bool gotArea = true;

        //
        // Geometry is given as XY, but we snap it to IL and XL.
        //
        if (modelSettings->getSnapGridToSeismicData()) {
          SegyGeometry * templateGeometry = NULL;
          if (areaSpecification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM) {
            templateGeometry = modelSettings->getAreaParameters();
          }
          else if (areaSpecification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE) {
            RotatedSurface surf(inputFiles->getAreaSurfaceFile());
            templateGeometry = new SegyGeometry(surf);
          }
          else {
            errText += "CRAVA has been asked to identify a proper ILXL inversion area based\n";
            errText += "on XY input information, but no UTM coordinates or surface have\n";
            errText += "been specified in model file.\n";
            gotArea = false;
          }
          if (gotArea) {
            areaILXL = fullGeometry->findAreaILXL(templateGeometry);
          }
        }
        else {
          areaILXL = modelSettings->getAreaILXL();
        }

        if (gotArea) {
          try {
            bool interpolated, aligned;
            geometry = fullGeometry->GetILXLSubGeometry(areaILXL, interpolated, aligned);

            std::string text;
            if(interpolated == true) {
              if(aligned == true) {
                text  = "Check IL/XL specification: Specified IL- or XL-step is not an integer multiple\n";
                text += "   of those found in the seismic data. Furthermore, the distance between first\n";
                text += "   and last XL and/or IL does not match the step size.\n";
                TaskList::addTask(text);
              }
              else {
                text  = "Check IL/XL specification: Specified IL- or XL-step is not an integer multiple\n";
                text += "   of those found in the seismic data.\n";
                TaskList::addTask(text);
              }
            }
            else if(aligned == true) {
              text  = "Check IL/XL specification: Either start or end of IL and/or XL interval does not\n";
              text += "   align with IL/XL in seismic data, or end IL and/or XL is not an integer multiple\n";
              text += "   of steps away from the start.\n";
              TaskList::addTask(text);
            }
          }
          catch (NRLib::Exception & e) {
            errText += "Error: "+std::string(e.what());
            geometry->WriteILXL(true);
            geometry = NULL;
            failed = true;
          }
        }
        delete fullGeometry;
      }
      else {
        geometry->WriteILXL();
      }
      if(!failed) {
        modelSettings->setAreaParameters(geometry);
        ILXLGeometry = geometry;
      }
    }
    else {
      errText += tmpErrText;
      failed = true;
    }
  }
  if(!failed)
  {
    const SegyGeometry * areaParams = modelSettings->getAreaParameters();
    failed = timeSimbox->setArea(areaParams, errText);

    if(failed)
    {
      writeAreas(areaParams,timeSimbox,areaType);
      errText += "The specified AREA extends outside the surface(s).\n";
    }
    else
    {
      LogKit::LogFormatted(LogKit::Low,"\nResolution                x0           y0            lx         ly     azimuth         dx      dy\n");
      LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
      double azimuth = (-1)*timeSimbox->getAngle()*(180.0/M_PI);
      if (azimuth < 0)
        azimuth += 360.0;
      LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                           areaType.c_str(),
                           timeSimbox->getx0(), timeSimbox->gety0(),
                           timeSimbox->getlx(), timeSimbox->getly(), azimuth,
                           timeSimbox->getdx(), timeSimbox->getdy());
    }

    float minHorRes = modelSettings->getMinHorizontalRes();
    if (timeSimbox->getdx() < minHorRes || timeSimbox->getdy() < minHorRes){
      failed = true;
      errText += "The horizontal resolution in dx and dy should normally be above "+NRLib::ToString(minHorRes)
        +" m. If you need a denser\n sampling, please specify a new <advanced-settings><minimum-horizontal-resolution>\n";
    }

    if(!failed)
    {
      //
      // Set IL/XL information in geometry
      // ---------------------------------
      //
      // Skip for estimation mode if possible:
      //   a) For speed
      //   b) Grid data may not be available.
      if (modelSettings->getEstimationMode() == false || estimationModeNeedILXL == true) {
        if(ILXLGeometry == NULL) {
          int gridType = IO::findGridType(gridFile);
          bool ilxl_info_available = ((gridType == IO::SEGY) || (gridType == IO::CRAVA));
          if (ilxl_info_available) {
            LogKit::LogFormatted(LogKit::High,"\nFinding IL/XL information from grid data file \'"+gridFile+"\'\n");
            std::string tmpErrText;
            getGeometryFromGridOnFile(gridFile,
                                      modelSettings->getTraceHeaderFormat(0,0), //Trace header format is the same for all time lapses
                                      ILXLGeometry,
                                      tmpErrText);
            if(ILXLGeometry == NULL) {
              errText += tmpErrText;
              failed = true;
            }
          }
          else {
            LogKit::LogFormatted(LogKit::High,"\nCannot extract IL/XL information from non-SEGY grid data file \'"+gridFile+"\'\n");
          }
        }
        if(ILXLGeometry != NULL) {
          if(timeSimbox->isAligned(ILXLGeometry))
            timeSimbox->setILXL(ILXLGeometry);
          delete ILXLGeometry;
        }
      }

      // Rotate variograms relative to simbox
      modelSettings->rotateVariograms(static_cast<float> (timeSimbox_->getAngle()));

      //
      // Set SURFACES
      //

      setSimboxSurfaces(timeSimbox,
                        inputFiles->getTimeSurfFiles(),
                        modelSettings,
                        errText,
                        failed);

      if(!failed)
      {
        if(modelSettings->getUseLocalWavelet() && timeSimbox->getIsConstantThick())
        {
          LogKit::LogFormatted(LogKit::Warning,"\nWARNING: LOCALWAVELET is ignored when using constant thickness in DEPTH.\n");
          TaskList::addTask("If local wavelet is to be used, constant thickness in depth should be removed.");
        }


        int status = timeSimbox->calculateDz(modelSettings->getLzLimit(),errText);
        estimateZPaddingSize(timeSimbox, modelSettings);

        float minSampDens = modelSettings->getMinSamplingDensity();
        if (timeSimbox->getdz()*timeSimbox->getMinRelThick() < minSampDens){
          failed   = true;
          errText += "We normally discourage denser sampling than "+NRLib::ToString(minSampDens);
          errText += "ms in the time grid. If you really need\nthis, please use ";
          errText += "<project-settings><advanced-settings><minimum-sampling-density>\n";
        }

        if(status == Simbox::BOXOK)
        {
          logIntervalInformation(timeSimbox, "Time output interval:","Two-way-time");
          //
          // Make extended time simbox
          //
          if(inputFiles->getCorrDirFile() != "") {
            //
            // Get correlation direction
            //
            try {
              Surface tmpSurf(inputFiles->getCorrDirFile());
              if(timeSimbox->CheckSurface(tmpSurf) == true)
                correlationDirection = new Surface(tmpSurf);
              else {
                errText += "Error: Correlation surface does not cover volume.\n";
                failed = true;
              }
            }
            catch (NRLib::Exception & e) {
              errText += e.what();
              failed = true;
            }

            if(failed == false && modelSettings->getForwardModeling() == false) {
              //Extends timeSimbox for correlation coverage. Original stored in timeCutSimbox
              setupExtendedTimeSimbox(timeSimbox, correlationDirection,
                                      timeCutSimbox,
                                      modelSettings->getOutputGridFormat(),
                                      modelSettings->getOutputGridDomain(),
                                      modelSettings->getOtherOutputFlag());
            }

            estimateZPaddingSize(timeSimbox, modelSettings);

            status = timeSimbox->calculateDz(modelSettings->getLzLimit(),errText);

            if(status == Simbox::BOXOK)
              logIntervalInformation(timeSimbox, "Time inversion interval (extended relative to output interval due to correlation):","Two-way-time");
            else
            {
              errText += "Could not make the time simulation grid.\n";
              failed = true;
            }

            if(modelSettings->getForwardModeling() == false && failed == false) {
              setupExtendedBackgroundSimbox(timeSimbox, correlationDirection, timeBGSimbox,
                                            modelSettings->getOutputGridFormat(),
                                            modelSettings->getOutputGridDomain(),
                                            modelSettings->getOtherOutputFlag());
              status = timeBGSimbox->calculateDz(modelSettings->getLzLimit(),errText);
              if(status == Simbox::BOXOK)
                logIntervalInformation(timeBGSimbox, "Time interval used for background modelling:","Two-way-time");
              else
              {
                errText += "Could not make the grid for background model.\n";
                failed = true;
              }
            }
          }

          if(failed == false) {
            estimateXYPaddingSizes(timeSimbox, modelSettings);

            unsigned long long int gridsize = static_cast<unsigned long long int>(modelSettings->getNXpad())*modelSettings->getNYpad()*modelSettings->getNZpad();

            if(gridsize > std::numeric_limits<unsigned int>::max()) {
              float fsize = 4.0f*static_cast<float>(gridsize)/static_cast<float>(1024*1024*1024);
              float fmax  = 4.0f*static_cast<float>(std::numeric_limits<unsigned int>::max()/static_cast<float>(1024*1024*1024));
              errText += "Grids as large as "+NRLib::ToString(fsize,1)+"GB cannot be handled. The largest accepted grid size\n";
              errText += "is "+NRLib::ToString(fmax)+"GB. Please reduce the number of layers or the lateral resolution.\n";
              failed = true;
            }

            LogKit::LogFormatted(LogKit::Low,"\nTime simulation grids:\n");
            LogKit::LogFormatted(LogKit::Low,"  Output grid         %4i * %4i * %4i   : %10llu\n",
                                 timeSimbox->getnx(),timeSimbox->getny(),timeSimbox->getnz(),
                                 static_cast<unsigned long long int>(timeSimbox->getnx())*timeSimbox->getny()*timeSimbox->getnz());
            LogKit::LogFormatted(LogKit::Low,"  FFT grid            %4i * %4i * %4i   :%11llu\n",
                                 timeSimbox->GetNXpad(),timeSimbox->GetNYpad(),timeSimbox->GetNZpad(),
                                 static_cast<unsigned long long int>(modelSettings->getNXpad())*modelSettings->getNYpad()*modelSettings->getNZpad());
          }

          //
          // Make time simbox with constant thicknesses (needed for log filtering and facies probabilities)
          //
          timeSimboxConstThick = new Simbox(timeSimbox);
          Surface tsurf(dynamic_cast<const Surface &> (timeSimbox->GetTopSurface()));
          timeSimboxConstThick->setDepth(tsurf, 0, timeSimbox->getlz(), timeSimbox->getdz());

          if((modelSettings->getOtherOutputFlag() & IO::EXTRA_SURFACES) > 0 && (modelSettings->getOutputGridDomain() & IO::TIMEDOMAIN) > 0) {
            std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_ConstThick";
            std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_ConstThick";
            timeSimboxConstThick->writeTopBotGrids(topSurf,
                                                   baseSurf,
                                                   IO::PathToInversionResults(),
                                                   modelSettings->getOutputGridFormat());
          }

          if(timeSimbox->getdz() >= 10.0 && modelSettings->getFaciesProbFromRockPhysics() == true) {
            errText += "dz is too large to generate synthetic well data when estimating facies probabilities using rock physics models. Need dz < 10.";
            failed = true;
          }
        }
        else
        {
          errText += "Could not make time simulation grid.\n";
          failed = true;
        }
      }
      else
      {
        timeSimbox->externalFailure();
        failed = true;
      }
    }
  }

}



void
ModelGeneral::logIntervalInformation(const Simbox      * simbox,
                                     const std::string & header_text1,
                                     const std::string & header_text2)
{
  LogKit::LogFormatted(LogKit::Low,"\n"+header_text1+"\n");
  double zmin, zmax;
  simbox->getMinMaxZ(zmin,zmax);
  LogKit::LogFormatted(LogKit::Low," %13s          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       header_text2.c_str(),
                       zmin+simbox->getlz()*simbox->getAvgRelThick()*0.5,
                       zmin,zmax);
  LogKit::LogFormatted(LogKit::Low,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       simbox->getlz()*simbox->getAvgRelThick(),
                       simbox->getlz()*simbox->getMinRelThick(),
                       simbox->getlz());
  LogKit::LogFormatted(LogKit::Low,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n",
                       simbox->getdz()*simbox->getAvgRelThick(),
                       simbox->getdz(),
                       simbox->getdz()*simbox->getMinRelThick());
}
*/
void ModelGeneral::setSimboxSurfaces(Simbox                        *& simbox,
                                     const std::vector<std::string> & surfFile,
                                     ModelSettings                  * model_settings,
                                     std::string                    & errText,
                                     bool                           & failed)
{
  const std::string & topName = surfFile[0];

  bool   generateSeismic    = model_settings->getForwardModeling();
  bool   estimationMode     = model_settings->getEstimationMode();
  bool   generateBackground = model_settings->getGenerateBackground();
  bool   parallelSurfaces   = model_settings->getParallelTimeSurfaces();
  int    nz                 = model_settings->getTimeNz();
  int    outputFormat       = model_settings->getOutputGridFormat();
  int    outputDomain       = model_settings->getOutputGridDomain();
  int    outputGridsElastic = model_settings->getOutputGridsElastic();
  int    outputGridsOther   = model_settings->getOutputGridsOther();
  int    outputGridsSeismic = model_settings->getOutputGridsSeismic();
  double dTop               = model_settings->getTimeDTop();
  double lz                 = model_settings->getTimeLz();
  double dz                 = model_settings->getTimeDz();

  Surface * z0Grid = NULL;
  Surface * z1Grid = NULL;
  try {
    if (NRLib::IsNumber(topName)) {
      // Find the smallest surface that covers the simbox. For simplicity
      // we use only four nodes (nx=ny=2).
      double xMin, xMax;
      double yMin, yMax;
      findSmallestSurfaceGeometry(simbox->getx0(), simbox->gety0(),
                                  simbox->getlx(), simbox->getly(),
                                  simbox->getAngle(),
                                  xMin,yMin,xMax,yMax);
      z0Grid = new Surface(xMin-100, yMin-100, xMax-xMin+200, yMax-yMin+200, 2, 2, atof(topName.c_str()));
    }
    else {
      Surface tmpSurf(topName);
      z0Grid = new Surface(tmpSurf);
    }
  }
  catch (NRLib::Exception & e) {
    errText += e.what();
    failed = true;
  }

  if (!failed) {
    if (parallelSurfaces) { //Only one reference surface
      simbox->setDepth(*z0Grid, dTop, lz, dz, model_settings->getRunFromPanel());
    }
    else {
      const std::string & baseName = surfFile[1];
      try {
        if (NRLib::IsNumber(baseName)) {
          // Find the smallest surface that covers the simbox. For simplicity
          // we use only four nodes (nx=ny=2).
          double xMin, xMax;
          double yMin, yMax;
          findSmallestSurfaceGeometry(simbox->getx0(), simbox->gety0(),
                                      simbox->getlx(), simbox->getly(),
                                      simbox->getAngle(),
                                      xMin,yMin,xMax,yMax);
          z1Grid = new Surface(xMin-100, yMin-100, xMax-xMin+200, yMax-yMin+200, 2, 2, atof(baseName.c_str()));
        }
        else {
          Surface tmpSurf(baseName);
          z1Grid = new Surface(tmpSurf);
        }
      }
      catch (NRLib::Exception & e) {
        errText += e.what();
        failed = true;
      }
      if (!failed) {
        try {
          simbox->setDepth(*z0Grid, *z1Grid, nz, model_settings->getRunFromPanel());
        }
        catch (NRLib::Exception & e) {
          errText += e.what();
          std::string text("Seismic data");
          writeAreas(model_settings->getAreaParameters(),simbox,text);
          failed = true;
        }
      }
    }
    if (!failed) {
      if ((outputDomain & IO::TIMEDOMAIN) > 0) {
        std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime();
        std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
        simbox->setTopBotName(topSurf,baseSurf,outputFormat);
        if (generateSeismic) {
          simbox->writeTopBotGrids(topSurf,
                                   baseSurf,
                                   IO::PathToSeismicData(),
                                   outputFormat);
        }
        else if (!estimationMode){
          if (outputGridsElastic > 0 || outputGridsOther > 0 || outputGridsSeismic > 0)
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToInversionResults(),
                                     outputFormat);
        }
        if ((outputFormat & IO::STORM) > 0) { // These copies are only needed with the STORM format
          if ((outputGridsElastic & IO::BACKGROUND) > 0 ||
              (outputGridsElastic & IO::BACKGROUND_TREND) > 0 ||
              (estimationMode && generateBackground)) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToBackground(),
                                     outputFormat);
          }
          if ((outputGridsOther & IO::CORRELATION) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToCorrelations(),
                                     outputFormat);
          }
          if ((outputGridsSeismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToSeismicData(),
                                     outputFormat);
          }
          if ((outputGridsOther & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToVelocity(),
                                     outputFormat);
          }
        }
      }
    }
  }
  delete z0Grid;
  delete z1Grid;
}

//void
//ModelGeneral::setupExtendedTimeSimbox(Simbox   * timeSimbox,
//                                      Surface  * corrSurf,
//                                      Simbox  *& timeCutSimbox,
//                                      int        outputFormat,
//                                      int        outputDomain,
//                                      int        otherOutput)
//{
//  timeCutSimbox = new Simbox(timeSimbox);
//  NRLib::Vector corrPlanePars = findPlane(corrSurf);
//  Surface * meanSurf;
//  if (corrSurf->GetNI() > 2)
//    meanSurf = new Surface(*corrSurf);
//  else {
//    meanSurf = new Surface(dynamic_cast<const Surface &>(timeSimbox->GetTopSurface()));
//    if (meanSurf->GetNI() == 2) { //Extend corrSurf to cover other surfaces.
//      double minX = meanSurf->GetXMin();
//      double maxX = meanSurf->GetXMax();
//      double minY = meanSurf->GetYMin();
//      double maxY = meanSurf->GetYMax();
//      if (minX > corrSurf->GetXMin())
//        minX = corrSurf->GetXMin();
//      if (maxX < corrSurf->GetXMax())
//        maxX = corrSurf->GetXMax();
//      if (minY > corrSurf->GetYMin())
//        minY = corrSurf->GetYMin();
//      if (maxY < corrSurf->GetYMax())
//        maxY = corrSurf->GetYMax();
//      corrSurf->SetDimensions(minX, minY, maxX-minX, maxY-minY);
//    }
//  }
//  int i;
//  for (i=0;i<static_cast<int>(meanSurf->GetN());i++)
//    (*meanSurf)(i) = 0;
//
//  meanSurf->AddNonConform(&(timeSimbox->GetTopSurface()));
//  meanSurf->AddNonConform(&(timeSimbox->GetBotSurface()));
//  meanSurf->Multiply(0.5);
//  NRLib::Vector refPlanePars = findPlane(meanSurf);
//
//  refPlanePars -= corrPlanePars;
//
//  gradX_ = refPlanePars(1);
//  gradY_ = refPlanePars(2);
//
//  Surface * refPlane = createPlaneSurface(refPlanePars, meanSurf);
//  delete meanSurf;
//  meanSurf = NULL;
//
//  std::string fileName = "Correlation_Rotation_Plane";
//  IO::writeSurfaceToFile(*refPlane, fileName, IO::PathToCorrelations(), outputFormat);
//
//  refPlane->AddNonConform(corrSurf);
//
//  Surface topSurf(*refPlane);
//  topSurf.SubtractNonConform(&(timeSimbox->GetTopSurface()));
//  double shiftTop = topSurf.Max();
//  shiftTop *= -1.0;
//  topSurf.Add(shiftTop);
//  topSurf.AddNonConform(&(timeSimbox->GetTopSurface()));
//
//  Surface botSurf(*refPlane);
//  botSurf.SubtractNonConform(&(timeSimbox->GetBotSurface()));
//  double shiftBot = botSurf.Min();
//  shiftBot *= -1.0;
//  double thick    = shiftBot-shiftTop;
//  double dz       = timeCutSimbox->getdz();
//  int    nz       = int(thick/dz);
//  double residual = thick - nz*dz;
//  if (residual > 0.0) {
//    shiftBot += dz-residual;
//    nz++;
//  }
//  if (nz != timeCutSimbox->getnz()) {
//    LogKit::LogFormatted(LogKit::High,"\nNumber of layers in inversion increased from %d",timeCutSimbox->getnz());
//    LogKit::LogFormatted(LogKit::High," to %d in grid created using correlation direction.\n",nz);
//  }
//  botSurf.Add(shiftBot);
//  botSurf.AddNonConform(&(timeSimbox->GetBotSurface()));
//
//  timeSimbox->setDepth(topSurf, botSurf, nz);
//
//  if ((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
//    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Extended";
//    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Extended";
//    timeSimbox->writeTopBotGrids(topSurf,
//                                 baseSurf,
//                                 IO::PathToInversionResults(),
//                                 outputFormat);
//  }
//
//  delete refPlane;
//}

//void
//ModelGeneral::setupExtendedBackgroundSimbox(Simbox   * timeSimbox,
//                                            Surface  * corrSurf,
//                                            Simbox  *& timeBGSimbox,
//                                            int        outputFormat,
//                                            int        outputDomain,
//                                            int        otherOutput)
//{
//  //
//  // Move correlation surface for easier handling.
//  //
//  Surface tmpSurf(*corrSurf);
//  double avg = tmpSurf.Avg();
//  if (avg > 0)
//    tmpSurf.Subtract(avg);
//  else
//    tmpSurf.Add(avg); // This situation is not very likely, but ...
//
//  //
//  // Find top surface of background simbox.
//  //
//  // The funny/strange dTop->Multiply(-1.0) is due to NRLIB's current
//  // inability to set dTop equal to Simbox top surface.
//  //
//  Surface dTop(tmpSurf);
//  dTop.SubtractNonConform(&(timeSimbox->GetTopSurface()));
//  dTop.Multiply(-1.0);
//  double shiftTop = dTop.Min();
//  Surface topSurf(tmpSurf);
//  topSurf.Add(shiftTop);
//
//  //
//  // Find base surface of background simbox
//  //
//  Surface dBot(tmpSurf);
//  dBot.SubtractNonConform(&(timeSimbox->GetBotSurface()));
//  dBot.Multiply(-1.0);
//  double shiftBot = dBot.Max();
//  Surface botSurf(tmpSurf);
//  botSurf.Add(shiftBot);
//
//  //
//  // Calculate number of layers of background simbox
//  //
//  tmpSurf.Assign(0.0);
//  tmpSurf.AddNonConform(&botSurf);
//  tmpSurf.SubtractNonConform(&topSurf);
//  double dMax = tmpSurf.Max();
//  double dt = timeSimbox->getdz();
//  int nz;
//  //
//  // NBNB-PAL: I think it is a good idea to use a maximum dt of 10ms.
//  //
//  //if (dt < 10.0) {
//  //  LogKit::LogFormatted(LogKit::High,"\nReducing sampling density for background",dt);
//  //  LogKit::LogFormatted(LogKit::High," modelling from %.2fms to 10.0ms\n");
//  //  dt = 10.0;  // A sampling density of 10.0ms is good enough for BG model
//  // }
//  nz = static_cast<int>(ceil(dMax/dt));
//
//  //
//  // Make new simbox
//  //
//  timeBGSimbox = new Simbox(timeSimbox);
//  timeBGSimbox->setDepth(topSurf, botSurf, nz);
//
//  if ((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
//    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_BG";
//    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_BG";
//    timeBGSimbox->writeTopBotGrids(topSurf,
//                                   baseSurf,
//                                   IO::PathToBackground(),
//                                   outputFormat);
//  }
//}

NRLib::Vector
ModelGeneral::findPlane(Surface * surf)
{
  NRLib::SymmetricMatrix A = NRLib::SymmetricZeroMatrix(3);
  NRLib::Vector b(3);
  NRLib::Vector x(3);

  b = 0;

  int nData = 0;

  for (int i=0 ; i<static_cast<int>(surf->GetN()) ; i++) {
    double x, y, z;
    surf->GetXY(i, x, y);
    z = (*surf)(i);
    if (!surf->IsMissing(z)) {
      nData++;
      A(0,1) += x;
      A(0,2) += y;
      A(1,1) += x*x;
      A(1,2) += x*y;
      A(2,2) += y*y;
      b(0)   += z;
      b(1)   += x*z;
      b(2)   += y*z;
    }
  }

  A(0,0) = nData;

  NRLib::CholeskySolve(A, b, x);

  return x;
}


Surface *
ModelGeneral::createPlaneSurface(const NRLib::Vector & planeParams,
                                 Surface             * templateSurf)
{
  Surface * result = new Surface(*templateSurf);
  for (int i=0;i<static_cast<int>(result->GetN());i++) {
    double x,y;
    result->GetXY(i,x,y);
    (*result)(i) = planeParams(0)+planeParams(1)*x+planeParams(2)*y;
  }
  return(result);
}


void
ModelGeneral::estimateXYPaddingSizes(Simbox         * timeSimbox,
                                     ModelSettings *& model_settings)
{
  double dx      = timeSimbox->getdx();
  double dy      = timeSimbox->getdy();
  double lx      = timeSimbox->getlx();
  double ly      = timeSimbox->getly();
  int    nx      = timeSimbox->getnx();
  int    ny      = timeSimbox->getny();
  int    nz      = timeSimbox->getnz();

  double xPadFac = model_settings->getXPadFac();
  double yPadFac = model_settings->getYPadFac();
  double xPad    = xPadFac*lx;
  double yPad    = yPadFac*ly;

  if (model_settings->getEstimateXYPadding())
  {
    float  range1 = model_settings->getLateralCorr()->getRange();
    float  range2 = model_settings->getLateralCorr()->getSubRange();
    float  angle  = model_settings->getLateralCorr()->getAngle();
    double factor = 0.5;  // Lateral correlation is not very important. Half a range is probably more than enough

    xPad          = factor * std::max(fabs(range1*cos(angle)), fabs(range2*sin(angle)));
    yPad          = factor * std::max(fabs(range1*sin(angle)), fabs(range2*cos(angle)));
    xPad          = std::max(xPad, dx);     // Always require at least on grid cell
    yPad          = std::max(yPad, dy);     // Always require at least one grid cell
    xPadFac       = std::min(1.0, xPad/lx); // A padding of more than 100% is insensible
    yPadFac       = std::min(1.0, yPad/ly);
  }

  int nxPad = setPaddingSize(nx, xPadFac);
  int nyPad = setPaddingSize(ny, yPadFac);
  int nzPad = timeSimbox->GetNZpad();

  double true_xPadFac = static_cast<double>(nxPad - nx)/static_cast<double>(nx);
  double true_yPadFac = static_cast<double>(nyPad - ny)/static_cast<double>(ny);
  double true_zPadFac = model_settings->getZPadFac();
  double true_xPad    = true_xPadFac*lx;
  double true_yPad    = true_yPadFac*ly;
  double true_zPad    = true_zPadFac*(timeSimbox->getlz()*timeSimbox->getMinRelThick());

  timeSimbox->SetNXpad(nxPad);
  timeSimbox->SetNYpad(nyPad);
  timeSimbox->SetXPadFactor(true_xPadFac);
  timeSimbox->SetYPadFactor(true_yPadFac);

  std::string text1;
  std::string text2;
  int log_level = LogKit::Medium;
  if (model_settings->getEstimateXYPadding()) {
    text1 = " estimated from lateral correlation ranges in internal grid";
    log_level = LogKit::Low;
  }
  if (model_settings->getEstimateZPadding()) {
    text2 = " estimated from an assumed wavelet length";
    log_level = LogKit::Low;
  }

  LogKit::LogFormatted(log_level,"\nPadding sizes"+text1+":\n");
  LogKit::LogFormatted(log_level,"  xPad, xPadFac, nx, nxPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_xPad, true_xPadFac, nx, nxPad);
  LogKit::LogFormatted(log_level,"  yPad, yPadFac, ny, nyPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_yPad, true_yPadFac, ny, nyPad);
  LogKit::LogFormatted(log_level,"\nPadding sizes"+text2+":\n");
  LogKit::LogFormatted(log_level,"  zPad, zPadFac, nz, nzPad                 : %5.fms, %5.3f, %5d, %4d\n",
                       true_zPad, true_zPadFac, nz, nzPad);
}

void
ModelGeneral::estimateZPaddingSize(Simbox         * timeSimbox,
                                   ModelSettings *& model_settings)
{
  int    nz          = timeSimbox->getnz();
  double minLz       = timeSimbox->getlz()*timeSimbox->getMinRelThick();
  double zPadFac     = model_settings->getZPadFac();
  double zPad        = zPadFac*minLz;

  if (model_settings->getEstimateZPadding())
  {
    double wLength = static_cast<double>(model_settings->getDefaultWaveletLength());
    double pfac    = 1.0;
    zPad           = wLength/pfac;                             // Use half a wavelet as padding
    zPadFac        = std::min(1.0, zPad/minLz);                // More than 100% padding is not sensible
  }
  int nzPad        = setPaddingSize(nz, zPadFac);
  zPadFac          = static_cast<double>(nzPad - nz)/static_cast<double>(nz);

  timeSimbox->SetNZpad(nzPad);
  timeSimbox->SetZPadFactor(zPadFac);
}

//void
//ModelGeneral::readGridFromFile(const std::string       & fileName,
//                               const std::string       & parName,
//                               const float               offset,
//                               FFTGrid                *& grid,
//                               const SegyGeometry     *& geometry,
//                               const TraceHeaderFormat * format,
//                               int                       gridType,
//                               const Simbox            * timeSimbox,
//                               const Simbox            * timeCutSimbox,
//                               const ModelSettings     * model_settings,
//                               std::string             & errText,
//                               bool                      nopadding)
//{
//  int fileType = IO::findGridType(fileName);
//
//  if (fileType == IO::CRAVA)
//  {
//    int nxPad, nyPad, nzPad;
//    if (nopadding)
//    {
//      nxPad = timeSimbox->getnx();
//      nyPad = timeSimbox->getny();
//      nzPad = timeSimbox->getnz();
//    }
//    else
//    {
//      nxPad = model_settings->getNXpad();
//      nyPad = model_settings->getNYpad();
//      nzPad = model_settings->getNZpad();
//    }
//    LogKit::LogFormatted(LogKit::Low,"\nReading grid \'"+parName+"\' from file "+fileName);
//    grid = createFFTGrid(timeSimbox->getnx(),
//                         timeSimbox->getny(),
//                         timeSimbox->getnz(),
//                         nxPad,
//                         nyPad,
//                         nzPad,
//                         model_settings->getFileGrid());
//
//    grid->setType(gridType);
//    grid->readCravaFile(fileName, errText, nopadding);
//  }
//  else if (fileType == IO::SEGY)
//    readSegyFile(fileName, grid, timeSimbox, timeCutSimbox, model_settings, geometry,
//                 gridType, parName, offset, format, errText, nopadding);
//  else if (fileType == IO::STORM)
//    readStormFile(fileName, grid, gridType, parName, timeSimbox, model_settings, errText, false, nopadding);
//  else if (fileType == IO::SGRI)
//    readStormFile(fileName, grid, gridType, parName, timeSimbox, model_settings, errText, true, nopadding);
//  else
//  {
//    errText += "\nReading of file \'"+fileName+"\' for grid type \'"+parName+"\'failed. File type not recognized.\n";
//  }
//
//}

void
ModelGeneral::printSettings(ModelSettings     * model_settings,
                            const InputFiles  * input_files)
{
  LogKit::WriteHeader("Model settings");

  LogKit::LogFormatted(LogKit::Low,"\nGeneral settings:\n");
  if (model_settings->getForwardModeling()==true)
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : forward\n");
  else if (model_settings->getEstimationMode()==true)
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : estimation\n");

  else if (model_settings->getNumberOfSimulations() == 0)
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : prediction\n");
  else
  {
    LogKit::LogFormatted(LogKit::Low,"  Modelling mode                           : simulation\n");
    if (input_files->getSeedFile()=="") {
      if (model_settings->getSeed() == 0)
        LogKit::LogFormatted(LogKit::Low,"  Seed                                     :          0 (default seed)\n");
      else
        LogKit::LogFormatted(LogKit::Low,"  Seed                                     : %10d\n",model_settings->getSeed());
    }
    else
      LogKit::LogFormatted(LogKit::Low,"  Seed read from file                      : %10s\n",input_files->getSeedFile().c_str());


    LogKit::LogFormatted(LogKit::Low,"  Number of realisations                   : %10d\n",model_settings->getNumberOfSimulations());
  }
  if (model_settings->getForwardModeling()==false)
  {
    LogKit::LogFormatted(LogKit::Low,"  Kriging                                  : %10s\n",(model_settings->getKrigingParameter()>0 ? "yes" : "no"));
    LogKit::LogFormatted(LogKit::Low,"  Facies probabilities                     : %10s\n",(model_settings->getEstimateFaciesProb() ? "yes" : "no"));
    LogKit::LogFormatted(LogKit::Low,"  Synthetic seismic                        : %10s\n",(model_settings->getGenerateSeismicAfterInv() ? "yes" : "no" ));
  }

  if (model_settings->getEstimateFaciesProb()) {
    LogKit::LogFormatted(LogKit::Low,"\nSettings for facies probability estimation:\n");
    LogKit::LogFormatted(LogKit::Low,"  Use elastic parameters relative to trend : %10s\n",(model_settings->getFaciesProbRelative()     ? "yes" : "no"));
    LogKit::LogFormatted(LogKit::Low,"  Include Vs information in estimation     : %10s\n",(model_settings->getNoVsFaciesProb()         ? "no"  : "yes"));
    LogKit::LogFormatted(LogKit::Low,"  Use filtered well logs for estimation    : %10s\n",(model_settings->getUseFilterForFaciesProb() ? "yes" : "no"));
  }

  LogKit::LogFormatted(LogKit::Low,"\nInput/Output settings:\n");
  std::string log_text("*NONE*");
  int log_level = model_settings->getLogLevel();
  if (log_level == LogKit::L_Error)
    log_text = "ERROR";
  else if (log_level == LogKit::L_Warning)
    log_text = "WARNING";
  else if (log_level == LogKit::L_Low)
    log_text = "LOW";
  else if (log_level == LogKit::L_Medium)
    log_text = "MEDIUM";
  else if (log_level == LogKit::L_High)
    log_text = "HIGH";
  else if (log_level == LogKit::L_DebugLow)
     log_text = "DEBUGLOW";
  else if (log_level == LogKit::L_DebugHigh)
    log_text = "DEBUGHIGH";
  LogKit::LogFormatted(LogKit::Low, "  Log level                                : %10s\n",log_text.c_str());
  if (input_files->getInputDirectory() != "")
    LogKit::LogFormatted(LogKit::High,"  Input directory                          : %10s\n",input_files->getInputDirectory().c_str());
  if (IO::getOutputPath() != "")
    LogKit::LogFormatted(LogKit::High,"  Output directory                         : %10s\n",IO::getOutputPath().c_str());

  int grid_format          = model_settings->getOutputGridFormat();
  int grid_domain          = model_settings->getOutputGridDomain();
  int output_grids_other   = model_settings->getOutputGridsOther();
  int output_grids_elastic = model_settings->getOutputGridsElastic();
  int output_grids_seismic = model_settings->getOutputGridsSeismic();

  if (output_grids_elastic > 0  || output_grids_seismic > 0  || output_grids_other > 0) {
    LogKit::LogFormatted(LogKit::Medium,"\nGrid output formats:\n");
    if (grid_format & IO::SEGY) {
      const std::string & format_name = model_settings->getTraceHeaderFormatOutput()->GetFormatName();
      LogKit::LogFormatted(LogKit::Medium,"  Segy - %-10s                        :        yes\n",format_name.c_str());
    }
    if (grid_format & IO::STORM)
      LogKit::LogFormatted(LogKit::Medium,"  Storm                                    :        yes\n");
    if (grid_format & IO::ASCII)
      LogKit::LogFormatted(LogKit::Medium,"  ASCII                                    :        yes\n");
    if (grid_format & IO::SGRI)
      LogKit::LogFormatted(LogKit::Medium,"  Norsar                                   :        yes\n");
    if (grid_format & IO::CRAVA)
      LogKit::LogFormatted(LogKit::Medium,"  Crava                                    :        yes\n");

    LogKit::LogFormatted(LogKit::Medium,"\nGrid output domains:\n");
    if (grid_domain & IO::TIMEDOMAIN)
      LogKit::LogFormatted(LogKit::Medium,"  Time                                     :        yes\n");
    if (grid_domain & IO::DEPTHDOMAIN)
      LogKit::LogFormatted(LogKit::Medium,"  Depth                                    :        yes\n");
  }

  if (output_grids_elastic > 0 &&
      model_settings->getForwardModeling() == false) {
    LogKit::LogFormatted(LogKit::Medium,"\nOutput of elastic parameters:\n");
    if ((output_grids_elastic & IO::VP) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Pressure-wave velocity  (Vp)             :        yes\n");
    if ((output_grids_elastic & IO::VS) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Shear-wave velocity  (Vs)                :        yes\n");
    if ((output_grids_elastic & IO::RHO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Density  (Rho)                           :        yes\n");
    if ((output_grids_elastic & IO::AI) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Acoustic impedance  (AI)                 :        yes\n");
    if ((output_grids_elastic & IO::VPVSRATIO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Vp/Vs ratio                              :        yes\n");
    if ((output_grids_elastic & IO::SI) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Shear impedance  (SI)                    :        yes\n");
    if ((output_grids_elastic & IO::MURHO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  MuRho  (SI*SI)                           :        yes\n");
    if ((output_grids_elastic & IO::LAMBDARHO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  LambdaRho  (AI*AI - 2*SI*SI)             :        yes\n");
    if ((output_grids_elastic & IO::LAMELAMBDA) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Lame's first parameter                   :        yes\n");
    if ((output_grids_elastic & IO::LAMEMU) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Lame's second parameter (shear modulus)  :        yes\n");
    if ((output_grids_elastic & IO::POISSONRATIO) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Poisson ratio  (X-1)/2(X-2), X=(Vp/Vs)^2 :        yes\n");
    if ((output_grids_elastic & IO::BACKGROUND) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Background (Vp, Vs, Rho)                 :        yes\n");
    if ((output_grids_elastic & IO::BACKGROUND_TREND) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Background trend (Vp, Vs, Rho)           :        yes\n");
  }

  if (model_settings->getForwardModeling() ||
      output_grids_seismic > 0) {
    LogKit::LogFormatted(LogKit::Medium,"\nOutput of seismic data:\n");
    if ((output_grids_seismic & IO::SYNTHETIC_SEISMIC_DATA) > 0 || model_settings->getForwardModeling())
      LogKit::LogFormatted(LogKit::Medium,"  Synthetic seismic data (forward modelled):        yes\n");
    if ((output_grids_seismic & IO::ORIGINAL_SEISMIC_DATA) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Original seismic data (in output grid)   :        yes\n");
    if ((output_grids_seismic & IO::RESIDUAL) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Seismic data residuals                   :        yes\n");
  }

  if (model_settings->getEstimateFaciesProb()) {
    LogKit::LogFormatted(LogKit::Medium,"\nOutput of facies probability volumes:\n");
    if ((output_grids_other & IO::FACIESPROB) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Facies probabilities                     :        yes\n");
    if ((output_grids_other & IO::FACIESPROB_WITH_UNDEF) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Facies probabilities with undefined value:        yes\n");
  }

  if ((output_grids_other & IO::CORRELATION)>0 ||
      (output_grids_other & IO::EXTRA_GRIDS)  >0 ||
      (output_grids_other & IO::TIME_TO_DEPTH_VELOCITY)>0) {
    LogKit::LogFormatted(LogKit::Medium,"\nOther grid output:\n");
    if ((output_grids_other & IO::CORRELATION) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Posterior correlations                   :        yes\n");
    if ((output_grids_other & IO::EXTRA_GRIDS) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Help grids (see use manual)              :        yes\n");
    if ((output_grids_other & IO::TIME_TO_DEPTH_VELOCITY) > 0)
      LogKit::LogFormatted(LogKit::Medium,"  Time-to-depth velocity                   :        yes\n");
  }

  if (model_settings->getFileGrid())
    LogKit::LogFormatted(LogKit::Medium,"\nAdvanced settings:\n");
  else
    LogKit::LogFormatted(LogKit::High,"\nAdvanced settings:\n");

  LogKit::LogFormatted(LogKit::Medium, "  Use intermediate disk storage for grids  : %10s\n", (model_settings->getFileGrid() ? "yes" : "no"));

  if (input_files->getReflMatrFile() != "")
    LogKit::LogFormatted(LogKit::Medium, "  Take reflection matrix from file         : %10s\n", input_files->getReflMatrFile().c_str());

  if (model_settings->getVpVsRatio() != RMISSING)
    LogKit::LogFormatted(LogKit::High ,"  Vp-Vs ratio used in reflection coef.     : %10.2f\n", model_settings->getVpVsRatio());

  LogKit::LogFormatted(LogKit::High, "  RMS panel mode                           : %10s\n"  , (model_settings->getRunFromPanel() ? "yes" : "no"));
  LogKit::LogFormatted(LogKit::High ,"  Smallest allowed length increment (dxy)  : %10.2f\n", model_settings->getMinHorizontalRes());
  LogKit::LogFormatted(LogKit::High ,"  Smallest allowed time increment (dt)     : %10.2f\n", model_settings->getMinSamplingDensity());

  if (model_settings->getKrigingParameter()>0) { // We are doing kriging
    LogKit::LogFormatted(LogKit::High ,"  Data in neighbourhood when doing kriging : %10.2f\n", model_settings->getKrigingParameter());
    LogKit::LogFormatted(LogKit::High, "  Smooth kriged parameters                 : %10s\n", (model_settings->getDoSmoothKriging() ? "yes" : "no"));
  }

  LogKit::LogFormatted(LogKit::High,"\nUnit settings/assumptions:\n");
  LogKit::LogFormatted(LogKit::High,"  Time                                     : %10s\n","ms TWT");
  LogKit::LogFormatted(LogKit::High,"  Frequency                                : %10s\n","Hz");
  LogKit::LogFormatted(LogKit::High,"  Length                                   : %10s\n","m");
  LogKit::LogFormatted(LogKit::High,"  Velocities                               : %10s\n","m/s");
  LogKit::LogFormatted(LogKit::High,"  Density                                  : %10s\n","g/cm3");
  LogKit::LogFormatted(LogKit::High,"  Angles                                   : %10s\n","   degrees (clockwise relative to north when applicable)");

  //
  // WELL PROCESSING
  //
  if (model_settings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::High,"\nSettings for well processing:\n");
    LogKit::LogFormatted(LogKit::High,"  Threshold for merging log entries        : %10.2f ms\n",model_settings->getMaxMergeDist());
    LogKit::LogFormatted(LogKit::High,"  Threshold for Vp-Vs rank correlation     : %10.2f\n",model_settings->getMaxRankCorr());
    LogKit::LogFormatted(LogKit::High,"  Threshold for deviation angle            : %10.1f (=%.2fm/ms TWT)\n",
                         model_settings->getMaxDevAngle(),tan(model_settings->getMaxDevAngle()*M_PI/180.0));
    LogKit::LogFormatted(LogKit::High,"  High cut for background modelling        : %10.1f\n",model_settings->getMaxHzBackground());
    LogKit::LogFormatted(LogKit::High,"  High cut for seismic resolution          : %10.1f\n",model_settings->getMaxHzSeismic());
    LogKit::LogFormatted(LogKit::High,"  Estimate Vp-Vs ratio from well data      : %10s\n", (model_settings->getVpVsRatioFromWells() ? "yes" : "no"));
  }
  LogKit::LogFormatted(LogKit::High,"\nRange of allowed parameter values:\n");
  LogKit::LogFormatted(LogKit::High,"  Vp  - min                                : %10.0f\n",model_settings->getAlphaMin());
  LogKit::LogFormatted(LogKit::High,"  Vp  - max                                : %10.0f\n",model_settings->getAlphaMax());
  LogKit::LogFormatted(LogKit::High,"  Vs  - min                                : %10.0f\n",model_settings->getBetaMin());
  LogKit::LogFormatted(LogKit::High,"  Vs  - max                                : %10.0f\n",model_settings->getBetaMax());
  LogKit::LogFormatted(LogKit::High,"  Rho - min                                : %10.1f\n",model_settings->getRhoMin());
  LogKit::LogFormatted(LogKit::High,"  Rho - max                                : %10.1f\n",model_settings->getRhoMax());

  //
  // WELL DATA
  //
  if (model_settings->getNumberOfWells() > 0)
  {
    LogKit::LogFormatted(LogKit::Low,"\nWell logs:\n");
    const std::vector<std::string> & log_names = model_settings->getLogNames();

    if (log_names.size() > 0)
    {
      LogKit::LogFormatted(LogKit::Low,"  Time                                     : %10s\n",  log_names[0].c_str());
      if (NRLib::Uppercase(log_names[1])=="VP" ||
         NRLib::Uppercase(log_names[1])=="LFP_VP")
        LogKit::LogFormatted(LogKit::Low,"  p-wave velocity                          : %10s\n",log_names[1].c_str());
      else
        LogKit::LogFormatted(LogKit::Low,"  Sonic                                    : %10s\n",log_names[1].c_str());
      if (NRLib::Uppercase(log_names[3])=="VS" ||
         NRLib::Uppercase(log_names[3])=="LFP_VS")
        LogKit::LogFormatted(LogKit::Low,"  s-wave velocity                          : %10s\n",log_names[3].c_str());
      else
        LogKit::LogFormatted(LogKit::Low,"  Shear sonic                              : %10s\n",log_names[3].c_str());
      LogKit::LogFormatted(LogKit::Low,"  Density                                  : %10s\n",  log_names[2].c_str());
      if (model_settings->getFaciesLogGiven())
        LogKit::LogFormatted(LogKit::Low,"  Facies                                   : %10s\n",log_names[4].c_str());
      if (model_settings->getPorosityLogGiven())
        LogKit::LogFormatted(LogKit::Low,"  Porosity                                  : %10s\n",  log_names[5].c_str());
    }
    else
    {
      LogKit::LogFormatted(LogKit::Low,"  Time                                     : %10s\n","TWT");
      LogKit::LogFormatted(LogKit::Low,"  Sonic                                    : %10s\n","DT");
      LogKit::LogFormatted(LogKit::Low,"  Shear sonic                              : %10s\n","DTS");
      LogKit::LogFormatted(LogKit::Low,"  Density                                  : %10s\n","RHOB");
      LogKit::LogFormatted(LogKit::Low,"  Facies                                   : %10s\n","FACIES");
    }
    LogKit::LogFormatted(LogKit::Low,"\nWell files:\n");
    for (int i = 0; i < model_settings->getNumberOfWells(); i++)
    {
      LogKit::LogFormatted(LogKit::Low,"  %-2d                                       : %s\n",i+1,input_files->getWellFile(i).c_str());
    }
    bool generate_background  = model_settings->getGenerateBackground();
    bool estimate_facies_prob = model_settings->getFaciesLogGiven();
    bool estimate_wavelet     = false;
    for (int i=0; i < model_settings->getNumberOfTimeLapses(); i++){
      std::vector<bool> estimate_wavelet_all_traces = model_settings->getEstimateWavelet(i);
      for (int j = 0; j < model_settings->getNumberOfAngles(i); j++)
        estimate_wavelet = estimate_wavelet || estimate_wavelet_all_traces[j];
    }
    if (generate_background || estimate_facies_prob || estimate_wavelet)
    {
      LogKit::LogFormatted(LogKit::Low,"\nUse well in estimation of:                   ");
      if (generate_background)  LogKit::LogFormatted(LogKit::Low,"BackgroundTrend  ");
      if (estimate_wavelet)     LogKit::LogFormatted(LogKit::Low,"WaveletEstimation  ");
      if (estimate_facies_prob) LogKit::LogFormatted(LogKit::Low,"FaciesProbabilities");
      LogKit::LogFormatted(LogKit::Low,"\n");
      for (int i = 0; i < model_settings->getNumberOfWells(); i++)
      {
        LogKit::LogFormatted(LogKit::Low,"  %-2d                                       : ",i+1);
        if (generate_background) {
          if (model_settings->getIndicatorBGTrend(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::Low,"    %-11s  ","yes");
          else if (model_settings->getIndicatorBGTrend(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::Low,"    %-11s  ","no");
          else
            LogKit::LogFormatted(LogKit::Low,"    %-11s  ","yes");
        }
        if (estimate_wavelet) {
          if (model_settings->getIndicatorWavelet(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::Low,"    %-13s  ","yes");
          else if (model_settings->getIndicatorWavelet(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::Low,"    %-13s  ","no");
          else
            LogKit::LogFormatted(LogKit::Low,"    %-12s  ","if possible");
        }
        if (estimate_facies_prob) {
          if (model_settings->getIndicatorFacies(i) == ModelSettings::YES)
            LogKit::LogFormatted(LogKit::Low,"    %-12s","yes");
          else if (model_settings->getIndicatorFacies(i) == ModelSettings::NO)
            LogKit::LogFormatted(LogKit::Low,"    %-12s","no");
          else
            LogKit::LogFormatted(LogKit::Low,"    %-12s","if possible");
        }
        LogKit::LogFormatted(LogKit::Low,"\n");
      }
    }
    if ( model_settings->getOptimizeWellLocation() )
    {
      LogKit::LogFormatted(LogKit::Low,"\nFor well, optimize position for            : Angle with Weight\n");
      for (int i = 0; i < model_settings->getNumberOfWells(); i++)
      {
        int n_move_angles = model_settings->getNumberOfWellAngles(i);
        if ( n_move_angles > 0 )
        {
          LogKit::LogFormatted(LogKit::Low," %2d %46.1f %10.1f\n",i+1,(model_settings->getWellMoveAngle(i,0)*180/M_PI),model_settings->getWellMoveWeight(i,0));
          for (int j=1; j < n_move_angles; j++)
            LogKit::LogFormatted(LogKit::Low," %49.1f %10.1f\n",(model_settings->getWellMoveAngle(i,j)*180/M_PI),model_settings->getWellMoveWeight(i,j));
        }
        LogKit::LogFormatted(LogKit::Low,"\n");
      }
    }
  }

  //
  // AREA
  //
  std::string grid_file;
  int area_specification = model_settings->getAreaSpecification();
  if (model_settings->getForwardModeling()) {
    LogKit::LogFormatted(LogKit::Low,"\nSeismic area:\n");
    grid_file = input_files->getBackFile(0);    // Get geometry from earth model (Vp)
  }
  else {
    LogKit::LogFormatted(LogKit::Low,"\nInversion area");
    if (area_specification == ModelSettings::AREA_FROM_GRID_DATA ||
       area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM ||
       area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE)
      grid_file = input_files->getSeismicFile(0,0); // Get area from first seismic data volume
  }
  if (area_specification == ModelSettings::AREA_FROM_GRID_DATA) {
    const std::vector<int> & area_ILXL = model_settings->getAreaILXL();
    LogKit::LogFormatted(LogKit::Low," taken from grid\n");
    LogKit::LogFormatted(LogKit::Low,"  Grid                                     : "+grid_file+"\n");
    if (area_ILXL.size() > 0) {
      if (area_ILXL[0] != IMISSING)
        LogKit::LogFormatted(LogKit::Low,"  In-line start                            : %10d\n", area_ILXL[0]);
      if (area_ILXL[1] != IMISSING)
        LogKit::LogFormatted(LogKit::Low,"  In-line end                              : %10d\n", area_ILXL[1]);
      if (area_ILXL[4] != IMISSING)
        LogKit::LogFormatted(LogKit::Low,"  In-line step                             : %10d\n", area_ILXL[4]);
      if (area_ILXL[2] != IMISSING)
        LogKit::LogFormatted(LogKit::Low,"  Cross-line start                         : %10d\n", area_ILXL[2]);
      if (area_ILXL[3] != IMISSING)
        LogKit::LogFormatted(LogKit::Low,"  Cross-line end                           : %10d\n", area_ILXL[3]);
      if (area_ILXL[5] != IMISSING)
        LogKit::LogFormatted(LogKit::Low,"  Cross-line step                          : %10d\n", area_ILXL[5]);
    }
  }
  else if (area_specification == ModelSettings::AREA_FROM_UTM ||
           area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM) {
    LogKit::LogFormatted(LogKit::Low," given as UTM coordinates\n");
    const SegyGeometry * geometry = model_settings->getAreaParameters();
    LogKit::LogFormatted(LogKit::Low,"  Reference point x                        : %10.1f\n", geometry->GetX0());
    LogKit::LogFormatted(LogKit::Low,"  Reference point y                        : %10.1f\n", geometry->GetY0());
    LogKit::LogFormatted(LogKit::Low,"  Length x                                 : %10.1f\n", geometry->Getlx());
    LogKit::LogFormatted(LogKit::Low,"  Length y                                 : %10.1f\n", geometry->Getly());
    if (area_specification == ModelSettings::AREA_FROM_UTM) {
      LogKit::LogFormatted(LogKit::Low,"  Sample density x                         : %10.1f\n", geometry->GetDx());
      LogKit::LogFormatted(LogKit::Low,"  Sample density y                         : %10.1f\n", geometry->GetDy());
    }
    LogKit::LogFormatted(LogKit::Low,"  Rotation                                 : %10.4f\n", geometry->GetAngle()*(180.0/NRLib::Pi)*(-1));
    if (area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM) {
      LogKit::LogFormatted(LogKit::Low,"and snapped to seismic data\n");
      LogKit::LogFormatted(LogKit::Low,"  Grid                                     : "+grid_file+"\n");
    }
  }
  else if (area_specification == ModelSettings::AREA_FROM_SURFACE) {
    LogKit::LogFormatted(LogKit::Low," taken from surface\n");
    LogKit::LogFormatted(LogKit::Low,"  Reference surface                        : "+input_files->getAreaSurfaceFile()+"\n");
    if (area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE) {
      LogKit::LogFormatted(LogKit::Low," and snapped to seismic data\n");
      LogKit::LogFormatted(LogKit::Low,"  Grid                                     : "+grid_file+"\n");
    }
  }

  //
  // SURFACES
  //
  LogKit::LogFormatted(LogKit::Low,"\nTime surfaces:\n");
  if (model_settings->getParallelTimeSurfaces())
  {
    LogKit::LogFormatted(LogKit::Low,"  Reference surface                        : "+input_files->getTimeSurfFile(0)+"\n");
    LogKit::LogFormatted(LogKit::Low,"  Shift to top surface                     : %10.1f\n", model_settings->getTimeDTop());
    LogKit::LogFormatted(LogKit::Low,"  Time slice                               : %10.1f\n", model_settings->getTimeLz());
    LogKit::LogFormatted(LogKit::Low,"  Sampling density                         : %10.1f\n", model_settings->getTimeDz());
    LogKit::LogFormatted(LogKit::Low,"  Number of layers                         : %10d\n",   int(model_settings->getTimeLz()/model_settings->getTimeDz()+0.5));
  }
  else
  {
    const std::string & top_name  = input_files->getTimeSurfFile(0);
    const std::string & base_name = input_files->getTimeSurfFile(1);

    if (NRLib::IsNumber(top_name))
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : %10.2f\n",atof(top_name.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Top surface                              : "+top_name+"\n");

    if (NRLib::IsNumber(base_name))
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : %10.2f\n", atof(base_name.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Base surface                             : "+base_name+"\n");
      LogKit::LogFormatted(LogKit::Low,"  Number of layers                         : %10d\n", model_settings->getTimeNz());

    LogKit::LogFormatted(LogKit::Low,"  Minimum allowed value for lmin/lmax      : %10.2f\n", model_settings->getLzLimit());
  }
  if (input_files->getCorrDirFile() != "")
    LogKit::LogFormatted(LogKit::Low,"\n  Correlation direction                    : "+input_files->getCorrDirFile()+"\n");

  if (model_settings->getDoDepthConversion())
  {
    LogKit::LogFormatted(LogKit::Low,"\nDepth conversion:\n");
    if (input_files->getDepthSurfFile(0) != "")
      LogKit::LogFormatted(LogKit::Low,"  Top depth surface                        : "+input_files->getDepthSurfFile(0)+"\n");
    else
      LogKit::LogFormatted(LogKit::Low,"  Top depth surface                        : %s\n", "Made from base depth surface and velocity field");
    if (input_files->getDepthSurfFile(1) != "")
      LogKit::LogFormatted(LogKit::Low,"  Base depth surface                       : "+input_files->getDepthSurfFile(1)+"\n");
    else
      LogKit::LogFormatted(LogKit::Low,"  Base depth surface                       : %s\n", "Made from top depth surface and velocity field");
    std::string velocityField = input_files->getVelocityField();
    if (model_settings->getVelocityFromInversion()) {
      velocityField = "Use Vp from inversion";
    }
     LogKit::LogFormatted(LogKit::Low,"  Velocity field                           : "+velocityField+"\n");
  }

  const std::string & top_WEI  = input_files->getWaveletEstIntFileTop(0);
  const std::string & base_WEI = input_files->getWaveletEstIntFileBase(0);

  if (top_WEI != "" || base_WEI != "") {
    LogKit::LogFormatted(LogKit::Low,"\nWavelet estimation interval:\n");
    if (NRLib::IsNumber(top_WEI))
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : %10.2f\n",atof(top_WEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : "+top_WEI+"\n");

    if (NRLib::IsNumber(base_WEI))
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : %10.2f\n",atof(base_WEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : "+base_WEI+"\n");
  }

  const std::string & top_FEI  = input_files->getFaciesEstIntFile(0);
  const std::string & base_FEI = input_files->getFaciesEstIntFile(1);

  if (top_FEI != "" || base_FEI != "") {
    LogKit::LogFormatted(LogKit::Low,"\nFacies estimation interval:\n");
    if (NRLib::IsNumber(top_FEI))
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : %10.2f\n",atof(top_FEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Start time                               : "+top_FEI+"\n");

    if (NRLib::IsNumber(base_FEI))
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : %10.2f\n",atof(base_FEI.c_str()));
    else
      LogKit::LogFormatted(LogKit::Low,"  Stop time                                : "+base_FEI+"\n");
  }

  //
  // BACKGROUND
  //
  if (model_settings->getGenerateBackground())
  {
    LogKit::LogFormatted(LogKit::Low,"\nBackground model (estimated):\n");
    if (input_files->getBackVelFile() != "")
      LogKit::LogFormatted(LogKit::Low,"  Trend for p-wave velocity                : "+input_files->getBackVelFile()+"\n");
    Vario       * vario  = model_settings->getBackgroundVario();
    GenExpVario * pVario = dynamic_cast<GenExpVario*>(vario);
    LogKit::LogFormatted(LogKit::Low,"  Variogram\n");
    LogKit::LogFormatted(LogKit::Low,"    Model                                  : %10s\n",(vario->getType()).c_str());
    if (pVario != NULL)
    LogKit::LogFormatted(LogKit::Low,"    Power                                  : %10.1f\n",pVario->getPower());
    LogKit::LogFormatted(LogKit::Low,"    Range                                  : %10.1f\n",vario->getRange());
    if (vario->getAnisotropic())
    {
      LogKit::LogFormatted(LogKit::Low,"    Subrange                               : %10.1f\n",vario->getSubRange());
      LogKit::LogFormatted(LogKit::Low,"    Azimuth                                : %10.1f\n",90.0 - vario->getAngle()*(180/M_PI));
    }
    LogKit::LogFormatted(LogKit::Low,"  High cut frequency for well logs         : %10.1f\n",model_settings->getMaxHzBackground());
    if (model_settings->getMultizoneBackground() == true) {
      std::vector<std::string> surface_files = input_files->getMultizoneSurfaceFiles();
      std::vector<int> erosion               = model_settings->getErosionPriority();
      std::vector<double> uncertainty        = model_settings->getSurfaceUncertainty();
      std::vector<int> structure             = model_settings->getCorrelationStructure();
      int nZones = static_cast<int>(surface_files.size()-1);
      LogKit::LogFormatted(LogKit::Low,"\n  Multizone background model:\n");
      LogKit::LogFormatted(LogKit::Low,"    Top surface file                       : "+surface_files[0]+"\n");
      LogKit::LogFormatted(LogKit::Low,"    Top surface erosion priority           : %10d\n",erosion[0]);
      for (int i=0; i<nZones; i++) {
        LogKit::LogFormatted(LogKit::Low,"\n    Zone%2d\n",i+1);
        LogKit::LogFormatted(LogKit::Low,"      Base surface file                    : "+surface_files[i+1]+"\n");
        LogKit::LogFormatted(LogKit::Low,"      Base surface erosion priority        : %10d\n",erosion[i+1]);
        LogKit::LogFormatted(LogKit::Low,"      Base surface Beta uncertainty        : %10.1f\n",uncertainty[i+1]);
        if (structure[i+1] == ModelSettings::TOP)
          LogKit::LogFormatted(LogKit::Low,"      Correlation structure                :        Top\n");
        else if (structure[i+1] == ModelSettings::BASE)
          LogKit::LogFormatted(LogKit::Low,"      Correlation structure                :       Base\n");
        else if (structure[i+1] == ModelSettings::COMPACTION)
          LogKit::LogFormatted(LogKit::Low,"      Correlation structure                : Compaction\n");
      }
    }
  }
  else
  {
    if (model_settings->getForwardModeling()==true)
      LogKit::LogFormatted(LogKit::Low,"\nEarth model:\n");
    else
      LogKit::LogFormatted(LogKit::Low,"\nBackground model:\n");

    if (model_settings->getUseAIBackground()) {
      if (model_settings->getConstBackValue(0) > 0)
        LogKit::LogFormatted(LogKit::Low,"  Acoustic impedance                       : %10.1f\n",model_settings->getConstBackValue(0));
      else
        LogKit::LogFormatted(LogKit::Low,"  Acoustic impedance read from file        : "+input_files->getBackFile(0)+"\n");
    }
    else {
      if (model_settings->getConstBackValue(0) > 0)
        LogKit::LogFormatted(LogKit::Low,"  P-wave velocity                          : %10.1f\n",model_settings->getConstBackValue(0));
      else
        LogKit::LogFormatted(LogKit::Low,"  P-wave velocity read from file           : "+input_files->getBackFile(0)+"\n");
    }

    if (model_settings->getUseSIBackground()) {
      if (model_settings->getConstBackValue(1) > 0)
        LogKit::LogFormatted(LogKit::Low,"  Shear impedance                          : %10.1f\n",model_settings->getConstBackValue(1));
      else
        LogKit::LogFormatted(LogKit::Low,"  Shear impedance read from file           : "+input_files->getBackFile(1)+"\n");
    }
    else if (model_settings->getUseVpVsBackground()) {
      if (model_settings->getConstBackValue(1) > 0)
        LogKit::LogFormatted(LogKit::Low,"  Vp/Vs                                    : %10.1f\n",model_settings->getConstBackValue(1));
      else
        LogKit::LogFormatted(LogKit::Low,"  Vp/Vs  read from file                    : "+input_files->getBackFile(1)+"\n");
    }
    else {
      if (model_settings->getConstBackValue(1) > 0)
        LogKit::LogFormatted(LogKit::Low,"  S-wave velocity                          : %10.1f\n",model_settings->getConstBackValue(1));
      else
        LogKit::LogFormatted(LogKit::Low,"  S-wave velocity read from file           : "+input_files->getBackFile(1)+"\n");
    }

    if (model_settings->getConstBackValue(2) > 0)
      LogKit::LogFormatted(LogKit::Low,"  Density                                  : %10.1f\n",model_settings->getConstBackValue(2));
    else
      LogKit::LogFormatted(LogKit::Low,"  Density read from file                   : "+input_files->getBackFile(2)+"\n");
  }

  TraceHeaderFormat * thf_old = model_settings->getTraceHeaderFormat();
  if (thf_old != NULL)
  {
    LogKit::LogFormatted(LogKit::Low,"\nAdditional SegY trace header format:\n");
    if (thf_old != NULL) {
      LogKit::LogFormatted(LogKit::Low,"  Format name                              : "+thf_old->GetFormatName()+"\n");
      if (thf_old->GetBypassCoordScaling())
        LogKit::LogFormatted(LogKit::Low,"  Bypass coordinate scaling                :        yes\n");
      if (!thf_old->GetStandardType())
      {
        LogKit::LogFormatted(LogKit::Low,"  Start pos coordinate scaling             : %10d\n",thf_old->GetScalCoLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos trace x coordinate             : %10d\n",thf_old->GetUtmxLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos trace y coordinate             : %10d\n",thf_old->GetUtmyLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos inline index                   : %10d\n",thf_old->GetInlineLoc());
        LogKit::LogFormatted(LogKit::Low,"  Start pos crossline index                : %10d\n",thf_old->GetCrosslineLoc());
        LogKit::LogFormatted(LogKit::Low,"  Coordinate system                        : %10s\n",thf_old->GetCoordSys()==0 ? "UTM" : "ILXL" );
      }
    }
  }

  if (model_settings->getForwardModeling())
  {
    //
    // SEISMIC
    //
    LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for seismic:\n");
    LogKit::LogFormatted(LogKit::Low,"  Generating seismic                       : %10s\n","yes");
    std::vector<float> angle = model_settings->getAngle(0);
    for (int i = 0; i < model_settings->getNumberOfAngles(0); i++) //Forward modeling can only be done for one time lapse
    {
      LogKit::LogFormatted(LogKit::Low,"\nSettings for AVO stack %d:\n",i+1);
      LogKit::LogFormatted(LogKit::Low,"  Angle                                    : %10.1f\n",(angle[i]*180/M_PI));
      LogKit::LogFormatted(LogKit::Low,"  Read wavelet from file                   : "+input_files->getWaveletFile(0,i)+"\n");
    }
  }
  else
  {
    //
    // PRIOR CORRELATION
    //
    Vario * corr = model_settings->getLateralCorr();
    if (corr != NULL) {
      GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
      LogKit::LogFormatted(LogKit::Low,"\nPrior correlation (of residuals):\n");
      LogKit::LogFormatted(LogKit::Low,"  Range of allowed parameter values:\n");
      LogKit::LogFormatted(LogKit::Low,"    Var{Vp}  - min                         : %10.1e\n",model_settings->getVarAlphaMin());
      LogKit::LogFormatted(LogKit::Low,"    Var{Vp}  - max                         : %10.1e\n",model_settings->getVarAlphaMax());
      LogKit::LogFormatted(LogKit::Low,"    Var{Vs}  - min                         : %10.1e\n",model_settings->getVarBetaMin());
      LogKit::LogFormatted(LogKit::Low,"    Var{Vs}  - max                         : %10.1e\n",model_settings->getVarBetaMax());
      LogKit::LogFormatted(LogKit::Low,"    Var{Rho} - min                         : %10.1e\n",model_settings->getVarRhoMin());
      LogKit::LogFormatted(LogKit::Low,"    Var{Rho} - max                         : %10.1e\n",model_settings->getVarRhoMax());
      LogKit::LogFormatted(LogKit::Low,"  Lateral correlation:\n");
      LogKit::LogFormatted(LogKit::Low,"    Model                                  : %10s\n",(corr->getType()).c_str());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::Low,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::Low,"    Range                                  : %10.1f\n",corr->getRange());
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::Low,"    Subrange                               : %10.1f\n",corr->getSubRange());
        LogKit::LogFormatted(LogKit::Low,"    Azimuth                                : %10.1f\n",90.0 - corr->getAngle()*(180/M_PI));
      }
    }
    //
    // PRIOR FACIES
    //
    if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE ||
        model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
        // Can not be written when FACIES_FROM_WELLS as this information not is extracted yet
    {
      LogKit::LogFormatted(LogKit::Low,"\nPrior facies probabilities:\n");
      if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
      {
        typedef std::map<std::string,float> mapType;
        mapType myMap = model_settings->getPriorFaciesProb();

        for (mapType::iterator i=myMap.begin();i!=myMap.end();i++)
          LogKit::LogFormatted(LogKit::Low,"   %-12s                            : %10.2f\n",(i->first).c_str(),i->second);
      }
      else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
      {
        typedef std::map<std::string,std::string> mapType;
        mapType myMap = input_files->getPriorFaciesProbFile();

        for (mapType::iterator i=myMap.begin();i!=myMap.end();i++)
          LogKit::LogFormatted(LogKit::Low,"   %-12s                            : %10s\n",(i->first).c_str(),(i->second).c_str());
      }
    }

    if (model_settings->getFaciesProbFromRockPhysics()) {
      LogKit::LogFormatted(LogKit::Low,"\nRock physics:\n");
      {
        const std::map<std::string, DistributionsRockStorage *>& rock_storage = model_settings->getRockStorage();
        std::map<std::string, DistributionsRockStorage *>::const_iterator it;
        for (it = rock_storage.begin() ; it != rock_storage.end() ; it++)
          LogKit::LogFormatted(LogKit::Low,"  Rock                                     : %10s\n",(it->first).c_str());
      }
      {
        const std::map<std::string, DistributionsDryRockStorage *>& dryrock_storage = model_settings->getDryRockStorage();
        std::map<std::string, DistributionsDryRockStorage *>::const_iterator it;
        for (it = dryrock_storage.begin() ; it != dryrock_storage.end() ; it++)
          LogKit::LogFormatted(LogKit::Low,"  Dry Rock                                 : %10s\n",(it->first).c_str());
      }
      {
        const std::map<std::string, DistributionsSolidStorage *>& solid_storage = model_settings->getSolidStorage();
        std::map<std::string, DistributionsSolidStorage *>::const_iterator it;
        for (it = solid_storage.begin() ; it != solid_storage.end() ; it++)
          LogKit::LogFormatted(LogKit::Low,"  Solid                                    : %10s\n",(it->first).c_str());
      }
      {
        const std::map<std::string, DistributionsFluidStorage *>& fluid_storage = model_settings->getFluidStorage();
        std::map<std::string, DistributionsFluidStorage *>::const_iterator it;
        for (it = fluid_storage.begin() ; it != fluid_storage.end() ; it++)
          LogKit::LogFormatted(LogKit::Low,"  Fluid                                    : %10s\n",(it->first).c_str());
      }
      {
        const std::map<std::string, std::vector<DistributionWithTrendStorage *> >& reservoir_storage = model_settings->getReservoirVariable();
        std::map<std::string, std::vector<DistributionWithTrendStorage *> >::const_iterator it;
        for (it = reservoir_storage.begin() ; it != reservoir_storage.end() ; it++)
          LogKit::LogFormatted(LogKit::Low,"  Reservoir variable                       : %10s\n",(it->first).c_str());
      }
      {
         std::vector<std::string> trend_param = model_settings->getTrendCubeParameters();
         std::vector<int>         cube_type   = model_settings->getTrendCubeType();
         if (trend_param.size() > 0) {
           LogKit::LogFormatted(LogKit::Low,"  \nTrend cubes:\n");
           for (size_t i=0; i<trend_param.size(); i++) {
             if (cube_type[i] == ModelSettings::CUBE_FROM_FILE)
               LogKit::LogFormatted(LogKit::Low,"   From file                               : %10s\n", (trend_param[i]).c_str());
             else if (cube_type[i] == ModelSettings::STRATIGRAPHIC_DEPTH)
               LogKit::LogFormatted(LogKit::Low,"   Stratigraphic                           : %10s\n", (trend_param[i]).c_str());
             else if (cube_type[i] == ModelSettings::TWT)
               LogKit::LogFormatted(LogKit::Low,"   TWT                                     : %10s\n", (trend_param[i]).c_str());
           }
         }
      }
    }

    //
    // SEISMIC
    //
    if (model_settings->getNoSeismicNeeded()==false)
    {
      LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for seismic:\n");
      LogKit::LogFormatted(LogKit::Low,"  White noise component                    : %10.2f\n",model_settings->getWNC());
      LogKit::LogFormatted(LogKit::Low,"  Low cut for inversion                    : %10.1f\n",model_settings->getLowCut());
      LogKit::LogFormatted(LogKit::Low,"  High cut for inversion                   : %10.1f\n",model_settings->getHighCut());
      LogKit::LogFormatted(LogKit::Low,"  Guard zone outside interval of interest  : %10.1f ms\n",model_settings->getGuardZone());
      LogKit::LogFormatted(LogKit::Low,"  Smoothing length in guard zone           : %10.1f ms\n",model_settings->getSmoothLength());
      LogKit::LogFormatted(LogKit::Low,"  Interpolation threshold                  : %10.1f ms\n",model_settings->getEnergyThreshold());

      if (model_settings->getDo4DInversion())
        LogKit::LogFormatted(LogKit::Low,"\n4D seismic data:\n");

      for (int i=0; i<model_settings->getNumberOfTimeLapses(); i++){
        if (model_settings->getDo4DInversion())
          LogKit::LogFormatted(LogKit::Low,"\nVintage:\n");
        if (model_settings->getVintageMonth(i)==IMISSING && model_settings->getVintageYear(i) != IMISSING)
          LogKit::LogFormatted(LogKit::Low,"    %-2d                                     : %10d\n", i+1, model_settings->getVintageYear(i));
        else if (model_settings->getVintageDay(i)==IMISSING && model_settings->getVintageMonth(i) != IMISSING && model_settings->getVintageYear(i) != IMISSING)
          LogKit::LogFormatted(LogKit::Low,"    %-2d                                     : %10d %4d\n", i+1, model_settings->getVintageYear(i), model_settings->getVintageMonth(i));
        else if (model_settings->getVintageDay(i)!=IMISSING && model_settings->getVintageMonth(i)!=IMISSING && model_settings->getVintageYear(i) != IMISSING)
          LogKit::LogFormatted(LogKit::Low,"    %-2d                                     : %10d %4d %4d\n", i+1, model_settings->getVintageYear(i), model_settings->getVintageMonth(i), model_settings->getVintageDay(i));

      corr  = model_settings->getAngularCorr(i);
      GenExpVario * pCorr = dynamic_cast<GenExpVario*>(corr);
      LogKit::LogFormatted(LogKit::Low,"  Angular correlation:\n");
      LogKit::LogFormatted(LogKit::Low,"    Model                                  : %10s\n",(corr->getType()).c_str());
      if (pCorr != NULL)
        LogKit::LogFormatted(LogKit::Low,"    Power                                  : %10.1f\n",pCorr->getPower());
      LogKit::LogFormatted(LogKit::Low,"    Range                                  : %10.1f\n",corr->getRange()*180.0/M_PI);
      if (corr->getAnisotropic())
      {
        LogKit::LogFormatted(LogKit::Low,"    Subrange                               : %10.1f\n",corr->getSubRange()*180.0/M_PI);
        LogKit::LogFormatted(LogKit::Low,"    Angle                                  : %10.1f\n",corr->getAngle());
      }
      bool estimate_noise = false;
      for (int j = 0; j < model_settings->getNumberOfAngles(i); j++) {
        estimate_noise = estimate_noise || model_settings->getEstimateSNRatio(i,j);
      }
      LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for wavelet:\n");
      if (estimate_noise)
        LogKit::LogFormatted(LogKit::Low,"  Maximum shift in noise estimation        : %10.1f\n",model_settings->getMaxWaveletShift());
      LogKit::LogFormatted(LogKit::High,  "  Minimum relative amplitude               : %10.3f\n",model_settings->getMinRelWaveletAmp());
      LogKit::LogFormatted(LogKit::High,  "  Wavelet tapering length                  : %10.1f\n",model_settings->getWaveletTaperingL());
      LogKit::LogFormatted(LogKit::High, "  Tuning factor for 3D wavelet estimation  : %10.1f\n", model_settings->getWavelet3DTuningFactor());
      LogKit::LogFormatted(LogKit::High, "  Smoothing range for gradient (3D wavelet): %10.1f\n", model_settings->getGradientSmoothingRange());
      LogKit::LogFormatted(LogKit::High, "  Estimate well gradient for seismic data  : %10s\n", (model_settings->getEstimateWellGradientFromSeismic() ? "yes" : "no"));

      if (model_settings->getOptimizeWellLocation()) {
        LogKit::LogFormatted(LogKit::Low,"\nGeneral settings for well locations:\n");
        LogKit::LogFormatted(LogKit::Low,"  Maximum offset                           : %10.1f\n",model_settings->getMaxWellOffset());
        LogKit::LogFormatted(LogKit::Low,"  Maximum vertical shift                   : %10.1f\n",model_settings->getMaxWellShift());
      }
        std::vector<float> angle            = model_settings->getAngle(i);
        std::vector<float> SN_ratio         = model_settings->getSNRatio(i);
        std::vector<bool>  estimate_wavelet = model_settings->getEstimateWavelet(i);
        std::vector<bool>  match_energies   = model_settings->getMatchEnergies(i);




        for (int j = 0; j < model_settings->getNumberOfAngles(i); j++)
        {
          LogKit::LogFormatted(LogKit::Low,"\nSettings for AVO stack %d:\n",j+1);
          LogKit::LogFormatted(LogKit::Low,"  Angle                                    : %10.1f\n",(angle[j]*180/M_PI));
          LogKit::LogFormatted(LogKit::Low,"  SegY start time                          : %10.1f\n",model_settings->getSegyOffset(i));
          TraceHeaderFormat * thf = model_settings->getTraceHeaderFormat(i,j);
          if (thf != NULL)
          {
            LogKit::LogFormatted(LogKit::Low,"  SegY trace header format:\n");
            LogKit::LogFormatted(LogKit::Low,"    Format name                            : "+thf->GetFormatName()+"\n");
            if (thf->GetBypassCoordScaling())
              LogKit::LogFormatted(LogKit::Low,"    Bypass coordinate scaling              :        yes\n");
            if (!thf->GetStandardType())
            {
              LogKit::LogFormatted(LogKit::Low,"    Start pos coordinate scaling           : %10d\n",thf->GetScalCoLoc());
              LogKit::LogFormatted(LogKit::Low,"    Start pos trace x coordinate           : %10d\n",thf->GetUtmxLoc());
              LogKit::LogFormatted(LogKit::Low,"    Start pos trace y coordinate           : %10d\n",thf->GetUtmyLoc());
              LogKit::LogFormatted(LogKit::Low,"    Start pos inline index                 : %10d\n",thf->GetInlineLoc());
              LogKit::LogFormatted(LogKit::Low,"    Start pos crossline index              : %10d\n",thf->GetCrosslineLoc());
              LogKit::LogFormatted(LogKit::Low,"    Coordinate system                      : %10s\n",thf->GetCoordSys()==0 ? "UTM" : "ILXL" );
            }
          }
          LogKit::LogFormatted(LogKit::Low,"  Data                                     : "+input_files->getSeismicFile(i,j)+"\n");
          if (estimate_wavelet[j])
            LogKit::LogFormatted(LogKit::Low,"  Estimate wavelet                         : %10s\n", "yes");
          else
            LogKit::LogFormatted(LogKit::Low,"  Read wavelet from file                   : "+input_files->getWaveletFile(i,j)+"\n");
          if (model_settings->getEstimateLocalShift(i,j))
           LogKit::LogFormatted(LogKit::Low,"  Estimate local shift map                 : %10s\n", "yes");
          else if (input_files->getShiftFile(i,j) != "")
            LogKit::LogFormatted(LogKit::Low,"  Local shift map                          : "+input_files->getShiftFile(i,j)+"\n");
          if (model_settings->getEstimateLocalScale(i,j))
            LogKit::LogFormatted(LogKit::Low,"  Estimate local scale map                 : %10s\n", "yes");
          else if (input_files->getScaleFile(i,j) != "")
            LogKit::LogFormatted(LogKit::Low,"  Local scale map                          : "+input_files->getScaleFile(i,j)+"\n");
          if (match_energies[j])
            LogKit::LogFormatted(LogKit::Low,"  Match empirical and theoretical energies : %10s\n", "yes");
          if (!estimate_wavelet[j] && !match_energies[j]){
            if (model_settings->getEstimateGlobalWaveletScale(i,j))
              LogKit::LogFormatted(LogKit::Low,"  Estimate global wavelet scale            : %10s\n","yes");
            else
              LogKit::LogFormatted(LogKit::Low,"  Global wavelet scale                     : %10.2f\n",model_settings->getWaveletScale(i,j));
          }
          if (model_settings->getEstimateSNRatio(i,j))
            LogKit::LogFormatted(LogKit::Low,"  Estimate signal-to-noise ratio           : %10s\n", "yes");
          else
            LogKit::LogFormatted(LogKit::Low,"  Signal-to-noise ratio                    : %10.1f\n",SN_ratio[j]);
          if (model_settings->getEstimateLocalNoise(i,j)) {
            if (input_files->getLocalNoiseFile(i,j) == "")
              LogKit::LogFormatted(LogKit::Low,"  Estimate local signal-to-noise ratio map : %10s\n", "yes");
            else
              LogKit::LogFormatted(LogKit::Low,"  Local signal-to-noise ratio map          : "+input_files->getLocalNoiseFile(i,j)+"\n");
          }
          if (model_settings->getEstimateLocalNoise(i,j))
            LogKit::LogFormatted(LogKit::Low,"  Estimate local noise                     : %10s\n", "yes");
          if (input_files->getLocalNoiseFile(i,j) != "")
            LogKit::LogFormatted(LogKit::Low,"  Local noise                              : "+input_files->getLocalNoiseFile(i,j)+"\n");
        }
      }
    }
  }
}

//void
//ModelGeneral::getCorrGradIJ(float & corrGradI, float &corrGradJ) const
//{
//  double angle  = simbox_->getAngle();
//  double cosrot = cos(angle);
//  double sinrot = sin(angle);
//  double dx     = simbox_->getdx();
//  double dy     = simbox_->getdy();
//
//  double cI =  dx*cosrot*gradX_ + dy*sinrot*gradY_;
//  double cJ = -dx*sinrot*gradX_ + dy*cosrot*gradY_;
//
//  corrGradI = float(cI/simbox_->getdz());
//  corrGradJ = float(cJ/simbox_->getdz());
//}

//void
//ModelGeneral::processDepthConversion(Simbox            * timeCutSimbox,
//                                     Simbox            * timeSimbox,
//                                     ModelSettings     * model_settings,
//                                     const InputFiles  * input_files,
//                                     std::string       & errText,
//                                     bool              & failed)
//{
//  FFTGrid * velocity = NULL;
//  if (timeCutSimbox != NULL)
//    loadVelocity(velocity, timeCutSimbox, timeCutSimbox, model_settings,
//                 input_files->getVelocityField(), velocityFromInversion_,
//                 errText, failed);
//  else
//    loadVelocity(velocity, timeSimbox, timeCutSimbox, model_settings,
//                 input_files->getVelocityField(), velocityFromInversion_,
//                 errText, failed);
//
//  if (!failed)
//  {
//    timeDepthMapping_ = new GridMapping();
//    timeDepthMapping_->setDepthSurfaces(input_files->getDepthSurfFiles(), failed, errText);
//    if (velocity != NULL)
//    {
//      velocity->setAccessMode(FFTGrid::RANDOMACCESS);
//      timeDepthMapping_->calculateSurfaceFromVelocity(velocity, timeSimbox);
//      timeDepthMapping_->setDepthSimbox(timeSimbox, timeSimbox->getnz(),
//                                        model_settings->getOutputGridFormat(),
//                                        failed, errText);            // NBNB-PAL: Er dettet riktig nz (timeCut vs time)?
//      timeDepthMapping_->makeTimeDepthMapping(velocity, timeSimbox);
//      velocity->endAccess();
//
//      if ((model_settings->getOutputGridsOther() & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
//        std::string baseName  = IO::FileTimeToDepthVelocity();
//        std::string sgriLabel = std::string("Time-to-depth velocity");
//        float       offset    = model_settings->getSegyOffset(0);//Only allow one segy offset for time lapse data
//        velocity->writeFile(baseName,
//                            IO::PathToVelocity(),
//                            timeSimbox,
//                            sgriLabel,
//                            offset,
//                            timeDepthMapping_,
//                            timeCutMapping_);
//      }
//    }
//    else if (velocity==NULL && velocityFromInversion_==false)
//    {
//      timeDepthMapping_->setDepthSimbox(timeSimbox,
//                                        timeSimbox->getnz(),
//                                        model_settings->getOutputGridFormat(),
//                                        failed,
//                                        errText);
//
//    }
//  }
//  if (velocity != NULL)
//    delete velocity;
//}

//void ModelGeneral::processRockPhysics(Simbox                       * timeSimbox,
//                                      Simbox                       * timeCutSimbox,
//                                      ModelSettings                * model_settings,
//                                      bool                         & failed,
//                                      std::string                  & errTxt,
//                                      const InputFiles             * input_files)
//{
//  if (model_settings->getFaciesProbFromRockPhysics()){
//
//    LogKit::WriteHeader("Processing Rock Physics");
//
//    trend_cubes_ = CravaTrend(timeSimbox,
//                              timeCutSimbox,
//                              model_settings,
//                              failed,
//                              errTxt,
//                              input_files);
//
//    int n_vintages = model_settings->getNumberOfVintages();
//
//    const std::string&                       path                               = input_files->getInputDirectory();
//    const std::vector<std::string>&          trend_cube_parameters              = model_settings->getTrendCubeParameters();
//    const std::vector<std::vector<double> >& trend_cube_sampling                = trend_cubes_.GetTrendCubeSampling();
//
//    const std::map<std::string, std::vector<DistributionWithTrendStorage *> >& reservoir_variable = model_settings->getReservoirVariable();
//    for (std::map<std::string, std::vector<DistributionWithTrendStorage *> >::const_iterator it = reservoir_variable.begin(); it != reservoir_variable.end(); it++) {
//
//      std::vector<DistributionWithTrendStorage *> storage = it->second;
//      std::vector<DistributionWithTrend *> dist_vector(storage.size());
//
//      for (size_t i=0; i<storage.size(); i++) {
//        dist_vector[i]                = storage[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
//      }
//
//      reservoir_variables_[it->first] = dist_vector;
//    }
//
//    if (errTxt == "") {
//
//      float alpha_min     = model_settings->getAlphaMin();
//      float alpha_max     = model_settings->getAlphaMax();
//      float beta_min      = model_settings->getBetaMin();
//      float beta_max      = model_settings->getBetaMax();
//      float rho_min       = model_settings->getRhoMin();
//      float rho_max       = model_settings->getRhoMax();
//      float var_alpha_min = model_settings->getVarAlphaMin();
//      float var_alpha_max = model_settings->getVarAlphaMax();
//      float var_beta_min  = model_settings->getVarBetaMin();
//      float var_beta_max  = model_settings->getVarBetaMax();
//      float var_rho_min   = model_settings->getVarRhoMin();
//      float var_rho_max   = model_settings->getVarRhoMax();
//
//      const std::map<std::string, DistributionsFluidStorage   *>& fluid_storage    = model_settings->getFluidStorage();
//      const std::map<std::string, DistributionsSolidStorage   *>& solid_storage    = model_settings->getSolidStorage();
//      const std::map<std::string, DistributionsDryRockStorage *>& dry_rock_storage = model_settings->getDryRockStorage();
//      const std::map<std::string, DistributionsRockStorage    *>& rock_storage     = model_settings->getRockStorage();
//
//      // Map out reservoir variables for using in rocks to access resampling trigger.
//      std::vector<std::vector<DistributionWithTrend *> > res_var_vintage(0);
//      if (reservoir_variables_.size() > 0) {
//        size_t nVintages = reservoir_variables_.begin()->second.size();
//        res_var_vintage.resize(nVintages);
//        for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator var_it = reservoir_variables_.begin();
//          var_it != reservoir_variables_.end();var_it++)
//        {
//          for (size_t vin_index=0;vin_index < var_it->second.size();vin_index++)
//            res_var_vintage[vin_index].push_back((var_it->second)[vin_index]);
//        }
//      }
//
//
//      std::map<std::string, float> facies_probabilities = model_settings->getPriorFaciesProb();
//      std::map<std::string, std::string> facies_cubes   = input_files->getPriorFaciesProbFile();
//      std::vector<std::string> all_facies_names         = facies_names_;
//
//      for (std::map<std::string, float>::iterator it_prob = facies_probabilities.begin(); it_prob != facies_probabilities.end(); it_prob++)
//        all_facies_names.push_back(it_prob->first);
//      for (std::map<std::string, std::string>::iterator it_cube = facies_cubes.begin(); it_cube != facies_cubes.end(); it_cube++)
//        all_facies_names.push_back(it_cube->first);
//
//      std::sort(all_facies_names.begin(), all_facies_names.end());
//
//      std::string prev_facies = "";
//      for (size_t it=0; it<all_facies_names.size(); it++) {
//        if (all_facies_names[it] != prev_facies) {
//          prev_facies = all_facies_names[it];
//
//          std::map<std::string, DistributionsRockStorage *>::const_iterator iter = rock_storage.find(all_facies_names[it]);
//          if (iter != rock_storage.end()) {
//
//            std::string rockErrTxt = "";
//
//            std::string name = iter->first;
//            LogKit::LogFormatted(LogKit::Low, "\nRock \'"+name+"\':\n");
//
//            DistributionsRockStorage * storage    = iter   ->second;
//            std::vector<DistributionsRock *> rock = storage->GenerateDistributionsRock(n_vintages,
//                                                                                       path,
//                                                                                       trend_cube_parameters,
//                                                                                       trend_cube_sampling,
//                                                                                       rock_storage,
//                                                                                       solid_storage,
//                                                                                       dry_rock_storage,
//                                                                                       fluid_storage,
//                                                                                       rockErrTxt);
//
//            if (rockErrTxt == "") {
//
//              int n_vintages = static_cast<int>(rock.size());
//              if (n_vintages > 1)
//                LogKit::LogFormatted(LogKit::Low, "Number of vintages: %4d\n", n_vintages);
//
//              for (int i=0; i<n_vintages; i++) {
//                if (n_vintages > 1)
//                  LogKit::LogFormatted(LogKit::Low, "\nVintage number: %4d\n", i+1);
//
//                //Completing the top level rocks, by setting access to reservoir variables and sampling distribution.
//                rock[i]->CompleteTopLevelObject(res_var_vintage[i]);
//
//                std::vector<bool> has_trends = rock[i]->HasTrend();
//                bool              has_trend = false;
//                for (size_t j=0; j<has_trends.size(); j++) {
//                  if (has_trends[j] == true) {
//                    has_trend = true;
//                    break;
//                  }
//                }
//
//                std::vector<double> expectation  = rock[i]->GetMeanLogExpectation();
//                NRLib::Grid2D<double> covariance = rock[i]->GetMeanLogCovariance();
//
//                printExpectationAndCovariance(expectation, covariance, has_trend);
//
//                std::string tmpErrTxt = "";
//                if (std::exp(expectation[0]) < alpha_min  || std::exp(expectation[0]) > alpha_max) {
//                  tmpErrTxt += "Vp value of "+NRLib::ToString(std::exp(expectation[0]))+" detected: ";
//                  tmpErrTxt += "Vp should be in the interval ("+NRLib::ToString(alpha_min)+", "+NRLib::ToString(alpha_max)+") m/s\n";
//                }
//                if (std::exp(expectation[1]) < beta_min  || std::exp(expectation[1]) > beta_max) {
//                  if (typeid(*(storage)) == typeid(ReussRockStorage))
//                    tmpErrTxt += "Vs value of 0 detected. Note that the Reuss model gives Vs=0; hence it can not be used to model a facies\n";
//                  else
//                    tmpErrTxt += "Vs value of "+NRLib::ToString(std::exp(expectation[1]))+" detected: ";
//                  tmpErrTxt += "Vs should be in the interval ("+NRLib::ToString(beta_min)+", "+NRLib::ToString(beta_max)+") m/s\n";
//                }
//                if (std::exp(expectation[2]) < rho_min  || std::exp(expectation[2]) > rho_max) {
//                  tmpErrTxt += "Rho value of "+NRLib::ToString(std::exp(expectation[2]))+" detected: ";
//                  tmpErrTxt += "Rho should be in the interval ("+NRLib::ToString(rho_min)+", "+NRLib::ToString(rho_max)+") g/cm^3\n";
//                }
//
//                if (tmpErrTxt != "") {
//                  errTxt += "\nToo high or low seismic properties calculated for rock '"+iter->first+"':\n";
//                  errTxt += tmpErrTxt;
//                }
//
//                std::string varErrTxt = "";
//                if (covariance(0,0) < var_alpha_min  || covariance(0,0) > var_alpha_max) {
//                  varErrTxt += "Var(log Vp) value of "+NRLib::ToString(covariance(0,0))+" detected: ";
//                  varErrTxt += "Var(log Vp) should be in the interval ("+NRLib::ToString(var_alpha_min)+", "+NRLib::ToString(var_alpha_max)+")\n";
//                }
//                if (covariance(1,1) < var_beta_min  || covariance(1,1) > var_beta_max) {
//                  varErrTxt += "Var(log Vs) value of "+NRLib::ToString(covariance(1,1))+" detected: ";
//                  varErrTxt += "Var(log Vs) should be in the interval ("+NRLib::ToString(var_beta_min)+", "+NRLib::ToString(var_beta_max)+")\n";
//                }
//                if (covariance(2,2) < var_rho_min  || covariance(2,2) > var_rho_max) {
//                  varErrTxt += "Var(log Rho) value of "+NRLib::ToString(covariance(2,2))+" detected: ";
//                  varErrTxt += "Var(log Rho) should be in the interval ("+NRLib::ToString(var_rho_min)+", "+NRLib::ToString(var_rho_max)+")\n";
//                }
//
//                if (varErrTxt != "") {
//                  errTxt += "\nToo high or low variance of seismic properties calculated for rock '"+iter->first+"':\n";
//                  errTxt += varErrTxt;
//                }
//
//                // Check correlations
//                float corr01 = static_cast<float>(covariance(0,1)/(sqrt(covariance(0,0)*covariance(1,1))));
//                float corr02 = static_cast<float>(covariance(0,2)/(sqrt(covariance(0,0)*covariance(2,2))));
//                float corr12 = static_cast<float>(covariance(1,2)/(sqrt(covariance(1,1)*covariance(2,2))));
//
//                NRLib::SymmetricMatrix corr(3);
//                corr(0,0) = 1;
//                corr(1,1) = 1;
//                corr(2,2) = 1;
//                corr(0,1) = corr01;
//                corr(0,2) = corr02;
//                corr(1,2) = corr12;
//
//                try {
//                  NRLib::CholeskyInvert(corr);
//                }
//                catch (NRLib::Exception & e) {
//                  errTxt += e.what();
//                  errTxt += " for rock '"+iter->first+"':\n";
//                  errTxt += "  The variabels in the rock model are probably linearly dependent\n";
//                }
//
//                if (varErrTxt != "" || tmpErrTxt != "")
//                  break;
//
//              }
//
//              rock_distributions_[all_facies_names[it]] = rock;
//            }
//            else
//              errTxt += rockErrTxt;
//          }
//          else
//            errTxt += "The facies "+all_facies_names[it]+" is not one of the rocks in the rock physics model\n";
//          //rock_distributions_[it->first] = rock;
//        }
//      }
//    }
//
//    if (errTxt != "")
//      failed = true;
//  }
//}

//void ModelGeneral::printExpectationAndCovariance(const std::vector<double>   & expectation,
//                                                 const NRLib::Grid2D<double> & covariance,
//                                                 const bool                  & has_trend) const
//{
//  if (has_trend == true)
//      LogKit::LogFormatted(LogKit::Low,"\nMean expectation and covariance estimated over all trend values:\n");
//  else
//    LogKit::LogFormatted(LogKit::Low,"\nEstimated expectation and covariance:\n");
//  LogKit::LogFormatted(LogKit::Low,"\n");
//  LogKit::LogFormatted(LogKit::Low,"Expectation            Vp        Vs       Rho\n");
//  LogKit::LogFormatted(LogKit::Low,"----------------------------------------------\n");
//  LogKit::LogFormatted(LogKit::Low,"                  %7.2f   %7.2f   %7.3f \n",
//                       std::exp(expectation[0]), std::exp(expectation[1]), std::exp(expectation[2]));
//
//  LogKit::LogFormatted(LogKit::Low,"\n");
//  LogKit::LogFormatted(LogKit::Low,"Variances           ln Vp     ln Vs    ln Rho\n");
//  LogKit::LogFormatted(LogKit::Low,"----------------------------------------------\n");
//  LogKit::LogFormatted(LogKit::Low,"                  %.1e   %.1e   %.1e\n", covariance(0,0), covariance(1,1), covariance(2,2));
//
//  float corr01 = static_cast<float>(covariance(0,1)/(sqrt(covariance(0,0)*covariance(1,1))));
//  float corr02 = static_cast<float>(covariance(0,2)/(sqrt(covariance(0,0)*covariance(2,2))));
//  float corr12 = static_cast<float>(covariance(1,2)/(sqrt(covariance(1,1)*covariance(2,2))));
//
//  LogKit::LogFormatted(LogKit::Low,"\n");
//  LogKit::LogFormatted(LogKit::Low,"Corr   | ln Vp     ln Vs    ln Rho \n");
//  LogKit::LogFormatted(LogKit::Low,"-------+---------------------------\n");
//  LogKit::LogFormatted(LogKit::Low,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
//  LogKit::LogFormatted(LogKit::Low,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
//  LogKit::LogFormatted(LogKit::Low,"ln Rho |                     %5.2f \n",1.0f);
//  LogKit::LogFormatted(LogKit::Low,"\n");
//}

//void
//ModelGeneral::loadVelocity(FFTGrid           *& velocity,
//                           const Simbox       * timeSimbox,
//                           const Simbox       * timeCutSimbox,
//                           const ModelSettings * model_settings,
//                           const std::string  & velocityField,
//                           bool               & velocityFromInversion,
//                           std::string        & errText,
//                           bool               & failed)
//{
//  LogKit::WriteHeader("Setup time-to-depth relationship");
//
//  if (model_settings->getVelocityFromInversion() == true)
//  {
//    velocityFromInversion = true;
//    velocity = NULL;
//  }
//  else if (velocityField == "")
//    velocity = NULL;
//  else
//  {
//    const SegyGeometry      * dummy1 = NULL;
//    const TraceHeaderFormat * dummy2 = NULL;
//    const float               offset = model_settings->getSegyOffset(0); //Segy offset needs to be the same for all time lapse data
//    std::string errorText("");
//    readGridFromFile(velocityField,
//                     "velocity field",
//                     offset,
//                     velocity,
//                     dummy1,
//                     dummy2,
//                     FFTGrid::PARAMETER,
//                     timeSimbox,
//                     timeCutSimbox,
//                     model_settings,
//                     errorText);
//
//    if (errorText == "") { // No errors
//      //
//      // Check that the velocity grid is veldefined.
//      //
//      float logMin = model_settings->getAlphaMin();
//      float logMax = model_settings->getAlphaMax();
//      const int nzp = velocity->getNzp();
//      const int nyp = velocity->getNyp();
//      const int nxp = velocity->getNxp();
//      const int nz = velocity->getNz();
//      const int ny = velocity->getNy();
//      const int nx = velocity->getNx();
//      int tooLow  = 0;
//      int tooHigh = 0;
//      velocity->setAccessMode(FFTGrid::READ);
//      int rnxp = 2*(nxp/2 + 1);
//      for (int k = 0; k < nzp; k++)
//        for (int j = 0; j < nyp; j++)
//          for (int i = 0; i < rnxp; i++) {
//            if (i < nx && j < ny && k < nz) {
//              float value = velocity->getNextReal();
//              if (value < logMin && value != RMISSING) {
//                tooLow++;
//              }
//              if (value > logMax && value != RMISSING)
//                tooHigh++;
//            }
//          }
//      velocity->endAccess();
//
//      if (tooLow+tooHigh > 0) {
//        std::string text;
//        text += "\nThe velocity grid used as trend in the background model of Vp";
//        text += "\ncontains too small and/or too high velocities:";
//        text += "\n  Minimum Vp = "+NRLib::ToString(logMin,2)+"    Number of too low values  : "+NRLib::ToString(tooLow);
//        text += "\n  Maximum Vp = "+NRLib::ToString(logMax,2)+"    Number of too high values : "+NRLib::ToString(tooHigh);
//        text += "\nThe range of allowed values can changed using the ALLOWED_PARAMETER_VALUES keyword\n";
//        text += "\naborting...\n";
//        errText += "Reading of file '"+velocityField+"' for background velocity field failed.\n";
//        errText += text;
//        failed = true;
//      }
//    }
//    else {
//      errorText += "Reading of file \'"+velocityField+"\' for background velocity field failed.\n";
//      errText += errorText;
//      failed = true;
//    }
//  }
//}

void
ModelGeneral::writeAreas(const SegyGeometry * areaParams,
                         Simbox             * timeSimbox,
                         std::string        & text)
{
  double areaX0   = areaParams->GetX0();
  double areaY0   = areaParams->GetY0();
  double areaLx   = areaParams->Getlx();
  double areaLy   = areaParams->Getly();
  double areaDx   = areaParams->GetDx();
  double areaDy   = areaParams->GetDy();
  double areaRot  = areaParams->GetAngle();
  double areaXmin = RMISSING;
  double areaXmax = RMISSING;
  double areaYmin = RMISSING;
  double areaYmax = RMISSING;

  findSmallestSurfaceGeometry(areaX0, areaY0, areaLx, areaLy, areaRot,
                              areaXmin, areaYmin, areaXmax, areaYmax);

  LogKit::LogFormatted(LogKit::Low,"\nThe top and/or base time surfaces do not cover the area specified by the "+text);
  LogKit::LogFormatted(LogKit::Low,"\nPlease extrapolate surfaces or specify a smaller AREA in the model file.\n");
  LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx        ly     azimuth          dx      dy\n");
  LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
  double azimuth = (-1)*areaRot*(180.0/M_PI);
  if (azimuth < 0)
    azimuth += 360.0;
  LogKit::LogFormatted(LogKit::Low,"Model area       %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n\n",
                       areaX0, areaY0, areaLx, areaLy, azimuth, areaDx, areaDy);

  LogKit::LogFormatted(LogKit::Low,"Area                    xmin         xmax           ymin        ymax\n");
  LogKit::LogFormatted(LogKit::Low,"--------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %11.2f %11.2f\n",
                       text.c_str(),areaXmin, areaXmax, areaYmin, areaYmax);
  const NRLib::Surface<double> & top  = timeSimbox->GetTopSurface();
  const NRLib::Surface<double> & base = timeSimbox->GetBotSurface();
  LogKit::LogFormatted(LogKit::Low,"Top surface      %11.2f  %11.2f    %11.2f %11.2f\n",
                       top.GetXMin(), top.GetXMax(), top.GetYMin(), top.GetYMax());
  LogKit::LogFormatted(LogKit::Low,"Base surface     %11.2f  %11.2f    %11.2f %11.2f\n",
                       base.GetXMin(), base.GetXMax(), base.GetYMin(), base.GetYMax());
}

void
ModelGeneral::findSmallestSurfaceGeometry(const double   x0,
                                          const double   y0,
                                          const double   lx,
                                          const double   ly,
                                          const double   rot,
                                          double       & xMin,
                                          double       & yMin,
                                          double       & xMax,
                                          double       & yMax)
{
  xMin = x0 - ly*sin(rot);
  xMax = x0 + lx*cos(rot);
  yMin = y0;
  yMax = y0 + lx*sin(rot) + ly*cos(rot);
  if (rot < 0) {
    xMin = x0;
    xMax = x0 + lx*cos(rot) - ly*sin(rot);
    yMin = y0 + lx*sin(rot);
    yMax = y0 + ly*cos(rot);
  }
}

void
ModelGeneral::getGeometryFromGridOnFile(const std::string         gridFile,
                                        const TraceHeaderFormat * thf,
                                        SegyGeometry           *& geometry,
                                        std::string             & errText)
{
  geometry = NULL;

  if (gridFile != "") { //May change the condition here, but need geometry if we want to set XL/IL
    int fileType = IO::findGridType(gridFile);
    if (fileType == IO::CRAVA) {
      geometry = geometryFromCravaFile(gridFile);
    }
    else if (fileType == IO::SEGY) {
      try
      {
        geometry = SegY::FindGridGeometry(gridFile, thf);
      }
      catch (NRLib::Exception & e)
      {
        errText = e.what();
      }
    }
    else if (fileType == IO::STORM)
      geometry = geometryFromStormFile(gridFile, errText);
    else if (fileType==IO::SGRI) {
      bool scale = true;
      geometry = geometryFromStormFile(gridFile, errText, scale);
    }
    else {
      errText = "Trying to read grid dimensions from unknown file format.\n";
    }
  }
  else {
    errText = "Cannot get geometry from file. The file name is empty.\n";
  }
}

SegyGeometry *
ModelGeneral::geometryFromCravaFile(const std::string & fileName)
{
  std::ifstream binFile;
  NRLib::OpenRead(binFile, fileName, std::ios::in | std::ios::binary);

  std::string fileType;
  getline(binFile,fileType);

  double x0      = NRLib::ReadBinaryDouble(binFile);
  double y0      = NRLib::ReadBinaryDouble(binFile);
  double dx      = NRLib::ReadBinaryDouble(binFile);
  double dy      = NRLib::ReadBinaryDouble(binFile);
  int    nx      = NRLib::ReadBinaryInt(binFile);
  int    ny      = NRLib::ReadBinaryInt(binFile);
  double IL0     = NRLib::ReadBinaryDouble(binFile);
  double XL0     = NRLib::ReadBinaryDouble(binFile);
  double ilStepX = NRLib::ReadBinaryDouble(binFile);
  double ilStepY = NRLib::ReadBinaryDouble(binFile);
  double xlStepX = NRLib::ReadBinaryDouble(binFile);
  double xlStepY = NRLib::ReadBinaryDouble(binFile);
  double rot     = NRLib::ReadBinaryDouble(binFile);

  binFile.close();

  SegyGeometry * geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                             IL0, XL0, ilStepX, ilStepY,
                                             xlStepX, xlStepY, rot);
  return(geometry);
}

SegyGeometry *
ModelGeneral::geometryFromStormFile(const std::string & fileName,
                                    std::string       & errText,
                                    bool scale)
{
  SegyGeometry  * geometry  = NULL;
  StormContGrid * stormgrid = NULL;
  std::string     tmpErrText;
  float scalehor;
  if (scale==false)
  {
    scalehor = 1.0;
  }
  else //from sgri file
  {
    LogKit::LogFormatted(LogKit::Low,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
    scalehor  = 1000.0;
  }
  try
  {
    stormgrid = new StormContGrid(0,0,0);
    stormgrid->ReadFromFile(fileName);
    stormgrid->SetMissingCode(RMISSING);
  }
  catch (NRLib::Exception & e)
  {
    tmpErrText = e.what();
  }

  if (tmpErrText == "") {
    double x0      = stormgrid->GetXMin()*scalehor;
    double y0      = stormgrid->GetYMin()*scalehor;
    double dx      = stormgrid->GetDX()*scalehor;
    double dy      = stormgrid->GetDY()*scalehor;
    int    nx      = static_cast<int>(stormgrid->GetNI());
    int    ny      = static_cast<int>(stormgrid->GetNJ());
    double rot     = stormgrid->GetAngle();
    double IL0     = 0.0;  ///< Dummy value since information is not contained in format
    double XL0     = 0.0;  ///< Dummy value since information is not contained in format
    double ilStepX =   1;  ///< Dummy value since information is not contained in format
    double ilStepY =   1;  ///< Dummy value since information is not contained in format
    double xlStepX =   1;  ///< Dummy value since information is not contained in format
    double xlStepY =   1;  ///< Dummy value since information is not contained in format
    geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                IL0, XL0, ilStepX, ilStepY,
                                xlStepX, xlStepY, rot);
  }
  else {
    errText += tmpErrText;
  }

  if (stormgrid != NULL)
    delete stormgrid;

  return(geometry);
}

std::map<std::string, DistributionsRock *>
ModelGeneral::getRockDistributionTime0() const
{
  std::map<std::string, DistributionsRock *> rock_dist_t0;

  for (std::map<std::string, std::vector<DistributionsRock *> >::const_iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
    std::string name = it->first;
    std::vector<DistributionsRock *> rock_dist = it->second;
    rock_dist_t0[name] = rock_dist[0];
  }

  return(rock_dist_t0);
}

FFTGrid*
ModelGeneral::createFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid)
{
  FFTGrid* fftGrid;

  if (fileGrid)
    fftGrid =  new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fftGrid =  new FFTGrid(nx, ny, nz, nxp, nyp, nzp);

  return(fftGrid);
}

int
ModelGeneral::computeTime(int year, int month, int day) const
{
  if (year == IMISSING)
    return(0);

  int deltaYear = year-1900; //Ok baseyear.
  int time = 365*deltaYear+deltaYear/4; //Leap years.
  if (month == IMISSING)
    time += 182;
  else {
    std::vector<int> accDays(12,0);
    accDays[1]  = accDays[0]  + 31;
    accDays[2]  = accDays[1]  + 28;
    accDays[3]  = accDays[2]  + 31;
    accDays[4]  = accDays[3]  + 30;
    accDays[5]  = accDays[4]  + 31;
    accDays[6]  = accDays[5]  + 30;
    accDays[7]  = accDays[6]  + 31;
    accDays[8]  = accDays[7]  + 31;
    accDays[9]  = accDays[8]  + 30;
    accDays[10] = accDays[9]  + 31;
    accDays[11] = accDays[10] + 30;

    time += accDays[month];

    if (day == IMISSING)
      time += 15;
    else
      time += day;
  }
  return(time);
}

//void
//ModelGeneral::generateRockPhysics3DBackground(const std::vector<DistributionsRock *>           & rock_distribution,
//                                              const std::vector<float>                         & probability,
//                                              FFTGrid                                          & vp,
//                                              FFTGrid                                          & vs,
//                                              FFTGrid                                          & rho)
//{
//  // Set up of expectations grids
//
//  // Variables for looping through FFTGrids
//  const int nz   = vp.getNz();
//  const int ny   = vp.getNy();
//  const int nx   = vp.getNx();
//  const int nzp  = vp.getNzp();
//  const int nyp  = vp.getNyp();
//  const int nxp = vp.getNxp();
//  const int rnxp = vp.getRNxp();
//
//  LogKit::LogFormatted(LogKit::Low,"\nGenerating background model from rock physics:\n");
//
//  float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
//  float nextMonitor = monitorSize;
//  std::cout
//    << "\n  0%       20%       40%       60%       80%      100%"
//    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
//    << "\n  ^";
//
//  const size_t number_of_facies = probability.size();
//
//  // Temporary grids for storing top and base values of (vp,vs,rho) for use in linear interpolation in the padding
//  NRLib::Grid2D<float> topVp  (nx, ny, 0.0);
//  NRLib::Grid2D<float> topVs  (nx, ny, 0.0);
//  NRLib::Grid2D<float> topRho (nx, ny, 0.0);
//  NRLib::Grid2D<float> baseVp (nx ,ny, 0.0);
//  NRLib::Grid2D<float> baseVs (nx, ny, 0.0);
//  NRLib::Grid2D<float> baseRho(nx, ny, 0.0);
//
//  vp.setAccessMode(FFTGrid::WRITE);
//  vs.setAccessMode(FFTGrid::WRITE);
//  rho.setAccessMode(FFTGrid::WRITE);
//
//  // Loop through all cells in the FFTGrids
//  for (int k = 0; k < nzp; k++) {
//    for (int j = 0; j < nyp; j++) {
//      for (int i = 0; i < rnxp; i++) {
//
//        // If outside/If in the padding in x- and y-direction,
//        // set expectation equal to something at right scale
//        // (top value for closest edge)
//        // NBNB OK Can be made better linear interoplation between first and last value in i an j direction as well
//        if (i >= nx || j >= ny) {
//          int indexI;
//          int indexJ;
//          indexI = i > (nx+nxp)/2 ? 0   : nx-1;
//          indexJ = j > (ny+nyp)/2 ? 0   : ny-1;
//          indexI = std::min(i,indexI);
//          indexJ = std::min(j,indexJ);
//
//          float vpVal  = topVp(indexI,indexJ);
//          float vsVal  = topVs(indexI,indexJ);
//          float rhoVal = topRho(indexI,indexJ);
//          vp.setNextReal(vpVal);
//          vs.setNextReal(vsVal);
//          rho.setNextReal(rhoVal);
//        }
//
//        // If outside in z-direction, use linear interpolation between top and base values of the expectations
//        else if (k >= nz) {
//          double t  = double(nzp-k+1)/(nzp-nz+1);
//          double vpVal =  topVp(i,j)*t  + baseVp(i,j)*(1-t);
//          double vsVal =  topVs(i,j)*t  + baseVs(i,j)*(1-t);
//          double rhoVal = topRho(i,j)*t + baseRho(i,j)*(1-t);
//
//          // Set interpolated values in expectation grids
//          vp.setNextReal(static_cast<float>(vpVal));
//          vs.setNextReal(static_cast<float>(vsVal));
//          rho.setNextReal(static_cast<float>(rhoVal));
//        }
//
//        // Otherwise use trend values to get expectation values for each facies from the rock
//        else {
//          std::vector<double> trend_position = trend_cubes_.GetTrendPosition(i,j,k);
//
//          std::vector<float> expectations(3, 0);  // Antar initialisert til 0.
//
//          std::vector<std::vector<double> > expectation_m(number_of_facies);
//          for (size_t f = 0; f < number_of_facies; f++)
//            expectation_m[f] = rock_distribution[f]->GetLogExpectation(trend_position);
//
//          // Sum up for all facies: probability for a facies multiplied with the expectations of (vp, vs, rho) given the facies
//          for (size_t f = 0; f < number_of_facies; f++){
//            expectations[0] += static_cast<float>(expectation_m[f][0] * probability[f]);
//            expectations[1] += static_cast<float>(expectation_m[f][1] * probability[f]);
//            expectations[2] += static_cast<float>(expectation_m[f][2] * probability[f]);
//          }
//
//          // Set values in expectation grids
//          vp.setNextReal(expectations[0]);
//          vs.setNextReal(expectations[1]);
//          rho.setNextReal(expectations[2]);
//
//          // Store top and base values of the expectations for later use in interpolation in the padded region.
//          if (k==0) {
//            topVp(i,j)  = expectations[0];
//            topVs(i,j)  = expectations[1];
//            topRho(i,j) = expectations[2];
//          }
//          else if (k==nz-1) {
//            baseVp(i,j)  = expectations[0];
//            baseVs(i,j)  = expectations[1];
//            baseRho(i,j) = expectations[2];
//          }
//        }
//      }
//    }
//
//    // Log progress
//    if (k+1 >= static_cast<int>(nextMonitor) && k < nz) {
//      nextMonitor += monitorSize;
//      std::cout << "^";
//      fflush(stdout);
//    }
//  }
//
//  vp.endAccess();
//  vs.endAccess();
//  rho.endAccess();
//}

void
ModelGeneral::calculateCovariancesFromRockPhysics(const std::vector<DistributionsRock *>           & rock_distribution,
                                                  const std::vector<float>                         & probability,
                                                  NRLib::Grid2D<double>                            & param_corr,
                                                  std::string                                      & errTxt)
{

  LogKit::LogFormatted(LogKit::Low,"\nGenerating covariances from rock physics\n");

  bool has_trend = false;
  for (size_t i=0; i<rock_distribution.size(); i++) {
    std::vector<bool> rock_has_trend = rock_distribution[i]->HasTrend();

    for (int j=0; j<2; j++) {
      if (rock_has_trend[j] == true)
        has_trend = true;
    }
  }

  if (has_trend == true) {

    std::vector<int> trend_cube_size = trend_cubes_.GetSizeTrendCubes();

    const int nx = trend_cube_size[0];
    const int ny = trend_cube_size[1];
    const int nz = trend_cube_size[2];

    int modulus = 100;

    float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
    float nextMonitor = monitorSize;
    std::cout
      << "\n  0%       20%       40%       60%       80%      100%"
      << "\n  |    |    |    |    |    |    |    |    |    |    |  "
      << "\n  ^";


    // Local storage for summed combined variances
    NRLib::Grid2D<double> sumVariance(3,3,0);

    int n_samples = 0;

    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {

          if ( ( (i+1)*(j+1)*(k+1) ) % modulus == 0) {

            std::vector<double> trend_position = trend_cubes_.GetTrendPosition(i,j,k);

            NRLib::Grid2D<double> sigma_sum(3,3,0);

            calculateCovarianceInTrendPosition(rock_distribution,
                                               probability,
                                               trend_position,
                                               sigma_sum);


            for (size_t a=0; a<3; a++){
              for (size_t b=0; b<3; b++)
                sumVariance(a,b) += sigma_sum(a,b);
            }

            n_samples++;
          }
        }
      }

      // Log progress
      if (k+1 >= static_cast<int>(nextMonitor) && k < nz) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }

    if (n_samples > 0) {
      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++)
          param_corr(i,j) = sumVariance(i,j)/n_samples;
      }
    }

    else
      errTxt += "Could not build a covariance structure from rock physics.\n";
  }

  else {
    std::vector<double> trend_position(2, 0.0);

    NRLib::Grid2D<double> sigma_sum(3,3,0);

    calculateCovarianceInTrendPosition(rock_distribution,
                                       probability,
                                       trend_position,
                                       sigma_sum);


    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++)
        param_corr(i,j) = sigma_sum(i,j);
    }

  }
}

void
ModelGeneral::calculateCovarianceInTrendPosition(const std::vector<DistributionsRock *> & rock_distribution,
                                                 const std::vector<float>               & probability,
                                                 const std::vector<double>              & trend_position,
                                                 NRLib::Grid2D<double>                  & sigma_sum) const
{
  int number_of_facies = static_cast<int>(rock_distribution.size());

  std::vector<std::vector<double> > expectation_m(number_of_facies);

  for (int f = 0; f < number_of_facies; f++)
    expectation_m[f] = rock_distribution[f]->GetLogExpectation(trend_position);

  // Sum up for all facies: probability for a facies multiplied with the expectations of (vp, vs, rho) given the facies
  std::vector<float> expectations(3, 0);
  for (int f = 0; f < number_of_facies; f++){
    expectations[0] += static_cast<float>(expectation_m[f][0] * probability[f]);
    expectations[1] += static_cast<float>(expectation_m[f][1] * probability[f]);
    expectations[2] += static_cast<float>(expectation_m[f][2] * probability[f]);
  }

  // Compute combined variance for all facies in the given grid cell.
  // Calculation of variance for a set of pdfs with given probability in the rock physics model:
  //
  // Var(X) = E([X-E(X)]^2) = E([X - E(X|facies)]^2) + E([E(X|facies) -E(X)]^2) = E(Var(X|facies)) + Var(E(X|facies))
  //        = Sum_{over all facies} (probability of facies * variance given facies) + sum_{over all facies} probability of facies * (expected value given facies - EX)^2,
  // where EX is the sum of probability of a facies multiplied with expectation of \mu given facies
  //
  // For all facies: Summing up expected value of variances and variance of expected values

  for (int f = 0; f < number_of_facies; f++) {
    NRLib::Grid2D<double> sigma = rock_distribution[f]->GetLogCovariance(trend_position);

    // For all elements in the 3x3 matrix of the combined variance
    for (size_t a=0; a<3; a++) {
      for (size_t b=0; b<3; b++) {
        double sigma_weigth = probability[f] * (sigma(a,b) + (expectation_m[f][a] - expectations[a])*(expectation_m[f][b] - expectations[b]));
        sigma_sum(a,b)     += sigma_weigth;
      }
    }
  }
}


//void
//ModelGeneral::setUp3DPartOf4DBackground(const std::vector<DistributionsRock *>           & rock,
//                                        const std::vector<float>                         & probability,
//                                        const Simbox                                     & timeSimbox,
//                                        const ModelSettings                              & model_settings,
//                                        SeismicParametersHolder                          & seismicParameters,
//                                        State4D                                          & state4d,
//                                        std::string                                      & /*errTxt*/)
//{
//  LogKit::WriteHeader("Prior Expectations / Background Model");
//
//  // Allocates the static mu grid: 3 grids.
//
//  // Static grids for 4D inversion, filled with 3D rock physics background
//  FFTGrid * vp_stat;
//  FFTGrid * vs_stat;
//  FFTGrid * rho_stat;
//
//  // Parameters for generating new FFTGrids
//  const int nx    = timeSimbox.getnx();
//  const int ny    = timeSimbox.getny();
//  const int nz    = timeSimbox.getnz();
//  const int nxPad = timeSimbox.GetNXpad();
//  const int nyPad = timeSimbox.GetNYpad();
//  const int nzPad = timeSimbox.GetNZpad();
//
//  // Creating grids for mu static
//  vp_stat  = createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, model_settings.getFileGrid());
//  vs_stat  = createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, model_settings.getFileGrid());
//  rho_stat = createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, model_settings.getFileGrid());
//
//  vp_stat ->createRealGrid();
//  vs_stat ->createRealGrid();
//  rho_stat->createRealGrid();
//
//  // For the static variables, generate expectation grids and variance coefficients from the 3D settings.
//  generateRockPhysics3DBackground(rock, probability, *vp_stat, *vs_stat, *rho_stat);
//
//  vp_stat ->setType(FFTGrid::PARAMETER);
//  vs_stat ->setType(FFTGrid::PARAMETER);
//  rho_stat->setType(FFTGrid::PARAMETER);
//
//  state4d.setStaticMu(vp_stat, vs_stat, rho_stat);
//
//  seismicParameters.copyBackgroundParameters(vp_stat, vs_stat, rho_stat);
//
//}

void
ModelGeneral::CopyCorrelationsTo4DState(SeismicParametersHolder                          & seismicParameters,
                                       State4D                                          & state4d)
{
  // Allocates the static sigma grids: 6 grids.

  // Static sigma
  FFTGrid * vp_vp_stat;
  FFTGrid * vp_vs_stat;
  FFTGrid * vp_rho_stat;
  FFTGrid * vs_vs_stat;
  FFTGrid * vs_rho_stat;
  FFTGrid * rho_rho_stat;

  // Copying grids for sigma static
  vp_vp_stat   = new FFTGrid( seismicParameters.GetCovVp());
  vp_vs_stat   = new FFTGrid( seismicParameters.GetCrCovVpVs());
  vp_rho_stat  = new FFTGrid( seismicParameters.GetCrCovVpRho());
  vs_vs_stat   = new FFTGrid( seismicParameters.GetCovVs());
  vs_rho_stat  = new FFTGrid( seismicParameters.GetCrCovVsRho());
  rho_rho_stat = new FFTGrid( seismicParameters.GetCovRho());

  state4d.setStaticSigma(vp_vp_stat, vp_vs_stat, vp_rho_stat, vs_vs_stat, vs_rho_stat, rho_rho_stat);
}

//void
//ModelGeneral::processWells(std::vector<WellData *> & wells,
//                             Simbox              * timeSimbox,
//                             ModelSettings      *& model_settings,
//                             const InputFiles    * input_files,
//                             std::string         & errText,
//                             bool                & failed)
//{
//  int     nWells         = model_settings->getNumberOfWells();
//
//  if (nWells > 0) {
//
//    double wall=0.0, cpu=0.0;
//    TimeKit::getTime(wall,cpu);
//
//    LogKit::WriteHeader("Reading and processing wells");
//
//    bool    faciesLogGiven = model_settings->getFaciesLogGiven();
//    int     nFacies        = 0;
//    int     error = 0;
//
//    std::string tmpErrText("");
//    wells.resize(nWells);
//    for (int i=0 ; i<nWells ; i++) {
//      wells[i] = new WellData(input_files->getWellFile(i),
//        model_settings->getLogNames(),
//        model_settings->getInverseVelocity(),
//        model_settings,
//        model_settings->getIndicatorFacies(i),
//        model_settings->getIndicatorFilter(i),
//        model_settings->getIndicatorWavelet(i),
//        model_settings->getIndicatorBGTrend(i),
//        model_settings->getIndicatorRealVs(i),
//        faciesLogGiven);
//      if (wells[i]->checkError(tmpErrText) != 0) {
//        errText += tmpErrText;
//        error = 1;
//      }
//    }
//
//
//    if (error == 0) {
//      if (model_settings->getFaciesLogGiven()) {
//        setFaciesNamesFromWells(wells, model_settings, tmpErrText, error);
//        nFacies = static_cast<int>(facies_names_.size()); // nFacies is set in setFaciesNames()
//      }
//      if (error>0)
//        errText += "Prior facies probabilities from wells failed.\n"+tmpErrText;
//
//      int   * validWells    = new int[nWells];
//      bool  * validIndex    = new bool[nWells];
//      int   * nMerges       = new int[nWells];
//      int   * nInvalidAlpha = new int[nWells];
//      int   * nInvalidBeta  = new int[nWells];
//      int   * nInvalidRho   = new int[nWells];
//      float * rankCorr      = new float[nWells];
//      float * devAngle      = new float[nWells];
//      int  ** faciesCount   = NULL;
//
//      if (nFacies > 0) {
//        faciesCount = new int * [nWells];
//        for (int i = 0 ; i < nWells ; i++)
//          faciesCount[i] = new int[nFacies];
//      }
//
//      int count = 0;
//      int nohit=0;
//      int empty=0;
//      int facieslognotok = 0;
//      int upwards=0;
//      LogKit::LogFormatted(LogKit::Low,"\n");
//      for (int i=0 ; i<nWells ; i++)
//      {
//        bool skip = false;
//        LogKit::LogFormatted(LogKit::Low,wells[i]->getWellname()+" : \n");
//        if (wells[i]!=NULL) {
//          if (wells[i]->checkSimbox(timeSimbox) == 1) {
//            skip = true;
//            nohit++;
//            TaskList::addTask("Consider increasing the inversion volume such that well "+wells[i]->getWellname()+ " can be included");
//          }
//          if (wells[i]->getNd() == 0) {
//            LogKit::LogFormatted(LogKit::Low,"  IGNORED (no log entries found)\n");
//            skip = true;
//            empty++;
//            TaskList::addTask("Check the log entries in well "+wells[i]->getWellname()+".");
//          }
//          if (wells[i]->isFaciesOk()==0) {
//            LogKit::LogFormatted(LogKit::Low,"   IGNORED (facies log has wrong entries)\n");
//            skip = true;
//            facieslognotok++;
//            TaskList::addTask("Check the facies logs in well "+wells[i]->getWellname()+".\n       The facies logs in this well are wrong and the well is ignored");
//          }
//          if (wells[i]->removeDuplicateLogEntries(timeSimbox, nMerges[i]) == false) {
//            LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is too far from monotonous in time)\n");
//            skip = true;
//            upwards++;
//            TaskList::addTask("Check the TWT log in well "+wells[i]->getWellname()+".\n       The well is moving too much upwards, and the well is ignored");
//          }
//          if (skip)
//            validIndex[i] = false;
//          else {
//            validIndex[i] = true;
//            wells[i]->setWrongLogEntriesUndefined(nInvalidAlpha[i], nInvalidBeta[i], nInvalidRho[i]);
//            wells[i]->filterLogs();
//            //wells[i]->findMeanVsVp(waveletEstimInterval_);
//            wells[i]->lookForSyntheticVsLog(rankCorr[i]);
//            wells[i]->calculateDeviation(devAngle[i], timeSimbox);
//
//            if (nFacies > 0)
//              wells[i]->countFacies(timeSimbox,faciesCount[i]);
//            validWells[count] = i;
//            count++;
//          }
//        }
//      }
//      //
//      // Write summary.
//      //
//      LogKit::LogFormatted(LogKit::Low,"\n");
//      LogKit::LogFormatted(LogKit::Low,"                                      Invalid                                    \n");
//      LogKit::LogFormatted(LogKit::Low,"Well                    Merges      Vp   Vs  Rho  synthVs/Corr    Deviated/Angle \n");
//      LogKit::LogFormatted(LogKit::Low,"---------------------------------------------------------------------------------\n");
//      for (int i=0 ; i<nWells ; i++) {
//        if (validIndex[i])
//          LogKit::LogFormatted(LogKit::Low,"%-23s %6d    %4d %4d %4d     %3s / %5.3f      %3s / %4.1f\n",
//          wells[i]->getWellname().c_str(),
//          nMerges[i],
//          nInvalidAlpha[i],
//          nInvalidBeta[i],
//          nInvalidRho[i],
//          (wells[i]->hasSyntheticVsLog() ? "yes" : " no"),
//          rankCorr[i],
//          (devAngle[i] > model_settings->getMaxDevAngle() ? "yes" : " no"),
//          devAngle[i]);
//        else
//          LogKit::LogFormatted(LogKit::Low,"%-23s      -       -    -    -       - /     -       -  /    -\n",
//          wells[i]->getWellname().c_str());
//      }
//
//      //
//      // Print facies count for each well
//      //
//      if (nFacies > 0) {
//        //
//        // Probabilities
//        //
//        LogKit::LogFormatted(LogKit::Low,"\nFacies distributions for each well: \n");
//        LogKit::LogFormatted(LogKit::Low,"\nWell                    ");
//        for (int i = 0 ; i < nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Low,"%12s ",facies_names_[i].c_str());
//        LogKit::LogFormatted(LogKit::Low,"\n");
//        for (int i = 0 ; i < 24+13*nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Low,"-");
//        LogKit::LogFormatted(LogKit::Low,"\n");
//        for (int i = 0 ; i < nWells ; i++) {
//          if (validIndex[i]) {
//            float tot = 0.0;
//            for (int f = 0 ; f < nFacies ; f++)
//              tot += static_cast<float>(faciesCount[i][f]);
//            LogKit::LogFormatted(LogKit::Low,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++) {
//              if (tot > 0) {
//                float faciesProb = static_cast<float>(faciesCount[i][f])/tot;
//                LogKit::LogFormatted(LogKit::Low,"%12.4f ",faciesProb);
//              }
//              else
//                LogKit::LogFormatted(LogKit::Low,"         -   ");
//            }
//            LogKit::LogFormatted(LogKit::Low,"\n");
//          }
//          else {
//            LogKit::LogFormatted(LogKit::Low,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++)
//              LogKit::LogFormatted(LogKit::Low,"         -   ");
//            LogKit::LogFormatted(LogKit::Low,"\n");
//
//          }
//        }
//        LogKit::LogFormatted(LogKit::Low,"\n");
//        //
//        // Counts
//        //
//        LogKit::LogFormatted(LogKit::Medium,"\nFacies counts for each well: \n");
//        LogKit::LogFormatted(LogKit::Medium,"\nWell                    ");
//        for (int i = 0 ; i < nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Medium,"%12s ",facies_names_[i].c_str());
//        LogKit::LogFormatted(LogKit::Medium,"\n");
//        for (int i = 0 ; i < 24+13*nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Medium,"-");
//        LogKit::LogFormatted(LogKit::Medium,"\n");
//        for (int i = 0 ; i < nWells ; i++) {
//          if (validIndex[i]) {
//            float tot = 0.0;
//            for (int f = 0 ; f < nFacies ; f++)
//              tot += static_cast<float>(faciesCount[i][f]);
//            LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++) {
//              LogKit::LogFormatted(LogKit::Medium,"%12d ",faciesCount[i][f]);
//            }
//            LogKit::LogFormatted(LogKit::Medium,"\n");
//          }
//          else {
//            LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++)
//              LogKit::LogFormatted(LogKit::Medium,"         -   ");
//            LogKit::LogFormatted(LogKit::Medium,"\n");
//
//          }
//        }
//        LogKit::LogFormatted(LogKit::Medium,"\n");
//      }
//
//      //
//      // Remove invalid wells
//      //
//      for (int i=0 ; i<nWells ; i++)
//        if (!validIndex[i])
//          delete wells[i];
//      for (int i=0 ; i<count ; i++)
//        wells[i] = wells[validWells[i]];
//      for (int i=count ; i<nWells ; i++)
//        wells[i] = NULL;
//      nWells = count;
//      model_settings->setNumberOfWells(nWells);
//
//      delete [] validWells;
//      delete [] validIndex;
//      delete [] nMerges;
//      delete [] nInvalidAlpha;
//      delete [] nInvalidBeta;
//      delete [] nInvalidRho;
//      delete [] rankCorr;
//      delete [] devAngle;
//
//      if (nohit>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) do not hit the inversion volume and will be ignored.\n",nohit);
//      if (empty>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) contain no log entries and will be ignored.\n",empty);
//      if (facieslognotok>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) have wrong facies logs and will be ignored.\n",facieslognotok);
//      if (upwards>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) are moving upwards in TWT and will be ignored.\n",upwards);
//      if (nWells==0 && model_settings->getNoWellNedded()==false) {
//        LogKit::LogFormatted(LogKit::Low,"\nERROR: There are no wells left for data analysis. Please check that the inversion area given");
//        LogKit::LogFormatted(LogKit::Low,"\n       below is correct. If it is not, you probably have problems with coordinate scaling.");
//        LogKit::LogFormatted(LogKit::Low,"\n                                   X0          Y0        DeltaX      DeltaY      Angle");
//        LogKit::LogFormatted(LogKit::Low,"\n       -------------------------------------------------------------------------------");
//        LogKit::LogFormatted(LogKit::Low,"\n       Inversion area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n",
//          timeSimbox->getx0(), timeSimbox->gety0(),
//          timeSimbox->getlx(), timeSimbox->getly(),
//          (timeSimbox->getAngle()*180)/M_PI);
//        errText += "No wells available for estimation.";
//        error = 1;
//      }
//
//      if (nFacies > 0) {
//        int fc;
//        for (int i = 0; i < nFacies; i++){
//          fc = 0;
//          for (int j = 0; j < nWells; j++){
//            fc+=faciesCount[j][i];
//          }
//          if (fc == 0){
//            LogKit::LogFormatted(LogKit::Low,"\nWARNING: Facies %s is not observed in any of the wells, and posterior facies probability can not be estimated for this facies.\n",facies_names_[i].c_str() );
//            TaskList::addTask("In order to estimate prior facies probability for facies "+ facies_names_[i] + " add wells which contain observations of this facies.\n");
//          }
//        }
//        for (int i = 0 ; i<nWells ; i++)
//          delete [] faciesCount[i];
//        delete [] faciesCount;
//      }
//
//    }
//    failed = error > 0;
//    Timings::setTimeWells(wall,cpu);
//  }
//}

//void
//ModelGeneral::setFaciesNamesFromWells(std::vector<WellData *>     wells,
//                                      ModelSettings            *& model_settings,
//                                      std::string               & tmpErrText,
//                                      int                       & error)
//{
//  int min,max;
//  int globalmin = 0;
//  int globalmax = 0;
//  bool first = true;
//  for (int w = 0; w < model_settings->getNumberOfWells(); w++) {
//    if (wells[w]->isFaciesLogDefined())
//    {
//      wells[w]->getMinMaxFnr(min,max);
//      if (first==true)
//      {
//        globalmin = min;
//        globalmax = max;
//        first = false;
//      }
//      else
//      {
//        if (min<globalmin)
//          globalmin = min;
//        if (max>globalmax)
//          globalmax = max;
//      }
//    }
//  }
//
//  int nnames = globalmax - globalmin + 1;
//  std::vector<std::string> names(nnames);
//
//  for (int w=0 ; w<model_settings->getNumberOfWells() ; w++)
//  {
//    if (wells[w]->isFaciesLogDefined())
//    {
//      for (int i=0 ; i < wells[w]->getNFacies() ; i++)
//      {
//        std::string name = wells[w]->getFaciesName(i);
//        int         fnr  = wells[w]->getFaciesNr(i) - globalmin;
//
//        if (names[fnr] == "") {
//          names[fnr] = name;
//        }
//        else if (names[fnr] != name)
//        {
//          tmpErrText += "Problem with facies logs. Facies names and numbers are not uniquely defined.\n";
//          error++;
//        }
//      }
//    }
//  }
//
//  LogKit::LogFormatted(LogKit::Low,"\nFaciesLabel      FaciesName           ");
//  LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------\n");
//  for (int i=0 ; i<nnames ; i++)
//    if (names[i] != "")
//      LogKit::LogFormatted(LogKit::Low,"    %2d           %-20s\n",i+globalmin,names[i].c_str());
//
//  int nFacies = 0;
//  for (int i=0 ; i<nnames ; i++)
//    if (names[i] != "")
//      nFacies++;
//
//  for (int i=0 ; i<nnames ; i++) {
//    if (names[i] != "") {
//      faciesLabels_.push_back(globalmin + i);
//      facies_names_.push_back(names[i]);
//    }
//  }
//}

//void
//ModelGeneral::processWellLocation(FFTGrid                       ** seisCube,
//                                    float                       ** reflectionMatrix,
//                                    ModelSettings                * model_settings,
//                                    const std::vector<Surface *> & interval)
//{
//  LogKit::WriteHeader("Estimating optimized well location");
//
//  double  deltaX, deltaY;
//  float   sum;
//  float   kMove;
//  float   moveAngle;
//  int     iMove;
//  int     jMove;
//  int     i,j,w;
//  int     iMaxOffset;
//  int     jMaxOffset;
//  int     nMoveAngles = 0;
//  int     nWells      = model_settings->getNumberOfWells();
//  int     nAngles     = model_settings->getNumberOfAngles(0);//Well location is not estimated when using time lapse data
//  float   maxShift    = model_settings->getMaxWellShift();
//  float   maxOffset   = model_settings->getMaxWellOffset();
//  double  angle       = simbox_->getAngle();
//  double  dx          = simbox_->getdx();
//  double  dy          = simbox_->getdy();
//  std::vector<float> seismicAngle = model_settings->getAngle(0); //Use first time lapse as this not is allowed in 4D
//
//  std::vector<float> angleWeight(nAngles);
//  LogKit::LogFormatted(LogKit::Low,"\n");
//  LogKit::LogFormatted(LogKit::Low,"  Well             Shift[ms]       DeltaI   DeltaX[m]   DeltaJ   DeltaY[m] \n");
//  LogKit::LogFormatted(LogKit::Low,"  ----------------------------------------------------------------------------------\n");
//
//  for (w = 0 ; w < nWells ; w++) {
//    if ( wells_[w]->isDeviated()==true )
//      continue;
//
//    BlockedLogs * bl = wells_[w]->getBlockedLogsOrigThick();
//    nMoveAngles = model_settings->getNumberOfWellAngles(w);
//
//    if ( nMoveAngles==0 )
//      continue;
//
//    for ( i=0; i<nAngles; i++ )
//      angleWeight[i] = 0;
//
//    for ( i=0; i<nMoveAngles; i++ ){
//      moveAngle   = model_settings->getWellMoveAngle(w,i);
//
//      for ( j=0; j<nAngles; j++ ){
//        if ( moveAngle == seismicAngle[j]){
//          angleWeight[j] = model_settings->getWellMoveWeight(w,i);
//          break;
//        }
//      }
//    }
//
//    sum = 0;
//    for ( i=0; i<nAngles; i++ )
//      sum += angleWeight[i];
//    if ( sum == 0 )
//      continue;
//
//    iMaxOffset = static_cast<int>(std::ceil(maxOffset/dx));
//    jMaxOffset = static_cast<int>(std::ceil(maxOffset/dy));
//
//    bl->findOptimalWellLocation(seisCube,simbox_,reflectionMatrix,nAngles,angleWeight,maxShift,iMaxOffset,jMaxOffset,interval,iMove,jMove,kMove);
//
//    deltaX = iMove*dx*cos(angle) - jMove*dy*sin(angle);
//    deltaY = iMove*dx*sin(angle) + jMove*dy*cos(angle);
//    wells_[w]->moveWell(simbox_,deltaX,deltaY,kMove);
//    wells_[w]->deleteBlockedLogsOrigThick();
//    wells_[w]->setBlockedLogsOrigThick( new BlockedLogs(wells_[w], simbox_, model_settings->getRunFromPanel()) );
//    LogKit::LogFormatted(LogKit::Low,"  %-13s %11.2f %12d %11.2f %8d %11.2f \n",
//    wells_[w]->getWellname().c_str(), kMove, iMove, deltaX, jMove, deltaY);
//  }
//
//   for (w = 0 ; w < nWells ; w++){
//     nMoveAngles = model_settings->getNumberOfWellAngles(w);
//
//    if ( wells_[w]->isDeviated()==true && nMoveAngles > 0 )
//    {
//      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Well %7s is treated as deviated and can not be moved.\n",
//          wells_[w]->getWellname().c_str());
//      TaskList::addTask("Well "+NRLib::ToString(wells_[w]->getWellname())+" can not be moved. Remove <optimize-location-to> for this well");
//    }
//   }
//}

//void
//ModelGeneral::processPriorCorrelations(Background                     * background,
//                                       std::vector<WellData *>          wells,
//                                       const Simbox                   * timeSimbox,
//                                       const ModelSettings            * model_settings,
//                                       const std::vector<float>       & priorFacies,
//                                       FFTGrid                       ** seisCube,
//                                       const InputFiles               * input_files,
//                                       SeismicParametersHolder        & seismicParameters,
//                                       std::string                    & errText,
//                                       bool                           & failed)
//{
//  bool printResult = ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0 ||
//                      model_settings->getEstimationMode() == true);
//  if (model_settings->getDoInversion() || printResult)
//  {
//    LogKit::WriteHeader("Prior Covariance");
//
//    double wall=0.0, cpu=0.0;
//    TimeKit::getTime(wall,cpu);
//
//    const std::string & paramCovFile = input_files->getParamCorrFile();
//    const std::string & corrTFile    = input_files->getTempCorrFile();
//
//    bool estimateParamCov = paramCovFile == "";
//    bool estimateTempCorr = (corrTFile    == "" && model_settings->getUseVerticalVariogram() == false);
//
//    //
//    // Read parameter covariance (Var0) from file or set from output from function generateRockPhysics3DBackground.
//    // Consistency check that only one option (file or rock physics) is possible, is done in XmlModelFile::checkInversionConsistency
//    //
//    float ** paramCorr = NULL;
//    bool failedParamCorr = false;
//    std::string tmpErrText("");
//    if (!estimateParamCov) {
//      paramCorr = ModelAVODynamic::readMatrix(paramCovFile, 3, 3, "parameter covariance", tmpErrText);
//      validateCorrelationMatrix(paramCorr, model_settings, tmpErrText);
//      if (paramCorr == NULL || tmpErrText != "") {
//        errText += "Reading of file "+paramCovFile+" for parameter covariance matrix failed\n";
//        errText += tmpErrText;
//        failedParamCorr = true;
//      }
//    }
//    else if (model_settings->getFaciesProbFromRockPhysics() == true) {
//      estimateParamCov = false;
//      paramCorr = new float * [3];
//      for (int i=0;i<3;i++) {
//        paramCorr[i] = new float[3];
//      }
//
//      int n_facies = static_cast<int>(facies_names_.size());
//
//      std::vector<DistributionsRock *> rock_distribution(n_facies);
//      typedef std::map<std::string, DistributionsRock *> rfMapType;
//      rfMapType rfMap = getRockDistributionTime0();
//
//      for (int i=0; i<n_facies; i++) {
//        rfMapType::iterator iter = rfMap.find(facies_names_[i]);
//        if (iter != rfMap.end())
//          rock_distribution[i] = iter->second;
//      }
//
//      NRLib::Grid2D<double> param_corr(3,3);
//      calculateCovariancesFromRockPhysics(rock_distribution,
//                                          priorFacies,
//                                          param_corr,
//                                          errText);
//
//
//      for (int i=0; i<3; i++) {
//        for (int j=0; j<3; j++)
//          paramCorr[i][j] = static_cast<float>(param_corr(i,j));
//      }
//
//      if (tmpErrText != "")
//      {
//        errText += "Parameter covariance matrix from rock physics failed\n";
//        errText += tmpErrText;
//        failedParamCorr = true;
//      }
//    }
//
//    //
//    // Estimate lateral correlation from seismic data
//    //
//    priorCorrXY_ = findCorrXYGrid(timeSimbox, model_settings);
//
//    if (model_settings->getLateralCorr()==NULL) // NBNB-PAL: this will never be true (default lateral corr)
//    {
//      int timelapse = 0; // Setting timelapse = 0 as this is the generation of prior model
//      estimateCorrXYFromSeismic(priorCorrXY_, seisCube, model_settings->getNumberOfAngles(timelapse));
//    }
//
//    int nCorrT = model_settings->getNZpad();
//    if ((nCorrT % 2) == 0)
//      nCorrT = nCorrT/2+1;
//    else
//      nCorrT = nCorrT/2;
//
//    std::vector<float> corrT;
//
//    bool failedTempCorr = false;
//    if (!estimateTempCorr)
//    {
//      if (model_settings->getUseVerticalVariogram() == true) {
//        corrT.resize(nCorrT+1);
//        float tempCorrRange = model_settings->getTemporalCorrelationRange();
//        float dz = static_cast<float>(timeSimbox->getdz());
//        for (int i=0; i<=nCorrT; i++){
//          //using an exponential variogram with a = 1/3 (Chiles and Delfiner 1999)
//          corrT[i] = exp(-3*dz*i/tempCorrRange);
//        }
//      }
//      else{
//        std::string tmpErrText("");
//        float ** corrMat = ModelAVODynamic::readMatrix(corrTFile, 1, nCorrT+1, "temporal correlation", tmpErrText);
//        if (corrMat == NULL)
//        {
//          errText += "Reading of file '"+corrTFile+"' for temporal correlation failed\n";
//          errText += tmpErrText;
//          failedTempCorr = true;
//        }
//        corrT.resize(nCorrT);
//        if (!failedTempCorr)
//        {
//          for (int i=0;i<nCorrT;i++)
//            corrT[i] = corrMat[0][i+1];
//          delete [] corrMat[0];
//          delete [] corrMat;
//        }
//      }
//    }
//
//    float ** pointVar0 = NULL;
//    if (estimateParamCov || estimateTempCorr) //Need well estimation
//    {
//      std::string tmpErrTxt;
//      Analyzelog * analyze = new Analyzelog(wells,
//                                            background,
//                                            timeSimbox,
//                                            model_settings,
//                                            tmpErrTxt);
//      if (tmpErrTxt != "") {
//        errText += tmpErrTxt;
//        failedParamCorr = true;
//      }
//
//      if (estimateParamCov)
//        paramCorr = analyze->getVar0();
//      else
//        delete [] analyze->getVar0();
//
//      pointVar0 = analyze->getPointVar0();
//
//      float * estCorrT = analyze->getCorrT();
//      if (estimateTempCorr) {
//        corrT.resize(nCorrT);
//        int nEst = analyze->getNumberOfLags();
//        int i, max = nEst;
//        if (max > nCorrT)
//          max = nCorrT;
//        for (i=0;i<max;i++)
//          corrT[i] = estCorrT[i];
//        if (i<nCorrT) {
//          LogKit::LogFormatted(LogKit::High,
//            "\nOnly able to estimate %d of %d lags needed in temporal correlation. The rest are set to 0.\n", nEst, nCorrT);
//          for (;i<nCorrT;i++)
//            corrT[i] = 0.0f;
//        }
//      }
//      delete [] estCorrT;
//
//      delete analyze;
//    }
//
//    if (failedParamCorr || failedTempCorr)
//      failed = true;
//
//    if (!failed) {
//
//      const int nx        = timeSimbox->getnx();
//      const int ny        = timeSimbox->getny();
//      const int nz        = timeSimbox->getnz();
//      const int nxPad     = model_settings->getNXpad();
//      const int nyPad     = model_settings->getNYpad();
//      const int nzPad     = model_settings->getNZpad();
//
//      float dt = static_cast<float>(timeSimbox->getdz());
//      float lowCut = model_settings->getLowCut();
//      int lowIntCut = int(floor(lowCut*(nzPad*0.001*dt))); // computes the integer whis corresponds to the low cut frequency.
//
//      float corrGradI;
//      float corrGradJ;
//      getCorrGradIJ(corrGradI, corrGradJ);
//      //makeCorr2DPositiveDefinite( priorCorrXY_);
//
//      seismicParameters.setCorrelationParameters(paramCorr,
//                                                 corrT,
//                                                 priorCorrXY_,
//                                                 lowIntCut,
//                                                 corrGradI,
//                                                 corrGradJ,
//                                                 nx,
//                                                 ny,
//                                                 nz,
//                                                 nxPad,
//                                                 nyPad,
//                                                 nzPad);
//
//      for (int i=0; i<3; i++)
//        delete [] paramCorr[i];
//      delete [] paramCorr;
//
//      if (printResult)
//        seismicParameters.writeFilePriorVariances(model_settings, corrT, priorCorrXY_, dt);
//      seismicParameters.printPriorVariances();
//    }
//
//
//    if (failedTempCorr == true || failedParamCorr == true)
//    {
//      errText += "Could not construct prior covariance. Unknown why...\n";
//      failed = true;
//    }
//
//    Timings::setTimePriorCorrelation(wall,cpu);
//  }
//}

//    bool    faciesLogGiven = modelSettings->getFaciesLogGiven();
//    int     nFacies        = 0;
//    int     error = 0;
//
//    std::string tmpErrText("");
//    wells.resize(nWells);
//    for(int i=0 ; i<nWells ; i++) {
//      wells[i] = new WellData(inputFiles->getWellFile(i),
//        modelSettings->getLogNames(),
//        modelSettings->getInverseVelocity(),
//        modelSettings,
//        modelSettings->getIndicatorFacies(i),
//        modelSettings->getIndicatorFilter(i),
//        modelSettings->getIndicatorWavelet(i),
//        modelSettings->getIndicatorBGTrend(i),
//        modelSettings->getIndicatorRealVs(i),
//        faciesLogGiven);
//      if(wells[i]->checkError(tmpErrText) != 0) {
//        errText += tmpErrText;
//        error = 1;
//      }
//    }
//
//
//    if (error == 0) {
//      if(modelSettings->getFaciesLogGiven()) {
//        setFaciesNamesFromWells(wells, modelSettings, tmpErrText, error);
//        nFacies = static_cast<int>(faciesNames_.size()); // nFacies is set in setFaciesNames()
//      }
//      if (error>0)
//        errText += "Prior facies probabilities from wells failed.\n"+tmpErrText;
//
//      int   * validWells    = new int[nWells];
//      bool  * validIndex    = new bool[nWells];
//      int   * nMerges       = new int[nWells];
//      int   * nInvalidAlpha = new int[nWells];
//      int   * nInvalidBeta  = new int[nWells];
//      int   * nInvalidRho   = new int[nWells];
//      float * rankCorr      = new float[nWells];
//      float * devAngle      = new float[nWells];
//      int  ** faciesCount   = NULL;
//
//      if(nFacies > 0) {
//        faciesCount = new int * [nWells];
//        for (int i = 0 ; i < nWells ; i++)
//          faciesCount[i] = new int[nFacies];
//      }
//
//      int count = 0;
//      int nohit=0;
//      int empty=0;
//      int facieslognotok = 0;
//      int upwards=0;
//      LogKit::LogFormatted(LogKit::Low,"\n");
//      for (int i=0 ; i<nWells ; i++)
//      {
//        bool skip = false;
//        LogKit::LogFormatted(LogKit::Low,wells[i]->getWellname()+" : \n");
//        if(wells[i]!=NULL) {
//          if(wells[i]->checkSimbox(timeSimbox) == 1) {
//            skip = true;
//            nohit++;
//            TaskList::addTask("Consider increasing the inversion volume such that well "+wells[i]->getWellname()+ " can be included");
//          }
//          if(wells[i]->getNd() == 0) {
//            LogKit::LogFormatted(LogKit::Low,"  IGNORED (no log entries found)\n");
//            skip = true;
//            empty++;
//            TaskList::addTask("Check the log entries in well "+wells[i]->getWellname()+".");
//          }
//          if(wells[i]->isFaciesOk()==0) {
//            LogKit::LogFormatted(LogKit::Low,"   IGNORED (facies log has wrong entries)\n");
//            skip = true;
//            facieslognotok++;
//            TaskList::addTask("Check the facies logs in well "+wells[i]->getWellname()+".\n       The facies logs in this well are wrong and the well is ignored");
//          }
//          if(wells[i]->removeDuplicateLogEntries(timeSimbox, nMerges[i]) == false) {
//            LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is too far from monotonous in time)\n");
//            skip = true;
//            upwards++;
//            TaskList::addTask("Check the TWT log in well "+wells[i]->getWellname()+".\n       The well is moving too much upwards, and the well is ignored");
//          }
//          if(skip)
//            validIndex[i] = false;
//          else {
//            validIndex[i] = true;
//            wells[i]->setWrongLogEntriesUndefined(nInvalidAlpha[i], nInvalidBeta[i], nInvalidRho[i]);
//            wells[i]->filterLogs();
//            //wells[i]->findMeanVsVp(waveletEstimInterval_);
//            wells[i]->lookForSyntheticVsLog(rankCorr[i]);
//            wells[i]->calculateDeviation(devAngle[i], timeSimbox);
//
//            if (nFacies > 0)
//              wells[i]->countFacies(timeSimbox,faciesCount[i]);
//            validWells[count] = i;
//            count++;
//          }
//        }
//      }
//      //
//      // Write summary.
//      //
//      LogKit::LogFormatted(LogKit::Low,"\n");
//      LogKit::LogFormatted(LogKit::Low,"                                      Invalid                                    \n");
//      LogKit::LogFormatted(LogKit::Low,"Well                    Merges      Vp   Vs  Rho  synthVs/Corr    Deviated/Angle \n");
//      LogKit::LogFormatted(LogKit::Low,"---------------------------------------------------------------------------------\n");
//      for(int i=0 ; i<nWells ; i++) {
//        if (validIndex[i])
//          LogKit::LogFormatted(LogKit::Low,"%-23s %6d    %4d %4d %4d     %3s / %5.3f      %3s / %4.1f\n",
//          wells[i]->getWellname().c_str(),
//          nMerges[i],
//          nInvalidAlpha[i],
//          nInvalidBeta[i],
//          nInvalidRho[i],
//          (wells[i]->hasSyntheticVsLog() ? "yes" : " no"),
//          rankCorr[i],
//          (devAngle[i] > modelSettings->getMaxDevAngle() ? "yes" : " no"),
//          devAngle[i]);
//        else
//          LogKit::LogFormatted(LogKit::Low,"%-23s      -       -    -    -       - /     -       -  /    -\n",
//          wells[i]->getWellname().c_str());
//      }
//
//      //
//      // Print facies count for each well
//      //
//      if(nFacies > 0) {
//        //
//        // Probabilities
//        //
//        LogKit::LogFormatted(LogKit::Low,"\nFacies distributions for each well: \n");
//        LogKit::LogFormatted(LogKit::Low,"\nWell                    ");
//        for (int i = 0 ; i < nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Low,"%12s ",faciesNames_[i].c_str());
//        LogKit::LogFormatted(LogKit::Low,"\n");
//        for (int i = 0 ; i < 24+13*nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Low,"-");
//        LogKit::LogFormatted(LogKit::Low,"\n");
//        for (int i = 0 ; i < nWells ; i++) {
//          if (validIndex[i]) {
//            float tot = 0.0;
//            for (int f = 0 ; f < nFacies ; f++)
//              tot += static_cast<float>(faciesCount[i][f]);
//            LogKit::LogFormatted(LogKit::Low,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++) {
//              if (tot > 0) {
//                float faciesProb = static_cast<float>(faciesCount[i][f])/tot;
//                LogKit::LogFormatted(LogKit::Low,"%12.4f ",faciesProb);
//              }
//              else
//                LogKit::LogFormatted(LogKit::Low,"         -   ");
//            }
//            LogKit::LogFormatted(LogKit::Low,"\n");
//          }
//          else {
//            LogKit::LogFormatted(LogKit::Low,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++)
//              LogKit::LogFormatted(LogKit::Low,"         -   ");
//            LogKit::LogFormatted(LogKit::Low,"\n");
//
//          }
//        }
//        LogKit::LogFormatted(LogKit::Low,"\n");
//        //
//        // Counts
//        //
//        LogKit::LogFormatted(LogKit::Medium,"\nFacies counts for each well: \n");
//        LogKit::LogFormatted(LogKit::Medium,"\nWell                    ");
//        for (int i = 0 ; i < nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Medium,"%12s ",faciesNames_[i].c_str());
//        LogKit::LogFormatted(LogKit::Medium,"\n");
//        for (int i = 0 ; i < 24+13*nFacies ; i++)
//          LogKit::LogFormatted(LogKit::Medium,"-");
//        LogKit::LogFormatted(LogKit::Medium,"\n");
//        for (int i = 0 ; i < nWells ; i++) {
//          if (validIndex[i]) {
//            float tot = 0.0;
//            for (int f = 0 ; f < nFacies ; f++)
//              tot += static_cast<float>(faciesCount[i][f]);
//            LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++) {
//              LogKit::LogFormatted(LogKit::Medium,"%12d ",faciesCount[i][f]);
//            }
//            LogKit::LogFormatted(LogKit::Medium,"\n");
//          }
//          else {
//            LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells[i]->getWellname().c_str());
//            for (int f = 0 ; f < nFacies ; f++)
//              LogKit::LogFormatted(LogKit::Medium,"         -   ");
//            LogKit::LogFormatted(LogKit::Medium,"\n");
//
//          }
//        }
//        LogKit::LogFormatted(LogKit::Medium,"\n");
//      }
//
//      //
//      // Remove invalid wells
//      //
//      for(int i=0 ; i<nWells ; i++)
//        if (!validIndex[i])
//          delete wells[i];
//      for(int i=0 ; i<count ; i++)
//        wells[i] = wells[validWells[i]];
//      for(int i=count ; i<nWells ; i++)
//        wells[i] = NULL;
//      nWells = count;
//      modelSettings->setNumberOfWells(nWells);
//
//      delete [] validWells;
//      delete [] validIndex;
//      delete [] nMerges;
//      delete [] nInvalidAlpha;
//      delete [] nInvalidBeta;
//      delete [] nInvalidRho;
//      delete [] rankCorr;
//      delete [] devAngle;
//
//      if (nohit>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) do not hit the inversion volume and will be ignored.\n",nohit);
//      if (empty>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) contain no log entries and will be ignored.\n",empty);
//      if(facieslognotok>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) have wrong facies logs and will be ignored.\n",facieslognotok);
//      if(upwards>0)
//        LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) are moving upwards in TWT and will be ignored.\n",upwards);
//      if (nWells==0 && modelSettings->getNoWellNedded()==false) {
//        LogKit::LogFormatted(LogKit::Low,"\nERROR: There are no wells left for data analysis. Please check that the inversion area given");
//        LogKit::LogFormatted(LogKit::Low,"\n       below is correct. If it is not, you probably have problems with coordinate scaling.");
//        LogKit::LogFormatted(LogKit::Low,"\n                                   X0          Y0        DeltaX      DeltaY      Angle");
//        LogKit::LogFormatted(LogKit::Low,"\n       -------------------------------------------------------------------------------");
//        LogKit::LogFormatted(LogKit::Low,"\n       Inversion area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n",
//          timeSimbox->getx0(), timeSimbox->gety0(),
//          timeSimbox->getlx(), timeSimbox->getly(),
//          (timeSimbox->getAngle()*180)/M_PI);
//        errText += "No wells available for estimation.";
//        error = 1;
//      }
//
//      if(nFacies > 0) {
//        int fc;
//        for(int i = 0; i < nFacies; i++){
//          fc = 0;
//          for(int j = 0; j < nWells; j++){
//            fc+=faciesCount[j][i];
//          }
//          if(fc == 0){
//            LogKit::LogFormatted(LogKit::Low,"\nWARNING: Facies %s is not observed in any of the wells, and posterior facies probability can not be estimated for this facies.\n",faciesNames_[i].c_str() );
//            TaskList::addTask("In order to estimate prior facies probability for facies "+ faciesNames_[i] + " add wells which contain observations of this facies.\n");
//          }
//        }
//        for (int i = 0 ; i<nWells ; i++)
//          delete [] faciesCount[i];
//        delete [] faciesCount;
//      }
//
//    }
//    failed = error > 0;
//    Timings::setTimeWells(wall,cpu);
//  }
//}

//void
//ModelGeneral::setFaciesNamesFromWells(std::vector<WellData *>     wells,
//                                      ModelSettings            *& modelSettings,
//                                      std::string               & tmpErrText,
//                                      int                       & error)
//{
//  int min,max;
//  int globalmin = 0;
//  int globalmax = 0;
//  bool first = true;
//  for (int w = 0; w < modelSettings->getNumberOfWells(); w++) {
//    if(wells[w]->isFaciesLogDefined())
//    {
//      wells[w]->getMinMaxFnr(min,max);
//      if(first==true)
//      {
//        globalmin = min;
//        globalmax = max;
//        first = false;
//      }
//      else
//      {
//        if(min<globalmin)
//          globalmin = min;
//        if(max>globalmax)
//          globalmax = max;
//      }
//    }
//  }
//
//  int nnames = globalmax - globalmin + 1;
//  std::vector<std::string> names(nnames);
//
//  for(int w=0 ; w<modelSettings->getNumberOfWells() ; w++)
//  {
//    if(wells[w]->isFaciesLogDefined())
//    {
//      for(int i=0 ; i < wells[w]->getNFacies() ; i++)
//      {
//        std::string name = wells[w]->getFaciesName(i);
//        int         fnr  = wells[w]->getFaciesNr(i) - globalmin;
//
//        if(names[fnr] == "") {
//          names[fnr] = name;
//        }
//        else if(names[fnr] != name)
//        {
//          tmpErrText += "Problem with facies logs. Facies names and numbers are not uniquely defined.\n";
//          error++;
//        }
//      }
//    }
//  }
//
//  LogKit::LogFormatted(LogKit::Low,"\nFaciesLabel      FaciesName           ");
//  LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------\n");
//  for(int i=0 ; i<nnames ; i++)
//    if(names[i] != "")
//      LogKit::LogFormatted(LogKit::Low,"    %2d           %-20s\n",i+globalmin,names[i].c_str());
//
//  int nFacies = 0;
//  for(int i=0 ; i<nnames ; i++)
//    if(names[i] != "")
//      nFacies++;
//
//  for(int i=0 ; i<nnames ; i++) {
//    if(names[i] != "") {
//      faciesLabels_.push_back(globalmin + i);
//      faciesNames_.push_back(names[i]);
//    }
//  }
//}
//
//void
//ModelGeneral::processWellLocation(FFTGrid                       ** seisCube,
//                                    float                       ** reflectionMatrix,
//                                    ModelSettings                * modelSettings,
//                                    const std::vector<Surface *> & interval)
//{
//  LogKit::WriteHeader("Estimating optimized well location");
//
//  double  deltaX, deltaY;
//  float   sum;
//  float   kMove;
//  float   moveAngle;
//  int     iMove;
//  int     jMove;
//  int     i,j,w;
//  int     iMaxOffset;
//  int     jMaxOffset;
//  int     nMoveAngles = 0;
//  int     nWells      = modelSettings->getNumberOfWells();
//  int     nAngles     = modelSettings->getNumberOfAngles(0);//Well location is not estimated when using time lapse data
//  float   maxShift    = modelSettings->getMaxWellShift();
//  float   maxOffset   = modelSettings->getMaxWellOffset();
//  double  angle       = timeSimbox_->getAngle();
//  double  dx          = timeSimbox_->getdx();
//  double  dy          = timeSimbox_->getdy();
//  std::vector<float> seismicAngle = modelSettings->getAngle(0); //Use first time lapse as this not is allowed in 4D
//
//  std::vector<float> angleWeight(nAngles);
//  LogKit::LogFormatted(LogKit::Low,"\n");
//  LogKit::LogFormatted(LogKit::Low,"  Well             Shift[ms]       DeltaI   DeltaX[m]   DeltaJ   DeltaY[m] \n");
//  LogKit::LogFormatted(LogKit::Low,"  ----------------------------------------------------------------------------------\n");
//
//  for (w = 0 ; w < nWells ; w++) {
//    if( wells_[w]->isDeviated()==true )
//      continue;
//
//    BlockedLogs * bl = wells_[w]->getBlockedLogsOrigThick();
//    nMoveAngles = modelSettings->getNumberOfWellAngles(w);
//
//    if( nMoveAngles==0 )
//      continue;
//
//    for( i=0; i<nAngles; i++ )
//      angleWeight[i] = 0;
//
//    for( i=0; i<nMoveAngles; i++ ){
//      moveAngle   = modelSettings->getWellMoveAngle(w,i);
//
//      for( j=0; j<nAngles; j++ ){
//        if( moveAngle == seismicAngle[j]){
//          angleWeight[j] = modelSettings->getWellMoveWeight(w,i);
//          break;
//        }
//      }
//    }
//
//    sum = 0;
//    for( i=0; i<nAngles; i++ )
//      sum += angleWeight[i];
//    if( sum == 0 )
//      continue;
//
//    iMaxOffset = static_cast<int>(std::ceil(maxOffset/dx));
//    jMaxOffset = static_cast<int>(std::ceil(maxOffset/dy));
//
//    bl->findOptimalWellLocation(seisCube,timeSimbox_,reflectionMatrix,nAngles,angleWeight,maxShift,iMaxOffset,jMaxOffset,interval,iMove,jMove,kMove);
//
//    deltaX = iMove*dx*cos(angle) - jMove*dy*sin(angle);
//    deltaY = iMove*dx*sin(angle) + jMove*dy*cos(angle);
//    wells_[w]->moveWell(timeSimbox_,deltaX,deltaY,kMove);
//    wells_[w]->deleteBlockedLogsOrigThick();
//    wells_[w]->setBlockedLogsOrigThick( new BlockedLogs(wells_[w], timeSimbox_, modelSettings->getRunFromPanel()) );
//    LogKit::LogFormatted(LogKit::Low,"  %-13s %11.2f %12d %11.2f %8d %11.2f \n",
//    wells_[w]->getWellname().c_str(), kMove, iMove, deltaX, jMove, deltaY);
//  }
//
//   for (w = 0 ; w < nWells ; w++){
//     nMoveAngles = modelSettings->getNumberOfWellAngles(w);
//
//    if( wells_[w]->isDeviated()==true && nMoveAngles > 0 )
//    {
//      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Well %7s is treated as deviated and can not be moved.\n",
//          wells_[w]->getWellname().c_str());
//      TaskList::addTask("Well "+NRLib::ToString(wells_[w]->getWellname())+" can not be moved. Remove <optimize-location-to> for this well");
//    }
//   }
//}
//void
//ModelGeneral::processPriorCorrelations(Background                     * background,
//                                       std::vector<WellData *>          wells,
//                                       const Simbox                   * timeSimbox,
//                                       const ModelSettings            * modelSettings,
//                                       const std::vector<float>       & priorFacies,
//                                       FFTGrid                       ** seisCube,
//                                       const InputFiles               * inputFiles,
//                                       SeismicParametersHolder        & seismicParameters,
//                                       std::string                    & errText,
//                                       bool                           & failed)
//{
//  bool printResult = ((modelSettings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0 ||
//                      modelSettings->getEstimationMode() == true);
//  if (modelSettings->getDoInversion() || printResult)
//  {
//    LogKit::WriteHeader("Prior Covariance");
//
//    double wall=0.0, cpu=0.0;
//    TimeKit::getTime(wall,cpu);
//
//    const std::string & paramCovFile = inputFiles->getParamCorrFile();
//    const std::string & corrTFile    = inputFiles->getTempCorrFile();
//
//    bool estimateParamCov = paramCovFile == "";
//    bool estimateTempCorr = (corrTFile    == "" && modelSettings->getUseVerticalVariogram() == false);
//
//    //
//    // Read parameter covariance (Var0) from file or set from output from function generateRockPhysics3DBackground.
//    // Consistency check that only one option (file or rock physics) is possible, is done in XmlModelFile::checkInversionConsistency
//    //
//    float ** paramCorr = NULL;
//    bool failedParamCorr = false;
//    std::string tmpErrText("");
//    if(!estimateParamCov) {
//      paramCorr = ModelAVODynamic::readMatrix(paramCovFile, 3, 3, "parameter covariance", tmpErrText);
//      validateCorrelationMatrix(paramCorr, modelSettings, tmpErrText);
//      if(paramCorr == NULL || tmpErrText != "") {
//        errText += "Reading of file "+paramCovFile+" for parameter covariance matrix failed\n";
//        errText += tmpErrText;
//        failedParamCorr = true;
//      }
//    }
//    else if (modelSettings->getFaciesProbFromRockPhysics() == true) {
//      estimateParamCov = false;
//      paramCorr = new float * [3];
//      for(int i=0;i<3;i++) {
//        paramCorr[i] = new float[3];
//      }
//
//      int n_facies = static_cast<int>(faciesNames_.size());
//
//      std::vector<DistributionsRock *> rock_distribution(n_facies);
//      typedef std::map<std::string, DistributionsRock *> rfMapType;
//      rfMapType rfMap = getRockDistributionTime0();
//
//      for(int i=0; i<n_facies; i++) {
//        rfMapType::iterator iter = rfMap.find(faciesNames_[i]);
//        if(iter != rfMap.end())
//          rock_distribution[i] = iter->second;
//      }
//
//      NRLib::Grid2D<double> param_corr(3,3);
//      calculateCovariancesFromRockPhysics(rock_distribution,
//                                          priorFacies,
//                                          param_corr,
//                                          errText);
//
//
//      for(int i=0; i<3; i++) {
//        for(int j=0; j<3; j++)
//          paramCorr[i][j] = static_cast<float>(param_corr(i,j));
//      }
//
//      if (tmpErrText != "")
//      {
//        errText += "Parameter covariance matrix from rock physics failed\n";
//        errText += tmpErrText;
//        failedParamCorr = true;
//      }
//    }
//
//    //
//    // Estimate lateral correlation from seismic data
//    //
//    priorCorrXY_ = findCorrXYGrid(timeSimbox, modelSettings);
//
//    if(modelSettings->getLateralCorr()==NULL) // NBNB-PAL: this will never be true (default lateral corr)
//    {
//      int timelapse = 0; // Setting timelapse = 0 as this is the generation of prior model
//      estimateCorrXYFromSeismic(priorCorrXY_, seisCube, modelSettings->getNumberOfAngles(timelapse));
//    }
//
//    int nCorrT = timeSimbox->GetNZpad();
//    if((nCorrT % 2) == 0)
//      nCorrT = nCorrT/2+1;
//    else
//      nCorrT = nCorrT/2;
//
//    std::vector<float> corrT;
//
//    bool failedTempCorr = false;
//    if(!estimateTempCorr)
//    {
//      if(modelSettings->getUseVerticalVariogram() == true) {
//        corrT.resize(nCorrT+1);
//        float tempCorrRange = modelSettings->getTemporalCorrelationRange();
//        float dz = static_cast<float>(timeSimbox->getdz());
//        for(int i=0; i<=nCorrT; i++){
//          //using an exponential variogram with a = 1/3 (Chiles and Delfiner 1999)
//          corrT[i] = exp(-3*dz*i/tempCorrRange);
//        }
//      }
//      else{
//        std::string tmpErrText("");
//        float ** corrMat = ModelAVODynamic::readMatrix(corrTFile, 1, nCorrT+1, "temporal correlation", tmpErrText);
//        if(corrMat == NULL)
//        {
//          errText += "Reading of file '"+corrTFile+"' for temporal correlation failed\n";
//          errText += tmpErrText;
//          failedTempCorr = true;
//        }
//        corrT.resize(nCorrT);
//        if (!failedTempCorr)
//        {
//          for(int i=0;i<nCorrT;i++)
//            corrT[i] = corrMat[0][i+1];
//          delete [] corrMat[0];
//          delete [] corrMat;
//        }
//      }
//    }
//
//    float ** pointVar0 = NULL;
//    if (estimateParamCov || estimateTempCorr) //Need well estimation
//    {
//      std::string tmpErrTxt;
//      Analyzelog * analyze = new Analyzelog(wells,
//                                            background,
//                                            timeSimbox,
//                                            modelSettings,
//                                            tmpErrTxt);
//      if (tmpErrTxt != "") {
//        errText += tmpErrTxt;
//        failedParamCorr = true;
//      }
//
//      if(estimateParamCov)
//        paramCorr = analyze->getVar0();
//      else
//        delete [] analyze->getVar0();
//
//      pointVar0 = analyze->getPointVar0();
//
//      float * estCorrT = analyze->getCorrT();
//      if(estimateTempCorr) {
//        corrT.resize(nCorrT);
//        int nEst = analyze->getNumberOfLags();
//        int i, max = nEst;
//        if(max > nCorrT)
//          max = nCorrT;
//        for(i=0;i<max;i++)
//          corrT[i] = estCorrT[i];
//        if(i<nCorrT) {
//          LogKit::LogFormatted(LogKit::High,
//            "\nOnly able to estimate %d of %d lags needed in temporal correlation. The rest are set to 0.\n", nEst, nCorrT);
//          for(;i<nCorrT;i++)
//            corrT[i] = 0.0f;
//        }
//      }
//      delete [] estCorrT;
//
//      delete analyze;
//    }
//
//    if (failedParamCorr || failedTempCorr)
//      failed = true;
//
//    if (!failed) {
//
//      const int nx        = timeSimbox->getnx();
//      const int ny        = timeSimbox->getny();
//      const int nz        = timeSimbox->getnz();
//      const int nxPad     = timeSimbox->GetNXpad();
//      const int nyPad     = timeSimbox->GetNYpad();
//      const int nzPad     = timeSimbox->GetNZpad();
//
//      float dt = static_cast<float>(timeSimbox->getdz());
//      float lowCut = modelSettings->getLowCut();
//      int lowIntCut = int(floor(lowCut*(nzPad*0.001*dt))); // computes the integer whis corresponds to the low cut frequency.
//
//      float corrGradI;
//      float corrGradJ;
//      getCorrGradIJ(corrGradI, corrGradJ);
//      //makeCorr2DPositiveDefinite( priorCorrXY_);
//
//      seismicParameters.setCorrelationParameters(paramCorr,
//                                                 corrT,
//                                                 priorCorrXY_,
//                                                 lowIntCut,
//                                                 corrGradI,
//                                                 corrGradJ,
//                                                 nx,
//                                                 ny,
//                                                 nz,
//                                                 nxPad,
//                                                 nyPad,
//                                                 nzPad);
//
//      for(int i=0; i<3; i++)
//        delete [] paramCorr[i];
//      delete [] paramCorr;
//
//      if(printResult)
//        seismicParameters.writeFilePriorVariances(modelSettings, corrT, priorCorrXY_, dt);
//      seismicParameters.printPriorVariances();
//    }
//
//
//    if(failedTempCorr == true || failedParamCorr == true)
//    {
//      errText += "Could not construct prior covariance. Unknown why...\n";
//      failed = true;
//    }
//
//    Timings::setTimePriorCorrelation(wall,cpu);
//  }
//}

void ModelGeneral::validateCorrelationMatrix(float               ** C,
                                             const ModelSettings *  model_settings,
                                             std::string         &  errTxt)
{
  float minAlpha = model_settings->getVarAlphaMin();
  float maxAlpha = model_settings->getVarAlphaMax();
  float minBeta  = model_settings->getVarBetaMin();
  float maxBeta  = model_settings->getVarBetaMax();
  float minRho   = model_settings->getVarRhoMin();
  float maxRho   = model_settings->getVarRhoMax();

  float C00      = C[0][0];
  float C11      = C[1][1];
  float C22      = C[2][2];
  float C01      = C[0][1];
  float C10      = C[1][0];
  float C02      = C[0][2];
  float C20      = C[2][0];
  float C12      = C[1][2];
  float C21      = C[2][1];

  if (C00 < minAlpha || C00 > maxAlpha) {
    errTxt += "The prior Vp variance is outside valid range:\n";
    errTxt += "  Given value   : " + NRLib::ToString(C00) + "\n";
    errTxt += "  Minimum value : " + NRLib::ToString(minAlpha) + "\n";
    errTxt += "  Maximum value : " + NRLib::ToString(maxAlpha) + "\n";
  }
  if (C11 < minBeta || C11 > maxBeta) {
    errTxt += "The prior Vs variance is outside valid range:\n";
    errTxt += "  Given value   : " + NRLib::ToString(C11) + "\n";
    errTxt += "  Minimum value : " + NRLib::ToString(minBeta) + "\n";
    errTxt += "  Maximum value : " + NRLib::ToString(maxBeta) + "\n";
  }
  if (C22 < minRho || C22 > maxRho) {
    errTxt += "The prior density variance is outside valid range:\n";
    errTxt += "  Given value   : " + NRLib::ToString(C22) + "\n";
    errTxt += "  Minimum value : " + NRLib::ToString(minRho) + "\n";
    errTxt += "  Maximum value : " + NRLib::ToString(maxRho) + "\n";
  }

  float corr01 = C01/(std::sqrt(C00)*std::sqrt(C11));
  float corr02 = C02/(std::sqrt(C00)*std::sqrt(C22));
  float corr12 = C12/(std::sqrt(C11)*std::sqrt(C22));

  if (corr01 < -1.0 || corr01 > 1.0) {
    errTxt += "The prior Vp-Vs correlation is illegal (" + NRLib::ToString(corr01) + ")\n";
  }
  if (corr02 < -1.0 || corr02 > 1.0) {
    errTxt += "The prior Vp-Rho correlation is illegal (" + NRLib::ToString(corr02) + ")\n";
  }
  if (corr12 < -1.0 || corr12 > 1.0) {
    errTxt += "The prior Vs-Rho correlation is illegal (" + NRLib::ToString(corr12) + ")\n";
  }

  if (std::abs(C01 - C10) > 0.0f) {
    errTxt += "The prior covariance matrix is not symmetric in Vp and Vs\n";
    errTxt += "  Corr(Vp,Vs) : " + NRLib::ToString(C01) + "\n";
    errTxt += "  Corr(Vs,Vp) : " + NRLib::ToString(C10) + "\n";
  }
  if (std::abs(C02 - C20) > 0.0f) {
    errTxt += "The prior covariance matrix is not symmetric in Vp and Rho\n";
    errTxt += "  Corr(Vp,Rho) : " + NRLib::ToString(C02) + "\n";
    errTxt += "  Corr(Rho,Vp) : " + NRLib::ToString(C20) + "\n";
  }
  if (std::abs(C12 - C21) > 0.0f) {
    errTxt += "The prior covariance matrix is not symmetric in Vs and Rho\n";
    errTxt += "  Corr(Vs,Rho) : " + NRLib::ToString(C12) + "\n";
    errTxt += "  Corr(Rho,Vs) : " + NRLib::ToString(C21) + "\n";
  }
}

//Surface *
//ModelGeneral::findCorrXYGrid(const Simbox * timeSimbox, const ModelSettings * model_settings)
//{
//  float dx  = static_cast<float>(timeSimbox->getdx());
//  float dy  = static_cast<float>(timeSimbox->getdy());
//
//  int   nx  = model_settings->getNXpad();
//  int   ny  = model_settings->getNYpad();
//
//  Surface * grid = new Surface(0, 0, dx*nx, dy*ny, nx, ny, RMISSING);
//
//  if (model_settings->getLateralCorr()!=NULL) // NBNB-PAL: Denne her blir aldri null etter at jeg la inn en default lateral correlation i modelsettings.
//  {
//    int refi,refj;
//    for (int j=0;j<ny;j++)
//    {
//      for (int i=0;i<nx;i++)
//      {
//        if (i<(nx/2+1))
//        {
//          refi = i;
//        }
//        else
//        {
//          refi = i-nx;
//        }
//        if (j< (ny/2+1))
//        {
//          refj = j;
//        }
//        else
//        {
//          refj = j-ny;
//        }
//        (*grid)(j*nx+i) = model_settings->getLateralCorr()->corr(refi*dx, refj*dy);
//      }
//    }
//  }
//  return(grid);
//}

//  int   nx  = timeSimbox->GetNXpad();
//  int   ny  = timeSimbox->GetNYpad();
//
//  Surface * grid = new Surface(0, 0, dx*nx, dy*ny, nx, ny, RMISSING);
//
//  if(modelSettings->getLateralCorr()!=NULL) // NBNB-PAL: Denne her blir aldri null etter at jeg la inn en default lateral correlation i modelsettings.
//  {
//    int refi,refj;
//    for(int j=0;j<ny;j++)
//    {
//      for(int i=0;i<nx;i++)
//      {
//        if(i<(nx/2+1))
//        {
//          refi = i;
//        }
//        else
//        {
//          refi = i-nx;
//        }
//        if(j< (ny/2+1))
//        {
//          refj = j;
//        }
//        else
//        {
//          refj = j-ny;
//        }
//        (*grid)(j*nx+i) = modelSettings->getLateralCorr()->corr(refi*dx, refj*dy);
//      }
//    }
//  }
//  return(grid);
//}

void
ModelGeneral::estimateCorrXYFromSeismic(Surface *& corrXY,
                                        FFTGrid ** seisCube,
                                        int numberOfAngles)
{
  FFTGrid * transf;
  float   * grid;

  int n = static_cast<int>(corrXY->GetNI()*corrXY->GetNJ());
  grid = new float[n];

  for (int i=0 ; i<n ; i++)
    grid[i] = 0.0;

  for (int i=0 ; i<numberOfAngles ; i++)
  {
    if (seisCube[i]->isFile())
      transf = new FFTFileGrid(static_cast<FFTFileGrid *>(seisCube[i])); //move new out of loop? Copy grid instead
    else
      transf = new FFTGrid(seisCube[i]); //move new out of loop? Copy grid instead

    transf->setAccessMode(FFTGrid::RANDOMACCESS);
    transf->fftInPlace();
    transf->square();
    transf->invFFTInPlace();
    transf->collapseAndAdd( grid ); //the result of the collapse (the result for z=0) is is added to grid
    transf->endAccess();
    delete transf;
  }
  float sill = grid[0];
  for (int i=0;i<n;i++)
    (*corrXY)(i) = grid[i]/sill;
  delete [] grid;
}

void ModelGeneral::checkFaciesNamesConsistency(ModelSettings     *& model_settings,
                                               const InputFiles   * input_files,
                                               std::string        & tmpErrText) const
{
  int nFacies = static_cast<int>(facies_names_.size());

  // Compare names in wells with names given in rock physics prior model
  if (rock_distributions_.size() > 0) {
    int nRocks  = static_cast<int>(rock_distributions_.size());
    if (nRocks > nFacies)
      tmpErrText += "Problem with facies logs. The number of rocks in the rock physics prior model is larger than the number of facies found in the wells.\n";
    for (int i=0; i<nFacies; i++) {
      if (rock_distributions_.find(facies_names_[i]) == rock_distributions_.end())
        tmpErrText += "Problem with facies logs. Facies "+facies_names_[i]+" found in a well is not one of the rocks given in rock physics prior model\n";
    }
  }

  // Compare names in wells with names given in .xml-file
  if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
  {
    typedef std::map<std::string,float> mapType;
    mapType myMap = model_settings->getPriorFaciesProb();

    for (int i=0;i<nFacies;i++)
    {
      mapType::iterator iter = myMap.find(facies_names_[i]);
      if (iter==myMap.end())
        tmpErrText += "Problem with facies logs. Facies "+facies_names_[i]+" is not one of the facies given in the xml-file.\n";
    }
  }

  // Compare names in wells with names given as input in proability cubes
  else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
  {
    typedef std::map<std::string,std::string> mapType;
    mapType myMap = input_files->getPriorFaciesProbFile();

    for (int i=0;i<nFacies;i++)
    {
      mapType::iterator iter = myMap.find(facies_names_[i]);
      if (iter==myMap.end())
        tmpErrText += "Problem with facies logs. Facies "+facies_names_[i]+" is not one of the facies given in the xml-file.\n";
    }
  }
}

//void
//ModelGeneral::setFaciesNamesFromRockPhysics()
//{
//  typedef std::map<std::string, DistributionsRock *> mapType;
//
//  int i = 0;
//  for (std::map<std::string, std::vector<DistributionsRock *> >::const_iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
//    facies_names_.push_back(it->first);
//    faciesLabels_.push_back(i);
//    i++;
//  }
//}

//void
//ModelGeneral::processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
//                                     std::vector<WellData *>        wells,
//                                     Simbox                       * timeSimbox,
//                                     Simbox                       * timeCutSimbox,
//                                     ModelSettings                * model_settings,
//                                     bool                         & failed,
//                                     std::string                  & errTxt,
//                                     const InputFiles             * input_files)
//{
//  if (model_settings->getEstimateFaciesProb() || model_settings->getDo4DInversion())
//  {
//    LogKit::WriteHeader("Prior Facies Probabilities");
//
//    if (facies_names_.size() == 0)
//      setFaciesNamesFromRockPhysics();
//
//    std::string tmpErrText = "";
//    checkFaciesNamesConsistency(model_settings,
//                                input_files,
//                                tmpErrText);
//    if (tmpErrText != "")
//      errTxt += "Prior facies probabilities failed.\n"+tmpErrText;
//
//    int nFacies = static_cast<int>(facies_names_.size());
//
//    if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_WELLS)
//    {
//      if (nFacies > 0)
//      {
//        int   nz      = timeSimbox->getnz();
//        float dz      = static_cast<float>(timeSimbox->getdz());
//        int   nWells  = model_settings->getNumberOfWells();
//        int   ndata   = nWells*nz;
//
//        int ** faciesCount = new int * [nWells];
//        for (int w = 0 ; w < nWells ; w++)
//          faciesCount[w] = new int[nFacies];
//
//        for (int w = 0 ; w < nWells ; w++)
//          for (int i = 0 ; i < nFacies ; i++)
//            faciesCount[w][i] = 0;
//
//        int * faciesLog = new int[ndata];   // NB! *internal* log numbering (0, 1, 2, ...)
//        for (int i = 0 ; i < ndata ; i++)
//          faciesLog[i] = IMISSING;
//
//        float * vtAlpha   = new float[nz];  // vt = vertical trend
//        float * vtBeta    = new float[nz];
//        float * vtRho     = new float[nz];
//        int   * vtFacies  = new int[nz];
//
//        int nUsedWells = 0;
//
//        for (int w = 0 ; w < nWells ; w++)
//        {
//          if (wells[w]->getNFacies() > 0) // Well has facies log
//          {
//            //
//            // Note that we use timeSimbox to calculate prior facies probabilities
//            // instead of the simbox with parallel top and base surfaces. This
//            // will make the prior probabilities slightly different, but that
//            // should not be a problem.
//            //
//            BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
//            int nBlocks = bl->getNumberOfBlocks();
//            //
//            // Set facies data outside facies estimation interval IMISSING
//            //
//            int * blFaciesLog = new int[nBlocks];
//            Utils::copyVector(bl->getFacies(), blFaciesLog, nBlocks);
//
//            if (faciesEstimInterval.size() > 0) {
//              const double * xPos  = bl->getXpos();
//              const double * yPos  = bl->getYpos();
//              const double * zPos  = bl->getZpos();
//              for (int i = 0 ; i < nBlocks ; i++) {
//                const double zTop  = faciesEstimInterval[0]->GetZ(xPos[i],yPos[i]);
//                const double zBase = faciesEstimInterval[1]->GetZ(xPos[i],yPos[i]);
//                if ( (zPos[i] - 0.5*dz) < zTop || (zPos[i] + 0.5*dz) > zBase)
//                  blFaciesLog[i] = IMISSING;
//              }
//            }
//
//            bl->getVerticalTrend(bl->getAlpha(),vtAlpha);
//            bl->getVerticalTrend(bl->getBeta(),vtBeta);
//            bl->getVerticalTrend(bl->getRho(),vtRho);
//            bl->getVerticalTrend(blFaciesLog,vtFacies);
//            delete [] blFaciesLog;
//
//            for (int i=0 ; i<nz ; i++)
//            {
//              int facies;
//              if (vtAlpha[i] != RMISSING && vtBeta[i] != RMISSING && vtRho[i] != RMISSING)
//                facies = vtFacies[i];
//              else
//                facies = IMISSING;
//
//              faciesLog[w*nz + i] = facies;
//              if (facies != IMISSING)
//                faciesCount[w][facies]++;
//            }
//            nUsedWells++;
//          }
//        }
//        delete [] vtAlpha;
//        delete [] vtBeta;
//        delete [] vtRho;
//        delete [] vtFacies;
//
//        if (nUsedWells > 0) {
//          //
//          // Probabilities
//          //
//          LogKit::LogFormatted(LogKit::Low,"\nFacies distributions for each blocked well: \n");
//          LogKit::LogFormatted(LogKit::Low,"\nBlockedWell              ");
//          for (int i = 0 ; i < nFacies ; i++)
//            LogKit::LogFormatted(LogKit::Low,"%12s ",facies_names_[i].c_str());
//          LogKit::LogFormatted(LogKit::Low,"\n");
//          for (int i = 0 ; i < 24+13*nFacies ; i++)
//            LogKit::LogFormatted(LogKit::Low,"-");
//          LogKit::LogFormatted(LogKit::Low,"\n");
//          for (int w = 0 ; w < nWells ; w++)
//          {
//            if (wells[w]->getNFacies() > 0) // Well has facies log
//            {
//              float tot = 0.0;
//              for (int i = 0 ; i < nFacies ; i++) {
//                tot += static_cast<float>(faciesCount[w][i]);
//              }
//
//              LogKit::LogFormatted(LogKit::Low,"%-23s ",wells[w]->getWellname().c_str());
//              for (int i = 0 ; i < nFacies ; i++) {
//                float faciesProb = static_cast<float>(faciesCount[w][i])/tot;
//                LogKit::LogFormatted(LogKit::Low," %12.4f",faciesProb);
//              }
//              LogKit::LogFormatted(LogKit::Low,"\n");
//            }
//          }
//          LogKit::LogFormatted(LogKit::Low,"\n");
//          //
//          // Counts
//          //
//          LogKit::LogFormatted(LogKit::Medium,"\nFacies counts for each blocked well: \n");
//
//          LogKit::LogFormatted(LogKit::Medium,"\nBlockedWell              ");
//          for (int i = 0 ; i < nFacies ; i++)
//            LogKit::LogFormatted(LogKit::Medium,"%12s ",facies_names_[i].c_str());
//          LogKit::LogFormatted(LogKit::Medium,"\n");
//          for (int i = 0 ; i < 24+13*nFacies ; i++)
//            LogKit::LogFormatted(LogKit::Medium,"-");
//          LogKit::LogFormatted(LogKit::Medium,"\n");
//          for (int w = 0 ; w < nWells ; w++)
//          {
//            if (wells[w]->getNFacies() > 0)
//            {
//              float tot = 0.0;
//              for (int i = 0 ; i < nFacies ; i++)
//                tot += static_cast<float>(faciesCount[w][i]);
//              LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells[w]->getWellname().c_str());
//              for (int i = 0 ; i < nFacies ; i++) {
//                LogKit::LogFormatted(LogKit::Medium," %12d",faciesCount[w][i]);
//              }
//              LogKit::LogFormatted(LogKit::Medium,"\n");
//            }
//          }
//          LogKit::LogFormatted(LogKit::Medium,"\n");
//
//          for (int w = 0 ; w < nWells ; w++)
//            delete [] faciesCount[w];
//          delete [] faciesCount;
//
//          //
//          // Make prior facies probabilities
//          //
//          float sum = 0.0f;
//          int * nData = new int[nFacies];
//          for (int i=0 ; i<nFacies ; i++)
//            nData[i] = 0;
//
//          for (int i=0 ; i<ndata ; i++) {
//            if (faciesLog[i] != IMISSING) {
//              nData[faciesLog[i]]++;
//            }
//          }
//          delete [] faciesLog;
//
//          for (int i=0 ; i<nFacies ; i++)
//            sum += nData[i];
//
//          if (sum > 0) {
//            LogKit::LogFormatted(LogKit::Low,"Facies probabilities based on all blocked wells:\n\n");
//            LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
//            LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
//            priorFacies_.resize(nFacies);
//            for (int i=0 ; i<nFacies ; i++) {
//              priorFacies_[i] = float(nData[i])/sum;
//              LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",facies_names_[i].c_str(),priorFacies_[i]);
//            }
//          }
//          else {
//            LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No valid facies log entries have been found\n");
//            model_settings->setEstimateFaciesProb(false);
//            TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
//
//          }
//          delete [] nData;
//        }
//        else
//        {
//          LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Estimation of facies probabilites have been requested, but there");
//          LogKit::LogFormatted(LogKit::Warning,"\n         are no wells with facies available and CRAVA will therefore not");
//          LogKit::LogFormatted(LogKit::Warning,"\n         be able to estimate these probabilities...\n");
//          model_settings->setEstimateFaciesProb(false);
//
//          TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
//        }
//      }
//      else
//      {
//        LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Estimation of facies probabilites have been requested, but no facies");
//        LogKit::LogFormatted(LogKit::Warning,"\n         have been found and CRAVA will therefore not be able to estimate");
//        LogKit::LogFormatted(LogKit::Warning,"\n         these probabilities...\n");
//        model_settings->setEstimateFaciesProb(false);
//        TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
//      }
//    }
//    else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
//    {
//      priorFacies_.resize(nFacies);
//      typedef std::map<std::string,float> mapType;
//      mapType myMap = model_settings->getPriorFaciesProb();
//
//      for (int i=0;i<nFacies;i++)
//      {
//        mapType::iterator iter = myMap.find(facies_names_[i]);
//        if (iter!=myMap.end())
//          priorFacies_[i] = iter->second;
//        else
//        {
//          LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",facies_names_[i].c_str());
//          model_settings->setEstimateFaciesProb(false);
//          TaskList::addTask("Check that facies " +NRLib::ToString(facies_names_[i].c_str())+" is given a prior probability in the xml-file");
//        }
//      }
//      LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
//      LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
//      for (int i=0 ; i<nFacies ; i++) {
//        LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",facies_names_[i].c_str(),priorFacies_[i]);
//      }
//
//    }
//    else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
//    {
//      readPriorFaciesProbCubes(input_files,
//                               model_settings,
//                               priorFaciesProbCubes_,
//                               timeSimbox,
//                               timeCutSimbox,
//                               errTxt,
//                               failed);
//
//       typedef std::map<std::string,std::string> mapType;
//       mapType myMap = input_files->getPriorFaciesProbFile();
//
//       LogKit::LogFormatted(LogKit::Low,"Facies         Probability in file\n");
//       LogKit::LogFormatted(LogKit::Low,"----------------------------------\n");
//       for (mapType::iterator it=myMap.begin();it!=myMap.end();it++)
//         LogKit::LogFormatted(LogKit::Low,"%-15s %10s\n",(it->first).c_str(),(it->second).c_str());
//
//    }
//  }
//  if (errTxt != "")
//    failed = true;
//}

//void
//ModelGeneral::readPriorFaciesProbCubes(const InputFiles        * input_files,
//                                       ModelSettings           * model_settings,
//                                       std::vector<FFTGrid *>  & priorFaciesProbCubes,
//                                       Simbox                  * timeSimbox,
//                                       Simbox                  * timeCutSimbox,
//                                       std::string             & errTxt,
//                                       bool                    & failed)
//{
//  int nFacies = static_cast<int>(facies_names_.size());
//  priorFaciesProbCubes.resize(nFacies);
//
//  typedef std::map<std::string,std::string> mapType;
//  mapType myMap = input_files->getPriorFaciesProbFile();
//  for (int i=0;i<nFacies;i++)
//  {
//    mapType::iterator iter = myMap.find(facies_names_[i]);
//
//    if (iter!=myMap.end())
//    {
//      const std::string & faciesProbFile = iter->second;
//      const SegyGeometry      * dummy1 = NULL;
//      const TraceHeaderFormat * dummy2 = NULL;
//      const float               offset = model_settings->getSegyOffset(0); //Facies estimation only allowed for one time lapse
//      std::string errorText("");
//      ModelGeneral::readGridFromFile(faciesProbFile,
//                                     "priorfaciesprob",
//                                     offset,
//                                     priorFaciesProbCubes[i],
//                                     dummy1,
//                                     dummy2,
//                                     FFTGrid::PARAMETER,
//                                     timeSimbox,
//                                     timeCutSimbox,
//                                     model_settings,
//                                     errorText,
//                                     true);
//      if (errorText != "")
//      {
//        errorText += "Reading of file \'"+faciesProbFile+"\' for prior facies probability for facies \'"
//                     +facies_names_[i]+"\' failed\n";
//        errTxt += errorText;
//        failed = true;
//      }
//    }
//    else
//    {
//      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",
//                           facies_names_[i].c_str());
//      TaskList::addTask("Check that facies "+NRLib::ToString(facies_names_[i].c_str())+" is given prior probability in the xml-file");
//      model_settings->setEstimateFaciesProb(false);
//      break;
//    }
//  }
//}

//bool
//ModelGeneral::process4DBackground(ModelSettings           *& model_settings,
//                                  const InputFiles         * input_files,
//                                  SeismicParametersHolder  & seismicParameters,
//                                  std::string              & errText,
//                                  bool                     & failed,
//                                  NRLib::Vector            & initialMean,
//                                  NRLib::Matrix            & initialCov)
//{
//  bool failedPriorCorr               = false;
//  bool failedRockPhysics4DBackground = false;
//
//  // Variables to be created in this function
//  Background * background = NULL;
//
//  FFTGrid **seisCube = NULL;  // vet ikke hva denne skal være
//
//  int n_facies = static_cast<int>(facies_names_.size());
//
//  std::vector<float> priorProbability(n_facies);
//  typedef std::map<std::string,float> mapType;
//  mapType myMap = model_settings->getPriorFaciesProb();
//
//  for (int i=0;i<n_facies;i++) {
//    mapType::iterator iter = myMap.find(facies_names_[i]);
//    if (iter!=myMap.end())
//      priorProbability[i] = static_cast<float>(iter->second);
//  }
//
//  std::vector<DistributionsRock *> rock_distribution(n_facies);
//  typedef std::map<std::string, DistributionsRock *> rfMapType;
//  rfMapType rfMap = getRockDistributionTime0();
//
//  for (int i=0; i<n_facies; i++) {
//    rfMapType::iterator iter = rfMap.find(facies_names_[i]);
//    if (iter != rfMap.end())
//      rock_distribution[i] = iter->second;
//  }
//
//  std::string tmpError = "";
//  setUp3DPartOf4DBackground(rock_distribution,
//                            priorProbability,
//                            simbox_,
//                            *model_settings,
//                            seismicParameters,
//                            state4d_,
//                            tmpError);
//
//  processPriorCorrelations(background,
//                           wells_,
//                           simbox_,
//                           model_settings,
//                           priorProbability,
//                           seisCube,
//                           input_files,
//                           seismicParameters,
//                           errText,
//                           failedPriorCorr);
//
//
//  copyCorrelationsTo4DState( seismicParameters, state4d_);
//  const int nx    = simbox_->getnx();
//  const int ny    = simbox_->getny();
//  const int nz    = simbox_->getnz();
//  const int nxPad = model_settings->getNXpad();
//  const int nyPad = model_settings->getNYpad();
//  const int nzPad = model_settings->getNZpad();
//
//  complete4DBackground(nx, ny, nz, nxPad, nyPad, nzPad,initialMean,initialCov);
//
//  if (tmpError != "") {
//    errText += tmpError;
//    failedRockPhysics4DBackground = true;
//  }
//
//  failed = failedPriorCorr || failedRockPhysics4DBackground;
//
//  return failed;
//
//}

//void
//ModelGeneral::readPriorFaciesProbCubes(const InputFiles        * inputFiles,
//                                       ModelSettings           * modelSettings,
//                                       std::vector<FFTGrid *>  & priorFaciesProbCubes,
//                                       Simbox                  * timeSimbox,
//                                       Simbox                  * timeCutSimbox,
//                                       std::string             & errTxt,
//                                       bool                    & failed)
//{
//  int nFacies = static_cast<int>(faciesNames_.size());
//  priorFaciesProbCubes.resize(nFacies);
//
//  typedef std::map<std::string,std::string> mapType;
//  mapType myMap = inputFiles->getPriorFaciesProbFile();
//  for(int i=0;i<nFacies;i++)
//  {
//    mapType::iterator iter = myMap.find(faciesNames_[i]);
//
//    if(iter!=myMap.end())
//    {
//      const std::string & faciesProbFile = iter->second;
//      const SegyGeometry      * dummy1 = NULL;
//      const TraceHeaderFormat * dummy2 = NULL;
//      const float               offset = modelSettings->getSegyOffset(0); //Facies estimation only allowed for one time lapse
//      std::string errorText("");
//      ModelGeneral::readGridFromFile(faciesProbFile,
//                                     "priorfaciesprob",
//                                     offset,
//                                     priorFaciesProbCubes[i],
//                                     dummy1,
//                                     dummy2,
//                                     FFTGrid::PARAMETER,
//                                     timeSimbox,
//                                     timeCutSimbox,
//                                     modelSettings,
//                                     errorText,
//                                     true);
//      if(errorText != "")
//      {
//        errorText += "Reading of file \'"+faciesProbFile+"\' for prior facies probability for facies \'"
//                     +faciesNames_[i]+"\' failed\n";
//        errTxt += errorText;
//        failed = true;
//      }
//    }
//    else
//    {
//      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",
//                           faciesNames_[i].c_str());
//      TaskList::addTask("Check that facies "+NRLib::ToString(faciesNames_[i].c_str())+" is given prior probability in the xml-file");
//      modelSettings->setEstimateFaciesProb(false);
//      break;
//    }
//  }
//}
//
//bool
//ModelGeneral::process4DBackground(ModelSettings           *& modelSettings,
//                                  const InputFiles         * inputFiles,
//                                  SeismicParametersHolder  & seismicParameters,
//                                  std::string              & errText,
//                                  bool                     & failed,
//                                  NRLib::Vector            & initialMean,
//                                  NRLib::Matrix            & initialCov)
//{
//  bool failedPriorCorr               = false;
//  bool failedRockPhysics4DBackground = false;
//
//  // Variables to be created in this function
//  Background * background = NULL;
//
//  FFTGrid **seisCube = NULL;  // vet ikke hva denne skal være
//
//  int n_facies = static_cast<int>(faciesNames_.size());
//
//  std::vector<float> priorProbability(n_facies);
//  typedef std::map<std::string,float> mapType;
//  mapType myMap = modelSettings->getPriorFaciesProb();
//
//  for(int i=0;i<n_facies;i++) {
//    mapType::iterator iter = myMap.find(faciesNames_[i]);
//    if(iter!=myMap.end())
//      priorProbability[i] = static_cast<float>(iter->second);
//  }
//
//  std::vector<DistributionsRock *> rock_distribution(n_facies);
//  typedef std::map<std::string, DistributionsRock *> rfMapType;
//  rfMapType rfMap = getRockDistributionTime0();
//
//  for(int i=0; i<n_facies; i++) {
//    rfMapType::iterator iter = rfMap.find(faciesNames_[i]);
//    if(iter != rfMap.end())
//      rock_distribution[i] = iter->second;
//  }
//
//  std::string tmpError = "";
//  setUp3DPartOf4DBackground(rock_distribution,
//                            priorProbability,
//                            timeSimbox_,
//                            *modelSettings,
//                            seismicParameters,
//                            state4d_,
//                            tmpError);
//
//  processPriorCorrelations(background,
//                           wells_,
//                           timeSimbox_,
//                           modelSettings,
//                           priorProbability,
//                           seisCube,
//                           inputFiles,
//                           seismicParameters,
//                           errText,
//                           failedPriorCorr);
//
//
//  copyCorrelationsTo4DState( seismicParameters, state4d_);
//  const int nx    = timeSimbox_->getnx();
//  const int ny    = timeSimbox_->getny();
//  const int nz    = timeSimbox_->getnz();
//  const int nxPad = timeSimbox_->GetNXpad();
//  const int nyPad = timeSimbox_->GetNYpad();
//  const int nzPad = timeSimbox_->GetNZpad();
//
//  complete4DBackground(nx, ny, nz, nxPad, nyPad, nzPad,initialMean,initialCov);
//
//  if(tmpError != "") {
//    errText += tmpError;
//    failedRockPhysics4DBackground = true;
//  }
//
//  failed = failedPriorCorr || failedRockPhysics4DBackground;
//
//  return failed;
//
//}

void
ModelGeneral::SetupState4D(ModelSettings           *& modelSettings,
                           SeismicParametersHolder  & seismicParameters,
                           NRLib::Vector            & initialMean,
                           NRLib::Matrix            & initialCov)
{
  //H Difference: Earlier a background was made from rockphysics3d, which was copied to seismicParameters and state4d.
  //Now, use background created in CommonData
  state4d_.setStaticMu(seismicParameters.GetMeanVp(), seismicParameters.GetMeanVs(), seismicParameters.GetMeanRho());

  CopyCorrelationsTo4DState(seismicParameters, state4d_);

  const int nx    = simbox_->getnx();
  const int ny    = simbox_->getny();
  const int nz    = simbox_->getnz();
  const int nxPad = simbox_->GetNXpad();
  const int nyPad = simbox_->GetNYpad();
  const int nzPad = simbox_->GetNZpad();

  Complete4DBackground(nx, ny, nz, nxPad, nyPad, nzPad, initialMean, initialCov);
}

void
ModelGeneral::Complete4DBackground(const int nx, const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad,NRLib::Vector &initial_mean,NRLib::Matrix &initial_cov)
{
  // Static grids (3 + 6) are set in process4DBackground.
  // Dynamic grids (3 + 6 + 9) are set here.

  FFTGrid * dynamicVp;
  FFTGrid * dynamicVs;
  FFTGrid * dynamicRho;
  FFTGrid * dynamicVpVp;

  FFTGrid *dynamicVpVs;
  FFTGrid *dynamicVpRho;
  FFTGrid *dynamicVsVs;
  FFTGrid *dynamicVsRho;
  FFTGrid *dynamicRhoRho;

  FFTGrid *staticDynamicVpVp;
  FFTGrid *staticDynamicVpVs;
  FFTGrid *staticDynamicVpRho;
  FFTGrid *staticDynamicVsVp;
  FFTGrid *staticDynamicVsVs;
  FFTGrid *staticDynamicVsRho;
  FFTGrid *staticDynamicRhoVp;
  FFTGrid *staticDynamicRhoVs;
  FFTGrid *staticDynamicRhoRho;

  dynamicVp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVp->fillInConstant(0.0);
  dynamicVp->setType(FFTGrid::PARAMETER);
  dynamicVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVs->fillInConstant(0.0);
  dynamicVs->setType(FFTGrid::PARAMETER);
  dynamicRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicRho->fillInConstant(0.0);
  dynamicRho->setType(FFTGrid::PARAMETER);

  state4d_.setDynamicMu(dynamicVp, dynamicVs, dynamicRho);
  initial_mean=state4d_.GetFullMean000();

  dynamicVpVp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVpVp->fillInConstant(0.0);
  dynamicVpVp->setType(FFTGrid::COVARIANCE);
  dynamicVpVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVpVs->fillInConstant(0.0);
  dynamicVpVs->setType(FFTGrid::COVARIANCE);
  dynamicVpRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVpRho->fillInConstant(0.0);
  dynamicVpRho->setType(FFTGrid::COVARIANCE);
  dynamicVsVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVsVs->fillInConstant(0.0);
  dynamicVsVs->setType(FFTGrid::COVARIANCE);
  dynamicVsRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVsRho->fillInConstant(0.0);
  dynamicVsRho->setType(FFTGrid::COVARIANCE);
  dynamicRhoRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicRhoRho->fillInConstant(0.0);
  dynamicRhoRho->setType(FFTGrid::COVARIANCE);

  state4d_.setDynamicSigma(dynamicVpVp, dynamicVpVs, dynamicVpRho,
                                        dynamicVsVs, dynamicVsRho,
                                                     dynamicRhoRho);

  staticDynamicVpVp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVpVp->fillInConstant(0.0);
  staticDynamicVpVp->setType(FFTGrid::COVARIANCE);
  staticDynamicVpVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVpVs->fillInConstant(0.0);
  staticDynamicVpVs->setType(FFTGrid::COVARIANCE);
  staticDynamicVpRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVpRho->fillInConstant(0.0);
  staticDynamicVpRho->setType(FFTGrid::COVARIANCE);
  staticDynamicVsVp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVsVp->fillInConstant(0.0);
  staticDynamicVsVp->setType(FFTGrid::COVARIANCE);
  staticDynamicVsVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVsVs->fillInConstant(0.0);
  staticDynamicVsVs->setType(FFTGrid::COVARIANCE);
  staticDynamicVsRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVsRho->fillInConstant(0.0);
  staticDynamicVsRho->setType(FFTGrid::COVARIANCE);
  staticDynamicRhoVp = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicRhoVp->fillInConstant(0.0);
  staticDynamicRhoVp->setType(FFTGrid::COVARIANCE);
  staticDynamicRhoVs = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicRhoVs->fillInConstant(0.0);
  staticDynamicRhoVs->setType(FFTGrid::COVARIANCE);
  staticDynamicRhoRho = ModelGeneral::createFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicRhoRho->fillInConstant(0.0);
  staticDynamicRhoRho->setType(FFTGrid::COVARIANCE);

  state4d_.setStaticDynamicSigma(staticDynamicVpVp,  staticDynamicVpVs,  staticDynamicVpRho,
                                 staticDynamicVsVp,  staticDynamicVsVs,  staticDynamicVsRho,
                                 staticDynamicRhoVp, staticDynamicRhoVs, staticDynamicRhoRho);

  initial_cov=state4d_.GetFullCov();

  state4d_.FFT();
}

void
ModelGeneral::advanceTime(int time_step, SeismicParametersHolder & seismicParameters,ModelSettings* model_settings)
{
  bool debug=false;
  if (debug) dump4Dparameters(model_settings, "_prior", time_step);  // note this prior should be equal to
                                                                    // next_prior in previous step
  if (debug) dumpSeismicParameters(model_settings,"_posterior", time_step,seismicParameters);
  state4d_.split(seismicParameters);
  if (debug) dump4Dparameters(model_settings, "_posterior", time_step);
  state4d_.evolve(time_step, time_evolution_); //NBNB grad I grad J
  //if (debug) dump4Dparameters(model_settings, "_next_prior", time_step+1);
  state4d_.merge(seismicParameters);
  if (debug) dumpSeismicParameters(model_settings,"_next_prior", time_step+1,seismicParameters);
  seismicParameters.invFFTAllGrids(); //merge gives FFT-transformed version, need the standard for now.
}


void
ModelGeneral::lastUpdateOfStaticAndDynamicParts(SeismicParametersHolder &  seismicParameters,ModelSettings* model_settings)
{
  bool debug=true;
  int time_step=time_evolution_.GetNTimSteps()-1;
  if (debug) dumpSeismicParameters(model_settings,"_posterior", time_step,seismicParameters);

  state4d_.split(seismicParameters);
  dump4Dparameters(model_settings, "_posterior", time_step);

}

bool
ModelGeneral::do4DRockPhysicsInversion(ModelSettings* model_settings)
{

  std::vector<FFTGrid*> predictions = state4d_.doRockPhysicsInversion(*time_line_, rock_distributions_.begin()->second,  time_evolution_);
  int nParamOut =predictions.size();

  std::vector<std::string> labels(nParamOut);

  int i=0;

  for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator it = reservoir_variables_.begin(); it != reservoir_variables_.end(); it++)
  {
    labels[i] = it->first;
    i++;
  }

  std::string  outPre =  "mu_";

  for (int i=0;i<nParamOut;i++)
  {
     std::string fileName;
     fileName= outPre + labels[i];
     ParameterOutput::writeToFile(simbox_,this, model_settings, predictions[i] , fileName, labels[i]);
  }

  return 0;
}


void
ModelGeneral::dumpSeismicParameters(ModelSettings* model_settings, std::string identifyer, int timestep,SeismicParametersHolder &  current_state)
{

  std::string  label;
  std::string fileName;
  std::stringstream tag;
  bool transformHere=false;

  if (current_state.GetMeanVp()->getIsTransformed())
  {
    transformHere=true;
    current_state.invFFTAllGrids();
  }

  // write mu current
  tag.str(std::string());tag.clear();label = "mean_vp_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings,  current_state.GetMeanVp() , fileName,  tag.str(),true);
  /*
  tag.str(std::string());tag.clear();label = "mean_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings,  current_state.GetMuBeta(), fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetMuRho() , fileName, tag.str() ,true);
  // */
  // write sigma current
  tag.str(std::string());tag.clear();label = "cov_vp_vp_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCovVp() , fileName,  tag.str(),true);
  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCrCovAlphaBeta() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCrCovAlphaRho() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCovBeta() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_current_step_"; tag << label << timestep << identifyer ; fileName= tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCrCovBetaRho() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCovRho() , fileName,  tag.str(),true);
  // */
  if (transformHere)
    current_state.FFTAllGrids();
}

void
ModelGeneral::dump4Dparameters(ModelSettings* model_settings, std::string identifyer, int timestep)
{
  state4d_.iFFT();

  std::string  outPath =  "";
  std::string  label;
  std::string fileName;
  std::stringstream tag;

  // write mu static
  tag.str(std::string());tag.clear();label = "mean_vp_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuVpStatic() , fileName,  tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "mean_vs_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuVsStatic() , fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuRhoStatic() , fileName, tag.str(),true);
  // */
  // write mu dynamic
  tag.str(std::string());tag.clear();label = "mean_vp_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuVpDynamic() , fileName, tag.str(),true);
  /*
  tag.str(std::string());tag.clear();label = "mean_vs_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuVsDynamic() , fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuRhoDynamic() , fileName,  tag.str(),true);
  // */


  // write sigma static - static
  tag.str(std::string());tag.clear();label = "cov_vp_vp_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVpStaticStatic() , fileName,  tag.str(),true);
  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVsStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpRhoStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVsStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsRhoStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoRhoStaticStatic() , fileName,  tag.str(),true);
     // */
  // write sigma dynamic - dynamic
  tag.str(std::string());tag.clear();label = "cov_vp_vp_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVpDynamicDynamic() , fileName,  tag.str(),true);
  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVsDynamicDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpRhoDynamicDynamic(), fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVsDynamicDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsRhoDynamicDynamic(), fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoRhoDynamicDynamic() , fileName,  tag.str(),true);
  // */
  // write sigma static - dynamic
  tag.str(std::string());tag.clear();label = "cov_vp_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVpStaticDynamic() , fileName,  tag.str(),true);
  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVsStaticDynamic() , fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpRhoStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVpStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVsStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsRhoStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoVpStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoVsStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoRhoStaticDynamic() , fileName,  tag.str(),true);
   // */
  state4d_.FFT();
}

void
ModelGeneral::makeCorr2DPositiveDefinite(Surface         * corrXY)
{
  int      nxp    = corrXY->GetNI();
  int      nyp    = corrXY->GetNJ();
  FFTGrid  helper = FFTGrid(nxp,nyp,1,nxp,nyp,1);
  helper.createRealGrid();
  helper.setType(FFTGrid::COVARIANCE);

  for (int i =0;i<nxp;i++)
    for (int j =0;j<nyp;j++)
    {
      float value = float((*corrXY)(i+j*nxp));
      helper.setRealValue(i,j,0,value);
    }

  helper.fftInPlace();
  int cnxp =helper.getCNxp();
  for (int i =0;i<cnxp;i++)
    for (int j =0;j<nyp;j++)
    {
      fftw_complex value;
      value=helper.getComplexValue(i,j,0);
      value.re=std::sqrt(value.re*value.re+value.im*value.im);
      value.im=0.0f;
      helper.setComplexValue(i,j,0,value);
    }

  helper.invFFTInPlace();
  double scale=1.0/double(helper.getRealValue(0,0,0));

  printf("\nFix in latteral correlation in CRAVA results in a variance increase of %f %% (of 100%%) \n",(scale-1.0)*100);

  for (int i =0;i<nxp;i++)
    for (int j =0;j<nyp;j++)
       (*corrXY)(i+j*nxp)=helper.getRealValue(i,j,0)*scale;
}
