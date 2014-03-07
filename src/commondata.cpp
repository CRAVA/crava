/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>
#define _USE_MATH_DEFINES

#include "src/commondata.h"
#include "src/analyzelog.h"
#include "src/simbox.h"
#include "src/timeline.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/seismicstorage.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/multiintervalgrid.h"
#include "src/cravatrend.h"
#include "src/background.h"
#include "src/blockedlogscommon.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/segy/segytrace.hpp"
#include "nrlib/well/norsarwell.hpp"

#include "rplib/distributionsrock.h"

#include "lib/timekit.hpp"
#include "src/timings.h"

CommonData::CommonData(ModelSettings * model_settings,
                       InputFiles    * input_files):
  outer_temp_simbox_(false),
  read_seismic_(false),
  read_wells_(false),
  block_wells_(false),
  setup_reflection_matrix_(false),
  temporary_wavelet_(false),
  optimize_well_location_(false),
  wavelet_handling_(false),
  setup_multigrid_(false),
  setup_trend_cubes_(false),
  setup_estimation_rock_physics_(false),
  setup_prior_facies_probabilities_(false),
  setup_depth_conversion_(false),
  setup_background_model_(false),
  setup_prior_correlation_(false),
  setup_timeline_(false),
  //prior_corr_estimation_(false),
  setup_gravity_inversion_(false),
  setup_traveltime_inversion_(false),
  //trend_cubes_(false),
  refmat_from_file_global_vpvs_(false),
  velocity_from_inversion_(false),
  multiple_interval_grid_(NULL),
  time_line_(NULL),
  estimation_simbox_()
{

  LogKit::WriteHeader("Common Data");
  std::string err_text = "";

  forward_modeling_ = model_settings->getForwardModeling();

  // 1. set up outer simbox
  outer_temp_simbox_ = CreateOuterTemporarySimbox(model_settings, input_files, estimation_simbox_, full_inversion_volume_, err_text);

  // 2. read seismic data
  read_seismic_ = ReadSeismicData(model_settings, input_files, err_text);

  // 3. read well data
  read_wells_ = ReadWellData(model_settings, &estimation_simbox_, input_files, log_names_, model_settings->getLogNames(),
                             model_settings->getInverseVelocity(), model_settings->getFaciesLogGiven(), err_text);

  // 8. Setup of multiple interval grid
  bool multi_failed = false;
  multiple_interval_grid_ = new MultiIntervalGrid(model_settings, input_files, &estimation_simbox_, err_text, multi_failed);
  setup_multigrid_ = !multi_failed;

  // 4. block wells for estimation
  // if well position is to be optimized or
  // if wavelet/noise should be estimated or
  // if correlations should be estimated
  if (model_settings->getOptimizeWellLocation() || model_settings->getEstimateWaveletNoise() || model_settings->getEstimateCorrelations())
    block_wells_ = BlockWellsForEstimation(model_settings, estimation_simbox_, multiple_interval_grid_,
                                            wells_, mapped_blocked_logs_, mapped_blocked_logs_for_correlation_, err_text);

  // 5. Reflection matrix and wavelet
  setup_reflection_matrix_ = SetupReflectionMatrix(model_settings, input_files, err_text);
  if (model_settings->getOptimizeWellLocation() && read_seismic_ && setup_reflection_matrix_)
    temporary_wavelet_ = SetupTemporaryWavelet(model_settings, err_text);

  // 6. Optimization of well location
  if (model_settings->getOptimizeWellLocation() && read_seismic_ && read_wells_ && setup_reflection_matrix_ && temporary_wavelet_)
    optimize_well_location_ = OptimizeWellLocations(model_settings, input_files, &estimation_simbox_, wells_, mapped_blocked_logs_, seismic_data_, reflection_matrix_, err_text);

  // 7. Wavelet Handling
  wavelet_handling_ = WaveletHandling(model_settings, input_files, err_text);

  // 9. Trend Cubes
  if (setup_multigrid_ && model_settings->getFaciesProbFromRockPhysics() && model_settings->getTrendCubeParameters().size() > 0) {
    setup_trend_cubes_ = SetupTrendCubes(model_settings, input_files, multiple_interval_grid_, err_text);
  }

  // 10. Rock Physics
  if (read_wells_ && setup_multigrid_ && model_settings->getFaciesProbFromRockPhysics()) {
    if (model_settings->getTrendCubeParameters().size() > 0) { // If trends are used, the setup of trend cubes must be ok as well
      if (setup_trend_cubes_) {
        setup_estimation_rock_physics_ = SetupRockPhysics(model_settings, input_files, multiple_interval_grid_, multiple_interval_grid_->GetTrendCubes(),
                                                          mapped_blocked_logs_, n_trend_cubes_, err_text);
      }
    }
    else {
        setup_estimation_rock_physics_ = SetupRockPhysics(model_settings, input_files, multiple_interval_grid_, multiple_interval_grid_->GetTrendCubes(),
                                                          mapped_blocked_logs_, n_trend_cubes_, err_text);
    }
  }

  // 11. Setup of prior facies probabilities
  if (setup_multigrid_) {
    if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_WELLS) {
      if (read_wells_)
        setup_prior_facies_probabilities_ = SetupPriorFaciesProb(model_settings, input_files, err_text);
    }
    else
      setup_prior_facies_probabilities_ = SetupPriorFaciesProb(model_settings, input_files, err_text);
  }

  // 12. Set up background model
  if (setup_multigrid_) {
    if (model_settings->getGenerateBackground() || model_settings->getEstimateBackground()) {
      if (model_settings->getGenerateBackgroundFromRockPhysics() == false && read_wells_)
        setup_background_model_ = SetupBackgroundModel(model_settings, input_files, err_text);
      else if (model_settings->getGenerateBackgroundFromRockPhysics() == true && setup_estimation_rock_physics_)
        setup_background_model_ = SetupBackgroundModel(model_settings, input_files, err_text);
    }
    else //Not estimation
      setup_background_model_ = SetupBackgroundModel(model_settings, input_files, err_text);
  }

  // 13. Setup of prior correlation
   //H-DEBUGGING
  //if(read_seismic_){
  //  if(model_settings->getEstimateCorrelations() == true){
  //    if(read_wells_ && setup_multigrid_){
  //        setup_prior_correlation_ = SetupPriorCorrelation(model_settings, input_files, wells_, mapped_blocked_logs_for_correlation_, multiple_interval_grid_->GetIntervalSimboxes(),
  //          prior_facies_, facies_names_, multiple_interval_grid_->GetTrendCubes(), seismic_data_, multiple_interval_grid_->GetBackgroundParameters(), err_text);
  //    }
  //    else{
  //      err_text += "Could not set up prior correlations in estimation mode, since this requires a correct setup of the grid and the wells.\n";
  //    }
  //  }
  //  else{
  //    setup_prior_correlation_ = SetupPriorCorrelation(model_settings, input_files, wells_, mapped_blocked_logs_for_correlation_, multiple_interval_grid_->GetIntervalSimboxes(),
  //          prior_facies_, facies_names_, multiple_interval_grid_->GetTrendCubes(), seismic_data_, multiple_interval_grid_->GetBackgroundParameters(), err_text);
  //  }
  //}
  //else{
  //  err_text += "Could not set up prior correlations since this requires seismic data.\n";
  //}


  // 14. Set up TimeLine class
  setup_timeline_ = SetupTimeLine(model_settings, input_files, err_text);

  // 15. Data for gravity inversion
  setup_gravity_inversion_ = SetupGravityInversion(model_settings, input_files, err_text);

  // 16. Data for Travel time Inversion
  //setup_traveltime_inversion_ = SetupTravelTimeInversion(model_settings, input_files, err_text);

  // 17. Depth Conversion
  if (multiple_interval_grid_->GetNIntervals() == 1 && model_settings->getDoDepthConversion())
    setup_depth_conversion_ = SetupDepthConversion(model_settings, input_files, err_text);

  //Punkt o: diverse:
  ReadAngularCorrelations(model_settings, err_text);

  //TODO: Handle if err_text != "".
  if(err_text != "") {
    LogKit::LogFormatted(LogKit::Error,"\n Error when loading and processing data: \n");
    LogKit::LogFormatted(LogKit::Error, err_text);
    exit(1);
  }

}

CommonData::~CommonData() {

  //if (multiple_interval_grid_ != NULL) //H Error in test_case 1
  //  delete multiple_interval_grid_;


  //delete estimation_simbox_;
  //delete full_inversion_volume_;
}

bool CommonData::CreateOuterTemporarySimbox(ModelSettings   * model_settings,
                                            InputFiles      * input_files,
                                            Simbox          & estimation_simbox,
                                            NRLib::Volume   & full_inversion_volume,
                                            std::string     & err_text_common) {

  // parameters
  std::string err_text = "";
  std::string grid_file("");
  int  area_specification  = model_settings->getAreaSpecification();

  LogKit::WriteHeader("Setting up outer modelling grid");

  // FIND SEGY GEOMETRY FOR INVERSION AREA -------------------------------------------------------------------------

  bool estimation_mode_need_ILXL = model_settings->getEstimationMode() &&
                                (area_specification == ModelSettings::AREA_FROM_GRID_DATA ||
                                 area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM ||
                                 area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE ||
                                (model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0 ||
                                (model_settings->getOutputGridFormat() & IO::SEGY) > 0);

  if (forward_modeling_)
    grid_file = input_files->getBackFile(0);    // Get geometry from earth model (Vp)
  else {
    if (model_settings->getEstimationMode() == false || estimation_mode_need_ILXL)
      grid_file = input_files->getSeismicFile(0,0); // Get area from first seismic data volume
  }
  SegyGeometry * ILXL_geometry = NULL; //Geometry with correct XL and IL settings.

  std::string area_type;

  // tilfelle i uten snap to seismic
  if (area_specification == ModelSettings::AREA_FROM_UTM)
  {
    // The geometry is already present in model_settings (geometry_ or geometry_full_ ? )
    area_type = "Model file";
    LogKit::LogFormatted(LogKit::High,"\nArea information has been taken from model file\n");
  }
  // tilfelle ii uten snap to seismic
  else if (area_specification == ModelSettings::AREA_FROM_SURFACE){
    LogKit::LogFormatted(LogKit::High,"\nFinding area information from surface \'"+input_files->getAreaSurfaceFile()+"\'\n");
    area_type = "Surface";
    RotatedSurface surf(input_files->getAreaSurfaceFile());
    SegyGeometry geometry(surf);
    model_settings->setAreaParameters(&geometry);
  }

  else if (area_specification == ModelSettings::AREA_FROM_GRID_DATA || // tilfelle iii
           area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM || // tilfelle i med snap to seismic
           area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE){ // tilfelle ii med snap to seismic
    LogKit::LogFormatted(LogKit::High,"\nFinding inversion area from grid data in file \'"+grid_file+"\'\n");
    area_type = "Grid data";
    std::string tmp_err_text;
    SegyGeometry * geometry;
    GetGeometryFromGridOnFile(grid_file,
                              model_settings->getTraceHeaderFormat(0,0), //Trace header format is the same for all time lapses
                              geometry,
                              tmp_err_text);
    if (geometry!=NULL){
      //tilfelle iv og v
      if (model_settings->getAreaILXL().size() > 0 || model_settings->getSnapGridToSeismicData()) {
        SegyGeometry * full_geometry = geometry;

        std::vector<int> area_ILXL;
        bool got_area = true;

        //
        // Geometry is given as XY, but we snap it to IL and XL.
        //
        if (model_settings->getSnapGridToSeismicData()) {
          SegyGeometry * template_geometry = NULL;
          if (area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM) {
            template_geometry = model_settings->getAreaParameters();
          }
          else if (area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE) {
            RotatedSurface surf(input_files->getAreaSurfaceFile());
            template_geometry = new SegyGeometry(surf);
          }
          else {
            err_text += "CRAVA has been asked to identify a proper ILXL inversion area based\n";
            err_text += "on XY input information, but no UTM coordinates or surface have\n";
            err_text += "been specified in model file.\n";
            got_area = false;
          }
          if (got_area) {
            area_ILXL = full_geometry->findAreaILXL(template_geometry);
          }
        }
        else {
          area_ILXL = model_settings->getAreaILXL();
        }

        if (got_area) {
          try {
            bool interpolated, aligned;
            geometry = full_geometry->GetILXLSubGeometry(area_ILXL, interpolated, aligned);

            std::string text;
            if (interpolated == true) {
              if (aligned == true) {
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
            else if (aligned == true) {
              text  = "Check IL/XL specification: Either start or end of IL and/or XL interval does not\n";
              text += "   align with IL/XL in seismic data, or end IL and/or XL is not an integer multiple\n";
              text += "   of steps away from the start.\n";
              TaskList::addTask(text);
            }
          }
          catch (NRLib::Exception & e) {
            err_text += "Error: "+std::string(e.what());
            geometry->WriteILXL(true);
            geometry = NULL;
          }
        }
        delete full_geometry;
      }
      else {
        geometry->WriteILXL();
      }
      if (err_text == "") {
        model_settings->setAreaParameters(geometry);
        ILXL_geometry = geometry;
      }
    }
    else {
      err_text += tmp_err_text;
    }
  }

  // SET THE AREA FOR THE ESTIMATION SIMBOX AND CHECK THAT IT IS OK -------------------------------------------------------------

  if (err_text == "") {
    bool failed = false;
    const SegyGeometry * area_params = model_settings->getAreaParameters();
    failed = estimation_simbox.setArea(area_params, err_text);
    if (failed)
    {
      WriteAreas(area_params,&estimation_simbox_,area_type);
      err_text += "The specified AREA extends outside the surface(s).\n";
    }
    else{
      LogKit::LogFormatted(LogKit::Low,"\nResolution                x0           y0            lx         ly     azimuth         dx      dy\n");
      LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
      double azimuth = (-1)*estimation_simbox.getAngle()*(180.0/M_PI);
      if (azimuth < 0)
        azimuth += 360.0;
      LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                           area_type.c_str(),
                           estimation_simbox.getx0(), estimation_simbox.gety0(),
                           estimation_simbox.getlx(), estimation_simbox.getly(), azimuth,
                           estimation_simbox.getdx(), estimation_simbox.getdy());
    }

    float min_hor_res = model_settings->getMinHorizontalRes();
    if (estimation_simbox.getdx() < min_hor_res || estimation_simbox.getdy() < min_hor_res){
      failed = true;
      err_text += "The horizontal resolution in dx and dy should normally be above "+NRLib::ToString(min_hor_res)
        +" m. If you need a denser\n sampling, please specify a new <advanced-settings><minimum-horizontal-resolution>\n";
    }

    // SET IL/XL GEOMETRY IN ESTIMATION SIMBOX -------------------------------------------------------------------------------------
    if (!failed) {

      //
      // Set IL/XL information in geometry
      // ---------------------------------
      //
      // Skip for estimation mode if possible:
      //   a) For speed
      //   b) Grid data may not be available.
      if (model_settings->getEstimationMode() == false || estimation_mode_need_ILXL == true) {
        if (ILXL_geometry == NULL) {
          int gridType = IO::findGridType(grid_file);
          bool ilxl_info_available = ((gridType == IO::SEGY) || (gridType == IO::CRAVA));
          if (ilxl_info_available) {
            LogKit::LogFormatted(LogKit::High,"\nFinding IL/XL information from grid data file \'"+grid_file+"\'\n");
            std::string tmp_err_text;
            GetGeometryFromGridOnFile(grid_file,
                                      model_settings->getTraceHeaderFormat(0,0), //Trace header format is the same for all time lapses
                                      ILXL_geometry,
                                      tmp_err_text);
            if (ILXL_geometry == NULL) {
              err_text += tmp_err_text;
            }
          }
          else {
            LogKit::LogFormatted(LogKit::High,"\nCannot extract IL/XL information from non-SEGY grid data file \'"+grid_file+"\'\n");
          }
        }
        if (ILXL_geometry != NULL) {
          if (estimation_simbox.isAligned(ILXL_geometry))
            estimation_simbox.setILXL(ILXL_geometry);
          delete ILXL_geometry;
        }
      }

      // Rotate variograms relative to simbox
      model_settings->rotateVariograms(static_cast<float> (estimation_simbox.getAngle()));

      // SET TOP AND BASE SURFACES FOR THE ESTIMATION SIMBOX -----------------------------------------------

      // Set number of layers in estimation_simbox equal to input file
      int nz = GetNzFromGridOnFile(model_settings, grid_file, err_text);

      // if multiple intervals
      if (model_settings->getIntervalNames().size() > 0) {
        SetSurfacesMultipleIntervals(model_settings, full_inversion_volume, estimation_simbox, nz, input_files, err_text);
      }
      // single interval described by either one or two surfaces
      else{
        SetSurfacesSingleInterval(model_settings, full_inversion_volume, estimation_simbox, nz, input_files->getTimeSurfFiles(), err_text);
      }
    }
  }

  // CALCULATE XY PADDING -------------------------------------------------------------------------
  if (err_text == "") {
    EstimateXYPaddingSizes(&estimation_simbox, model_settings);
  }

  else {
    err_text_common += err_text;
    return false;
  }

  return true;
}


bool CommonData::ReadSeismicData(ModelSettings  * model_settings,
                                 InputFiles     * input_files,
                                 std::string    & err_text_common) {

  std::string err_text = "";

  //Skip if there is no AVO-seismic.
  int timelapse_seismic_files = 0;

  for (size_t i = 0; i < input_files->getTimeLapseSeismicFiles().size(); i++) {
    for (size_t j = 0; j < input_files->getTimeLapseSeismicFiles()[i].size(); j++) {
      if(input_files->getTimeLapseSeismicFiles()[i][j] != "")
        timelapse_seismic_files++;
    }
  }

  if (timelapse_seismic_files == 0)
    return(true);

  //Skip if forward-modeling
  if (forward_modeling_ == true)
    return(true);

  //Skip if estimation mode and wavelet/noise is not set for estimation.
  if (model_settings->getEstimationMode() == true && model_settings->getEstimateWaveletNoise() == false)
    return(true);

  const std::vector<std::vector<std::string> > seismic_timelapse_files = input_files->getTimeLapseSeismicFiles();
  int n_timelapses = model_settings->getNumberOfTimeLapses();

  LogKit::WriteHeader("Reading seismic data");

  for (int this_timelapse = 0; this_timelapse < n_timelapses; this_timelapse++) {

    if (input_files->getNumberOfSeismicFiles(this_timelapse) > 0) {

      std::vector<float> angles = model_settings->getAngle(this_timelapse);
      std::vector<float> offset = model_settings->getLocalSegyOffset(this_timelapse);
      int n_angles              = model_settings->getNumberOfAngles(this_timelapse);

      std::vector<SeismicStorage> seismic_data_angle;

      for (int i = 0; i < n_angles; i++) {

        if (offset[i] < 0)
          offset[i] = model_settings->getSegyOffset(this_timelapse);

        std::string file_name = input_files->getSeismicFile(this_timelapse, i);
        int file_type = IO::findGridType(file_name);

        if (file_type == IO::SEGY) { //Copied out reading part from ModelGeneral::readSegyFile, no resampling

          SegY              * segy   = NULL;
          TraceHeaderFormat * format = model_settings->getTraceHeaderFormat(this_timelapse, i);

          if (format == NULL) { //Unknown format
            std::vector<TraceHeaderFormat*> traceHeaderFormats(0);

            if (model_settings->getTraceHeaderFormat() != NULL)
              traceHeaderFormats.push_back(model_settings->getTraceHeaderFormat());

            segy = new SegY(file_name,
                            offset[i],
                            traceHeaderFormats,
                            true); // Add standard formats to format search
          }
          else { //Known format, read directly.
            segy = new SegY(file_name, offset[i], *format);
          }

          float guard_zone = model_settings->getGuardZone();

          bool cover_ok = false;
          //H-TODO Move check that seismic data cover grid until interval_simboxes are made and before inversion.
          cover_ok = CheckThatDataCoverGrid(segy,
                                            offset[i],
                                            &estimation_simbox_,
                                            guard_zone,
                                            err_text);
          if (cover_ok) {
            float padding = 2*guard_zone;
            bool relative_padding = false;
            bool only_volume = true;

            segy->ReadAllTraces(&full_inversion_volume_,
                                padding,
                                only_volume,
                                relative_padding);
            segy->CreateRegularGrid(); //sets geometry

            SeismicStorage seismicdata(file_name, SeismicStorage::SEGY, angles[i], segy);
            seismic_data_angle.push_back(seismicdata);
          }
          else {
            err_text += "Data from segy-file " + file_name + " is not read.\n";
          }

        } //SEGY
        else if (file_type == IO::STORM || file_type == IO::SGRI) { //From ModelGeneral::readStormFile
          StormContGrid * stormgrid = NULL;

          try {
            stormgrid = new StormContGrid(0,0,0);
            stormgrid->ReadFromFile(file_name);
          }
          catch (NRLib::Exception & e) {
            err_text += "Error when reading storm-file " + file_name +": " + NRLib::ToString(e.what()) + "\n";
            break;
          }

          bool cover_ok = false;
          bool scale    = false;
          if (file_type == IO::SGRI)
            scale = true;
          std::string err_text_tmp = "";

          //H-TODO Move this check until after interval_simboxes are made and before inversion.
          //This cover_check failes in test_cast 15 (against estimation_simbox)
          cover_ok = true;
          //cover_ok = CheckThatDataCoverGrid(stormgrid, //Checks that stormgrid.z > simbox.z
          //                                  &estimation_simbox_,
          //                                  model_settings->getGuardZone(),
          //                                  err_text_tmp,
          //                                  scale);

          if (cover_ok == false) {
            err_text += "Data from storm file " + file_name + " is not read.\n";
            err_text += err_text_tmp;
          }

          if (cover_ok) {
            SeismicStorage seismicdata_tmp;

            if (file_type == IO::STORM)
              seismicdata_tmp = SeismicStorage(file_name, SeismicStorage::STORM, angles[i], stormgrid);
            else
              seismicdata_tmp = SeismicStorage(file_name, SeismicStorage::SGRI, angles[i], stormgrid);

            seismic_data_angle.push_back(seismicdata_tmp);
          }
        } //STORM / SGRI
        else if (file_type == IO::CRAVA) {

          int nx_pad = estimation_simbox_.GetNXpad();
          int ny_pad = estimation_simbox_.GetNYpad();
          int nz_pad = 0;  //Not set before ReadSeismicData. Get it from file.

          GetZPaddingFromCravaFile(file_name, err_text, nz_pad);

          FFTGrid *  grid = CreateFFTGrid(estimation_simbox_.getnx(),
                                          estimation_simbox_.getny(),
                                          estimation_simbox_.getnz(),
                                          nx_pad,
                                          ny_pad,
                                          nz_pad,
                                          model_settings->getFileGrid());

          std::string angle    = NRLib::ToString(angles[i]*(180/M_PI), 1);
          std::string par_name = "Seismic data angle stack "+angle;
          LogKit::LogFormatted(LogKit::Low,"\nReading grid \'"+par_name+"\' from file "+file_name);

          grid->setAccessMode(FFTGrid::RANDOMACCESS);
          grid->readCravaFile(file_name, err_text);

          SeismicStorage seismicdata(file_name, SeismicStorage::FFTGRID, angles[i], grid);
          seismic_data_angle.push_back(seismicdata);

        }
        else {
          err_text += "Error when reading file " + file_name +". File type not recognized.\n";
        }
      } //nAngles

      seismic_data_[this_timelapse] = seismic_data_angle;

    }//ifSeismicFiles
  } //n_timeLapses

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

bool
CommonData::CheckThatDataCoverGrid(const SegY   * segy,
                                   float          offset,
                                   const Simbox & time_simbox,
                                   float          guard_zone,
                                   std::string &  err_text) const {
  // Seismic data coverage (translate to CRAVA grid by adding half a grid cell)
  float dz = segy->GetDz();
  float z0 = offset + 0.5f*dz;
  float zn = z0 + (segy->GetNz() - 1)*dz;

  // Top and base of interval of interest
  float top_grid = static_cast<float>(time_simbox.getTopZMin());
  float bot_grid = static_cast<float>(time_simbox.getBotZMax());

  // Find guard zone
  float top_guard = top_grid - guard_zone;
  float bot_guard = bot_grid + guard_zone;

  if (top_guard < z0) {
    float z0_new = z0 - ceil((z0 - top_guard)/dz)*dz;
    err_text += "\nThere is not enough seismic data above the interval of interest. The seismic data\n";
    err_text += "must start at "+NRLib::ToString(z0_new)+"ms (in CRAVA grid) to allow for a ";
    err_text += NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n";
    err_text += "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n";
    err_text += "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n";
    err_text += "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n";
    err_text += "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n";
    err_text += "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n";
    err_text += "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n";

    return false;
  }
  if (bot_guard > zn) {
    float zn_new = zn + ceil((bot_guard - zn)/dz)*dz;
    err_text += "\nThere is not enough seismic data below the interval of interest. The seismic data\n";
    err_text += "must end at "+NRLib::ToString(zn_new)+"ms (in CRAVA grid) to allow for a ";
    err_text += NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n";
    err_text += "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n";
    err_text += "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n";
    err_text += "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n";
    err_text += "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n";
    err_text += "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n";
    err_text += "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n";

    return false;
  }

  return true;
}

bool
CommonData::CheckThatDataCoverGrid(StormContGrid * stormgrid,
                                   const Simbox  & simbox,
                                   float           guard_zone,
                                   std::string   & err_text,
                                   bool            scale) const {

  double storm_top = stormgrid->GetTopZMin(stormgrid->GetNI(), stormgrid->GetNJ());
  double storm_bot = stormgrid->GetBotZMax(stormgrid->GetNI(), stormgrid->GetNJ());

  if (scale) {
    storm_top /= 0.001;
    storm_bot /= 0.001;
  }

  // Top and base of interval of interest
  double top_grid = simbox.getTopZMin();
  double bot_grid = simbox.getBotZMax();

  // Find guard zone
  double top_guard = top_grid - guard_zone;
  double bot_guard = bot_grid + guard_zone;

  if (top_guard < storm_top) {
    err_text += "\nThere is not enough seismic data above the interval of interest. The seismic data\n";
    err_text += "must start at "+NRLib::ToString(top_guard)+"ms (in CRAVA grid) to allow for a ";
    err_text += NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n";
    err_text += "  Seismic data start (CRAVA grid) : "+NRLib::ToString(storm_top,1)+"\n";
    err_text += "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n";
    err_text += "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n";
    err_text += "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n";
    err_text += "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n";
    err_text += "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(storm_bot,1)+"\n";

    return false;
  }

  if (bot_guard > storm_bot) {
    err_text += "\nThere is not enough seismic data below the interval of interest. The seismic data\n";
    err_text += "must end at "+NRLib::ToString(bot_guard)+"ms (in CRAVA grid) to allow for a ";
    err_text += NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n";
    err_text += "  Seismic data start (CRAVA grid) : "+NRLib::ToString(storm_top,1)+"\n";
    err_text += "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n";
    err_text += "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n";
    err_text += "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n";
    err_text += "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n";
    err_text += "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(storm_bot,1)+"\n";

    return false;
  }

  return true;
}


bool CommonData::ReadWellData(ModelSettings                  * model_settings,
                              Simbox                         * estimation_simbox,
                              InputFiles                     * input_files,
                              std::vector<std::string>       & log_names,
                              const std::vector<std::string> & log_names_from_user,
                              const std::vector<bool>        & inverse_velocity,
                              bool                             facies_log_given,
                              std::string                    & err_text_common) {

  std::string err_text = "";

  // Get all log names given by user
  for (size_t i = 0; i < log_names_from_user.size(); i++) {
    if (log_names_from_user[i] != "")
      log_names.push_back(log_names_from_user[i]);
  }

  int n_wells  = model_settings->getNumberOfWells();

  if (n_wells == 0)
    return true;

  int n_facies = facies_names_.size();
  try {
    if (n_wells > 0)
      LogKit::WriteHeader("Reading wells");

    int no_hit            = 0;
    int empty             = 0;
    int facies_log_not_ok = 0;
    int upwards           = 0;
    std::vector<bool>  valid_index(n_wells);
    std::vector<int>   n_merges(n_wells);
    std::vector<int>   n_invalid_vp(n_wells);
    std::vector<int>   n_invalid_vs(n_wells);
    std::vector<int>   n_invalid_rho(n_wells);
    std::vector<float> rank_corr(n_wells);
    std::vector<float> dev_angle(n_wells);

    std::vector<std::string> well_names(n_wells);
    std::vector<bool> well_synthetic_vs_log(n_wells);

    std::vector<std::vector<int> > facies_count(n_wells);
    for (int i = 0; i < n_wells; i++) {
      facies_count[i].resize(n_facies);
    }

    for (int i = 0; i < n_wells; i++) {

      std::string well_file_name = input_files->getWellFile(i);
      bool read_ok = false;

      // If the facies log is given, it is always entry 4 in the log name list
      std::string facies_log = "FACIES";
      if (log_names_from_user.size() > 4 && log_names_from_user[4] != "")
        facies_log = log_names_from_user[4];

      NRLib::Well new_well(well_file_name, read_ok, facies_log);
      LogKit::LogFormatted(LogKit::Low, new_well.GetWellName()+" : \n");

      std::vector<int> facies_nr;
      std::vector<std::string> facies_names;

      //Process logs
      if (well_file_name.find(".nwh",0) != std::string::npos)
        ProcessLogsNorsarWell(new_well, log_names, inverse_velocity, facies_log_given, err_text);
      else if (well_file_name.find(".rms",0) != std::string::npos || well_file_name.find(".txt",0) != std::string::npos)
        ProcessLogsRMSWell(new_well, log_names, inverse_velocity, facies_log_given, err_text);

      //Cut wells against full_inversion_volume
      if (err_text == "")
        CutWell(well_file_name, new_well, full_inversion_volume_);

      //Store facies names.
      if (model_settings->getFaciesLogGiven()) {
        ReadFaciesNamesFromWellLogs(new_well, facies_nr, facies_names);
      }

      facies_nr_wells_.push_back(facies_nr);
      facies_names_wells_.push_back(facies_names);

      new_well.SetUseForFaciesProbabilities(model_settings->getIndicatorFacies(i));
      new_well.SetUseForFiltering(model_settings->getIndicatorFilter(i));
      new_well.SetRealVsLog(model_settings->getIndicatorRealVs(i));
      new_well.SetUseForBackgroundTrend(model_settings->getIndicatorBGTrend(i));
      new_well.SetUseForWaveletEstimation(model_settings->getIndicatorWavelet(i));

      //Check if well is valid
      bool well_valid = true;
      bool facies_ok  = true;

      if (new_well.CheckSimbox(estimation_simbox) == 1) {
        well_valid = false;
        no_hit++;
        TaskList::addTask("Consider increasing the inversion volume such that well "+new_well.GetWellName()+ " can be included");
      }
      if (new_well.GetNData() == 0) {
        LogKit::LogFormatted(LogKit::Low,"  IGNORED (no log entries found)\n");
        well_valid = false;
        empty++;
        TaskList::addTask("Check the log entries in well "+new_well.GetWellName()+".");
      }
      //Check well for valid facies
      if (new_well.HasFaciesLog() == true) {
        const std::vector<int> & facies_log = new_well.GetDiscLog("Facies");
        int n_facies                        = new_well.GetNFacies();
        bool facies_ok_tmp                  = false;
        for (size_t j = 0; j < facies_log.size(); j++) {
          if (facies_log[j] != IMISSING) {
            for (int k = 0; k < n_facies; k++) {
              if (facies_nr[k] == facies_log[j]) {
                facies_ok_tmp = true;
                break;
              }
            }
            if (facies_ok_tmp == false) { //Found a facies in facies_log that is not defined in facies_nr
              facies_ok = false;
              break;
            }
          }
        }
      }
      if (facies_ok == false) {
        LogKit::LogFormatted(LogKit::Low,"   IGNORED (facies log has wrong entries)\n");
        well_valid = false;
        facies_log_not_ok++;
        TaskList::addTask("Check the facies logs in well "+new_well.GetWellName()+".\n       The facies logs in this well are wrong and the well is ignored");
      }
      bool monotonous = RemoveDuplicateLogEntriesFromWell(new_well, model_settings, estimation_simbox, n_merges[i]);
      if (monotonous == false) {
        LogKit::LogFormatted(LogKit::Low,"   IGNORED (well is too far from monotonous in time)\n");
        well_valid = false;
        upwards++;
        TaskList::addTask("Check the TWT log in well "+new_well.GetWellName()+".\n       The well is moving too much upwards, and the well is ignored");
      }

      well_names[i]            = new_well.GetWellName();
      well_synthetic_vs_log[i] = new_well.HasSyntheticVsLog();

      if (read_ok == false) {
        err_text += "Well format of file " + well_file_name + " not recognized.\n";
      }
      else if (well_valid == true) {

        valid_index[i] = true;
        SetWrongLogEntriesInWellUndefined(new_well, model_settings, n_invalid_vp[i], n_invalid_vs[i], n_invalid_rho[i]);
        FilterLogs(new_well, model_settings);
        LookForSyntheticVsLog(new_well, model_settings, rank_corr[i]);
        CalculateDeviation(new_well, model_settings, dev_angle[i], estimation_simbox);

        if (n_facies > 0)
          CountFaciesInWell(new_well, estimation_simbox, n_facies, facies_nr, facies_count[i]);

        wells_.push_back(new_well);

        if (model_settings->getFaciesLogGiven())
          facies_log_wells_.push_back(true);
        else
          facies_log_wells_.push_back(false);
      }
      else
        valid_index[i] = false;

    } //n_wells

    //Combines facies information from wells
    if (model_settings->getFaciesLogGiven() && err_text == "")
      SetFaciesNamesFromWells(model_settings, err_text);

    //
    // Write summary.
    //
    LogKit::LogFormatted(LogKit::Low,"\n");
    LogKit::LogFormatted(LogKit::Low,"                                      Invalid                                    \n");
    LogKit::LogFormatted(LogKit::Low,"Well                    Merges      Vp   Vs  Rho  synthVs/Corr    Deviated/Angle \n");
    LogKit::LogFormatted(LogKit::Low,"---------------------------------------------------------------------------------\n");
    for (int i = 0; i < n_wells; i++) {
      if (valid_index[i])
        LogKit::LogFormatted(LogKit::Low,"%-23s %6d    %4d %4d %4d     %3s / %5.3f      %3s / %4.1f\n",
        well_names[i].c_str(),
        n_merges[i],
        n_invalid_vp[i],
        n_invalid_vs[i],
        n_invalid_rho[i],
        (well_synthetic_vs_log[i] ? "yes" : " no"),
        rank_corr[i],
        (dev_angle[i] > model_settings->getMaxDevAngle() ? "yes" : " no"),
        dev_angle[i]);
      else
        LogKit::LogFormatted(LogKit::Low,"%-23s      -       -    -    -       - /     -       -  /    -\n",
        well_names[i].c_str());
    }

    //
    // Print facies count for each well
    //
    if (n_facies > 0) {
      //
      // Probabilities
      //
      LogKit::LogFormatted(LogKit::Low,"\nFacies distributions for each well: \n");
      LogKit::LogFormatted(LogKit::Low,"\nWell                    ");
      for (int i = 0; i < n_facies; i++)
        LogKit::LogFormatted(LogKit::Low,"%12s ",facies_names_[i].c_str());
      LogKit::LogFormatted(LogKit::Low,"\n");
      for (int i = 0 ; i < 24+13*n_facies ; i++)
        LogKit::LogFormatted(LogKit::Low,"-");
      LogKit::LogFormatted(LogKit::Low,"\n");
      for (int i = 0; i < n_wells; i++) {
        if (valid_index[i]) {
          float tot = 0.0;
          for (int f = 0; f < n_facies; f++)
            tot += static_cast<float>(facies_count[i][f]);
          LogKit::LogFormatted(LogKit::Low,"%-23s ", well_names[i].c_str());
          for (int f = 0; f < n_facies; f++) {
            if (tot > 0) {
              float facies_prob = static_cast<float>(facies_count[i][f])/tot;
              LogKit::LogFormatted(LogKit::Low,"%12.4f ",facies_prob);
            }
            else
              LogKit::LogFormatted(LogKit::Low,"         -   ");
          }
          LogKit::LogFormatted(LogKit::Low,"\n");
        }
        else {
          LogKit::LogFormatted(LogKit::Low,"%-23s ", well_names[i].c_str());
          for (int f = 0; f < n_facies; f++)
            LogKit::LogFormatted(LogKit::Low,"         -   ");
          LogKit::LogFormatted(LogKit::Low,"\n");

        }
      }
      LogKit::LogFormatted(LogKit::Low,"\n");
      //
      // Counts
      //
      LogKit::LogFormatted(LogKit::Medium,"\nFacies counts for each well: \n");
      LogKit::LogFormatted(LogKit::Medium,"\nWell                    ");
      for (int i = 0; i < n_facies; i++)
        LogKit::LogFormatted(LogKit::Medium,"%12s ",facies_names_[i].c_str());
      LogKit::LogFormatted(LogKit::Medium,"\n");
      for (int i = 0; i < 24+13*n_facies; i++)
        LogKit::LogFormatted(LogKit::Medium,"-");
      LogKit::LogFormatted(LogKit::Medium,"\n");
      for (int i = 0; i < n_wells; i++) {
        if (valid_index[i]) {
          float tot = 0.0;
          for (int f = 0; f < n_facies; f++)
            tot += static_cast<float>(facies_count[i][f]);
          LogKit::LogFormatted(LogKit::Medium,"%-23s ", well_names[i].c_str());
          for (int f = 0; f < n_facies; f++) {
            LogKit::LogFormatted(LogKit::Medium,"%12d ",facies_count[i][f]);
          }
          LogKit::LogFormatted(LogKit::Medium,"\n");
        }
        else {
          LogKit::LogFormatted(LogKit::Medium,"%-23s ", well_names[i].c_str());
          for (int f = 0; f < n_facies; f++)
            LogKit::LogFormatted(LogKit::Medium,"         -   ");
          LogKit::LogFormatted(LogKit::Medium,"\n");
        }
      }
      LogKit::LogFormatted(LogKit::Medium,"\n");
    }

    //Update n_wells
    int n_wells_all = n_wells;
    n_wells = wells_.size();
    model_settings->setNumberOfWells(n_wells);

    if (no_hit > 0)
      LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) do not hit the inversion volume and will be ignored.\n", no_hit);
    if (empty > 0)
      LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) contain no log entries and will be ignored.\n", empty);
    if (facies_log_not_ok > 0)
      LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) have wrong facies logs and will be ignored.\n", facies_log_not_ok);
    if (upwards > 0)
      LogKit::LogFormatted(LogKit::Low,"\nWARNING: %d well(s) are moving upwards in TWT and will be ignored.\n", upwards);
    if (n_wells == 0 && model_settings->getNoWellNedded()==false) {
      LogKit::LogFormatted(LogKit::Low,"\nERROR: There are no wells left for data analysis. Please check that the inversion area given");
      LogKit::LogFormatted(LogKit::Low,"\n       below is correct. If it is not, you probably have problems with coordinate scaling.");
      LogKit::LogFormatted(LogKit::Low,"\n                                   X0          Y0        DeltaX      DeltaY      Angle");
      LogKit::LogFormatted(LogKit::Low,"\n       -------------------------------------------------------------------------------");
      LogKit::LogFormatted(LogKit::Low,"\n       Inversion area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n",
        estimation_simbox->getx0(), estimation_simbox->gety0(),
        estimation_simbox->getlx(), estimation_simbox->getly(),
        (estimation_simbox->getAngle()*180)/M_PI);
      err_text += "No wells available for estimation.";
    }

    if (n_facies > 0) {
      int fc;
      for (int i = 0; i < n_facies; i++) {
        fc = 0;
        for (int j = 0; j < n_wells_all; j++) {
          fc+=facies_count[j][i];
        }
        if (fc == 0) {
          LogKit::LogFormatted(LogKit::Low,"\nWARNING: Facies %s is not observed in any of the wells, and posterior facies probability can not be estimated for this facies.\n",facies_names_[i].c_str() );
          TaskList::addTask("In order to estimate prior facies probability for facies "+ facies_names_[i] + " add wells which contain observations of this facies.\n");
        }
      }
    }
  }
  catch (NRLib::Exception & e) {
    err_text += "Error: " + NRLib::ToString(e.what());
  }

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

bool CommonData::RemoveDuplicateLogEntriesFromWell(NRLib::Well   & well,
                                                   ModelSettings * model_settings,
                                                   const Simbox  * simbox,
                                                   int           & n_merges)
{
  bool debug      = false;
  bool monotonous = true; //Check that well does not move unreasonably much upwards.
  float max_upward_dist = -3.0f*static_cast<float>(simbox->getdz()); //Do not allow upwards movement of more than three cells. Calibrate.

  float min_merge_dist = model_settings->getMaxMergeDist();
  int   n_data         = well.GetNData();

  std::vector<double> x_pos_resampled;
  std::vector<double> y_pos_resampled;
  std::vector<double> z_pos_resampled;  // time step
  std::vector<double> vp_resampled;
  std::vector<double> vs_resampled;
  std::vector<double> rho_resampled;
  std::vector<int>    facies_resampled; // Always included (for convenience)

  const std::vector<double> & x_pos = well.GetContLog("X_pos");
  const std::vector<double> & y_pos = well.GetContLog("Y_pos");
  const std::vector<double> & z_pos = well.GetContLog("Z_pos");
  const std::vector<double> & vp    = well.GetContLog("Vp");
  const std::vector<double> & vs    = well.GetContLog("Vs");
  const std::vector<double> & rho   = well.GetContLog("Rho");

  int ii = 0;
  int istart = 0;                                // First element in merge
  int iend;                                      // Last element in merge
  while (istart < n_data) {                      // Loop over elements in original log
    if (istart == n_data - 1)                    // Are we at last element?
      iend = istart;
    else {
      iend = istart + 1;                         // Start looking one element ahead
      while (z_pos[iend] - z_pos[istart] < min_merge_dist && iend < n_data - 1) {
        iend++;
        if(monotonous == true && z_pos[iend] - z_pos[istart] < max_upward_dist)
          monotonous = false;
      }
      iend--;
    }

    bool print_to_screen = debug && iend != istart;
    if (print_to_screen)
      LogKit::LogFormatted(LogKit::Low,"Merge log entries %d -> %d into new entry %d\n",istart,iend,ii);

    MergeCells("time", z_pos_resampled, z_pos, ii, istart, iend, print_to_screen);
    MergeCells("x   ", x_pos_resampled, x_pos, ii, istart, iend, print_to_screen);
    MergeCells("y   ", y_pos_resampled, y_pos, ii, istart, iend, print_to_screen);
    MergeCells("Vp  ", vp_resampled,    vp,    ii, istart, iend, print_to_screen);
    MergeCells("Vs  ", vs_resampled,    vs,    ii, istart, iend, print_to_screen);
    MergeCells("Rho ", rho_resampled,   rho,   ii, istart, iend, print_to_screen);
    if(well.HasFaciesLog()) {
      const std::vector<int>  & facies = well.GetDiscLog("Facies");
      MergeCellsDiscrete("Facies ", facies_resampled, facies, ii, istart, iend, print_to_screen);
    }

    if (print_to_screen)
      LogKit::LogFormatted(LogKit::Low,"\n");

    istart = iend + 1;
    ii++;
  }

  n_merges = 0;
  if (ii != n_data) {
    n_merges = n_data - ii;
    LogKit::LogFormatted(LogKit::Low,"   Duplicate log entries merged with neighbour : %d\n",n_merges);
    n_data = ii;
  }

  well.SetNumberOfData(n_data);

  //Store new logs
  well.RemoveContLog("X_pos");
  well.AddContLog("X_pos", x_pos_resampled);

  well.RemoveContLog("Y_pos");
  well.AddContLog("Y_pos", y_pos_resampled);

  well.RemoveContLog("Z_pos");
  well.AddContLog("Z_pos", z_pos_resampled);

  well.RemoveContLog("Vp");
  well.AddContLog("Vp", vp_resampled);

  well.RemoveContLog("Vs");
  well.AddContLog("Vs", vs_resampled);

  well.RemoveContLog("Rho");
  well.AddContLog("Rho", rho_resampled);

  if (well.HasFaciesLog()) {
    well.RemoveDiscLog("Facies");
    well.AddDiscLog("Facies", facies_resampled);
  }

  return(monotonous);
}

//----------------------------------------------------------------------------
void CommonData::MergeCells(const std::string         & name,
                            std::vector<double>       & pos_resampled,
                            const std::vector<double> & pos,
                            int                         ii,
                            int                         istart,
                            int                         iend,
                            bool                        print_to_screen)
{
  int n_sample = 0;
  pos_resampled.push_back(0.0);

  for (int i = istart; i < iend + 1; i++) {
    if (pos[i] != RMISSING) {
      pos_resampled[ii] += pos[i];
      n_sample++;
      if (print_to_screen)
        LogKit::LogFormatted(LogKit::Low,"%s     Old:%d   pos = %.3f\n",name.c_str(),i,pos[i]);
    }
  }
  if (n_sample == 0)
    pos_resampled[ii] = RMISSING;
  else
    pos_resampled[ii] /= n_sample;
  if (print_to_screen)
    LogKit::LogFormatted(LogKit::Low,"%s     New:%d   pos = %.3f\n",name.c_str(),ii,pos_resampled[ii]);
}

void CommonData::MergeCellsDiscrete(const std::string      & name,
                                    std::vector<int>       & log_resampled,
                                    const std::vector<int> & log,
                                    int                      ii,
                                    int                      istart,
                                    int                      iend,
                                    bool                     print_to_screen)
{
  int n_sample = 0;
  for (int i = istart; i < iend + 1; i++) {
    if (log[i] != RMISSING) {
      //log_resampled[ii] += log[i];
      n_sample++;
      if (print_to_screen)
        LogKit::LogFormatted(LogKit::Low,"%s     Old:%d   log = %d\n",name.c_str(),i,log[i]);
    }
  }

  log_resampled.push_back(log[istart]);
  if (print_to_screen)
    LogKit::LogFormatted(LogKit::Low,"%s     New:%d   log = %d\n",name.c_str(),ii,log_resampled[ii]);

}

//----------------------------------------------------------------------------
void CommonData::SetWrongLogEntriesInWellUndefined(NRLib::Well   & well,
                                                   ModelSettings * model_settings,
                                                   int           & count_vp,
                                                   int           & count_vs,
                                                   int           & count_rho)
{
  //
  // Log values outside the minimum and maximum limits given below
  // are set to RMISSING as such values are very unlikely.
  //
  bool debug = true;

  float vp_min  = model_settings->getAlphaMin();
  float vp_max  = model_settings->getAlphaMax();
  float vs_min  = model_settings->getBetaMin();
  float vs_max  = model_settings->getBetaMax();
  float rho_min = model_settings->getRhoMin();
  float rho_max = model_settings->getRhoMax();

  count_vp  = 0;
  count_vs  = 0;
  count_rho = 0;

  std::vector<double> & vp  = well.GetContLog("Vp");
  std::vector<double> & vs  = well.GetContLog("Vs");
  std::vector<double> & rho = well.GetContLog("Rho");
  const std::vector<double> & z_pos = well.GetContLog("Z_pos");

  int n_data = well.GetNData();

  for (int i = 0; i < n_data; i++) {
    if (vp[i] != RMISSING && (vp[i] < vp_min || vp[i] > vp_max)) {
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   Set undefined:   time = %.2f   Vp = %.2f\n",z_pos[i],vp[i]);
      vp[i] = RMISSING;
      count_vp++;
    }
    if (vs[i] != RMISSING && (vs[i] < vs_min || vs[i] > vs_max)) {
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   Set undefined:   time = %.2f   Vs = %.2f\n",z_pos[i],vs[i]);
      vs[i] = RMISSING;
      count_vs++;
    }
    if (rho[i] != RMISSING && (rho[i] < rho_min || rho[i] > rho_max)) {
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   Set undefined:   time = %.2f   Rho = %.2f\n",z_pos[i],rho[i]);
      rho[i] = RMISSING;
      count_rho++;
    }
  }

  if (count_vp > 0 || count_vs || count_rho> 0)
    LogKit::LogFormatted(LogKit::Low,"   Log entries have been set undefined (Vp:%d Vs:%d Rho:%d)\n",
                         count_vp, count_vs, count_rho);
}

//----------------------------------------------
bool compare(const std::pair<int, float>& i1, const std::pair<int, float>& i2)
{
  return (i1.second < i2.second);
}

void CommonData::LookForSyntheticVsLog(NRLib::Well   & well,
                                       ModelSettings * model_settings,
                                       float         & rank_correlation)
{
  float corr_threshold = model_settings->getMaxRankCorr();

  //
  // Estimate the correlation between Vp and Vs logs. To be able to identify
  // nonlinear relationships between the logs we use rank correlation.
  //
  typedef std::pair<int, float> Item;
  std::vector<Item> sorted_vp;
  std::vector<Item> sorted_vs;

  //
  // Store Vp and Vs in sortable structs. Note that we can only use nonmissing values.
  //
  int n_data = well.GetNData();
  const std::vector<double> & vp = well.GetContLog("Vp");
  const std::vector<double> & vs = well.GetContLog("Vs");

  int real_vs_log = well.GetRealVsLog();

  for (int i = 0; i < n_data; i++) {
    if (vp[i] != RMISSING && vs[i] != RMISSING) {
      sorted_vp.push_back(Item(i, vp[i]));
      sorted_vs.push_back(Item(i, vs[i]));
    }
  }

  int n = static_cast<int>(sorted_vp.size());
  if (n > 0) {
    //
    // Sort Vp and Vs logs.
    //
    std::sort(sorted_vp.begin(), sorted_vp.end(), compare);
    std::sort(sorted_vs.begin(), sorted_vs.end(), compare);

    //
    // Estimate correlation between sorted alpha and beta sorted by alpha
    //
    float mean = float(n)/2.0f; // We start the indexing at 0 rather than 1
    float cov_rank = 0.0;       // Covariance between ranks of alpha and beta (sorted on alpha)
    float var_rank = 0.0;       // Variance in ranks of alpha and beta (which are equal)

    for (int i = 0; i < n; i++) {
      var_rank +=(i - mean)*(i - mean);
      for (int j = 0; j < n; j++) {

        if (sorted_vs[j].first == sorted_vp[i].first) {
          cov_rank += (j - mean)*(i - mean);
        }
      }
    }
    rank_correlation = cov_rank/var_rank; // Skip division by n-1 in both nominator and denominator

    if (rank_correlation > corr_threshold) {
      if (real_vs_log == ModelSettings::NOTSET) {
        LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. Treating Vs log as synthetic.\n",rank_correlation);
        real_vs_log = ModelSettings::NO;
      }
      else {
        if (real_vs_log == ModelSettings::YES)
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f, but well log is defined as real.\n",rank_correlation);
        else
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as synthetic.)\n",rank_correlation);
      }
    }
    else {
      switch(real_vs_log) {
        case ModelSettings::YES :
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as real.)\n",rank_correlation);
          break;
        case ModelSettings::NO :
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is defined as synthetic.)\n",rank_correlation);
          break;
        default :
          LogKit::LogFormatted(LogKit::Low,"   Vp-Vs rank correlation is %5.3f. (Well log is treated as real.)\n",rank_correlation);
          break;
      }
    }

  }
  else
  {
    LogKit::LogFormatted(LogKit::Low,"   Cannot calculate Vp-Vs rank correlation. One or both logs are empty.\n");
  }

  well.SetRealVsLog(real_vs_log);

  bool use_filter    = model_settings->getUseFilterForFaciesProb();
  bool use_vp_vs_rho = model_settings->getNoVsFaciesProb() == false;
  int  use_for_facies_probabilities = well.GetUseForFaciesProbabilities();

  if(use_filter && use_vp_vs_rho && real_vs_log == ModelSettings::NO && use_for_facies_probabilities == ModelSettings::NOTSET)
    use_for_facies_probabilities = ModelSettings::NO;

  well.SetUseForFaciesProbabilities(use_for_facies_probabilities);
}

//----------------------------------------------------------------------------
void CommonData::CountFaciesInWell(NRLib::Well            & well,
                                   Simbox                 * simbox,
                                   int                      n_facies,
                                   const std::vector<int> & facies_nr,
                                   std::vector<int>       & facies_count)
{
  for (int i = 0; i < n_facies; i++)
    facies_count[i] = 0;

  int n_data = well.GetNData();
  const std::vector<double> & x_pos  = well.GetContLog("X_pos");
  const std::vector<double> & y_pos  = well.GetContLog("Y_pos");
  const std::vector<double> & z_pos  = well.GetContLog("Z_pos");
  const std::vector<double> & vp     = well.GetContLog("Vp");
  const std::vector<double> & vs     = well.GetContLog("Vs");
  const std::vector<double> & rho    = well.GetContLog("Rho");
  const std::vector<int>    & facies = well.GetDiscLog("Facies");

  if (n_facies > 0) {
    for(int i = 0; i < n_data; i++) {
      if (facies[i] != IMISSING) {
        //
        // Count facies only when logs
        //
        if (vp[i] != RMISSING && vs[i] != RMISSING && rho[i] != RMISSING) {
          //
          // Count facies only for part of well inside simbox
          //
          if(simbox->isInside(x_pos[i], y_pos[i])) {
            if(z_pos[i] > simbox->getTop(x_pos[i], y_pos[i]) && z_pos[i] < simbox->getBot(x_pos[i], y_pos[i]))  {
              //
              // Find facies index...
              //
              for (int j = 0; j < n_facies; j++) {
                if (facies_nr[j] == facies[i]) {
                  facies_count[j]++;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void CommonData::FilterLogs(NRLib::Well   & well,
                            ModelSettings * model_settings)
{
  float max_hz_background = model_settings->getMaxHzBackground();
  float max_hz_seismic    = model_settings->getMaxHzSeismic();

  int n_data = well.GetNData();

  std::vector<double> vp_interpolated(n_data);
  std::vector<double> vs_interpolated(n_data);
  std::vector<double> rho_interpolated(n_data);

  std::vector<double> vp_resampled(n_data);
  std::vector<double> vs_resampled(n_data);
  std::vector<double> rho_resampled(n_data);

  std::vector<double> vp_filtered(n_data);
  std::vector<double> vs_filtered(n_data);
  std::vector<double> rho_filtered(n_data);

  std::vector<double> time_resampled(n_data);
  double dt = 0.0;

  std::vector<double> vp_background_resolution(n_data,  RMISSING);
  std::vector<double> vs_background_resolution(n_data,  RMISSING);
  std::vector<double> rho_background_resolution(n_data, RMISSING);

  std::vector<double> vp_seismic_resolution(n_data,  RMISSING);
  std::vector<double> vs_seismic_resolution(n_data,  RMISSING);
  std::vector<double> rho_seismic_resolution(n_data, RMISSING);

  const std::vector<double> & vp  = well.GetContLog("Vp"); //_raw_logs  = GetVpRawLogs();
  const std::vector<double> & vs  = well.GetContLog("Vs");
  const std::vector<double> & rho = well.GetContLog("Rho");

  const std::vector<double> & z_pos = well.GetContLog("Z_pos");
  //
  // Time
  //
  bool filtered = ResampleTime(time_resampled, z_pos, n_data, dt); //False if well not monotonous in time.

  if (filtered) {
    //
    // Vp
    //
    ResampleLog(vp_resampled, vp, z_pos, time_resampled, n_data, dt);         // May generate missing values
    InterpolateLog(vp_interpolated, vp_resampled, n_data);                    // Interpolate missing values

    ApplyFilter(vp_filtered, vp_interpolated, n_data, dt, max_hz_background);
    ResampleLog(vp_resampled, vp_filtered, time_resampled, z_pos, n_data, dt);
    InterpolateLog(vp_background_resolution, vp_resampled, n_data);

    ApplyFilter(vp_filtered, vp_interpolated, n_data, dt, max_hz_seismic);
    ResampleLog(vp_resampled, vp_filtered, time_resampled, z_pos, n_data, dt);
    InterpolateLog(vp_seismic_resolution, vp_resampled, n_data);

    //
    // Vs
    //
    ResampleLog(vs_resampled, vs, z_pos, time_resampled, n_data, dt);
    InterpolateLog(vs_interpolated, vs_resampled, n_data);

    ApplyFilter(vs_filtered, vs_interpolated, n_data, dt, max_hz_background);
    ResampleLog(vs_resampled, vs_filtered, time_resampled, z_pos, n_data, dt);
    InterpolateLog(vs_background_resolution, vs_resampled, n_data);

    ApplyFilter(vs_filtered, vs_interpolated, n_data, dt, max_hz_seismic);
    ResampleLog(vs_resampled, vs_filtered, time_resampled, z_pos, n_data, dt);
    InterpolateLog(vs_seismic_resolution, vs_resampled, n_data);

    //
    // Rho
    //
    ResampleLog(rho_resampled, rho, z_pos, time_resampled, n_data, dt);
    InterpolateLog(rho_interpolated, rho_resampled, n_data);

    ApplyFilter(rho_filtered, rho_interpolated, n_data, dt, max_hz_background);
    ResampleLog(rho_resampled, rho_filtered, time_resampled, z_pos, n_data, dt);
    InterpolateLog(rho_background_resolution, rho_resampled, n_data);

    ApplyFilter(rho_filtered, rho_interpolated, n_data, dt, max_hz_seismic);
    ResampleLog(rho_resampled, rho_filtered, time_resampled, z_pos, n_data, dt);
    InterpolateLog(rho_seismic_resolution, rho_resampled, n_data);
  }

  well.AddContLogSeismicResolution("Vp",  vp_seismic_resolution);
  well.AddContLogSeismicResolution("Vs",  vs_seismic_resolution);
  well.AddContLogSeismicResolution("Rho", rho_seismic_resolution);

  well.AddContLogBackgroundResolution("Vp",  vp_background_resolution);
  well.AddContLogBackgroundResolution("Vs",  vs_background_resolution);
  well.AddContLogBackgroundResolution("Rho", rho_background_resolution);

}

//----------------------------------------------------------------------------
bool CommonData::ResampleTime(std::vector<double>       & time_resampled,
                              const std::vector<double> & z_pos,
                              int                         nd,
                              double                    & dt) {
  //Only resample if monotonous increasing in time.
  double time_begin = z_pos[0];
  double time_end   = z_pos[nd - 1];
  bool   monotonous = true;
  for (int i = 1 ; (i < nd && monotonous == true); i++)
    if (z_pos[i] < z_pos[i-1])
      monotonous = false;

  if (monotonous == false)
    return(false);

  //
  // Make new time scale with constant sampling density
  //

  if (time_begin != RMISSING && time_end != RMISSING) {
    dt = (time_end - time_begin)/(nd - 1);            // average sampling density
    for (int i = 0; i < nd; i++)
      time_resampled[i] = time_begin + i*dt;
  }
  else {
    LogKit::LogFormatted(LogKit::Warning,"WARNING: First or last time sample is undefined. Cannot estimate average sampling density.\n");
    LogKit::LogFormatted(LogKit::Warning,"         time[first] = %12.2f\n",time_begin);
    LogKit::LogFormatted(LogKit::Warning,"         time[last]  = %12.2f\n",time_end);
    return(false);
  }

  return(true);

  //printf("i time[i] dt  time_resampled[i] dt   %d  %7.3f         %7.3f\n",0,time[0],time_resampled[0]);
  //for (unsigned int i = 1 ; i < nd ; i++)
  //{
  //  printf("i time[i] dt  time_resampled[i] dt   %d  %7.3f %.3f   %7.3f %.3f\n",i,time[i],time[i]-time[i-1],time_resampled[i],dt);
  //}
}

//----------------------------------------------------------------------------
void CommonData::ResampleLog(std::vector<double>       & log_resampled,
                             const std::vector<double> & log,
                             const std::vector<double> & time,
                             const std::vector<double> & time_resampled,
                             int                         nd,
                             double                      dt) {
  bool resample_log = true;

  if (!resample_log) {
    for (int i = 0; i < nd; i++) {
      log_resampled[i] = log[i];
    }
  }

  //
  // Initialise as undefined
  //
  for (int i = 0 ; i < nd ; i++) {
    log_resampled[i] = RMISSING;
  }

  //
  // Set end points equal to original log. Set the rest undefined
  //
  log_resampled[0]    = log[0];
  log_resampled[nd-1] = log[nd-1];

  //
  // Paals first-shot interpolation...
  //
  // Aritmethic mean of raw values in interval (i-0.5)*dt < t < (i+0.5)*dt
  //
  int j = 0;
  while (time[j] < time_resampled[1] - 0.5*dt) // Find starting position
    j++;

  for (int i = 1; i < nd - 1; i++) { // End points are already set
    // Start gathering values
    double       value = 0.0;
    unsigned int count = 0;

    while ((time[j] < time_resampled[i] + 0.5*dt) && (j < nd-1)) {
      if (log[j] != RMISSING) {
        value += log[j];
        count++;
        //printf("i j value   %d %d  %7.3f\n",i,j,value);

      }
      j++;
    }
    if (count > 0)
      log_resampled[i] = value/count;
  }

  //for (unsigned int i = 0 ; i < nd ; i++)
  //{
  //  printf("i log[i] log_resampled[i]    %d  %7.3f   %7.3f\n",i,log[i],log_resampled[i]);
  //}
}

//----------------------------------------------------------------------------
void CommonData::InterpolateLog(std::vector<double>       & log_interpolated,
                                const std::vector<double> & log_resampled,
                                int                         nd) {
  int i;
  for (i = 0; i < nd; i++) {
    log_interpolated[i] = log_resampled[i];
  }

  // Skip leading RMISSING

  i = 0;
  while (i < nd && log_resampled[i] == RMISSING)
    i++;

  // When the log has intermediate RMISSING... use linear interpolation. Skip trailing RMISSING

  while (i < nd) {
    if (log_resampled[i] == RMISSING) {
      int last_nonmissing = i - 1;                // last defined value (e.g. ..., 2.31, -99999, -99999, ...)
      while (i < nd && log_resampled[i] == RMISSING)  //                               ^^^^
        i++;
      if (i < nd) {
        int   first_nonmissing = i;               // first defined value (e.g. ..., -99999, -99999, 2.29, ...)
        int   j0 = last_nonmissing;               //                                                ^^^^
        int   j1 = first_nonmissing;
        double l0 = log_resampled[j0];
        double l1 = log_resampled[j1];

        double a = (l1 - l0)/double(j1 - j0);

        for (int j = j0 + 1; j < j1; j++) {
          log_interpolated[j] = a*(j - j0) + l0;
        }
      }
    }
    else
      i++;
  }
  //for (int i=0 ; i<nd ; i++) {
  //  printf("i log_resampled[i] log_interpolated[i]  %d  %.3f  %.3f\n",i,log_resampled[i],log_interpolated[i]);
  //}
}

//----------------------------------------------------------------------------
void CommonData::ApplyFilter(std::vector<double> & log_filtered,
                             std::vector<double> & log_interpolated,
                             int                   n_time_samples,
                             double                dt_milliseconds,
                             float                 max_hz) {
  //
  // Extract nonmissing part of log
  //
  int i=0;
  while (i < n_time_samples && log_interpolated[i] == RMISSING)
    i++;
  int first_nonmissing = i;
  i = n_time_samples - 1;
  while (i > 0 && log_interpolated[i] == RMISSING)
    i--;
  int last_nonmissing = i;
  int n_time_samples_defined = last_nonmissing - first_nonmissing + 1;

  for (i = 0; i < n_time_samples; i++) {            // Initialise with RMISSING
    log_filtered[i] = RMISSING;
  }

  if (n_time_samples_defined > 0) {
    //
    // Setup FFT
    //
    int   nt  = 2*n_time_samples_defined;
    int   cnt = nt/2 + 1;
    int   rnt = 2*cnt;

    fftw_real*    rAmp = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
    fftw_complex* cAmp = reinterpret_cast<fftw_complex*>(rAmp);

    for (i = 0; i < n_time_samples_defined; i++) {          // Array to filter is made symmetric
      rAmp[i]      = static_cast<fftw_real>(log_interpolated[first_nonmissing + i]);
      rAmp[nt-i-1] = rAmp[i];
    }

    //for (int i=0 ; i<nt ; i++) {
    //  printf("i=%d, log_interpolated[i]=%7.4f\n",i,rAmp[i]);
    //}

    //
    // Transform to Fourier domain
    //
    rfftwnd_plan p1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
    rfftwnd_one_real_to_complex(p1, rAmp, cAmp);
    fftwnd_destroy_plan(p1);

    //for (int i=0 ; i<cnt ; i++) {
    //  printf("i=%2d, cAmp.re[i]=%11.4f  cAmp.im[i]=%11.4f\n",i,cAmp[i].re,cAmp[i].im);
    //}

    //
    // Filter using Odd's magic vector...
    //
    float dt  = static_cast<float> (dt_milliseconds/1000.0); // Sampling density in seconds
    float T   = (nt - 1)*dt;                                 // Time sample
    float w   = 1/T;                                         // Lowest frequency that can be extracted from log
    int   N   = int(max_hz/w + 0.5f);                        // Number of elements of Fourier vector to keep

    if (cnt < N+1) {
      LogKit::LogMessage(LogKit::Warning, "Warning: The vertical resolution is too low to allow filtering of well logs to %3.1f Hz.\n");
    }

    float * magic_vector = new float[cnt];
    for (i = 0; ((i < N+1) && (i < cnt)); i++) {
      magic_vector[i] = 1.0;
    }
    for (;i < cnt; i++) {
      magic_vector[i] = 0.0;
    }
    for (i = 0; i < cnt; i++) {
      cAmp[i].re *= magic_vector[i];
      cAmp[i].im *= magic_vector[i];
    }

    //for (int i=0 ; i<cnt ; i++) {
    //  printf("i=%2d, cAmp.re[i]=%11.4f  cAmp.im[i]=%11.4f\n",i,cAmp[i].re,cAmp[i].im);
    //}

    //
    // Backtransform to time domain
    //
    rfftwnd_plan p2 = rfftwnd_create_plan(1, &nt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
    rfftwnd_one_complex_to_real(p2, cAmp, rAmp);
    fftwnd_destroy_plan(p2);

    float scale= float(1.0/nt);
    for (i = 0; i < rnt; i++) {
      rAmp[i] *= scale;
    }

    //
    // Fill log_filtered[]
    //
    for (i = 0; i < n_time_samples_defined; i++) {
      log_filtered[first_nonmissing + i] = rAmp[i];      // Fill with values where defined
    }
    delete [] magic_vector;
    fftw_free(rAmp);

    //for (int i=0 ; i<n_time_samples ; i++) {
    //  printf("i log_interpolated[i] log_filtered[i]  %d  %.3f  %.3f\n",i,log_interpolated[i],log_filtered[i]);
    //}
  }
}


void CommonData::CutWell(std::string           well_file_name,
                         NRLib::Well         & well,
                         const NRLib::Volume & full_inversion_volume) {

  //This is run after ProcessLogsNorsarWell and ProcessLogsRMSWell so log names should be equal
  //Possible Logs: X_pos, Y_pos, Z_pos, Dt, Vp, Dts, Vs, Rho, Facies

  const std::vector<double> & x_old = well.GetContLog("X_pos");
  const std::vector<double> & y_old = well.GetContLog("Y_pos");
  const std::vector<double> & z_old = well.GetContLog("Z_pos");

  std::vector<double> x_new;
  std::vector<double> y_new;
  std::vector<double> z_new;
  std::vector<double> dt_new;
  std::vector<double> vp_new;
  std::vector<double> dts_new;
  std::vector<double> vs_new;
  std::vector<double> rho_new;
  std::vector<int>    facies_new;

  const NRLib::Surface<double> & top_surf = full_inversion_volume.GetTopSurface();
  const NRLib::Surface<double> & bot_surf = full_inversion_volume.GetBotSurface();

  double top_surf_min = top_surf.Min();
  double top_surf_max = top_surf.Max();
  double bot_surf_min = bot_surf.Min();
  double bot_surf_max = bot_surf.Max();

  bool need_cutting = false;
  //if (z_old[0] < top_surf.GetZ(x_old[0], y_old[0]) || z_old[z_old.size()-1] > bot_surf.GetZ(x_old[z_old.size()-1], y_old[z_old.size()-1]))
  if (z_old[0] < top_surf_max || z_old[z_old.size()-1] > bot_surf_min)
    need_cutting = true;

  if (need_cutting == true) { //H Cut missing values?

    for (size_t i = 0; i < z_old.size(); i++) {

      if (z_old[i] != RMISSING) {

        double top_tmp = top_surf.GetZ(x_old[i], y_old[i]);
        if (top_tmp = WELLMISSING)
          top_tmp = top_surf_min;

        double bot_tmp = bot_surf.GetZ(x_old[i], y_old[i]);
        if (bot_tmp = WELLMISSING)
          bot_tmp = bot_surf_max;

        //if ( !(z_old[i] < top_surf.GetZ(x_old[i], y_old[i]) || z_old[i] > bot_surf.GetZ(x_old[i], y_old[i])) ) { //Inside, keep
        if ( !(z_old[i] < top_tmp || z_old[i] > bot_tmp) ) { //Inside, keep

          x_new.push_back(x_old[i]);
          y_new.push_back(y_old[i]);
          z_new.push_back(z_old[i]);

          if (well.HasContLog("Dt"))
            dt_new.push_back(well.GetContLog("Dt")[i]);
          if (well.HasContLog("Vp"))
            vp_new.push_back(well.GetContLog("Vp")[i]);
          if (well.HasContLog("Dts"))
            dts_new.push_back(well.GetContLog("Dts")[i]);
          if (well.HasContLog("Vs"))
            vs_new.push_back(well.GetContLog("Vs")[i]);
          if (well.HasContLog("Rho"))
            rho_new.push_back(well.GetContLog("Rho")[i]);
          if (well.HasDiscLog("Facies"))
            facies_new.push_back(well.GetDiscLog("Facies")[i]);
        }
      }
    }

    //Replace logs
    well.RemoveContLog("X_pos");
    well.AddContLog("X_pos", x_new);
    well.RemoveContLog("Y_pos");
    well.AddContLog("Y_pos", y_new);
    well.RemoveContLog("Z_pos");
    well.AddContLog("Z_pos", z_new);
    if (well.HasContLog("Dt")) {
      well.RemoveContLog("Dt");
      well.AddContLog("Dt", dt_new);
    }
    if (well.HasContLog("Vp")) {
      well.RemoveContLog("Vp");
      well.AddContLog("Vp", vp_new);
    }
    if (well.HasContLog("Dts")) {
      well.RemoveContLog("Dts");
      well.AddContLog("Dts", dts_new);
    }
    if (well.HasContLog("Vs")) {
      well.RemoveContLog("Vs");
      well.AddContLog("Vs", vs_new);
    }
    if (well.HasContLog("Rho")) {
      well.RemoveContLog("Rho");
      well.AddContLog("Rho", rho_new);
    }
    if (well.HasDiscLog("Facies")) {
      well.RemoveDiscLog("Facies");
      well.AddDiscLog("Facies", facies_new);
    }

    well.SetNumberOfData(z_new.size());

    int nonmissing_data = 0;
    for (size_t i = 0; i < z_new.size(); i++) {
      if (z_new[i] != WELLMISSING)
        nonmissing_data++;
    }
    well.SetNumberOfNonMissingData(nonmissing_data);

  }
}


void CommonData::ProcessLogsNorsarWell(NRLib::Well              & new_well,
                                       std::vector<std::string> & log_names_from_user,
                                       const std::vector<bool>  & inverse_velocity,
                                       bool                       facies_log_given,
                                       std::string              & err_text) {

  const int factor_kilometer = 1000;

  if (log_names_from_user.size() == 0) {
    log_names_from_user.push_back("TWT");
    log_names_from_user.push_back("DT");
    log_names_from_user.push_back("RHOB");
    log_names_from_user.push_back("DTS");
    log_names_from_user.push_back("FACIES");
  }

  for (size_t i = 0; i < log_names_from_user.size(); i++) {
    // If the well does not contain the log specified by the user, return an error message
    if (!new_well.HasContLog(log_names_from_user[i]) && !new_well.HasDiscLog(log_names_from_user[i])){
      err_text+="Could not find log \'" + log_names_from_user[i] + "\' in well file \'"+new_well.GetWellName()+"\'.\n";
    }
  }

  // Norsar wells must have a UTMX log
  if (new_well.HasContLog("UTMX")) {
    std::vector<double> x_pos_temp = new_well.GetContLog("UTMX");
    for (unsigned int i=0; i < x_pos_temp.size(); i++) {
      x_pos_temp[i] = x_pos_temp[i]*factor_kilometer;
    }
    new_well.AddContLog("X_pos", x_pos_temp);
    new_well.RemoveContLog("UTMX");
  }
  else {
    err_text += "Could not find log 'UTMX' in well file "+new_well.GetWellName()+".\n";
  }

  // Norsar wells must have a UTMY log
  if (new_well.HasContLog("UTMY")) {
    std::vector<double> y_pos_temp = new_well.GetContLog("UTMY");
    for (unsigned int i = 0; i < y_pos_temp.size(); i++) {
      y_pos_temp[i] = y_pos_temp[i]*factor_kilometer;
    }
    new_well.AddContLog("Y_pos", y_pos_temp);
    new_well.RemoveContLog("UTMY");
  }
  else {
    err_text += "Could not find log 'UTMY' in well file "+new_well.GetWellName()+".\n";
  }

  int nonmissing_data = 0; //Count number of data not Missing (nd_ in welldata.h)

  //Store TWT as Z_pos
  if (new_well.HasContLog("TWT")) {
    std::vector<double> twt_temp = new_well.GetContLog("TWT");
    for (unsigned int i = 0; i < twt_temp.size(); i++) {
      twt_temp[i] = twt_temp[i]*factor_kilometer;

      if (twt_temp[i] != WELLMISSING)
        nonmissing_data++;

    }
    new_well.AddContLog("Z_pos", twt_temp);

    new_well.SetNumberOfNonMissingData(nonmissing_data);
  }
  else { // Process MD log if TVD is not available?
    err_text += "Could not find log 'TWT' in well file "+new_well.GetWellName()+".\n";
  }

  // Time is always entry 0 in the log name list and is always called TWT
  //if (new_well.HasContLog(log_names_from_user[0])){

  double well_missing = new_well.GetContMissing();

  // Vp/Dt is always entry 1 in the log name list
  if (new_well.HasContLog(log_names_from_user[1])) {

    //Change missing values. Norsar wells has its own missing value from the well-log
    std::vector<double> vp_dt_tmp = new_well.GetContLog(log_names_from_user[1]);
    for (size_t i = 0; i < vp_dt_tmp.size(); i++) {
      if (vp_dt_tmp[i] == well_missing)
        vp_dt_tmp[i] = RMISSING;
    }

    if (inverse_velocity[0] == true) {
      new_well.RemoveContLog(log_names_from_user[1]);
      new_well.AddContLog("Dt", vp_dt_tmp); //new_well.GetContLog(log_names_from_user[1])
    }
    else {
      new_well.RemoveContLog(log_names_from_user[1]);
      new_well.AddContLog("Vp", vp_dt_tmp);
    }
  }

  // Rho is always entry 2 in the log name list
  if (new_well.HasContLog(log_names_from_user[2])) {

    std::vector<double> rho_tmp = new_well.GetContLog(log_names_from_user[2]);
    for (size_t i = 0; i < rho_tmp.size(); i++) {
      if (rho_tmp[i] == well_missing)
        rho_tmp[i] = RMISSING;
    }

    new_well.RemoveContLog(log_names_from_user[2]);
    new_well.AddContLog("Rho", rho_tmp);
  }

  // Vs/Dts is always entry 3 in the log name list
  if (new_well.HasContLog(log_names_from_user[3])) {

    std::vector<double> vs_dts_tmp = new_well.GetContLog(log_names_from_user[3]);
    for (size_t i = 0; i < vs_dts_tmp.size(); i++) {
      if (vs_dts_tmp[i] == well_missing)
        vs_dts_tmp[i] = RMISSING;
    }

    if (inverse_velocity[1] == true) {
      new_well.RemoveContLog(log_names_from_user[3]);
      new_well.AddContLog("Dts", vs_dts_tmp);
    }
    else {
      new_well.RemoveContLog(log_names_from_user[3]);
      new_well.AddContLog("Vs", vs_dts_tmp);
    }
  }

  // If the facies log is given, it is always entry 4 in the log name list
  if (facies_log_given) {
    int well_missing_int = new_well.GetIntMissing();

    if (new_well.HasDiscLog(log_names_from_user[4])) {

      std::vector<int> & facies_tmp = new_well.GetDiscLog(log_names_from_user[4]);
      for (size_t i = 0; i < facies_tmp.size(); i++) {
        if (facies_tmp[i] == well_missing_int)
          facies_tmp[i] = IMISSING;
      }

      new_well.RemoveDiscLog(log_names_from_user[4]);
      new_well.AddDiscLog("Facies", facies_tmp);
    }
  }

}

void CommonData::ProcessLogsRMSWell(NRLib::Well                     & new_well,
                                    std::vector<std::string>        & log_names_from_user,
                                    const std::vector<bool>         & inverse_velocity,
                                    bool                              facies_log_given,
                                    std::string                     & error_text) {

  const double factor_usfeet_to_meters = 304800.0;

  if (log_names_from_user.size() == 0) {
    log_names_from_user.push_back("TWT");
    log_names_from_user.push_back("DT");
    log_names_from_user.push_back("RHOB");
    log_names_from_user.push_back("DTS");
    log_names_from_user.push_back("FACIES");
  }

  for (size_t i = 0; i < log_names_from_user.size(); i++) {
    // If the well does not contain the log specified by the user, return an error message
    if (!new_well.HasContLog(log_names_from_user[i]) && !new_well.HasDiscLog(log_names_from_user[i])) {
      error_text+="Could not find log \'" + log_names_from_user[i] + "\' in well file \'"+new_well.GetWellName()+"\'.\n";
    }
  }

  // RMS wells must have an x log
  if (new_well.HasContLog("X")) {
    new_well.AddContLog("X_pos", new_well.GetContLog("X"));
    new_well.RemoveContLog("X");
  }
  else {
    error_text += "Could not find log 'X' in well file "+new_well.GetWellName()+".\n";
  }

  // RMS wells must have a y log
  if (new_well.HasContLog("Y")) {
    new_well.AddContLog("Y_pos", new_well.GetContLog("Y"));
    new_well.RemoveContLog("Y");
  }
  else {
    error_text += "Could not find log 'Y' in well file "+new_well.GetWellName()+".\n";
  }

  //Store TWT as Z_pos
  if (new_well.HasContLog("Z"))
    new_well.RemoveContLog("Z");
  if (new_well.HasContLog("TWT")) {
    new_well.AddContLog("Z_pos", new_well.GetContLog("TWT"));
  }
  else {
    error_text += "Could not find log 'TWT' in well file "+new_well.GetWellName()+".\n";
  }

  int nonmissing_data = 0; //Count number of data not Missing (nd_ in welldata.h)
  const std::vector<double> & z_tmp = new_well.GetContLog("Z_pos");
  for (size_t i = 0; i < z_tmp.size(); i++) {
    if (z_tmp[i] != WELLMISSING)
      nonmissing_data++;
  }
  new_well.SetNumberOfNonMissingData(nonmissing_data);

  // Time is always entry 0 in the log name list
  // Time is always called 'TWT', so no need to rename it
  //if (new_well.HasContLog(log_names_from_user[0]))

  // Vp is always entry 1 in the log name list and entry 0 in inverse_velocity
  if (new_well.HasContLog(log_names_from_user[1])) {
    const std::vector<double> & vp_temp = new_well.GetContLog(log_names_from_user[1]);
    std::vector<double> vp(vp_temp.size());
    if (inverse_velocity[0]) {
      for (unsigned int i=0; i<vp_temp.size(); i++) {
        if (vp_temp[i] != WELLMISSING) //H Added check if vp-log had missing. Later it is checked against RMISSING
          vp[i] = static_cast<double>(factor_usfeet_to_meters/vp_temp[i]);
        else
          vp[i] = RMISSING;
      }
    }
    else {
      for (unsigned int i=0; i<vp_temp.size(); i++) {
        if (vp_temp[i] != WELLMISSING)
          vp[i] = static_cast<double>(vp_temp[i]);
        else
          vp[i] = RMISSING;
      }
    }
    new_well.RemoveContLog(log_names_from_user[1]);
    new_well.AddContLog("Vp", vp);
  }

  // Vs is always entry 3 in the log name list
  if (new_well.HasContLog(log_names_from_user[3])) {
    const std::vector<double> & vs_temp = new_well.GetContLog(log_names_from_user[3]);
    std::vector<double> vs(vs_temp.size());
    if (inverse_velocity[1]) {
      for (unsigned int i = 0; i < vs_temp.size(); i++) {
        if (vs_temp[i] != WELLMISSING)
          vs[i] = static_cast<double>(factor_usfeet_to_meters/vs_temp[i]);
        else
          vs[i] = RMISSING;
      }
    }
    else {
      for (unsigned int i=0; i<vs_temp.size(); i++) {
        if (vs_temp[i] != WELLMISSING)
          vs[i] = static_cast<double>(vs_temp[i]);
        else
          vs[i] = RMISSING;
      }
    }
    new_well.RemoveContLog(log_names_from_user[3]);
    new_well.AddContLog("Vs", vs);
  }

  // Rho is always entry 2 in the log name list
  if (new_well.HasContLog(log_names_from_user[2])) {
    new_well.AddContLog("Rho", new_well.GetContLog(log_names_from_user[2]));
    new_well.RemoveContLog(log_names_from_user[2]);
  }

  // If defined, Facies is always entry 4 in the log name list
  if (facies_log_given) {
    if (new_well.HasDiscLog(log_names_from_user[4])) {
      new_well.AddDiscLog("Facies", new_well.GetDiscLog(log_names_from_user[4]));
      new_well.RemoveDiscLog(log_names_from_user[4]);
    }
  }
}

//void
//CommonData::ReadFaciesNamesFromWellFile(ModelSettings            * model_settings,
//                                        std::string                well_file_name,
//                                        std::vector<int>         & facies_nr,
//                                        std::vector<std::string> & facies_names,
//                                        std::string              & err_txt)
//{
//  std::ifstream file;
//  std::string token;
//  std::string dummy_str;
//  std::string well_name;
//  NRLib::OpenRead(file, well_file_name);
//
//  int nlog; // number of logs in file
//  int line = 0;
//  NRLib::DiscardRestOfLine(file,line,false); //First two lines contain info we do not need.
//  NRLib::DiscardRestOfLine(file,line,false);
//  NRLib::ReadNextToken(file, token, line);
//  well_name = token;
//  NRLib::DiscardRestOfLine(file,line,false); //xpos, ypos.
//
//  //facies_wells_.push_back(well_name);
//
//  nlog   = NRLib::ReadNext<int>(file, line);
//
//  int nVar = 5;       // z,alpha,beta,rho, and facies
//  std::vector<std::string> log_names = model_settings->getLogNames();
//
//  std::vector<std::string> parameterList(5);
//  if (log_names[0] != "") // Assume that all lognames are filled present if first is.
//  {
//    parameterList = log_names;
//  }
//  else
//  {
//  parameterList[0] = "TWT";
//  parameterList[1] = "DT";
//  parameterList[2] = "RHOB";
//  parameterList[3] = "DTS";
//  parameterList[4] = "FACIES";
//  }
//
//  int * pos = new int[nVar];
//  for (int i=0;i<nVar;i++)
//    pos[i] = IMISSING;
//
//  std::string facies_log_name;
//
//  //Find number of facies
//  int n_facies = 0;
//  for (int i=0; i < nlog; i++)
//  {
//    NRLib::ReadNextToken(file,token,line);
//    for (int j=0; j < nVar; j++)
//    {
//      if ( NRLib::Uppercase(token)==parameterList[j])
//      {
//        pos[j] = i + 4;
//        if (j==4)
//        {
//          facies_log_name = parameterList[4];
//          // facies log - save names
//          NRLib::ReadNextToken(file,token,line); // read code word DISC
//          if (token != "DISC")
//          {
//            err_txt += "Facies log must be discrete for well " + well_file_name + ".\n";
//            return;
//            //LogKit::LogFormatted(LogKit::Error,"ERROR: Facies log must be discrete.\n");
//          }
//          // Find number of facies
//          std::getline(file, dummy_str);
//          std::vector<std::string> tokenLine = NRLib::GetTokens(dummy_str);
//          n_facies = static_cast<int>(tokenLine.size())/2;
//        }
//      }
//    }
//    if (token != "DISC")
//      NRLib::DiscardRestOfLine(file,line,false);
//  }
//
//  file.close();
//  file.clear();
//
//  int k;
//  facies_nr.resize(n_facies);
//  facies_names.resize(n_facies);
//  //std::vector<int> facies_nr_tmp(n_facies);
//  //std::vector<std::string> facies_names_tmp(n_facies);
//
//  NRLib::OpenRead(file, well_file_name);
//  line = 0;
//  for (int i=0; i < 4+nlog; i++)
//  {
//    NRLib::ReadNextToken(file,token,line);
//    if (NRLib::Uppercase(token) == parameterList[4])
//    {
//      NRLib::ReadNextToken(file,token,line); // read code word DISC
//      // facies types given here
//      for (k=0; k < n_facies; k++)
//      {
//        NRLib::ReadNextToken(file,token,line);
//        facies_nr[k] = NRLib::ParseType<int>(token);
//        NRLib::ReadNextToken(file,token,line);
//        facies_names[i] = token;
//      }
//    }
//    NRLib::DiscardRestOfLine(file,line,false);
//  }
//
//  //facies_nr_wells_.push_back(facies_nr_tmp);
//  //facies_names_wells_.push_back(facies_names_tmp);
//}

void
CommonData::ReadFaciesNamesFromWellLogs(NRLib::Well              & well,
                                        std::vector<int>         & facies_nr,
                                        std::vector<std::string> & facies_names)
                                        //std::string              & err_txt)
{
  const std::map<int, std::string> & facies_map = well.GetFaciesMap();
  int n_facies = facies_map.size();

  facies_nr.resize(n_facies);
  facies_names.resize(n_facies);

  for (size_t i = 0; i < facies_map.size(); i++) {
    facies_nr[i] = i;
    facies_names[i] = facies_map.find(i)->second;
  }
}

void CommonData::SetFaciesNamesFromWells(ModelSettings            *& model_settings,
                                         std::string               & err_text) {

  int min, max;
  int globalmin = 0;
  int globalmax = 0;
  int n_facies = 0;
  std::vector<int> facies_nr;
  bool first = true;
  for (int w = 0; w < model_settings->getNumberOfWells(); w++) {
    n_facies = facies_names_wells_[w].size();
    facies_nr = facies_nr_wells_[w];

    if (facies_log_wells_[w] == true) {
      GetMinMaxFnr(min,max, n_facies, facies_nr);

      if (first==true) {
        globalmin = min;
        globalmax = max;
        first = false;
      }
      else {
        if (min<globalmin)
          globalmin = min;
        if (max>globalmax)
          globalmax = max;
      }
    }
  }

  int n_names = globalmax - globalmin + 1;
  std::vector<std::string> names(n_names);

  for (int w = 0; w < model_settings->getNumberOfWells(); w++) {

    if (facies_log_wells_[w] == true) {

      n_facies = facies_names_wells_[w].size();

      for (int i = 0; i < n_facies; i++) {

        std::string name = facies_names_wells_[w][i];
        int         fnr  = facies_nr_wells_[w][i] - globalmin;

        if (names[fnr] == "") {
          names[fnr] = name;
        }
        else if (names[fnr] != name) {
          err_text += "Problem with facies logs. Facies names and numbers are not uniquely defined.\n";
        }
      }
    }
  }

  LogKit::LogFormatted(LogKit::Low,"\nFaciesLabel      FaciesName           ");
  LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------\n");
  for (int i=0 ; i < n_names ; i++) {
    if (names[i] != "")
      LogKit::LogFormatted(LogKit::Low,"    %2d           %-20s\n",i+globalmin,names[i].c_str());
  }

  //n_facies = 0;
  //for (int i=0 ; i < n_names ; i++)
  //  if (names[i] != "")
  //    n_facies++;

  for (int i=0 ; i < n_names ; i++) {
    if (names[i] != "") {
      facies_nr_.push_back(globalmin + i);
      facies_names_.push_back(names[i]);
    }
  }
}

void
CommonData::GetMinMaxFnr(int            & min,
                         int            & max,
                         const int        n_facies,
                         std::vector<int> facies_nr)
{
  int i;
  //int premin, premax;
  min = facies_nr[0];
  max = facies_nr[0];
  for (i=1;i<n_facies;i++)
  {
    if (facies_nr[i]<min)
      min = facies_nr[i];
    if (facies_nr[i]>max)
      max = facies_nr[i];
  }
}

bool CommonData::SetupReflectionMatrix(ModelSettings * model_settings,
                                       InputFiles    * input_files,
                                       std::string   & err_text_common) {

  LogKit::WriteHeader("Setting up reflection matrix");

  //If Vp/Vs is given per interval: A default reflection matrix is set up here, and then altered later per interval.

  std::string err_text = "";
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model or wells.
  //
  const std::string & refl_matr_file = input_files->getReflMatrFile();
  const double vpvs                  = model_settings->getVpVsRatio();
  float                             ** reflection_matrix;

  int n_timelapses = model_settings->getNumberOfTimeLapses(); //Returnerer timeLapseAngle_.size()
  for (int i = 0; i < n_timelapses; i++) {

    int n_angles = model_settings->getNumberOfAngles(i);

    if (refl_matr_file != "") { //File should have one line for each seismic data file. Check: if (input_files->getNumberOfSeismicFiles(thisTimeLapse) > 0 ) ?
      std::string tmp_err_text("");
      reflection_matrix = ReadMatrix(refl_matr_file, n_angles, 3, "reflection matrix", tmp_err_text);
      if (reflection_matrix == NULL) {
        err_text += "Reading of file "+refl_matr_file+ " for reflection matrix failed\n";
        err_text += tmp_err_text;
      }
      else
        LogKit::LogFormatted(LogKit::Low,"\nReflection parameters read from file.\n\n");

      refmat_from_file_global_vpvs_ = true;
    }

    else if (vpvs != RMISSING) {
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs ratio specified in model file.\n");
      double vsvp = 1.0/vpvs;
      SetupDefaultReflectionMatrix(reflection_matrix, vsvp, model_settings, n_angles, i);
      refmat_from_file_global_vpvs_ = true;
    }
    else {
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs equal to 2\n");
      double vsvp = 0.5;
      SetupDefaultReflectionMatrix(reflection_matrix, vsvp, model_settings, n_angles, i);
      refmat_from_file_global_vpvs_ = false;
    }

    reflection_matrix_[i] = reflection_matrix;

  } //nTimeLapses

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

bool CommonData::SetupTemporaryWavelet(ModelSettings * model_settings,
                                       //InputFiles    * input_files,
                                       std::string   & err_text_common) {
  //Set up temporary wavelet
  LogKit::WriteHeader("Setting up temporary wavelet");

  std::string err_text = "";
  int n_timelapses     = model_settings->getNumberOfTimeLapses();

  //3. Use Ricker - wavelet.
  //4. 1 wavelet per angle
  //5 Frequency per ange: Take 100 traces from first AVO-vintage on this angle. Find peak-frequency for these.
  int n_angles = model_settings->getNumberOfAngles(0);
  //int error = 0;

  for (int j = 0 ; j < n_angles; j++) {
    //Check all timelapses for this angle, choose the first one;
    int this_timelapse = 0;
    int vintage_year   = model_settings->getVintageYear(0);
    int vintage_month  = model_settings->getVintageMonth(0);
    int vintage_day    = model_settings->getVintageDay(0);
    for (int k = 1; k < n_timelapses; k++) {
      if (model_settings->getVintageYear(k) <= vintage_year && model_settings->getVintageMonth(k) <= vintage_month && model_settings->getVintageDay(k) <= vintage_day) {
        vintage_year   = model_settings->getVintageYear(k);
        vintage_month  = model_settings->getVintageMonth(k);
        vintage_day    = model_settings->getVintageDay(k);
        this_timelapse = k;
      }
    }
    std::vector<float> angles = model_settings->getAngle(this_timelapse);

    std::vector<float> frequency_peaks;
    std::vector<std::vector<float> > trace_data(100);
    std::vector<float> trace_length;
    seismic_data_[this_timelapse][j].GetSparseTraceData(trace_data, trace_length, 100);

    //FFT to find peak-frequency.
    for (int k = 0; k < 100; k++) {
      int n_trace     = trace_data[k].size();
      int n_trace_fft = ((n_trace / 2) + 1)*2;

      fftw_real * seis_r    = new fftw_real[n_trace_fft];
      fftw_complex * seis_c = reinterpret_cast<fftw_complex*>(seis_r);

      for (int kk=0; kk < n_trace; kk++)
        seis_r[kk] = trace_data[k][kk];

      Utils::fft(seis_r, seis_c, n_trace);

      float peak_tmp = 0.0;
      float length_tmp = 0.0;
      for (int kk = 0; kk < n_trace; kk++) {
        //fftw_complex value = seis_c[kk];
        length_tmp = std::sqrt( (seis_c[kk].re)*(seis_c[kk].re) + (seis_c[kk].im)*(seis_c[kk].im) );
        if (length_tmp > peak_tmp && length_tmp != std::numeric_limits<float>::infinity())
            peak_tmp = length_tmp;
      }

      peak_tmp /= trace_length[k];
      frequency_peaks.push_back(peak_tmp);

      if (seis_r != NULL)
        delete seis_r;
      //if (seis_c != NULL)
      // delete seis_c;
    }

    float mean_frequency   = 0.0;
    float mean_grid_height = 0.0;
    for (size_t k = 0; k < frequency_peaks.size(); k++) {
      mean_frequency   += frequency_peaks[k];
      mean_grid_height += trace_length[k];
    }
    mean_frequency   /= frequency_peaks.size();
    mean_grid_height /= frequency_peaks.size(); //Will be 1 if SEGY is used.

    mean_frequency   /= mean_grid_height;

    int tmp_error = 0;
    Wavelet * wavelet_tmp = new Wavelet1D(model_settings, reflection_matrix_[this_timelapse][j], angles[j], mean_frequency, tmp_error);

    if (tmp_error == 0)
      temporary_wavelets_.push_back(wavelet_tmp);
    else
      err_text += "Error setting up a temporary wavelet for angle " + NRLib::ToString(angles[j]) + ".\n";
  }

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

float ** CommonData::ReadMatrix(const std::string & file_name,
                                int                 n1,
                                int                 n2,
                                const std::string & read_reason,
                                std::string       & err_text)
{
  float * tmp_res = new float[n1*n2+1];
  std::ifstream in_file;
  NRLib::OpenRead(in_file, file_name);
  std::string text = "Reading "+read_reason+" from file "+file_name+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
  std::string storage;
  int index   = 0;
  bool failed = false;

  while(failed == false && in_file >> storage) {
    if (index < n1*n2) {
      try {
        tmp_res[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        err_text += "Error in "+file_name+"\n";
        err_text += e.what();
        failed = true;
      }
    }
    index++;
  }
  if (failed == false) {
    if (index != n1*n2) {
      failed = true;
      err_text += "Found "+NRLib::ToString(index)+" in file "+file_name+", expected "+NRLib::ToString(n1*n2)+".\n";
    }
  }

  float ** result = NULL;
  if (failed == false) {
    LogKit::LogFormatted(LogKit::Low,"ok.\n");
    result = new float * [n1];
    int i, j;
    index = 0;
    for (i=0;i<n1;i++) {
      result[i] = new float[n2];
      for (j=0;j<n2;j++) {
        result[i][j] = tmp_res[index];
        index++;
      }
    }
  }
  else
    err_text += "Reading matrix " + file_name + "failed.\n";
    //LogKit::LogFormatted(LogKit::Low,"failed.\n");

  delete [] tmp_res;
  return(result);
}


void
CommonData::SetupDefaultReflectionMatrix(float             **& reflection_matrix,
                                         double                vsvp,
                                         const ModelSettings * model_settings,
                                         int                   n_angles,
                                         int                   this_timelapse)
{
  int i;
  float ** A = new float * [n_angles];

  double vsvp2 = vsvp*vsvp;
  std::vector<int> seismic_type = model_settings->getSeismicType(this_timelapse);
  std::vector<float> angles = model_settings->getAngle(this_timelapse);

  for (i = 0; i < n_angles; i++)
  {
    double angle = static_cast<double>(angles[i]);
    A[i] = new float[3];
    double sint = sin(angle);
    double sint2 = sint*sint;
    if (seismic_type[i] == ModelSettings::STANDARDSEIS) { //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float( (1.0 +tan2t )/2.0 ) ;
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0 );
    }
    else if (seismic_type[i] == ModelSettings::PSSEIS) {
      double cost = cos(angle);
      double cosp = sqrt(1-vsvp2*sint2);
      double fac = 0.5*sint/cosp;

      A[i][0] = 0;
      A[i][1] = float(4.0*fac*(vsvp2*sint2-vsvp*cost*cosp));
      A[i][2] = float(fac*(-1.0+2*vsvp2*sint2+2*vsvp*cost*cosp));
    }
  }
  reflection_matrix = A;
  double vpvs = 1.0f/vsvp;
  LogKit::LogFormatted(LogKit::Low,"\nMaking reflection parameters using a Vp/Vs ratio of %4.2f\n",vpvs);
  std::string text;
  if (vpvs < model_settings->getVpVsRatioMin()) {
    LogKit::LogFormatted(LogKit::Warning,"\nA very small Vp/Vs-ratio has been detected. Values below %.2f are regarded unlikely. \n",model_settings->getVpVsRatioMin());
    text = "Check the Vp/Vs-ratio. A small value has been found. If the value is acceptable,\n";
    text += " you can remove this task using the <minimim-vp-vs-ratio> keyword.\n";
    TaskList::addTask(text);
  }
  else if (vpvs > model_settings->getVpVsRatioMax()) {
    LogKit::LogFormatted(LogKit::Warning,"\nA very large Vp/Vs-ratio has been detected. Values above %.2f are regarded unlikely. \n",model_settings->getVpVsRatioMax());
    text = "Check the Vp/Vs-ratio. A large value has been found. If the value is acceptable,\n";
    text += " you can remove this task using the <maximum-vp-vs-ratio> keyword.\n";
    TaskList::addTask(text);
  }
}

bool CommonData::WaveletHandling(ModelSettings * model_settings,
                                 InputFiles    * input_files,
                                 std::string   & err_text_common) {

  int n_timeLapses     = model_settings->getNumberOfTimeLapses();
  int error            = 0;
  std::string err_text = "";

  std::string err_text_tmp("");
  std::vector<Surface *> wavelet_estim_interval;
  FindWaveletEstimationInterval(input_files, wavelet_estim_interval, err_text_tmp);

  if (err_text_tmp != "")
    err_text += "Error when finding wavelet estimation interval: " + err_text_tmp + "\n";

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  for (int i = 0; i < n_timeLapses; i++) {

    int n_angles = model_settings->getNumberOfAngles(i);

    std::vector<float> sn_ratio = model_settings->getSNRatio(i);
    std::vector<float> angles   = model_settings->getAngle(i);
    use_local_noise_            = model_settings->getUseLocalNoise(i);

    //Fra ModelAvoDynamic::processSeismic:
    std::vector<bool> estimate_wavelets = model_settings->getEstimateWavelet(i);

    //Estimation of a wavelet requires the reading of seismic, reading of wells and reflection matrix to be ok. Check blocking of wells since they are used in estimation.
    bool estimate_failed = false;
    for (size_t j = 0; j < estimate_wavelets.size(); j++) {
      if (estimate_wavelets[j]) {
        if (read_seismic_ == false || read_wells_ == false || setup_reflection_matrix_ == false || block_wells_ == false) {
          estimate_failed = true;
          return false;
        }
      }
    }

    std::vector<bool> use_ricker_wavelet = model_settings->getUseRickerWavelet(i);

    std::vector<Wavelet *> wavelet(n_angles);

    std::vector<Grid2D *> local_noise_scale; ///< Scale factors for local noise
    std::vector<Grid2D *> local_shift;
    std::vector<Grid2D *> local_scale;
    local_noise_scale.resize(n_angles);
    local_shift.resize(n_angles);
    local_scale.resize(n_angles);
    bool has_3D_wavelet = false;

    for (int j = 0; j < n_angles; j++) {

      local_noise_scale[j] = NULL;
      local_shift[j] = NULL;
      local_scale[j] = NULL;

      if (model_settings->getWaveletDim(j) == Wavelet::THREE_D)
        has_3D_wavelet = true;
      if (estimate_wavelets[j] == true)
        model_settings->setWaveletScale(i,j,1.0);
    }

    unsigned int n_wells = model_settings->getNumberOfWells();

    //Store ref_time_grad_x, ref_time_grad_y, t_grad_x and t_grad_y. They are needed to reestimate SNRatio of a 3DWavelet in ModelAVODynamic
    t_grad_x_.resize(n_wells);
    t_grad_y_.resize(n_wells);
    NRLib::Grid2D<float> structure_depth_grad_x; ///< Depth gradient in x-direction for structure ( correlationDirection-t0)*v0/2
    NRLib::Grid2D<float> structure_depth_grad_y; ///< Depth gradient in y-direction for structure ( correlationDirection-t0)*v0/2

    bool failed = false;

    if (has_3D_wavelet) {
      if (input_files->getRefSurfaceFile() != "") {
        Surface t0_surf;
        try {
          t0_surf =Surface(input_files->getRefSurfaceFile());
        }
        catch (NRLib::Exception & e) {
          err_text += e.what();
          failed = true;
        }
        if (!failed) {
          //Correlation direction from ModelGeneral::makeTimeSimboxes
          Surface * correlation_direction = NULL;

          if(input_files->getCorrDirFile() != "") { //H can be given as top and base?
            try {
              Surface tmp_surf(input_files->getCorrDirFile());
              if (estimation_simbox_.CheckSurface(tmp_surf) == true)
                correlation_direction = new Surface(tmp_surf);
              else {
                err_text += "Error: Correlation surface does not cover volume.\n";
              }
            }
            catch (NRLib::Exception & e) {
              err_text += e.what();
            }
          }

          double v0=model_settings->getAverageVelocity();
          ComputeStructureDepthGradient(v0,
                                        model_settings->getGradientSmoothingRange(),
                                        &t0_surf,
                                        correlation_direction,
                                        structure_depth_grad_x,
                                        structure_depth_grad_y);
          Wavelet3D::setGradientMaps(structure_depth_grad_x,
                                     structure_depth_grad_y);
          ComputeReferenceTimeGradient(&t0_surf,
                                       ref_time_grad_x_,
                                       ref_time_grad_y_);
        }
        else {
          err_text += "Problems reading reference time surface in (x,y).\n";
        }
      }
      bool estimate_well_gradient = model_settings->getEstimateWellGradientFromSeismic();
      float distance, sigma_m;
      model_settings->getTimeGradientSettings(distance, sigma_m, i);
      std::vector<std::vector<double> > SigmaXY;

      for (size_t w = 0; w < n_wells; w++) {
        if (!estimate_well_gradient & ((structure_depth_grad_x.GetN()> 0) & (structure_depth_grad_y.GetN()>0))) {
          double v0=model_settings->getAverageVelocity();
          mapped_blocked_logs_.find(wells_[w].GetWellName())->second->SetSeismicGradient(v0, structure_depth_grad_x, structure_depth_grad_y, ref_time_grad_x_, ref_time_grad_y_, t_grad_x_[w], t_grad_y_[w]);
        }
        else {
          mapped_blocked_logs_.find(wells_[w].GetWellName())->second->SetTimeGradientSettings(distance, sigma_m);
          mapped_blocked_logs_.find(wells_[w].GetWellName())->second->FindSeismicGradient(seismic_data_[i], &estimation_simbox_, n_angles, t_grad_x_[w], t_grad_y_[w], SigmaXY);
        }
      }
    }

    if (estimation_simbox_.getdz() > 4.01f && model_settings->getEstimateNumberOfWavelets(i) > 0) { // Require this density for wavelet estimation
      LogKit::LogFormatted(LogKit::Low,"\n\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
      LogKit::LogFormatted(LogKit::Low," CRAVA are not reliable and the output results should be treated accordingly.\n");
      LogKit::LogFormatted(LogKit::Low," The number of layers must be increased. \n");
      std::string text("");
      text += "Increase the number of layers to improve the quality of the wavelet estimation.\n";
      text += " The minimum sampling density is "+NRLib::ToString(estimation_simbox_.getdz())+", and it should be ";
      text += "lower than 4.0.\n To obtain the desired density, the number of layers should be at least ";
      text += NRLib::ToString(static_cast<int>(estimation_simbox_.GetLZ()/4.0))+"\n";
      TaskList::addTask(text);
    }

    // check if local noise is set for some angles.
    bool local_noise_set = false;
    std::vector<std::vector<double> > synt_seis_angles(n_angles);

    for (int j = 0; j < n_angles; j++) {

      float angle = float(angles[j]*180.0/M_PI);
      LogKit::LogFormatted(LogKit::Low,"\nAngle stack : %.1f deg",angle);

      SeismicStorage * seismic_data = NULL;
      if (forward_modeling_ == false)
        seismic_data = &seismic_data_[i][j];

      if (model_settings->getWaveletDim(j) == Wavelet::ONE_D)
        error += Process1DWavelet(model_settings,
                                  input_files,
                                  seismic_data,
                                  mapped_blocked_logs_,
                                  wavelet_estim_interval,
                                  *reflection_matrix_[i],
                                  synt_seis_angles[j],
                                  err_text,
                                  wavelet[j],
                                  local_noise_scale[j],
                                  local_shift[j],
                                  local_scale[j],
                                  i, //Timelapse
                                  j, //Angle
                                  angles[j],
                                  sn_ratio[j],
                                  estimate_wavelets[j],
                                  use_ricker_wavelet[j],
                                  use_local_noise_);
      else
        error += Process3DWavelet(model_settings,
                                  input_files,
                                  seismic_data,
                                  mapped_blocked_logs_,
                                  wavelet_estim_interval,
                                  *reflection_matrix_[i],
                                  err_text,
                                  wavelet[j],
                                  i, //Timelapse
                                  j, //Angle
                                  angles[j],
                                  sn_ratio[j],
                                  ref_time_grad_x_,
                                  ref_time_grad_y_,
                                  t_grad_x_,
                                  t_grad_y_,
                                  estimate_wavelets[j]);

      if (local_noise_scale[j] != NULL)
        local_noise_set = true;

    } //angle

    if (local_noise_set == true) {
      for (int k=0; k < n_angles; k++) {
        if (local_noise_scale[k]==NULL)
          local_noise_scale[k] = new Grid2D(estimation_simbox_.getnx(),
                                            estimation_simbox_.getny(),
                                            1.0);
      }
    }

    wavelets_[i]              = wavelet;
    local_noise_scale_[i]     = local_noise_scale;
    local_shift_[i]           = local_shift;
    local_scale_[i]           = local_scale;
    global_noise_estimate_[i] = sn_ratio;
    sn_ratio_[i]              = sn_ratio;
    synt_seis_[i]             = synt_seis_angles;
  } //timelapse



  Timings::setTimeWavelets(wall,cpu);

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }
  else if (error > 0) { //Should be covered by err_text
    err_text_common += err_text;
    return false;
  }

  return true;
}

void
CommonData::FindWaveletEstimationInterval(InputFiles             * input_files,
                                          std::vector<Surface *> & wavelet_estim_interval,
                                          std::string            & err_text)

{
  const double x0 = estimation_simbox_.getx0();
  const double y0 = estimation_simbox_.gety0();
  const double lx = estimation_simbox_.getlx();
  const double ly = estimation_simbox_.getly();
  const int nx    = estimation_simbox_.getnx();
  const int ny    = estimation_simbox_.getny();

  //
  // Get wavelet estimation interval
  //
  const std::string & topWEI  = input_files->getWaveletEstIntFileTop(0); //Same for all time lapses
  const std::string & baseWEI = input_files->getWaveletEstIntFileBase(0);//Same for all time lapses

  if (topWEI != "" && baseWEI != "") {
    wavelet_estim_interval.resize(2);
    try {
      if (NRLib::IsNumber(topWEI))
        wavelet_estim_interval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWEI.c_str()));
      else {
        Surface tmpSurf(topWEI);
        wavelet_estim_interval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      err_text += e.what();
    }

    try {
      if (NRLib::IsNumber(baseWEI))
        wavelet_estim_interval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWEI.c_str()));
      else {
        Surface tmpSurf(baseWEI);
        wavelet_estim_interval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      err_text += e.what();
    }
  }

}

int
CommonData::Process1DWavelet(const ModelSettings                      * model_settings,
                             const InputFiles                         * input_files,
                             const SeismicStorage                     * seismic_data,
                             std::map<std::string, BlockedLogsCommon *> mapped_blocked_logs,
                             const std::vector<Surface *>             & wavelet_estim_interval,
                             const float                              * reflection_matrix,
                             std::vector<double>                      & synt_seis,
                             std::string                              & err_text,
                             Wavelet                                 *& wavelet,
                             Grid2D                                  *& local_noise_scale, //local noise estimates?
                             Grid2D                                  *& local_shift,
                             Grid2D                                  *& local_scale,
                             unsigned int                               i_timelapse,
                             unsigned int                               j_angle,
                             const float                                angle,
                             float                                    & sn_ratio,
                             bool                                       estimate_wavelet,
                             bool                                       use_ricker_wavelet,
                             bool                                       use_local_noise)
{

  int error = 0;
  //Grid2D * shiftGrid(NULL);
  //Grid2D * gainGrid(NULL);
  if (model_settings->getUseLocalWavelet() && input_files->getScaleFile(i_timelapse,j_angle) != "") {
      Surface help(input_files->getScaleFile(i_timelapse, j_angle));
      //gainGrid = new Grid2D(estimation_simbox_.getnx(),estimation_simbox_.getny(), 0.0);
      local_scale = new Grid2D(estimation_simbox_.getnx(),estimation_simbox_.getny(), 0.0);
      ResampleSurfaceToGrid2D(&help, local_scale);
  }
  if (model_settings->getUseLocalWavelet() && input_files->getShiftFile(i_timelapse,j_angle) != "") {
    Surface helpShift(input_files->getShiftFile(i_timelapse, j_angle));
    local_shift = new Grid2D(estimation_simbox_.getnx(),estimation_simbox_.getny(), 0.0);
    ResampleSurfaceToGrid2D(&helpShift, local_shift);
  }
  if (use_local_noise && input_files->getLocalNoiseFile(i_timelapse,j_angle) != "") {
    Surface helpNoise(input_files->getLocalNoiseFile(i_timelapse, j_angle));
    local_noise_scale = new Grid2D(estimation_simbox_.getnx(), estimation_simbox_.getny(), 0.0);
    ResampleSurfaceToGrid2D(&helpNoise, local_noise_scale);
  }

  Wavelet * wavelet_pre_resampling = NULL;

  if (estimate_wavelet) {
    wavelet = new Wavelet1D(&estimation_simbox_,
                            seismic_data,
                            mapped_blocked_logs,
                            wavelet_estim_interval,
                            model_settings,
                            reflection_matrix,
                            synt_seis,
                            j_angle,
                            error,
                            err_text);
  }
  else { //Not estimation modus
    if (use_ricker_wavelet)
        wavelet = new Wavelet1D(model_settings,
                                reflection_matrix,
                                angle,
                                model_settings->getRickerPeakFrequency(i_timelapse,j_angle),
                                error);
    else {
      const std::string & wavelet_file = input_files->getWaveletFile(i_timelapse,j_angle);
      int file_format = GetWaveletFileFormat(wavelet_file, err_text);
      if (file_format < 0) {
        err_text += "Unknown file format of file '"+wavelet_file+"'.\n";
        error++;
      }
      else
        wavelet = new Wavelet1D(wavelet_file,
                                file_format,
                                model_settings,
                                reflection_matrix,
                                angle,
                                error,
                                err_text);
    }
    // Calculate a preliminary scale factor to see if wavelet is in the same size order as the data. A large or small value might cause problems.
    if (seismic_data != NULL) { // If forward modeling, we have no seismic, can not prescale wavelet.
      float       prescale  = wavelet->findGlobalScaleForGivenWavelet(model_settings, &estimation_simbox_, seismic_data, mapped_blocked_logs);
      const float lim_high  = 3.0f;
      const float lim_low   = 0.33f;

      if (model_settings->getEstimateGlobalWaveletScale(i_timelapse,j_angle)) // prescale, then we have correct size order, and later scale estimation will be ok.
          wavelet->multiplyRAmpByConstant(prescale);
      else {
        if (model_settings->getWaveletScale(i_timelapse,j_angle)!= 1.0f && (prescale>lim_high || prescale<lim_low)) {
            std::string text = "The wavelet given for angle no "+NRLib::ToString(j_angle)+" is badly scaled. Ask Crava to estimate global wavelet scale.\n";
          if (model_settings->getEstimateLocalScale(i_timelapse,j_angle)) {
            err_text += text;
            error++;
          }
          else {
            LogKit::LogFormatted(LogKit::Warning,"\nWARNING: "+text);
            TaskList::addTask("The wavelet is badly scaled. Consider having CRAVA estimate global wavelet scale");
          }
        }
      }
    }

    //If wavelet isn't estimated it needs to be resampled in order to calculate SN-Ratio
    //The wavelet also needs to be resampled in modelAVODynamic to the different intervals.
    //Here we make a different copy to be sent to modelAVODynamic, so the same wavelet isn't resampled twice.
    wavelet_pre_resampling = new Wavelet1D(wavelet);

    if (error == 0) {

      wavelet->resample(static_cast<float>(estimation_simbox_.getdz()),
                        estimation_simbox_.getnz(),
                        estimation_simbox_.getnz());
    }

  }

  if (error == 0) {
    wavelet->scale(model_settings->getWaveletScale(i_timelapse,j_angle));

    if (forward_modeling_ == false && model_settings->getNumberOfWells() > 0) {

      std::vector<std::vector<double> > seis_logs(mapped_blocked_logs.size());
      int w = 0;
        for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = mapped_blocked_logs.begin(); it != mapped_blocked_logs.end(); it++) {
        std::map<std::string, BlockedLogsCommon *>::const_iterator iter = mapped_blocked_logs.find(it->first);
        BlockedLogsCommon * blocked_log = iter->second;

        seis_logs[w].resize(blocked_log->GetNumberOfBlocks());
        blocked_log->GetBlockedGrid(seismic_data, &estimation_simbox_, seis_logs[w]);

        w++;
      }

      float SNRatio_tmp = wavelet->calculateSNRatioAndLocalWavelet(&estimation_simbox_,
                                                                   seis_logs,
                                                                   mapped_blocked_logs,
                                                                   model_settings,
                                                                   err_text,
                                                                   error,
                                                                   j_angle,
                                                                   local_noise_scale,
                                                                   //shiftGrid,
                                                                   local_shift,
                                                                   //gainGrid,
                                                                   local_scale,
                                                                   sn_ratio,
                                                                   model_settings->getWaveletScale(i_timelapse, j_angle),
                                                                   model_settings->getEstimateSNRatio(i_timelapse, j_angle),
                                                                   model_settings->getEstimateGlobalWaveletScale(i_timelapse, j_angle),
                                                                   model_settings->getEstimateLocalNoise(i_timelapse, j_angle),
                                                                   model_settings->getEstimateLocalShift(i_timelapse, j_angle),
                                                                   model_settings->getEstimateLocalScale(i_timelapse, j_angle),
                                                                   estimate_wavelet);
      if (model_settings->getEstimateSNRatio(i_timelapse,j_angle))
        sn_ratio = SNRatio_tmp;
    }

    if (error == 0) {

      if ((model_settings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (model_settings->getEstimationMode() && estimate_wavelet)) {
        std::string type;
        if (estimate_wavelet) {
          type = "Estimated_";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,true); // dt_max = 1.0;
        }
        else if (model_settings->getWaveletScale(i_timelapse,j_angle) == 1.00) {
          type = "";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
        else {
          type = "Scaled_";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
      }
      const float SNLow  = 1.0;
      const float SNHigh = 10.0;
      if ((sn_ratio <=SNLow  || sn_ratio > SNHigh) && forward_modeling_==false) {
        err_text += "Illegal signal-to-noise ratio of "+NRLib::ToString(sn_ratio)+" for cube "+NRLib::ToString(j_angle+1)+".\n";
        err_text += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
        error++;
      }

      bool useLocalNoise = model_settings->getEstimateLocalNoise(i_timelapse,j_angle) || input_files->getLocalNoiseFile(i_timelapse,j_angle) != "";
      bool useLocalShift = model_settings->getEstimateLocalShift(i_timelapse,j_angle) || input_files->getShiftFile(i_timelapse,j_angle)      != "";
      bool useLocalGain  = model_settings->getEstimateLocalScale(i_timelapse,j_angle) || input_files->getScaleFile(i_timelapse,j_angle)      != "";

      if (useLocalNoise)
        ReadAndWriteLocalGridsToFile(input_files->getLocalNoiseFile(i_timelapse,j_angle),
                                     IO::PrefixLocalNoise(),
                                     1.0,  // Scale map with this factor before writing to disk
                                     model_settings,
                                     local_noise_scale,
                                     angle);

      if (useLocalShift) {
        ReadAndWriteLocalGridsToFile(input_files->getShiftFile(i_timelapse,j_angle),
                                     IO::PrefixLocalWaveletShift(),
                                     1.0,
                                     model_settings,
                                     local_shift,
                                     angle);
        wavelet->setShiftGrid(local_shift);
      }

      if (useLocalGain) {
        ReadAndWriteLocalGridsToFile(input_files->getScaleFile(i_timelapse,j_angle),
                                     IO::PrefixLocalWaveletGain(),
                                     1.0,
                                     model_settings,
                                     local_scale,
                                     angle);
        wavelet->setGainGrid(local_scale);
      }
    }
  }

  if (wavelet_pre_resampling != NULL)
    wavelet = new Wavelet1D(wavelet_pre_resampling);

  return error;
}

int
CommonData::Process3DWavelet(const ModelSettings                      * model_settings,
                             const InputFiles                         * input_files,
                             const SeismicStorage                     * seismic_data,
                             std::map<std::string, BlockedLogsCommon *> mapped_blocked_logs,
                             const std::vector<Surface *>             & wavelet_estim_interval,
                             const float                              * reflection_matrix,
                             std::string                              & err_text,
                             Wavelet                                 *& wavelet,
                             unsigned int                               i_timelapse,
                             unsigned int                               j_angle,
                             float                                      angle,
                             float                                      sn_ratio,
                             const NRLib::Grid2D<float>               & ref_time_grad_x,
                             const NRLib::Grid2D<float>               & ref_time_grad_y,
                             const std::vector<std::vector<double> >  & t_grad_x,
                             const std::vector<std::vector<double> >  & t_grad_y,
                             bool                                       estimate_wavelet)
{
  int error = 0;
  if (estimate_wavelet) {
    wavelet = new Wavelet3D(input_files->getWaveletFilterFile(j_angle),
                            wavelet_estim_interval,
                            ref_time_grad_x,
                            ref_time_grad_y,
                            t_grad_x,
                            t_grad_y,
                            seismic_data,
                            model_settings,
                            mapped_blocked_logs,
                            &estimation_simbox_,
                            reflection_matrix,
                            j_angle,
                            error,
                            err_text);
  }
  else { //Not estimation modus
    const std::string & wavelet_file = input_files->getWaveletFile(i_timelapse,j_angle);
    int file_format = GetWaveletFileFormat(wavelet_file, err_text);
    if (file_format < 0) {
      err_text += "Unknown file format of file '"+wavelet_file+"'.\n";
      error++;
    }
    else {
      wavelet = new Wavelet3D(wavelet_file,
                              file_format,
                              model_settings,
                              reflection_matrix,
                              angle,
                              error,
                              err_text,
                              input_files->getWaveletFilterFile(j_angle));
      if (error == 0)
        wavelet->resample(static_cast<float>(estimation_simbox_.getdz()),
                          estimation_simbox_.getnz(),
                          estimation_simbox_.GetNZpad());
    }
  }
  if ((model_settings->getEstimationMode() == false) && estimation_simbox_.getIsConstantThick()) {
    err_text += "Simbox with constant thicknessis not implemented for modelling or inversion when 3D wavelet.\n";
    error++;
  }
  if (error == 0) {
    wavelet->scale(model_settings->getWaveletScale(i_timelapse,j_angle));
    bool localEst = (model_settings->getEstimateLocalScale(i_timelapse,j_angle) ||
                     model_settings->getEstimateLocalShift(i_timelapse,j_angle) ||
                     model_settings->getEstimateLocalNoise(i_timelapse,j_angle) ||
                     model_settings->getEstimateGlobalWaveletScale(i_timelapse,j_angle) ||
                     model_settings->getEstimateSNRatio(i_timelapse,j_angle));

    if (localEst && forward_modeling_ == false) {
      float sn_ratio_tmp = wavelet->calculateSNRatio(&estimation_simbox_,
                                                     seismic_data,
                                                     mapped_blocked_logs,
                                                     model_settings,
                                                     err_text,
                                                     error,
                                                     ref_time_grad_x,
                                                     ref_time_grad_y,
                                                     t_grad_x,
                                                     t_grad_y,
                                                     j_angle,
                                                     sn_ratio,
                                                     model_settings->getEstimateSNRatio(i_timelapse,j_angle),
                                                     estimate_wavelet);
      if (model_settings->getEstimateSNRatio(i_timelapse,j_angle))
        sn_ratio = sn_ratio_tmp;

    }
    if (error == 0) {
      if ((model_settings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (model_settings->getEstimationMode() && estimate_wavelet)) {
        std::string type;
        if (estimate_wavelet) {
          type = "Estimated_";
          if (wavelet->getDim()==1)
            wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,true); // dt_max = 1.0;
          else
            wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
        else if (model_settings->getWaveletScale(i_timelapse,j_angle) == 1.00) {
          type = "";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
        else {
          type = "Scaled_";
          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
        }
      }

      const float SNLow = 1.0;
      const float SNHigh = 10.0;
      if ((sn_ratio <=SNLow || sn_ratio > SNHigh) && forward_modeling_==false) {
        err_text += "Illegal signal-to-noise ratio of "+NRLib::ToString(sn_ratio)+" for cube "+NRLib::ToString(j_angle+1)+".\n";
        err_text += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
        error++;
      }
    }
  }

  return error;
}

void
CommonData::ComputeStructureDepthGradient(double                 v0,
                                          double                 radius,
                                          const Surface        * t0_surf,
                                          const Surface        * correlation_direction,
                                          NRLib::Grid2D<float> & structure_depth_grad_x,
                                          NRLib::Grid2D<float> & structure_depth_grad_y)
 {
   double ds = 12.5;

   int nx = estimation_simbox_.getnx();
   int ny = estimation_simbox_.getny();
   structure_depth_grad_x.Resize(nx,ny);
   structure_depth_grad_y.Resize(nx,ny);
   double mp=v0*0.001*0.5; // 0.001 is due to s vs ms convension

   for (int i = 0; i < nx; i++) {
     for (int j = 0; j < ny; j++) {
       double x,y;
       double gx,gy,gxTmp,gyTmp;
       gx=0.0;
       gy=0.0;
       estimation_simbox_.getXYCoord(i,j,x,y);
       CalculateSmoothGrad(t0_surf, x, y, radius, ds,gxTmp, gyTmp);
       gx=-gxTmp;
       gy=-gyTmp;
       if (correlation_direction != NULL) {
         CalculateSmoothGrad(correlation_direction, x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp;
         gy+=gyTmp;
       }
       else {
         CalculateSmoothGrad( &(dynamic_cast<const Surface &> (estimation_simbox_.GetTopSurface())), x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp*0.5;
         gy+=gyTmp*0.5;
         CalculateSmoothGrad( &(dynamic_cast<const Surface &> (estimation_simbox_.GetBotSurface())), x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp*0.5;
         gy+=gyTmp*0.5;
       }

       gx*=mp;
       gy*=mp;
       structure_depth_grad_x(i,j) =float(gx);
       structure_depth_grad_y(i,j) =float(gy);
     }
   }
 }


void
CommonData::CalculateSmoothGrad(const Surface * surf, double x, double y, double radius, double ds, double& gx, double& gy)
{
  /// Return smoothed Gradient. Computes the gradient as a regression
  /// among points within a given distance from the central point.
  // Returns missing if central point is outside the grid.
  /// Returns otherwise it is ok,
  int i,j,k,l;
  int disc_radius = int(floor(radius/ds));
  int n_points = (2*disc_radius+1)*(2*disc_radius+1);
  std::vector<double> Z(3*n_points,0.0);
  std::vector<double> Y(n_points,0.0);
  std::vector<double> cov(9,0.0);
  std::vector<double> invCov(9,0.0);
  std::vector<double> proj(3,0.0);
  double z0 = 0.0;
  int    cy = 0;
  int    cz = 0;

  double baseDepth =surf->GetZ(x,y);

  for (i=- disc_radius; i <= disc_radius; i++) {
    for (j=- disc_radius; j <= disc_radius; j++) {
      double dx= i*ds;
      double dy= j*ds;
      double foo=surf->GetZ(x+dx,y+dy);
      Y[cy] = foo;
      Z[cz] = 1.0;
      Z[cz + 1] = dx;
      Z[cz + 2] = dy;
      cy++;
      cz += 3;
    }
  }

  int nData=0;
  for (i = 0; i < n_points; i++) {
    if (!surf->IsMissing(Y[i])) {
      nData++;
      for (k=0; k < 3; k++) {
        for (l=0;l<3;l++)
          cov[k+l*3]+=Z[k + 3*i] * Z[l + 3*i];

        proj[k]+=Z[k + 3*i]*(Y[i]-baseDepth);
      }
    }
  }

  double det = cov[0]*(cov[4]*cov[8] - cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8] - cov[5]*cov[6])
                  + cov[2]*(cov[3]*cov[7] - cov[4]*cov[6]);

  if (det != 0) {
      invCov[0] = (cov[4]*cov[8] - cov[5]*cov[7]) / det;
      invCov[1] = (cov[2]*cov[7] - cov[1]*cov[8]) / det;
      invCov[2] = (cov[1]*cov[5] - cov[2]*cov[4]) / det;
      invCov[3] = (cov[5]*cov[6] - cov[3]*cov[8]) / det;
      invCov[4] = (cov[0]*cov[8] - cov[2]*cov[6]) / det;
      invCov[5] = (cov[2]*cov[3] - cov[0]*cov[5]) / det;
      invCov[6] = (cov[3]*cov[7] - cov[4]*cov[6]) / det;
      invCov[7] = (cov[1]*cov[6] - cov[0]*cov[7]) / det;
      invCov[8] = (cov[0]*cov[4] - cov[1]*cov[3]) / det;

      z0 = baseDepth;
      gx = 0.0;
      gy = 0.0;
      for (k=0; k < 3; k++) {
        z0 += invCov[k]*proj[k]; //NBNB check
        gx += invCov[3+k]*proj[k];
        gy += invCov[6+k]*proj[k];
      }
  }
  else {
   gx = RMISSING;
   gy = RMISSING;
  }
}

void
CommonData::ComputeReferenceTimeGradient(const Surface       * t0_surf,
                                         NRLib::Grid2D<float> &ref_time_grad_x,
                                         NRLib::Grid2D<float> &ref_time_grad_y)
 {
   double radius = 50.0;
   double ds     = 12.5;
   int nx        = estimation_simbox_.getnx();
   int ny        = estimation_simbox_.getny();
   ref_time_grad_x.Resize(nx,ny);
   ref_time_grad_y.Resize(nx,ny);
   for (int i = 0; i < nx; i++) {
     for (int j = 0; j < ny; j++) {
       double x,y;
       double gx,gy;
       estimation_simbox_.getXYCoord(i,j,x,y);
       CalculateSmoothGrad(t0_surf, x, y, radius, ds,gx, gy);
       ref_time_grad_x(i,j) =float(gx);
       ref_time_grad_y(i,j) =float(gy);
     }
   }
 }

void
CommonData::ResampleSurfaceToGrid2D(const Surface * surface,
                                    Grid2D        * outgrid)
{
  for (int i=0; i < estimation_simbox_.getnx(); i++) {
    for (int j=0; j < estimation_simbox_.getny(); j++) {
      double x, y, z;
      estimation_simbox_.getCoord(i, j, 0, x, y, z);
      (*outgrid)(i,j) = static_cast<float>(surface->GetZ(x,y));
    }
  }
}

int
CommonData::GetWaveletFileFormat(const std::string & file_name, std::string & err_text)
{
  int fileformat = -1;
  int line       = 0;
  int pos;
  std::string dummyStr;
  std::string targetString;

  std::ifstream file;
  NRLib::OpenRead(file, file_name);

  std::getline(file,dummyStr);
  line++;
  targetString = "pulse file-3";
  pos = Utils::findEnd(dummyStr, 0, targetString);
  if (pos >= 0)
    fileformat = Wavelet::NORSAR;
  file.close();
  file.clear();

  if (fileformat<0) { // not norsar format
      // test for jason file format
    NRLib::OpenRead(file,file_name);
    line         = 0;
    int thisLine = 0;
    bool lineIsComment = true;
    while (lineIsComment == true) {
      NRLib::ReadNextToken(file,dummyStr,line);
      if (NRLib::CheckEndOfFile(file)) {
        err_text += "End of wavelet file "+file_name+" is premature\n";
        return 0;
      }
      else {
        if (thisLine == line) {
          NRLib::DiscardRestOfLine(file,line,false);
          thisLine = line;
        }
        if ((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
          lineIsComment = false;
      }
    }
    file.close();
    if (NRLib::IsNumber(dummyStr)) // not convertable number
      fileformat= Wavelet::JASON;
  }
  return fileformat;
}

void
CommonData::ReadAndWriteLocalGridsToFile(const std::string   & file_name,
                                         const std::string   & type,
                                         const float           scale_factor,
                                         const ModelSettings * model_settings,
                                         const Grid2D        * grid,
                                         const float           angle)
{
  bool   estimationMode   = model_settings->getEstimationMode();
  int    outputFormat     = model_settings->getOutputGridFormat();
  double angle_tmp        = angle*180.0/M_PI;

  Surface * help = NULL;

  if (file_name != "") {
    std::string toPath = NRLib::RemovePath(file_name);

    if (type == IO::PrefixLocalNoise())
      toPath = NRLib::PrependDir(IO::PathToNoise(), toPath);
    else
      toPath = NRLib::PrependDir(IO::PathToWavelets(), toPath);

    NRLib::CreateDirIfNotExists(toPath);
    NRLib::CopyFile(file_name, toPath, true);
  }
  else {
    if (grid != NULL) {
      ResampleGrid2DToSurface(&estimation_simbox_, grid, help);
    }
  }
  if ((estimationMode ||
    ((type==IO::PrefixLocalWaveletGain() || type==IO::PrefixLocalWaveletShift()) && (model_settings->getWaveletOutputFlag() & IO::LOCAL_WAVELETS)>0) ||
    (type==IO::PrefixLocalNoise() && (model_settings->getOtherOutputFlag() & IO::LOCAL_NOISE)>0)) &&
     help != NULL)
  {
    std::string baseName = type + NRLib::ToString(angle_tmp, 1);
    help->Multiply(scale_factor);
    if (type==IO::PrefixLocalNoise())
      IO::writeSurfaceToFile(*help, baseName, IO::PathToNoise(), outputFormat);
    else
      IO::writeSurfaceToFile(*help, baseName, IO::PathToWavelets(), outputFormat);
  }
  if (help != NULL)
    delete help;
}

void
CommonData::ResampleGrid2DToSurface(const Simbox   * simbox,
                                    const Grid2D   * grid,
                                    Surface       *& surface)
{
  double xmin,xmax,ymin,ymax;
  simbox->getMinAndMaxXY(xmin,xmax,ymin,ymax);
  int nx,ny;
  double angle = simbox->getAngle()*180.0/M_PI;
  if (angle > -45 || angle < 45)
  {
    nx = static_cast<int>(floor(simbox->getnx()*1.0/std::cos(simbox->getAngle())+0.5)) * 2;
    ny = static_cast<int>(floor(simbox->getny()*1.0/std::cos(simbox->getAngle())+0.5)) * 2;
  }
  else
  {
    nx = static_cast<int>(floor(simbox->getnx()*1.0/std::sin(simbox->getAngle())+0.5)) * 2;
    ny = static_cast<int>(floor(simbox->getny()*1.0/std::sin(simbox->getAngle())+0.5)) * 2;
  }
  surface = new Surface(xmin,ymin,xmax-xmin,ymax-ymin,nx,ny,0.0);
  double x,y;
  int i1,j1;
  for (int i=0;i<nx;i++) {
    for (int j=0;j<ny;j++) {
      surface->GetXY(i,j,x,y);
      simbox->getIndexes(x,y,i1,j1);
      if (i1==IMISSING || j1== IMISSING)
        surface->SetMissing(i,j);
      else
        (*surface)(i,j) = (*grid)(i1,j1);
    }
  }
}

//bool CommonData::optimizeWellLocations() {
//  return true;
//}

//bool CommonData::SetupEstimationRockPhysics(){
//  return true;
//}

int
CommonData::ComputeTime(int year, int month, int day) const
{
  if (year == IMISSING)
    return(0);

  int delta_year = year-1900; //Ok baseyear.
  int time = 365*delta_year+delta_year/4; //Leap years.
  if (month == IMISSING)
    time += 182;
  else {
    std::vector<int> acc_days(12,0);
    acc_days[1]  = acc_days[0]  + 31;
    acc_days[2]  = acc_days[1]  + 28;
    acc_days[3]  = acc_days[2]  + 31;
    acc_days[4]  = acc_days[3]  + 30;
    acc_days[5]  = acc_days[4]  + 31;
    acc_days[6]  = acc_days[5]  + 30;
    acc_days[7]  = acc_days[6]  + 31;
    acc_days[8]  = acc_days[7]  + 31;
    acc_days[9]  = acc_days[8]  + 30;
    acc_days[10] = acc_days[9]  + 31;
    acc_days[11] = acc_days[10] + 30;

    time += acc_days[month];

    if (day == IMISSING)
      time += 15;
    else
      time += day;
  }
  return(time);
}

void
CommonData::GetGeometryFromGridOnFile(const std::string           grid_file,
                                      const TraceHeaderFormat   * thf,
                                      SegyGeometry             *& geometry,
                                      std::string               & err_text)
{
  geometry = NULL;

  if (grid_file != "") { //May change the condition here, but need geometry if we want to set XL/IL
    int file_type = IO::findGridType(grid_file);
    if (file_type == IO::CRAVA) {
      geometry = GetGeometryFromCravaFile(grid_file);
    }
    else if (file_type == IO::SEGY) {
      try
      {
        geometry = SegY::FindGridGeometry(grid_file, thf);
      }
      catch (NRLib::Exception & e)
      {
        err_text = e.what();
      }
    }
    else if (file_type == IO::STORM)
      geometry = GetGeometryFromStormFile(grid_file, err_text);
    else if (file_type==IO::SGRI) {
      bool scale = true;
      geometry = GetGeometryFromStormFile(grid_file, err_text, scale);
    }
    else {
      err_text = "Trying to read grid dimensions from unknown file format.\n";
    }
  }
  else {
    err_text = "Cannot get geometry from file. The file name is empty.\n";
  }
}

SegyGeometry * CommonData::GetGeometryFromCravaFile(const std::string & file_name)
{
  std::ifstream bin_file;
  NRLib::OpenRead(bin_file, file_name, std::ios::in | std::ios::binary);

  std::string file_type;
  getline(bin_file,file_type);

  double x0        = NRLib::ReadBinaryDouble(bin_file);
  double y0        = NRLib::ReadBinaryDouble(bin_file);
  double dx        = NRLib::ReadBinaryDouble(bin_file);
  double dy        = NRLib::ReadBinaryDouble(bin_file);
  int    nx        = NRLib::ReadBinaryInt(bin_file);
  int    ny        = NRLib::ReadBinaryInt(bin_file);
  double IL0       = NRLib::ReadBinaryDouble(bin_file);
  double XL0       = NRLib::ReadBinaryDouble(bin_file);
  double il_step_X = NRLib::ReadBinaryDouble(bin_file);
  double il_step_Y = NRLib::ReadBinaryDouble(bin_file);
  double xl_step_X = NRLib::ReadBinaryDouble(bin_file);
  double xl_step_Y = NRLib::ReadBinaryDouble(bin_file);
  double rot       = NRLib::ReadBinaryDouble(bin_file);

  bin_file.close();

  SegyGeometry * geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                             IL0, XL0, il_step_X, il_step_Y,
                                             xl_step_X, xl_step_Y, rot);
  return(geometry);
}

SegyGeometry * CommonData::GetGeometryFromStormFile(const std::string & file_name,
                                                    std::string       & err_text,
                                                    bool                scale)
{
  SegyGeometry  * geometry  = NULL;
  StormContGrid * storm_grid = NULL;
  std::string     tmp_err_text;
  float scale_hor;
  if (scale==false)
  {
    scale_hor = 1.0;
  }
  else //from sgri file
  {
    LogKit::LogFormatted(LogKit::Low,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
    scale_hor  = 1000.0;
  }
  try
  {
    storm_grid = new StormContGrid(0,0,0);
    storm_grid->ReadFromFile(file_name);
    storm_grid->SetMissingCode(RMISSING);
  }
  catch (NRLib::Exception & e)
  {
    tmp_err_text = e.what();
  }

  if (tmp_err_text == "") {
    double x0        = storm_grid->GetXMin()*scale_hor;
    double y0        = storm_grid->GetYMin()*scale_hor;
    double dx        = storm_grid->GetDX()*scale_hor;
    double dy        = storm_grid->GetDY()*scale_hor;
    int    nx        = static_cast<int>(storm_grid->GetNI());
    int    ny        = static_cast<int>(storm_grid->GetNJ());
    double rot       = storm_grid->GetAngle();
    double IL0       = 0.0;  ///< Dummy value since information is not contained in format
    double XL0       = 0.0;  ///< Dummy value since information is not contained in format
    double IL_step_X =   1;  ///< Dummy value since information is not contained in format
    double IL_step_Y =   1;  ///< Dummy value since information is not contained in format
    double XL_step_X =   1;  ///< Dummy value since information is not contained in format
    double XL_step_Y =   1;  ///< Dummy value since information is not contained in format
    geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, ///< When XL, IL is available.
                                IL0, XL0, IL_step_X, IL_step_Y,
                                XL_step_X, XL_step_Y, rot);
  }
  else {
    err_text += tmp_err_text;
  }

  if (storm_grid != NULL)
    delete storm_grid;

  return(geometry);
}

void CommonData::WriteAreas(const SegyGeometry * area_params,
                            Simbox             * time_simbox,
                            std::string        & text)
{
  double area_x0   = area_params->GetX0();
  double area_y0   = area_params->GetY0();
  double area_lx   = area_params->Getlx();
  double area_ly   = area_params->Getly();
  double area_dx   = area_params->GetDx();
  double area_dy   = area_params->GetDy();
  double area_rot  = area_params->GetAngle();
  double area_xmin = RMISSING;
  double area_xmax = RMISSING;
  double area_ymin = RMISSING;
  double area_ymax = RMISSING;

  FindSmallestSurfaceGeometry(area_x0, area_y0, area_lx, area_ly, area_rot,
                              area_xmin, area_ymin, area_xmax, area_ymax);

  LogKit::LogFormatted(LogKit::Low,"\nThe top and/or base time surfaces do not cover the area specified by the "+text);
  LogKit::LogFormatted(LogKit::Low,"\nPlease extrapolate surfaces or specify a smaller AREA in the model file.\n");
  LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx        ly     azimuth          dx      dy\n");
  LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");
  double azimuth = (-1)*area_rot*(180.0/NRLib::Pi);
  if (azimuth < 0)
    azimuth += 360.0;
  LogKit::LogFormatted(LogKit::Low,"Model area       %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n\n",
                       area_x0, area_y0, area_lx, area_ly, azimuth, area_dx, area_dy);

  LogKit::LogFormatted(LogKit::Low,"Area                    xmin         xmax           ymin        ymax\n");
  LogKit::LogFormatted(LogKit::Low,"--------------------------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"%-12s     %11.2f  %11.2f    %11.2f %11.2f\n",
                       text.c_str(),area_xmin, area_xmax, area_ymin, area_ymax);
  const NRLib::Surface<double> & top  = time_simbox->GetTopSurface();
  const NRLib::Surface<double> & base = time_simbox->GetBotSurface();
  LogKit::LogFormatted(LogKit::Low,"Top surface      %11.2f  %11.2f    %11.2f %11.2f\n",
                       top.GetXMin(), top.GetXMax(), top.GetYMin(), top.GetYMax());
  LogKit::LogFormatted(LogKit::Low,"Base surface     %11.2f  %11.2f    %11.2f %11.2f\n",
                       base.GetXMin(), base.GetXMax(), base.GetYMin(), base.GetYMax());
}


void CommonData::FindSmallestSurfaceGeometry(const double   x0,
                                             const double   y0,
                                             const double   lx,
                                             const double   ly,
                                             const double   rot,
                                             double       & x_min,
                                             double       & y_min,
                                             double       & x_max,
                                             double       & y_max)
{
  x_min = x0 - ly*sin(rot);
  x_max = x0 + lx*cos(rot);
  y_min = y0;
  y_max = y0 + lx*sin(rot) + ly*cos(rot);
  if (rot < 0) {
    x_min = x0;
    x_max = x0 + lx*cos(rot) - ly*sin(rot);
    y_min = y0 + lx*sin(rot);
    y_max = y0 + ly*cos(rot);
  }
}

int CommonData::GetNzFromGridOnFile(ModelSettings     * model_settings,
                                    const std::string & grid_file,
                                    std::string       & err_text) {

  int nz = 0;
  if (grid_file != "") {
    int file_type = IO::findGridType(grid_file);

    if (file_type == IO::CRAVA) {
      GetZPaddingFromCravaFile(grid_file, err_text, nz);
    }

    else if (file_type == IO::SEGY) {

      SegY              * segy   = NULL;
      TraceHeaderFormat * format = model_settings->getTraceHeaderFormat(0, 0);
      float offset               = model_settings->getSegyOffset(0);

      if (format == NULL) { //Unknown format
        std::vector<TraceHeaderFormat*> traceHeaderFormats(0);

        if (model_settings->getTraceHeaderFormat() != NULL)
          traceHeaderFormats.push_back(model_settings->getTraceHeaderFormat());

        segy = new SegY(grid_file,
                        offset,
                        traceHeaderFormats,
                        true); // Add standard formats to format search
      }
      else { //Known format, read directly.
        segy = new SegY(grid_file, offset, *format);
      }

      nz = segy->GetNz();
      if (segy != NULL)
        delete segy;
    }

    else if (file_type == IO::STORM || file_type == IO::SGRI) {

      StormContGrid * stormgrid = NULL;
      stormgrid = new StormContGrid(0,0,0);
      stormgrid->ReadFromFile(grid_file);
      nz = stormgrid->GetNK();

      if (stormgrid != NULL)
        delete stormgrid;

    }
    else {
      err_text = "Trying to read grid dimensions from unknown file format.\n";
    }
  }
  else {
    err_text = "Cannot read grid from file. The file name is empty.\n";
  }

  return nz;
}

void CommonData::SetSurfacesSingleInterval(const ModelSettings              * const model_settings,
                                           NRLib::Volume                    & full_inversion_volume,
                                           Simbox                           & estimation_simbox,
                                           int                                nz,
                                           const std::vector<std::string>   & surf_file,
                                           std::string                      & err_text) {
  bool failed = false;

  const std::string       & top_surface_file_name         = surf_file[0];
  int                       nx                            = model_settings->getAreaParameters()->GetNx();
  int                       ny                            = model_settings->getAreaParameters()->GetNy();
  bool                      generate_seismic              = forward_modeling_;
  bool                      estimation_mode               = model_settings->getEstimationMode();
  bool                      generate_background           = model_settings->getGenerateBackground();
  bool                      parallel_surfaces             = model_settings->getParallelTimeSurfaces();
  //int                       nz                            = model_settings->getTimeNz();
  int                       output_format                 = model_settings->getOutputGridFormat();
  int                       output_domain                 = model_settings->getOutputGridDomain();
  int                       output_grids_elastic          = model_settings->getOutputGridsElastic();
  int                       output_grids_other            = model_settings->getOutputGridsOther();
  int                       output_grids_seismic          = model_settings->getOutputGridsSeismic();
  double                    d_top                         = model_settings->getTimeDTop();
  double                    lz                            = model_settings->getTimeLz();
  double                    dz                            = model_settings->getTimeDz();

  Surface * top_surface  = NULL;
  Surface * base_surface = NULL;

  Surface * top_surface_flat  = NULL;
  Surface * base_surface_flat = NULL;

  try {
    double x_min, x_max;
    double y_min, y_max;
    FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                estimation_simbox.getlx(), estimation_simbox.getly(),
                                estimation_simbox.getAngle(),
                                x_min,y_min,x_max,y_max);
    if (NRLib::IsNumber(top_surface_file_name)) {
      LogKit::LogFormatted(LogKit::Low,"Top surface: Flat surface at depth %11.2f \n", atof(top_surface_file_name.c_str()));
      // Find the smallest surface that covers the simbox.

      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, nx, ny, atof(top_surface_file_name.c_str()));
    }
    else {
      LogKit::LogFormatted(LogKit::Low,"Top surface file name: " + top_surface_file_name +" \n");
      top_surface = new Surface(top_surface_file_name);
    }
  }
  catch (NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  if (!failed) {
    if (parallel_surfaces) { //Only one reference surface
      //estimation_simbox.setDepth(*top_surface_flat, d_top, lz, dz, model_settings->getRunFromPanel()); //H Removed. Top_surface_flat not created yet.
      top_surface->Add(d_top);
      base_surface = new Surface(*top_surface);
      base_surface->Add(lz);
      LogKit::LogFormatted(LogKit::Low,"Base surface: parallel to the top surface, shifted %11.2f down.\n", lz);
      //full_inversion_volume.SetSurfaces(*top_surface, *base_surface, model_settings->getRunFromPanel()); //H
    }
    else { //Two reference surfaces
      const std::string & base_surface_file_name = surf_file[1];
      try {
        if (NRLib::IsNumber(base_surface_file_name)) {
          LogKit::LogFormatted(LogKit::Low,"Base surface: Flat surface at depth %11.2f \n", atof(base_surface_file_name.c_str()));
          // Find the smallest surface that covers the simbox.
          double x_min, x_max;
          double y_min, y_max;
          FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                      estimation_simbox.getlx(), estimation_simbox.getly(),
                                      estimation_simbox.getAngle(),
                                      x_min, y_min, x_max, y_max);
          base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, nx, ny, atof(base_surface_file_name.c_str()));
        }
        else {
          LogKit::LogFormatted(LogKit::Low,"Base surface file name: " + base_surface_file_name +" \n");
          base_surface = new Surface(base_surface_file_name);
        }
      }
      catch (NRLib::Exception & e) {
        err_text += e.what();
        failed = true;
      }

    }
    //if (!failed) {
    //  if ((output_domain & IO::TIMEDOMAIN) > 0) {
    //    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime();
    //    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
    //    estimation_simbox.setTopBotName(topSurf,baseSurf,output_format);
    //    if (generate_seismic) {
    //      estimation_simbox.writeTopBotGrids(topSurf,
    //                                         baseSurf,
    //                                         IO::PathToSeismicData(),
    //                                         output_format);
    //    }
    //    else if (!estimation_mode){
    //      if (output_grids_elastic > 0 || output_grids_other > 0 || output_grids_seismic > 0)
    //        estimation_simbox.writeTopBotGrids(topSurf,
    //                                           baseSurf,
    //                                           IO::PathToInversionResults(),
    //                                           output_format);
    //    }
    //    if ((output_format & IO::STORM) > 0) { // These copies are only needed with the STORM format
    //      if ((output_grids_elastic & IO::BACKGROUND) > 0 ||
    //          (output_grids_elastic & IO::BACKGROUND_TREND) > 0 ||
    //          (estimation_mode && generate_background)) {
    //        estimation_simbox.writeTopBotGrids(topSurf,
    //                                           baseSurf,
    //                                           IO::PathToBackground(),
    //                                           output_format);
    //      }
    //      if ((output_grids_other & IO::CORRELATION) > 0) {
    //        estimation_simbox.writeTopBotGrids(topSurf,
    //                                           baseSurf,
    //                                           IO::PathToCorrelations(),
    //                                           output_format);
    //      }
    //      if ((output_grids_seismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
    //        estimation_simbox.writeTopBotGrids(topSurf,
    //                                           baseSurf,
    //                                           IO::PathToSeismicData(),
    //                                           output_format);
    //      }
    //      if ((output_grids_other & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
    //        estimation_simbox.writeTopBotGrids(topSurf,
    //                                           baseSurf,
    //                                           IO::PathToVelocity(),
    //                                           output_format);
    //      }
    //    }
    //  }
    //}
  }
  if (!failed){
    try{ //initialize full_inversion_volume and set the flat top and base surfaces of the simbox

      full_inversion_volume = NRLib::Volume(estimation_simbox.getx0(), estimation_simbox.gety0(),
        estimation_simbox.getlx(), estimation_simbox.getly(),
        *top_surface, *base_surface, estimation_simbox.getAngle());

      // find the top and bottom points on the full inversion volume
      double z_min_top_surface = full_inversion_volume.GetTopZMin(nx,ny);
      double z_max_base_surface = full_inversion_volume.GetBotZMax(nx,ny);
      // round outwards from the center of the volume to nearest multiple 4
      double z_min = floor(z_min_top_surface/4)*4;
      double z_max = ceil(z_max_base_surface/4)*4;
      top_surface_flat = new Surface(top_surface->GetXMin(), top_surface->GetYMin(),
        top_surface->GetXMax()-top_surface->GetXMin(), top_surface->GetYMax()-top_surface->GetYMin(),
        nx, ny, z_min);
      base_surface_flat = new Surface(base_surface->GetXMin(), base_surface->GetYMin(),
        base_surface->GetXMax()-base_surface->GetXMin(), base_surface->GetYMax()-base_surface->GetYMin(),
        nx, ny, z_max);

      estimation_simbox.setDepth(*top_surface_flat, *base_surface_flat, nz, model_settings->getRunFromPanel());
      estimation_simbox.calculateDz(model_settings->getLzLimit(), err_text);
    }
    catch(NRLib::Exception & e) {
      err_text += e.what();
    }
  }
  if (!failed) {
    if ((output_domain & IO::TIMEDOMAIN) > 0) {
      std::string topSurf = IO::PrefixSurface() + IO::PrefixTop() + IO::PrefixTime();
      std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
      estimation_simbox.setTopBotName(topSurf,baseSurf,output_format);
      if (generate_seismic) {
        estimation_simbox.writeTopBotGrids(topSurf,
                                  baseSurf,
                                  IO::PathToSeismicData(),
                                  output_format);
      }
      else if (!estimation_mode){
        if (output_grids_elastic > 0 || output_grids_other > 0 || output_grids_seismic > 0)
          estimation_simbox.writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToInversionResults(),
                                    output_format);
      }
      if ((output_format & IO::STORM) > 0) { // These copies are only needed with the STORM format
        if ((output_grids_elastic & IO::BACKGROUND) > 0 ||
            (output_grids_elastic & IO::BACKGROUND_TREND) > 0 ||
            (estimation_mode && generate_background)) {
          estimation_simbox.writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToBackground(),
                                    output_format);
        }
        if ((output_grids_other & IO::CORRELATION) > 0) {
          estimation_simbox.writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToCorrelations(),
                                    output_format);
        }
        if ((output_grids_seismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
          estimation_simbox.writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToSeismicData(),
                                    output_format);
        }
        if ((output_grids_other & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
          estimation_simbox.writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToVelocity(),
                                    output_format);
        }
      }
    }
  }



  delete top_surface;
  delete base_surface;
  delete top_surface_flat;
  delete base_surface_flat;
}

void CommonData::SetSurfacesMultipleIntervals(const ModelSettings             * const model_settings,
                                              NRLib::Volume                   & full_inversion_volume,
                                              Simbox                          & estimation_simbox,
                                              int                               nz,
                                              const InputFiles                * input_files,
                                              std::string                     & err_text) {

  // Get interval surface data ------------------------------------------------------------------------------

  LogKit::LogFormatted(LogKit::Low,"\nSetting top and base surfaces for the entire inversion volume:\n");
  bool failed = false;

  int                                       nx                            =  model_settings->getAreaParameters()->GetNx();
  int                                       ny                            =  model_settings->getAreaParameters()->GetNy();
  const std::vector<std::string>            interval_names                =  model_settings->getIntervalNames();
  unsigned int                              n_intervals                   =  interval_names.size();
  const std::map<std::string, std::string>  interval_base_time_surfaces   =  input_files->getIntervalBaseTimeSurfaces();
  // implicit assumption that the top surface is the first one given and the base surface is the last one
  const std::string                         top_surface_file_name         =  input_files->getTimeSurfFile(0);
  const std::string                         base_surface_file_name        =  interval_base_time_surfaces.find( interval_names[n_intervals-1])->second;

  bool                                      generate_seismic              = forward_modeling_;
  bool                                      estimation_mode               = model_settings->getEstimationMode();
  bool                                      generate_background           = model_settings->getGenerateBackground();
  bool                                      parallel_surfaces             = model_settings->getParallelTimeSurfaces();
  //int                                       nz                            = model_settings->getTimeNz();
  int                                       output_format                 = model_settings->getOutputGridFormat();
  int                                       output_domain                 = model_settings->getOutputGridDomain();
  int                                       output_grids_elastic          = model_settings->getOutputGridsElastic();
  int                                       output_grids_other            = model_settings->getOutputGridsOther();
  int                                       output_grids_seismic          = model_settings->getOutputGridsSeismic();
  double                                    d_top                         = model_settings->getTimeDTop();

  Surface * top_surface  = NULL;
  Surface * base_surface = NULL;

  Surface * top_surface_flat  = NULL;
  Surface * base_surface_flat = NULL;

  try{
    if (NRLib::IsNumber(top_surface_file_name)){
      LogKit::LogFormatted(LogKit::Low,"Top surface: Flat surface at depth %11.2f \n", atof(top_surface_file_name.c_str()));
      double x_min, x_max;
      double y_min, y_max;
      FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                  estimation_simbox.getlx(), estimation_simbox.getly(),
                                  estimation_simbox.getAngle(), x_min,y_min,x_max,y_max);
      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, nx, ny, atof(top_surface_file_name.c_str()));
    }
    else{
      LogKit::LogFormatted(LogKit::Low,"Top surface file name: " + top_surface_file_name +" \n");
      top_surface = new Surface(top_surface_file_name);
    }

  }
  catch(NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  if (!failed){
    try{
      if (NRLib::IsNumber(base_surface_file_name)){
        LogKit::LogFormatted(LogKit::Low,"Base surface: Flat surface at depth %11.2f \n", atof(base_surface_file_name.c_str()));
        double x_min, x_max;
        double y_min, y_max;
        FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                    estimation_simbox.getlx(), estimation_simbox.getly(),
                                    estimation_simbox.getAngle(), x_min,y_min,x_max,y_max);
        base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, nx, ny, atof(base_surface_file_name.c_str()));
      }
      else{
        LogKit::LogFormatted(LogKit::Low,"Base surface file name: " + base_surface_file_name +" \n");
        base_surface = new Surface(base_surface_file_name);
      }
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }
  if (!failed){
    try{ // initialize full_inversion_volume and set the flat top and base surfaces of the simbox
      full_inversion_volume =  NRLib::Volume(estimation_simbox.getx0(), estimation_simbox.gety0(),
        estimation_simbox.getlx(), estimation_simbox.getly(),
        *top_surface, *base_surface, estimation_simbox.getAngle());

      // find the top and bottom points on the full inversion volume
      double z_min_top_surface = full_inversion_volume.GetTopZMin(nx,ny);
      double z_max_base_surface = full_inversion_volume.GetBotZMax(nx,ny);
      // round outwards from the center of the volume to nearest multiple 4
      double z_min = floor(z_min_top_surface/4)*4;
      double z_max = ceil(z_max_base_surface/4)*4;
      top_surface_flat = new Surface(top_surface->GetXMin(), top_surface->GetYMin(),
        top_surface->GetXMax()-top_surface->GetXMin(), top_surface->GetYMax()-top_surface->GetYMin(),
        nx, ny, z_min);
      base_surface_flat = new Surface(base_surface->GetXMin(), base_surface->GetYMin(),
        base_surface->GetXMax()-base_surface->GetXMin(), base_surface->GetYMax()-base_surface->GetYMin(),
        model_settings->getAreaParameters()->GetNx(), model_settings->getAreaParameters()->GetNy(), z_max);

      estimation_simbox.setDepth(*top_surface_flat, *base_surface_flat, nz, model_settings->getRunFromPanel());
      estimation_simbox.calculateDz(model_settings->getLzLimit(), err_text);
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }
  if (!failed) {
      if ((output_domain & IO::TIMEDOMAIN) > 0) {
        std::string topSurf = IO::PrefixSurface() + IO::PrefixTop() + IO::PrefixTime();
        std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
        estimation_simbox.setTopBotName(topSurf,baseSurf,output_format);
        if (generate_seismic) {
          estimation_simbox.writeTopBotGrids(topSurf,
                                   baseSurf,
                                   IO::PathToSeismicData(),
                                   output_format);
        }
        else if (!estimation_mode){
          if (output_grids_elastic > 0 || output_grids_other > 0 || output_grids_seismic > 0)
            estimation_simbox.writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToInversionResults(),
                                     output_format);
        }
        if ((output_format & IO::STORM) > 0) { // These copies are only needed with the STORM format
          if ((output_grids_elastic & IO::BACKGROUND) > 0 ||
              (output_grids_elastic & IO::BACKGROUND_TREND) > 0 ||
              (estimation_mode && generate_background)) {
            estimation_simbox.writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToBackground(),
                                     output_format);
          }
          if ((output_grids_other & IO::CORRELATION) > 0) {
            estimation_simbox.writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToCorrelations(),
                                     output_format);
          }
          if ((output_grids_seismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
            estimation_simbox.writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToSeismicData(),
                                     output_format);
          }
          if ((output_grids_other & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
            estimation_simbox.writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToVelocity(),
                                     output_format);
          }
        }
      }
    }

  delete top_surface;
  delete base_surface;
  delete top_surface_flat;
  delete base_surface_flat;

}


bool CommonData::BlockWellsForEstimation(const ModelSettings                            * const model_settings,
                                         const Simbox                                   & estimation_simbox,
                                         const MultiIntervalGrid                        * multiple_interval_grid,
                                         std::vector<NRLib::Well>                       & wells,
                                         std::map<std::string, BlockedLogsCommon *>     & mapped_blocked_logs_common,
                                         std::map<std::string, BlockedLogsCommon *>     & mapped_blocked_logs_for_correlation,
                                         std::string                                    & err_text_common) {

  std::string err_text = "";
  const std::vector<Simbox> interval_simboxes = multiple_interval_grid->GetIntervalSimboxes();
  LogKit::WriteHeader("Blocking wells for estimation");

  // Continuous parameters that are to be used in BlockedLogs
  continuous_logs_to_be_blocked_.push_back("Vp");
  continuous_logs_to_be_blocked_.push_back("Vs");
  continuous_logs_to_be_blocked_.push_back("Rho");
  continuous_logs_to_be_blocked_.push_back("MD");

  // Discrete parameters that are to be used in BlockedLogs


  // Block logs to surrounding estimation simbox
  try{
    for (unsigned int i=0; i<wells.size(); i++) {
      BlockedLogsCommon * blocked_log = new BlockedLogsCommon(&wells[i], continuous_logs_to_be_blocked_, discrete_logs_to_be_blocked_,
                                                              &estimation_simbox, model_settings->getRunFromPanel(), err_text);
      mapped_blocked_logs_common.insert(std::pair<std::string, BlockedLogsCommon *>(wells[i].GetWellName(), blocked_log));


      //blocked_log->FilterLogs(model_settings->getMaxHzBackground(),
      //                        model_settings->getMaxHzSeismic());

    }
  }
  catch(NRLib::Exception & e){
    err_text += e.what();
  }

  // Block logs to interval simboxes for estimation
  if(err_text == ""){
    try{
      for (unsigned int i=0; i<wells.size(); i++){
        mapped_blocked_logs_for_correlation.insert(std::pair<std::string, BlockedLogsCommon *>(wells[i].GetWellName(),
          new BlockedLogsCommon(&wells[i], continuous_logs_to_be_blocked_, discrete_logs_to_be_blocked_, multiple_interval_grid,
                                model_settings->getRunFromPanel(), err_text)));
      }
    }
    catch(NRLib::Exception & e){
    err_text += e.what();
    }
  }

  if(err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}


bool  CommonData::OptimizeWellLocations(ModelSettings                                 * model_settings,
                                        InputFiles                                    * input_files,
                                        const Simbox                                  * estimation_simbox,
                                        //const NRLib::Volume                           & volume,
                                        std::vector<NRLib::Well>                      & wells,
                                        std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs,
                                        std::map<int, std::vector<SeismicStorage> >   & seismic_data,
                                        std::map<int, float **>                       & reflection_matrix,
                                        std::string                                   & err_text_common) {

  std::string err_text = "";

  LogKit::WriteHeader("Estimating optimized well location");

  std::vector<Surface *>      well_move_interval;
  LoadWellMoveInterval(input_files, estimation_simbox, well_move_interval, err_text);

  size_t n_wells = wells.size();

  double  delta_X, delta_Y;
  float   sum;
  float   k_move;
  float   move_angle;
  int     i_move;
  int     j_move;
  int     i_max_offset;
  int     j_max_offset;
  int     n_move_angles = 0;
  int     n_angles     = model_settings->getNumberOfAngles(0);//Well location is not estimated when using time lapse data
  float   max_shift    = model_settings->getMaxWellShift();
  float   max_offset   = model_settings->getMaxWellOffset();
  double  angle       = estimation_simbox->getAngle();
  double  dx          = estimation_simbox->getdx();
  double  dy          = estimation_simbox->getdy();
  std::vector<float> seismicAngle = model_settings->getAngle(0); //Use first time lapse as this not is allowed in 4D

  std::vector<float> angle_weight(n_angles);
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"  Well             Shift[ms]       DeltaI   DeltaX[m]   DeltaJ   DeltaY[m] \n");
  LogKit::LogFormatted(LogKit::Low,"  ----------------------------------------------------------------------------------\n");

  for (int w = 0 ; w < static_cast<int>(n_wells) ; w++) {
    if (wells[w].IsDeviated())
      continue;
    std::string well_name = wells[w].GetWellName();
    std::map<std::string, BlockedLogsCommon * >::iterator it = mapped_blocked_logs.find(well_name);
    if (it == mapped_blocked_logs.end()){
      err_text += "Blocked log not found for well  " + wells[w].GetWellName() + "\n";
      break;
    }
    BlockedLogsCommon * bl = it->second;
    n_move_angles = model_settings->getNumberOfWellAngles(w);

    if ( n_move_angles==0 )
      continue;

    for (int i=0; i<n_angles; i++ )
      angle_weight[i] = 0;

    for (int i=0; i<n_move_angles; i++ ){
      move_angle   = model_settings->getWellMoveAngle(w,i);

      for (int j=0; j<n_angles; j++ ){
        if ( move_angle == seismicAngle[j]){
          angle_weight[j] = model_settings->getWellMoveWeight(w,i);
          break;
        }
      }
    }

    sum = 0;
    for (int i=0; i<n_angles; i++ )
      sum += angle_weight[i];
    if ( sum == 0 )
      continue;

    i_max_offset = static_cast<int>(std::ceil(max_offset/dx));
    j_max_offset = static_cast<int>(std::ceil(max_offset/dy));

    bl->FindOptimalWellLocation(seismic_data[0], estimation_simbox, reflection_matrix[0], n_angles,angle_weight,
                                max_shift,i_max_offset,j_max_offset, well_move_interval,i_move,j_move,k_move);

    delta_X = i_move*dx*cos(angle) - j_move*dy*sin(angle);
    delta_Y = i_move*dx*sin(angle) + j_move*dy*cos(angle);
    MoveWell(wells_[w], estimation_simbox,delta_X,delta_Y,k_move);
    // delete old blocked well and create new
    delete bl;
    mapped_blocked_logs.erase(it);
    mapped_blocked_logs.insert(std::pair<std::string, BlockedLogsCommon *>(well_name, new BlockedLogsCommon(&wells[w], continuous_logs_to_be_blocked_, discrete_logs_to_be_blocked_,
                                                                                                            estimation_simbox, model_settings->getRunFromPanel(), err_text) ) );
    LogKit::LogFormatted(LogKit::Low,"  %-13s %11.2f %12d %11.2f %8d %11.2f \n",
    wells[w].GetWellName().c_str(), k_move, i_move, delta_X, j_move, delta_Y);
  }

  for (int w = 0 ; w < static_cast<int>(n_wells) ; w++){
    n_move_angles = model_settings->getNumberOfWellAngles(w);

    if ( wells[w].IsDeviated()==true && n_move_angles > 0 ) {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Well %7s is treated as deviated and can not be moved.\n",
                                           wells[w].GetWellName().c_str());
      TaskList::addTask("Well "+NRLib::ToString(wells[w].GetWellName())+" can not be moved. Remove <optimize-location-to> for this well");
    }
  }

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

void CommonData::MoveWell(NRLib::Well  & well,
                          const Simbox * simbox,
                          double         delta_X,
                          double         delta_Y,
                          double         k_move){

  double delta_Z;
  double top_old, top_new;

  if (well.HasContLog("X_pos") && well.HasContLog("Y_pos")){
    std::vector<double> & x_pos = well.GetContLog("X_pos");
    std::vector<double> & y_pos = well.GetContLog("Y_pos");
    std::vector<double> & z_pos = well.GetContLog("Z_pos");

    top_old = simbox->getTop(x_pos[0], y_pos[0]);

    for (unsigned int i=0; i<x_pos.size(); i++)
    {
      if (x_pos[i] != RMISSING)
        x_pos[i] = x_pos[i]+delta_X;
      if (y_pos[i] != RMISSING)
        y_pos[i] = y_pos[i]+delta_Y;
    }

    top_new = simbox->getTop(x_pos[0], y_pos[0]);

    delta_Z = top_new - top_old + k_move;

    for (unsigned int i=0; i<z_pos.size(); i++)
      z_pos[i] = z_pos[i]+delta_Z;
  }
}

void  CommonData::CalculateDeviation(NRLib::Well            & new_well,
                                     const ModelSettings    * const model_settings,
                                     float                  & dev_angle,
                                     Simbox                 * simbox)
{
  float maxDevAngle   = model_settings->getMaxDevAngle();
  float thr_deviation = float(tan(maxDevAngle*NRLib::Pi/180.0));  // Largest allowed deviation
  float max_deviation =  0.0f;
  float max_dz        = 10.0f;                      // Calculate slope each 10'th millisecond.

  const std::vector<double> & x_pos = new_well.GetContLog("X_pos");
  const std::vector<double> & y_pos = new_well.GetContLog("Y_pos");
  const std::vector<double> & z_pos = new_well.GetContLog("Z_pos");

  //
  // Find first log entry in simbox
  //
  int iFirst = IMISSING;
  for (unsigned int i=0 ; i < x_pos.size() ; i++)
  {
    if (simbox->isInside(x_pos[i], y_pos[i]))
    {
      if (z_pos[i] > simbox->getTop(x_pos[i], y_pos[i]))
      {
        iFirst = i;
        break;
      }
    }
  }

  if (iFirst != IMISSING) {
    //
    // Find last log entry in simbox
    //
    int iLast = iFirst;
    for (unsigned int i = iFirst + 1 ; i < x_pos.size() ; i++)
    {
      if (simbox->isInside(x_pos[i], y_pos[i]))
        {
          if (z_pos[i] > simbox->getBot(x_pos[i], y_pos[i]))
            break;
        }
      else
        break;
      iLast = i;
    }

    if (iLast > iFirst) {
      double x0 = x_pos[iFirst];
      double y0 = y_pos[iFirst];
      double z0 = z_pos[iFirst];
      for (int i = iFirst+1 ; i < iLast+1 ; i++) {
        double x1 = x_pos[i];
        double y1 = y_pos[i];
        double z1 = z_pos[i];
        float dz = static_cast<float>(z1 - z0);

        if (dz > max_dz || i == iLast) {
          float deviation = static_cast<float>(sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0))/dz);
          if (deviation > max_deviation) {
            x0 = x1;
            y0 = y1;
            z0 = z1;
            max_deviation = deviation;
          }
        }
      }
    }
    dev_angle = static_cast<float>(atan(max_deviation)*180.0/NRLib::Pi);
    LogKit::LogFormatted(LogKit::Low,"   Maximum local deviation is %.1f degrees.",dev_angle);

    if (max_deviation > thr_deviation)
    {
      if (new_well.GetUseForWaveletEstimation() == ModelSettings::NOTSET) {
        new_well.SetUseForWaveletEstimation(ModelSettings::NO);
      }
      //if (use_for_wavelet_estimation == ModelSettings::NOTSET) {
      //  use_for_wavelet_estimation = ModelSettings::NO;
      //}
      new_well.SetDeviated(true);
      LogKit::LogFormatted(LogKit::Low," Well is treated as deviated.\n");
    }
    else
    {
      new_well.SetDeviated(false);
      LogKit::LogFormatted(LogKit::Low,"\n");
    }
  }
  else {
    new_well.SetDeviated(false);
    LogKit::LogFormatted(LogKit::Low,"Well is outside inversion interval. Cannot calculate deviation.\n");
  }

}

void CommonData::LoadWellMoveInterval(const InputFiles             * input_files,
                                      const Simbox                 * estimation_simbox,
                                      std::vector<Surface *>       & well_move_interval,
                                      std::string                  & err_text) {

  const double x0 = estimation_simbox->getx0();
  const double y0 = estimation_simbox->gety0();
  const double lx = estimation_simbox->getlx();
  const double ly = estimation_simbox->getly();
  const int    nx = estimation_simbox->getnx();
  const int    ny = estimation_simbox->getny();
  //
  // Get well move interval
  //
  const std::string & topWMI  = input_files->getWellMoveIntFile(0);
  const std::string & baseWMI = input_files->getWellMoveIntFile(1);

  if (topWMI != "" && baseWMI != "") {
    well_move_interval.resize(2);
    try {
      if (NRLib::IsNumber(topWMI))
        well_move_interval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topWMI.c_str()));
      else {
        Surface tmpSurf(topWMI);
        well_move_interval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      err_text += e.what();
    }

    try {
      if (NRLib::IsNumber(baseWMI))
        well_move_interval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseWMI.c_str()));
      else {
        Surface tmpSurf(baseWMI);
        well_move_interval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      err_text += e.what();
    }
  }
}

bool CommonData::SetupTrendCubes(ModelSettings                  * model_settings,
                                 InputFiles                     * input_files,
                                 MultiIntervalGrid              * multiple_interval_grid,
                                 std::string                    & err_text_common) {

  LogKit::WriteHeader("Setting up Trend Cubes");

  std::string err_text = "";
  //std::vector<CravaTrend> trend_cubes;

  // Get trend variables from model settings
  const std::vector<std::string>  & trend_cube_parameters     = model_settings->getTrendCubeParameters();
  const std::vector<int>          & trend_cube_type           = model_settings->getTrendCubeType();
  //trend_cubes_.resize(multiple_interval_grid->GetNIntervals());
  const std::vector<std::string>  & interval_names            =  model_settings->getIntervalNames();

  // Initialize values
  std::vector<std::vector<NRLib::Grid<double> *> > trend_cubes_tmp;  //Vector(trends) vector(intervals)
  trend_cubes_tmp.resize(trend_cube_parameters.size());

  for (size_t i = 0; i < trend_cube_parameters.size(); i++) {
    trend_cubes_tmp[i].resize(interval_names.size());
    for (size_t j = 0; j < interval_names.size(); j++) {
      trend_cubes_tmp[i][j] = new NRLib::Grid<double>();
    }
  }

  try{

    for (size_t i = 0; i< trend_cube_parameters.size(); i++) {
      if (trend_cube_type[i] == ModelSettings::CUBE_FROM_FILE) {
        // 1. Read the file into an FFTGrid
        const std::string         log_name   = "Trend cube '"+trend_cube_parameters[i]+"'";
        const SegyGeometry      * dummy1     = NULL;
        const TraceHeaderFormat * dummy2     = NULL;
        const float               offset     = model_settings->getSegyOffset(0); //Facies estimation only allowed for one time lapse
        std::string err_text_tmp             = "";

        ReadGridFromFile(input_files->getTrendCube(i),
                         log_name,
                         offset,
                         trend_cubes_tmp[i],
                         dummy1,
                         dummy2,
                         PARAMETER,
                         multiple_interval_grid->GetIntervalSimboxes(),
                         estimation_simbox_,
                         model_settings,
                         err_text_tmp,
                         true);

        if (err_text_tmp != "") {
          err_text += "Reading of file \'"+input_files->getTrendCube(i)+"\' failed\n";
          err_text += err_text_tmp;
        }
      }
    }

    std::vector<CravaTrend> trend_cubes;
    trend_cubes.resize(multiple_interval_grid->GetNIntervals());

    for (int i = 0; i<multiple_interval_grid->GetNIntervals(); i++) {

      std::vector<NRLib::Grid<double> *> trend_cubes_interval;
      for(size_t j = 0; j < trend_cube_parameters.size(); j++) {
        trend_cubes_interval.push_back(trend_cubes_tmp[j][i]);
      }

      trend_cubes[i] = CravaTrend(multiple_interval_grid->GetIntervalSimbox(i),
                                  model_settings,
                                  input_files,
                                  interval_names[i],
                                  trend_cube_type,
                                  trend_cube_parameters,
                                  trend_cubes_interval,
                                  err_text);

    }

    //Trend cubes stored in multiple interval grid
    multiple_interval_grid->AddTrendCubes(trend_cubes);
    n_trend_cubes_ = trend_cubes.size(); //Is this needed?

  }catch(NRLib::Exception & e){
    err_text += e.what();
   }

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

bool CommonData::SetupRockPhysics(const ModelSettings                               * model_settings,
                                  const InputFiles                                  * input_files,
                                  const MultiIntervalGrid                           * multiple_interval_grid,
                                  const std::vector<CravaTrend>                     & trend_cubes,
                                  const std::map<std::string, BlockedLogsCommon *>  & mapped_blocked_logs,
                                  int                                                 n_trend_cubes,
                                  std::string                                       & err_text_common) {

  LogKit::WriteHeader("Processing Rock Physics");

  (void) multiple_interval_grid;
  (void) mapped_blocked_logs;

  std::string err_text = "";

  // rock physics data
  const std::vector<std::string>                    interval_names          = model_settings->getIntervalNames();
  //int                                               n_vintages              = model_settings->getNumberOfVintages();
  const std::string                                 path                    = input_files->getInputDirectory();
  const std::vector<std::string>                    trend_cube_parameters   = model_settings->getTrendCubeParameters();
  std::vector<std::vector<std::vector<double> > >   trend_cube_sampling(n_trend_cubes_);
  const std::vector<std::vector<float> >            dummy_blocked_logs;
  const std::map<std::string,
    std::vector<DistributionWithTrendStorage *> >   reservoir_variable      = model_settings->getReservoirVariable();

  //
  //reservoir_variables_.resize(n_trend_cubes);
  //rock_distributions_.resize(n_trend_cubes);

  // generate distribution for each reservoir variable
  for (int i=0; i < n_trend_cubes; i++) { // the number of trend cubes is the same as the number of intervals
    trend_cube_sampling[i]                                                                     = trend_cubes[i].GetTrendCubeSampling();
    for (std::map<std::string, std::vector<DistributionWithTrendStorage *> >::const_iterator it = reservoir_variable.begin(); it != reservoir_variable.end(); it++) {

      std::vector<DistributionWithTrendStorage *>   storage = it->second;
      std::vector<DistributionWithTrend *>          dist_vector(storage.size());

      for (size_t j=0; j<storage.size(); j++) {
        dist_vector[j] = storage[j]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling[i], err_text);
      }

      reservoir_variables_[it->first]  = dist_vector; //H Removed intervals for reservoir_variables_ -> Still use trend_cube_sampling per interval/trend_cube?
    }
  }


  if (err_text == "") {

    // elastic min/max data
    float alpha_min     = model_settings->getAlphaMin();
    float alpha_max     = model_settings->getAlphaMax();
    float beta_min      = model_settings->getBetaMin();
    float beta_max      = model_settings->getBetaMax();
    float rho_min       = model_settings->getRhoMin();
    float rho_max       = model_settings->getRhoMax();
    float var_alpha_min = model_settings->getVarAlphaMin();
    float var_alpha_max = model_settings->getVarAlphaMax();
    float var_beta_min  = model_settings->getVarBetaMin();
    float var_beta_max  = model_settings->getVarBetaMax();
    float var_rho_min   = model_settings->getVarRhoMin();
    float var_rho_max   = model_settings->getVarRhoMax();
    //int   n_wells       = model_settings->getNumberOfWells();

    // Block logs
    /*
    std::vector<BlockedLogsCommon *> blocked_logs_rock_physics(n_wells, NULL);
    if (n_wells > 0) {
      for (int i=0 ; i<n_wells ; i++)
        blocked_logs_rock_physics[i] = new BlockedLogsCommon();
    }
    */

    // map between reservoir variables and storage classes
    const std::map<std::string, DistributionsFluidStorage   *>  fluid_storage    = model_settings->getFluidStorage();
    const std::map<std::string, DistributionsSolidStorage   *>  solid_storage    = model_settings->getSolidStorage();
    const std::map<std::string, DistributionsDryRockStorage *>  dry_rock_storage = model_settings->getDryRockStorage();
    const std::map<std::string, DistributionsRockStorage    *>  rock_storage     = model_settings->getRockStorage();

    // Map reservoir variables for use in rocks to access resampling trigger.
    std::vector<std::vector<DistributionWithTrend *> > res_var_vintage;
    //for (int i=0; i<n_trend_cubes_; i++){
      if (reservoir_variables_.size() > 0) {
        size_t n_vintages = reservoir_variables_.begin()->second.size();
        res_var_vintage.resize(n_vintages);
        for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator var_it = reservoir_variables_.begin();
          var_it != reservoir_variables_.end();var_it++)
        {
          for (size_t vin_index=0; vin_index < var_it->second.size();vin_index++)
            res_var_vintage[vin_index].push_back((var_it->second)[vin_index]);
        }
      }
    //}

    // Facies names
    std::map<std::string, std::map<std::string, float> > prior_facies_prob_interval = model_settings->getPriorFaciesProbIntervals();
    //std::map<std::string, float> facies_probabilities = model_settings->getPriorFaciesProb();
    std::map<std::string, std::string> facies_cubes   = input_files->getPriorFaciesProbFile();
    std::vector<std::string> all_facies_names         = facies_names_;

    // EN hva skjer her?
    /*
    for (std::map<std::string, float>::iterator it_prob = facies_probabilities.begin(); it_prob != facies_probabilities.end(); it_prob++)
      all_facies_names.push_back(it_prob->first);
    for (std::map<std::string, std::string>::iterator it_cube = facies_cubes.begin(); it_cube != facies_cubes.end(); it_cube++)
      all_facies_names.push_back(it_cube->first);
    */

    std::sort(all_facies_names.begin(), all_facies_names.end());

    std::string prev_facies = "";

    for (size_t i = 0; i<interval_names.size(); i++){
      for (size_t f=0; f<all_facies_names.size(); f++) {
        if (all_facies_names[f] != prev_facies) {
          prev_facies = all_facies_names[f];

          std::map<std::string, DistributionsRockStorage *>::const_iterator iter = rock_storage.find(all_facies_names[f]);
          if (iter != rock_storage.end()) {

            std::string rock_err_txt = "";

            std::string name = iter->first;
            LogKit::LogFormatted(LogKit::Low, "\nRock \'"+name+"\':\n");

            DistributionsRockStorage * storage    = iter->second;
            // rock: vector of size equal to n_vintages
            /*std::vector<DistributionsRock *> rock = storage->GenerateDistributionsRock(n_vintages,
                                                                                        path,
                                                                                        trend_cube_parameters,
                                                                                        trend_cube_sampling[i],
                                                                                        blocked_logs_rock_physics[i],
                                                                                        rock_storage,
                                                                                        solid_storage,
                                                                                        dry_rock_storage,
                                                                                        fluid_storage,
                                                                                        rock_err_txt);*/

            if (rock_err_txt == "") {

              int n_vintages =  1;//static_cast<int>(rock.size());
              if (n_vintages > 1)
                LogKit::LogFormatted(LogKit::Low, "Number of vintages: %4d\n", n_vintages);

              for (int t=0; t<n_vintages; t++) {
                if (n_vintages > 1)
                  LogKit::LogFormatted(LogKit::Low, "\nVintage number: %4d\n", t+1);

                //Completing the top level rocks, by setting access to reservoir variables and sampling distribution.
                //rock[t]->CompleteTopLevelObject(res_var_vintage[i][t]);

                std::vector<bool> has_trends (1,false);// = rock[t]->HasTrend();
                bool              has_trend = false;
                for (size_t j=0; j<has_trends.size(); j++) {
                  if (has_trends[j] == true) {
                    has_trend = true;
                    break;
                  }
                }

                std::vector<double> expectation  (1,0);//rock[t]->GetMeanLogExpectation();
                NRLib::Grid2D<double> covariance  (1,0);//rock[t]->GetMeanLogCovariance();

                PrintExpectationAndCovariance(expectation, covariance, has_trend);

                std::string tmp_err_txt = "";
                if (std::exp(expectation[0]) < alpha_min  || std::exp(expectation[0]) > alpha_max) {
                  tmp_err_txt += "Vp value of "+NRLib::ToString(std::exp(expectation[0]))+" detected: ";
                  tmp_err_txt += "Vp should be in the interval ("+NRLib::ToString(alpha_min)+", "+NRLib::ToString(alpha_max)+") m/s\n";
                }
                if (std::exp(expectation[1]) < beta_min  || std::exp(expectation[1]) > beta_max) {
                  if (typeid(*(storage)) == typeid(ReussRockStorage))
                    tmp_err_txt += "Vs value of 0 detected. Note that the Reuss model gives Vs=0; hence it can not be used to model a facies\n";
                  else
                    tmp_err_txt += "Vs value of "+NRLib::ToString(std::exp(expectation[1]))+" detected: ";
                  tmp_err_txt += "Vs should be in the interval ("+NRLib::ToString(beta_min)+", "+NRLib::ToString(beta_max)+") m/s\n";
                }
                if (std::exp(expectation[2]) < rho_min  || std::exp(expectation[2]) > rho_max) {
                  tmp_err_txt += "Rho value of "+NRLib::ToString(std::exp(expectation[2]))+" detected: ";
                  tmp_err_txt += "Rho should be in the interval ("+NRLib::ToString(rho_min)+", "+NRLib::ToString(rho_max)+") g/cm^3\n";
                }

                if (tmp_err_txt != "") {
                  err_text += "\nToo high or low seismic properties calculated for rock '"+iter->first+"':\n";
                  err_text += tmp_err_txt;
                }

                std::string var_err_txt = "";
                if (covariance(0,0) < var_alpha_min  || covariance(0,0) > var_alpha_max) {
                  var_err_txt += "Var(log Vp) value of "+NRLib::ToString(covariance(0,0))+" detected: ";
                  var_err_txt += "Var(log Vp) should be in the interval ("+NRLib::ToString(var_alpha_min)+", "+NRLib::ToString(var_alpha_max)+")\n";
                }
                if (covariance(1,1) < var_beta_min  || covariance(1,1) > var_beta_max) {
                  var_err_txt += "Var(log Vs) value of "+NRLib::ToString(covariance(1,1))+" detected: ";
                  var_err_txt += "Var(log Vs) should be in the interval ("+NRLib::ToString(var_beta_min)+", "+NRLib::ToString(var_beta_max)+")\n";
                }
                if (covariance(2,2) < var_rho_min  || covariance(2,2) > var_rho_max) {
                  var_err_txt += "Var(log Rho) value of "+NRLib::ToString(covariance(2,2))+" detected: ";
                  var_err_txt += "Var(log Rho) should be in the interval ("+NRLib::ToString(var_rho_min)+", "+NRLib::ToString(var_rho_max)+")\n";
                }

                if (var_err_txt != "") {
                  err_text += "\nToo high or low variance of seismic properties calculated for rock '"+iter->first+"':\n";
                  err_text += var_err_txt;
                }

                // Check correlations
                float corr01 = static_cast<float>(covariance(0,1)/(sqrt(covariance(0,0)*covariance(1,1))));
                float corr02 = static_cast<float>(covariance(0,2)/(sqrt(covariance(0,0)*covariance(2,2))));
                float corr12 = static_cast<float>(covariance(1,2)/(sqrt(covariance(1,1)*covariance(2,2))));

                NRLib::SymmetricMatrix corr(3);
                corr(0,0) = 1;
                corr(1,1) = 1;
                corr(2,2) = 1;
                corr(0,1) = corr01;
                corr(0,2) = corr02;
                corr(1,2) = corr12;

                try {
                  NRLib::CholeskyInvert(corr);
                }
                catch (NRLib::Exception & e) {
                  err_text += e.what();
                  err_text += " for rock '"+iter->first+"':\n";
                  err_text += "  The variables in the rock model are probably linearly dependent\n";
                }

                if (var_err_txt != "" || tmp_err_txt != "")
                  break;

              }
              // set rock distribution for interval i and facies f
              //rock_distributions_[interval_names[i]][all_facies_names[f]] = rock;
            }
            else
              err_text += rock_err_txt;
          }
          else
            err_text += "The facies "+all_facies_names[f]+" is not one of the rocks in the rock physics model\n";
          //rock_distributions_[i][all_facies_names[f]] = rock;
        }
      }
    }
    /*
    for (int i=0; i<n_wells; i++)
        delete blocked_logs_rock_physics[i];
        */
  }


  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

bool CommonData::SetupPriorFaciesProb(ModelSettings * model_settings,
                                      InputFiles    * input_files,
                                      std::string   & err_text_common) {

  //std::vector<Surface *> outer_facies_estim_interval; // Whole facies interval for all zones/intervals
  std::string err_text = "";

  if (model_settings->getEstimateFaciesProb() || model_settings->getDo4DInversion()) {
    LogKit::WriteHeader("Prior Facies Probabilities");

    if (facies_names_.size() == 0)
      SetFaciesNamesFromRockPhysics();

    std::string tmp_err_text = "";
    CheckFaciesNamesConsistency(model_settings,
                                input_files,
                                tmp_err_text);

    if (tmp_err_text != "") {
      err_text += "Prior facies probabilities failed.\n"+tmp_err_text;
    }

    int n_facies = static_cast<int>(facies_names_.size());

    tmp_err_text = "";
    FindFaciesEstimationInterval(input_files, facies_estim_interval_, tmp_err_text);

    if (tmp_err_text != "") {
      err_text += "Reading facies estimation interval failed.\n"+tmp_err_text;
    }

    if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_WELLS) {
      if (n_facies > 0) {
        prior_facies_.resize(multiple_interval_grid_->GetNIntervals());

        for (int i_interval = 0; i_interval < multiple_interval_grid_->GetNIntervals(); i_interval++) {

          int   nz      = estimation_simbox_.getnz();
          float dz      = static_cast<float>(estimation_simbox_.getdz());
          int   n_wells = model_settings->getNumberOfWells();
          int   n_data  = n_wells*nz;

          int ** facies_count = new int * [n_wells];
          for (int w = 0 ; w < n_wells ; w++)
            facies_count[w] = new int[n_facies];

          for (int w = 0 ; w < n_wells ; w++)
            for (int i = 0 ; i < n_facies ; i++)
              facies_count[w][i] = 0;

          int * facies_log = new int[n_data];   // NB! *internal* log numbering (0, 1, 2, ...)
          for (int i = 0 ; i < n_data ; i++)
            facies_log[i] = IMISSING;

          std::vector<double> vt_vp(nz); // vt = vertical trend
          std::vector<double> vt_vs(nz);
          std::vector<double> vt_rho(nz);
          std::vector<int>    vt_facies(nz);

          int n_used_wells = 0;

          int w_well = 0;
          for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = mapped_blocked_logs_.begin(); it != mapped_blocked_logs_.end(); it++) {
            std::map<std::string, BlockedLogsCommon *>::const_iterator iter = mapped_blocked_logs_.find(it->first);

            if (facies_log_wells_[w_well] == true) { // Well has facies log //if (wells[w]->getNFacies() > 0)
              //
              // Note that we use timeSimbox to calculate prior facies probabilities
              // instead of the simbox with parallel top and base surfaces. This
              // will make the prior probabilities slightly different, but that
              // should not be a problem.
              //

              BlockedLogsCommon * blocked_log = iter->second;
              int n_blocks = blocked_log->GetNumberOfBlocks();

              //
              // Set facies data outside facies estimation interval IMISSING
              //
              // If multiinterval: Set facies data outside current interval IMISSING

              std::vector<int> bl_facies_log = blocked_log->GetFaciesBlocked();

              // Outside this interval
              if (multiple_interval_grid_->GetNIntervals() > 1) {
                const std::vector<double> & x_pos = blocked_log->GetXposBlocked();
                const std::vector<double> & y_pos = blocked_log->GetYposBlocked();
                const std::vector<double> & z_pos = blocked_log->GetZposBlocked();

                for (int i = 0 ; i < n_blocks ; i++) {
                  const Simbox * interval_simbox = multiple_interval_grid_->GetIntervalSimbox(i_interval);

                  const double z_top  = interval_simbox->GetTopSurface().GetZ(x_pos[i], y_pos[i]);
                  const double z_base = interval_simbox->GetBotSurface().GetZ(x_pos[i], y_pos[i]);
                  if ( (z_pos[i] - 0.5*dz) < z_top || (z_pos[i] + 0.5*dz) > z_base)
                    bl_facies_log[i] = IMISSING;
                }
              }

              if (facies_estim_interval_.size() > 0) {
                const std::vector<double> & x_pos = blocked_log->GetXposBlocked();
                const std::vector<double> & y_pos = blocked_log->GetYposBlocked();
                const std::vector<double> & z_pos = blocked_log->GetZposBlocked();

                for (int i = 0 ; i < n_blocks ; i++) {
                  const double z_top  = facies_estim_interval_[0]->GetZ(x_pos[i], y_pos[i]);
                  const double z_base = facies_estim_interval_[1]->GetZ(x_pos[i], y_pos[i]);
                  if ( (z_pos[i] - 0.5*dz) < z_top || (z_pos[i] + 0.5*dz) > z_base)
                    bl_facies_log[i] = IMISSING;
                }
              }

              blocked_log->GetVerticalTrend(blocked_log->GetVpBlocked(), vt_vp);
              blocked_log->GetVerticalTrend(blocked_log->GetVsBlocked(), vt_vs);
              blocked_log->GetVerticalTrend(blocked_log->GetRhoBlocked(), vt_rho);
              blocked_log->GetVerticalTrend(&bl_facies_log[0], vt_facies);

              for (int i = 0; i < nz; i++) {
                int facies = 0;
                if (vt_vp[i] != RMISSING && vt_vs[i] != RMISSING && vt_rho[i] != RMISSING)
                  facies = vt_facies[i];
                else
                  facies = IMISSING;

                facies_log[w_well*nz + i] = facies;
                if (facies != IMISSING)
                  facies_count[w_well][facies]++;
              }
              n_used_wells++;
            }
            w_well++;
          }

          if (n_used_wells > 0) {
            //
            // Probabilities
            //
            LogKit::LogFormatted(LogKit::Low,"\nFacies distributions for each blocked well: \n");
            LogKit::LogFormatted(LogKit::Low,"\nBlockedWell              ");
            for (int i = 0 ; i < n_facies ; i++)
              LogKit::LogFormatted(LogKit::Low,"%12s ",facies_names_[i].c_str());
            LogKit::LogFormatted(LogKit::Low,"\n");
            for (int i = 0 ; i < 24+13*n_facies ; i++)
              LogKit::LogFormatted(LogKit::Low,"-");
            LogKit::LogFormatted(LogKit::Low,"\n");
            for (int w = 0 ; w < n_wells ; w++) {

              if (facies_log_wells_[w] == true) { // Well has facies log
                float tot = 0.0;
                for (int i = 0 ; i < n_facies ; i++) {
                  tot += static_cast<float>(facies_count[w][i]);
                }

                LogKit::LogFormatted(LogKit::Low,"%-23s ",wells_[w].GetWellName().c_str());
                for (int i = 0; i < n_facies; i++) {
                  float facies_prob = static_cast<float>(facies_count[w][i])/tot;
                  LogKit::LogFormatted(LogKit::Low," %12.4f",facies_prob);
                }
                LogKit::LogFormatted(LogKit::Low,"\n");
              }
            }
            LogKit::LogFormatted(LogKit::Low,"\n");
            //
            // Counts
            //
            LogKit::LogFormatted(LogKit::Medium,"\nFacies counts for each blocked well: \n");

            LogKit::LogFormatted(LogKit::Medium,"\nBlockedWell              ");
            for (int i = 0 ; i < n_facies ; i++)
              LogKit::LogFormatted(LogKit::Medium,"%12s ",facies_names_[i].c_str());
            LogKit::LogFormatted(LogKit::Medium,"\n");
            for (int i = 0 ; i < 24+13*n_facies ; i++)
              LogKit::LogFormatted(LogKit::Medium,"-");
            LogKit::LogFormatted(LogKit::Medium,"\n");
            for (int w = 0 ; w < n_wells ; w++) {
              //
              if (facies_log_wells_[w] == true) {
                float tot = 0.0;
                for (int i = 0 ; i < n_facies ; i++)
                  tot += static_cast<float>(facies_count[w][i]);
                LogKit::LogFormatted(LogKit::Medium,"%-23s ",wells_[w].GetWellName().c_str());
                for (int i = 0 ; i < n_facies ; i++) {
                  LogKit::LogFormatted(LogKit::Medium," %12d",facies_count[w][i]);
                }
                LogKit::LogFormatted(LogKit::Medium,"\n");
              }
            }
            LogKit::LogFormatted(LogKit::Medium,"\n");

            for (int w = 0 ; w < n_wells ; w++)
              delete [] facies_count[w];
            delete [] facies_count;

            //
            // Make prior facies probabilities
            //
            float sum = 0.0f;
            int * n_data_facies = new int[n_facies];
            for (int i=0; i < n_facies; i++)
              n_data_facies[i] = 0;

            for (int i = 0; i < n_data; i++) {
              if (facies_log[i] != IMISSING) {
                n_data_facies[facies_log[i]]++;
              }
            }
            delete [] facies_log;

            for (int i = 0; i < n_facies; i++)
              sum += n_data_facies[i];

            if (sum > 0) {
              LogKit::LogFormatted(LogKit::Low,"Facies probabilities based on all blocked wells:\n\n");
              LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
              LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
              prior_facies_[i_interval].resize(n_facies);
              for (int i = 0; i < n_facies; i++) {
                prior_facies_[i_interval][i] = float(n_data_facies[i])/sum;
                LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",facies_names_[i].c_str(),prior_facies_[i_interval][i]);
              }
            }
            else {
              LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No valid facies log entries have been found\n");
              model_settings->setEstimateFaciesProb(false);
              TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");

            }
            delete [] n_data_facies;
          }
          else {
            LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Estimation of facies probabilites have been requested, but there");
            LogKit::LogFormatted(LogKit::Warning,"\n         are no wells with facies available and CRAVA will therefore not");
            LogKit::LogFormatted(LogKit::Warning,"\n         be able to estimate these probabilities...\n");
            model_settings->setEstimateFaciesProb(false);

            TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
          }
        }//i_intervals
      }
      else {
        LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Estimation of facies probabilites have been requested, but no facies");
        LogKit::LogFormatted(LogKit::Warning,"\n         have been found and CRAVA will therefore not be able to estimate");
        LogKit::LogFormatted(LogKit::Warning,"\n         these probabilities...\n");
        model_settings->setEstimateFaciesProb(false);
        TaskList::addTask("Consider using a well containing facies log entries to be able to estimate facies probabilities.");
      }
    }
    else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE) {
      for (int i_interval = 0; i_interval < multiple_interval_grid_->GetNIntervals(); i_interval++) {

        prior_facies_.resize(multiple_interval_grid_->GetNIntervals());
        typedef std::map<std::string,float> map_type;
        map_type my_map;

        if (model_settings->getIntervalNames().size() > 0) {
          std::string interval_name = model_settings->getIntervalName(i_interval);
          my_map = model_settings->getPriorFaciesProbInterval(interval_name);
        }
        else
          my_map = model_settings->getPriorFaciesProb();

        for (int i = 0; i < n_facies; i++) {
          map_type::iterator iter = my_map.find(facies_names_[i]);
          if (iter!=my_map.end())
            prior_facies_[i_interval][i] = iter->second;
          else {
            LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",facies_names_[i].c_str());
            model_settings->setEstimateFaciesProb(false);
            TaskList::addTask("Check that facies " +NRLib::ToString(facies_names_[i].c_str())+" is given a prior probability in the xml-file");
          }
        }
        LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
        LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
        for (int i=0; i < n_facies; i++) {
          LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",facies_names_[i].c_str(),prior_facies_[i]);
        }
      }//i_interval
    }
    else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES) {
      //If intervals, send in vector of simboxes and grids to ReadPriorFaciesProbCubes. Splitting of intervals is done in ReadGridFromFile.
      std::vector<Simbox> & interval_simboxes = multiple_interval_grid_->GetIntervalSimboxes();
      std::vector<std::vector<NRLib::Grid<double> *> > prior_facies_prob_cubes;

      if (model_settings->getIntervalNames().size() == 0)
        interval_simboxes[0] = &estimation_simbox_;

      std::string err_text_tmp = "";
      ReadPriorFaciesProbCubes(input_files,
                               model_settings,
                               prior_facies_prob_cubes,  //Vector(facies) vector(intervals)
                               interval_simboxes,
                               err_text_tmp);

      if (err_text_tmp == "") {

        //Store as vector(intervals) vector(facies)
        int n_facies = facies_names_.size();
        for (size_t i = 0; i < model_settings->getIntervalNames().size(); i++) {
          multiple_interval_grid_->GetPriorFaciesProbCubesInterval(i).resize(n_facies);

          for (int j = 0; j < n_facies; j++) {
            multiple_interval_grid_->AddPriorFaciesCube(i, j, prior_facies_prob_cubes[j][i]);
          }
        }

        typedef std::map<std::string,std::string> mapType;
        mapType myMap = input_files->getPriorFaciesProbFile();

        LogKit::LogFormatted(LogKit::Low,"Facies         Probability in file\n");
        LogKit::LogFormatted(LogKit::Low,"----------------------------------\n");
        for (mapType::iterator it=myMap.begin(); it!=myMap.end(); it++)
          LogKit::LogFormatted(LogKit::Low,"%-15s %10s\n",(it->first).c_str(),(it->second).c_str());
      }
      else {
        err_text += "Error when reading facies from cubes: \n";
        err_text += err_text_tmp;
      }
    }
  }

  // facies_estim_interval to be read in ModelAVOStatic
  //if (outer_facies_estim_interval.size() == 2) {
  //  if (outer_facies_estim_interval[0] != NULL)
  //    delete outer_facies_estim_interval[0];
  //  if (outer_facies_estim_interval[1] != NULL)
  //    delete outer_facies_estim_interval[1];
  //}

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

void CommonData::FindFaciesEstimationInterval(InputFiles             * input_files,
                                              std::vector<Surface *> & facies_estim_interval,
                                              std::string            & err_text)
{
  //
  // Get facies estimation interval
  //
  const std::string & topFEI  = input_files->getFaciesEstIntFile(0);
  const std::string & baseFEI = input_files->getFaciesEstIntFile(1);
  const double x0             = estimation_simbox_.getx0();
  const double y0             = estimation_simbox_.gety0();
  const double lx             = estimation_simbox_.getlx();
  const double ly             = estimation_simbox_.getly();
  const int    nx             = estimation_simbox_.getnx();
  const int    ny             = estimation_simbox_.getny();

  if (topFEI != "" && baseFEI != "") {
    facies_estim_interval.resize(2);
    try {
      if (NRLib::IsNumber(topFEI))
        facies_estim_interval[0] = new Surface(x0,y0,lx,ly,nx,ny,atof(topFEI.c_str()));
      else {
        Surface tmpSurf(topFEI);
        facies_estim_interval[0] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      err_text += e.what();
    }

    try {
      if (NRLib::IsNumber(baseFEI))
        facies_estim_interval[1] = new Surface(x0,y0,lx,ly,nx,ny,atof(baseFEI.c_str()));
      else {
        Surface tmpSurf(baseFEI);
        facies_estim_interval[1] = new Surface(tmpSurf);
      }
    }
    catch (NRLib::Exception & e) {
      err_text += e.what();
    }
  }
}

void
CommonData::ReadPriorFaciesProbCubes(const InputFiles                                 * input_files,
                                     ModelSettings                                    * model_settings,
                                     std::vector<std::vector<NRLib::Grid<double> *> > & prior_facies_prob_cubes, //Vector(facies) vector(intervals)
                                     const std::vector<Simbox>                        & interval_simboxes,
                                     std::string                                      & err_text)
{
  int n_facies = static_cast<int>(facies_names_.size());
  prior_facies_prob_cubes.resize(n_facies);

  typedef std::map<std::string,std::string> mapType;
  mapType myMap = input_files->getPriorFaciesProbFile();
  for (int i=0; i < n_facies; i++) {
    mapType::iterator iter = myMap.find(facies_names_[i]);

    if (iter!=myMap.end()) {
      const std::string & facies_prob_file = iter->second;
      const SegyGeometry      * dummy1 = NULL;
      const TraceHeaderFormat * dummy2 = NULL;
      const float               offset = model_settings->getSegyOffset(0); //Facies estimation only allowed for one time lapse
      std::string error_text("");

      if (model_settings->getIntervalNames().size() > 0)
        prior_facies_prob_cubes[i].resize(model_settings->getIntervalNames().size());

      for (size_t j = 0; j < model_settings->getIntervalNames().size(); j++)
        prior_facies_prob_cubes[i][j] = new NRLib::Grid<double>();

      ReadGridFromFile(facies_prob_file,
                       "priorfaciesprob",
                       offset,
                       prior_facies_prob_cubes[i],
                       dummy1,
                       dummy2,
                       PARAMETER,
                       interval_simboxes,
                       estimation_simbox_,
                       model_settings,
                       error_text);

      if (error_text != "") {
        error_text += "Reading of file \'"+facies_prob_file+"\' for prior facies probability for facies \'"
                     +facies_names_[i]+"\' failed\n";
        err_text += error_text;
      }
    }
    else
    {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: No prior facies probability found for facies %12s\n",
                           facies_names_[i].c_str());
      TaskList::addTask("Check that facies "+NRLib::ToString(facies_names_[i].c_str())+" is given prior probability in the xml-file");
      model_settings->setEstimateFaciesProb(false);
      break;
    }
  }
}

void
CommonData::CheckFaciesNamesConsistency(ModelSettings     *& model_settings,
                                        const InputFiles   * input_files,
                                        std::string        & tmp_err_text) const
{
  int n_facies = static_cast<int>(facies_names_.size());

  // Compare names in wells with names given in rock physics prior model
  if (rock_distributions_.size() > 0) {
    int n_rocks  = static_cast<int>(rock_distributions_.size());
    if (n_rocks > n_facies)
      tmp_err_text += "Problem with facies logs. The number of rocks in the rock physics prior model is larger than the number of facies found in the wells.\n";
    for (int i=0; i<n_facies; i++) {
      if (rock_distributions_.find(facies_names_[i]) == rock_distributions_.end())
        tmp_err_text += "Problem with facies logs. Facies "+facies_names_[i]+" found in a well is not one of the rocks given in rock physics prior model\n";
    }
  }

  // Compare names in wells with names given in .xml-file
  if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_MODEL_FILE)
  {
    typedef std::map<std::string,float> mapType;
    mapType myMap = model_settings->getPriorFaciesProb();

    for (int i=0;i<n_facies;i++)
    {
      mapType::iterator iter = myMap.find(facies_names_[i]);
      if (iter==myMap.end())
        tmp_err_text += "Problem with facies logs. Facies "+facies_names_[i]+" is not one of the facies given in the xml-file.\n";
    }
  }

  // Compare names in wells with names given as input in proability cubes
  else if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES)
  {
    typedef std::map<std::string,std::string> mapType;
    mapType myMap = input_files->getPriorFaciesProbFile();

    for (int i=0;i<n_facies;i++)
    {
      mapType::iterator iter = myMap.find(facies_names_[i]);
      if (iter==myMap.end())
        tmp_err_text += "Problem with facies logs. Facies "+facies_names_[i]+" is not one of the facies given in the xml-file.\n";
    }
  }
}

void
CommonData::SetFaciesNamesFromRockPhysics()
{
  typedef std::map<std::string, DistributionsRock *> mapType;

  int i = 0;
  for (std::map<std::string, std::vector<DistributionsRock *> >::const_iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
    facies_names_.push_back(it->first);
    facies_nr_.push_back(i);
    i++;
  }
}

void
CommonData::ReadGridFromFile(const std::string                  & file_name,
                             const std::string                  & par_name,
                             const float                          offset,
                             std::vector<NRLib::Grid<double> *> & interval_grids,
                             const SegyGeometry                *& geometry,
                             const TraceHeaderFormat            * format,
                             int                                  grid_type,
                             const std::vector<Simbox>          & interval_simboxes,
                             const Simbox                       & estimation_simbox,
                             const ModelSettings                * model_settings,
                             std::string                        & err_text,
                             bool                                 nopadding) {

  int fileType = IO::findGridType(file_name);

  LogKit::LogFormatted(LogKit::Low,"\nReading grid \'"+par_name+"\' from file "+file_name);

  if (fileType == IO::CRAVA) {

    int nx_pad = estimation_simbox_.GetNXpad();
    int ny_pad = estimation_simbox_.GetNYpad();
    int nz_pad = 0;

    GetZPaddingFromCravaFile(file_name, err_text, nz_pad);

    FFTGrid *  crava_grid = CreateFFTGrid(estimation_simbox_.getnx(),
                                          estimation_simbox_.getny(),
                                          estimation_simbox_.getnz(),
                                          nx_pad,
                                          ny_pad,
                                          nz_pad,
                                          model_settings->getFileGrid());

    crava_grid->setAccessMode(FFTGrid::RANDOMACCESS);
    crava_grid->readCravaFile(file_name, err_text);

    //Intervals not used, no need to resample
    // Alternative: create a readCravaFile for NRLibGrid. Padding?
    if (model_settings->getIntervalNames().size() == 0) {
      int nx = interval_simboxes[0].getnx();
      int ny = interval_simboxes[0].getny();
      int nz = interval_simboxes[0].getnz();

      interval_grids[0] = new NRLib::Grid<double>(nx, ny, nz);

      //Copy elements from fft_grid
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            interval_grids[0]->SetValue(i, j, k, crava_grid->getRealValue(i, j, k));
          }
        }
      }

    }
    else {
      for (size_t i_interval = 0; i_interval < interval_simboxes.size(); i_interval++) {

        int nx = interval_simboxes[i_interval].getnx();
        int ny = interval_simboxes[i_interval].getny();
        int nz = interval_simboxes[i_interval].getnz();

        interval_grids[i_interval]->Resize(nx, ny, nz);

        FFTGrid       * fft_grid_tmp = NULL;
        SegY          * segy_tmp     = NULL;
        StormContGrid * storm_tmp    = NULL;

        int missing_traces_simbox  = 0;
        int missing_traces_padding = 0;
        int dead_traces_simbox     = 0;

        int file_type = IO::findGridType(file_name);
        bool is_segy  = false;
        bool is_storm = false;
        bool scale    = false;

        if (file_type == IO::SEGY)
          is_segy = true;
        else if (file_type == IO::STORM)
          is_storm = true;
        else if (file_type == IO::SGRI) {
          is_storm = true;
          scale    = true;
        }

        FillInData(interval_grids[i_interval],
                   fft_grid_tmp,
                   interval_simboxes[i_interval],
                   storm_tmp,
                   segy_tmp,
                   crava_grid,
                   model_settings->getSmoothLength(),
                   missing_traces_simbox,
                   missing_traces_padding,
                   dead_traces_simbox,
                   grid_type,
                   scale,
                   is_segy,
                   is_storm);

        if (fft_grid_tmp != NULL)
          delete fft_grid_tmp;
        if (storm_tmp != NULL)
          delete storm_tmp;
        if (segy_tmp != NULL)
          delete segy_tmp;

      }
    }

    if (crava_grid != NULL)
      delete crava_grid;

  }
  else if (fileType == IO::SEGY)
    ReadSegyFile(file_name,
                  interval_grids,
                  interval_simboxes,
                  estimation_simbox,
                  model_settings,
                  geometry,
                  grid_type,
                  par_name,
                  offset,
                  format,
                  err_text);
  else if (fileType == IO::STORM)
    ReadStormFile(file_name,
                  interval_grids,
                  grid_type,
                  par_name,
                  interval_simboxes,
                  model_settings,
                  err_text,
                  false);
  else if (fileType == IO::SGRI)
    ReadStormFile(file_name,
                  interval_grids,
                  grid_type,
                  par_name,
                  interval_simboxes,
                  model_settings,
                  err_text,
                  true);
  else {
    err_text += "\nReading of file \'"+file_name+"\' for grid type \'"+par_name+"\'failed. File type not recognized.\n";
  }
}

//void CommonData::ReadCravaFile(NRLib::Grid<double> & grid,
//                               const std::string   & file_name,
//                               std::string         & err_text) {
//  std::string error;
//  try {
//    std::ifstream bin_file;
//    NRLib::OpenRead(bin_file, file_name, std::ios::in | std::ios::binary);
//
//    std::string file_type;
//    getline(bin_file, file_type);
//
//    double dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryInt(bin_file);
//    dummy = NRLib::ReadBinaryInt(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    dummy = NRLib::ReadBinaryDouble(bin_file);
//    //int rnxp = NRLib::ReadBinaryInt(bin_file);
//    //int nyp  = NRLib::ReadBinaryInt(bin_file);
//    //int nzp  = NRLib::ReadBinaryInt(bin_file);
//
//    //if (rnxp != rnxp_ || nyp != nyp_ || nzp != nzp_) {
//    //  LogKit::LogFormatted(LogKit::Low,"\n\nERROR: The grid has different dimensions than the model grid. Check the padding settings");
//    //  LogKit::LogFormatted(LogKit::Low,"\n                rnxp   nyp   nzp");
//    //  LogKit::LogFormatted(LogKit::Low,"\n--------------------------------");
//    //  LogKit::LogFormatted(LogKit::Low,"\nModel grid  :   %4d  %4d  %4d",rnxp_,nyp_,nzp_);
//    //  LogKit::LogFormatted(LogKit::Low,"\nGrid on file:   %4d  %4d  %4d\n",rnxp ,nyp ,nzp );
//    //  binFile.close();
//    //  throw(NRLib::Exception("Grid dimension is wrong for file '"+fileName+"'."));
//    //}
//
//    //createRealGrid(!nopadding);
//    //add_ = !nopadding;
//    //int i;
//
//    int ni = grid.GetNI();
//    int nj = grid.GetNJ();
//    int nk = grid.GetNK();
//    size_t index = 0;
//
//    //int size = ni*nj*nk;
//
//    for (int i = 0; i < ni; i++) {
//      for (int j = 0; j < nj; j++) {
//        for (int k = 0; k < nk; k++) {
//          index = grid.GetIndex(i,j,k);
//          grid(index) = NRLib::ReadBinaryFloat(bin_file);
//        }
//      }
//    }
//
//    //for (i=0;i<rsize_;i++)
//    //  rvalue_[i] = NRLib::ReadBinaryFloat(bin_file);
//
//    bin_file.close();
//  }
//  catch (NRLib::Exception & e) {
//    error = std::string("Error: ") + e.what() + "\n";
//  }
//
//  err_text += error;
//}

void CommonData::GetZPaddingFromCravaFile(const std::string & file_name,
                                          std::string       & err_text,
                                          int               & nz_pad) {

  std::string err_text_tmp;
  try {
    std::ifstream bin_file;
    NRLib::OpenRead(bin_file, file_name, std::ios::in | std::ios::binary);

    std::string file_type;
    getline(bin_file, file_type);

    double dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryInt(bin_file);
    dummy = NRLib::ReadBinaryInt(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryDouble(bin_file);
    dummy = NRLib::ReadBinaryInt(bin_file); //rnxp
    dummy = NRLib::ReadBinaryInt(bin_file);
    nz_pad = NRLib::ReadBinaryInt(bin_file);

   bin_file.close();
  }
  catch (NRLib::Exception & e) {
    err_text_tmp = std::string("Error: ") + e.what() + "\n";
  }

  err_text += err_text_tmp;

}

void
CommonData::ReadSegyFile(const std::string                  & file_name,
                         std::vector<NRLib::Grid<double> *> & interval_grids,
                         const std::vector<Simbox>          & interval_simboxes,
                         const Simbox                       & estimation_simbox, //Simbox to check against SegY-data = estimation_simbox
                         const ModelSettings                * model_settings,
                         const SegyGeometry                *& geometry,
                         int                                  grid_type,
                         const std::string                  & par_name,
                         float                                offset,
                         const TraceHeaderFormat            * format,
                         std::string                        & err_text,
                         bool                                 nopadding) {

  (void) nopadding;
  (void) par_name;

  SegY * segy = NULL;
  bool failed = false;

  try
  {
    //
    // Currently we have only one optional TraceHeaderFormat, but this can
    // be augmented to a list with several formats ...
    //
    if (format == NULL) { //Unknown format
      std::vector<TraceHeaderFormat*> traceHeaderFormats(0);
      if (model_settings->getTraceHeaderFormat() != NULL) {
        traceHeaderFormats.push_back(model_settings->getTraceHeaderFormat());
      }
      segy = new SegY(file_name,
                      offset,
                      traceHeaderFormats,
                      true); // Add standard formats to format search
    }
    else //Known format, read directly.
      segy = new SegY(file_name, offset, *format);

    float guard_zone = model_settings->getGuardZone();

    bool cover_ok = false;
    cover_ok = CheckThatDataCoverGrid(segy,
                                      offset,
                                      estimation_simbox,
                                      guard_zone,
                                      err_text);

    if (cover_ok) {
      bool  only_volume      = true;
      // This is *not* the same as FFT-grid padding. If the padding
      // size is changed from 2*guard_zone, the smoothing done in
      // FFTGrid::smoothTraceInGuardZone() will become incorrect.
      float padding         = 2*guard_zone;
      bool  relative_padding = false;

      segy->ReadAllTraces(&estimation_simbox,
                          padding,
                          only_volume,
                          relative_padding);
      segy->CreateRegularGrid();
    }
    else {
      failed = true;
    }
  }
  catch (NRLib::Exception & e)
  {
    err_text += e.what();
    failed = true;
  }

  if (!failed)
  {

    for (size_t i_interval = 0; i_interval < interval_simboxes.size(); i_interval++) {

      int missing_traces_simbox  = 0;
      int missing_traces_padding = 0;
      int dead_traces_simbox     = 0;

      const SegyGeometry * geo;
      geo = segy->GetGeometry();
      geo->WriteGeometry();
      if (grid_type == DATA)
        geometry = new SegyGeometry(geo);

      int xpad = interval_simboxes[i_interval].getnx();
      int ypad = interval_simboxes[i_interval].getny();
      int zpad = interval_simboxes[i_interval].getnz();

      interval_grids[i_interval]->Resize(xpad, ypad, zpad);

      StormContGrid * stormgrid_tmp = NULL;
      FFTGrid * fft_grid_tmp        = NULL;
      FFTGrid * fft_grid_tmp_2      = NULL;
      FillInData(interval_grids[i_interval],
                 fft_grid_tmp,
                 interval_simboxes[i_interval],
                 stormgrid_tmp,
                 segy,
                 fft_grid_tmp_2,
                 model_settings->getSmoothLength(),
                 missing_traces_simbox,
                 missing_traces_padding,
                 dead_traces_simbox,
                 grid_type);
      if (stormgrid_tmp != NULL)
       delete stormgrid_tmp;
      if (fft_grid_tmp != NULL)
        delete fft_grid_tmp;
      if (fft_grid_tmp_2 != NULL)
        delete fft_grid_tmp_2;

      if (missing_traces_simbox > 0) {
        if (missing_traces_simbox == interval_simboxes[i_interval].getnx()*interval_simboxes[i_interval].getny()) {
          err_text += "Error: Data in file "+file_name+" was completely outside the inversion area.\n";
          failed = true;
        }
        else {
          if (grid_type == PARAMETER) {
            err_text += "Grid in file "+file_name+" does not cover the inversion area.\n";
          }
          else {
            LogKit::LogMessage(LogKit::Warning, "WARNING: "+NRLib::ToString(missing_traces_simbox)
                               +" grid columns are outside the area defined by the seismic data.\n");
            std::string text;
            text += "Check seismic volumes and inversion area: A part of the inversion area is outside\n";
            text += "   the seismic data specified in file \'"+file_name+"\'.";
            TaskList::addTask(text);
          }
        }
      }
      //if (missingTracesPadding > 0) {
      //  int nx     = simbox->getnx();
      //  int ny     = simbox->getny();
      //  int nxpad  = xpad - nx;
      //  int nypad  = ypad - ny;
      //  int nxypad = nxpad*ny + nx*nypad - nxpad*nypad;
      //  LogKit::LogMessage(LogKit::High, "Number of grid columns in padding that are outside area defined by seismic data : "
      //                     +NRLib::ToString(missingTracesPadding)+" of "+NRLib::ToString(nxypad)+"\n");
      //}
      if (dead_traces_simbox > 0) {
        LogKit::LogMessage(LogKit::High, "Number of grid columns with no seismic data (nearest trace is dead) : "
                           +NRLib::ToString(dead_traces_simbox)+" of "+NRLib::ToString(interval_simboxes[i_interval].getnx()*interval_simboxes[i_interval].getny())+"\n");
      }
    } //i_interval
  }
  if (segy != NULL)
    delete segy;
}

void CommonData::FillInData(NRLib::Grid<double> * grid_new,
                            FFTGrid             * fft_grid_new,
                            const Simbox        & simbox,
                            StormContGrid       * storm_grid,
                            const SegY          * segy,
                            FFTGrid             * fft_grid_old,
                            float                 smooth_length,
                            int                 & missing_traces_simbox,
                            int                 & missing_traces_padding,
                            int                 & dead_traces_simbox,
                            int                   grid_type,
                            bool                  scale,
                            bool                  is_segy,
                            bool                  is_storm)
{
  assert(grid_type != CTMISSING);

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  float scalevert = 1.0f;
  float scalehor  = 1.0f;

  if (scale == true) { //from sgri file
    LogKit::LogFormatted(LogKit::Low,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
    scalevert = 0.001f; //1000.0;
    scalehor  = 0.001f; //1000.0;
  }

  bool is_nrlib_grid = true;
  if (grid_new->GetN() == 0) //Send in an empty NRLib::Grid if we want fft-grid
    is_nrlib_grid = false;

  int nx   = grid_new->GetNI();
  int ny   = grid_new->GetNJ();
  int nz   = grid_new->GetNK();
  int nxp  = grid_new->GetNI();
  int nyp  = grid_new->GetNJ();
  int nzp  = grid_new->GetNK();
  int rnxp = grid_new->GetNI(); //2*(nxp/2+1);

  if (is_nrlib_grid == false) {
    nx   = fft_grid_new->getNx();
    ny   = fft_grid_new->getNy();
    nz   = fft_grid_new->getNz();
    nxp  = fft_grid_new->getNxp();
    nyp  = fft_grid_new->getNyp();
    nzp  = fft_grid_new->getNzp();
    rnxp = fft_grid_new->getRNxp();
  }

  LogKit::LogFormatted(LogKit::Low,"\nResampling data into %dx%dx%d grid:", nxp, nyp, nzp);

  float monitorSize = std::max(1.0f, static_cast<float>(nyp*rnxp)*0.02f);
  float nextMonitor = monitorSize;
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
  printf("\n  ^");

  //
  // Find proper length of time samples to get N*log(N) performance in FFT.
  //
  size_t n_samples = 0;
  float  dz_data   = 0.0f;
  float  dz_min    = 0.0f;
  if (is_segy) {
    n_samples = segy->FindNumberOfSamplesInLongestTrace();
    dz_data   = segy->GetDz();
    dz_min    = dz_data/4.0f;
  }
  else if (is_storm) {
    n_samples = storm_grid->GetNK();
  }
  else
    n_samples = fft_grid_old->getNzp();

  int nt = FindClosestFactorableNumber(static_cast<int>(n_samples));
  int mt = 4*nt;           // Use four times the sampling density for the fine-meshed data

  //
  // Create FFT plans
  //
  rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  //
  // Do resampling
  //
  missing_traces_simbox  = 0; // Part of simbox is outside seismic data
  missing_traces_padding = 0; // Part of padding is outside seismic data
  dead_traces_simbox     = 0; // Simbox is inside seismic data but trace is missing

  for (int j = 0 ; j < nyp ; j++) {
    for (int i = 0 ; i < rnxp ; i++) {

      int refi = GetFillNumber(i, nx, nxp); // Find index (special treatment for padding)
      int refj = GetFillNumber(j, ny, nyp); // Find index (special treatment for padding)
      int refk = 0;

      double x  = 0.0;
      double y  = 0.0;
      double z0 = 0.0;
      simbox.getCoord(refi, refj, refk, x, y, z0);  // Get lateral position and z-start (z0)
      x  *= scalehor;
      y  *= scalehor;
      z0 *= scalevert;

      double dz = simbox.getdz(refi, refj)*scalevert;
      float  xf = static_cast<float>(x);
      float  yf = static_cast<float>(y);

      bool is_inside = false;
      if (is_segy)
        is_inside = segy->GetGeometry()->IsInside(xf, yf);
      else if (is_storm) {
        if (storm_grid->IsInside(xf, yf) == 1)
          is_inside = true;
      }
      else {
        if (fft_grid_old->getRealValue(refi, refj, 0) != RMISSING) //H Returns missing if outside grid or rvalue_[index] is missing. How to check against xf and yf?
          is_inside = true;
      }

      if (is_inside == true) {
        bool  missing = true;
        float z0_data = RMISSING;

        std::vector<float> data_trace;
        size_t grid_i = 0;
        size_t grid_j = 0;
        double grid_x = 0.0;
        double grid_y = 0.0;
        double grid_z = 0.0;
        float  value  = 0.0f;

        //Get data_trace for this i and j.
        if (is_segy) {
          segy->GetNearestTrace(data_trace, missing, z0_data, xf, yf);
        }
        else if (is_storm) {
          float z_min = 0.0f;
          float z_max = 0.0f;

          storm_grid->FindXYIndex(xf, yf, grid_i, grid_j);
          for (size_t k = 0; k < grid_new->GetNK(); k++) {
            storm_grid->FindCenterOfCell(grid_i, grid_j, k, grid_x, grid_y, grid_z);
            value = storm_grid->GetValueZInterpolated(grid_x, grid_y, grid_z);
            data_trace.push_back(value);

            if (k == 0)
              z_min = static_cast<float>(grid_z);
            if (k == grid_new->GetNK()-1)
              z_max = static_cast<float>(grid_z);
          }

          dz_data = (z_max- z_min) / storm_grid->GetNK();
          dz_min = dz_data/4.0f;
          z0_data = z_min;
        }
        else {
          data_trace = fft_grid_old->getRealTrace2(refi, refj);

          double z_min = 0.0;
          double z_max = 0.0;

          simbox.getZCoord(0, xf, yf, z_min); //H FFTGrid doesnt have top/bot surfaces. Correct to use simbox?
          simbox.getZCoord(fft_grid_old->getNz(), xf, yf, z_max);

          dz_data = (z_max- z_min) / fft_grid_old->getNz();
          dz_min = dz_data/4.0f;
          z0_data = static_cast<float>(z_min);
        }

        size_t n_trace = data_trace.size();
        float trend_first = 0.0f;
        float trend_last  = 0.0f;

        if (grid_type != DATA) {
          //Remove zeroes. F.ex. background on segy-format with a non-constant top-surface, the vector is filled with zeroes at the beginning.
          if (data_trace[0] == 0) {
            std::vector<float> data_trace_new;
            for (size_t k_trace = 0; k_trace < n_trace; k_trace++) {
              if (data_trace[k_trace] != 0)
                data_trace_new.push_back(data_trace[k_trace]);
            }
            data_trace = data_trace_new;
            n_trace = data_trace.size();
          }

          n_samples = data_trace.size();
          nt = FindClosestFactorableNumber(static_cast<int>(n_samples));
          mt = 4*nt;

          fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
          fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

          //Remove trend from trace
          trend_first = data_trace[0];
          trend_last = data_trace[n_trace - 1];
          float trend_inc = (trend_last - trend_first) / (n_trace - 1);
          for (size_t k_trace = 0; k_trace < data_trace.size(); k_trace++) {
            data_trace[k_trace] -= trend_first + k_trace * trend_inc;
          }
        }

        if (is_segy == false || (is_segy == true && !missing)) {
          int         cnt      = nt/2 + 1;
          int         rnt      = 2*cnt;
          int         cmt      = mt/2 + 1;
          int         rmt      = 2*cmt;

          fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
          fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));

          float       dz_grid  = static_cast<float>(dz);
          float       z0_grid  = static_cast<float>(z0);

          std::vector<float> grid_trace(nzp);

          smooth_length*=scalevert;
          if (grid_type == DATA) {
            SmoothTraceInGuardZone(data_trace,
                                   dz_data,
                                   smooth_length);
          }

          ResampleTrace(data_trace,
                        fftplan1,
                        fftplan2,
                        rAmpData,
                        rAmpFine,
                        cnt,
                        rnt,
                        cmt,
                        rmt);

          std::vector<float> data_trace_trend_long;
          if (grid_type != DATA) {
            float trend_inc = (trend_last - trend_first) / (rmt - 1);

            data_trace_trend_long.resize(rmt);
            for (int k_trace = 0; k_trace < rmt; k_trace++) {
              data_trace_trend_long[k_trace] = trend_first + k_trace * trend_inc;
            }
          }

          //Includes a shift
          InterpolateGridValues(grid_trace,
                                z0_grid,     // Centre of first cell
                                dz_grid,
                                rAmpFine,
                                z0_data,     // Time of first data sample
                                dz_min,
                                rmt,
                                nz,
                                nzp);

          //Interpolate and shift trend before adding to grid_trace.
          //Alternative: add trend before interpolating and change values under l2 < 0 || l1 > n_fine
          if (grid_type != DATA) {
            std::vector<float> trend_interpolated(nzp);
            InterpolateAndShiftTrend(trend_interpolated,
                                     z0_grid,     // Centre of first cell
                                     dz_grid,
                                     data_trace_trend_long,
                                     z0_data,     // Time of first data sample
                                     dz_min,
                                     rmt,
                                     nz,
                                     nzp);
                                     //grid.GetNK());

            //Add trend
            for (size_t k_trace = 0; k_trace < grid_trace.size(); k_trace++)
              grid_trace[k_trace] += trend_interpolated[k_trace];
          }

          //err_text += err_text_tmp;

          fftw_free(rAmpData);
          fftw_free(rAmpFine);

          if (is_nrlib_grid)
            SetTrace(grid_trace, grid_new, i, j);
          else
            SetTrace(grid_trace, fft_grid_new, i, j);
        }
        else {
          if (is_nrlib_grid)
            SetTrace(0.0f, grid_new, i, j); // Dead traces (in case we allow them)
          else
            SetTrace(0.0f, fft_grid_new, i, j);

          dead_traces_simbox++;
        }
      }
      else {
        if (is_nrlib_grid)
          SetTrace(0.0f, grid_new, i, j);   // Outside seismic data grid
        else
          SetTrace(0.0f, fft_grid_new, i, j);

        if (i < nx && j < ny)
          missing_traces_simbox++;
        else
          missing_traces_padding++; //Won't happen with NRLib::Grid
      }

      if (rnxp*j + i + 1 >= static_cast<int>(nextMonitor)) {
        nextMonitor += monitorSize;
        printf("^");
        fflush(stdout);
      }

    }
  }
  LogKit::LogFormatted(LogKit::Low,"\n");

  fftwnd_destroy_plan(fftplan1);
  fftwnd_destroy_plan(fftplan2);

  Timings::setTimeResamplingSeismic(wall,cpu);
}

int CommonData::GetFillNumber(int i, int n, int np) {

  //  for the series                 i = 0,1,2,3,4,5,6,7
  //  GetFillNumber(i, 5 , 8)  returns   0,1,2,3,4,4,1,0 (cut middle, i.e 3,2)
  //  GetFillNumber(i, 4 , 8)  returns   0,1,2,3,3,2,1,0 (copy)
  //  GetFillNumber(i, 3 , 8)  returns   0,1,2,2,1,1,1,0 (drag middle out, i.e. 1)

  int refi     =  0;
  int BeloWnp, AbovEn;

  if (i< np)
  {
    if (i<n)
      // then it is in the main cube
      refi  =  i;
    else
    {
      // Get cyclic extention
      BeloWnp  = np-i-1;
      AbovEn   = i-n+1;
      if (AbovEn < BeloWnp)
      {
        // Then the index is closer to end than start.
        refi=std::max(n-AbovEn,n/2);
      } else {
        // The it is closer to  start than the end
        refi=std::min(BeloWnp,n/2);
      }//endif
    }//endif
  }//endif
  else
  {
    // This happens when the index is larger than the padding size
    // this happens in some cases because rnxp_ is larger than nxp_
    // and the x cycle is of length rnxp_
    refi=IMISSING;
  }//endif
  return(refi);
}

//void CommonData::FillInData(NRLib::Grid<double> & grid,
//                            const Simbox        & simbox,
//                            StormContGrid       * storm_grid,
//                            const SegY          * segy,
//                            float                 smooth_length,
//                            int                 & missing_traces_simbox,
//                            int                 & missing_traces_padding,
//                            int                 & dead_traces_simbox,
//                            std::string         & err_text,
//                            int                   grid_type,
//                            bool                  scale,
//                            bool                  is_segy)
//{
//  assert(grid_type != CTMISSING);
//
//  double wall=0.0, cpu=0.0;
//  TimeKit::getTime(wall,cpu);
//
//  float scalevert, scalehor;
//  if (scale == false || is_segy == true)
//  {
//    scalevert = 1.0;
//    scalehor  = 1.0;
//  }
//  else //from sgri file
//  {
//    LogKit::LogFormatted(LogKit::Low,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
//    scalevert = 0.001f; //1000.0;
//    scalehor  = 0.001f; //1000.0;
//  }
//
//  int ni = grid.GetNI();
//  int nj = grid.GetNJ();
//  int nk = grid.GetNK();
//
//  LogKit::LogFormatted(LogKit::Low,"\nResampling seismic data into %dx%dx%d grid:",ni,nj,nk);
//
//  float monitorSize = std::max(1.0f, static_cast<float>(nj*ni)*0.02f); //nyp_*rnxp_
//  float nextMonitor = monitorSize;
//  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
//  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
//  printf("\n  ^");
//
//  //
//  // Find proper length of time samples to get N*log(N) performance in FFT.
//  //
//  size_t n_samples = 0;
//  float  dz_data   = 0.0f;
//  float  dz_min    = 0.0f;
//  if (is_segy) {
//    n_samples = segy->FindNumberOfSamplesInLongestTrace();
//    dz_data   = segy->GetDz();
//    dz_min    = dz_data/4.0f;
//  }
//  else
//    n_samples = storm_grid->GetNK();
//
//  int nt = FindClosestFactorableNumber(static_cast<int>(n_samples));
//  int mt = 4*nt;           // Use four times the sampling density for the fine-meshed data
//
//  //
//  // Create FFT plans
//  //
//  rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
//  rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
//
//  //
//  // Do resampling
//  //
//  missing_traces_simbox  = 0; // Part of simbox is outside seismic data
//  missing_traces_padding = 0; // Part of padding is outside seismic data
//  dead_traces_simbox     = 0; // Simbox is inside seismic data but trace is missing
//
//  for (int j = 0 ; j < nj ; j++) {
//    for (int i = 0 ; i < ni ; i++) {
//
//      int refi  = i; //getFillNumber(i, nx_, nxp_ ); // Find index (special treatment for padding)
//      int refj  = j; //getFillNumber(j, ny_, nyp_ ); // Find index (special treatment for padding)
//      int refk  = 0;
//
//      double x, y, z0;
//      simbox.getCoord(refi, refj, refk, x, y, z0);  // Get lateral position and z-start (z0)
//      x*=  scalehor;
//      y*=  scalehor;
//      z0*= scalevert;
//
//      double dz = simbox.getdz(refi, refj)*scalevert;
//      float  xf = static_cast<float>(x);
//      float  yf = static_cast<float>(y);
//
//      bool is_inside = false;
//      if (is_segy)
//        is_inside = segy->GetGeometry()->IsInside(xf, yf);
//      else {
//        if (storm_grid->IsInside(xf, yf) == 1)
//          is_inside = true;
//      }
//
//      if (is_inside == true) {
//        bool  missing = true;
//        float z0_data = RMISSING;
//
//        std::vector<float> data_trace;
//        size_t grid_i = 0;
//        size_t grid_j = 0;
//        double grid_x = 0.0;
//        double grid_y = 0.0;
//        double grid_z = 0.0;
//        float value =   0.0f;
//
//        float z_min = 0.0f;
//        float z_max = 0.0f;
//
//        //Get data_trace for this i and j.
//        if (is_segy) {
//          segy->GetNearestTrace(data_trace, missing, z0_data, xf, yf);
//        }
//        else {
//          storm_grid->FindXYIndex(xf, yf, grid_i, grid_j);
//          for (size_t k = 0; k < grid.GetNK(); k++) {
//            storm_grid->FindCenterOfCell(grid_i, grid_j, k, grid_x, grid_y, grid_z);
//            value = storm_grid->GetValueZInterpolated(grid_x, grid_y, grid_z);
//            data_trace.push_back(value);
//
//            if (k == 0)
//              z_min = static_cast<float>(grid_z);
//            if (k == grid.GetNK()-1)
//              z_max = static_cast<float>(grid_z);
//          }
//
//          dz_data = (z_max- z_min) / grid.GetNK();
//          dz_min = dz_data/4.0f;
//          z0_data = z_min;
//        }
//
//        size_t n_trace = data_trace.size();
//        float trend_first = 0.0f;
//        float trend_last  = 0.0f;
//
//        if (grid_type != DATA) {
//          //Remove zeroes. F.ex. background on segy-format with a non-constant top-surface, the vector is filled with zeroes at the beginning.
//          if (data_trace[0] == 0) {
//            std::vector<float> data_trace_new;
//            for (size_t k_trace = 0; k_trace < n_trace; k_trace++) {
//              if (data_trace[k_trace] != 0)
//                data_trace_new.push_back(data_trace[k_trace]);
//            }
//            data_trace = data_trace_new;
//            n_trace = data_trace.size();
//          }
//
//          n_samples = data_trace.size();
//          nt = FindClosestFactorableNumber(static_cast<int>(n_samples));
//          mt = 4*nt;
//
//          fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
//          fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
//
//          //Remove trend from trace
//          trend_first = data_trace[0];
//          trend_last = data_trace[n_trace - 1];
//          float trend_inc = (trend_last - trend_first) / (n_trace - 1);
//          for (size_t k_trace = 0; k_trace < data_trace.size(); k_trace++) {
//            data_trace[k_trace] -= trend_first + k_trace * trend_inc;
//          }
//        }
//
//        if (is_segy == false || (is_segy == true && !missing)) {
//          int         cnt      = nt/2 + 1;
//          int         rnt      = 2*cnt;
//          int         cmt      = mt/2 + 1;
//          int         rmt      = 2*cmt;
//
//          fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
//          fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));
//
//          float       dz_grid  = static_cast<float>(dz);
//          float       z0_grid  = static_cast<float>(z0);
//
//          std::vector<float> grid_trace(nk);
//
//          std::string err_text_tmp = "";
//          smooth_length*=scalevert;
//          if (grid_type == DATA) {
//            SmoothTraceInGuardZone(data_trace,
//                                   dz_data,
//                                   smooth_length);
//          }
//
//          ResampleTrace(data_trace,
//                        fftplan1,
//                        fftplan2,
//                        rAmpData,
//                        rAmpFine,
//                        cnt,
//                        rnt,
//                        cmt,
//                        rmt);
//
//          std::vector<float> data_trace_trend_long;
//          if (grid_type != DATA) {
//            float trend_inc = (trend_last - trend_first) / (rmt - 1);
//
//            data_trace_trend_long.resize(rmt);
//            for (int k_trace = 0; k_trace < rmt; k_trace++) {
//              data_trace_trend_long[k_trace] = trend_first + k_trace * trend_inc;
//            }
//          }
//
//          //Includes a shift
//          InterpolateGridValues(grid_trace,
//                                z0_grid,     // Centre of first cell
//                                dz_grid,
//                                rAmpFine,
//                                z0_data,     // Time of first data sample
//                                dz_min,
//                                rmt,
//                                grid.GetNK());
//
//          //Interpolate and shift trend before adding to grid_trace.
//          //Alternative: add trend before interpolating and change values under l2 < 0 || l1 > n_fine
//          if (grid_type != DATA) {
//            std::vector<float> trend_interpolated(nk);
//            InterpolateAndShiftTrend(trend_interpolated,
//                                     z0_grid,     // Centre of first cell
//                                     dz_grid,
//                                     data_trace_trend_long,
//                                     z0_data,     // Time of first data sample
//                                     dz_min,
//                                     rmt,
//                                     grid.GetNK());
//
//            //Add trend
//            for (size_t k_trace = 0; k_trace < grid_trace.size(); k_trace++)
//              grid_trace[k_trace] += trend_interpolated[k_trace];
//          }
//
//          err_text += err_text_tmp;
//
//          fftw_free(rAmpData);
//          fftw_free(rAmpFine);
//
//          SetTrace(grid_trace, grid, i, j);
//        }
//        else {
//          SetTrace(0.0f, grid, i, j); // Dead traces (in case we allow them)
//          dead_traces_simbox++;
//        }
//      }
//      else {
//        SetTrace(0.0f, grid, i, j);   // Outside seismic data grid
//        if (i < ni && j < nj )
//          missing_traces_simbox++;
//        else
//          missing_traces_padding++; //H Won't happen with NRLibGrid
//      }
//
//      if (ni*j + i + 1 >= static_cast<int>(nextMonitor)) { //rnxp_
//        nextMonitor += monitorSize;
//        printf("^");
//        fflush(stdout);
//      }
//
//    }
//  }
//  LogKit::LogFormatted(LogKit::Low,"\n");
//
//  fftwnd_destroy_plan(fftplan1);
//  fftwnd_destroy_plan(fftplan2);
//
//  Timings::setTimeResamplingSeismic(wall,cpu);
//}

int CommonData::FindClosestFactorableNumber(int leastint)
{
  int i,j,k,l,m,n;
  int factor   =       1;

  int maxant2    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(2.0f) ));
  int maxant3    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(3.0f) ));
  int maxant5    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(5.0f) ));
  int maxant7    = static_cast<int>(ceil(static_cast<double>(log(static_cast<float>(leastint))) / log(7.0f) ));
  int maxant11   = 0;
  int maxant13   = 0;
  int closestprod= static_cast<int>(pow(2.0f,maxant2));

  /* kan forbedres ved aa trekke fra i endepunktene.i for lokkene*/
  for (i=0;i<maxant2+1;i++)
    for (j=0;j<maxant3+1;j++)
      for (k=0;k<maxant5+1;k++)
        for (l=0;l<maxant7+1;l++)
          for (m=0;m<maxant11+1;m++)
            for (n=maxant11;n<maxant13+1;n++)
            {
              factor = static_cast<int>(pow(2.0f,i)*pow(3.0f,j)*pow(5.0f,k)*
                pow(7.0f,l)*pow(11.0f,m)*pow(13.0f,n));
              if ((factor >=  leastint) &&  (factor <  closestprod))
              {
                closestprod=factor;
              }
            }
            return closestprod;
}

void CommonData::SmoothTraceInGuardZone(std::vector<float> & data_trace,
                                        float                dz_data,
                                        float                smooth_length)
{
  // We recommend a guard zone of at least half a wavelet on each side of
  // the target zone and that half a wavelet of the guard zone is smoothed.
  //
  // By default, we have: guard_zone = smooth_length = 0.5*wavelet = 100ms
  //
  // Originally, we checked here that the guard zone was larger than or equal
  // to the smooth length. However, as a trace is generally not located in
  // the center of the grid cell, we made an invalid comparison of z-values
  // from different reference systems. Only when the top/base surface is flat
  // such a comparison becomes valid. Instead, we hope and assume that the
  // SegY volume made in ModelGeneral::readSegyFile() has been made with
  // the correct guard zone in method.

 //
  // k=n_smooth is the first sample within the simbox. For this sample
  // the smoothing factor (if it had been applied) should be one.
  //
  int n_smooth = static_cast<int>(floor(smooth_length/dz_data));

  for (int k = 0 ; k < n_smooth ; k++) {
    double theta   = static_cast<double>(k)/static_cast<double>(n_smooth);
    float  sinT    = static_cast<float>(std::sin(NRLib::PiHalf*theta));
    data_trace[k] *= sinT*sinT;
  }

  int n_data = data_trace.size();
  int kstart = n_data - n_smooth;

  for (int k = 0 ; k < n_smooth ; k++) {
    double theta   = static_cast<double>(k + 1)/static_cast<double>(n_smooth);
    float  cosT    = static_cast<float>(std::cos(NRLib::PiHalf*theta));
    data_trace[kstart + k] *= cosT*cosT;
  }
}

void CommonData::ResampleTrace(const std::vector<float> & data_trace,
                               const rfftwnd_plan       & fftplan1,
                               const rfftwnd_plan       & fftplan2,
                               fftw_real                * rAmpData,
                               fftw_real                * rAmpFine,
                               int                        cnt,
                               int                        rnt,
                               int                        cmt,
                               int                        rmt)
{
  fftw_complex * cAmpData = reinterpret_cast<fftw_complex*>(rAmpData);
  fftw_complex * cAmpFine = reinterpret_cast<fftw_complex*>(rAmpFine);

  //
  // Fill vector to be FFT'ed
  //
  int n_data = static_cast<int>(data_trace.size());

  for (int i = 0 ; i < n_data ; i++) {
    rAmpData[i] = data_trace[i];
  }
  // Pad with zeros
  for (int i = n_data ; i < rnt ; i++) {
    rAmpData[i] = 0.0f;
  }

  //
  // Transform to Fourier domain
  //
  rfftwnd_one_real_to_complex(fftplan1, rAmpData, cAmpData);

  //
  // Fill fine-sampled grid
  //
  for (int i = 0 ; i < cnt ; i++) {
    cAmpFine[i].re = cAmpData[i].re;
    cAmpFine[i].im = cAmpData[i].im;
  }
  // Pad with zeros (cmt is always greater than cnt)
  for (int i = cnt ; i < cmt ; i++) {
    cAmpFine[i].re = 0.0f;
    cAmpFine[i].im = 0.0f;
  }

  //
  // Fine-sampled grid: Fourier --> Time
  //
  rfftwnd_one_complex_to_real(fftplan2, cAmpFine, rAmpFine);

  //
  // Scale and fill grid_trace
  //
  float scale = 1/static_cast<float>(rnt);
  for (int i = 0 ; i < rmt ; i++) {
    rAmpFine[i] = scale*rAmpFine[i];
  }
}

void CommonData::InterpolateGridValues(std::vector<float> & grid_trace,
                                       float                z0_grid,
                                       float                dz_grid,
                                       fftw_real          * rAmpFine,
                                       float                z0_data,
                                       float                dz_fine,
                                       int                  n_fine,
                                       int                  nz,
                                       int                  nzp)
{
  //
  // Bilinear interpolation
  //
  // refk establishes link between traces order and grid order
  // In trace:    A A A B B B B B B C C C     (A and C are values in padding)
  // In grid :    B B B B B B C C C A A A

  float z0_shift    = z0_grid - z0_data;
  float inv_dz_fine = 1.0f/dz_fine;

  int n_grid = static_cast<int>(grid_trace.size());

  for (int k = 0 ; k < n_grid ; k++) {
    int refk = GetZSimboxIndex(k, nz, nzp);
    float dl = (z0_shift + static_cast<float>(refk)*dz_grid)*inv_dz_fine;
    int   l1 = static_cast<int>(floor(dl));
    int   l2 = static_cast<int>(ceil(dl));

    if (l2 < 0 || l1 > n_fine - 1) {
      grid_trace[k] = 0.0f;
    }
    else {
      if (l1 < 0) {
        grid_trace[k] = rAmpFine[l2];
      }
      else if (l2 > n_fine - 1) {
        grid_trace[k] = rAmpFine[l1];
      }
      else if (l1 == l2) {
        grid_trace[k] = rAmpFine[l1];
      }
      else {
        float w1 = ceil(dl) - dl;
        float w2 = dl - floor(dl);
        grid_trace[k] = w1*rAmpFine[l1] + w2*rAmpFine[l2];
      }
    }
  }
}

void CommonData::InterpolateAndShiftTrend(std::vector<float>       & interpolated_trend,
                                          float                      z0_grid,
                                          float                      dz_grid,
                                          const std::vector<float> & trend_long,
                                          float                      z0_data,
                                          float                      dz_fine,
                                          int                        n_fine,
                                          int                        nz,
                                          int                        nzp)
{
  //
  // Bilinear interpolation
  //
  // refk establishes link between traces order and grid order
  // In trace:    A A A B B B B B B C C C     (A and C are values in padding)
  // In grid :    B B B B B B C C C A A A

  float z0_shift    = z0_grid - z0_data;
  float inv_dz_fine = 1.0f/dz_fine;

  int n_grid = static_cast<int>(interpolated_trend.size());

  for (int k = 0 ; k < n_grid ; k++) {
    int refk = GetZSimboxIndex(k, nz, nzp);
    float dl = (z0_shift + static_cast<float>(refk)*dz_grid)*inv_dz_fine;
    int   l1 = static_cast<int>(floor(dl));
    int   l2 = static_cast<int>(ceil(dl));

    if (l2 < 0 || l1 > n_fine - 1) {
      if (l2 < 0)
        interpolated_trend[k] = trend_long[0];
      if (l1 > n_fine -1)
        interpolated_trend[k] = trend_long[n_fine - 1];
    }
    else {
      if (l1 < 0) {
        interpolated_trend[k] = trend_long[l2];
      }
      else if (l2 > n_fine - 1) {
        interpolated_trend[k] = trend_long[l1];
      }
      else if (l1 == l2) {
        interpolated_trend[k] = trend_long[l1];
      }
      else {
        float w1 = ceil(dl) - dl;
        float w2 = dl - floor(dl);
        interpolated_trend[k] = w1*trend_long[l1] + w2*trend_long[l2];
      }
    }
  }
}

int CommonData::GetZSimboxIndex(int k,
                                int nz,
                                int nzp)
{
  int refk;

  if (k < (nz + nzp)/2)
    refk=k;
  else
    refk=k-nzp;

  return refk;
}

//void CommonData::SetTrace(const std::vector<float> & trace,
//                          NRLib::Grid<double>      & grid,
//                          size_t                     i,
//                          size_t                     j)
//{
//  for (size_t k = 0; k < grid.GetNK(); k++) {
//    grid(i, j, k) = static_cast<double>(trace[k]);
//  }
//}

void CommonData::SetTrace(const std::vector<float> & trace,
                          NRLib::Grid<double>      * grid,
                          size_t                     i,
                          size_t                     j)
{
  for (size_t k = 0; k < grid->GetNK(); k++) {
    double value_tmp = static_cast<double>(trace[k]);
    grid->SetValue(i, j, k, value_tmp);
  }
}

//void CommonData::SetTrace(float                 value,
//                          NRLib::Grid<double> & grid,
//                          size_t                i,
//                          size_t                j)
//{
//  for (size_t k = 0; k < grid.GetNK(); k++) {
//    grid(i, j, k) = static_cast<double>(value);
//  }
//}

void CommonData::SetTrace(float                 value,
                          NRLib::Grid<double> * grid,
                          size_t                i,
                          size_t                j)
{
  for (size_t k = 0; k < grid->GetNK(); k++) {
    double value_tmp = static_cast<double>(value);
    grid->SetValue(i, j, k, value_tmp);
  }
}

void CommonData::SetTrace(const std::vector<float> & trace,
                          FFTGrid                  * grid,
                          size_t                     i,
                          size_t                     j)
{
  for (int k = 0; k < grid->getNzp(); k++) {
    grid->setRealValue(i, j, k, trace[k], true);
  }
}

void CommonData::SetTrace(float     value,
                          FFTGrid * grid,
                          size_t    i,
                          size_t    j)
{
  for (int k = 0; k < grid->getNzp(); k++) {
    grid->setRealValue(i, j, k, value, true);
  }
}

void
CommonData::ReadStormFile(const std::string                  & file_name,
                          std::vector<NRLib::Grid<double> *> & interval_grids,
                          const int                            grid_type,
                          const std::string                  & par_name,
                          const std::vector<Simbox>          & interval_simboxes,
                          const ModelSettings                * model_settings,
                          std::string                        & err_text,
                          bool                                 scale,
                          bool                                 nopadding) {
 (void) nopadding;
 (void) par_name;
  StormContGrid * stormgrid = NULL;
  bool failed = false;

  try {
    stormgrid = new StormContGrid(0,0,0);
    stormgrid->ReadFromFile(file_name);
  }
  catch (NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  int missing_traces_simbox  = 0;
  int missing_traces_padding = 0;
  int dead_traces_simbox     = 0;

  if (failed == false) {

    for (size_t i_interval = 0; i_interval < interval_simboxes.size(); i_interval++) {

      int nx = interval_simboxes[i_interval].getnx();
      int ny = interval_simboxes[i_interval].getny();
      int nz = interval_simboxes[i_interval].getnz();

      interval_grids[i_interval]->Resize(nx, ny, nz, 0.0);

      try {

        SegY * segy_tmp          = NULL;
        FFTGrid * fft_grid_tmp   = NULL;
        FFTGrid * fft_grid_tmp_2 = NULL;
        FillInData(interval_grids[i_interval],
                   fft_grid_tmp,
                   interval_simboxes[i_interval],
                   stormgrid,
                   segy_tmp,
                   fft_grid_tmp_2,
                   model_settings->getSmoothLength(),
                   missing_traces_simbox,
                   missing_traces_padding,
                   dead_traces_simbox,
                   grid_type,
                   scale,
                   false, //is_segy
                   true); //is_storm

        if (segy_tmp != NULL)
         delete segy_tmp;
        if (fft_grid_tmp != NULL)
          delete fft_grid_tmp;
        if (fft_grid_tmp_2 != NULL)
          delete fft_grid_tmp_2;

      }
      catch (NRLib::Exception & e) {
        err_text += std::string(e.what());
      }

      if (missing_traces_simbox > 0) { //outsideTraces
        if (missing_traces_simbox == interval_simboxes[i_interval].getnx()*interval_simboxes[i_interval].getny()) {
          err_text += "Error: Data in file \'"+file_name+"\' was completely outside the inversion area.\n";
          failed = true;
        }
        else {
          if (grid_type == FFTGrid::PARAMETER) {
            err_text += "Error: Data read from file \'"+file_name+"\' does not cover the inversion area.\n";
          }
          else {
            LogKit::LogMessage(LogKit::Warning, "WARNING: "+NRLib::ToString(missing_traces_simbox)
                               + " grid columns were outside the seismic data in file \'"+file_name+"\'.\n");
            TaskList::addTask("Check seismic data and inversion area: One or volumes did not have data enough to cover entire grid.\n");
         }
        }
      }
      //if (missing_traces_padding > 0) {
      //  int nx     = interval_simboxes[i_interval].getnx();
      //  int ny     = interval_simboxes[i_interval].getny();
      //  int nxpad  = xpad - nx;
      //  int nypad  = ypad - ny;
      //  int nxypad = nxpad*ny + nx*nypad - nxpad*nypad;
      //  LogKit::LogMessage(LogKit::High, "Number of grid columns in padding that are outside area defined by seismic data : "
      //                     +NRLib::ToString(missing_traces_padding)+" of "+NRLib::ToString(nxypad)+"\n");
      //}
      if (dead_traces_simbox > 0) {
        LogKit::LogMessage(LogKit::High, "Number of grid columns with no seismic data (nearest trace is dead) : "
                           +NRLib::ToString(dead_traces_simbox)+" of "+NRLib::ToString(interval_simboxes[i_interval].getnx()*interval_simboxes[i_interval].getny())+"\n");
      }

      //if (outsideTraces > 0) {
      //  if (outsideTraces == interval_simboxes[i_interval].getnx()*interval_simboxes[i_interval].getny()) {
      //    err_text += "Error: Data in file \'"+file_name+"\' was completely outside the inversion area.\n";
      //    failed = true;
      //  }
      //  else {
      //    if (grid_type == PARAMETER) {
      //      err_text += "Error: Data read from file \'"+file_name+"\' does not cover the inversion area.\n";
      //    }
      //    else {
      //      LogKit::LogMessage(LogKit::Warning, "WARNING: "+NRLib::ToString(outsideTraces)
      //                         + " grid columns were outside the seismic data in file \'"+file_name+"\'.\n");
      //      TaskList::addTask("Check seismic data and inversion area: One or volumes did not have data enough to cover entire grid.\n");
      //    }
      //  }
      //}
    } //i_interval

    if (stormgrid != NULL)
      delete stormgrid;

  }
}

bool CommonData::SetupDepthConversion(ModelSettings  * model_settings,
                                      InputFiles     * input_files,
                                      std::string    & err_text_common) {


  //From ModelGeneral::ProcessDepthConversion
  std::string err_text = "";

  NRLib::Grid<double> * velocity = NULL;
  std::string velocity_field     = input_files->getVelocityField();

  LoadVelocity(velocity,
               &estimation_simbox_,
               model_settings,
               velocity_field,
               velocity_from_inversion_,
               err_text);

  if (err_text == "") {
    bool failed_dummy = false;

    time_depth_mapping_ = new GridMapping();
    time_depth_mapping_->setDepthSurfaces(input_files->getDepthSurfFiles(), failed_dummy, err_text);

    if (velocity != NULL) {
      time_depth_mapping_->calculateSurfaceFromVelocity(velocity, &estimation_simbox_);
      time_depth_mapping_->setDepthSimbox(&estimation_simbox_, estimation_simbox_.getnz(),
                                        model_settings->getOutputGridFormat(),
                                        failed_dummy, err_text);            // NBNB-PAL: Er dettet riktig nz (timeCut vs time)?
      time_depth_mapping_->makeTimeDepthMapping(velocity, &estimation_simbox_);

      if ((model_settings->getOutputGridsOther() & IO::TIME_TO_DEPTH_VELOCITY) > 0) { //H Currently create a temporary fft_grid and use the writing in fftgrid.cpp. Add writing-functions to NRLib::Grid?
        std::string base_name  = IO::FileTimeToDepthVelocity();
        std::string sgri_label = std::string("Time-to-depth velocity");
        float       offset    = model_settings->getSegyOffset(0);//Only allow one segy offset for time lapse data

        FFTGrid * velocity_fft = new FFTGrid(velocity, velocity->GetNI(), velocity->GetNJ(), velocity->GetNK());
        velocity_fft->writeFile(base_name,
                                IO::PathToVelocity(),
                                &estimation_simbox_,
                                sgri_label,
                                offset,
                                time_depth_mapping_);
        if (velocity_fft != NULL)
          delete velocity_fft;
      }
    }
    else if (velocity == NULL && velocity_from_inversion_ ==false) {
      time_depth_mapping_->setDepthSimbox(&estimation_simbox_,
                                          estimation_simbox_.getnz(),
                                          model_settings->getOutputGridFormat(),
                                          failed_dummy,
                                          err_text);
    }
  }

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}


bool CommonData::SetupBackgroundModel(ModelSettings  * model_settings,
                                      InputFiles     * input_files,
                                      std::string    & err_text_common) {

  std::string err_text = "";

  if (forward_modeling_)
    LogKit::WriteHeader("Earth Model");
  else
    LogKit::WriteHeader("Prior Expectations / Background Model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  //Skip if background is specified as false in estimation mode
  if(model_settings->getEstimationMode()     == true  && model_settings->getEstimateWaveletNoise() == true &&
     model_settings->getEstimateBackground() == false && model_settings->getEstimateCorrelations() == false)
     return true;

  if (model_settings->getGenerateBackground()) {

    if (model_settings->getGenerateBackgroundFromRockPhysics() == false) {

      NRLib::Grid<double> * velocity = new NRLib::Grid<double>();
      std::string back_vel_file = input_files->getBackVelFile();
      if (back_vel_file != "") {
        bool dummy;
        LoadVelocity(velocity,
                     &estimation_simbox_,
                     model_settings,
                     back_vel_file,
                     dummy,
                     err_text);
      }
      if (err_text == "") {

        if (model_settings->getBackgroundVario() == NULL) {
          err_text += "There is no variogram available for the background modelling.\n";
        }

        Surface           * correlation_direction = NULL;
        Simbox            * bg_simbox             = NULL;
        BlockedLogsCommon * bl_bg                 = NULL;
        std::map<std::string, BlockedLogsCommon *> mapped_bg_bl;

        if (input_files->getCorrDirFile() != "") {

          Surface tmpSurf(input_files->getCorrDirFile());
          if (estimation_simbox_.CheckSurface(tmpSurf) == true)
            correlation_direction = new Surface(tmpSurf);
          else {
            err_text += "Error: Correlation surface does not cover volume.\n";
          }

          //Create a sepeate background simbox based on correlation_direction
          //H In ModelGeneral setupExtendedTimeSimbox and setupExtendedBackgroundSimbox were created seperatly
          //  Could extended_background_simbox we covered/included by the new simbox-format?
          SetupExtendedBackgroundSimbox(&estimation_simbox_, correlation_direction, bg_simbox, model_settings->getOutputGridFormat(), model_settings->getOutputGridDomain(), model_settings->getOtherOutputFlag());
          std::string err_text_tmp = "";
          int status = bg_simbox->calculateDz(model_settings->getLzLimit(), err_text_tmp);
          if(status != Simbox::BOXOK) {
            err_text += "Could not make the grid for background model.\n";
            err_text += err_text_tmp;
          }

        }
        else if (input_files->getCorrDirBaseFile() != "" || model_settings->getCorrDirBaseConform()) {
          //H Create a backgroundsimbox based on the different options for correlation direction? top/base conform.
          //Similiar to MultiIntervalGrid::SetupIntervalSimbox

          //SetupExtendedBackgroundSimbox( two correlations surfaces);

        }

        for (size_t i = 0; i < wells_.size(); i++) {
          bl_bg = NULL;

          // Get all continuous and discrete logs
          std::vector<std::string> cont_logs_to_be_blocked;
          std::vector<std::string> disc_logs_to_be_blocked;

          const std::map<std::string,std::vector<double> > & cont_logs = wells_[i].GetContLog();
          const std::map<std::string,std::vector<int> >    & disc_logs = wells_[i].GetDiscLog();

          for (std::map<std::string,std::vector<double> >::const_iterator it = cont_logs.begin(); it!=cont_logs.end(); it++) {
            cont_logs_to_be_blocked.push_back(it->first);
          }
          for (std::map<std::string,std::vector<int> >::const_iterator it = disc_logs.begin(); it!=disc_logs.end(); it++) {
            disc_logs_to_be_blocked.push_back(it->first);
          }

          if (bg_simbox == NULL)
            bl_bg = new BlockedLogsCommon(&wells_[i],
                                          cont_logs_to_be_blocked,
                                          disc_logs_to_be_blocked,
                                          &estimation_simbox_,
                                          false,
                                          err_text);
          else
            bl_bg = new BlockedLogsCommon(&wells_[i],
                                          cont_logs_to_be_blocked,
                                          disc_logs_to_be_blocked,
                                          bg_simbox,
                                          false,
                                          err_text);

          mapped_bg_bl.insert(std::pair<std::string, BlockedLogsCommon *>(wells_[i].GetWellName(), bl_bg));

        }


        if (model_settings->getIntervalNames().size() == 0) {
          std::vector<std::vector<NRLib::Grid<double> *> > & parameters = multiple_interval_grid_->GetBackgroundParameters(); //interval vector(parameter)

          if (model_settings->getMultizoneBackground() == true)
            Background(parameters[0], wells_, &estimation_simbox_, model_settings, input_files->getMultizoneSurfaceFiles(), err_text); //Multizone background model (not multiple intervals)
          else //Neither multizone or multiinterval
            Background(parameters[0], wells_, velocity, &estimation_simbox_, bg_simbox, mapped_blocked_logs_, mapped_bg_bl, model_settings, err_text);

          double vs_vp_ratio = FindMeanVsVp(parameters[0][0], parameters[0][1]);

          for (int j = 0; j < 3; j++)
            multiple_interval_grid_->AddBackgroundParameterForInterval(0, j, parameters[0][j]);

          multiple_interval_grid_->SetBackgroundVsVpRatio(0, vs_vp_ratio);
        }
        else { //Multiple intervals, not multizone background
          std::vector<std::vector<NRLib::Grid<double> *> > & parameters = multiple_interval_grid_->GetBackgroundParameters(); //vector(intervals) vector(parameters)
          std::vector<double> vs_vp_ratios(model_settings->getIntervalNames().size());

          Background(parameters, wells_, multiple_interval_grid_, model_settings, err_text);

          for (size_t i = 0; i < model_settings->getIntervalNames().size(); i++) {
            //for (int j = 0; j < 3; j++)
            //  multiple_interval_grid_->AddBackgroundParameterForInterval(i, j, &parameters[j][i]);

            vs_vp_ratios[i] = FindMeanVsVp(parameters[i][0], parameters[i][1]);
          }

          multiple_interval_grid_->SetBackgroundVsVpRatios(vs_vp_ratios);
        }
      }
    }
    else {

      std::vector<std::vector<NRLib::Grid<double> *> > & parameters = multiple_interval_grid_->GetBackgroundParameters();

      for (int i_interval = 0; i_interval < multiple_interval_grid_->GetNIntervals(); i_interval++) {

        // Get prior probabilities for the facies in a vector
        std::vector<std::string> facies_names = facies_names_;
        int                      n_facies     = static_cast<int>(facies_names.size());

        std::vector<float> prior_probability = prior_facies_[i_interval];

        std::vector<DistributionsRock *> rock_distribution(n_facies);  //Get from rock_distributions_

        typedef std::map<std::string, DistributionsRock *> rf_map_type;
        rf_map_type rf_map = GetRockDistributionTime0();

        for (int j = 0; j < n_facies; j++) {
          rf_map_type::iterator iter = rf_map.find(facies_names[j]);
          if (iter != rf_map.end())
            rock_distribution[j] = iter->second;
        }

        // filling in the backModel in Background
        GenerateRockPhysics3DBackground(rock_distribution,
                                        prior_probability,
                                        parameters[i_interval][0],
                                        parameters[i_interval][1],
                                        parameters[i_interval][2],
                                        estimation_simbox_,
                                        multiple_interval_grid_->GetTrendCube(i_interval));
        double vs_vp_ratio = FindMeanVsVp(parameters[i_interval][0], parameters[i_interval][1]);

        //for (int j = 0; j < 3; j++)
        //  multiple_interval_grid_->AddBackgroundParameterForInterval(i_interval, j, parameters[j]);

        multiple_interval_grid_->SetBackgroundVsVpRatio(i_interval, vs_vp_ratio);

      }
    }
  }
  else {

    int n_intervals = multiple_interval_grid_->GetNIntervals();

    std::vector<std::vector<NRLib::Grid<double> *> > parameters; //vector(parameters) vector(intervals). Parameters on the outside due to reading grid from file.
    parameters.resize(3);
    for (int i = 0; i < 3; i++) {
      parameters[i].resize(n_intervals);
      for (int j = 0; j < n_intervals; j++) {
        parameters[i][j] = new NRLib::Grid<double>();
      }
    }

    std::vector<std::string> par_name;
    if (model_settings->getUseAIBackground())
      par_name.push_back("AI "+model_settings->getBackgroundType());
    else
      par_name.push_back("Vp "+model_settings->getBackgroundType());
    if (model_settings->getUseSIBackground())
      par_name.push_back("SI "+model_settings->getBackgroundType());
    else if (model_settings->getUseVpVsBackground())
      par_name.push_back("Vp/Vs "+model_settings->getBackgroundType());
    else
      par_name.push_back("Vs "+model_settings->getBackgroundType());
    par_name.push_back("Rho "+model_settings->getBackgroundType());

    for (int i = 0; i < 3; i++) {
      float const_back_value        = model_settings->getConstBackValue(i);
      const std::string & back_file = input_files->getBackFile(i);
      //parameters_tmp[i].resize(n_intervals);

      if (const_back_value < 0) {
        if (back_file.size() > 0) {
          const SegyGeometry      * dummy1 = NULL;
          const TraceHeaderFormat * dummy2 = NULL;
          const float               offset = model_settings->getSegyOffset(0); //H Currently set to 0. In ModelAVODynamic getSegyOffset(thisTimeLapse) was used. Create loop over timelapses?
          std::string err_text_tmp = "";

          ReadGridFromFile(back_file,
                           par_name[i],
                           offset,
                           parameters[i], //This is split in intervals inside.
                           dummy1,
                           dummy2,
                           PARAMETER,
                           multiple_interval_grid_->GetIntervalSimboxes(),
                           &estimation_simbox_,
                           model_settings,
                           err_text_tmp);

          if (err_text_tmp != "") {
            err_text += err_text_tmp;
            err_text += "Reading of file '"+back_file+"' for parameter '"+par_name[i]+"' failed\n\n";
          }
          else {
            for (int j = 0; j < n_intervals; j++) {
              double min = 0.0;
              double max = 0.0;
              double avg = 0.0;
              parameters[i][j]->GetAvgMinMax(avg, min, max);
              SetUndefinedCellsToGlobalAverageGrid(parameters[i][j], avg);
              parameters[i][j]->LogTransform(RMISSING);
            }
          }
        }
        else {
          err_text += "Reading of file for parameter "+par_name[i]+" failed. No file name is given.\n";
        }
      }
      else if (const_back_value > 0) {

        for (int j = 0; j < n_intervals; j++) {

          int nx = 0;
          int ny = 0;
          int nz = 0;

          if (model_settings->getIntervalNames().size() > 0) {
            const Simbox * simbox = multiple_interval_grid_->GetIntervalSimbox(j);
            nx = simbox->getnx();
            ny = simbox->getny();
            nz = simbox->getnz();
          }
          else {
            nx = estimation_simbox_.getnx();
            ny = estimation_simbox_.getny();
            nz = estimation_simbox_.getnz();
          }

          double log_value = log(const_back_value);
          parameters[i][j]->Resize(nx, ny, nz, log_value);

          //back_model[i] = CreateFFTGrid(nx, ny, nz, nx_pad, ny_pad, nz_pad, model_settings->getFileGrid());
          //back_model[i]->setType(FFTGrid::PARAMETER);
          //back_model[i]->fillInConstant(float( log( const_back_value )));
          //back_model[i]->calculateStatistics();
        }
      }
      else {
        err_text += "Trying to set background model to 0 for parameter "+par_name[i]+"\n";
      }
    } //for i = 0,1,2

    if (err_text == "") {

      for (int j = 0; j < n_intervals; j++) {

        if (multiple_interval_grid_->GetIntervalNames().size() > 0)
          LogKit::LogFormatted(LogKit::Low, "\nInterval " + multiple_interval_grid_->GetIntervalName(j) + "\n");
        LogKit::LogFormatted(LogKit::Low, "\nSummary                Average   Minimum   Maximum\n");
        LogKit::LogFormatted(LogKit::Low, "--------------------------------------------------\n");

        for (int i = 0; i < 3; i++) {

          double avg = 0.0;
          double min = 0.0;
          double max = 0.0;
          parameters[i][j]->GetAvgMinMax(avg, min, max);

          LogKit::LogFormatted(LogKit::Low, "%-20s %9.2f %9.2f %9.2f\n",
                               par_name[i].c_str(),
                               avg,
                               min,
                               max);

        }
        if (model_settings->getUseAIBackground()) { // Vp = AI/Rho     ==> lnVp = lnAI - lnRho
          LogKit::LogMessage(LogKit::Low, "\nMaking Vp background from AI and Rho\n");
          SubtractGrid(parameters[0][j], parameters[2][j]);
          //back_model[0]->subtract(back_model[2]);
        }
        if (model_settings->getUseSIBackground()) { // Vs = SI/Rho     ==> lnVs = lnSI - lnRho
          LogKit::LogMessage(LogKit::Low, "\nMaking Vs background from SI and Rho\n");
          SubtractGrid(parameters[1][j], parameters[2][j]);
          //back_model[1]->subtract(back_model[2]);
        }
        else if (model_settings->getUseVpVsBackground()) { // Vs = Vp/(Vp/Vs) ==> lnVs = lnVp - ln(Vp/Vs)
          LogKit::LogMessage(LogKit::Low, "\nMaking Vs background from Vp and Vp/Vs\n");
          SubtractGrid(parameters[1][j], parameters[0][j]);
          ChangeSignGrid(parameters[1][j]);
          //back_model[1]->subtract(back_model[0]);
          //back_model[1]->changeSign();
        }
      }
    }

    //Add background model to multiintervalgrid
    std::vector<double> vs_vp_ratios(n_intervals);

    for (int i = 0; i < n_intervals; i++) {
     for (int j = 0; j < 3; j++) {
       multiple_interval_grid_->AddBackgroundParameterForInterval(i, j, parameters[j][i]);
     }

     vs_vp_ratios[i] = FindMeanVsVp(parameters[0][i], parameters[1][i]);
    }
    multiple_interval_grid_->SetBackgroundVsVpRatios(vs_vp_ratios);

  }

    //if (failed == false) { //H Writing of background models missing.
    //  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
    //    background->writeBackgrounds(timeSimbox,
    //                                 timeDepthMapping,
    //                                 timeCutMapping,
    //                                 model_settings->getFileGrid(),
    //                                 *model_settings->getTraceHeaderFormatOutput());
    //  }
    //}

  //} //i_interval

  Timings::setTimePriorExpectation(wall,cpu);

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

double CommonData::FindMeanVsVp(const NRLib::Grid<double> * vp,
                                const NRLib::Grid<double> * vs) {

  double mean = 0;
  int ni  = vp->GetNI();
  int nj  = vp->GetNJ();
  int nk  = vp->GetNK();
  for (int k=0; k < nk; k++) {
    for (int j=0; j < nj; j++) {
      for (int i=0; i < ni; i++) {
        double v1 = vp->GetValue(i,j,k);
        double v2 = vs->GetValue(i,j,k);
        mean += exp(v2-v1);
      }
    }
  }
  mean = mean/double(ni*nj*nk);

  return(mean);
}

void CommonData::SetUndefinedCellsToGlobalAverageGrid(NRLib::Grid<double> * grid,
                                                      const double          avg) {

  //float globalAvg = getAvgReal();

  //if (globalAvg == RMISSING) {
  //  calculateStatistics();
  //}

  long long count = 0;

  //setAccessMode(RANDOMACCESS);
  for (size_t k = 0 ; k < grid->GetNK(); k++) {
    for (size_t j = 0 ; j < grid->GetNJ(); j++) {
      for (size_t i = 0 ; i < grid->GetNI(); i++) {
        double value = grid->GetValue(i, j, k);
        if (value == RMISSING) {
          grid->SetValue(i, j, k, avg);
          count ++;
          //setRealValue(i,j,k,globalAvg,true);
          //if (i < nx_ && j < ny_ && k < nz_) { // Do not count undefined in padding (confusing for user)
          //  count++;
          //}
        }
      }
    }
  }
  //endAccess();

  if (count > 0) {
    long long nxyz = static_cast<long>(grid->GetNI())*static_cast<long>(grid->GetNJ())*static_cast<long>(grid->GetNK());
    LogKit::LogFormatted(LogKit::Medium, "\nThe grid contains %ld undefined grid cells (%.2f%). Setting these to global average\n",
                         count, 100.0f*static_cast<float>(count)/(static_cast<float>(count) + static_cast<float>(nxyz)),
                         avg);
  }

}

void CommonData::SubtractGrid(NRLib::Grid<double>       * to_grid,
                              const NRLib::Grid<double> * from_grid) {

  int ni = to_grid->GetNI();
  int nj = to_grid->GetNJ();
  int nk = to_grid->GetNK();

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      for (int k = 0; k < nk; k++) {
        double value_tmp = from_grid->GetValue(i, j, k);
        to_grid->SetValue(i, j, k, value_tmp);
        //to_grid(i,j,k) -= from_grid(i,j,k);
      }
    }
  }
}

void CommonData::ChangeSignGrid(NRLib::Grid<double> * grid) {

  int ni = grid->GetNI();
  int nj = grid->GetNJ();
  int nk = grid->GetNK();

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      for (int k = 0; k < nk; k++) {
        double value_tmp = grid->GetValue(i, j, k);
        value_tmp *= -1;
        grid->SetValue(i, j, k, value_tmp);
        //grid(i,j,k) *= -1;
      }
    }
  }

}


//NRLib::Grid<double>
//CommonData::FFTGridRealToGrid(const FFTGrid * fft_grid) {
//
//  NRLib::Grid<double> grid(0,0,0);
//
//  if (fft_grid != NULL) {
//    int nx = fft_grid->getNx();
//    int ny = fft_grid->getNy();
//    int nz = fft_grid->getNz();
//    //int nxp = fft_grid->getNxp();
//    int nyp = fft_grid->getNyp();
//    int nzp = fft_grid->getNzp();
//
//    int rnxp = fft_grid->getRNxp();
//
//    NRLib::Grid<double> grid(nx,ny,nz);
//
//    for (int k=0; k<nzp; k++) {
//      for (int j=0; j<nyp; j++) {
//        for (int i=0; i<rnxp; i++) {
//          if (i < nx && j < ny && k < nz)
//            grid(i,j,k) = fft_grid->getRealValue(i,j,k);
//        }
//      }
//    }
//  }
//
//  return grid;
//}


void CommonData::LoadVelocity(NRLib::Grid<double>  * velocity,
                              const Simbox         & estimation_simbox,
                              const ModelSettings  * model_settings,
                              const std::string    & velocity_field,
                              bool                 & velocity_from_inversion,
                              std::string          & err_text) {

  LogKit::WriteHeader("Setup time-to-depth relationship");

  if (model_settings->getVelocityFromInversion() == true) {
    velocity_from_inversion = true;
    return;
  }
  else if (velocity_field == "")
    return;
  else {
    const SegyGeometry      * dummy1 = NULL;
    const TraceHeaderFormat * dummy2 = NULL;
    const float               offset = model_settings->getSegyOffset(0); //Segy offset needs to be the same for all time lapse data
    std::string err_text_tmp         = "";

    std::vector<NRLib::Grid<double> *> grids(1);
    //grids[0] = velocity;
    std::vector<Simbox> simboxes;
    simboxes.push_back(estimation_simbox);

    ReadGridFromFile(velocity_field,
                     "velocity field",
                     offset,
                     grids, //velocity,
                     dummy1,
                     dummy2,
                     PARAMETER,
                     simboxes, //estimation_simbox,
                     estimation_simbox,
                     model_settings,
                     err_text_tmp);

    velocity = grids[0];

    if (err_text_tmp == "") { // No errors
      //
      // Check that the velocity grid is veldefined.
      //
      float log_min = model_settings->getAlphaMin();
      float log_max = model_settings->getAlphaMax();

      const int nz = velocity->GetNK();
      const int ny = velocity->GetNJ();
      const int nx = velocity->GetNI();
      int too_low  = 0;
      int too_high = 0;

      for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
          for (int i = 0; i < nx; i++) {
            double value = velocity->GetValue(i,j,k);
            if (value < log_min && value != RMISSING)
              too_low++;
            if (value > log_max && value != RMISSING)
              too_high++;
          }
        }
      }

      if (too_low+too_high > 0) {
        std::string text;
        text += "\nThe velocity grid used as trend in the background model of Vp";
        text += "\ncontains too small and/or too high velocities:";
        text += "\n  Minimum Vp = "+NRLib::ToString(log_min,2)+"    Number of too low values  : "+NRLib::ToString(too_low);
        text += "\n  Maximum Vp = "+NRLib::ToString(log_max,2)+"    Number of too high values : "+NRLib::ToString(too_high);
        text += "\nThe range of allowed values can changed using the ALLOWED_PARAMETER_VALUES keyword\n";
        text += "\naborting...\n";
        err_text += "Reading of file '"+velocity_field+"' for background velocity field failed.\n";
        err_text += text;
      }
    }
    else {
      err_text_tmp += "Reading of file \'"+velocity_field+"\' for background velocity field failed.\n";
      err_text += err_text_tmp;
    }
  }
}

std::map<std::string, DistributionsRock *>
CommonData::GetRockDistributionTime0() const
{
  //std::vector<std::map<std::string, DistributionsRock *> > rock_dist_t0;
  std::map<std::string, DistributionsRock *> rock_dist_t0;

  //for (size_t i = 0; i < rock_distributions_.size(); i++){
  //
  //  for (std::map<std::string, std::vector<DistributionsRock *> >::const_iterator it = rock_distributions_[i].begin(); it != rock_distributions_[i].end(); it++) {
  //    std::string name = it->first;
  //    std::vector<DistributionsRock *> rock_dist = it->second;
  //    rock_dist_t0[i][name] = rock_dist[0];
  //  }
  //}
  for (std::map<std::string, std::vector<DistributionsRock *> >::const_iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
    std::string name = it->first;
    std::vector<DistributionsRock *> rock_dist = it->second;
    rock_dist_t0[name] = rock_dist[0];
  }

  return rock_dist_t0;
}

void CommonData::GenerateRockPhysics3DBackground(const std::vector<DistributionsRock *> & rock_distribution,
                                                 const std::vector<float>               & probability,
                                                 //std::vector<NRLib::Grid<double> >      & parameters,
                                                 NRLib::Grid<double>                    * vp,
                                                 NRLib::Grid<double>                    * vs,
                                                 NRLib::Grid<double>                    * rho,
                                                 Simbox                                 & simbox,
                                                 const CravaTrend                       & trend_cube)
                                                 //int                                      i_interval)
{
  // Set up of expectations grids

  // Variables for looping through FFTGrids
  const int nz   = simbox.getnz();
  const int ny   = simbox.getny();
  const int nx   = simbox.getnx();
  //const int nzp  = vp.getNzp();
  //const int nyp  = vp.getNyp();
  //const int nxp = vp.getNxp();
  //const int rnxp = vp.getRNxp();

  LogKit::LogFormatted(LogKit::Low,"\nGenerating background model from rock physics:\n");

  float monitor_size = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float next_monitor = monitor_size;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  const size_t number_of_facies = probability.size();

  // Temporary grids for storing top and base values of (vp,vs,rho) for use in linear interpolation in the padding
  NRLib::Grid2D<float> top_vp  (nx, ny, 0.0);
  NRLib::Grid2D<float> top_vs  (nx, ny, 0.0);
  NRLib::Grid2D<float> top_rho (nx, ny, 0.0);
  NRLib::Grid2D<float> base_vp (nx ,ny, 0.0);
  NRLib::Grid2D<float> base_vs (nx, ny, 0.0);
  NRLib::Grid2D<float> base_rho(nx, ny, 0.0);

  //vp.setAccessMode(FFTGrid::WRITE);
  //vs.setAccessMode(FFTGrid::WRITE);
  //rho.setAccessMode(FFTGrid::WRITE);

  // Loop through all cells in the FFTGrids
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {

        // If outside/If in the padding in x- and y-direction,
        // set expectation equal to something at right scale
        // (top value for closest edge)
        // NBNB OK Can be made better linear interoplation between first and last value in i an j direction as well
        //if (i >= nx || j >= ny) {
        //  int indexI;
        //  int indexJ;
        //  indexI = (nx+nxp)/2 ? 0   : nx-1;
        //  indexJ = (ny+nyp)/2 ? 0   : ny-1;
        //  indexI = std::min(i,indexI);
        //  indexJ = std::min(j,indexJ);

        //  float vp_val  = top_vp(indexI,indexJ);
        //  float vs_val  = top_vs(indexI,indexJ);
        //  float rho_val = top_rho(indexI,indexJ);

        //  vp.setNextReal(vpVal);
        //  vs.setNextReal(vsVal);
        //  rho.setNextReal(rhoVal);
        //}

        // If outside in z-direction, use linear interpolation between top and base values of the expectations
        //else if (k >= nz) {
        //  double t  = double(nzp-k+1)/(nzp-nz+1);
        //  double vp_val =  top_vp(i,j)*t  + base_vp(i,j)*(1-t);
        //  double vs_val =  top_vs(i,j)*t  + base_vs(i,j)*(1-t);
        //  double rho_val = top_rho(i,j)*t + base_rho(i,j)*(1-t);

        //  // Set interpolated values in expectation grids
        //  vp.setNextReal(static_cast<float>(vp_val));
        //  vs.setNextReal(static_cast<float>(vs_val));
        //  rho.setNextReal(static_cast<float>(rho_val));
        //}

        // Otherwise use trend values to get expectation values for each facies from the rock
        //else {
        //std::vector<double> trend_position = trend_cubes_[i_interval].GetTrendPosition(i,j,k);
        std::vector<double> trend_position = trend_cube.GetTrendPosition(i,j,k);

        std::vector<float> expectations(3, 0);  // Antar initialisert til 0.

        std::vector<std::vector<double> > expectation_m(number_of_facies);
        for (size_t f = 0; f < number_of_facies; f++)
          expectation_m[f] = rock_distribution[f]->GetLogExpectation(trend_position);

        // Sum up for all facies: probability for a facies multiplied with the expectations of (vp, vs, rho) given the facies
        for (size_t f = 0; f < number_of_facies; f++){
          expectations[0] += static_cast<float>(expectation_m[f][0] * probability[f]);
          expectations[1] += static_cast<float>(expectation_m[f][1] * probability[f]);
          expectations[2] += static_cast<float>(expectation_m[f][2] * probability[f]);
        }

        // Set values in expectation grids
        //vp.setNextReal(expectations[0]);
        //vs.setNextReal(expectations[1]);
        //rho.setNextReal(expectations[2]);
        double value_tmp = static_cast<double>(expectations[0]);
        vp->SetValue(i, j, k, value_tmp);
        value_tmp = static_cast<double>(expectations[1]);
        vs->SetValue(i, j, k, value_tmp);
        value_tmp = static_cast<double>(expectations[2]);
        rho->SetValue(i, j, k, value_tmp);

        // Store top and base values of the expectations for later use in interpolation in the padded region.
        if (k==0) {
          top_vp(i,j)  = expectations[0];
          top_vs(i,j)  = expectations[1];
          top_rho(i,j) = expectations[2];
        }
        else if (k==nz-1) {
          base_vp(i,j)  = expectations[0];
          base_vs(i,j)  = expectations[1];
          base_rho(i,j) = expectations[2];
        }
        //}
      }
    }

    // Log progress
    if (k+1 >= static_cast<int>(next_monitor) && k < nz) {
      next_monitor += monitor_size;
      std::cout << "^";
      fflush(stdout);
    }
  }

  //vp.endAccess();
  //vs.endAccess();
  //rho.endAccess();
}

//-----------------------------------------------------------------------------------------------------

void CommonData::PrintExpectationAndCovariance(const std::vector<double>   & expectation,
                                               const NRLib::Grid2D<double> & covariance,
                                               const bool                  & has_trend) const{
  if (has_trend == true)
      LogKit::LogFormatted(LogKit::Low,"\nMean expectation and covariance estimated over all trend values:\n");
  else
    LogKit::LogFormatted(LogKit::Low,"\nEstimated expectation and covariance:\n");
  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Expectation            Vp        Vs       Rho\n");
  LogKit::LogFormatted(LogKit::Low,"----------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"                  %7.2f   %7.2f   %7.3f \n",
                       std::exp(expectation[0]), std::exp(expectation[1]), std::exp(expectation[2]));

  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Variances           ln Vp     ln Vs    ln Rho\n");
  LogKit::LogFormatted(LogKit::Low,"----------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"                  %.1e   %.1e   %.1e\n", covariance(0,0), covariance(1,1), covariance(2,2));

  float corr01 = static_cast<float>(covariance(0,1)/(sqrt(covariance(0,0)*covariance(1,1))));
  float corr02 = static_cast<float>(covariance(0,2)/(sqrt(covariance(0,0)*covariance(2,2))));
  float corr12 = static_cast<float>(covariance(1,2)/(sqrt(covariance(1,1)*covariance(2,2))));

  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"Corr   | ln Vp     ln Vs    ln Rho \n");
  LogKit::LogFormatted(LogKit::Low,"-------+---------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"ln Vp  | %5.2f     %5.2f     %5.2f \n",1.0f, corr01, corr02);
  LogKit::LogFormatted(LogKit::Low,"ln Vs  |           %5.2f     %5.2f \n",1.0f, corr12);
  LogKit::LogFormatted(LogKit::Low,"ln Rho |                     %5.2f \n",1.0f);
  LogKit::LogFormatted(LogKit::Low,"\n");
}

void CommonData::SetupExtendedBackgroundSimbox(Simbox   * simbox,
                                               Surface  * corr_surf,
                                               Simbox  *& bg_simbox,
                                               int        output_format,
                                               int        output_domain,
                                               int        other_output)
{
  //
  // Move correlation surface for easier handling.
  //
  Surface tmp_surf(*corr_surf);
  double avg = tmp_surf.Avg();
  if (avg > 0)
    tmp_surf.Subtract(avg);
  else
    tmp_surf.Add(avg); // This situation is not very likely, but ...

  //
  // Find top surface of background simbox.
  //
  // The funny/strange dTop->Multiply(-1.0) is due to NRLIB's current
  // inability to set dTop equal to Simbox top surface.
  //
  Surface d_top(tmp_surf);
  d_top.SubtractNonConform(&(simbox->GetTopSurface()));
  d_top.Multiply(-1.0);
  double shift_top = d_top.Min();
  Surface top_surf(tmp_surf);
  top_surf.Add(shift_top);

  //
  // Find base surface of background simbox
  //
  Surface d_bot(tmp_surf);
  d_bot.SubtractNonConform(&(simbox->GetBotSurface()));
  d_bot.Multiply(-1.0);
  double shift_bot = d_bot.Max();
  Surface bot_surf(tmp_surf);
  bot_surf.Add(shift_bot);

  //
  // Calculate number of layers of background simbox
  //
  tmp_surf.Assign(0.0);
  tmp_surf.AddNonConform(&bot_surf);
  tmp_surf.SubtractNonConform(&top_surf);
  double d_max = tmp_surf.Max();
  double dt = simbox->getdz();
  int nz;
  //
  // NBNB-PAL: I think it is a good idea to use a maximum dt of 10ms.
  //
  //if (dt < 10.0) {
  //  LogKit::LogFormatted(LogKit::High,"\nReducing sampling density for background",dt);
  //  LogKit::LogFormatted(LogKit::High," modelling from %.2fms to 10.0ms\n");
  //  dt = 10.0;  // A sampling density of 10.0ms is good enough for BG model
  // }
  nz = static_cast<int>(ceil(d_max/dt));

  //
  // Make new simbox
  //
  bg_simbox = new Simbox(simbox);
  bg_simbox->setDepth(top_surf, bot_surf, nz);

  if ((other_output & IO::EXTRA_SURFACES) > 0 && (output_domain & IO::TIMEDOMAIN) > 0) {
    std::string top_surf_name  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_BG";
    std::string base_surf_name = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_BG";
    bg_simbox->writeTopBotGrids(top_surf_name,
                                base_surf_name,
                                IO::PathToBackground(),
                                output_format);
  }
}

void CommonData::SetupExtendedBackgroundSimbox(Simbox   * simbox,
                                               Surface  * top_corr_surf,
                                               Surface  * base_corr_surf,
                                               Simbox  *& bg_simbox,
                                               int        output_format,
                                               int        output_domain,
                                               int        other_output)
{
  (void) output_format;
  (void) output_domain;
  (void) other_output;
  (void) bg_simbox;
  (void) simbox;
  (void) top_corr_surf, base_corr_surf;
  //Setup background simbox if there are two correlation surfaces
  //Similar to MultiIntervalGrid::SetupIntervalSimbox or Simbox::Simbox
  /*


  //
  // Move correlation surface for easier handling.
  //
  Surface tmp_surf(*corr_surf);
  double avg = tmp_surf.Avg();
  if (avg > 0)
    tmp_surf.Subtract(avg);
  else
    tmp_surf.Add(avg); // This situation is not very likely, but ...

  //
  // Find top surface of background simbox.
  //
  // The funny/strange dTop->Multiply(-1.0) is due to NRLIB's current
  // inability to set dTop equal to Simbox top surface.
  //
  Surface dTop(tmp_surf);
  dTop.SubtractNonConform(&(simbox->GetTopSurface()));
  dTop.Multiply(-1.0);
  double shiftTop = dTop.Min();
  Surface topSurf(tmp_surf);
  topSurf.Add(shiftTop);

  //
  // Find base surface of background simbox
  //
  Surface dBot(tmp_surf);
  dBot.SubtractNonConform(&(simbox->GetBotSurface()));
  dBot.Multiply(-1.0);
  double shiftBot = dBot.Max();
  Surface botSurf(tmp_surf);
  botSurf.Add(shiftBot);

  //
  // Calculate number of layers of background simbox
  //
  tmp_surf.Assign(0.0);
  tmp_surf.AddNonConform(&botSurf);
  tmp_surf.SubtractNonConform(&topSurf);
  double dMax = tmp_surf.Max();
  double dt = simbox->getdz();
  int nz;
  //
  // NBNB-PAL: I think it is a good idea to use a maximum dt of 10ms.
  //
  //if (dt < 10.0) {
  //  LogKit::LogFormatted(LogKit::High,"\nReducing sampling density for background",dt);
  //  LogKit::LogFormatted(LogKit::High," modelling from %.2fms to 10.0ms\n");
  //  dt = 10.0;  // A sampling density of 10.0ms is good enough for BG model
  // }
  nz = static_cast<int>(ceil(dMax/dt));

  //
  // Make new simbox
  //
  bg_simbox = new Simbox(simbox);
  bg_simbox->setDepth(topSurf, botSurf, nz);

  if ((other_output & IO::EXTRA_SURFACES) > 0 && (output_domain & IO::TIMEDOMAIN) > 0) {
    std::string top_surf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_BG";
    std::string base_surf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_BG";
    bg_simbox->writeTopBotGrids(top_surf,
                                base_surf,
                                IO::PathToBackground(),
                                output_format);
  }
  */
}
//}

bool CommonData::SetupPriorCorrelation(const ModelSettings * model_settings,
                                       const InputFiles * input_files,
                                       const std::vector<NRLib::Well> & wells,
                                       const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs_for_correlation,
                                       const std::vector<Simbox> & interval_simboxes,
                                       const std::vector<std::vector<float> > & prior_facies_prob,
                                       const std::vector<std::string> & facies_names,
                                       const std::vector<CravaTrend> & trend_cubes,
                                       const std::map<int, std::vector<SeismicStorage> > & seismic_data,
                                       const std::vector<std::vector<NRLib::Grid<double> *> > & background,
                                       std::string & err_text_common){

  (void) seismic_data, facies_names, prior_facies_prob;
  LogKit::WriteHeader("Setup of Prior Correlation");

  bool failed = false;
  std::string err_text = "";
  Analyzelog * analyze_all = NULL;
  int min_data_for_estimation = 100; ///< Minimum number of blocks for estimation of time and param correlation
  bool print_result = ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0 ||
                                              model_settings->getEstimationMode() == true);
  size_t n_intervals = interval_simboxes.size();
  size_t n_facies = facies_names.size(); //static_cast<int>(prior_facies_[0].size()); //H Changed since prior_facies_ isn't necessarily set.
  std::vector<std::string> interval_names = model_settings->getIntervalNames();

  // Local parameters to be set in this function ----------------------------------------------------
  prior_corr_T_.resize(n_intervals);
  prior_param_cov_.resize(n_intervals);
  prior_corr_XY_.resize(n_intervals);

  if (model_settings->getDoInversion() || print_result)
  {
    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    // A Parameter covariance -----------------------------------------------------------------------

    // 1. From file
    // 2. From rock physics
    // 3. Estimate from data

    // Get parameter correlation file
    const std::string & param_cov_file = input_files->getParamCovFile();
    // If the parameter correlation file is empty and rock physics is being used
    bool param_cov_from_rock_physics = (param_cov_file == "" && model_settings->getFaciesProbFromRockPhysics());
    // If there is no parameter corr file and no rock physics, estimate parameter correlation
    bool estimate_param_cov = (param_cov_file == "" && !param_cov_from_rock_physics);

    //
    // Read parameter covariance (Var0) from file or set from output from function generateRockPhysics3DBackground.
    // Consistency check that only one option (file or rock physics) is possible, is done in XmlModelFile::checkInversionConsistency
    //
    //std::vector<NRLib::Grid2D<double> > param_corr_array;
    float ** temp_array;
    bool failed_param_cov = false;
    std::string tmp_err_text = "";

    // 1. If there is a param corr file defined, use this
    if(param_cov_file != "") {
      //prior_param_corr_.resize(n_intervals);
      temp_array = ReadMatrix(param_cov_file, 3, 3, "parameter covariance", tmp_err_text);
      ValidateCovarianceMatrix(temp_array, model_settings, tmp_err_text);
      if(temp_array == NULL || tmp_err_text != "") {
        err_text += "Reading of file " + param_cov_file + " for parameter covariance matrix failed\n";
        err_text += tmp_err_text;
        failed_param_cov = true;
      }
      // copy
      else{
        for(size_t i = 0; i<n_intervals; i++){
          prior_param_cov_[i].resize(3,3);
          for(size_t j = 0; j<3; j++){
            for(size_t k = 0; k<3; k++){
              // parameter covariance is the same for all intervals
              prior_param_cov_[i](j,k) = temp_array[j][k];
            }
          }
        }
      }
    }
    // 2. If there is no param correlation file and rock physics is used,
    // get param corr at lag 0 from the rock physics distribution.
    // (One 3x3 correlation matrix per interval)
    else if (param_cov_from_rock_physics) {

      temp_array = new float * [3];
      for(size_t i=0;i<3;i++) {
          temp_array[i] = new float[3];
      }

      for (size_t i = 0; i<n_intervals; i++){
        prior_param_cov_[i].resize(3,3);
      }

      std::vector<DistributionsRock *> rock_distribution(n_intervals);
      typedef std::map<std::string, DistributionsRock *> Rp_Map_Type;
      Rp_Map_Type rp_map = GetRockDistributionTime0();

      for(size_t i=0; i<n_facies; i++) {
        std::map<std::string, DistributionsRock *>::iterator iter = rp_map.find(facies_names_[i]);
        if(iter != rp_map.end()){
          rock_distribution[i] = iter->second;
        }
        else{
          err_text += "Rock distribution not found for facies '"+facies_names_[i]+"'\n";
        }
      }

      // if multiple interval settings
      if (interval_names.size()>0){
        for(size_t i = 0; i<n_intervals; i++){
          CalculateCovarianceFromRockPhysics(rock_distribution,
                                               model_settings->getPriorFaciesProbInterval(interval_names[i]),
                                               facies_names_,
                                               trend_cubes[i],
                                               prior_param_cov_[i],
                                               err_text);
        }
      }
      // if multiple interval settings is not being used
      else{
        CalculateCovarianceFromRockPhysics(rock_distribution,
                                             model_settings->getPriorFaciesProb(),
                                             facies_names_,
                                             trend_cubes[0],
                                             prior_param_cov_[0],
                                             err_text);
      }

      if (tmp_err_text != "")
      {
        err_text += "Parameter covariance matrix from rock physics failed\n";
        err_text += tmp_err_text;
        failed_param_cov = true;
      }
      // 3 Parameter estimation from well data -> further down
    }

    //
    // B Lateral correlation -------------------------------------------------------------------
    //

    for (size_t i = 0; i < n_intervals; i++){
      prior_corr_XY_[i] = FindCorrXYGrid(&(interval_simboxes[i]), model_settings);
    }

    //
    // C Temporal correlation -------------------------------------------------------------------

    // 1. From file
    // 2. Correlation range
    // 3. Estimate from data

    const std::string & corr_time_file = input_files->getTempCorrFile();
    bool temporal_corr_range_given = (input_files->getTempCorrFile() == "" && model_settings->getUseVerticalVariogram());
    bool estimate_temp_corr = (corr_time_file == "" && model_settings->getUseVerticalVariogram() == false);
    bool failed_temp_corr = false;
    failed_param_cov = false;

    for (size_t i = 0; i<n_intervals; i++){
      if(!failed_temp_corr && !failed_param_cov){

        // Number of bins
        int n_corr_T = interval_simboxes[i].GetNZpad();
        if((n_corr_T % 2) == 0)
          n_corr_T = n_corr_T/2+1;
        else
          n_corr_T = n_corr_T/2;

        // 2. Use variogram with correlation range
        if(temporal_corr_range_given == true) {
          prior_corr_T_[i].resize(n_corr_T+1,0);
          float temp_corr_range = model_settings->getTemporalCorrelationRange();
          float dz = static_cast<float>(interval_simboxes[i].getdz());
          for(int j=0; j<=n_corr_T; j++){
            //using an exponential variogram with a = 1/3 (Chiles and Delfiner 1999)
            prior_corr_T_[i][j] = exp(-3*dz*j/temp_corr_range);
          }
        }
        // 1. Defined in text file
        else{
          std::string tmp_err_text("");
          float ** corr_mat = ReadMatrix(corr_time_file, 1, n_corr_T+1, "temporal correlation", tmp_err_text);
          if(corr_mat == NULL)
          {
            err_text += "Reading of file '"+corr_time_file+"' for temporal correlation failed\n";
            err_text += tmp_err_text;
            failed_temp_corr = true;
          }
          prior_corr_T_[i].resize(n_corr_T,0);
          if (!failed_temp_corr)
          {
            for(int j=0;j<n_corr_T;j++)
              prior_corr_T_[i][j] = corr_mat[0][j+1];
            delete [] corr_mat[0];
            delete [] corr_mat;
          }
        }

        //
        // Estimation of parameter and temporal correlation ----------------------------------------
        //

        // A.3 Estimation of parameter correlation from data
        // C.3 Estimation of temporal correlation from data


        //float ** point_var_0 = NULL;
        if (estimate_param_cov || estimate_temp_corr){

          std::string tmp_err_txt;
          // Option 1: Estimate within this interval
          Analyzelog * analyze = new Analyzelog(wells, mapped_blocked_logs_for_correlation, background[i], &interval_simboxes[i], model_settings, tmp_err_txt);
          if (tmp_err_txt != "") {
            err_text += tmp_err_txt;
            failed_param_cov = true;
          }
          // Option 2: Estimate over all intervals if the multiple interval setting is being used
          else if(analyze->GetEnoughData() == false && interval_names.size() > 0 && analyze_all == NULL){
            analyze_all = new Analyzelog(wells, mapped_blocked_logs_for_correlation, background, interval_simboxes, model_settings, tmp_err_txt);
            if(analyze_all->GetEnoughData() == false){
              err_text += "There are not enough layers in the inversion intervals to estimate prior correlations.\n";
              failed_param_cov = true;
              failed_temp_corr = true;
            }
          }

          if(err_text == ""){
            NRLib::Matrix param_cov_array_temp;
            if(estimate_param_cov){
              if (analyze->GetEnoughData() == true){
                param_cov_array_temp = analyze->GetVar0();
              }
              else{
                param_cov_array_temp = analyze_all->GetVar0();
              }
            }

            for (int j=0; j<3; j++){
              for (int k=0; k<3; k++)
                prior_param_cov_[i](j,k) = param_cov_array_temp(j,k);
            }

            //prior_var_0_[i] = analyze->getPointVar0();

            std::vector<float> est_corr_T;
            if (analyze->GetEnoughData() == true){
              est_corr_T = analyze->GetCorrT();
            }
            else{
              est_corr_T = analyze_all->GetCorrT();
            }
            if(estimate_temp_corr) {
              //corr_T.resize(n_corr_T);
              int n_est = analyze->GetNumberOfLags();
              int j, max = n_est;
              if(max > n_corr_T)
                max = n_corr_T;
              for(j = 0; j < max; j++)
                prior_corr_T_[i][j] = est_corr_T[j];
              if(j<n_corr_T) {
                LogKit::LogFormatted(LogKit::High,
                  "\nOnly able to estimate %d of %d lags needed in temporal correlation. The rest are set to 0.\n", n_est, n_corr_T);
                for(;j<n_corr_T;j++)
                  prior_corr_T_[i][j] = 0.0f;
              }
            }
          }
          delete analyze;

        }

        if (failed_param_cov || failed_temp_corr)
          failed = true;


        if (!failed) {

          CheckCovarianceParameters(prior_param_cov_[i]);

          //if(print_result)
          // WriteFilePriorVariances(model_settings, corr_T, prior_corr_XY_[i], dt);
          //PrintPriorVariances();

        }
      }
    }

    // Delete parameter correlation matrices

    for(int j=0; j<3; j++)
      delete [] temp_array[j];
    delete [] temp_array;


    if(failed_temp_corr == true || failed_param_cov == true)
    {
      err_text += "Could not construct prior correlation. Unknown why...\n";
      failed = true;
    }

    Timings::setTimePriorCorrelation(wall,cpu);
  }

  delete analyze_all;

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}
 /*
  const std::vector<NRLib::Grid<double> > & CommonData::GetCovParametersInterval(int i_interval){
    return cov_params_interval_;
  }

  const std::vector<NRLib::Grid<double> > & CommonData::GetCorrParametersInterval(int i_interval){
    return corr_params_interval_;
  }

  const NRLib::Matrix                     & CommonData::GetPriorVar0(int i_interval){
    return prior_var_0_;
  }
  */

void CommonData::CheckCovarianceParameters(NRLib::Matrix & param_cov)
{
  // check that covariance is well conditioned and robustify
  NRLib::Vector e_vals(3);
  NRLib::Matrix e_vec(3,3);
  NRLib::Matrix tmp(3,3);
  tmp = param_cov;
  NRLib::ComputeEigenVectors(tmp, e_vals,e_vec);

  double max_val = e_vals(0);

  for(int i=1;i<3;i++)
    if(max_val<e_vals(i))
      max_val=e_vals(i);

  for(int k=0;k<3;k++){
    if(e_vals(k)<max_val*0.001){
      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
          param_cov(i,j)+= e_vec(k,i)*e_vec(k,j)*(max_val*0.0011-e_vals(k));
    }
  }

}

// -------------------------------------------------------------------------------------------------
void  CommonData::GetCorrGradIJ(float         & corr_grad_I,
                                float         & corr_grad_J,
                                const Simbox  * simbox) const{
  double angle  = simbox->getAngle();
  double cosrot = cos(angle);
  double sinrot = sin(angle);
  double dx     = simbox->getdx();
  double dy     = simbox->getdy();
  double grad_x = simbox->GetGradX();
  double grad_y = simbox->GetGradY();

  double cI =  dx*cosrot*grad_x + dy*sinrot*grad_y;
  double cJ = -dx*sinrot*grad_y + dy*cosrot*grad_y;

  corr_grad_I = float(cI/simbox->getdz());
  corr_grad_J = float(cJ/simbox->getdz());
}

void CommonData::ValidateCovarianceMatrix(float               ** C,
                                           const ModelSettings * model_settings,
                                           std::string         & err_txt){
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
    err_txt += "The prior Vp variance is outside valid range:\n";
    err_txt += "  Given value   : " + NRLib::ToString(C00) + "\n";
    err_txt += "  Minimum value : " + NRLib::ToString(minAlpha) + "\n";
    err_txt += "  Maximum value : " + NRLib::ToString(maxAlpha) + "\n";
  }
  if (C11 < minBeta || C11 > maxBeta) {
    err_txt += "The prior Vs variance is outside valid range:\n";
    err_txt += "  Given value   : " + NRLib::ToString(C11) + "\n";
    err_txt += "  Minimum value : " + NRLib::ToString(minBeta) + "\n";
    err_txt += "  Maximum value : " + NRLib::ToString(maxBeta) + "\n";
  }
  if (C22 < minRho || C22 > maxRho) {
    err_txt += "The prior density variance is outside valid range:\n";
    err_txt += "  Given value   : " + NRLib::ToString(C22) + "\n";
    err_txt += "  Minimum value : " + NRLib::ToString(minRho) + "\n";
    err_txt += "  Maximum value : " + NRLib::ToString(maxRho) + "\n";
  }

  float corr01 = C01/(std::sqrt(C00)*std::sqrt(C11));
  float corr02 = C02/(std::sqrt(C00)*std::sqrt(C22));
  float corr12 = C12/(std::sqrt(C11)*std::sqrt(C22));

  if (corr01 < -1.0 || corr01 > 1.0) {
    err_txt += "The prior Vp-Vs correlation is illegal (" + NRLib::ToString(corr01) + ")\n";
  }
  if (corr02 < -1.0 || corr02 > 1.0) {
    err_txt += "The prior Vp-Rho correlation is illegal (" + NRLib::ToString(corr02) + ")\n";
  }
  if (corr12 < -1.0 || corr12 > 1.0) {
    err_txt += "The prior Vs-Rho correlation is illegal (" + NRLib::ToString(corr12) + ")\n";
  }

  if (std::abs(C01 - C10) > 0.0f) {
    err_txt += "The prior covariance matrix is not symmetric in Vp and Vs\n";
    err_txt += "  Corr(Vp,Vs) : " + NRLib::ToString(C01) + "\n";
    err_txt += "  Corr(Vs,Vp) : " + NRLib::ToString(C10) + "\n";
  }
  if (std::abs(C02 - C20) > 0.0f) {
    err_txt += "The prior covariance matrix is not symmetric in Vp and Rho\n";
    err_txt += "  Corr(Vp,Rho) : " + NRLib::ToString(C02) + "\n";
    err_txt += "  Corr(Rho,Vp) : " + NRLib::ToString(C20) + "\n";
  }
  if (std::abs(C12 - C21) > 0.0f) {
    err_txt += "The prior covariance matrix is not symmetric in Vs and Rho\n";
    err_txt += "  Corr(Vs,Rho) : " + NRLib::ToString(C12) + "\n";
    err_txt += "  Corr(Rho,Vs) : " + NRLib::ToString(C21) + "\n";
  }
}

void  CommonData::CalculateCovarianceFromRockPhysics(const std::vector<DistributionsRock *>           & rock_distribution,
                                                     const std::map<std::string, float>               & probability,
                                                     const std::vector<std::string>                   & facies_names,
                                                     const CravaTrend                                 & trend_cubes,
                                                     NRLib::Matrix                                    & param_cov,
                                                     std::string                                      & err_txt){

  LogKit::LogFormatted(LogKit::Low,"\nGenerating covariances from rock physics\n");

  // order the facies prob vector after the facies_names vector
  std::vector<float> facies_prob(facies_names.size());
  for(unsigned int i=0; i<facies_names.size(); i++)
    facies_prob[i] = probability.find(facies_names[i])->second;

  bool has_trend = false;
  for (size_t i=0; i<rock_distribution.size(); i++) {
    std::vector<bool> rock_has_trend = rock_distribution[i]->HasTrend();

    for (int j=0; j<2; j++) {
      if (rock_has_trend[j] == true)
        has_trend = true;
    }
  }

  if (has_trend == true) {

    std::vector<int> trend_cube_size = trend_cubes.GetSizeTrendCubes();

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

            std::vector<double> trend_position = trend_cubes.GetTrendPosition(i,j,k);

            NRLib::Grid2D<double> sigma_sum(3,3,0);

            CalculateCovarianceInTrendPosition(rock_distribution,
                                               facies_prob,
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
          param_cov(i,j) = sumVariance(i,j)/n_samples;
      }
    }

    else
      err_txt += "Could not build a covariance structure from rock physics.\n";
  }

  else {
    std::vector<double> trend_position(2, 0.0);

    NRLib::Grid2D<double> sigma_sum(3,3,0);

    CalculateCovarianceInTrendPosition(rock_distribution,
                                       facies_prob,
                                       trend_position,
                                       sigma_sum);


    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++)
        param_cov(i,j) = sigma_sum(i,j);
    }

  }
}

void  CommonData::CalculateCovarianceInTrendPosition(const std::vector<DistributionsRock *> & rock_distribution,
                                                     const std::vector<float>               & probability,
                                                     const std::vector<double>              & trend_position,
                                                     NRLib::Grid2D<double>                  & sigma_sum) const{
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

Surface * CommonData::FindCorrXYGrid(const Simbox           * time_simbox,
                                     const ModelSettings    * model_settings) const{

  float dx  = static_cast<float>(time_simbox->getdx());
  float dy  = static_cast<float>(time_simbox->getdy());

  int   nx  = time_simbox->GetNXpad();
  int   ny  = time_simbox->GetNYpad();

  Surface * grid = new Surface(0, 0, dx*nx, dy*ny, nx, ny, RMISSING);

  if (model_settings->getLateralCorr()!=NULL) // NBNB-PAL: Denne her blir aldri null etter at jeg la inn en default lateral correlation i modelsettings.
  {
    int refi,refj;
    for (int j=0;j<ny;j++)
    {
      for (int i=0;i<nx;i++)
      {
        if (i<(nx/2+1))
        {
          refi = i;
        }
        else
        {
          refi = i-nx;
        }
        if (j< (ny/2+1))
        {
          refj = j;
        }
        else
        {
          refj = j-ny;
        }
        (*grid)(j*nx+i) = model_settings->getLateralCorr()->corr(refi*dx, refj*dy);
      }
    }
  }
  return(grid);
}

//void CommonData::SetCorrelationParameters(NRLib::Grid2D<double>     & param_corr,
//                                          const std::vector<float>  & prior_corr_t,
//                                          Surface                   * prior_corr_XY,
//                                          const int                 & min_int_fq,
//                                          const float               & corr_grad_i,
//                                          const float               & corr_grad_j,
//                                          const int                 & nx,
//                                          const int                 & ny,
//                                          const int                 & nz,
//                                          const int                 & nx_pad,
//                                          const int                 & ny_pad,
//                                          const int                 & nz_pad,
//                                          int                         i_interval)
//{
//  prior_var_0_[i_interval].resize(3,3);
//
//  prior_var_0_[i_interval](0,0) = static_cast<double>(param_corr(0,0));
//  prior_var_0_[i_interval](1,0) = static_cast<double>(param_corr(1,0));
//  prior_var_0_[i_interval](2,0) = static_cast<double>(param_corr(2,0));
//  prior_var_0_[i_interval](0,1) = static_cast<double>(param_corr(0,1));
//  prior_var_0_[i_interval](1,1) = static_cast<double>(param_corr(1,1));
//  prior_var_0_[i_interval](2,1) = static_cast<double>(param_corr(2,1));
//  prior_var_0_[i_interval](0,2) = static_cast<double>(param_corr(0,2));
//  prior_var_0_[i_interval](1,2) = static_cast<double>(param_corr(1,2));
//  prior_var_0_[i_interval](2,2) = static_cast<double>(param_corr(2,2));
//
//
//  // check if covariance is well conditioned and robustify
//  NRLib::Vector e_vals(3);
//  NRLib::Matrix e_vec(3,3);
//  NRLib::Matrix tmp(3,3);
//  tmp=prior_var_0_;
//  NRLib::ComputeEigenVectors(tmp, e_vals,e_vec);
//
//  double max_val = e_vals(0);
//
//  for(int i=1;i<3;i++)
//    if(max_val<e_vals(i))
//      max_val=e_vals(i);
//
//  for(int k=0;k<3;k++){
//    if(e_vals(k)<max_val*0.001){
//      for(int i=0;i<3;i++)
//        for(int j=0;j<3;j++)
//          prior_var_0_[i_interval](i,j)+= e_vec(k,i)*e_vec(k,j)*(max_val*0.0011-e_vals(k));
//    }
//  }
//  //tmp=priorVar0_;
//  //NRLib::ComputeEigenVectors(tmp,eVals,eVec);
//
//  createCorrGrids(nx, ny, nz, nxPad, nyPad, nzPad, false);
//
//  initializeCorrelations(priorCorrXY,
//                         priorCorrT,
//                         corrGradI,
//                         corrGradJ,
//                         minIntFq,
//                         nzPad);
//}


bool CommonData::SetupTimeLine(ModelSettings * model_settings,
                               InputFiles    * input_files,
                               std::string   & err_text_common) {

  //Set up timeline.
  time_line_ = new TimeLine();
  std::string err_text = "";
  //Activate below when gravity data are ready.
  //Do gravity first.
  //for (int i=0;i<modelSettings->getNumberOfGravityData();i++) {
  //  int time = computeTime(modelSettings->getGravityYear[i],
  //                         modelSettings->getGravityMonth[i],
  //                         modelSettings->getGravityDay[i]);
  //  timeLine_->AddEvent(time, TimeLine::GRAVITY, i);

  bool first_gravimetric_event = true;
  for (int i=0; i < model_settings->getNumberOfVintages(); i++) {
    //Vintages may have both travel time and AVO

    try {
      int time = ComputeTime(model_settings->getVintageYear(i),
                              model_settings->getVintageMonth(i),
                              model_settings->getVintageDay(i));
        // Do gravity first
        if (model_settings->getGravityTimeLapse(i)){
          if (first_gravimetric_event){
            // Do not save first gravity event in timeline
            first_gravimetric_event = false;
          }
          else{
            time_line_->AddEvent(time, TimeLine::GRAVITY, i);
          }
      }

      //Activate below when travel time is ready.
      //Travel time ebefore AVO for same vintage.
      //if (travel time for this vintage)
      //timeLine_->AddEvent(time, TimeLine::TRAVEL_TIME, i);
      //Travel time before AVO for same vintage.

      if (model_settings->getNumberOfAngles(i) > 0) //Check for AVO data, could be pure travel time.
        time_line_->AddEvent(time, TimeLine::AVO, i);
    }
    catch(NRLib::Exception & e) {
      err_text += "Error setting up TimeLine: " + std::string(e.what());
    }

  }

  //H 4D-part set ut later (3 b) 10)

  if (err_text != "") {
    err_text_common += err_text;
    return false;
  }

  return true;
}

bool CommonData::SetupGravityInversion(ModelSettings * model_settings,
                                       InputFiles    * input_files,
                                       std::string   & err_text_common) {

  //H Intervals? It is set up per timelapse, divide it it per intevals?
  //  Seperate based on depth and simox-intevals?
  //  Alternative: Read all here -> split up in modelGravityStatic/modelGravityDynamic where these are used.

  std::string err_text = "";

  bool failed                 = false;
  bool before_injection_start = false; // When do we know what this should be??

  bool do_gravity_inversion = true;
  int number_gravity_files  = 0;
  for (int i = 0; i < model_settings->getNumberOfVintages(); i++) {
    if (model_settings->getGravityTimeLapse(i))
      number_gravity_files++;
  }

  if (number_gravity_files == 0) {
    // Everything is ok - we do not need gravity inversion
    do_gravity_inversion = false;
  }

  if (number_gravity_files == 1) {
    do_gravity_inversion = false;
    err_text+="Need at least two gravity surveys for inversion.";
  }

  int n_timelapses = model_settings->getNumberOfTimeLapses();

  observation_location_utmx_.resize(n_timelapses);
  observation_location_utmy_.resize(n_timelapses);
  observation_location_depth_.resize(n_timelapses);
  gravity_response_.resize(n_timelapses);
  gravity_std_dev_.resize(n_timelapses);

  // Set up gravimetric baseline
  if (do_gravity_inversion) {

    for (int i_timelapse = 0; i_timelapse < n_timelapses; i_timelapse++) {

      // Find gravity data file for this timelapse
      std::string fileName = input_files->getGravimetricData(i_timelapse);

      int n_obs = 30;     //user input
      int n_columns = 5;  // We require data files to have five columns

      std::vector<float> observation_location_utmx_tmp(n_obs);
      std::vector<float> observation_location_utmy_tmp(n_obs);
      std::vector<float> observation_location_depth_tmp(n_obs);
      std::vector<float> gravity_response_tmp(n_obs);
      std::vector<float> gravity_std_dev_tmp(n_obs);

      if (i_timelapse == 0) {

        ReadGravityDataFile(fileName, "gravimetric base survey",
                            n_obs, n_columns,
                            observation_location_utmx_tmp,
                            observation_location_utmy_tmp,
                            observation_location_depth_tmp,
                            gravity_response_tmp,
                            gravity_std_dev_tmp,
                            err_text);
      }
      else if (i_timelapse <= model_settings->getNumberOfVintages()) {

        ReadGravityDataFile(fileName,
                            "gravimetric survey ", // +thisTimeLapse_,
                            n_obs, n_columns,
                            observation_location_utmx_tmp,
                            observation_location_utmy_tmp,
                            observation_location_depth_tmp,
                            gravity_response_tmp,
                            gravity_std_dev_tmp,
                            err_text);
      }

      observation_location_utmx_.push_back(observation_location_utmx_tmp);
      observation_location_utmy_.push_back(observation_location_utmy_tmp);
      observation_location_depth_.push_back(observation_location_depth_tmp);
      gravity_response_.push_back(gravity_response_tmp);
      gravity_std_dev_.push_back(gravity_std_dev_tmp);

    }
  }

  if (err_text != "") {
    err_text_common += "Error(s) with gravimetric surveys";
    err_text_common += err_text;
    return false;
  }

  return true;
}

void CommonData::ReadGravityDataFile(const std::string   & file_name,
                                     const std::string   & read_reason,
                                     int                   n_obs,
                                     int                   n_columns,
                                     std::vector <float> & obs_loc_utmx,
                                     std::vector <float> & obs_loc_utmy,
                                     std::vector <float> & obs_loc_depth,
                                     std::vector <float> & gravity_response,
                                     std::vector <float> & gravity_std_dev,
                                     std::string         & err_text) {

  float * tmpRes = new float[n_obs*n_columns+1];
  std::ifstream in_file;
  NRLib::OpenRead(in_file, file_name);
  std::string text = "Reading "+read_reason+" from file "+file_name+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
  std::string storage;
  int index = 0;
  bool failed = false;

  while(failed == false && in_file >> storage) {
    if (index < n_obs*n_columns) {
      try {
        tmpRes[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        err_text += "Error in "+file_name+"\n";
        err_text += e.what();
        failed = true;
      }
    }
    index++;
  }
  if (failed == false) {
    if (index != n_obs*n_columns) {
      failed = true;
      err_text += "Found "+NRLib::ToString(index)+" in file "+file_name+", expected "+NRLib::ToString(n_obs*n_columns)+".\n";
    }
  }

  if (failed == false) {
    LogKit::LogFormatted(LogKit::Low,"ok.\n");
    index = 0;
    for (int i=0; i < n_obs; i++) {
      obs_loc_utmx[i] = tmpRes[index];
      index++;
      obs_loc_utmy[i] = tmpRes[index];
      index++;
      obs_loc_depth[i] = tmpRes[index];
      index++;
      gravity_response[i] = tmpRes[index];
      index++;
      gravity_std_dev[i] = tmpRes[index];
      index++;
    }
  }
  else{
    failed = true;
    LogKit::LogFormatted(LogKit::Low,"failed.\n");
  }
  delete [] tmpRes;
}

bool CommonData::SetupTravelTimeInversion(ModelSettings * model_settings,
                                          InputFiles    * input_files,
                                          std::string   & err_text_common) {

  //This is from ModelTravelTimeDynamic. Need to add from ModelTravelTimeStatic when it is added.
  std::string err_text = "";

  bool failed_surfaces = false;
  int n_timelapses     = model_settings->getNumberOfTimeLapses();
  horizons_.resize(n_timelapses);

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall, cpu);

  for (int i_timelapse = 0; i_timelapse < n_timelapses; i_timelapse++) {

    ProcessHorizons(horizons_[i_timelapse],
                    input_files,
                    err_text,
                    failed_surfaces,
                    i_timelapse);

    //From ModelTravelTimeDynamic::processRmsData

    LogKit::WriteHeader("Reading RMS travel time data for timelapse " + NRLib::ToString(i_timelapse));
    bool failed = false;
    const SegyGeometry * geometry = new const SegyGeometry;
    //FFTGrid * rms_data            = new FFTGrid;
    NRLib::Grid<double> * rms_data = NULL;

    float offset = model_settings->getTravelTimeSegyOffset(i_timelapse);
    if (offset < 0)
      offset = model_settings->getSegyOffset(i_timelapse);

    const std::string & file_name = input_files->getRmsVelocities(i_timelapse);
    std::string         data_name = "RMS data";
    std::string         tmp_err_text = "";

    std::vector<NRLib::Grid<double> *> grids(1);
    grids[0] = rms_data;
    std::vector<Simbox> simboxes;
    simboxes.push_back(estimation_simbox_);

    ReadGridFromFile(file_name,
                     data_name,
                     offset,
                     grids, //rms_data,
                     geometry,
                     model_settings->getTravelTimeTraceHeaderFormat(i_timelapse),
                     DATA,
                     simboxes,
                     estimation_simbox_,
                     model_settings,
                     tmp_err_text);

    //rms_data = grids[0];

    if (tmp_err_text != "") {
      tmp_err_text += "\nReading of file \'"+file_name+"\' for "+data_name+" for timelapse " + NRLib::ToString(i_timelapse) + "failed.\n";
      err_text += tmp_err_text;
      failed = true;
    }

    if (failed == false)
      rms_data_.push_back(rms_data);

    LogKit::LogFormatted(LogKit::Low,"\n");

    if (failed == false) {
      bool segy_volumes_read = false;

      if (geometry != NULL)
        segy_volumes_read = true;

      if (segy_volumes_read == true) {
        LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx         ly     azimuth         dx      dy\n");
        LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");

        if (geometry != NULL) {
          double geo_angle = (-1)*estimation_simbox_.getAngle()*(180/M_PI);
          if (geo_angle < 0)
            geo_angle += 360.0;
          LogKit::LogFormatted(LogKit::Low,"RMS travel time data   %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                               geometry->GetX0(), geometry->GetY0(),
                               geometry->Getlx(), geometry->Getly(), geo_angle,
                               geometry->GetDx(), geometry->GetDy());
        }
      }

      //if ((modelSettings->getOutputGridsOther() & IO::RMS_VELOCITIES) > 0) {
      //  std::string baseName = IO::PrefixTravelTimeData();
      //  std::string sgriLabel = std::string("RMS travel time data");

      //  rms_data->writeFile(baseName,
      //                      IO::PathToTravelTimeData(),
      //                      timeSimbox,
      //                      sgriLabel,
      //                      offset,
      //                      timeDepthMapping,
      //                      timeCutMapping,
      //                      *modelSettings->getTraceHeaderFormatOutput());
      //}

      if (geometry != NULL)
        delete geometry;
    }

    Timings::setTimeSeismic(wall, cpu);
  } //i_timelapse

  if (err_text != "") {
    err_text_common += "Error(s) while loadting travel time data";
    err_text_common += err_text;
    return false;
  }

  return true;
}

void CommonData::ProcessHorizons(std::vector<Surface>   & horizons,
                                 const InputFiles       * input_files,
                                 std::string            & err_text,
                                 bool                   & failed,
                                 int                      i_timelapse)
{
  const std::vector<std::string> & travel_time_horizons = input_files->getTravelTimeHorizons(i_timelapse);

  int n_horizons = static_cast<int>(travel_time_horizons.size());

  if (n_horizons == 1) {
    if (travel_time_horizons[0] != "") {
      err_text += "Only one surface is given for inversion of the horizons in the travel time data. At least two surfaces should be given\n";
      failed = true;
    }
  }

  else {
    horizons.resize(n_horizons);
    for (int i=0; i<n_horizons; i++)
      horizons[i] = Surface(travel_time_horizons[i]);
  }

}

// --------------------------------------------------------------------------------
void  CommonData::EstimateXYPaddingSizes(Simbox          * interval_simbox,
                                         ModelSettings   * model_settings) const{
  double dx      = interval_simbox->getdx();
  double dy      = interval_simbox->getdy();
  double lx      = interval_simbox->getlx();
  double ly      = interval_simbox->getly();
  int    nx      = interval_simbox->getnx();
  int    ny      = interval_simbox->getny();
  //int    nz      = interval_simbox->getnz();

  // The padding factor in model settings is the user specified padding factor
  double x_pad_factor = model_settings->getXPadFac();
  double y_pad_factor = model_settings->getYPadFac();
  double x_pad    = x_pad_factor*lx;
  double y_pad    = y_pad_factor*ly;

  if (model_settings->getEstimateXYPadding())
  {
    float  range1 = model_settings->getLateralCorr()->getRange();
    float  range2 = model_settings->getLateralCorr()->getSubRange();
    float  angle  = model_settings->getLateralCorr()->getAngle();
    double factor = 0.5;  // Lateral correlation is not very important. Half a range is probably more than enough

    x_pad           = factor * std::max(fabs(range1*cos(angle)), fabs(range2*sin(angle)));
    y_pad           = factor * std::max(fabs(range1*sin(angle)), fabs(range2*cos(angle)));
    x_pad           = std::max(x_pad, dx);     // Always require at least one grid cell
    y_pad           = std::max(y_pad, dy);     // Always require at least one grid cell
    x_pad_factor    = std::min(1.0, x_pad/lx); // A padding of more than 100% is insensible
    y_pad_factor    = std::min(1.0, y_pad/ly);
  }

  int nx_pad = MultiIntervalGrid::FindPaddingSize(nx, x_pad_factor);
  int ny_pad = MultiIntervalGrid::FindPaddingSize(ny, y_pad_factor);
  //int nz_pad = interval_simbox->GetNZpad();

  double true_x_pad_factor = static_cast<double>(nx_pad - nx)/static_cast<double>(nx);
  double true_y_pad_factor = static_cast<double>(ny_pad - ny)/static_cast<double>(ny);
  double true_z_pad_factor = interval_simbox->GetZPadFactor();
  double true_xPad    = true_x_pad_factor*lx;
  double true_yPad    = true_y_pad_factor*ly;
  //double true_zPad    = true_z_pad_factor*(interval_simbox->getlz()*interval_simbox->getMinRelThick());

  interval_simbox->SetNXpad(nx_pad);
  interval_simbox->SetNYpad(ny_pad);
  interval_simbox->SetXPadFactor(true_x_pad_factor);
  interval_simbox->SetYPadFactor(true_y_pad_factor);

  std::string text1;
  std::string text2;
  int logLevel = LogKit::Medium;
  if (model_settings->getEstimateXYPadding()) {
    text1 = " estimated from lateral correlation ranges in internal grid";
    logLevel = LogKit::Low;
  }
  if (model_settings->getEstimateZPadding()) {
    text2 = " estimated from an assumed wavelet length";
    logLevel = LogKit::Low;
  }

  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text1+":\n");
  LogKit::LogFormatted(logLevel,"  xPad, xPadFac, nx, nxPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_xPad, true_x_pad_factor, nx, nx_pad);
  LogKit::LogFormatted(logLevel,"  yPad, yPadFac, ny, nyPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_yPad, true_y_pad_factor, ny, ny_pad);
  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text2+":\n");
  //LogKit::LogFormatted(logLevel,"  zPad, zPadFac, nz, nzPad                 : %5.fms, %5.3f, %5d, %4d\n",
  //                     true_zPad, true_z_pad_factor, nz, nz_pad);
}

FFTGrid*
CommonData::CreateFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool file_grid)
{

  FFTGrid* fft_grid;

  if (file_grid)
    fft_grid =  new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fft_grid =  new FFTGrid(nx, ny, nz, nxp, nyp, nzp);

  return(fft_grid);

}

void CommonData::ReadAngularCorrelations(ModelSettings * model_settings,
                                         std::string   & err_text) {

  for (int t = 0; t < model_settings->getNumberOfTimeLapses(); t++) {

    Vario * vario                     = model_settings->getAngularCorr(t);
    int n_angles                      = model_settings->getNumberOfAngles(t);
    const std::vector<float> & angles = model_settings->getAngle(t);

    std::vector<std::vector<float> > angle_corr(n_angles);
    for (int i = 0; i < n_angles; i++) {
      angle_corr[i].resize(n_angles);
      for (int j = 0; j < n_angles; j++) {
        float d_angle = angles[i] - angles[j];
        angle_corr[i][i] = vario->corr(d_angle, 0);
      }
    }
    angular_correlations_.push_back(angle_corr);
  }
}


/*
void CommonData::WriteFilePriorVariances(const ModelSettings * model_settings,
const std::vector<float> & prior_corr_T,
const Surface * prior_corr_XY,
const float & dt) const
{
std::string baseName1 = IO::PrefixPrior() + IO::FileParameterCov() + IO::SuffixCrava();
std::string baseName2 = IO::PrefixPrior() + IO::FileTemporalCorr() + IO::SuffixCrava();
std::string baseName3 = IO::PrefixPrior() + IO::FileLateralCorr();
std::string fileName1 = IO::makeFullFileName(IO::PathToCorrelations(), baseName1);
std::string fileName2 = IO::makeFullFileName(IO::PathToCorrelations(), baseName2);

std::ofstream file;
NRLib::OpenWrite(file, fileName1);
file << std::fixed
<< std::right
<< std::setprecision(10);
for(int i=0 ; i<3 ; i++) {
for(int j=0 ; j<3 ; j++) {
file << std::setw(13) << prior_var0_(i,j) << " ";
}
file << "\n";
}
file.close();

NRLib::OpenWrite(file, fileName2);
file << std::fixed
<< std::right
<< std::setprecision(8)
<< dt << "\n";
for(int i=0 ; i<static_cast<int>(prior_corr_T.size()); i++) {
file << std::setw(11) << prior_corr_T[i] << "\n";
}
file.close();

IO::writeSurfaceToFile(*prior_corr_XY, baseName3, IO::PathToCorrelations(), model_settings->getOutputGridFormat());
}
*/

//--------------------------------------------------------------------
/*
void CommonData::PrintPriorVariances() const{

LogKit::LogFormatted(LogKit::Low,"\nVariances and correlations for parameter residuals:\n");
LogKit::LogFormatted(LogKit::Low,"\n");
LogKit::LogFormatted(LogKit::Low,"Variances ln Vp ln Vs ln Rho \n");
LogKit::LogFormatted(LogKit::Low,"---------------------------------------------------------------\n");
LogKit::LogFormatted(LogKit::Low,"Inversion grid: %.1e %.1e %.1e (used by program)\n",prior_var_0(0,0),prior_var_0(1,1),prior_var_0(2,2));

float corr01 = static_cast<float>(priorVar0_(0,1)/(sqrt(priorVar0_(0,0)*priorVar0_(1,1))));
float corr02 = static_cast<float>(priorVar0_(0,2)/(sqrt(priorVar0_(0,0)*priorVar0_(2,2))));
float corr12 = static_cast<float>(priorVar0_(1,2)/(sqrt(priorVar0_(1,1)*priorVar0_(2,2))));
LogKit::LogFormatted(LogKit::Low,"\n");
LogKit::LogFormatted(LogKit::Low,"Corr | ln Vp ln Vs ln Rho \n");
LogKit::LogFormatted(LogKit::Low,"-------+---------------------------\n");
LogKit::LogFormatted(LogKit::Low,"ln Vp | %5.2f %5.2f %5.2f \n",1.0f, corr01, corr02);
LogKit::LogFormatted(LogKit::Low,"ln Vs | %5.2f %5.2f \n",1.0f, corr12);
LogKit::LogFormatted(LogKit::Low,"ln Rho | %5.2f \n",1.0f);
LogKit::LogFormatted(LogKit::Low,"\n");

if (std::abs(corr01) > 1.0) {
LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vp-Vs correlation is wrong (%.2f).\n",corr01);
TaskList::addTask("Check your prior correlations. Corr(Vp,Vs) is out of bounds.");
}
if (std::abs(corr02) > 1.0) {
LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vp-Rho correlation is wrong (%.2f).\n",corr02);
TaskList::addTask("Check your prior correlations. Corr(Vp,Rho) is out of bounds.");
}
if (std::abs(corr12) > 1.0) {
LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The Vs-Rho correlation is wrong (%.2f).\n",corr12);
TaskList::addTask("Check your prior correlations. Corr(Vs,Rho) is out of bounds.");
}
}
*/