/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/commondata.h"
#include "src/simbox.h"
#include "src/timeline.h"
#include "src/fftgrid.h"
#include "src/seismicstorage.h"
#include "src/wavelet1D.h"
#include "src/multiintervalgrid.h"
#include "src/cravatrend.h"

#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/segy/segytrace.hpp"
#include "nrlib/well/norsarwell.hpp"


CommonData::CommonData(ModelSettings  * model_settings,
                       InputFiles     * input_files):
  outer_temp_simbox_(false),
  read_seismic_(false),
  read_wells_(false),
  block_wells_(false),
  setup_reflection_matrix_(false),
  optimize_well_location_(false),
  wavelet_estimation_shape_(false),
  prior_corr_estimation_(false),
  setup_estimation_rock_physics_(false),
  trend_cubes_(false),
  setup_multigrid_(false),
  multiple_interval_grid_(NULL)
{
  bool failed = false;
  std::string err_text = "";

  LogKit::WriteHeader("Common Data");

  //if(readSeismicData(model_settings,
  //                   input_files) == true)
  //  read_seismic_ = true; //True or false if there is no seismic data?

  // 1. set up outer simbox
  outer_temp_simbox_ = CreateOuterTemporarySimbox(model_settings, input_files, estimation_simbox_, full_inversion_volume_, err_text);
  failed = !outer_temp_simbox_;

  // 2. read seismic data
  //if(ReadSeismicData(model_settings, input_files) == true)
    //read_seismic_ = true; //True or false if there is no seismic data?

  // 3. read well data
  if(!failed){
    read_wells_ = ReadWellData(model_settings, &estimation_simbox_, input_files, log_names_, model_settings->getLogNames(),
                                model_settings->getInverseVelocity(), model_settings->getFaciesLogGiven(), err_text);
    failed = !read_wells_;
  }

  // 4. block wells for estimation
  // if well position is to be optimized or
  // if wavelet/noise should be estimated or
  // if correlations should be estimated

  if(!failed){
    if (model_settings->getOptimizeWellLocation() || model_settings->getEstimateWaveletNoise() || model_settings->getEstimateCorrelations()){
      block_wells_ = BlockWellsForEstimation(model_settings, estimation_simbox_, wells_, mapped_blocked_logs_, err_text);
      failed = !block_wells_;
    }
  }

  // 5. Reflection matrix and wavelet
  if(!failed){
    SetupReflectionMatrixAndTempWavelet(model_settings,
                                      input_files);
  }

  // 6. Optimization of well location
  if(!failed){
    if(model_settings->getOptimizeWellLocation()){
      OptimizeWellLocations(model_settings, input_files, &estimation_simbox_, wells_, mapped_blocked_logs_, seismic_data_, reflection_matrix_, err_text, failed);
      optimize_well_location_ = !failed;
    }
  }

  // 7. Setup of multiple interval grid
  if(!failed){
    multiple_interval_grid_ = new MultiIntervalGrid(model_settings, input_files, &estimation_simbox_, err_text, failed);
    setup_multigrid_ = !failed;
  }

  // 8. Trend Cubes
  if(setup_multigrid_ && model_settings->getFaciesProbFromRockPhysics() 
      && model_settings->getTrendCubeParameters().size() > 0){
    SetupTrendCubes(model_settings, input_files, multiple_interval_grid_, err_text, failed);
    setup_trend_cubes_ = !failed;
  }

  // 9. Rock Physics
  if(read_wells_ && setup_multigrid_ && model_settings->getFaciesProbFromRockPhysics()){
    if(model_settings->getTrendCubeParameters().size() > 0){ // If trends are used, the setup of trend cubes must be ok as well
      if(setup_trend_cubes_){
        SetupRockPhysics(model_settings, input_files, multiple_interval_grid_, trend_cubes_,
                         mapped_blocked_logs_, n_trend_cubes_, err_text,failed);
        setup_estimation_rock_physics_ = !failed;
      }
    }
    else{
      SetupRockPhysics(model_settings, input_files, multiple_interval_grid_, trend_cubes_,
                         mapped_blocked_logs_, n_trend_cubes_, err_text,failed);
      setup_estimation_rock_physics_ = !failed;
    }
  }
  
}

CommonData::~CommonData() {
  //delete estimation_simbox_;
  //delete full_inversion_volume_;
}

bool CommonData::CreateOuterTemporarySimbox(ModelSettings   * model_settings,
                                            InputFiles      * input_files,
                                            Simbox          & estimation_simbox,
                                            NRLib::Volume   & full_inversion_volume,
                                            std::string     & err_text){

  // parameters
  bool failed = false;
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

  if(model_settings->getForwardModeling())
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
  else if(area_specification == ModelSettings::AREA_FROM_SURFACE){
    LogKit::LogFormatted(LogKit::High,"\nFinding area information from surface \'"+input_files->getAreaSurfaceFile()+"\'\n");
    area_type = "Surface";
    RotatedSurface surf(input_files->getAreaSurfaceFile());
    SegyGeometry geometry(surf);
    model_settings->setAreaParameters(&geometry);
  }

  else if(area_specification == ModelSettings::AREA_FROM_GRID_DATA || // tilfelle iii
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
    if(geometry!=NULL){
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
            err_text += "Error: "+std::string(e.what());
            geometry->WriteILXL(true);
            geometry = NULL;
            failed = true;
          }
        }
        delete full_geometry;
      }
      else {
        geometry->WriteILXL();
      }
      if(!failed) {
        model_settings->setAreaParameters(geometry);
        ILXL_geometry = geometry;
      }
    }
    else {
      err_text += tmp_err_text;
      failed = true;
    }
  }

  // SET THE AREA FOR THE ESTIMATION SIMBOX

  if(!failed){
    const SegyGeometry * area_params = model_settings->getAreaParameters();
    failed = estimation_simbox.setArea(area_params, err_text);
    if(failed)
    {
      WriteAreas(area_params,&estimation_simbox_,area_type);
      err_text += "The specified AREA extends outside the surface(s).\n";
    }
  }

  // SET TOP AND BASE SURFACES FOR THE FORWARD/INVERSION INTERVALS -----------------------------------------------

  // if multiple intervals
  if(model_settings->getIntervalNames().size() > 0){
    SetSurfacesMultipleIntervals(model_settings, full_inversion_volume, estimation_simbox, input_files, err_text, failed);
  }
  // single interval described by either one or two surfaces
  else{
    SetSurfacesSingleInterval(model_settings, full_inversion_volume, estimation_simbox, input_files->getTimeSurfFiles(), err_text, failed);
  }

  return (!failed);
}


bool CommonData::ReadSeismicData(ModelSettings  * model_settings,
                                 InputFiles     * input_files)
{
  //Skip if there is no AVO-seismic.
  int timelaps_seismic_files = 0;

  for(size_t i = 0; i < input_files->getTimeLapseSeismicFiles().size(); i++) {
    timelaps_seismic_files += input_files->getNumberOfSeismicFiles(i);
  }

  if(timelaps_seismic_files == 0 && input_files->getSeismicFiles().size() == 0)
    return(false);

  //Skip if mode is estimation and wavelet/noise is not set for estimation.

  if(model_settings->getEstimationMode() == true && model_settings->getEstimateWaveletNoise() == false)
    return(false);


  const std::vector<std::vector<std::string> > seismic_timelapse_files = input_files->getTimeLapseSeismicFiles();

  int nTimeLapses = model_settings->getNumberOfTimeLapses();
  int error = 0;

  for(int thisTimeLapse = 0; thisTimeLapse < nTimeLapses; thisTimeLapse++) {

    if(input_files->getNumberOfSeismicFiles(thisTimeLapse) > 0 ) {

      LogKit::WriteHeader("Reading seismic data");


      std::vector<float> angles = model_settings->getAngle(thisTimeLapse);
      std::vector<float> offset = model_settings->getLocalSegyOffset(thisTimeLapse);

      int numberOfAngles = model_settings->getNumberOfAngles(thisTimeLapse);
      std::vector<SeismicStorage> seismic_data_angle;

      for (int i = 0 ; i < numberOfAngles ; i++) {

        std::string filename = input_files->getSeismicFile(thisTimeLapse,i);
        int fileType = IO::findGridType(filename);

        if(fileType == IO::SEGY) { //From ModelGeneral::readSegyFile

          SegY * segy = NULL;

          TraceHeaderFormat * format = model_settings->getTraceHeaderFormat(thisTimeLapse,i);

          if(format == NULL) { //Unknown format
            std::vector<TraceHeaderFormat*> traceHeaderFormats(0);

            if (model_settings->getTraceHeaderFormat() != NULL)
            {

              traceHeaderFormats.push_back(model_settings->getTraceHeaderFormat());
            }
            segy = new SegY(filename,
                            offset[i],
                            traceHeaderFormats,
                            true); // Add standard formats to format search
          }
          else //Known format, read directly.
            segy = new SegY(filename, offset[i], *format);


          float guard_zone = model_settings->getGuardZone();


          if(CheckThatDataCoverGrid(segy, offset[i], &estimation_simbox_, guard_zone) == true) { //Change this to full_inversion_volume_?
            float padding         = 2*guard_zone;
            bool  relativePadding = false;
            bool  onlyVolume      = true;


            segy->ReadAllTraces(&full_inversion_volume_, //timeCutSimBox  //estimation_simbox_ eller full_inversion_volume_?
                                padding,
                                onlyVolume,
                                relativePadding);
            segy->CreateRegularGrid();


            SeismicStorage seismicdata(filename, SeismicStorage::SEGY, angles[i], segy);
            seismic_data_angle.push_back(seismicdata);
          }
          else {
            LogKit::LogFormatted(LogKit::Warning, "Data from segy-file " + filename + " is not read.\n");
            error = 1;
          }


        } //SEGY
        else if(fileType == IO::STORM || fileType == IO::SGRI) { //From ModelGeneral::readStormFile
          StormContGrid * stormgrid = NULL;
          bool failed = false;

          try
          {
            stormgrid = new StormContGrid(0,0,0);
            stormgrid->ReadFromFile(filename);
          }
          catch (NRLib::Exception & e)
          {
            LogKit::LogFormatted(LogKit::Warning, "Error when reading storm-file " + filename +": " + NRLib::ToString(e.what()) + "\n");
            failed = true;
            error = 1;
          }

          if(failed == false) {
            SeismicStorage seismicdata_tmp;

            if(fileType == IO::STORM)

              seismicdata_tmp = SeismicStorage(filename, SeismicStorage::STORM, angles[i], stormgrid);
            else

              seismicdata_tmp = SeismicStorage(filename, SeismicStorage::SGRI, angles[i], stormgrid);


            seismic_data_angle.push_back(seismicdata_tmp);
          }

        } //STORM / SGRI
        else {

          LogKit::LogFormatted(LogKit::Warning, "Error when reading file " + filename +". File type not recognized.\n");
          error = 1;
        }
      } //nAngles

      seismic_data_[thisTimeLapse] = seismic_data_angle;

    }//ifSeismicFiles
  } //nTimeLapses

  if(error == 0)
    return true;
  else
    return false;
}

bool
CommonData::CheckThatDataCoverGrid(const SegY   * segy,
                                   float         offset,
                                   const Simbox * timeCutSimbox,
                                   float         guard_zone)
{
  // Seismic data coverage (translate to CRAVA grid by adding half a grid cell)
  float dz = segy->GetDz();
  float z0 = offset + 0.5f*dz;
  float zn = z0 + (segy->GetNz() - 1)*dz;

  // Top and base of interval of interest
  float top_grid = static_cast<float>(timeCutSimbox->getTopZMin());
  float bot_grid = static_cast<float>(timeCutSimbox->getBotZMax());

  //Use this instead? input for timeCutSimbox is estimation_simbox_
  //int nx = estimation_simbox_.getnx();
  //int ny = estimation_simbox_.getny();
  //float top_grid = static_cast<float>(full_inversion_volume_.GetTopZMin(nx, ny));
  //float bot_grid = static_cast<float>(full_inversion_volume_.GetBotZMax(nx, ny));

  // Find guard zone
  float top_guard = top_grid - guard_zone;
  float bot_guard = bot_grid + guard_zone;

  if (top_guard < z0) {
    float z0_new = z0 - ceil((z0 - top_guard)/dz)*dz;
    LogKit::LogFormatted(LogKit::Warning, "\nThere is not enough seismic data above the interval of interest. The seismic data\n"
                                          " must start at "+NRLib::ToString(z0_new)+"ms (in CRAVA grid) to allow for a "
                                          + NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n"
                                          + "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n"
                                          + "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n"
                                          + "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n"
                                          + "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n"
                                          + "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n"
                                          + "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n");
    return false;
  }
  if (bot_guard > zn) {
    float zn_new = zn + ceil((bot_guard - zn)/dz)*dz;
    LogKit::LogFormatted(LogKit::Warning, "\nThere is not enough seismic data below the interval of interest. The seismic data\n"
                                          "must end at "+NRLib::ToString(zn_new)+"ms (in CRAVA grid) to allow for a "
                                          + NRLib::ToString(guard_zone)+"ms FFT guard zone:\n\n"
                                          + "  Seismic data start (CRAVA grid) : "+NRLib::ToString(z0,1)+"\n"
                                          + "  Top of upper guard zone         : "+NRLib::ToString(top_guard,1)+"\n"
                                          + "  Top of interval-of-interest     : "+NRLib::ToString(top_grid,1)+"\n\n"
                                          + "  Base of interval-of-interest    : "+NRLib::ToString(bot_grid,1)+"\n"
                                          + "  Base of lower guard zone        : "+NRLib::ToString(bot_guard,1)+"\n"
                                          + "  Seismic data end (CRAVA grid)   : "+NRLib::ToString(zn,1)+"\n");
    return false;
  }

  return true;
}


bool CommonData::ReadWellData(ModelSettings                   * model_settings,
                              Simbox                          * estimation_simbox,
                              InputFiles                      * input_files,
                              std::vector<std::string>        & log_names,
                              const std::vector<std::string>  & log_names_from_user,
                              const std::vector<bool>         & inverse_velocity,
                              bool                              facies_log_given,
                              std::string                     & err_text){
  bool failed = false;

  for(size_t i = 0; i<log_names_from_user.size(); i++){
    log_names.push_back(log_names_from_user[i]);
  }

  int nWells = model_settings->getNumberOfWells();
  std::vector<float> dev_angle(nWells);
  try{
    if(nWells > 0)
      LogKit::WriteHeader("Reading wells");


    //std::vector<std::string> logNames = model_settings->getLogNames();
    //std::vector<bool> inverseVelocity = model_settings->getInverseVelocity();
    //bool faciesLogGiven = model_settings->getFaciesLogGiven();

    for(int i=0 ; i<nWells; i++) {

      std::string wellFileName = input_files->getWellFile(i);
      bool read_ok = false;
      NRLib::Well new_well(wellFileName, read_ok);

      if(wellFileName.find(".nwh",0) != std::string::npos)
        ProcessLogsNorsarWell(new_well, log_names, inverse_velocity, facies_log_given, err_text, failed);
      else if(wellFileName.find(".rms",0) != std::string::npos)
        ProcessLogsRMSWell(new_well, log_names, inverse_velocity, facies_log_given, err_text, failed);

      if(read_ok == true){
        wells_.push_back(new_well);
      }
      else{
        LogKit::LogFormatted(LogKit::Error, "Well format of file " + wellFileName + " not recognized.");
      }

      CalculateDeviation(new_well, model_settings, dev_angle[i], estimation_simbox, model_settings->getIndicatorWavelet(i));
    }
  }catch (NRLib::Exception & e) {
    err_text += "Error: " + NRLib::ToString(e.what());
    failed = true;
  }

  return !failed;
}

void CommonData::ProcessLogsNorsarWell(NRLib::Well                     & new_well,
                                       std::vector<std::string>        & log_names_from_user,
                                       const std::vector<bool>         & inverse_velocity,
                                       bool                              facies_log_given,
                                       std::string                     & error_text,
                                       bool                            & failed){
  const int factor_kilometer = 1000;

  if(log_names_from_user.size() == 0){
    log_names_from_user.push_back("TWT");
    log_names_from_user.push_back("DT");
    log_names_from_user.push_back("RHOB");
    log_names_from_user.push_back("DTS");
    log_names_from_user.push_back("FACIES");
  }

  for(size_t i = 0; i<log_names_from_user.size(); i++){
    // If the well does not contain the log specified by the user, return an error message
    if(!new_well.HasContLog(log_names_from_user[i]) && !new_well.HasDiscLog(log_names_from_user[i])){
      error_text+="Could not find log \'" + log_names_from_user[i] + "\' in well file \'"+new_well.GetWellName()+"\'.\n";
      failed = true;
    }
  }
  

  // Norsar wells must have a UTMX log
  if(new_well.HasContLog("UTMX")){
    std::vector<double> x_pos_temp = new_well.GetContLog("UTMX");
    for(unsigned int i=0; i<x_pos_temp.size(); i++){
      x_pos_temp[i] = x_pos_temp[i]*factor_kilometer;
    }
    new_well.AddContLog("X_pos", x_pos_temp);
    new_well.RemoveContLog("UTMX");
  }else{
    failed = true;
    error_text += "Could not find log 'UTMX' in well file "+new_well.GetWellName()+".\n";
  }

  // Norsar wells must have a UTMY log
  if(new_well.HasContLog("UTMY")){
    std::vector<double> y_pos_temp = new_well.GetContLog("UTMY");
    for(unsigned int i=0; i<y_pos_temp.size(); i++){
      y_pos_temp[i] = y_pos_temp[i]*factor_kilometer;
    }
    new_well.AddContLog("Y_pos", y_pos_temp);
    new_well.RemoveContLog("UTMY");
  }else{
    failed = true;
    error_text += "Could not find log 'UTMY' in well file "+new_well.GetWellName()+".\n";
  }

  // Norsar wells must have a TVD log
  if(new_well.HasContLog("TVD")){
    std::vector<double> tvd_temp = new_well.GetContLog("TVD");
    for(unsigned int i=0; i<tvd_temp.size(); i++){
      tvd_temp[i] = tvd_temp[i]*factor_kilometer;
    }
    new_well.RemoveContLog("TVD");
    new_well.AddContLog("TVD", tvd_temp);
  }else{ // Process MD log if TVD is not available?
    failed = true;
    error_text += "Could not find log 'TVD' in well file "+new_well.GetWellName()+".\n";
  }

  // Time is always entry 0 in the log name list and is always called TWT
  //if(new_well.HasContLog(log_names_from_user[0])){

  // Vp/Dt is always entry 1 in the log name list
  if(new_well.HasContLog(log_names_from_user[1])){
    if(inverse_velocity[0] == true){
      new_well.AddContLog("Dt", new_well.GetContLog(log_names_from_user[0]));
      new_well.RemoveContLog(log_names_from_user[0]);
    }
    else{
      new_well.AddContLog("Vp", new_well.GetContLog(log_names_from_user[0]));
      new_well.RemoveContLog(log_names_from_user[0]);
    }
  }

  // Rho is always entry 2 in the log name list
  if(new_well.HasContLog(log_names_from_user[2])){
    new_well.AddContLog("Rho", new_well.GetContLog(log_names_from_user[2]));
    new_well.RemoveContLog(log_names_from_user[2]);
  }

  // Vs/Dts is always entry 3 in the log name list
  if(new_well.HasContLog(log_names_from_user[3])){
    if(inverse_velocity[1] == true){
      new_well.AddContLog("Dts", new_well.GetContLog(log_names_from_user[3]));
      new_well.RemoveContLog(log_names_from_user[3]);
    }
    else{
      new_well.AddContLog("Vp", new_well.GetContLog(log_names_from_user[3]));
      new_well.RemoveContLog(log_names_from_user[3]);
    }
  }

  // If the facies log is given, it is always entry 4 in the log name list
  if(facies_log_given){
    if(new_well.HasDiscLog(log_names_from_user[4])){
      new_well.AddDiscLog("Facies", new_well.GetDiscLog(log_names_from_user[4]));
      new_well.RemoveDiscLog(log_names_from_user[4]);
    }
  }


}

void CommonData::ProcessLogsRMSWell(NRLib::Well                     & new_well,
                                    std::vector<std::string>        & log_names_from_user,
                                    const std::vector<bool>         & inverse_velocity,
                                    bool                              facies_log_given,
                                    std::string                     & error_text,
                                    bool                            & failed){

  const double factor_usfeet_to_meters = 304800.0;

  if(log_names_from_user.size() == 0){
    log_names_from_user.push_back("TWT");
    log_names_from_user.push_back("DT");
    log_names_from_user.push_back("RHOB");
    log_names_from_user.push_back("DTS");
    log_names_from_user.push_back("FACIES");
  }

  for(size_t i = 0; i<log_names_from_user.size(); i++){
    // If the well does not contain the log specified by the user, return an error message
    if(!new_well.HasContLog(log_names_from_user[i]) && !new_well.HasDiscLog(log_names_from_user[i])){
      error_text+="Could not find log \'" + log_names_from_user[i] + "\' in well file \'"+new_well.GetWellName()+"\'.\n";
      failed = true;
    }
  }

  // RMS wells must have an x log
  if(new_well.HasContLog("x")){
    new_well.AddContLog("X_pos", new_well.GetContLog("x"));
    new_well.RemoveContLog("x");
  }else{
    failed = true;
    error_text += "Could not find log 'x' in well file "+new_well.GetWellName()+".\n";
  }

  // RMS wells must have a y log
  if(new_well.HasContLog("y")){
    new_well.AddContLog("Y_pos", new_well.GetContLog("y"));
    new_well.RemoveContLog("y");
  }else{
    failed = true;
    error_text += "Could not find log 'y' in well file "+new_well.GetWellName()+".\n";
  }

  // RMS wells must have a z log
  if(new_well.HasContLog("z")){
    new_well.AddContLog("TVD", new_well.GetContLog("z"));
    new_well.RemoveContLog("z");
  }else{
    failed = true;
    error_text += "Could not find log 'z' in well file "+new_well.GetWellName()+".\n";
  }

  // Time is always entry 0 in the log name list
  // Time is always called 'TWT', so no need to rename it
  //if(new_well.HasContLog(log_names_from_user[0]))

  // Vp is always entry 1 in the log name list and entry 0 in inverse_velocity
  if(new_well.HasContLog(log_names_from_user[1])){
    std::vector<double> vp_temp = new_well.GetContLog(log_names_from_user[1]);
    std::vector<double> vp(vp_temp.size());
    if(inverse_velocity[0]){
      for(unsigned int i=0; i<vp_temp.size(); i++){
        vp[i] = static_cast<double>(factor_usfeet_to_meters/vp_temp[i]);
      }
    }
    else{
      for(unsigned int i=0; i<vp_temp.size(); i++){
        vp[i] = static_cast<double>(vp_temp[i]);
      }
    }
    new_well.RemoveContLog(log_names_from_user[1]);
    new_well.AddContLog("Vp", vp);
  }

  // Vs is always entry 3 in the log name list
  if(new_well.HasContLog(log_names_from_user[3])){
    std::vector<double> vs_temp = new_well.GetContLog(log_names_from_user[3]);
    std::vector<double> vs(vs_temp.size());
    if(inverse_velocity[1]){
      for(unsigned int i=0; i<vs_temp.size(); i++){
        vs[i] = static_cast<double>(factor_usfeet_to_meters/vs_temp[i]);
      }
    }
    else{
      for(unsigned int i=0; i<vs_temp.size(); i++){
        vs[i] = static_cast<double>(vs_temp[i]);
      }
    }
    new_well.RemoveContLog(log_names_from_user[1]);
    new_well.AddContLog("Vs", vs);
  }

  // Rho is always entry 2 in the log name list
  if(new_well.HasContLog(log_names_from_user[2])){
    new_well.AddContLog("Rho", new_well.GetContLog(log_names_from_user[2]));
    new_well.RemoveContLog(log_names_from_user[3]);
  }

  // If defined, Facies is always entry 4 in the log name list
  if(facies_log_given){
    if(new_well.HasDiscLog(log_names_from_user[4])){
      new_well.AddDiscLog("Facies", new_well.GetDiscLog(log_names_from_user[4]));
      new_well.RemoveDiscLog(log_names_from_user[4]);
    }
  }
}

//void
//CommonData::readNorsarWell(const std::string              & wellFileName,
//                           NRLib::Well                    & new_well,
//                           const std::vector<std::string> & logNames,
//                           const std::vector<bool>        & inverseVelocity,
//                           bool                             faciesLogGiven,
//                           std::string                    & error)
//{
//  std::string name = NRLib::RemovePath(wellFileName);
//  name = NRLib::ReplaceExtension(name, "");
//  std::string well_name = name;
//  int faciesok = 0; // all faciesnumbers read are present in header //1 ok
//  int timemissing = 0;
//  double xpos0 = 0.0;
//  double ypos0 = 0;
//  int nFacies = 0;
//
//  std::vector<double> xpos;
//  std::vector<double> ypos;
//  std::vector<double> zpos;
//  std::vector<double> alpha;
//  std::vector<double> beta;
//  std::vector<double> rho;
//  std::vector<double> md;
//  std::vector<int> facies; // Always allocate a facies log (for code simplicity). Keep this ??
//
//  try
//  {
//    NRLib::NorsarWell well(wellFileName);
//
//    int nVar = 5;       // z,alpha,beta,rho, and facies
//
//    std::vector<std::string> parameterList;
//
//    bool vpLog = false;
//    bool vsLog = false;
//
//    if(logNames[0] != "") // Assume that all lognames are filled present if first is.
//    {
//      parameterList = logNames;
//      if (!faciesLogGiven)
//        nVar = 4;
//      vpLog = !inverseVelocity[0];
//      vsLog = !inverseVelocity[1];
//    }
//    else
//    {
//      parameterList[0] = "TWT";
//      parameterList[1] = "DT";
//      parameterList[2] = "RHOB";
//      parameterList[3] = "DTS";
//      parameterList[4] = "FACIES";
//    }
//
//    int nLogs  = 2+nVar;
//    int nExtra = 1; //MD log, needed for writing.
//    std::vector<double> * filler = NULL; //to eliminate warning.
//    std::vector<std::vector<double> *> logs(nLogs+nExtra, filler);
//    if(well.HasContLog("UTMX") == false) {
//      //error = 1;
//      //LogKit::LogFormatted(LogKit::Error,"Could not find log 'UTMX' in well file "+wellFileName+".\n");
//      error += "Could not find log 'UTMX' in well file "+wellFileName+".\n";
//      logs[0] = NULL;
//    }
//    else
//      logs[0] = &(well.GetContLog("UTMX"));
//    if(well.HasContLog("UTMY") == false) {
//      //error = 1;
//      //LogKit::LogFormatted(LogKit::Error,"Could not find log 'UTMY' in well file "+wellFileName+".\n");
//      error += "Could not find log 'UTMY' in well file "+wellFileName+".\n";
//      logs[1] = NULL;
//    }
//    logs[1] = &(well.GetContLog("UTMY"));
//    for(int i=0;i<nVar;i++) {
//      if(well.HasContLog(parameterList[i]) == false) {
//        logs[2+i] = NULL;
//        if(i != 4 || logNames[0] != "") {
//          //error = 1;
//          //LogKit::LogFormatted(LogKit::Error,"Could not find log "+parameterList[i]+" in well file "+wellFileName+".\n");
//          error += "Could not find log "+parameterList[i]+" in well file "+wellFileName+".\n";
//        }
//        else if(i==4)
//          nLogs = nLogs-1;
//      }
//      else
//        logs[2+i] = &(well.GetContLog(parameterList[i]));
//    }
//
//    //Added MD log.
//    int mdLog = nLogs;
//    if(well.HasContLog("MD") == false) {
//      //error = 1;
//      //LogKit::LogFormatted(LogKit::Error,"Could not find log 'MD' in well file "+wellFileName+".\n");
//      error += "Could not find log 'MD' in well file "+wellFileName+".\n";
//      logs[mdLog] = NULL;
//    }
//    logs[mdLog] = &(well.GetContLog("MD"));
//    nLogs++;
//
//    if(logs[2] == NULL)
//      timemissing = 1;
//    else
//      timemissing =0;
//
//    if(error == "") {
//      faciesok = 1;
//      std::vector<int> facCodes;
//      int nd = 0;
//      for(size_t i=0;i<logs[2]->size();i++)
//        if(well.IsMissing((*logs[2])[i]) == false)
//          nd++;
//
//      for(size_t i=0;i<logs[0]->size();i++) {
//        if(well.IsMissing((*logs[2])[i]) == false) {
//          xpos.push_back((*logs[0])[i]*1000);
//          ypos.push_back((*logs[1])[i]*1000);
//          zpos.push_back((*logs[2])[i]*1000);
//          if(!well.IsMissing((*logs[3])[i]))
//            alpha.push_back(static_cast<double>((*logs[3])[i]));
//          else
//            alpha.push_back(RMISSING);
//          if(!well.IsMissing((*logs[5])[i]))
//            beta.push_back(static_cast<float>((*logs[5])[i]));
//          else
//            beta.push_back(RMISSING);
//          if(!well.IsMissing((*logs[4])[i]))
//            rho.push_back(static_cast<float>((*logs[4])[i]));
//          else
//            rho.push_back(RMISSING);
//          if(mdLog != 6 && nLogs > 6 && logs[6] != NULL && !well.IsMissing((*logs[6])[i])) {
//            facies.push_back(static_cast<int>((*logs[6])[i]));
//
//            if(find(facCodes.begin(), facCodes.end(), facies[facies.size()]) == facCodes.end())
//              facCodes.push_back(facies[facies.size()]);
//          }
//          else
//            facies.push_back(IMISSING);
//          if(!well.IsMissing((*logs[mdLog])[i]))
//            md.push_back(static_cast<float>((*logs[mdLog])[i]));
//          else
//            md.push_back(RMISSING);
//        }
//      }
//      nFacies = static_cast<int>(facCodes.size());
//
//    }
//    xpos0 = well.GetXPosOrigin()*1000;
//    ypos0 = well.GetYPosOrigin()*1000;
//
//  }
//  catch (NRLib::Exception & e) {
//    //LogKit::LogFormatted(LogKit::Error,"Error: " + NRLib::ToString(e.what()));
//    //error = 1;
//    error += "Error: " + NRLib::ToString(e.what());
//  }
//
//  //Add well logs
//  new_well.SetWellName(well_name);
//  new_well.AddContLog("xpos", xpos);
//  new_well.AddContLog("ypos", ypos);
//  new_well.AddContLog("zpos", zpos);
//  new_well.AddContLog("alpha", alpha);
//  new_well.AddContLog("beta", beta);
//  new_well.AddContLog("rho", rho);
//  if(md.size() > 0)
//    new_well.AddContLog("md", md);
//
//  new_well.AddDiscLog("facies", facies);
//  new_well.SetMissing(RMISSING);
//
//  timemissing_[well_name] = timemissing;
//  xpos0_[well_name] = xpos0;
//  ypos0_[well_name] = ypos0;
//  faciesok_[well_name] = faciesok;
//  nFacies_[well_name] = nFacies;
//
//}


//void
//CommonData::readRMSWell(const std::string              & wellFileName,
//                        NRLib::Well                    & new_well,
//                        const std::vector<std::string> & logNames,
//                        const std::vector<bool>        & inverseVelocity,
//                        bool                             faciesLogGiven,
//                        std::string                    & error)
//{
//  double dummy = RMISSING;
//  int j,k;
//  //int timemissing = 0;
//  std::ifstream file;
//  NRLib::OpenRead(file, wellFileName);
//  std::string token, dummyStr;
//  std::vector<std::string> tokenLine;
//
//  int timemissing;
//  //From WellData::ReadRMSWell
//
//  std::vector<std::string>  faciesNames;
//  std::vector<double> xpos;
//  std::vector<double> ypos;
//  std::vector<double> zpos;
//  std::vector<double> alpha;
//  std::vector<double> beta;
//  std::vector<double> rho;
//  std::vector<int> facies; // Always allocate a facies log (for code simplicity). Keep this ??
//
//  int nlog; // number of logs in file
//  int line = 0;
//  NRLib::DiscardRestOfLine(file,line,false); //First two lines contain info we do not need.
//  NRLib::DiscardRestOfLine(file,line,false);
//  NRLib::ReadNextToken(file, token, line);
//  std::string well_name = token;
//  double xpos0 = NRLib::ReadNext<double>(file, line); //Needed?
//  double ypos0 = NRLib::ReadNext<double>(file, line); //Needed
//  NRLib::DiscardRestOfLine(file,line,false);
//  nlog   = NRLib::ReadNext<int>(file, line);
//
//  //int faciesok = 1; // all faciesnumbers read are present in header
//
//  //Start searching for key words.
//
//  int nVar = 5;       // z,alpha,beta,rho, and facies
//
//  std::vector<std::string> parameterList(5);
//
//  bool vpLog = false;
//  bool vsLog = false;
//
//  if(logNames[0] != "") // Assume that all lognames are filled present if first is.
//  {
//    parameterList = logNames;
//    if (!faciesLogGiven)
//      nVar = 4;
//    vpLog = !inverseVelocity[0];
//    vsLog = !inverseVelocity[1];
//  }
//  else
//  {
//    parameterList[0] = "TWT";
//    parameterList[1] = "DT";
//    parameterList[2] = "RHOB"; //DTS before RHOB?
//    parameterList[3] = "DTS";
//    parameterList[4] = "FACIES";
//  }
//
//  int * pos = new int[nVar];
//  for(k=0;k<nVar;k++)
//    pos[k] = IMISSING;
//
//  int nFacies = 0;
//  std::string faciesLogName;
//  for(k=0;k<nlog;k++)
//  {
//    NRLib::ReadNextToken(file,token,line);
//    for(j=0;j<nVar;j++)
//    {
//      if( NRLib::Uppercase(token)==parameterList[j])
//      {
//        pos[j] = k + 4;
//        if(j==4)
//        {
//          faciesLogName = parameterList[4];
//          // facies log - save names
//          NRLib::ReadNextToken(file,token,line); // read code word DISC
//          if (token != "DISC")
//          {
//            LogKit::LogFormatted(LogKit::Error,"ERROR: Facies log must be discrete.\n");
//            exit(1);
//          }
//          // Find number of facies
//          std::getline(file,dummyStr);
//          tokenLine = NRLib::GetTokens(dummyStr);
//          nFacies = static_cast<int>(tokenLine.size())/2;
//        }
//      }
//    }
//    if (token != "DISC")
//      NRLib::DiscardRestOfLine(file,line,false);
//  }
//
//  std::string missVar = "";
//  bool missing_error = false;
//  for(k=0 ; k<nVar ; k++)
//  {
//    if(pos[k]==IMISSING)
//    {
//      missVar += parameterList[k];
//      //error = 1;
//      missing_error = true;
//    }
//  }
//  if(missing_error == true)
//    error += "Cannot find log(s) "+missVar+" in well file "+wellFileName+".\n";
//    //LogKit::LogFormatted(LogKit::Error,"Cannot find log(s) "+missVar+" in well file "+wellFileName+".\n");
//
//  if(pos[0]==IMISSING)
//    timemissing = 1;
//  else
//    timemissing = 0;
//
//  // Find nd, the number of observations in well.
//  // Count the number of time observations which is not missing values.
//
//  int nData = 0;
//  int legalData = 0;
//  int nd;
//  while (NRLib::CheckEndOfFile(file)==false)
//  {
//    nData++;
//    try {
//      dummy = NRLib::ReadNext<double>(file,line); // Read x which we do not need yet.
//      dummy = NRLib::ReadNext<double>(file,line); // Read y which we do not need yet.
//      dummy = NRLib::ReadNext<double>(file,line); // Read z which we do not need.
//
//      for(j=4;j<=nlog+3;j++) {
//        dummy = NRLib::ReadNext<double>(file,line); // Read z which we do not need.
//        if(j==pos[0] && dummy != WELLMISSING) {
//          legalData++;   // Found legal TIME variable
//        }
//      }
//    }
//    catch (NRLib::IOError e) {
//      std::string text;
//      text += std::string("\nERROR: Reading of well \'") + wellFileName + "\' failed for log record ";
//      text += NRLib::ToString(nData) + " (not counting header lines).\n";
//      text += std::string("\nERROR message is \'") + e.what() + "\'";
//      LogKit::LogMessage(LogKit::Error,text);
//      exit(1);
//    }
//  }
//  file.close();
//  file.clear();
//  nd = legalData;
//
//  //
//  // Check that the number of logs found for each log entry agrees
//  // with the number of logs specified in header.
//  //
//  // A nicer and faster implementation for this is requested...
//  //
//
//  int logEntry = 0;
//  NRLib::OpenRead(file, wellFileName);
//  line = 0;
//  for(k=0;k<4+nlog;k++)
//    NRLib::DiscardRestOfLine(file,line,false);
//  while (NRLib::CheckEndOfFile(file)==false && error!="")
//  {
//    logEntry++;
//    int  elements = 0;
//    bool lastIsBlank = true;
//    std::getline(file,token);
//
//    int l = static_cast<int>(token.length());
//    for (k=0; k<l; k++)
//    {
//      if (token[k] != ' ' && token[k] != '\t' && token[k] != '\r' && token[k] != '\0')
//      {
//        if (lastIsBlank)
//          elements++;
//        lastIsBlank = false;
//      }
//      else
//        lastIsBlank = true;
//    }
//    if(elements != nlog+3)
//    {
//      //error = 1;
//      //LogKit::LogFormatted(LogKit::Error,"ERROR for well "+well_name +": The number of log elements (nlogs="+NRLib::ToString(elements-3)+")\n in line "+NRLib::ToString(logEntry)+" does not match header specifications (nlogs="+NRLib::ToString(nlog)+").\n");
//      error += "ERROR for well "+well_name +": The number of log elements (nlogs="+NRLib::ToString(elements-3)+")\n in line "+NRLib::ToString(logEntry)+" does not match header specifications (nlogs="+NRLib::ToString(nlog)+").\n";
//    }
//  }
//  file.close();
//  file.clear();
//
//  //
//  // Read logs
//  //
//  std::vector<int> faciesNr;
//  //if (nFacies > 0)
//  //  faciesNr.resize(nFacies);
//
//  NRLib::OpenRead(file, wellFileName);
//  line = 0;
//  for(k=0;k<4+nlog;k++)
//  {
//    NRLib::ReadNextToken(file,token,line);
//    if (NRLib::Uppercase(token) == parameterList[4])
//    {
//      NRLib::ReadNextToken(file,token,line); // read code word DISC
//      // facies types given here
//      for(int kk=0;kk<nFacies;kk++)
//      {
//        NRLib::ReadNextToken(file,token,line);
//        //faciesNr[kk] = NRLib::ParseType<int>(token);
//        faciesNr.push_back(NRLib::ParseType<int>(token));
//        NRLib::ReadNextToken(file,token,line);
//        faciesNames.push_back(token);
//      }
//    }
//    NRLib::DiscardRestOfLine(file,line,false);
//  }
//  double OPENWORKS_MISSING = -999.25;
//  bool wrongMissingValues = false;
//
//  double xpos_data = RMISSING;
//  double ypos_data = RMISSING;
//  double zpos_data = RMISSING;
//  double alpha_data = RMISSING;
//  double beta_data = RMISSING;
//  double rho_data = RMISSING;
//  int facies_data = IMISSING;
//  k         = -1;
//  int faciesok = 1;
//  int legal = 0;
//  for(k=0;k<nData;k++)
//  {
//    xpos_data  = NRLib::ReadNext<double>(file,line);
//    ypos_data  = NRLib::ReadNext<double>(file,line);
//    dummy = NRLib::ReadNext<double>(file,line);
//    for(j=4;j<=nlog+3;j++)
//    {
//      dummy = NRLib::ReadNext<double>(file,line);
//      if(j==pos[0])
//      {
//        //Found TIME variable
//        if(dummy != WELLMISSING && dummy != OPENWORKS_MISSING)
//        {
//          zpos_data = dummy;
//          k++;
//          legal = 1;
//        }
//        else
//          legal = 0;
//      }
//      else if(j==pos[1])
//      {
//        // Found ALPHA variable
//        if(dummy != WELLMISSING)
//          if (vpLog)
//            alpha_data = static_cast<double>(dummy);
//          else
//            alpha_data = static_cast<double>(304800.0/dummy);
//        else
//          alpha_data = RMISSING;
//      }
//      else if(j==pos[3])
//      {
//        // Found BETA variable
//        if(dummy != WELLMISSING)
//          if (vsLog)
//            beta_data = static_cast<double>(dummy);
//          else
//            beta_data = static_cast<double>(304800.0/dummy);
//        else
//          beta_data = RMISSING;
//      }
//      else if(j==pos[2])
//      {
//        //Found RHO variable
//        if(dummy != WELLMISSING)
//          rho_data = static_cast<double>(dummy);
//        else
//          rho_data = RMISSING;
//      }
//      else if(nVar > 4 && j==pos[4])
//      {
//        //Found facies variable
//        if(dummy != WELLMISSING)
//          facies_data = static_cast<int>(dummy);
//        else
//          facies_data = IMISSING;
//        if(facies_data!=IMISSING)
//        {
//          int faciesok_check = 0;
//          for(int kk=0 ; kk<nFacies ; kk++)
//          {
//            if(facies_data == faciesNr[kk]) {
//              faciesok_check = 1;
//              break;
//            }
//          }
//          if(faciesok_check == 0)
//            faciesok = 0;
//        }
//      }
//    }
//    if(legal == 1)
//    {
//      //Cut against full_inversion_volum
//      double z_top = full_inversion_volume_.GetTopSurface().GetZ(xpos_data, ypos_data);
//      double z_bot = full_inversion_volume_.GetBotSurface().GetZ(xpos_data, ypos_data);
//
//      if(zpos_data < z_top && zpos_data > z_bot) {
//        xpos.push_back(xpos_data);
//        ypos.push_back(ypos_data);
//        zpos.push_back(zpos_data);
//        alpha.push_back(alpha_data);
//        beta.push_back(beta_data);
//        rho.push_back(rho_data);
//        facies.push_back(facies_data);
//        if (alpha_data == OPENWORKS_MISSING)
//          wrongMissingValues = true;
//        if(beta_data == OPENWORKS_MISSING)
//          wrongMissingValues = true;
//        if(rho_data == OPENWORKS_MISSING)
//          wrongMissingValues = true;
//      }
//    }
//  }
//  file.close();
//
//  //Add well logs
//  new_well.SetWellName(well_name);
//  new_well.AddContLog("xpos", xpos);
//  new_well.AddContLog("ypos", ypos);
//  new_well.AddContLog("zpos", zpos);
//  new_well.AddContLog("alpha", alpha);
//  new_well.AddContLog("beta", beta);
//  new_well.AddContLog("rho", rho);
//
//  new_well.AddDiscLog("facies", facies);
//  new_well.SetMissing(RMISSING);
//
//  timemissing_[well_name] = timemissing;
//  xpos0_[well_name] = xpos0;
//  ypos0_[well_name] = ypos0;
//  faciesnames_[well_name] = faciesNames;
//  faciesok_[well_name] = faciesok;
//  faciesNr_[well_name] = faciesNr;
//  nFacies_[well_name] = nFacies;
//}



bool CommonData::SetupReflectionMatrixAndTempWavelet(ModelSettings  * model_settings,
                                                     InputFiles     * input_files) {
  LogKit::WriteHeader("Reflection matrix");
  //
  // About to process wavelets and energy information. Needs the a-matrix, so create
  // if not already made. A-matrix may need Vp/Vs-ratio from background model or wells.
  //
  const std::string & reflMatrFile = input_files->getReflMatrFile();
  const double        vpvs         = model_settings->getVpVsRatio();
  float                  ** reflectionMatrix;


  unsigned int nTimeLapses = model_settings->getNumberOfTimeLapses(); //Returnerer timeLapseAngle_.size()
  for(size_t thisTimeLapse = 0; thisTimeLapse < nTimeLapses; thisTimeLapse++) {


    std::vector<float> angles = model_settings->getAngle(thisTimeLapse);
    std::vector<float> offset = model_settings->getLocalSegyOffset(thisTimeLapse);


    int numberOfAngles = model_settings->getNumberOfAngles(thisTimeLapse);


    if (reflMatrFile != "") {  //File should have one line for each seismic data file. Check: if(input_files->getNumberOfSeismicFiles(thisTimeLapse) > 0 ) ?
      std::string tmpErrText("");
      reflectionMatrix = ReadMatrix(reflMatrFile, numberOfAngles, 3, "reflection matrix", tmpErrText);
      if(reflectionMatrix == NULL) {
        LogKit::LogFormatted(LogKit::Error, "Reading of file "+reflMatrFile+ " for reflection matrix failed\n");
        LogKit::LogFormatted(LogKit::Error, tmpErrText);
        return false;
      }

      LogKit::LogFormatted(LogKit::Low,"\nReflection parameters read from file.\n\n");
    }

    else if(vpvs != RMISSING) {
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs ratio specified in model file.\n");
      double vsvp = 1.0/vpvs;
      SetupDefaultReflectionMatrix(reflectionMatrix, vsvp, model_settings, numberOfAngles, thisTimeLapse);
    }
    else {
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs equal to 2\n");
      double vsvp = 1/2;
      SetupDefaultReflectionMatrix(reflectionMatrix, vsvp, model_settings, numberOfAngles, thisTimeLapse);
    }


      reflection_matrix_[thisTimeLapse] = reflectionMatrix;
  } //nTimeLapses


  setup_reflection_matrix_ = true;


  //Set up temporary wavelet
  LogKit::WriteHeader("Setting up temporary wavelet");


  //1. Check if optimize welllocation //Point f) Comes from xml-model file
    //2. Check if read seismic ok | read_seismic_ok_
  read_seismic_ = true;
  if(model_settings->getOptimizeWellLocation() == false && read_seismic_ == true) {

    //3. Use Ricker - wavelet.
    //4. 1 wavelet per angle
    //5  Frequency per ange: Take 100 traces from first AVO-vintage on this angle. Find peak-frequency for these.
    int numberOfAngles = model_settings->getNumberOfAngles(0);
    int error = 0;

    for (int i = 0 ; i < numberOfAngles ; i++) {
      //Check all timelapses for this angle, choose the first one;
      int thisTimeLapse = 0;
      int vintageyear = model_settings->getVintageYear(0);
      int vintagemonth = model_settings->getVintageMonth(0);
      int vintageday = model_settings->getVintageDay(0);
      for(unsigned int j = 1; j < nTimeLapses; j++) {
        if(model_settings->getVintageYear(j) <= vintageyear && model_settings->getVintageMonth(j) <= vintagemonth && model_settings->getVintageDay(j) <= vintageday) {
          vintageyear = model_settings->getVintageYear(j);
          vintagemonth = model_settings->getVintageMonth(j);
          vintageday = model_settings->getVintageDay(j);
          thisTimeLapse = j;
        }
      }
      std::vector<float> angles = model_settings->getAngle(thisTimeLapse);

      //int tmp_type = seismic_data_[thisTimeLapse][i].getSeismicType();
      int n_traces;
      std::vector<std::vector<float> > trace_data(100);
      std::vector<double> trace_length;
      std::vector<float> frequency_peaks;

      if(seismic_data_[thisTimeLapse][i].GetSeismicType() == SeismicStorage::SEGY) {
        SegY * segy = seismic_data_[thisTimeLapse][i].GetSegY();
        n_traces = segy->GetNTraces();
        trace_length.resize(100, 1); //Dummy, to avoid checking later if trace_data came from SEGY or STORM

        for(size_t j = 0; j < 100; j++) {
          int trace_index = j*(static_cast<int>(n_traces / 100));
          NRLib::SegYTrace * segy_tmp = segy->getTrace(trace_index);

          if(segy_tmp != NULL) {
            size_t start = segy_tmp->GetStart();
            size_t end = segy_tmp->GetEnd();
            for(size_t k = start; k < end; k++) {
              trace_data[j].push_back(segy_tmp->GetValue(k));
            }
          }
        }
      }
      else {
        StormContGrid * stormg = seismic_data_[thisTimeLapse][i].GetStorm();

        double x_tmp = 0.0;
        double y_tmp = 0.0;
        double z_tmp = 0.0;

        unsigned int index_i = 0;
        unsigned int index_j = 0;
        int trace_index = 0;

        for(int ii = 0; ii < 10; ii++) {

          index_i = ii*static_cast<int>(stormg->GetNI()/10);
          if(index_i >= stormg->GetNI())
            index_i = stormg->GetNI() -1;

          for(int jj = 0; jj < 10; jj++) {

            index_j = jj*static_cast<int>(stormg->GetNJ()/10);
            if(index_j >= stormg->GetNJ())
              index_j = stormg->GetNJ()-1;

            for(size_t kk = 0; kk < stormg->GetNK(); kk++) {
              stormg->FindCenterOfCell(index_i, index_j, kk, x_tmp, y_tmp, z_tmp);
              trace_data[trace_index].push_back(stormg->GetValueClosestInZ(x_tmp, y_tmp, z_tmp));
            }

            //Store length
            double top = stormg->GetTopSurface().GetZ(x_tmp, y_tmp);
            double bot = stormg->GetBotSurface().GetZ(x_tmp, y_tmp);
            trace_length.push_back(std::abs(bot-top));

            trace_index++;
          }
        }
      }

      //FFT to find peak-frequency.
      for(int j = 0; j < 100; j++) {
        int n_trace = trace_data[j].size();
        fftw_real    * seis_r = new fftw_real;
        fftw_complex * seis_c = reinterpret_cast<fftw_complex*>(seis_r);
        seis_r = new fftw_real[n_trace];

        for(int k=0; k < n_trace; k++)
          seis_r[k] = trace_data[j][k];

        Utils::fft(seis_r, seis_c, n_trace);

        //std::vector<float> seis_tmp;

        float peak_tmp = 0.0;
        for(int k = 0; k < n_trace; k++) {
          //seis_tmp.push_back(seis_r[k]);
          if(seis_r[k] > peak_tmp)
            peak_tmp = seis_r[k];
        }

        peak_tmp *= static_cast<float>(trace_length[j]);
        frequency_peaks.push_back(peak_tmp);
      }

      float mean_frequency = 0.0;
      float mean_grid_height = 0.0;
      for(size_t j = 0; j < frequency_peaks.size(); j++) {
        mean_frequency += frequency_peaks[j];
        mean_grid_height += static_cast<float>(trace_length[j]);
      }
      mean_frequency /= frequency_peaks.size();
      mean_grid_height /= frequency_peaks.size(); //Will be 1 if SEGY is used.

      mean_frequency /= mean_grid_height;

      int tmp_error = 0;
      Wavelet  * wavelet_tmp = new Wavelet1D(model_settings, reflectionMatrix[i], angles[i], mean_frequency, tmp_error);

      error += tmp_error;

      if(tmp_error == 0)
        temporary_wavelets_.push_back(wavelet_tmp);
      else
        LogKit::LogFormatted(LogKit::Error, "Error setting up a temporary wavelet for angle " + NRLib::ToString(angles[i]) + ".\n");
    }

    if(error == 0)
      temporary_wavelet_ = true;
  }


  return true;


}


float **  CommonData::ReadMatrix(const std::string & fileName, int n1, int n2,
                                 const std::string & readReason,
                                 std::string       & errText)
{
  float * tmpRes = new float[n1*n2+1];
  std::ifstream inFile;
  NRLib::OpenRead(inFile,fileName);
  std::string text = "Reading "+readReason+" from file "+fileName+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
  std::string storage;
  int index = 0;
  int error = 0;


  while(error == 0 && inFile >> storage) {
    if(index < n1*n2) {
      try {
        tmpRes[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        errText += "Error in "+fileName+"\n";
        errText += e.what();
        error = 1;
      }
    }
    index++;
  }
  if(error == 0) {
    if(index != n1*n2) {
      error = 1;
      errText += "Found "+NRLib::ToString(index)+" in file "+fileName+", expected "+NRLib::ToString(n1*n2)+".\n";
    }
  }


  float ** result = NULL;
  if(error == 0) {
    LogKit::LogFormatted(LogKit::Low,"ok.\n");
    result = new float * [n1];
    int i, j;
    index = 0;
    for(i=0;i<n1;i++) {
      result[i] = new float[n2];
      for(j=0;j<n2;j++) {
        result[i][j] = tmpRes[index];
        index++;
      }
    }
  }
  else

    LogKit::LogFormatted(LogKit::Low,"failed.\n");
  delete [] tmpRes;
  return(result);
}


void
CommonData::SetupDefaultReflectionMatrix(float             **& reflectionMatrix,
                                         double                vsvp,
                                         const ModelSettings * model_settings,
                                         int                   numberOfAngles,
                                         int                   thisTimeLapse)
{
  int      i;
  float ** A      = new float * [numberOfAngles];


  double           vsvp2       = vsvp*vsvp;
  std::vector<int> seismicType = model_settings->getSeismicType(thisTimeLapse);
  std::vector<float>        angles = model_settings->getAngle(thisTimeLapse);


  for(i = 0; i < numberOfAngles; i++)
  {
    double angle = static_cast<double>(angles[i]);
    A[i] = new float[3];
    double sint  = sin(angle);
    double sint2 = sint*sint;
    if(seismicType[i] == ModelSettings::STANDARDSEIS) {  //PP
      double tan2t=tan(angle)*tan(angle);

      A[i][0] = float( (1.0 +tan2t )/2.0 ) ;
      A[i][1] = float( -4*vsvp2 * sint2 );
      A[i][2] = float( (1.0-4.0*vsvp2*sint2)/2.0 );
    }
    else if(seismicType[i] == ModelSettings::PSSEIS) {
      double cost  = cos(angle);
      double cosp  = sqrt(1-vsvp2*sint2);
      double fac   = 0.5*sint/cosp;

      A[i][0] = 0;
      A[i][1] = float(4.0*fac*(vsvp2*sint2-vsvp*cost*cosp));
      A[i][2] = float(fac*(-1.0+2*vsvp2*sint2+2*vsvp*cost*cosp));
    }
  }
  reflectionMatrix = A;
  double vpvs = 1.0f/vsvp;
  LogKit::LogFormatted(LogKit::Low,"\nMaking reflection parameters using a Vp/Vs ratio of %4.2f\n",vpvs);
  std::string text;
  if (vpvs < model_settings->getVpVsRatioMin()) {
    LogKit::LogFormatted(LogKit::Warning,"\nA very small Vp/Vs-ratio has been detected. Values below %.2f are regarded unlikely. \n",model_settings->getVpVsRatioMin());
    text  = "Check the Vp/Vs-ratio. A small value has been found. If the value is acceptable,\n";
    text += "   you can remove this task using the <minimim-vp-vs-ratio> keyword.\n";
    TaskList::addTask(text);
  }
  else if (vpvs > model_settings->getVpVsRatioMax()) {
    LogKit::LogFormatted(LogKit::Warning,"\nA very large Vp/Vs-ratio has been detected. Values above %.2f are regarded unlikely. \n",model_settings->getVpVsRatioMax());
    text  = "Check the Vp/Vs-ratio. A large value has been found. If the value is acceptable,\n";
    text += "   you can remove this task using the <maximum-vp-vs-ratio> keyword.\n";
    TaskList::addTask(text);
  }
}

bool CommonData::WaveletHandling(ModelSettings * model_settings,
                                 InputFiles * input_files)
{
  (void) input_files;

  int nTimeLapses = model_settings->getNumberOfTimeLapses();
  //int error = 0;

  for(int i = 0; i < nTimeLapses; i++) {

    Wavelet ** wavelet;               ///< Wavelet for angle MOVE TO CLASS VARIABLE

    int nAngles = model_settings->getNumberOfAngles(i);

    //Fra ModelAvoDynamic::processSeismic:
    std::vector<bool> estimateWavelets = model_settings->getEstimateWavelet(i);
    std::vector<bool> useRickerWavelet = model_settings->getUseRickerWavelet(i);

    //double wall=0.0, cpu=0.0;
    //TimeKit::getTime(wall,cpu);

    wavelet = new Wavelet * [nAngles];
    std::vector<Grid2D *> localNoiseScale;       ///< Scale factors for local noise
    localNoiseScale.resize(nAngles);

    bool has3Dwavelet = false;
    for(int j=0; j <nAngles; j++) {
      localNoiseScale[i] = NULL;
      if (model_settings->getWaveletDim(i) == Wavelet::THREE_D)
        has3Dwavelet = true;
      if(estimateWavelets[i] == true)
        model_settings->setWaveletScale(i,j,1.0);
    }

    unsigned int                      nWells = model_settings->getNumberOfWells();

    std::vector<std::vector<double> > tGradX(nWells);
    std::vector<std::vector<double> > tGradY(nWells);

    NRLib::Grid2D<float>              refTimeGradX;         ///< Time gradient in x-direction for reference time surface (t0)
    NRLib::Grid2D<float>              refTimeGradY;         ///< Time gradient in x-direction for reference time surface (t0)
    NRLib::Grid2D<float>              structureDepthGradX;  ///< Depth gradient in x-direction for structure ( correlationDirection-t0)*v0/2
    NRLib::Grid2D<float>              structureDepthGradY;  ///< Depth gradient in y-direction for structure ( correlationDirection-t0)*v0/2

    std::string errText("");
    //bool failed = false;

    //if (has3Dwavelet) {
    //  if (input_files->getRefSurfaceFile() != "") {
    //    Surface  t0Surf;
    //    try {
    //      t0Surf =Surface(input_files->getRefSurfaceFile());
    //    }
    //    catch (NRLib::Exception & e) {
    //      errText += e.what();
    //      failed = true;
    //    }
    //     if(!failed)
    //    {
    //      double v0=model_settings->getAverageVelocity();
    //      computeStructureDepthGradient(v0,
    //                                    model_settings->getGradientSmoothingRange(),
    //                                    timeSimbox,
    //                                    &t0Surf,
    //                                    correlationDirection,
    //                                    structureDepthGradX,
    //                                    structureDepthGradY);
    //      Wavelet3D::setGradientMaps(structureDepthGradX,
    //                                 structureDepthGradY);
    //      computeReferenceTimeGradient(timeSimbox,
    //                                   &t0Surf,
    //                                   refTimeGradX,
    //                                   refTimeGradY);
    //    }
    //    else{
    //      errText += "Problems reading reference time surface in (x,y).\n";
    //      error = 1;
    //    }
    //  }
    //  bool estimateWellGradient = model_settings->getEstimateWellGradientFromSeismic();
    //  float distance, sigma_m;
    //  model_settings->getTimeGradientSettings(distance, sigma_m, i);
    //  std::vector<std::vector<double> > SigmaXY;
    //  for (unsigned int w=0; w<nWells; w++) {
    //    BlockedLogs *bl    = wells[w]->getBlockedLogsOrigThick();
    //    if(!estimateWellGradient & ((structureDepthGradX.GetN()> 0) & (structureDepthGradY.GetN()>0))){ // then we allready have the
    //      double v0=model_settings->getAverageVelocity();
    //      bl->setSeismicGradient(v0,structureDepthGradX,structureDepthGradY,refTimeGradX, refTimeGradY,tGradX[w], tGradY[w]);
    //    }else{
    //      bl->setTimeGradientSettings(distance, sigma_m);
    //      bl->findSeismicGradient(seisCube, timeSimbox, nAngles,tGradX[w], tGradY[w],SigmaXY);
    //    }
    //  }
    //}

    if (estimation_simbox_.getdz() > 4.01f && model_settings->getEstimateNumberOfWavelets(i) > 0)
    { // Require this density for wavelet estimation
      LogKit::LogFormatted(LogKit::Low,"\n\nWARNING: The minimum sampling density is lower than 4.0. The WAVELETS generated by \n");
      LogKit::LogFormatted(LogKit::Low,"         CRAVA are not reliable and the output results should be treated accordingly.\n");
      LogKit::LogFormatted(LogKit::Low,"         The number of layers must be increased.                                  \n");
      std::string text("");
      text += "Increase the number of layers to improve the quality of the wavelet estimation.\n";
      text += "   The minimum sampling density is "+NRLib::ToString(estimation_simbox_.getdz())+", and it should be ";
      text += "lower than 4.0.\n   To obtain the desired density, the number of layers should be at least ";
      text += NRLib::ToString(static_cast<int>(estimation_simbox_.GetLZ()/4.0))+"\n";
      TaskList::addTask(text);
    }

    // check if local noise is set for some angles.
    //bool localNoiseSet = false;
    std::vector<float> angles = model_settings->getAngle(i);

    for (int j = 0; j < nAngles; j++) {

      float angle = float(angles[i]*180.0/NRLib::Pi);
    LogKit::LogFormatted(LogKit::Low,"\nAngle stack : %.1f deg",angle);
    //if(model_settings->getForwardModeling()==false)
    //  seisCube[i]->setAccessMode(FFTGrid::RANDOMACCESS);
    //if (model_settings->getWaveletDim(j) == Wavelet::ONE_D)
    //  error += process1DWavelet(model_settings,
    //                            input_files,
    //                            timeSimbox,
    //                            seisCube,
    //                            wells,
    //                            waveletEstimInterval,
    //                            reflectionMatrix_[j],
    //                            errText,
    //                            wavelet[j],
    //                            i,
    //                            j,
    //                            useRickerWavelet[i]);
    //else
    //  error += process3DWavelet(model_settings,
    //                            input_files,
    //                            timeSimbox,
    //                            seisCube,
    //                            wells,
    //                            waveletEstimInterval,
    //                            reflectionMatrix_[j],
    //                            errText,
    //                            wavelet[j],
    //                            i,
    //                            j,
    //                            refTimeGradX,
    //                            refTimeGradY,
    //                            tGradX,
    //                            tGradY);



      //if(estimateWavelets[i]) {
      //  if(read_seismic_ == true && read_wells_ == true && setup_reflection_matrix_ == true) {
      //    //Block wells as in optimize well-location

      //    //use ProcessWavelets (modelavodynamic.cpp), but with blocked logs instead of well-object.

      //    //Store Noise-estimate (globally and localy), local shift and scale for each wavelet.

      //  }

      //}
      //else { //From file or ricker

      //  if(useRickerWavelet[i]) {

      //  }
      //  else { //From file

      //  }

      //  if(input_files->getWaveletFile(i, j) != "")
      //    ;

      //}




    }
  }
  return true;

}

//int
//CommonData::process1DWavelet(const ModelSettings          * model_settings,
//                                  const InputFiles             * input_files,
//                                  const Simbox                 * timeSimbox,
//                                  const FFTGrid        * const * seisCube,
//                                  std::vector<WellData *>        wells,
//                                  const std::vector<Surface *> & waveletEstimInterval,
//                                  const float                  * reflectionMatrix,
//                                  std::string                  & errText,
//                                  Wavelet                     *& wavelet,
//                                  //unsigned int                   i, //Timelapse
//                                  unsigned int                   j, //Angle
//                                  bool                           useRickerWavelet)
//{
//  int error = 0;
//  Grid2D * shiftGrid(NULL);
//  Grid2D * gainGrid(NULL);
//  if(model_settings->getUseLocalWavelet() && input_files->getScaleFile(timelapse,j) != "") {
//      Surface help(input_files->getScaleFile(thisTimeLapse_,j));
//      gainGrid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
//      resampleSurfaceToGrid2D(timeSimbox, &help, gainGrid);
//  }
//  if (model_settings->getUseLocalWavelet() && input_files->getShiftFile(timelapse,j) != ""){
//    Surface helpShift(input_files->getShiftFile(thisTimeLapse_,j));
//    shiftGrid = new Grid2D(timeSimbox->getnx(),timeSimbox->getny(), 0.0);
//    resampleSurfaceToGrid2D(timeSimbox, &helpShift, shiftGrid);
//  }
//  if (useLocalNoise_ && input_files->getLocalNoiseFile(thisTimeLapse_,j) != ""){
//    Surface helpNoise(input_files->getLocalNoiseFile(thisTimeLapse_,j));
//    localNoiseScale_[i] = new Grid2D(timeSimbox->getnx(), timeSimbox->getny(), 0.0);
//    resampleSurfaceToGrid2D(timeSimbox, &helpNoise, localNoiseScale_[i]);
//  }
//
//  if (estimateWavelet_[j])
//    wavelet = new Wavelet1D(timeSimbox,
//                            seisCube[j],
//                            wells,
//                            waveletEstimInterval,
//                            model_settings,
//                            reflectionMatrix,
//                            i,
//                            error,
//                            errText);
//
//  else { //Not estimation modus
//    if(useRickerWavelet)
//        wavelet = new Wavelet1D(model_settings,
//                                reflectionMatrix,
//                                angle_[j],
//                                model_settings->getRickerPeakFrequency(thisTimeLapse_,j),
//                                error);
//    else {
//      const std::string & waveletFile = input_files->getWaveletFile(thisTimeLapse_,j);
//      int fileFormat = getWaveletFileFormat(waveletFile,errText);
//      if(fileFormat < 0) {
//        errText += "Unknown file format of file '"+waveletFile+"'.\n";
//        error++;
//      }
//      else
//        wavelet = new Wavelet1D(waveletFile,
//                                fileFormat,
//                                model_settings,
//                                reflectionMatrix,
//                                angle_[j],
//                                error,
//                                errText);
//    }
//      // Calculate a preliminary scale factor to see if wavelet is in the same size order as the data. A large or small value might cause problems.
//      if(seisCube!=NULL) {// If forward modeling, we have no seismic, can not prescale wavelet.
//        float       prescale = wavelet->findGlobalScaleForGivenWavelet(model_settings, timeSimbox, seisCube[i], wells);
//        const float limHigh  = 3.0f;
//        const float limLow   = 0.33f;
//
//        if(model_settings->getEstimateGlobalWaveletScale(thisTimeLapse_,i)) // prescale, then we have correct size order, and later scale estimation will be ok.
//           wavelet->multiplyRAmpByConstant(prescale);
//        else {
//          if(model_settings->getWaveletScale(thisTimeLapse_,i)!= 1.0f && (prescale>limHigh || prescale<limLow)) {
//             std::string text = "The wavelet given for angle no "+NRLib::ToString(i)+" is badly scaled. Ask Crava to estimate global wavelet scale.\n";
//            if(model_settings->getEstimateLocalScale(thisTimeLapse_,i)) {
//              errText += text;
//              error++;
//            }
//            else {
//              LogKit::LogFormatted(LogKit::Warning,"\nWARNING: "+text);
//              TaskList::addTask("The wavelet is badly scaled. Consider having CRAVA estimate global wavelet scale");
//            }
//          }
//        }
//      }
//      if (error == 0)
//        wavelet->resample(static_cast<float>(timeSimbox->getdz()),
//                          timeSimbox->getnz(),
//                          model_settings->getNZpad());
//  }
//
//  if (error == 0) {
//    wavelet->scale(model_settings->getWaveletScale(thisTimeLapse_,j));
//
//    if (model_settings->getForwardModeling() == false && model_settings->getNumberOfWells() > 0) {
//      float SNRatio = wavelet->calculateSNRatioAndLocalWavelet(timeSimbox,
//                                                               seisCube[j],
//                                                               wells,
//                                                               model_settings,
//                                                               errText,
//                                                               error,
//                                                               j,
//                                                               localNoiseScale_[j],
//                                                               shiftGrid,
//                                                               gainGrid,
//                                                               SNRatio_[j],
//                                                               model_settings->getWaveletScale(thisTimeLapse_,j),
//                                                               model_settings->getEstimateSNRatio(thisTimeLapse_,j),
//                                                               model_settings->getEstimateGlobalWaveletScale(thisTimeLapse_,j),
//                                                               model_settings->getEstimateLocalNoise(thisTimeLapse_,j),
//                                                               model_settings->getEstimateLocalShift(thisTimeLapse_,j),
//                                                               model_settings->getEstimateLocalScale(thisTimeLapse_,j),
//                                                               estimateWavelet_[j]);
//
//      if(model_settings->getEstimateSNRatio(thisTimeLapse_,j))
//        SNRatio_[i] = SNRatio;
//    }
//
//    if (error == 0) {
//      if((model_settings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
//         (model_settings->getEstimationMode() && estimateWavelet_[j])) {
//        std::string type;
//        if (estimateWavelet_[j]) {
//          type = "Estimated_";
//          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,true); // dt_max = 1.0;
//        }
//        else if (model_settings->getWaveletScale(thisTimeLapse_,j) == 1.00) {
//          type = "";
//          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
//        }
//        else {
//          type = "Scaled_";
//          wavelet->writeWaveletToFile(IO::PrefixWavelet()+type, 1.0,false); // dt_max = 1.0;
//        }
//      }
//      const float SNLow  = 1.0;
//      const float SNHigh = 10.0;
//      if ((SNRatio_[j] <=SNLow  || SNRatio_[j] > SNHigh) && model_settings->getForwardModeling()==false) {
//        errText += "Illegal signal-to-noise ratio of "+NRLib::ToString(SNRatio_[j])+" for cube "+NRLib::ToString(j+1)+".\n";
//        errText += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
//        error++;
//      }
//
//      bool useLocalNoise = model_settings->getEstimateLocalNoise(thisTimeLapse_,j) || input_files->getLocalNoiseFile(thisTimeLapse_,j) != "";
//      bool useLocalShift = model_settings->getEstimateLocalShift(thisTimeLapse_,j) || input_files->getShiftFile(thisTimeLapse_,j)      != "";
//      bool useLocalGain  = model_settings->getEstimateLocalScale(thisTimeLapse_,j) || input_files->getScaleFile(thisTimeLapse_,j)      != "";
//
//      if (useLocalNoise)
//        readAndWriteLocalGridsToFile(input_files->getLocalNoiseFile(thisTimeLapse_,j),
//                                     IO::PrefixLocalNoise(),
//                                     1.0,  // Scale map with this factor before writing to disk
//                                     model_settings,
//                                     j,
//                                     timeSimbox,
//                                     localNoiseScale_[j]);
//
//      if (useLocalShift) {
//        readAndWriteLocalGridsToFile(input_files->getShiftFile(thisTimeLapse_,j),
//                                     IO::PrefixLocalWaveletShift(),
//                                     1.0,
//                                     model_settings,
//                                     j,
//                                     timeSimbox,
//                                     shiftGrid);
//        wavelet->setShiftGrid(shiftGrid);
//      }
//
//      if (useLocalGain) {
//        readAndWriteLocalGridsToFile(input_files->getScaleFile(thisTimeLapse_,j),
//                                     IO::PrefixLocalWaveletGain(),
//                                     1.0,
//                                     model_settings,
//                                     j,
//                                     timeSimbox,
//                                     gainGrid);
//        wavelet->setGainGrid(gainGrid);
//      }
//    }
//  }
//  return error;
//}


//bool CommonData::EstimateWaveletShape() {
//  return true;
//}
//
//bool CommonData::SetupEstimationRockPhysics(){
//  return true;
//}

int
CommonData::ComputeTime(int year, int month, int day) const
{
  if(year == IMISSING)
    return(0);

  int delta_year = year-1900; //Ok baseyear.
  int time = 365*delta_year+delta_year/4; //Leap years.
  if(month == IMISSING)
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

    if(day == IMISSING)
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

  if(grid_file != "") { //May change the condition here, but need geometry if we want to set XL/IL
    int file_type = IO::findGridType(grid_file);
    if(file_type == IO::CRAVA) {
      geometry = GetGeometryFromCravaFile(grid_file);
    }
    else if(file_type == IO::SEGY) {
      try
      {
        geometry = SegY::FindGridGeometry(grid_file, thf);
      }
      catch (NRLib::Exception & e)
      {
        err_text = e.what();
      }
    }
    else if(file_type == IO::STORM)
      geometry = GetGeometryFromStormFile(grid_file, err_text);
    else if(file_type==IO::SGRI) {
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

  double x0      = NRLib::ReadBinaryDouble(bin_file);
  double y0      = NRLib::ReadBinaryDouble(bin_file);
  double dx      = NRLib::ReadBinaryDouble(bin_file);
  double dy      = NRLib::ReadBinaryDouble(bin_file);
  int    nx      = NRLib::ReadBinaryInt(bin_file);
  int    ny      = NRLib::ReadBinaryInt(bin_file);
  double IL0     = NRLib::ReadBinaryDouble(bin_file);
  double XL0     = NRLib::ReadBinaryDouble(bin_file);
  double il_step_X = NRLib::ReadBinaryDouble(bin_file);
  double il_step_Y = NRLib::ReadBinaryDouble(bin_file);
  double xl_step_X = NRLib::ReadBinaryDouble(bin_file);
  double xl_step_Y = NRLib::ReadBinaryDouble(bin_file);
  double rot     = NRLib::ReadBinaryDouble(bin_file);

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
  if(scale==false)
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
    double x0      = storm_grid->GetXMin()*scale_hor;
    double y0      = storm_grid->GetYMin()*scale_hor;
    double dx      = storm_grid->GetDX()*scale_hor;
    double dy      = storm_grid->GetDY()*scale_hor;
    int    nx      = static_cast<int>(storm_grid->GetNI());
    int    ny      = static_cast<int>(storm_grid->GetNJ());
    double rot     = storm_grid->GetAngle();
    double IL0     = 0.0;  ///< Dummy value since information is not contained in format
    double XL0     = 0.0;  ///< Dummy value since information is not contained in format
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

void CommonData::SetSurfacesSingleInterval(const ModelSettings              * const model_settings,
                                           NRLib::Volume                    & full_inversion_volume,
                                           Simbox                           & estimation_simbox,
                                           const std::vector<std::string>   & surf_file,
                                           std::string                      & err_text,
                                           bool                             & failed)
{
  const std::string & top_surface_file_name = surf_file[0];

  //bool   generate_seismic    = model_settings->getForwardModeling();
  //bool   estimation_mode     = model_settings->getEstimationMode();
  //bool   generate_background = model_settings->getGenerateBackground();
  bool   parallel_surfaces   = model_settings->getParallelTimeSurfaces();
  //int    nz                 = model_settings->getTimeNz();
  //int    output_format       = model_settings->getOutputGridFormat();
  //int    output_domain       = model_settings->getOutputGridDomain();
  //int    output_grids_elastic = model_settings->getOutputGridsElastic();
  //int    output_grids_other   = model_settings->getOutputGridsOther();
  //int    output_grids_seismic = model_settings->getOutputGridsSeismic();
  double d_top               = model_settings->getTimeDTop();
  double lz                 = model_settings->getTimeLz();
  //double dz                 = model_settings->getTimeDz();

  Surface * top_surface = NULL;
  Surface * base_surface = NULL;

  Surface * top_surface_flat = NULL;
  Surface * base_surface_flat = NULL;

  try {
    double x_min, x_max;
    double y_min, y_max;
    FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                  estimation_simbox.getlx(), estimation_simbox.getly(),
                                  estimation_simbox.getAngle(),
                                  x_min,y_min,x_max,y_max);
    if (NRLib::IsNumber(top_surface_file_name)) {
      // Find the smallest surface that covers the simbox. For simplicity
      // we use only four nodes (nx=ny=2).
      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(top_surface_file_name.c_str()));
    }
    else {
      top_surface = new Surface(top_surface_file_name);
    }
  }
  catch (NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  if(!failed) {
    if(parallel_surfaces) { //Only one reference surface
      //simbox->setDepth(*top_surface_flat, d_top, lz, dz, model_settings->getRunFromPanel());
      top_surface->Add(d_top);
      base_surface = new Surface(*top_surface);
      base_surface->Add(lz);
      //full_inversion_volume->SetSurfaces(*top_surface, *base_surface, model_settings->getRunFromPanel());
    }
    else { //Two reference surfaces
      const std::string & base_surface_file_name = surf_file[1];
      try {
        if (NRLib::IsNumber(base_surface_file_name)) {
          // Find the smallest surface that covers the simbox. For simplicity
          // we use only four nodes (nx=ny=2).
          double x_min, x_max;
          double y_min, y_max;
          FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                      estimation_simbox.getlx(), estimation_simbox.getly(),
                                      estimation_simbox.getAngle(),
                                      x_min, y_min, x_max, y_max);
          base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(base_surface_file_name.c_str()));
        }
        else {
          base_surface = new Surface(base_surface_file_name);
        }
      }
      catch (NRLib::Exception & e) {
        err_text += e.what();
        failed = true;
      }

    }
    /*
    if (!failed) {
      if((output_domain & IO::TIMEDOMAIN) > 0) {
        std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime();
        std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime();
        simbox->setTopBotName(topSurf,baseSurf,output_format);
        if (generate_seismic) {
          simbox->writeTopBotGrids(topSurf,
                                   baseSurf,
                                   IO::PathToSeismicData(),
                                   output_format);
        }
        else if (!estimation_mode){
          if (output_grids_elastic > 0 || output_grids_other > 0 || output_grids_seismic > 0)
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToInversionResults(),
                                     output_format);
        }
        if((output_format & IO::STORM) > 0) { // These copies are only needed with the STORM format
          if ((output_grids_elastic & IO::BACKGROUND) > 0 ||
              (output_grids_elastic & IO::BACKGROUND_TREND) > 0 ||
              (estimation_mode && generate_background)) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToBackground(),
                                     output_format);
          }
          if ((output_grids_other & IO::CORRELATION) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToCorrelations(),
                                     output_format);
          }
          if ((output_grids_seismic & (IO::ORIGINAL_SEISMIC_DATA | IO::SYNTHETIC_SEISMIC_DATA)) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToSeismicData(),
                                     output_format);
          }
          if ((output_grids_other & IO::TIME_TO_DEPTH_VELOCITY) > 0) {
            simbox->writeTopBotGrids(topSurf,
                                     baseSurf,
                                     IO::PathToVelocity(),
                                     output_format);
          }
        }
      }
    } */
  }
  if(!failed){
    try{ //initialize full_inversion_volume and set the flat top and base surfaces of the simbox

      full_inversion_volume = NRLib::Volume(estimation_simbox.getx0(), estimation_simbox.gety0(),
        estimation_simbox.getlx(), estimation_simbox.getly(),
        *top_surface, *base_surface, estimation_simbox.getAngle());

      // find the top and bottom points on the full inversion volume
      double z_min_top_surface = full_inversion_volume.GetTopZMin(2,2);
      double z_max_base_surface = full_inversion_volume.GetBotZMax(2,2);
      // round outwards from the center of the volume to nearest multiple 4
      double z_min = floor(z_min_top_surface/4)*4;
      double z_max = ceil(z_max_base_surface/4)*4;
      top_surface_flat = new Surface(top_surface->GetXMin(), top_surface->GetYMin(),
        top_surface->GetXMax()-top_surface->GetXMin(), top_surface->GetYMax()-top_surface->GetYMin(),
        2, 2, z_min);
      base_surface_flat = new Surface(base_surface->GetXMin(), base_surface->GetYMin(),
        base_surface->GetXMax()-base_surface->GetXMin(), base_surface->GetYMax()-base_surface->GetYMin(),
        2, 2, z_max);
      estimation_simbox.setDepth(*top_surface_flat, *base_surface_flat, 2, model_settings->getRunFromPanel());
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
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
                                              const InputFiles                * input_files,
                                              std::string                     & err_text,
                                              bool                            & failed){

  // Get interval surface data ------------------------------------------------------------------------------

  const std::vector<std::string>            interval_names                         =  model_settings->getIntervalNames();
  unsigned int                              n_intervals                            =  interval_names.size();
  const std::map<std::string, std::string>  interval_base_time_surfaces            =  input_files->getIntervalBaseTimeSurfaces();
  // implicit assumption that the top surface is the first one given and the base surface is the last one
  const std::string                         top_surface_file_name                  =  input_files->getTimeSurfFile(0);
  const std::string                         base_surface_file_name                 =  interval_base_time_surfaces.find( interval_names[n_intervals-1])->second;

  Surface * top_surface = NULL;
  Surface * base_surface = NULL;

  Surface * top_surface_flat = NULL;
  Surface * base_surface_flat = NULL;

  try{
    if(NRLib::IsNumber(top_surface_file_name)){
      double x_min, x_max;
      double y_min, y_max;
      FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                  estimation_simbox.getlx(), estimation_simbox.getly(),
                                  estimation_simbox.getAngle(), x_min,y_min,x_max,y_max);
      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(top_surface_file_name.c_str()));
    }
    else{
      top_surface = new Surface(top_surface_file_name);
    }

  }
  catch(NRLib::Exception & e) {
    err_text += e.what();
    failed = true;
  }

  if(!failed){
    try{
      if(NRLib::IsNumber(base_surface_file_name)){
        double x_min, x_max;
        double y_min, y_max;
        FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                    estimation_simbox.getlx(), estimation_simbox.getly(),
                                    estimation_simbox.getAngle(), x_min,y_min,x_max,y_max);
        base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(base_surface_file_name.c_str()));
      }
      else{
        base_surface = new Surface(base_surface_file_name);
      }
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }
  if(!failed){
    try{ // initialize full_inversion_volume and set the flat top and base surfaces of the simbox
      full_inversion_volume =  NRLib::Volume(estimation_simbox.getx0(), estimation_simbox.gety0(),
        estimation_simbox.getlx(), estimation_simbox.getly(),
        *top_surface, *base_surface, estimation_simbox.getAngle());

      // find the top and bottom points on the full inversion volume
      double z_min_top_surface = full_inversion_volume.GetTopZMin(2,2);
      double z_max_base_surface = full_inversion_volume.GetBotZMax(2,2);
      // round outwards from the center of the volume to nearest multiple 4
      double z_min = floor(z_min_top_surface/4)*4;
      double z_max = ceil(z_max_base_surface/4)*4;
      top_surface_flat = new Surface(top_surface->GetXMin(), top_surface->GetYMin(),
        top_surface->GetXMax()-top_surface->GetXMin(), top_surface->GetYMax()-top_surface->GetYMin(),
        2, 2, z_min);
      base_surface_flat = new Surface(base_surface->GetXMin(), base_surface->GetYMin(),
        base_surface->GetXMax()-base_surface->GetXMin(), base_surface->GetYMax()-base_surface->GetYMin(),
        2, 2, z_max);
      estimation_simbox.setDepth(*top_surface_flat, *base_surface_flat, 2, model_settings->getRunFromPanel());
    }
    catch(NRLib::Exception & e){
      err_text += e.what();
      failed = true;
    }
  }

  delete top_surface;
  delete base_surface;
  delete top_surface_flat;
  delete base_surface_flat;

}


bool CommonData::BlockWellsForEstimation(const ModelSettings                            * const model_settings,
                                         const Simbox                                   & estimation_simbox,
                                         std::vector<NRLib::Well>                       & wells,
                                         std::map<std::string, BlockedLogsCommon *>     & mapped_blocked_logs_common,
                                         std::string                                    & err_text){
  bool failed = false;

  LogKit::WriteHeader("Blocking wells for estimation");

  // These are the continuous parameters that are to be used in BlockedLogs
  continuous_logs_to_be_blocked_.push_back("Vp");
  continuous_logs_to_be_blocked_.push_back("Vs");
  continuous_logs_to_be_blocked_.push_back("Rho");
  continuous_logs_to_be_blocked_.push_back("MD");

  // These are the discrete parameters that are to be used in BlockedLogs
  

  try{
    for (unsigned int i=0; i<wells.size(); i++){
      BlockedLogsCommon * blocked_log = new BlockedLogsCommon(&wells[i], continuous_logs_to_be_blocked_, discrete_logs_to_be_blocked_, 
        &estimation_simbox, model_settings->getRunFromPanel(), failed, err_text);
      mapped_blocked_logs_common.insert(std::pair<std::string, BlockedLogsCommon *>(wells[i].GetWellName(), blocked_log));
    }
  }catch(NRLib::Exception & e){
    err_text += e.what();
    failed = true;
  }

  return !failed;
}


void  CommonData::OptimizeWellLocations(ModelSettings                                 * model_settings,
                                        InputFiles                                    * input_files,
                                        const Simbox                                  * estimation_simbox,
                                        //const NRLib::Volume                           & volume,
                                        std::vector<NRLib::Well>                      & wells,
                                        std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs,
                                        std::map<int, std::vector<SeismicStorage> >   & seismic_data,
                                        std::map<int, float **>                       & reflection_matrix,
                                        std::string                                   & err_text,
                                        bool                                          & failed){

  LogKit::WriteHeader("Estimating optimized well location");

  std::vector<Surface *>      well_move_interval;
  LoadWellMoveInterval(input_files, estimation_simbox, well_move_interval, err_text, failed);

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
    if( wells[w].IsDeviated())
      continue;
    std::string well_name = wells[w].GetWellName();
    std::map<std::string, BlockedLogsCommon * >::iterator it = mapped_blocked_logs.find(well_name);
    if(it == mapped_blocked_logs.end()){
      failed = true;
      err_text += "Blocked log not found for well  " + wells[w].GetWellName() + "\n";
      break;
    }
    BlockedLogsCommon * bl = it->second;
    n_move_angles = model_settings->getNumberOfWellAngles(w);

    if( n_move_angles==0 )
      continue;

    for(int i=0; i<n_angles; i++ )
      angle_weight[i] = 0;

    for(int i=0; i<n_move_angles; i++ ){
      move_angle   = model_settings->getWellMoveAngle(w,i);

      for(int j=0; j<n_angles; j++ ){
        if( move_angle == seismicAngle[j]){
          angle_weight[j] = model_settings->getWellMoveWeight(w,i);
          break;
        }
      }
    }

    sum = 0;
    for(int i=0; i<n_angles; i++ )
      sum += angle_weight[i];
    if( sum == 0 )
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
                                                                              estimation_simbox, model_settings->getRunFromPanel(), failed, err_text) ) );
    LogKit::LogFormatted(LogKit::Low,"  %-13s %11.2f %12d %11.2f %8d %11.2f \n",
    wells[w].GetWellName().c_str(), k_move, i_move, delta_X, j_move, delta_Y);
  }

   for (int w = 0 ; w < static_cast<int>(n_wells) ; w++){
     n_move_angles = model_settings->getNumberOfWellAngles(w);

    if( wells[w].IsDeviated()==true && n_move_angles > 0 )
    {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: Well %7s is treated as deviated and can not be moved.\n",
          wells[w].GetWellName().c_str());
      TaskList::addTask("Well "+NRLib::ToString(wells[w].GetWellName())+" can not be moved. Remove <optimize-location-to> for this well");
    }
   }

}

void CommonData::MoveWell(const NRLib::Well & well,
                          const Simbox      * simbox, 
                          double              delta_X, 
                          double              delta_Y, 
                          double              k_move){

  double delta_Z;
  double top_old, top_new;

  if(well.HasContLog("X_pos") && well.HasContLog("Y_pos")){
    std::vector<double> x_pos = well.GetContLog("X_pos");
    std::vector<double> y_pos = well.GetContLog("Y_pos");
    std::vector<double> z_pos = well.GetContLog("TVD");
    top_old = simbox->getTop(x_pos[0], y_pos[0]);

    for(unsigned int i=0; i<x_pos.size(); i++)
    {
      if(x_pos[i] != RMISSING)
        x_pos[i] = x_pos[i]+delta_X;
      if(y_pos[i] != RMISSING)
        y_pos[i] = y_pos[i]+delta_Y;
    }

    top_new = simbox->getTop(x_pos[0], y_pos[0]);

    delta_Z = top_new - top_old + k_move;

    for(unsigned int i=0; i<z_pos.size(); i++)
      z_pos[i] = z_pos[i]+delta_Z;
  }
}

void  CommonData::CalculateDeviation(NRLib::Well            & new_well,
                                     const ModelSettings    * const model_settings,
                                     float                  & dev_angle,
                                     Simbox                 * simbox,
                                     int                      use_for_wavelet_estimation)
{
  float maxDevAngle   = model_settings->getMaxDevAngle();
  float thr_deviation = float(tan(maxDevAngle*NRLib::Pi/180.0));  // Largest allowed deviation
  float max_deviation =  0.0f;
  float max_dz        = 10.0f;                      // Calculate slope each 10'th millisecond.

  std::vector<double> x_pos = new_well.GetContLog("X_pos");
  std::vector<double> y_pos = new_well.GetContLog("Y_pos");
  std::vector<double> z_pos = new_well.GetContLog("TVD");

  //
  // Find first log entry in simbox
  //
  int iFirst = IMISSING;
  for(unsigned int i=0 ; i < x_pos.size() ; i++)
  {
    if(simbox->isInside(x_pos[i], y_pos[i]))
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
    for(unsigned int i = iFirst + 1 ; i < x_pos.size() ; i++)
    {
      if(simbox->isInside(x_pos[i], y_pos[i]))
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
      if(use_for_wavelet_estimation == ModelSettings::NOTSET) {
        use_for_wavelet_estimation = ModelSettings::NO;
      }
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
                                      std::string                  & err_text,
                                      bool                         & failed){

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
      failed = true;
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
      failed = true;
    }
  }
}

void   CommonData::SetupTrendCubes(ModelSettings                  * model_settings, 
                                   InputFiles                     * input_files, 
                                   MultiIntervalGrid              * multiple_interval_grid,
                                   std::string                    & error_text,
                                   bool                           & failed){

  const std::vector<std::string> trend_cube_parameters = model_settings->getTrendCubeParameters();
  const std::vector<int>         trend_cube_type       = model_settings->getTrendCubeType();
  trend_cubes_.resize(multiple_interval_grid->GetNIntervals());
  const std::vector<std::string>            interval_names                         =  model_settings->getIntervalNames();
  try{
    for (int i = 0; i<multiple_interval_grid->GetNIntervals(); i++){

      trend_cubes_[i] = CravaTrend(multiple_interval_grid->GetIntervalSimbox(i),
                                   multiple_interval_grid->GetSimbox(i),
                                   model_settings,
                                   input_files,
                                   interval_names[i],
                                   trend_cube_type,
                                   trend_cube_parameters,
                                   failed,
                                   error_text);

    }
  }catch(NRLib::Exception & e){
    error_text+= e.what();
    failed = true;
  }
}

void CommonData::SetupRockPhysics(const ModelSettings                               * model_settings,
                                  const InputFiles                                  * input_files,
                                  const MultiIntervalGrid                           * multiple_interval_grid,
                                  const std::vector<CravaTrend>                     & trend_cubes,
                                  const std::map<std::string, BlockedLogsCommon *>  & mapped_blocked_logs,
                                  int                                                 n_trend_cubes,
                                  std::string                                       & error_text,
                                  bool                                              & failed)
{
  LogKit::WriteHeader("Processing Rock Physics");

  // rock physics data
  const std::vector<std::string>                    interval_names          = model_settings->getIntervalNames();
  int                                               n_vintages              = model_settings->getNumberOfVintages();
  const std::string                                 path                    = input_files->getInputDirectory();
  const std::vector<std::string>                    trend_cube_parameters   = model_settings->getTrendCubeParameters();
  std::vector<std::vector<std::vector<double> > >   trend_cube_sampling(n_trend_cubes_);
  const std::vector<std::vector<float> >            dummy_blocked_logs;
  const std::map<std::string, 
    std::vector<DistributionWithTrendStorage *> >   reservoir_variable      = model_settings->getReservoirVariable();
  
  //
  reservoir_variables_.resize(n_trend_cubes);
  rock_distributions_.resize(n_trend_cubes);

  // generate distribution for each reservoir variable
  for(int i=0; i<n_trend_cubes; i++){ // the number of trend cubes is the same as the number of intervals
    trend_cube_sampling[i]                                                                     = trend_cubes[i].GetTrendCubeSampling();
    for(std::map<std::string, std::vector<DistributionWithTrendStorage *> >::const_iterator it = reservoir_variable.begin(); it != reservoir_variable.end(); it++) {

      std::vector<DistributionWithTrendStorage *>   storage = it->second;
      std::vector<DistributionWithTrend *>          dist_vector(storage.size());

      for(size_t j=0; j<storage.size(); j++) {
        dist_vector[j]                    = storage[j]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling[i], error_text);
      }

      reservoir_variables_[i][it->first]  = dist_vector;
    }
  }
  

  if(error_text == "") {

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
    int   n_wells       = model_settings->getNumberOfWells();

    // Block logs
    /*
    std::vector<BlockedLogsCommon *> blocked_logs_rock_physics(n_wells, NULL);
    if(n_wells > 0) {
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
    std::vector<std::vector<std::vector<DistributionWithTrend *> > > res_var_vintage(n_trend_cubes_);
    for(int i=0; i<n_trend_cubes_; i++){
      if(reservoir_variables_[i].size() > 0) {
        size_t n_vintages = reservoir_variables_[i].begin()->second.size();
        res_var_vintage[i].resize(n_vintages);
        for(std::map<std::string, std::vector<DistributionWithTrend *> >::iterator var_it = reservoir_variables_[i].begin();
          var_it != reservoir_variables_[i].end();var_it++)
        {
          for(size_t vin_index=0; vin_index < var_it->second.size();vin_index++)
            res_var_vintage[i][vin_index].push_back((var_it->second)[vin_index]);
        }
      }
    }

    // Facies names
    std::map<std::string, std::map<std::string, float> > prior_facies_prob_interval = model_settings->getPriorFaciesProbInterval();
    //std::map<std::string, float> facies_probabilities = model_settings->getPriorFaciesProb();
    std::map<std::string, std::string> facies_cubes   = input_files->getPriorFaciesProbFile();
    std::vector<std::string> all_facies_names         = facies_names_;

    // EN hva skjer her?
    /*
    for(std::map<std::string, float>::iterator it_prob = facies_probabilities.begin(); it_prob != facies_probabilities.end(); it_prob++)
      all_facies_names.push_back(it_prob->first);
    for(std::map<std::string, std::string>::iterator it_cube = facies_cubes.begin(); it_cube != facies_cubes.end(); it_cube++)
      all_facies_names.push_back(it_cube->first);
    */

    std::sort(all_facies_names.begin(), all_facies_names.end());

    std::string prev_facies = "";

    for(size_t i = 0; i<interval_names.size(); i++){
      for(size_t f=0; f<all_facies_names.size(); f++) {
        if(all_facies_names[f] != prev_facies) {
          prev_facies = all_facies_names[f];

          std::map<std::string, DistributionsRockStorage *>::const_iterator iter = rock_storage.find(all_facies_names[f]);
          if(iter != rock_storage.end()) {

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

            if(rock_err_txt == "") {

              int n_vintages =  1;//static_cast<int>(rock.size());
              if(n_vintages > 1)
                LogKit::LogFormatted(LogKit::Low, "Number of vintages: %4d\n", n_vintages);

              for(int t=0; t<n_vintages; t++) {
                if(n_vintages > 1)
                  LogKit::LogFormatted(LogKit::Low, "\nVintage number: %4d\n", t+1);

                //Completing the top level rocks, by setting access to reservoir variables and sampling distribution.
                //rock[t]->CompleteTopLevelObject(res_var_vintage[i][t]);

                std::vector<bool> has_trends (1,false);// = rock[t]->HasTrend();
                bool              has_trend = false;
                for(size_t j=0; j<has_trends.size(); j++) {
                  if(has_trends[j] == true) {
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
                  if(typeid(*(storage)) == typeid(ReussRockStorage))
                    tmp_err_txt += "Vs value of 0 detected. Note that the Reuss model gives Vs=0; hence it can not be used to model a facies\n";
                  else
                    tmp_err_txt += "Vs value of "+NRLib::ToString(std::exp(expectation[1]))+" detected: ";
                  tmp_err_txt += "Vs should be in the interval ("+NRLib::ToString(beta_min)+", "+NRLib::ToString(beta_max)+") m/s\n";
                }
                if (std::exp(expectation[2]) < rho_min  || std::exp(expectation[2]) > rho_max) {
                  tmp_err_txt += "Rho value of "+NRLib::ToString(std::exp(expectation[2]))+" detected: ";
                  tmp_err_txt += "Rho should be in the interval ("+NRLib::ToString(rho_min)+", "+NRLib::ToString(rho_max)+") g/cm^3\n";
                }

                if(tmp_err_txt != "") {
                  error_text += "\nToo high or low seismic properties calculated for rock '"+iter->first+"':\n";
                  error_text += tmp_err_txt;
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

                if(var_err_txt != "") {
                  error_text += "\nToo high or low variance of seismic properties calculated for rock '"+iter->first+"':\n";
                  error_text += var_err_txt;
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
                  error_text += e.what();
                  error_text += " for rock '"+iter->first+"':\n";
                  error_text += "  The variables in the rock model are probably linearly dependent\n";
                }

                if(var_err_txt != "" || tmp_err_txt != "")
                  break;

              }
              // set rock distribution for interval i and facies f
              //rock_distributions_[interval_names[i]][all_facies_names[f]] = rock;
            }
            else
              error_text += rock_err_txt;
          }
          else
            error_text += "The facies "+all_facies_names[f]+" is not one of the rocks in the rock physics model\n";
          //rock_distributions_[i][all_facies_names[f]] = rock;
        }
      }
    }
    /*
    for(int i=0; i<n_wells; i++)
        delete blocked_logs_rock_physics[i];
        */
  }

  if(error_text != "")
    failed = true;
}

//-----------------------------------------------------------------------------------------------------

void CommonData::PrintExpectationAndCovariance(const std::vector<double>   & expectation,
                                               const NRLib::Grid2D<double> & covariance,
                                               const bool                  & has_trend) const{
  if(has_trend == true)
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