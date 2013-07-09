/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/commondata.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/timeline.h"
#include "src/fftgrid.h"
#include "src/modelgeneral.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"
#include "src/seismicstorage.h"
#include "nrlib/well/norsarwell.hpp"

#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "nrlib/segy/segytrace.hpp"

#include "lib/timekit.hpp"
#include "src/timings.h"

CommonData::CommonData(ModelSettings  * model_settings,
                       InputFiles     * input_files)
:
  outer_temp_simbox_(false),
  read_seismic_(false),
  read_wells_(false),
  block_wells_(false),
  setup_reflection_matrix_(false),
  optimize_well_location_(false),
  wavelet_estimation_shape_(false),
  prior_corr_estimation_(false),
  setup_estimation_rock_physics_(false)
{
  bool failed = false;
  std::string err_text = "";
  //if(readSeismicData(modelSettings,
  //                   inputFiles) == true)
  //  read_seismic_ = true; //True or false if there is no seismic data?

  // 1. set up outer simbox
  outer_temp_simbox_ = CreateOuterTemporarySimbox(model_settings, input_files, estimation_simbox_, full_inversion_volume_, err_text);
  failed = !outer_temp_simbox_;

  // 2. read seismic data
  if(ReadSeismicData(model_settings, input_files) == true)
    read_seismic_ = true; //True or false if there is no seismic data?

  // 3. read well data
  if(!failed){
    read_wells_ = ReadWellData(model_settings, input_files, err_text);
    failed = !read_wells_;
  }

  // 4. block wells for estimation
  // if well position is to be optimized or
  // if wavelet/noise should be estimated or
  // if correlations should be estimated
  if(!failed){
    if (model_settings->getOptimizeWellLocation() || model_settings->getEstimateWaveletNoise() || model_settings->getEstimateCorrelations()){
      block_wells_ = BlockWellsForEstimation(model_settings, estimation_simbox_, wells_, blocked_logs_common_, err_text);
      failed = !block_wells_;
    }
  }

  if(SetupReflectionMatrixAndTempWavelet(model_settings,
                                      input_files) == true) {
    setup_reflection_matrix_ = true; //Set inside function?
    temporary_wavelet_ = true;
  }

  if(WaveletHandling(model_settings, input_files) == true)
    wavelet_handling_ = true;

  //if(readWellData(model_settings,
  //                input_files) == true)
  //  read_wells_ = true;

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
    // The geometry is already present in modelSettings (geometry_ or geometry_full_ ? )
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
    SetSurfacesMultipleIntervals(estimation_simbox, full_inversion_volume, input_files, model_settings, err_text, failed);
  }
  // single interval described by either one or two surfaces
  else{
    SetSurfacesSingleInterval(estimation_simbox, full_inversion_volume, input_files->getTimeSurfFiles(), model_settings, err_text, failed);
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

  int n_timeLapses = model_settings->getNumberOfTimeLapses();
  int error = 0;

  for(int thisTimeLapse = 0; thisTimeLapse < n_timeLapses; thisTimeLapse++) {

    if(input_files->getNumberOfSeismicFiles(thisTimeLapse) > 0 ) {

      LogKit::WriteHeader("Reading seismic data");


      std::vector<float> angles = model_settings->getAngle(thisTimeLapse);
      std::vector<float> offset = model_settings->getLocalSegyOffset(thisTimeLapse);

      int n_angles = model_settings->getNumberOfAngles(thisTimeLapse);
      std::vector<SeismicStorage> seismic_data_angle;

      for (int i = 0 ; i < n_angles; i++) {

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

          ///H
          //CheckThatDataCoverGrid used a TimeCutSimBox. Changed to full_inversion_volume.
          //"Bruk full_inversion_volume i kallet til som klippevolum i kallet til SegY-leseren."
          if(CheckThatDataCoverGrid(segy, offset[i], guard_zone) == true) {

            float padding         = 2*guard_zone;
            bool  relativePadding = false;
            bool  onlyVolume      = true;

            segy->ReadAllTraces(&full_inversion_volume_, //was timeCutSimBox.
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
  } //n_timeLapses

  if(error > 0)
    return false;

  return true;

}

bool
CommonData::CheckThatDataCoverGrid(const SegY   * segy,
                                   float          offset,
                                   //const Simbox * timeCutSimbox,
                                   //const NRLib::Volume   full_inversion_volume,
                                   float          guard_zone)
{
  // Seismic data coverage (translate to CRAVA grid by adding half a grid cell)
  float dz = segy->GetDz();
  float z0 = offset + 0.5f*dz;
  float zn = z0 + (segy->GetNz() - 1)*dz;

  // Top and base of interval of interest
  //float top_grid = static_cast<float>(timeCutSimbox->getTopZMin());
  //float bot_grid = static_cast<float>(timeCutSimbox->getBotZMax());

  //Use this instead? input for timeCutSimbox is estimation_simbox_
  int nx = estimation_simbox_.getnx();
  int ny = estimation_simbox_.getny();
  float top_grid = static_cast<float>(full_inversion_volume_.GetTopZMin(nx, ny));
  float bot_grid = static_cast<float>(full_inversion_volume_.GetBotZMax(nx, ny));

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


bool CommonData::ReadWellData(ModelSettings  * model_settings,
                              InputFiles     * input_files,
                              std::string    & err_text)
{
  bool failed = false;

  int nWells = model_settings->getNumberOfWells();
  try{
    if(nWells > 0)
      LogKit::WriteHeader("Reading wells");

    for(int i=0 ; i<nWells; i++) {

      std::string well_file_name = input_files->getWellFile(i);
      bool read_ok = false;
      NRLib::Well new_well(well_file_name, read_ok);

      if(well_file_name.find(".nwh",0) != std::string::npos)
        ProcessLogsNorsarWell(new_well, err_text, failed);
      else if(well_file_name.find(".rms",0) != std::string::npos)
        ProcessLogsRMSWell(new_well, err_text, failed);

      if(read_ok == true)
        wells_.push_back(new_well);
      else
        LogKit::LogFormatted(LogKit::Error, "Well format of file " + well_file_name + " not recognized.\n");

    }
  }catch (NRLib::Exception & e) {
    err_text += "Error: " + NRLib::ToString(e.what());
    failed = true;
  }

  return !failed;
}

void CommonData::ProcessLogsNorsarWell(NRLib::Well   & new_well,
                                       std::string   & error_text,
                                       bool          & failed){
  const int factor_kilometer = 1000;

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

  if(!new_well.HasContLog("TWT")){
    failed = true;
    error_text += "Could not find log 'TWT' in well file "+new_well.GetWellName()+".\n";
  }

  if(new_well.HasContLog("VP")){
    new_well.AddContLog("Vp", new_well.GetContLog("VP"));
    new_well.RemoveContLog("VP");
  }

  if(new_well.HasContLog("VS")){
    new_well.AddContLog("Vs", new_well.GetContLog("VS"));
    new_well.RemoveContLog("VS");
  }

  if(new_well.HasContLog("RHO")){
    new_well.AddContLog("Rho", new_well.GetContLog("RHO"));
    new_well.RemoveContLog("RHO");
  }

}

void CommonData::ProcessLogsRMSWell(NRLib::Well   & new_well,
                                    std::string   & error_text,
                                    bool          & failed){

  const double factor_usfeet_to_meters = 304800.0;

  if(new_well.HasContLog("x")){
    new_well.AddContLog("X_pos", new_well.GetContLog("x"));
    new_well.RemoveContLog("x");
  }else{
    failed = true;
    error_text += "Could not find log 'x' in well file "+new_well.GetWellName()+".\n";
  }

  if(new_well.HasContLog("y")){
    new_well.AddContLog("Y_pos", new_well.GetContLog("y"));
    new_well.RemoveContLog("y");
  }else{
    failed = true;
    error_text += "Could not find log 'y' in well file "+new_well.GetWellName()+".\n";
  }

  if(new_well.HasContLog("z")){
    new_well.AddContLog("TVD", new_well.GetContLog("z"));
    new_well.RemoveContLog("z");
  }else{
    failed = true;
    error_text += "Could not find log 'z' in well file "+new_well.GetWellName()+".\n";
  }

  if(new_well.HasContLog("DT")){
    std::vector<double> vp_temp = new_well.GetContLog("DT");
    std::vector<double> vp(vp_temp.size());
    for(unsigned int i=0; i<vp_temp.size(); i++){
      vp[i] = static_cast<double>(factor_usfeet_to_meters/vp_temp[i]);
    }
    new_well.RemoveContLog("DT");
    new_well.AddContLog("Vp", vp);
  }

  if(new_well.HasContLog("DTS")){
    std::vector<double> vs_temp = new_well.GetContLog("DTS");
    std::vector<double> vs(vs_temp.size());
    for(unsigned int i=0; i<vs_temp.size(); i++){
      vs[i] = static_cast<double>(factor_usfeet_to_meters/vs_temp[i]);
    }
    new_well.RemoveContLog("DTS");
    new_well.AddContLog("Vs", vs);
  }

  if(new_well.HasContLog("RHOB")){
    new_well.AddContLog("Rho", new_well.GetContLog("RHOB"));
    new_well.RemoveContLog("RHOB");
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
  const std::string & refl_matr_file = input_files->getReflMatrFile();
  const double        vpvs         = model_settings->getVpVsRatio();
  float                  ** reflection_matrix;

  int n_timelapses = model_settings->getNumberOfTimeLapses(); //Returnerer timeLapseAngle_.size()
  for(int i = 0; i < n_timelapses; i++) {

    std::vector<float> angles = model_settings->getAngle(i);
    std::vector<float> offset = model_settings->getLocalSegyOffset(i);

    int n_angles = model_settings->getNumberOfAngles(i);

    if (refl_matr_file != "") {  //File should have one line for each seismic data file. Check: if(input_files->getNumberOfSeismicFiles(thisTimeLapse) > 0 ) ?
      std::string tmpErrText("");
      reflection_matrix = ReadMatrix(refl_matr_file, n_angles, 3, "reflection matrix", tmpErrText);
      if(reflection_matrix == NULL) {
        LogKit::LogFormatted(LogKit::Error, "Reading of file "+refl_matr_file+ " for reflection matrix failed\n");
        LogKit::LogFormatted(LogKit::Error, tmpErrText);
        return false;
      }

      LogKit::LogFormatted(LogKit::Low,"\nReflection parameters read from file.\n\n");
    }

    else if(vpvs != RMISSING) {
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs ratio specified in model file.\n");
      double vsvp = 1.0/vpvs;
      SetupDefaultReflectionMatrix(reflection_matrix, vsvp, model_settings, n_angles, i);
    }
    else {
      LogKit::LogFormatted(LogKit::Low,"\nMaking reflection matrix with Vp/Vs equal to 2\n");
      double vsvp = 1/2;
      SetupDefaultReflectionMatrix(reflection_matrix, vsvp, model_settings, n_angles, i);
    }

      reflection_matrix_[i] = reflection_matrix;
  } //nTimeLapses

  setup_reflection_matrix_ = true;

  //Set up temporary wavelet
  LogKit::WriteHeader("Setting up temporary wavelet");

  //1. Check if optimize welllocation //Point f) Comes from xml-model file
    //2. Check if read seismic ok | read_seismic_ok_
  //read_seismic_ = true;
  if(model_settings->getOptimizeWellLocation() == false && read_seismic_ == true) {

    //3. Use Ricker - wavelet.
    //4. 1 wavelet per angle
    //5  Frequency per ange: Take 100 traces from first AVO-vintage on this angle. Find peak-frequency for these.
    int n_angles = model_settings->getNumberOfAngles(0);
    int error = 0;

    for (int j = 0 ; j < n_angles; j++) {
      //Check all timelapses for this angle, choose the first one;
      int this_timelapse = 0;
      int vintage_year = model_settings->getVintageYear(0);
      int vintage_month = model_settings->getVintageMonth(0);
      int vintage_day = model_settings->getVintageDay(0);
      for(int k = 1; k < n_timelapses; k++) {
        if(model_settings->getVintageYear(k) <= vintage_year && model_settings->getVintageMonth(k) <= vintage_month && model_settings->getVintageDay(k) <= vintage_day) {
          vintage_year = model_settings->getVintageYear(k);
          vintage_month = model_settings->getVintageMonth(k);
          vintage_day = model_settings->getVintageDay(k);
          this_timelapse = k;
        }
      }
      std::vector<float> angles = model_settings->getAngle(this_timelapse);

      std::vector<float> frequency_peaks;
      std::vector<std::vector<float> > trace_data;
      std::vector<float> trace_length;
      seismic_data_[this_timelapse][j].GetSparseTraceData(trace_data, trace_length, 100);

      //FFT to find peak-frequency.
      for(int k = 0; k < 100; k++) {
        int n_trace = trace_data[k].size();
        int n_trace_fft = ((n_trace / 2) + 1)*2;

        fftw_real * seis_r = new fftw_real[n_trace_fft];
        fftw_complex * seis_c = reinterpret_cast<fftw_complex*>(seis_r);

        for(int kk=0; kk < n_trace; kk++)
          seis_r[kk] = trace_data[k][kk];

        Utils::fft(seis_r, seis_c, n_trace);

        std::vector<float> tmp_re_test;
        std::vector<float> tmp_re;
        std::vector<float> tmp_im;

        //for(int iii = 0; iii < n_trace_fft; iii++) {
        //  tmp_re_test.push_back(seis_r[iii]);
        //  tmp_re.push_back(seis_c[iii].re);
        //  tmp_im.push_back(seis_c[iii].im);
        //}

        float peak_tmp = 0.0;
        float length_tmp = 0.0;
        for(int kk = 0; kk < n_trace; kk++) {
          //fftw_complex value = seis_c[kk];
          length_tmp = std::sqrt( (seis_c[kk].re)*(seis_c[kk].re) + (seis_c[kk].im)*(seis_c[kk].im) );
          if(length_tmp > peak_tmp && length_tmp != std::numeric_limits<float>::infinity())
              peak_tmp = length_tmp;
        }

        peak_tmp /= trace_length[k];
        frequency_peaks.push_back(peak_tmp);

        if(seis_r != NULL)
          delete seis_r;
        //if(seis_c != NULL)
        //  delete seis_c;
      }

      float mean_frequency = 0.0;
      float mean_grid_height = 0.0;
      for(size_t k = 0; k < frequency_peaks.size(); k++) {
        mean_frequency += frequency_peaks[k];
        mean_grid_height += trace_length[k];
      }
      mean_frequency /= frequency_peaks.size();
      mean_grid_height /= frequency_peaks.size(); //Will be 1 if SEGY is used.

      mean_frequency /= mean_grid_height;

      int tmp_error = 0;
      Wavelet  * wavelet_tmp = new Wavelet1D(model_settings, reflection_matrix_[this_timelapse][j], angles[j], mean_frequency, tmp_error);

      error += tmp_error;

      if(tmp_error == 0)
        temporary_wavelets_.push_back(wavelet_tmp);
      else
        LogKit::LogFormatted(LogKit::Error, "Error setting up a temporary wavelet for angle " + NRLib::ToString(angles[j]) + ".\n");

    }

    if(error == 0)
      temporary_wavelet_ = true;
    else
      return false;
  }

  return true;
}


float **
CommonData::ReadMatrix(const std::string & file_name,
                       int                 n1,
                       int                 n2,
                       const std::string & read_reason,
                       std::string       & err_text)
{
  float * tmp_res = new float[n1*n2+1];
  std::ifstream inFile;
  NRLib::OpenRead(inFile,file_name);
  std::string text = "Reading "+read_reason+" from file "+file_name+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
  std::string storage;
  int index = 0;
  int error = 0;


  while(error == 0 && inFile >> storage) {
    if(index < n1*n2) {
      try {
        tmp_res[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        err_text += "Error in "+file_name+"\n";
        err_text += e.what();
        error = 1;
      }
    }
    index++;
  }
  if(error == 0) {
    if(index != n1*n2) {
      error = 1;
      err_text += "Found "+NRLib::ToString(index)+" in file "+file_name+", expected "+NRLib::ToString(n1*n2)+".\n";
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
        result[i][j] = tmp_res[index];
        index++;
      }
    }
  }
  else
    LogKit::LogFormatted(LogKit::Low,"failed.\n");
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
  int      i;
  float ** A      = new float * [n_angles];


  double           vsvp2       = vsvp*vsvp;
  std::vector<int> seismicType = model_settings->getSeismicType(this_timelapse);
  std::vector<float>        angles = model_settings->getAngle(this_timelapse);


  for(i = 0; i < n_angles; i++)
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
  reflection_matrix = A;
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
                                 InputFiles    * input_files)
{
  int n_timeLapses = model_settings->getNumberOfTimeLapses();
  int error = 0;

  std::string err_text("");
  std::vector<Surface *> wavelet_estim_interval;
  FindWaveletEstimationInterval(input_files, wavelet_estim_interval, err_text);

  if(err_text != "")
    LogKit::LogFormatted(LogKit::Error, "Error when finding wavelet estimation interval: " + err_text + "\n");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  for(int i = 0; i < n_timeLapses; i++) {

    Wavelet ** wavelet;               ///< Wavelet for angle

    int n_angles = model_settings->getNumberOfAngles(i);

    std::vector<float> sn_ratio = model_settings->getSNRatio(i);
    std::vector<float> angles = model_settings->getAngle(i);
    bool use_local_noise = model_settings->getUseLocalNoise(i);

    //Fra ModelAvoDynamic::processSeismic:
    std::vector<bool> estimate_wavelets = model_settings->getEstimateWavelet(i);

    //Estimation of a wavelet requires the reading of seismic, reading of wells and reflection matrix to be ok. Check blocking of wells since they are used in estimation.
    bool estimate_failed = false;
    for(size_t j = 0; j < estimate_wavelets.size(); j++) {
      if(estimate_wavelets[j])
        if(read_seismic_ == false || read_wells_ == false || setup_reflection_matrix_ == false || block_wells_ == false) {
          estimate_failed = true;
          return false; //Utskrift av feilmelding?
        }
    }

    std::vector<bool> use_ricker_wavelet = model_settings->getUseRickerWavelet(i);

    wavelet = new Wavelet * [n_angles];
    std::vector<Grid2D *> local_noise_scale;       ///< Scale factors for local noise
    std::vector<Grid2D *> local_shift;
    std::vector<Grid2D *> local_scale;
    local_noise_scale.resize(n_angles);
    local_shift.resize(n_angles);
    local_scale.resize(n_angles);
    bool has_3D_wavelet = false;

    for(int j=0; j <n_angles; j++) {

      local_noise_scale[j] = NULL;
      local_shift[j] = NULL;
      local_scale[j] = NULL;

      if (model_settings->getWaveletDim(j) == Wavelet::THREE_D)
        has_3D_wavelet = true;
      if(estimate_wavelets[j] == true)
        model_settings->setWaveletScale(i,j,1.0);
    }

    unsigned int                      n_wells = model_settings->getNumberOfWells();

    std::vector<std::vector<double> > t_grad_x(n_wells);
    std::vector<std::vector<double> > t_grad_y(n_wells);

    NRLib::Grid2D<float>              ref_time_grad_x;         ///< Time gradient in x-direction for reference time surface (t0)
    NRLib::Grid2D<float>              ref_time_grad_y;         ///< Time gradient in x-direction for reference time surface (t0)
    NRLib::Grid2D<float>              structure_depth_grad_x;  ///< Depth gradient in x-direction for structure ( correlationDirection-t0)*v0/2
    NRLib::Grid2D<float>              structure_depth_grad_y;  ///< Depth gradient in y-direction for structure ( correlationDirection-t0)*v0/2

    bool failed = false;

    if (has_3D_wavelet) {
      if (input_files->getRefSurfaceFile() != "") {
        Surface  t0_surf;
        try {
          t0_surf =Surface(input_files->getRefSurfaceFile());
        }
        catch (NRLib::Exception & e) {
          err_text += e.what();
          failed = true;
        }
         if(!failed)
        {
          //Correlation direction from ModelGeneral::makeTimeSimboxes
          Surface * correlation_direction = NULL;
          try {
            Surface tmpSurf(input_files->getCorrDirFile()); //Top and base?
            if(estimation_simbox_.CheckSurface(tmpSurf) == true)
              correlation_direction = new Surface(tmpSurf);
            else {
              err_text += "Error: Correlation surface does not cover volume.\n"; //?
              failed = true;
            }
          }
          catch (NRLib::Exception & e) {
            err_text += e.what();
            failed = true;
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
                                       ref_time_grad_x,
                                       ref_time_grad_y);
        }
        else{
          err_text += "Problems reading reference time surface in (x,y).\n";
          error = 1;
        }
      }
      bool estimateWellGradient = model_settings->getEstimateWellGradientFromSeismic();
      float distance, sigma_m;
      model_settings->getTimeGradientSettings(distance, sigma_m, i);
      std::vector<std::vector<double> > SigmaXY;

      for(size_t w = 0; w < n_wells; w++) {
        if(!estimateWellGradient & ((structure_depth_grad_x.GetN()> 0) & (structure_depth_grad_y.GetN()>0))){
          double v0=model_settings->getAverageVelocity();
          blocked_logs_common_[w]->SetSeismicGradient(v0, structure_depth_grad_x, structure_depth_grad_y, ref_time_grad_x, ref_time_grad_y, t_grad_x[w], t_grad_y[w]);
        }else{
          blocked_logs_common_[w]->SetTimeGradientSettings(distance, sigma_m);
          blocked_logs_common_[w]->FindSeismicGradient(seismic_data_[i], &estimation_simbox_, n_angles, t_grad_x[w], t_grad_y[w], SigmaXY);
        }
      }
    }

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
    bool local_noise_set = false;

    for (int j = 0; j < n_angles; j++) {
      seismic_data_[i][j];

      float angle = float(angles[i]*180.0/M_PI);
      LogKit::LogFormatted(LogKit::Low,"\nAngle stack : %.1f deg",angle);
      //if(model_settings->getForwardModeling()==false)
      //  seisCube[i]->setAccessMode(FFTGrid::RANDOMACCESS);
      if (model_settings->getWaveletDim(j) == Wavelet::ONE_D)
        error += Process1DWavelet(model_settings,
                                  input_files,
                                  &seismic_data_[i][j],
                                  blocked_logs_common_,
                                  wavelet_estim_interval,
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
                                  use_local_noise);
      else
        error += Process3DWavelet(model_settings,
                                  input_files,
                                  &seismic_data_[i][j],
                                  blocked_logs_common_,
                                  wavelet_estim_interval,
                                  err_text,
                                  wavelet[j],
                                  i, //Timelapse
                                  j, //Angle
                                  angles[j],
                                  sn_ratio[j],
                                  ref_time_grad_x,
                                  ref_time_grad_y,
                                  t_grad_x,
                                  t_grad_y,
                                  estimate_wavelets[j]);

      if(local_noise_scale[j] != NULL)
        local_noise_set = true;
      //if(modelSettings->getForwardModeling()==false) // else, no seismic data
      //  seisCube[i]->endAccess();

    } //angle

    if(local_noise_set == true) {
      for(int i=0;i<n_angles;i++)
        if(local_noise_scale[i]==NULL)
          local_noise_scale[i] = new Grid2D(estimation_simbox_.getnx(),
                                            estimation_simbox_.getny(),
                                            1.0);
    }

    wavelets_[i] = wavelet;
    local_noise_scale_[i] = local_noise_scale;
    local_shift_[i] = local_shift;
    local_scale_[i] = local_scale;
    global_noise_estimate_[i] = sn_ratio;
    sn_ratio_[i] = sn_ratio;
  }  //timelapse

  Timings::setTimeWavelets(wall,cpu);

  return true;

}

int
CommonData::Process1DWavelet(const ModelSettings          * model_settings,
                             const InputFiles             * input_files,
                             const SeismicStorage           * seismic_data,
                             std::vector<BlockedLogsCommon *> blocked_logs,
                             const std::vector<Surface *> & wavelet_estim_interval,
                             std::string                  & err_text,
                             Wavelet                     *& wavelet,
                             Grid2D                      *& local_noise_scale, //local noise estimates?
                             Grid2D                      *& local_shift,
                             Grid2D                      *& local_scale,
                             unsigned int                   i_timelapse,
                             unsigned int                   j_angle,
                             const float                    angle,
                             float                          sn_ratio,
                             bool                           estimate_wavelet,
                             bool                           use_ricker_wavelet,
                             bool                           use_local_noise)
{

  int error = 0;
  //Grid2D * shiftGrid(NULL);
  //Grid2D * gainGrid(NULL);
  if(model_settings->getUseLocalWavelet() && input_files->getScaleFile(i_timelapse,j_angle) != "") {
      Surface help(input_files->getScaleFile(i_timelapse,j_angle));
      //gainGrid = new Grid2D(estimation_simbox_.getnx(),estimation_simbox_.getny(), 0.0);
      local_scale = new Grid2D(estimation_simbox_.getnx(),estimation_simbox_.getny(), 0.0);
      ResampleSurfaceToGrid2D(&help, local_scale);
  }
  if (model_settings->getUseLocalWavelet() && input_files->getShiftFile(i_timelapse,j_angle) != "") {
    Surface helpShift(input_files->getShiftFile(i_timelapse,j_angle));
    local_shift = new Grid2D(estimation_simbox_.getnx(),estimation_simbox_.getny(), 0.0);
    ResampleSurfaceToGrid2D(&helpShift, local_shift);
  }
  if (use_local_noise && input_files->getLocalNoiseFile(i_timelapse,j_angle) != "") {
    Surface helpNoise(input_files->getLocalNoiseFile(i_timelapse,j_angle));
    local_noise_scale = new Grid2D(estimation_simbox_.getnx(), estimation_simbox_.getny(), 0.0);
    ResampleSurfaceToGrid2D(&helpNoise, local_noise_scale);
  }

  if (estimate_wavelet)
    wavelet = new Wavelet1D(&estimation_simbox_,
                            seismic_data,
                            blocked_logs,
                            wavelet_estim_interval,
                            model_settings,
                            *reflection_matrix_[i_timelapse],
                            j_angle,
                            error,
                            err_text);

  else { //Not estimation modus
    if(use_ricker_wavelet)
        wavelet = new Wavelet1D(model_settings,
                                *reflection_matrix_[i_timelapse],
                                angle,
                                model_settings->getRickerPeakFrequency(i_timelapse,j_angle),
                                error);
    else {
      const std::string & waveletFile = input_files->getWaveletFile(i_timelapse,j_angle);
      int fileFormat = GetWaveletFileFormat(waveletFile, err_text);
      if(fileFormat < 0) {
        err_text += "Unknown file format of file '"+waveletFile+"'.\n";
        error++;
      }
      else
        wavelet = new Wavelet1D(waveletFile,
                                fileFormat,
                                model_settings,
                                *reflection_matrix_[i_timelapse],
                                angle,
                                error,
                                err_text);
    }
      // Calculate a preliminary scale factor to see if wavelet is in the same size order as the data. A large or small value might cause problems.
      if(seismic_data->GetFileName() != "") {
      //if(seisCube!=NULL) {// If forward modeling, we have no seismic, can not prescale wavelet.
        float       prescale = wavelet->findGlobalScaleForGivenWavelet(model_settings, &estimation_simbox_, seismic_data, blocked_logs);
        const float limHigh  = 3.0f;
        const float limLow   = 0.33f;

        if(model_settings->getEstimateGlobalWaveletScale(i_timelapse,j_angle)) // prescale, then we have correct size order, and later scale estimation will be ok.
           wavelet->multiplyRAmpByConstant(prescale);
        else {
          if(model_settings->getWaveletScale(i_timelapse,j_angle)!= 1.0f && (prescale>limHigh || prescale<limLow)) {
             std::string text = "The wavelet given for angle no "+NRLib::ToString(j_angle)+" is badly scaled. Ask Crava to estimate global wavelet scale.\n";
            if(model_settings->getEstimateLocalScale(i_timelapse,j_angle)) {
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
      if (error == 0)
        wavelet->resample(static_cast<float>(estimation_simbox_.getdz()),
                          estimation_simbox_.getnz(),
                          model_settings->getNZpad());
  }

  if (error == 0) {
    wavelet->scale(model_settings->getWaveletScale(i_timelapse,j_angle));

    if (model_settings->getForwardModeling() == false && model_settings->getNumberOfWells() > 0) {
      float SNRatio_tmp = wavelet->calculateSNRatioAndLocalWavelet(&estimation_simbox_,
                                                                   seismic_data,
                                                                   blocked_logs,
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
                                                                   model_settings->getWaveletScale(i_timelapse,j_angle),
                                                                   model_settings->getEstimateSNRatio(i_timelapse,j_angle),
                                                                   model_settings->getEstimateGlobalWaveletScale(i_timelapse,j_angle),
                                                                   model_settings->getEstimateLocalNoise(i_timelapse,j_angle),
                                                                   model_settings->getEstimateLocalShift(i_timelapse,j_angle),
                                                                   model_settings->getEstimateLocalScale(i_timelapse,j_angle),
                                                                   estimate_wavelet);
      if(model_settings->getEstimateSNRatio(i_timelapse,j_angle))
        sn_ratio = SNRatio_tmp;
        //SNRatio_[i_timelapse] = SNRatio;
    }

    if (error == 0) {
      if((model_settings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
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
      if ((sn_ratio <=SNLow  || sn_ratio > SNHigh) && model_settings->getForwardModeling()==false) {
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
  return error;
}

int
CommonData::Process3DWavelet(const ModelSettings                     * model_settings,
                             const InputFiles                        * input_files,
                             const SeismicStorage                    * seismic_data,
                             std::vector<BlockedLogsCommon *>          blocked_logs,
                             const std::vector<Surface *>            & wavelet_estim_interval,
                             std::string                             & err_text,
                             Wavelet                                *& wavelet,
                             unsigned int                              i_timelapse,
                             unsigned int                              j_angle,
                             float                                     angle,
                             float                                     sn_ratio,
                             const NRLib::Grid2D<float>              & ref_time_grad_x,
                             const NRLib::Grid2D<float>              & ref_time_grad_y,
                             const std::vector<std::vector<double> > & t_grad_x,
                             const std::vector<std::vector<double> > & t_grad_y,
                             bool                                      estimate_wavelet)
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
                            blocked_logs,
                            &estimation_simbox_,
                            *reflection_matrix_[i_timelapse],
                            j_angle,
                            error,
                            err_text);
  }
  else { //Not estimation modus
    const std::string & wavelet_file = input_files->getWaveletFile(i_timelapse,j_angle);
    int file_format = GetWaveletFileFormat(wavelet_file, err_text);
    if(file_format < 0) {
      err_text += "Unknown file format of file '"+wavelet_file+"'.\n";
      error++;
    }
    else {
      wavelet = new Wavelet3D(wavelet_file,
                              file_format,
                              model_settings,
                              *reflection_matrix_[i_timelapse],
                              angle,
                              error,
                              err_text,
                              input_files->getWaveletFilterFile(j_angle));
      if (error == 0)
        wavelet->resample(static_cast<float>(estimation_simbox_.getdz()),
                          estimation_simbox_.getnz(),
                          model_settings->getNZpad());
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

    if (localEst && model_settings->getForwardModeling() == false) {
      float sn_ratio_tmp = wavelet->calculateSNRatio(&estimation_simbox_,
                                                     seismic_data,
                                                     blocked_logs,
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
      if(model_settings->getEstimateSNRatio(i_timelapse,j_angle))
        sn_ratio = sn_ratio_tmp;

    }
    if (error == 0) {
      if((model_settings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
         (model_settings->getEstimationMode() && estimate_wavelet)) {
        std::string type;
        if (estimate_wavelet) {
          type = "Estimated_";
          if(wavelet->getDim()==1)
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

      const float SNLow  = 1.0;
      const float SNHigh = 10.0;
      if ((sn_ratio <=SNLow  || sn_ratio > SNHigh) && model_settings->getForwardModeling()==false) {
        err_text += "Illegal signal-to-noise ratio of "+NRLib::ToString(sn_ratio)+" for cube "+NRLib::ToString(j_angle+1)+".\n";
        err_text += "Ratio must be in interval "+NRLib::ToString(SNLow)+" < S/N ratio < "+NRLib::ToString(SNHigh)+"\n";
        error++;
      }
    }
  }

  return error;
}

void
CommonData::FindWaveletEstimationInterval(InputFiles             * input_files,
                                          std::vector<Surface *> & wavelet_estim_interval,
                                          std::string &            err_text)

{
  const double x0 = estimation_simbox_.getx0();
  const double y0 = estimation_simbox_.gety0();
  const double lx = estimation_simbox_.getlx();
  const double ly = estimation_simbox_.getly();
  const int    nx = estimation_simbox_.getnx();
  const int    ny = estimation_simbox_.getny();

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
      //failed = true;
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
      //failed = true;
    }
  }

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

   for(int i=0;i<nx;i++)
     for(int j=0;j<ny;j++)
     {
       double x,y;
       double gx,gy,gxTmp,gyTmp;
       gx=0.0;
       gy=0.0;
       estimation_simbox_.getXYCoord(i,j,x,y);
       CalculateSmoothGrad( t0_surf, x, y, radius, ds,gxTmp, gyTmp);
       gx=-gxTmp;
       gy=-gyTmp;
       if( correlation_direction !=NULL){
         CalculateSmoothGrad( correlation_direction, x, y, radius, ds,gxTmp, gyTmp);
         gx+=gxTmp;
         gy+=gyTmp;
       }else{
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

void
CommonData::ComputeReferenceTimeGradient(const Surface * t0_surf,
                                         NRLib::Grid2D<float> &ref_time_grad_x,
                                         NRLib::Grid2D<float> &ref_time_grad_y)
 {
   double radius = 50.0;
   double ds = 12.5;
   int nx = estimation_simbox_.getnx();
   int ny = estimation_simbox_.getny();
   ref_time_grad_x.Resize(nx,ny);
   ref_time_grad_y.Resize(nx,ny);
   for(int i=0;i<nx;i++)
     for(int j=0;j<ny;j++)
     {
       double x,y;
       double gx,gy;
       estimation_simbox_.getXYCoord(i,j,x,y);
       CalculateSmoothGrad( t0_surf, x, y, radius, ds,gx, gy);
       ref_time_grad_x(i,j) =float(gx);
       ref_time_grad_y(i,j) =float(gy);
     }
 }

void
CommonData::CalculateSmoothGrad(const Surface * surf, double x, double y, double radius, double ds,  double& gx, double& gy)
{
  /// Return smoothed Gradient. Computes the gradient as a regression
  /// among points within a given distance from the central point.
  //  Returns missing if central point is outside the grid.
  /// Returns  otherwise it is ok,
  int i,j,k,l;
  int disc_radius = int(floor(radius/ds));
  int n_points    = (2*disc_radius+1)*(2*disc_radius+1);
  std::vector<double> Z(3*n_points,0.0);
  std::vector<double> Y(n_points,0.0);
  std::vector<double> cov(9,0.0);
  std::vector<double> invCov(9,0.0);
  std::vector<double> proj(3,0.0);
  double z0;
  int cy=0;
  int cz=0;

  double baseDepth =surf->GetZ(x,y);


  for(i=- disc_radius;i<=disc_radius;i++)
    for( j=- disc_radius;j<=disc_radius;j++)
    {
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

  int nData=0;
  for(i=0; i < n_points; i++)
  {
    if(!surf->IsMissing(Y[i]))
    {
      nData++;
      for(k=0;k<3;k++)
      {
        for(l=0;l<3;l++)
          cov[k+l*3]+=Z[k + 3*i] * Z[l + 3*i];

        proj[k]+=Z[k + 3*i]*(Y[i]-baseDepth);
      }
    }
  }

  double det = cov[0]*(cov[4]*cov[8] - cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8] - cov[5]*cov[6])
                  +   cov[2]*(cov[3]*cov[7] - cov[4]*cov[6]);

  if(det != 0)
  {
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
      for(k=0;k<3;k++)
      {
        z0 += invCov[k]*proj[k]; //NBNB check
        gx += invCov[3+k]*proj[k];
        gy += invCov[6+k]*proj[k];
      }
  }
  else
  {
   gx = RMISSING;
   gy = RMISSING;
  }
}





void
CommonData::ResampleSurfaceToGrid2D(const Surface * surface,
                                    Grid2D        * outgrid)
{
  for(int i=0;i<estimation_simbox_.getnx();i++) {
    for(int j=0;j<estimation_simbox_.getny();j++) {
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

  if(fileformat<0) { // not norsar format
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
        if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
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
                                              //const unsigned int    i,
                                              //const Simbox        * timeSimbox,
                                              const Grid2D        * grid,
                                              const float           angle)
{
  bool   estimationMode   = model_settings->getEstimationMode();
  int    outputFormat     = model_settings->getOutputGridFormat();
  double angle_tmp        = angle*180.0/M_PI;

  Surface * help = NULL;

  if(file_name != "") {
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
  if(angle > -45 || angle < 45)
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
  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) {
      surface->GetXY(i,j,x,y);
      simbox->getIndexes(x,y,i1,j1);
      if(i1==IMISSING || j1== IMISSING)
        surface->SetMissing(i,j);
      else
        (*surface)(i,j) = (*grid)(i1,j1);
    }
  }
}




bool CommonData::optimizeWellLocations() {
  return true;
}

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
  double azimuth = (-1)*area_rot*(180.0/M_PI);
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

void CommonData::SetSurfacesSingleInterval(Simbox                           & estimation_simbox,
                                           NRLib::Volume                    & full_inversion_volume,
                                           const std::vector<std::string>   & surf_file,
                                           ModelSettings                    * model_settings,
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
  int    output_grids_seismic = model_settings->getOutputGridsSeismic();
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

void CommonData::SetSurfacesMultipleIntervals(Simbox                         & estimation_simbox,
                                              NRLib::Volume                  & full_inversion_volume,
                                              const InputFiles               * input_files,
                                              const ModelSettings            * model_settings,
                                              std::string                    & err_text,
                                              bool                           & failed){

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
                                         //const InputFiles                               * const input_files,
                                         const Simbox                                   & estimation_simbox,
                                         const std::vector<NRLib::Well>                 & wells,
                                         std::vector<BlockedLogsCommon *>               & blocked_logs_common,
                                         std::string                                    & err_text){
  bool failed = false;

  try{
    unsigned int n_wells = wells.size();

    for (unsigned int i=0; i<n_wells; i++){
      BlockedLogsCommon * blocked_log = new BlockedLogsCommon(&wells[i], &estimation_simbox, model_settings->getRunFromPanel());
      blocked_logs_common.push_back(blocked_log);
    }
  }catch(NRLib::Exception & e){
    err_text += e.what();
    failed = true;
  }

  return failed;
}