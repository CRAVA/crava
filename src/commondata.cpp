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
#include "src/simbox.h"
#include "src/timeline.h"
#include "src/tasklist.h"
#include "src/modelgeneral.h"

#include "nrlib/segy/segy.hpp"

CommonData::CommonData(ModelSettings  * model_settings,
                       InputFiles     * input_files){
  estimation_simbox_ = Simbox();
  std::string err_text = "";
  createOuterTemporarySimbox(model_settings, input_files, estimation_simbox_, full_inversion_volume_, err_text);

}

CommonData::~CommonData(){
  //delete estimation_simbox_;
  //delete full_inversion_volume_;
}

bool CommonData::createOuterTemporarySimbox(ModelSettings   * model_settings,
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

  else if(area_specification == ModelSettings::AREA_FROM_GRID_DATA         || // tilfelle iii 
          area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_UTM || // tilfelle i med snap to seismic
          area_specification == ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE){ // tilfelle ii med snap to seismic
    LogKit::LogFormatted(LogKit::High,"\nFinding inversion area from grid data in file \'"+grid_file+"\'\n");
    area_type = "Grid data";
    std::string tmp_err_text;
    SegyGeometry * geometry;
    getGeometryFromGridOnFile(grid_file,
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
      writeAreas(area_params,&estimation_simbox_,area_type);
      err_text += "The specified AREA extends outside the surface(s).\n";
    }
  }

  // SET TOP AND BASE SURFACES FOR THE FORWARD/INVERSION INTERVALS -----------------------------------------------

  // if multiple intervals
  if(model_settings->getIntervalNames().size() > 0){
    setSurfacesMultipleIntervals(estimation_simbox, full_inversion_volume, input_files, model_settings, err_text, failed);
  }
  // single interval described by either one or two surfaces
  else{
    setSurfacesSingleInterval(estimation_simbox, full_inversion_volume, input_files->getTimeSurfFiles(), model_settings, err_text, failed);
  }

  return (!failed);
}

bool CommonData::readSeismicData(){
  return true;

}

bool CommonData::readWellData(){
  return true;

}

bool CommonData::blockWellsForEstimation(){
  return true;

}

bool CommonData::setupReflectionMatrixAndTempWavelet(){
  return true;

}

bool CommonData::optimizeWellLocations(){
  return true;
}

bool CommonData::estimateWaveletShape(){
  return true;
}
  
bool CommonData::estimatePriorCorrelation(){
  return true;
}
  
bool CommonData::setupEstimationRockPhysics(){
  return true;
}

int
CommonData::computeTime(int year, int month, int day) const
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
CommonData::getGeometryFromGridOnFile(const std::string           grid_file,
                                      const TraceHeaderFormat   * thf,
                                      SegyGeometry             *& geometry,
                                      std::string               & err_text)
{
  geometry = NULL;

  if(grid_file != "") { //May change the condition here, but need geometry if we want to set XL/IL
    int file_type = IO::findGridType(grid_file);
    if(file_type == IO::CRAVA) {
      geometry = geometryFromCravaFile(grid_file);
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
      geometry = geometryFromStormFile(grid_file, err_text);
    else if(file_type==IO::SGRI) {
      bool scale = true;
      geometry = geometryFromStormFile(grid_file, err_text, scale);
    }
    else {
      err_text = "Trying to read grid dimensions from unknown file format.\n";
    }
  }
  else {
    err_text = "Cannot get geometry from file. The file name is empty.\n";
  }
}

SegyGeometry * CommonData::geometryFromCravaFile(const std::string & file_name)
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

SegyGeometry * CommonData::geometryFromStormFile(const std::string & file_name,
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

void CommonData::writeAreas(const SegyGeometry * area_params,
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

  findSmallestSurfaceGeometry(area_x0, area_y0, area_lx, area_ly, area_rot,
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


void CommonData::findSmallestSurfaceGeometry(const double   x0,
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

void CommonData::setSurfacesSingleInterval(Simbox                           & estimation_simbox,
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
    findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
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
          findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
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

void CommonData::setSurfacesMultipleIntervals(Simbox                         & estimation_simbox,
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
      findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
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
        findSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
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
