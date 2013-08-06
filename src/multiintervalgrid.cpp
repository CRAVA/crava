/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#include "src/multiintervalgrid.h"
#include "src/definitions.h"
#include "src/simbox.h"

#include "nrlib/grid/grid.hpp"

MultiIntervalGrid::MultiIntervalGrid(ModelSettings  * model_settings,
                                     InputFiles     * input_files,
                                     const Simbox   * estimation_simbox,
                                     std::string    & err_text,
                                     bool           & failed){

  std::vector<std::string> interval_names = model_settings->getIntervalNames();
  n_intervals_ = static_cast<int>(interval_names.size());
  int erosion_priority_top_surface                                = model_settings->getErosionPriorityTopSurface();
  const std::map<std::string,int> erosion_priority_base_surfaces  = model_settings->getErosionPriorityBaseSurfaces();

  Surface               * top_surface = NULL;
  Surface               * base_surface = NULL;
  std::vector<Surface>    eroded_surfaces(n_intervals_+1);
  std::string             previous_interval_name("");
  std::string             top_surface_file_name_temp("");
  std::string             base_surface_file_name_temp("");
  std::vector<Surface>    surfaces;
  std::vector<int>        erosion_priorities;

  // if there are multiple intervals
  if(n_intervals_ > 0){
    surfaces.resize(n_intervals_+1); //Store surfaces.
    erosion_priorities.resize(n_intervals_+1);
    interval_simboxes_.resize(n_intervals_);
  }
  // if there is only one interval
  else{
    surfaces.resize(1);
    interval_simboxes_.resize(1);
  }

  // 1. ERODE SURFACES -------------------------------------------------------------------

  try{
    top_surface_file_name_temp = input_files->getTimeSurfFile(0);
    erosion_priorities[0] = erosion_priority_top_surface;

    top_surface = MakeSurfaceFromFileName(top_surface_file_name_temp, *estimation_simbox);
    surfaces[0] = *top_surface;

    for(size_t i = 0; i < n_intervals_; i++) {

      std::string interval_name = model_settings->getIntervalName(i);
      base_surface_file_name_temp = input_files->getIntervalBaseTimeSurface(interval_name);
      erosion_priorities[i+1] = erosion_priority_base_surfaces.find(interval_name)->second;

      base_surface = MakeSurfaceFromFileName(base_surface_file_name_temp, *estimation_simbox);
      surfaces[i] =  *base_surface;
    }

    ErodeAllSurfaces(eroded_surfaces,
                     erosion_priorities,
                     surfaces,
                     *estimation_simbox);
  }
  catch(NRLib::Exception & e){
    failed = true;
    err_text += e.what();
  }

  // 2 SET UP INTERVAL_SIMBOXES ----------------------------------------------------------

  //Set up a vector of simboxes, one per interval.

  simboxes_.resize(interval_names.size());
  interval_simboxes_.resize(interval_names.size());

  if(!failed){
    try{
      const std::map<std::string, std::string> corr_dir_single_surfaces = input_files->getCorrDirIntervalFiles();
      const std::map<std::string, std::string> corr_dir_top_surfaces    = input_files->getCorrDirIntervalTopSurfaceFiles();
      const std::map<std::string, std::string> corr_dir_base_surfaces   = input_files->getCorrDirIntervalBaseSurfaceFiles();
      const std::map<std::string, bool> corr_dir_top_conform            = model_settings->getCorrDirIntervalTopConforms();
      const std::map<std::string, bool> corr_dir_base_conform           = model_settings->getCorrDirIntervalBaseConforms();

      SetUpIntervalSimboxes(model_settings,
                            estimation_simbox,
                            interval_names,
                            eroded_surfaces,
                            interval_simboxes_,
                            simboxes_,
                            corr_dir_single_surfaces,
                            corr_dir_top_surfaces,
                            corr_dir_base_surfaces,
                            corr_dir_top_conform,
                            corr_dir_base_conform,
                            err_text,
                            failed);
    }
    catch(NRLib::Exception & e){
      failed = true;
      err_text += e.what();
    }
  }


  // 3. SET UP BACKGROUND MODEL ----------------------------------------------------------

  std::vector<NRLib::Grid<float> > vp_intervals(n_intervals_);
  std::vector<NRLib::Grid<float> > vs_intervals(n_intervals_);
  std::vector<NRLib::Grid<float> > rho_intervals(n_intervals_);

  if(!failed){
    try{

      //dz to vector
      //std::vector<float> dz;
      //for(size_t i = 0; i < n_intervals_; i++)
      //  dz.push_back( interval_simboxes_[i].GetDz()*interval_simboxes_[i].GetAvgRelThick() * 4 );
      //float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

      if(n_intervals_ > 0){ // It is possible to give only one interval in multiple-intervals
        BuildSeismicPropertyIntervals(vp_intervals,
                                      vs_intervals,
                                      rho_intervals,
                                      interval_simboxes_,
                                      relative_grid_resolution_);
      }else{
        // what ?
      }
    }
    catch(NRLib::Exception & e){
    failed = true;
    err_text += e.what();
    }
  }
}

MultiIntervalGrid::~MultiIntervalGrid(){

}

void   MultiIntervalGrid::SetUpIntervalSimboxes(const ModelSettings                       * model_settings,
                                                const Simbox                              * estimation_simbox,
                                                const std::vector<std::string>            & interval_names,
                                                const std::vector<Surface>                & eroded_surfaces,
                                                std::vector<IntervalSimbox>               & interval_simboxes,
                                                std::vector<Simbox>                       & simboxes,
                                                const std::map<std::string, std::string>  & corr_dir_single_surfaces,
                                                const std::map<std::string, std::string>  & corr_dir_top_surfaces,
                                                const std::map<std::string, std::string>  & corr_dir_base_surfaces,
                                                const std::map<std::string, bool>         & corr_dir_top_conform,
                                                const std::map<std::string, bool>         & corr_dir_base_conform,
                                                std::string                               & err_text,
                                                bool                                      & failed) const{

  for (size_t i = 0; i< interval_names.size(); i++){

    std::string            interval_name                          = interval_names[i];
    Surface                top_surface                            = eroded_surfaces[i];
    Surface                base_surface                           = eroded_surfaces[i+1];
    int                    n_layers                               = model_settings->getTimeNzInterval(interval_name);
    std::map<std::string, std::string>::const_iterator it_single  = corr_dir_single_surfaces.find(interval_name);
    std::map<std::string, std::string>::const_iterator it_top     = corr_dir_top_surfaces.find(interval_name);
    std::map<std::string, std::string>::const_iterator it_base    = corr_dir_base_surfaces.find(interval_name);
    std::map<std::string, bool>::const_iterator it_top_conform    = corr_dir_top_conform.find(interval_name);
    std::map<std::string, bool>::const_iterator it_base_conform   = corr_dir_base_conform.find(interval_name);

    // Make a simbox for the original interval --------------------------------------------
    //SegyGeometry * geometry = model_settings->getAreaParameters();
    simboxes[i] = Simbox(estimation_simbox);
    simboxes[i].SetSurfaces(top_surface, base_surface);

    // Make extended interval_simbox for the inversion interval ---------------------------

    // Case 1: Single correlation surface
    if(it_single != corr_dir_single_surfaces.end() && it_top == corr_dir_top_surfaces.end() && it_base == corr_dir_base_surfaces.end()){
      Surface * corr_surf = MakeSurfaceFromFileName(it_single->second,  estimation_simbox);
      interval_simboxes[i] = IntervalSimbox(&simboxes[i], interval_names[i], n_layers,  top_surface, base_surface, corr_surf,
                                            err_text, failed);
    }
    // Case 2: Top and base correlation surfaces
    else if(it_single == corr_dir_single_surfaces.end() && it_top != corr_dir_top_surfaces.end() && it_base != corr_dir_base_surfaces.end()){
      Surface * corr_surf_top = MakeSurfaceFromFileName(it_top->second,  estimation_simbox);
      Surface * corr_surf_base = MakeSurfaceFromFileName(it_base->second,  estimation_simbox);
      interval_simboxes[i] = IntervalSimbox(&simboxes[i], interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                              corr_surf_top, corr_surf_base);
    }
    // Case 3: Top conform and base conform
    else if(it_top_conform == corr_dir_top_conform.end() && it_base_conform == corr_dir_base_conform.end()){
      interval_simboxes[i] = IntervalSimbox(&simboxes[i], interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                              NULL, NULL);
    }
    // Case 4: Top correlation surface and base conform
    else if(it_single == corr_dir_single_surfaces.end() && it_top != corr_dir_top_surfaces.end() && it_base_conform == corr_dir_base_conform.end()){
      Surface * corr_surf_top = MakeSurfaceFromFileName(it_top->second,  estimation_simbox);
      interval_simboxes[i] = IntervalSimbox(&simboxes[i], interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                              corr_surf_top, NULL);
    }
    // Case 5: Top conform and base correlation surface
    else if(it_top_conform == corr_dir_top_conform.end() && it_base_conform != corr_dir_base_conform.end()){
      Surface * corr_surf_base = MakeSurfaceFromFileName(it_base->second,  estimation_simbox);
      interval_simboxes[i] = IntervalSimbox(&simboxes[i], interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                              NULL, corr_surf_base);
    }
    // else something is wrong
    else{
      err_text += "\nCorrelation directions are not set correctly for interval " + interval_name[i];
      err_text += ".\n";
      failed = true;
    }
  }
}

Surface * MultiIntervalGrid::MakeSurfaceFromFileName(const std::string    & file_name,
                                                     const Simbox         & estimation_simbox) const{

  Surface * new_surface = NULL;

  if (!NRLib::IsNumber(file_name)) { // If the file name is a string
    new_surface = new Surface(file_name);
  }
  else { //If the file name is a value

    double x_min, x_max, y_min, y_max;

    FindSmallestSurfaceGeometry(estimation_simbox.getx0(), estimation_simbox.gety0(),
                                estimation_simbox.getlx(), estimation_simbox.getly(),
                                estimation_simbox.getAngle(), x_min, y_min, x_max, y_max);

    new_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(file_name.c_str()));
  }

  return new_surface;
}

void MultiIntervalGrid::FindSmallestSurfaceGeometry(const double   x0,
                                                    const double   y0,
                                                    const double   lx,
                                                    const double   ly,
                                                    const double   rot,
                                                    double       & x_min,
                                                    double       & y_min,
                                                    double       & x_max,
                                                    double       & y_max) const
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


void  MultiIntervalGrid::ErodeAllSurfaces(std::vector<Surface>            & eroded_surfaces,
                                          const std::vector<int>          & erosion_priorities,
                                          const std::vector<Surface>      & surfaces,
                                          const Simbox                    & simbox) const{
  int    n_surf     = static_cast<int>(eroded_surfaces.size());

  for(int i=0; i<n_surf; i++) {
    int l=0;
    while(i+1 != erosion_priorities[l])
      l++;

    Surface temp_surface = Surface(surfaces[l]);

    //Find closest eroded surface downward
    for(int k=l+1; k<n_surf; k++) {
      if(eroded_surfaces[k].GetN() > 0) {
        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, false);
        break;
      }
    }
    //Find closest eroded surface upward
    for(int k=l-1; k>=0; k--) {
      if(eroded_surfaces[k].GetN() > 0) {
        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, true);
        break;
      }
    }
    eroded_surfaces[l] = temp_surface;
  }
}

void  MultiIntervalGrid::ErodeSurface(Surface       &  surface,
                                      const Surface &  priority_surface,
                                      const Simbox  &  simbox,
                                      const bool    &  compare_upward) const{
  int nx    = simbox.getnx();
  int ny    = simbox.getny();
  double x0 = simbox.GetXMin();
  double y0 = simbox.GetYMin();
  double lx = simbox.GetLX();
  double ly = simbox.GetLY();

  NRLib::Grid2D<double> eroded_surface(nx,ny,0);
  double x;
  double y;
  double z;
  double z_priority;

  double missing = surface.GetMissingValue();
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      simbox.getXYCoord(i,j,x,y);

      z_priority = priority_surface.GetZ(x,y);
      z          = surface.GetZ(x,y);

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


  surface = Surface(x0, y0, lx, ly, eroded_surface);
}

void MultiIntervalGrid::BuildSeismicPropertyIntervals(std::vector<NRLib::Grid<float> >          & vp_interval,
                                                      std::vector<NRLib::Grid<float> >          & vs_interval,
                                                      std::vector<NRLib::Grid<float> >          & rho_interval,
                                                      const std::vector<IntervalSimbox>        & interval_simboxes,
                                                      std::vector<double>                      & relative_grid_resolution) const{
  (void) relative_grid_resolution;

  for(size_t i=0; i<n_intervals_; i++) {
    int    nx        = interval_simboxes[i].GetNx();
    int    ny        = interval_simboxes[i].GetNy();
    //double x_min     = interval_simboxes[i].GetXMin();
    //double y_min     = interval_simboxes[i].GetYMin();
    //double lx        = interval_simboxes[i].GetLX();
    //double ly        = interval_simboxes[i].GetLY();
    //double angle     = interval_simboxes[i].GetAngle();

    double  x;
    double  y;
    double  z_top;
    double  z_base;

    //Find maximum distance between the surfaces
    double max_distance = 0;

    for(int j=0; j<nx; j++) {
      for(int k=0; k<ny; k++) {
        interval_simboxes_[i].GetXYCoord(j,k,x,y);

        z_top  = interval_simboxes[i].GetTop(x,y);
        z_base = interval_simboxes[i].GetBot(x,y);

        if(z_top == RMISSING) {
          LogKit::LogFormatted(LogKit::Low,"ERROR: The top surface for interval \'"+interval_simboxes_[i].GetIntervalName()+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }
        else if(z_base == RMISSING) {
          LogKit::LogFormatted(LogKit::Low,"ERROR: The base surface for interval \'"+interval_simboxes_[i].GetIntervalName()+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }

        if(z_base-z_top > max_distance) {
          if(z_top != RMISSING && z_base != RMISSING)
            max_distance = z_base-z_top;
        }
      }
    }

    if(max_distance == 0) {
      LogKit::LogFormatted(LogKit::Low,"ERROR: Interval \'"+interval_simboxes_[i].GetIntervalName()+"\' has size zero. Check that its top surface is above the base surface.\n");
      exit(1);
    }

    //NRLib::Volume volume(x_min, y_min, lx, ly, interval_simboxes[i].GetTopSurface(), interval_simboxes[i].GetBotSurface(), angle);
    int nz_zone = 1;// static_cast<int>(std::ceil(max_distance/dz[i-1]));

    vp_interval[i]  = NRLib::Grid<float>(nx, ny, nz_zone, 0);
    vs_interval[i]  = NRLib::Grid<float>(nx, ny, nz_zone, 0);
    rho_interval[i] = NRLib::Grid<float>(nx, ny, nz_zone, 0);

    //For each interval, store the actual vertical resolution relative to the wanted.
    //relative_grid_resolution_.push_back();

  }
}

/*
void MultiIntervalGrid::BuildSeismicPropertyZones(InputFiles                               * input_files,
                                                  std::vector<NRLib::Grid<float>>          & alpha_zones,
                                                  std::vector<NRLib::Grid<float>>          & beta_zones,
                                                  std::vector<NRLib::Grid<float>>          & rho_zones,
                                                  const std::vector<Surface>               & surfaces,
                                                  const std::map<std::string, std::string> & interval_corr_dir_files,
                                                  const std::map<std::string, std::string> & interval_corr_dir_top_files,
                                                  const std::map<std::string, std::string> & interval_corr_dir_base_files,
                                                  const std::map<std::string, bool>        & interval_top_conforms,
                                                  const std::map<std::string, bool>        & interval_base_conforms,
                                                  const std::vector<std::string>           & interval_names,
                                                  const std::vector<float>                 & dz) const
{
  for(int i=1; i<n_intervals_+1; i++) {
    int    nx        = simboxes_[i].getnx();
    int    ny        = simboxes_[i].getny();
    double x_min     = simboxes_[i].GetXMin();
    double y_min     = simboxes_[i].GetYMin();
    double lx        = simboxes_[i].GetLX();
    double ly        = simboxes_[i].GetLY();
    double angle     = simboxes_[i].getAngle();

    Surface temp_top;
    Surface temp_base;
    double  x;
    double  y;
    double  z_top;
    double  z_base;

    Surface top  = surfaces[i-1];
    Surface base = surfaces[i];

    double top_missing  = top.GetMissingValue();
    double base_missing = base.GetMissingValue();

    //Find maximum distance between the surfaces
    double max_distance = 0;

    for(int j=0; j<nx; j++) {
      for(int k=0; k<ny; k++) {
        simboxes_[i].getXYCoord(j,k,x,y);

        z_top  = top.GetZ(x,y);
        z_base = base.GetZ(x,y);

        if(z_top == top_missing) {
          const std::string name = top.GetName();

          LogKit::LogFormatted(LogKit::Low,"ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }
        else if(z_base == base_missing) {
          const std::string name = base.GetName();

          LogKit::LogFormatted(LogKit::Low,"ERROR: Surface \'"+name+"\' does not cover the inversion grid, or it contains missing values.\n");
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

    //Make new top and base surfaces

    //Now:
    //1 Correlation direction file File
    //2 Correlation directions from top and base.
      //Top and base file
      //Top file and base conform
      //Top conform and base file

    if(interval_corr_dir_files.find(interval_names[i-1])->second != "") { //Ekstra flate
      temp_top  = Surface(interval_corr_dir_files.find(interval_names[i-1])->second);
      temp_base = temp_top;
      temp_base.Add(max_distance);
    }
    else {
      if(interval_top_conforms.find(interval_names[i-1])->second == true)
        temp_top = top;
      else
        temp_top = Surface(interval_corr_dir_top_files.find(interval_names[i-1])->second);

      if(interval_base_conforms.find(interval_names[i-1])->second == true)
        temp_base = base;
      else
        temp_base = Surface(interval_corr_dir_base_files.find(interval_names[i-1])->second);
    }

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

    NRLib::Volume volume(x_min, y_min, lx, ly, temp_top, temp_base, angle);

    int nz_zone = static_cast<int>(std::ceil(max_distance/dz[i-1]));

    alpha_zones[i-1] = StormContGrid(volume, nx, ny, nz_zone);
    beta_zones[i-1]  = StormContGrid(volume, nx, ny, nz_zone);
    rho_zones[i-1]   = StormContGrid(volume, nx, ny, nz_zone);

    //For each interval, store the actual vertical resolution relative to the wanted.
    //relative_grid_resolution_.push_back();

  }
}
*/
