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
                                     bool           & failed) {

  std::vector<std::string> interval_names_                        = model_settings->getIntervalNames();
  n_intervals_                                                    = static_cast<int>(interval_names_.size());
  int erosion_priority_top_surface                                = model_settings->getErosionPriorityTopSurface();
  const std::map<std::string,int> erosion_priority_base_surfaces  = model_settings->getErosionPriorityBaseSurfaces();

  Surface               * top_surface = NULL;
  Surface               * base_surface = NULL;
  //std::vector<Surface>    eroded_surfaces(n_intervals_+1);
  eroded_surfaces_.resize(n_intervals_+1);
  std::string             previous_interval_name("");
  std::string             top_surface_file_name_temp("");
  std::string             base_surface_file_name_temp("");
  std::vector<Surface>    surfaces;
  //std::vector<int>        erosion_priorities_;
  trend_cubes_.resize(n_intervals_);

  // if there are multiple intervals (there can potentially be 1 interval as well)
  if (interval_names_.size() > 0) {
    multiple_interval_setting_ = true;
    desired_grid_resolution_.resize(interval_names_.size());
    relative_grid_resolution_.resize(interval_names_.size());
    n_intervals_ = static_cast<int>(interval_names_.size());
    eroded_surfaces_.resize(n_intervals_ + 1);
    LogKit::WriteHeader("Setting up Multiple Interval Grid");
    surfaces.resize(n_intervals_+1); //Store surfaces.
    erosion_priorities_.resize(n_intervals_+1);
    interval_simboxes_.resize(n_intervals_);
    background_parameters_.resize(n_intervals_);
    for(int i = 0; i < n_intervals_; i++)
      background_parameters_[i].resize(3);
    background_vs_vp_ratios_.resize(n_intervals_);
  }
  // if there is only one interval
  else {
    multiple_interval_setting_ = false;
    desired_grid_resolution_.resize(1);
    relative_grid_resolution_.resize(1);
    eroded_surfaces_.resize(2);
    n_intervals_ = 1;
    LogKit::WriteHeader("Setting up Grid");
    surfaces.resize(1);
    interval_simboxes_.resize(1);
    background_parameters_.resize(1);
    background_parameters_[0].resize(3);
    background_vs_vp_ratios_.resize(1);
  }

  // 1. ERODE SURFACES AND SET SURFACES OF SIMBOXES -----------------------------------------

  try{
    // if multiple-intervals keyword is used in model settings
    if (multiple_interval_setting_){
      top_surface_file_name_temp = input_files->getTimeSurfFile(0);
      erosion_priorities_[0] = erosion_priority_top_surface;

      top_surface = MakeSurfaceFromFileName(top_surface_file_name_temp, *estimation_simbox);
      surfaces[0] = *top_surface;

      for (size_t i = 0; i < n_intervals_; i++) {

        std::string interval_name = model_settings->getIntervalName(i);
        base_surface_file_name_temp = input_files->getIntervalBaseTimeSurface(interval_name);
        erosion_priorities_[i+1] = erosion_priority_base_surfaces.find(interval_name)->second;

        base_surface = MakeSurfaceFromFileName(base_surface_file_name_temp, *estimation_simbox);
        surfaces[i+1] =  *base_surface;
      }



      if (!failed){

        for (size_t i = 0; i<n_intervals_; i++){
          desired_grid_resolution_[i] = FindResolution(&surfaces[i], &surfaces[i+1], estimation_simbox,
                                                     model_settings->getTimeNzInterval(interval_names_[i]));
        }

        ErodeAllSurfaces(eroded_surfaces_,
                         erosion_priorities_,
                         surfaces,
                         *estimation_simbox);
      }
      else{
        err_text += "Erosion of surfaces failed because the interval surfaces could not be set up correctly.\n";
      }
    }
    // if multiple-intervals is NOT used in model settings
    else{
      top_surface_file_name_temp = input_files->getTimeSurfFile(0);
      top_surface = MakeSurfaceFromFileName(top_surface_file_name_temp, *estimation_simbox);
      eroded_surfaces_[0] = *top_surface;

      base_surface_file_name_temp = input_files->getTimeSurfFile(1);
      base_surface = MakeSurfaceFromFileName(base_surface_file_name_temp, *estimation_simbox);
      eroded_surfaces_[1] = *base_surface;

      desired_grid_resolution_[0] = FindResolution(top_surface, base_surface, estimation_simbox,
                                                    model_settings->getTimeNz());
    }
  }
  catch(NRLib::Exception & e){
    failed = true;
    err_text += e.what();
  }

  delete top_surface;
  delete base_surface;

  // 2 SET UP INTERVAL_SIMBOXES ------------------------------------------------------------

  //Set up a vector of simboxes, one per interval.

  //interval_simboxes_.resize(interval_names_.size()); //H Removed

  //H Testing while SetupIntervalSimbox is uncomplete.
  //interval_names_.push_back("test");
  //multiple_interval_setting_ = true;

  if (!failed){
    try{
      // if multiple-intervals keyword is used in model settings
      if (multiple_interval_setting_){

        SetupIntervalSimboxes(model_settings,
                              estimation_simbox,
                              interval_names_,
                              eroded_surfaces_,
                              interval_simboxes_,
                              input_files->getCorrDirIntervalFiles(),
                              input_files->getCorrDirIntervalTopSurfaceFiles(),
                              input_files->getCorrDirIntervalBaseSurfaceFiles(),
                              model_settings->getCorrDirIntervalTopConform(),
                              model_settings->getCorrDirIntervalBaseConform(),
                              desired_grid_resolution_,
                              relative_grid_resolution_,
                              err_text,
                              failed);
      }
      // if multiple-intervals is NOT used in model settings
      else{

        SetupIntervalSimbox(model_settings,
                            estimation_simbox,
                            interval_simboxes_[0],
                            input_files->getCorrDirFile(),
                            input_files->getCorrDirTopFile(),
                            input_files->getCorrDirBaseFile(),
                            model_settings->getCorrDirTopConform(),
                            model_settings->getCorrDirBaseConform(),
                            desired_grid_resolution_[0],
                            relative_grid_resolution_[0],
                            err_text,
                            failed);
      }
    }
    catch(NRLib::Exception & e){
      failed = true;
      err_text += e.what();
    }
  }


  // 3. SET UP BACKGROUND MODEL ----------------------------------------------------------

  //if (model_settings->getIntervalNames().size() > 0) {
  //  parameters_.resize(model_settings->getIntervalNames().size());
  //  background_vs_vp_ratios_.resize(model_settings->getIntervalNames().size());
  //}

  //std::vector<NRLib::Grid<double> > vp_intervals(n_intervals_);
  //std::vector<NRLib::Grid<double> > vs_intervals(n_intervals_);
  //std::vector<NRLib::Grid<double> > rho_intervals(n_intervals_);

  //if (!failed){
  //  try{

  //    //dz to vector
  //    //std::vector<float> dz;
  //    //for (size_t i = 0; i < n_intervals_; i++)
  //    //  dz.push_back( interval_simboxes_[i].GetDz()*interval_simboxes_[i].GetAvgRelThick() * 4 );
  //    //float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

  //    if (n_intervals_ > 0){ // It is possible to give only one interval in multiple-intervals
  //      BuildSeismicPropertyIntervals(vp_intervals,
  //                                    vs_intervals,
  //                                    rho_intervals,
  //                                    interval_simboxes_,
  //                                    relative_grid_resolution_);
  //    }else{
  //      // what ?
  //    }
  //  }
  //  catch(NRLib::Exception & e){
  //  failed = true;
  //  err_text += e.what();
  //  }
  //}

  // Add inn surface files.
  surface_files_.push_back(input_files->getTimeSurfFile(0));

  const std::map<std::string, std::string> & interval_base_time_surfaces = input_files->getIntervalBaseTimeSurfaces();
  for (size_t i = 0; i > n_intervals_; i++) {
    surface_files_.push_back(interval_base_time_surfaces.find(interval_names_[i])->second);
  }

}

MultiIntervalGrid::~MultiIntervalGrid(){

}

// ---------------------------------------------------------------------------------------------------------------
void  MultiIntervalGrid::SetupIntervalSimbox(ModelSettings                               * model_settings,
                                             const Simbox                                * estimation_simbox,
                                             Simbox                                      & interval_simbox,
                                             const std::string                           & corr_dir_single_surf,
                                             const std::string                           & corr_dir_top_surf,
                                             const std::string                           & corr_dir_base_surf,
                                             bool                                          corr_dir_top_conform,
                                             bool                                          corr_dir_base_conform,
                                             double                                      & desired_grid_resolution,
                                             double                                      & relative_grid_resolution,
                                             std::string                                 & err_text,
                                             bool                                        & failed) const{


  // H TESTING WITH ONE CORRELATION DIRECTION
  Surface * corr_surf     = MakeSurfaceFromFileName(corr_dir_single_surf,  estimation_simbox);
  int n_layers = model_settings->getTimeNz();
  Surface top_surface = eroded_surfaces_[0];
  Surface base_surface = eroded_surfaces_[1];
  interval_simbox = Simbox(estimation_simbox, "test_interval", n_layers, top_surface, base_surface, corr_surf, err_text, failed);

    const SegyGeometry * area_params = model_settings->getAreaParameters();
    failed = interval_simbox.setArea(area_params, err_text);

  interval_simbox.SetTopBotErodedSurfaces(top_surface, base_surface);

  EstimateXYPaddingSizes(&interval_simbox, model_settings);





  // Calculate Z padding ----------------------------------------------------------------
  int status = interval_simbox.calculateDz(model_settings->getLzLimit(),err_text);
  EstimateZPaddingSize(&interval_simbox, model_settings);
}

// ---------------------------------------------------------------------------------------------------------------
void   MultiIntervalGrid::SetupIntervalSimboxes(ModelSettings                             * model_settings,
                                                const Simbox                              * estimation_simbox,
                                                const std::vector<std::string>            & interval_names,
                                                const std::vector<Surface>                & eroded_surfaces,
                                                std::vector<Simbox>                       & interval_simboxes,
                                                const std::map<std::string, std::string>  & corr_dir_single_surfaces,
                                                const std::map<std::string, std::string>  & corr_dir_top_surfaces,
                                                const std::map<std::string, std::string>  & corr_dir_base_surfaces,
                                                const std::map<std::string, bool>         & corr_dir_top_conform,
                                                const std::map<std::string, bool>         & corr_dir_base_conform,
                                                std::vector<double>                       & desired_grid_resolution,
                                                std::vector<double>                       & relative_grid_resolution,
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

    float min_samp_dens = model_settings->getMinSamplingDensity();
    double tmp = interval_simboxes[i].getMinRelThick();

    //if (interval_simboxes[i].getdz()*interval_simboxes[i].getMinRelThick() < min_samp_dens){ ///H Commented out. Move this?
    //  failed   = true;
    //  err_text += "We normally discourage denser sampling than "+NRLib::ToString(min_samp_dens);
    //  err_text += "ms in the time grid. If you really need\nthis, please use ";
    //  err_text += "<project-settings><advanced-settings><minimum-sampling-density>\n";
    //}


    // Make extended interval_simbox for the inversion interval ---------------------------

    if (interval_simboxes[i].status() == Simbox::EMPTY){
      LogIntervalInformation(interval_simboxes[i], interval_names[i], "Time output interval:","Two-way-time");

      // Case 1: Single correlation surface
      if (it_single != corr_dir_single_surfaces.end() && it_top == corr_dir_top_surfaces.end() && it_base == corr_dir_base_surfaces.end()){
        Surface * corr_surf     = MakeSurfaceFromFileName(it_single->second,  estimation_simbox);
        interval_simboxes[i]    = Simbox(estimation_simbox, interval_names[i], n_layers, top_surface, base_surface, corr_surf, err_text, failed);
        interval_simboxes[i].SetTopBotErodedSurfaces(top_surface, base_surface);
      }
      // Case 2: Top and base correlation surfaces
      else if (it_single == corr_dir_single_surfaces.end() && it_top != corr_dir_top_surfaces.end() && it_base != corr_dir_base_surfaces.end()){
        Surface * corr_surf_top = MakeSurfaceFromFileName(it_top->second,  estimation_simbox);
        Surface * corr_surf_base = MakeSurfaceFromFileName(it_base->second,  estimation_simbox);
        interval_simboxes[i] = Simbox(estimation_simbox, interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                                corr_surf_top, corr_surf_base);
        interval_simboxes[i].SetTopBotErodedSurfaces(top_surface, base_surface);
      }
      // Case 3: Top conform and base conform if (i) both are set conform or (ii) if no other corr surfaces have been defined
      else if ((it_top_conform->second == true && it_base_conform->second == true) ||
               (it_single == corr_dir_single_surfaces.end() && it_top == corr_dir_top_surfaces.end() && it_base == corr_dir_base_surfaces.end())){
        interval_simboxes[i] = Simbox(estimation_simbox, interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                                NULL, NULL);
        interval_simboxes[i].SetTopBotErodedSurfaces(top_surface, base_surface);
      }
      // Case 4: Top correlation surface and base conform
      else if (it_top != corr_dir_top_surfaces.end() && it_base_conform->second == true){
        Surface * corr_surf_top = MakeSurfaceFromFileName(it_top->second,  estimation_simbox);
        interval_simboxes[i] = Simbox(estimation_simbox, interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                                corr_surf_top, NULL);
        interval_simboxes[i].SetTopBotErodedSurfaces(top_surface, base_surface);
      }
      // Case 5: Top conform and base correlation surface
      else if (it_top_conform == corr_dir_top_conform.end() && it_base_conform != corr_dir_base_conform.end()){
        Surface * corr_surf_base = MakeSurfaceFromFileName(it_base->second,  estimation_simbox);
        interval_simboxes[i] = Simbox(estimation_simbox, interval_names[i], n_layers, top_surface, base_surface, err_text, failed,
                                                NULL, corr_surf_base);
        interval_simboxes[i].SetTopBotErodedSurfaces(top_surface, base_surface);
      }
      // else something is wrong
      else{
        err_text += "\nCorrelation directions are not set correctly for interval " + interval_name[i];
        err_text += ".\n";
        failed = true;
      }
    }
    // Calculate Z padding ----------------------------------------------------------------

    if (!failed){
      int status = interval_simboxes[i].calculateDz(model_settings->getLzLimit(),err_text);
      EstimateZPaddingSize(&interval_simboxes[i], model_settings);
      relative_grid_resolution[i] = interval_simboxes[i].getdz() / desired_grid_resolution[i];


      if (status == Simbox::BOXOK)
        LogIntervalInformation(&interval_simboxes[i], "Time inversion interval (extended relative to output interval due to correlation):","Two-way-time");
      else
      {
        err_text += "Could not make the time simulation grid.\n";
        failed = true;
      }

    // Calculate XY padding ---------------------------------------------------------------

    if (!failed) {

      EstimateXYPaddingSizes(&interval_simboxes[i], model_settings);

            unsigned long long int grid_size = static_cast<unsigned long long int>(model_settings->getNXpad())*model_settings->getNYpad()*model_settings->getNZpad();

            if (grid_size > std::numeric_limits<unsigned int>::max()) {
              float fsize = 4.0f*static_cast<float>(grid_size)/static_cast<float>(1024*1024*1024);
              float fmax  = 4.0f*static_cast<float>(std::numeric_limits<unsigned int>::max()/static_cast<float>(1024*1024*1024));
              err_text += "Grids as large as "+NRLib::ToString(fsize,1)+"GB cannot be handled. The largest accepted grid size\n";
              err_text += "is "+NRLib::ToString(fmax)+"GB. Please reduce the number of layers or the lateral resolution.\n";
              failed = true;
            }

            LogKit::LogFormatted(LogKit::Low,"\nTime simulation grids for interval \'"+interval_names[i]+"\':\n");
            LogKit::LogFormatted(LogKit::Low,"  Output grid        %4i * %4i * %4i   : %10llu\n",
                                 interval_simboxes[i].getnx(),interval_simboxes[i].getny(),interval_simboxes[i].getnz(),
                                 static_cast<unsigned long long int>(interval_simboxes[i].getnx())*interval_simboxes[i].getny()*interval_simboxes[i].getnz());
            LogKit::LogFormatted(LogKit::Low,"  FFT grid            %4i * %4i * %4i   :%11llu\n",
                                 model_settings->getNXpad(),model_settings->getNYpad(),model_settings->getNZpad(),
                                 static_cast<unsigned long long int>(model_settings->getNXpad())*model_settings->getNYpad()*model_settings->getNZpad());
          }
    }

    // Check consistency ------------------------------------------------------------------
    if (interval_simboxes[i].getdz() >= 10.0 && model_settings->getFaciesProbFromRockPhysics() == true) {
      err_text += "dz for interval \'" + interval_names[i] + "\' is too large to generate synthetic well data when estimating facies probabilities using rock physics models. Need dz < 10.";
      failed = true;
    }
  }
}

// --------------------------------------------------------------------------------
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

// --------------------------------------------------------------------------------
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

// --------------------------------------------------------------------------------
void  MultiIntervalGrid::ErodeAllSurfaces(std::vector<Surface>       & eroded_surfaces,
                                          const std::vector<int>     & erosion_priorities,
                                          const std::vector<Surface> & surfaces,
                                          const Simbox               & simbox) const{
  int    n_surf     = static_cast<int>(eroded_surfaces.size());

  for (int i=0; i<n_surf; i++) {
    int l=0;
    while(i+1 != erosion_priorities[l])
      l++;

    Surface temp_surface = Surface(surfaces[l]);

    //Find closest eroded surface downward
    for (int k=l+1; k<n_surf; k++) {
      if (eroded_surfaces[k].GetN() > 0) {
        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, false);
        break;
      }
    }
    //Find closest eroded surface upward
    for (int k=l-1; k>=0; k--) {
      if (eroded_surfaces[k].GetN() > 0) {
        ErodeSurface(temp_surface, eroded_surfaces[k], simbox, true);
        break;
      }
    }
    eroded_surfaces[l] = temp_surface;
  }
}

// --------------------------------------------------------------------------------
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
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      simbox.getXYCoord(i,j,x,y);

      z_priority = priority_surface.GetZ(x,y);
      z          = surface.GetZ(x,y);

      if (compare_upward) {
        if (z < z_priority && z != missing)
          eroded_surface(i,j) = z_priority;
        else
          eroded_surface(i,j) = z;
      }

      else {
        if (z > z_priority && z_priority != missing)
          eroded_surface(i,j) = z_priority;
        else
          eroded_surface(i,j) = z;
      }
    }
  }


  surface = Surface(x0, y0, lx, ly, eroded_surface);
}

void MultiIntervalGrid::BuildSeismicPropertyIntervals(std::vector<NRLib::Grid<double> >          & vp_interval,
                                                      std::vector<NRLib::Grid<double> >          & vs_interval,
                                                      std::vector<NRLib::Grid<double> >          & rho_interval,
                                                      const std::vector<Simbox>                 & interval_simboxes,
                                                      std::vector<double>                       & relative_grid_resolution) const{
  (void) relative_grid_resolution;

  for (size_t i=0; i<n_intervals_; i++) {
    int    nx        = interval_simboxes[i].getnx();
    int    ny        = interval_simboxes[i].getny();
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

    for (int j=0; j<nx; j++) {
      for (int k=0; k<ny; k++) {
        interval_simboxes[i].getXYCoord(j,k,x,y);

        z_top  = interval_simboxes[i].getTop(x,y);
        z_base = interval_simboxes[i].getBot(x,y);

        if (z_top == RMISSING) {
          LogKit::LogFormatted(LogKit::Low,"ERROR: The top surface for interval \'"+interval_simboxes_[i].GetIntervalName()+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }
        else if (z_base == RMISSING) {
          LogKit::LogFormatted(LogKit::Low,"ERROR: The base surface for interval \'"+interval_simboxes_[i].GetIntervalName()+"\' does not cover the inversion grid, or it contains missing values.\n");
          exit(1);
        }

        if (z_base-z_top > max_distance) {
          if (z_top != RMISSING && z_base != RMISSING)
            max_distance = z_base-z_top;
        }
      }
    }

    if (max_distance == 0) {
      LogKit::LogFormatted(LogKit::Low,"ERROR: Interval \'"+interval_simboxes_[i].GetIntervalName()+"\' has size zero. Check that its top surface is above the base surface.\n");
      exit(1);
    }

    //NRLib::Volume volume(x_min, y_min, lx, ly, interval_simboxes[i].GetTopSurface(), interval_simboxes[i].GetBotSurface(), angle);
    int nz_zone = 1;// static_cast<int>(std::ceil(max_distance/dz[i-1]));

    vp_interval[i]  = NRLib::Grid<double>(nx, ny, nz_zone, 0);
    vs_interval[i]  = NRLib::Grid<double>(nx, ny, nz_zone, 0);
    rho_interval[i] = NRLib::Grid<double>(nx, ny, nz_zone, 0);

    //For each interval, store the actual vertical resolution relative to the wanted.
    //relative_grid_resolution_.push_back();

  }
}

// --------------------------------------------------------------------------------
void MultiIntervalGrid::EstimateZPaddingSize(Simbox          * simbox,
                                             ModelSettings   * model_settings) const{
  int    nz             = simbox->getnz();
  double min_lz         = simbox->getlz()*simbox->getMinRelThick();
  double z_pad_fac      = model_settings->getZPadFac();
  double z_pad          = z_pad_fac*min_lz;

  if (model_settings->getEstimateZPadding())
  {
    double w_length    = static_cast<double>(model_settings->getDefaultWaveletLength());
    double p_fac      = 1.0;
    z_pad             = w_length/p_fac;                             // Use half a wavelet as padding
    z_pad_fac         = std::min(1.0, z_pad/min_lz);                  // More than 100% padding is not sensible
  }
  int nz_pad          = SetPaddingSize(nz, z_pad_fac);
  z_pad_fac           = static_cast<double>(nz_pad - nz)/static_cast<double>(nz);

  model_settings->setNZpad(nz_pad);
  model_settings->setZPadFac(z_pad_fac);
}

// --------------------------------------------------------------------------------
int MultiIntervalGrid::SetPaddingSize(int     nx,
                                      double  px) const{

  int leastint    = static_cast<int>(ceil(nx*(1.0f+px)));
  int closestprod = FindClosestFactorableNumber(leastint);
  return(closestprod);
}

// --------------------------------------------------------------------------------
int MultiIntervalGrid::FindClosestFactorableNumber(int leastint) const{
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

// --------------------------------------------------------------------------------
void  MultiIntervalGrid::LogIntervalInformation(const Simbox      & simbox,
                                                const std::string & interval_name,
                                                const std::string & header_text1,
                                                const std::string & header_text2) const{
  LogKit::LogFormatted(LogKit::Low,"\n"+header_text1+"\n");
  double zmin, zmax;
  simbox.getMinMaxZ(zmin,zmax);
  LogKit::LogFormatted(LogKit::Low," Interval name: "+interval_name +"\n");
  LogKit::LogFormatted(LogKit::Low," %13s          avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       header_text2.c_str(),
                       zmin+simbox.getlz()*simbox.getAvgRelThick()*0.5,
                       zmin,zmax);
  LogKit::LogFormatted(LogKit::Low,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       simbox.getlz()*simbox.getAvgRelThick(),
                       simbox.getlz()*simbox.getMinRelThick(),
                       simbox.getlz());
  LogKit::LogFormatted(LogKit::Low,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n",
                       simbox.getdz()*simbox.getAvgRelThick(),
                       simbox.getdz(),
                       simbox.getdz()*simbox.getMinRelThick());
}

void MultiIntervalGrid::LogIntervalInformation(const Simbox      * simbox,
                                               const std::string & header_text1,
                                               const std::string & header_text2) const{
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

// --------------------------------------------------------------------------------
double  MultiIntervalGrid::FindResolution(const Surface * top_surface,
                                          const Surface * base_surface,
                                          const Simbox  * estimation_simbox,
                                          int             n_layers) const{
  size_t nx = top_surface->GetNI();
  size_t ny = base_surface->GetNJ();

  double max_resolution = 0;

  for (size_t i = 0; i<nx; i++){
    for (size_t j = 0; j<ny; j++){
      double x,y;
      estimation_simbox->getXYCoord(i,j,x,y);
      double z_top  = top_surface->GetZ(x,y);
      double z_base = base_surface->GetZ(x,y);
      //double resolution = (z_top - z_base) / n_layers ;
      double resolution = (z_base - z_top) / n_layers; //H Changed since z_top < z_base
      if (resolution > max_resolution)
        max_resolution = resolution;
    }
  }

  return max_resolution;
}

// --------------------------------------------------------------------------------
void  MultiIntervalGrid::EstimateXYPaddingSizes(Simbox          * interval_simbox,
                                                ModelSettings   * model_settings) const{
  double dx      = interval_simbox->getdx();
  double dy      = interval_simbox->getdy();
  double lx      = interval_simbox->getlx();
  double ly      = interval_simbox->getly();
  int    nx      = interval_simbox->getnx();
  int    ny      = interval_simbox->getny();
  int    nz      = interval_simbox->getnz();

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

  int nxPad = SetPaddingSize(nx, xPadFac);
  int nyPad = SetPaddingSize(ny, yPadFac);
  int nzPad = model_settings->getNZpad();

  double true_xPadFac = static_cast<double>(nxPad - nx)/static_cast<double>(nx);
  double true_yPadFac = static_cast<double>(nyPad - ny)/static_cast<double>(ny);
  double true_zPadFac = model_settings->getZPadFac();
  double true_xPad    = true_xPadFac*lx;
  double true_yPad    = true_yPadFac*ly;
  double true_zPad    = true_zPadFac*(interval_simbox->getlz()*interval_simbox->getMinRelThick());

  model_settings->setNXpad(nxPad);
  model_settings->setNYpad(nyPad);
  model_settings->setXPadFac(true_xPadFac);
  model_settings->setYPadFac(true_yPadFac);

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
                       true_xPad, true_xPadFac, nx, nxPad);
  LogKit::LogFormatted(logLevel,"  yPad, yPadFac, ny, nyPad                 : %6.fm, %5.3f, %5d, %4d\n",
                       true_yPad, true_yPadFac, ny, nyPad);
  LogKit::LogFormatted(logLevel,"\nPadding sizes"+text2+":\n");
  LogKit::LogFormatted(logLevel,"  zPad, zPadFac, nz, nzPad                 : %5.fms, %5.3f, %5d, %4d\n",
                       true_zPad, true_zPadFac, nz, nzPad);
}
