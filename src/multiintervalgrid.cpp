/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/modelsettings.h"
#include "src/inputfiles.h"

#include "src/multiintervalgrid.h"
#include "src/definitions.h"
#include "src/simbox.h"

#include "nrlib/grid/grid.hpp"

MultiIntervalGrid::MultiIntervalGrid()
{
}

MultiIntervalGrid::MultiIntervalGrid(ModelSettings  * model_settings,
                                     InputFiles     * input_files)
{
  n_intervals_ = model_settings->getIntervalNames().size();
  std::vector<std::string> interval_names = model_settings->getIntervalNames();

  Surface * top_surface = NULL;
  Surface * base_surface = NULL;

  std::string last_interval_name("");
  std::string top_surface_file_name("");
  std::string base_surface_file_name("");

  std::vector<Surface> surfaces(n_intervals_+1); //Store surfaces.

  //Set up a vector of simboxes, one per inteval.
  for(size_t i = 0; i < n_intervals_; i++) {

    Simbox interval_simbox;
    std::string interval_name = model_settings->getIntervalName(i);

    if(i == 0) { //First interval
      top_surface_file_name = input_files->getTimeSurfFile(0);
      base_surface_file_name = input_files->getIntervalBaseTimeSurface(interval_name);
    }
    else {
      last_interval_name = model_settings->getIntervalName(i-1); //Base surface of last interval is this interval's top surface
      top_surface_file_name = input_files->getIntervalBaseTimeSurface(last_interval_name);
      base_surface_file_name = input_files->getIntervalBaseTimeSurface(interval_name);
    }

    if (!NRLib::IsNumber(top_surface_file_name)) { //Time file
      top_surface = new Surface(top_surface_file_name);

    }
    else { //Hvis Time er gitt som en value og ikke fil (fra Commondata.cpp:)
      const SegyGeometry * geometry = model_settings->getAreaParameters();

      std::string err_txt("");
      interval_simbox.setArea(geometry, err_txt);
      double x_min, x_max;
      double y_min, y_max;
      FindSmallestSurfaceGeometry(interval_simbox.getx0(), interval_simbox.gety0(),
                                  interval_simbox.getlx(), interval_simbox.getly(),
                                  interval_simbox.getAngle(), x_min, y_min, x_max, y_max);
      top_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(top_surface_file_name.c_str()));

      ///H Evt. Kopiere forrige flate og flytte den ned?
    }
    if (!NRLib::IsNumber(base_surface_file_name)) {
      base_surface = new Surface(base_surface_file_name);
    }
    else {
      const SegyGeometry * geometry = model_settings->getAreaParameters();

      std::string err_txt("");
      interval_simbox.setArea(geometry, err_txt);
      double x_min, x_max;
      double y_min, y_max;
      FindSmallestSurfaceGeometry(interval_simbox.getx0(), interval_simbox.gety0(),
                                  interval_simbox.getlx(), interval_simbox.getly(),
                                  interval_simbox.getAngle(), x_min,y_min,x_max,y_max);
      base_surface = new Surface(x_min-100, y_min-100, x_max-x_min+200, y_max-y_min+200, 2, 2, atof(base_surface_file_name.c_str()));
    }

    if(i == 0)
      surfaces.push_back(*top_surface);
    else
      surfaces.push_back(*base_surface);

    interval_simbox.SetSurfaces(*top_surface, *base_surface, model_settings->getRunFromPanel());

    simboxes_.push_back(interval_simbox);
  }

  //std::vector<int> correlation_structure = model_settings->getCorrelationStructure();
  std::vector<int> erosion_priority      = model_settings->getErosionPriority();
  int erosion_priority_top_surface       = model_settings->getErosionPriorityTopSurface();
  const std::map<std::string,int> & erosion_priority_base_surfaces = model_settings->getErosionPriorityBaseSurfaces();

  //int    nWells    = modelSettings->getNumberOfWells();
  //int    nZones    = static_cast<int>(correlation_structure.size()) - 1;
  //dz to vector
  std::vector<float> dz;
  for(size_t i = 0; i < n_intervals_; i++)
    dz.push_back( simboxes_[i].getdz()*simboxes_[i].getAvgRelThick() * 4 );
  //float  dz        = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick()) * 4; //NBNB Marit: Multiply by 4 to save memory

  std::vector<NRLib::Grid<float>> alpha(n_intervals_);
  std::vector<NRLib::Grid<float>> beta(n_intervals_);
  std::vector<NRLib::Grid<float>> rho(n_intervals_);

  BuildSeismicPropertyZones(input_files,
                            alpha,
                            beta,
                            rho,
                            surfaces,
                            input_files->getCorrDirIntervalFiles(),
                            input_files->getCorrDirIntervalTopSurfaceFiles(),
                            input_files->getCorrDirIntervalBaseSurfaceFiles(),
                            model_settings->getCorrDirIntervalTopConforms(),
                            model_settings->getCorrDirIntervalBaseConforms(),
                            interval_names,
                            dz);

  //std::vector<Surface *> eroded_surfaces(nZones+1);
  //for(int i=0; i<nZones+1; i++)
  //  eroded_surfaces[i] = NULL;

  //ErodeAllSurfaces(eroded_surfaces,
  //                 erosion_priority,
  //                 surface,
  //                 simbox);

  //Vektor av parametere, hver en vektor med lengde lik antall simboxer. Hvert element et 3DGrid = NRLIb::Grid

  //std::vector<NRLib::Grid<float>> alpha(simboxes_.size());
  //std::vector<NRLib::Grid<float>> beta(simboxes_.size());
  //std::vector<NRLib::Grid<float>> rho(simboxes_.size());


}

void MultiIntervalGrid::FindSmallestSurfaceGeometry(const double   x0,
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



MultiIntervalGrid::~MultiIntervalGrid()
{
}
