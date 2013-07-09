/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MULTIINTERVALGRID_H
#define MULTIINTERVALGRID_H

#include <math.h>
#include <string>

#include "src/definitions.h"

class ModelSettings;
class InputFiles;

class MultiIntervalGrid
{
public:
  MultiIntervalGrid();

  MultiIntervalGrid(ModelSettings  * model_settings,
                    InputFiles     * input_files);

  ~MultiIntervalGrid();

  //GET FUNCTIONS
  const std::vector<Simbox> & GetSimboxes() const { return simboxes_ ; }


private:

  void  FindSmallestSurfaceGeometry(const double   x0,
                                    const double   y0,
                                    const double   lx,
                                    const double   ly,
                                    const double   rot,
                                    double       & x_min,
                                    double       & y_min,
                                    double       & x_max,
                                    double       & y_max);

  void  BuildSeismicPropertyZones(InputFiles                               * input_files,
                                  std::vector<NRLib::Grid<float>>          & alpha_zones,
                                  std::vector<NRLib::Grid<float>>          & beta_zones,
                                  std::vector<NRLib::Grid<float>>          & rho_zones,
                                  const std::vector<Surface>               & surfaces,
                                  const std::map<std::string, std::string> & interval_corr_dir_files,
                                  const std::map<std::string, std::string> & interval_corr_dir_top_files,
                                  const std::map<std::string, std::string> & interval_corr_dir_base_files,
                                  const std::map<std::string, bool>        & top_conforms,
                                  const std::map<std::string, bool>        & base_conforms,
                                  const std::vector<std::string>           & interval_names,
                                  const std::vector<float>                 & dz) const;

  int n_intervals_;

  std::vector<Simbox> simboxes_;

  std::vector<double> relative_grid_resolution_; //Actual grid resolution relative to the wanted grid resolution.

};

#endif
