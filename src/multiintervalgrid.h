/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MULTIINTERVALGRID_H
#define MULTIINTERVALGRID_H

#include "src/definitions.h"
#include "src/intervalsimbox.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"

class MultiIntervalGrid
{
public:
  MultiIntervalGrid();

  MultiIntervalGrid(ModelSettings  * model_settings,
                    InputFiles     * input_files,
                    const Simbox   * estimation_simbox,
                    std::string    & err_text,
                    bool           & failed);

  ~MultiIntervalGrid();

  //GET FUNCTIONS
  int                                                      GetNIntervals()       const           { return n_intervals_          ;}
  const std::vector<Simbox>                              & GetSimboxes() const                   { return simboxes_             ;}
  const std::vector<IntervalSimbox>                      & GetIntervalSimboxes() const           { return interval_simboxes_    ;}
  const Simbox                                           * GetSimbox(int i) const                { return &simboxes_[i]          ;}
  const IntervalSimbox                                   * GetIntervalSimbox(int i) const        { return &interval_simboxes_[i] ;}
  const std::vector<std::vector<NRLib::Grid<double> > >  & GetParametersAllIntervals() const     { return parameters_               ;}
  const std::vector<NRLib::Grid<double> >                & GetParametersForInterval(int i) const;

  //SET FUNCTIONS


private:

  void  SetUpIntervalSimboxes(const ModelSettings                       * model_settings,
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
                              bool                                      & failed) const;

  Surface * MakeSurfaceFromFileName(const std::string        & file_name,
                                    const Simbox             & estimation_simbox) const;

  void ErodeSurface(Surface       &  surface,
                    const Surface &  priority_surface,
                    const Simbox  &  simbox,
                    const bool    &  compare_upward) const;

  void ErodeAllSurfaces(std::vector<Surface>            & eroded_surfaces,
                        const std::vector<int>          & erosion_priorities,
                        const std::vector<Surface>      & surface,
                        const Simbox                    & simbox) const;

  void  FindSmallestSurfaceGeometry(const double   x0,
                                    const double   y0,
                                    const double   lx,
                                    const double   ly,
                                    const double   rot,
                                    double       & x_min,
                                    double       & y_min,
                                    double       & x_max,
                                    double       & y_max) const;

  void  BuildSeismicPropertyIntervals(std::vector<NRLib::Grid<double> >          & vp_interval,
                                      std::vector<NRLib::Grid<double> >          & vs_interval,
                                      std::vector<NRLib::Grid<double> >          & rho_interval,
                                      const std::vector<IntervalSimbox>         & interval_simboxes,
                                      std::vector<double>                       & relative_grid_resolution) const;

  //void  BuildSeismicPropertyZones(InputFiles                               * input_files,
  //                                std::vector<NRLib::Grid<float>>          & alpha_zones,
  //                                std::vector<NRLib::Grid<float>>          & beta_zones,
  //                                std::vector<NRLib::Grid<float>>          & rho_zones,
  //                                const std::vector<Surface>               & surfaces,
  //                                const std::map<std::string, std::string> & interval_corr_dir_files,
  //                                const std::map<std::string, std::string> & interval_corr_dir_top_files,
  //                                const std::map<std::string, std::string> & interval_corr_dir_base_files,
  //                                const std::map<std::string, bool>        & top_conforms,
  //                                const std::map<std::string, bool>        & base_conforms,
  //                                const std::vector<std::string>           & interval_names,
  //                                const std::vector<float>                 & dz) const;

  // CLASS VARIABLES

  size_t n_intervals_;

  std::vector<Simbox>                                  simboxes_;                 // original inversion interval without correlation directions
  std::vector<IntervalSimbox>                          interval_simboxes_;        // extended simbox with correlation direction, must have same size as the parameters vector
  std::vector<std::vector<NRLib::Grid<double> > >      parameters_;               // must have same size as the simbox vector

  std::vector<double>                                  relative_grid_resolution_; //Actual grid resolution relative to the wanted grid resolution.

};

#endif
