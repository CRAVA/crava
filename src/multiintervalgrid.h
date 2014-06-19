/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MULTIINTERVALGRID_H
#define MULTIINTERVALGRID_H

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "nrlib/flens/nrlib_flens.hpp"

class MultiIntervalGrid
{
public:
  MultiIntervalGrid();

  MultiIntervalGrid(ModelSettings  * model_settings,
                    InputFiles     * input_files,
                    const Simbox   * estimation_simbox,
                    std::string    & err_text,
                    bool           & failed);

  MultiIntervalGrid(const MultiIntervalGrid * multi_interval_grid);

  ~MultiIntervalGrid();

  // MIXED FUNCTIONS
  int                                                      WhichSimbox(double x, double y, double z) const;
  static int                                               FindPaddingSize(int nx, double px);
  static int                                               FindClosestFactorableNumber(int leastint);

  // GET FUNCTIONS
  int                                                      GetNIntervals()                        const { return n_intervals_                 ;}
  const std::vector<Simbox *>                            & GetIntervalSimboxes()                  const { return interval_simboxes_           ;}
  std::vector<Simbox *>                                  & GetIntervalSimboxes()                        { return interval_simboxes_           ;}
  const Simbox                                           * GetIntervalSimbox(int i)               const { return interval_simboxes_[i]        ;}
  Simbox                                                 * GetIntervalSimbox(int i)                     { return interval_simboxes_[i]        ;}
  const std::string                                      & GetIntervalName(int i)                 const { return interval_names_[i]           ;}
  const std::vector<std::string>                         & GetIntervalNames()                     const { return interval_names_              ;}
  const std::vector<int>                                 & GetErosionPriorities()                 const { return erosion_priorities_          ;}
  const std::vector<std::string>                         & GetSurfaceFiles()                      const { return surface_files_               ;}
  const std::vector<double>                              & GetRelativeGridResolution()            const { return relative_grid_resolution_    ;}
  const std::vector<double>                              & GetDesiredGridResolution()             const { return desired_grid_resolution_     ;}
  double                                                   GetDzMin()                             const { return dz_min_                      ;}
  const std::vector<double>                              & GetDzRel()                             const { return dz_rel_                      ;}

  // SET FUNCTIONS
  void SetDzRel(std::vector<double> & dz_rel)                                                           { dz_rel_ = dz_rel                    ;}

private:

  void  SetupIntervalSimboxes(ModelSettings                             * model_settings,
                              const Simbox                              * estimation_simbox,
                              const std::vector<std::string>            & interval_names,
                              const std::vector<Surface>                & eroded_surfaces,
                              std::vector<Simbox *>                       & interval_simboxes,
                              const std::map<std::string, std::string>  & corr_dir_single_surfaces,
                              const std::map<std::string, std::string>  & corr_dir_top_surfaces,
                              const std::map<std::string, std::string>  & corr_dir_base_surfaces,
                              const std::map<std::string, bool>         & corr_dir_top_conform,
                              const std::map<std::string, bool>         & corr_dir_base_conform,
                              std::vector<double>                       & desired_grid_resolution,
                              std::vector<double>                       & relative_grid_resolution,
                              double                                    & dz_min,
                              std::vector<double>                       & dz_rel,
                              std::string                               & err_text,
                              bool                                      & failed) const;

  /*
  void  SetupIntervalSimbox(ModelSettings                               * model_settings,
                            const Simbox                                * estimation_simbox,
                            Simbox                                      & interval_simboxes,
                            const std::vector<Surface>                  & eroded_surfaces,
                            const std::string                           & corr_dir_single_surf,
                            const std::string                           & corr_dir_top_surf,
                            const std::string                           & corr_dir_base_surf,
                            bool                                          corr_dir_top_conform,
                            bool                                          corr_dir_base_conform,
                            double                                      & desired_grid_resolution,
                            double                                      & relative_grid_resolution,
                            std::string                                 & err_text,
                            bool                                        & failed) const;
 */

  Surface * MakeSurfaceFromFileName(const std::string        & file_name,
                                    const Simbox             & estimation_simbox) const;

  void ErodeSurface(Surface       &  surface,
                    const Surface &  priority_surface,
                    const Surface &  resolution_surface,
                    const bool    &  compare_upward) const;

  void ErodeAllSurfaces(std::vector<Surface>            & eroded_surfaces,
                        const std::vector<int>          & erosion_priorities,
                        const std::vector<Surface>      & surface,
                        const Simbox                    & simbox,
                        std::string                     & err_text) const;

  void  FindSmallestSurfaceGeometry(const double   x0,
                                    const double   y0,
                                    const double   lx,
                                    const double   ly,
                                    const double   rot,
                                    double       & x_min,
                                    double       & y_min,
                                    double       & x_max,
                                    double       & y_max) const;

  void  EstimateZPaddingSize(Simbox          * simbox,
                             ModelSettings   * model_settings) const;

  void  LogIntervalInformation(const Simbox         * simbox,
                               const std::string    & interval_name,
                               const std::string    & header_text1,
                               const std::string    & header_text2) const;

  void  LogIntervalInformation(const Simbox      * simbox,
                               const std::string & header_text1,
                               const std::string & header_text2) const;

  double  FindResolution(const Surface * top_surface,
                         const Surface * base_surface,
                         const Simbox  * estimation_simbox,
                         int             n_layers) const;

  // CLASS VARIABLES
  double                                               dz_min_;                    // Highest vertical resolution in the interval simboxes
  int                                                  n_intervals_;
  bool                                                 multiple_interval_setting_;
  std::vector<std::string>                             interval_names_;
  std::vector<int>                                     erosion_priorities_;
  std::vector<std::string>                             surface_files_;

  std::vector<Simbox *>                                interval_simboxes_;        // Extended simbox with padding and correlation direction, must have same size as the parameters vector


  std::vector<double>                                  desired_grid_resolution_;  // Max vertical distance between original interval surfaces divided by number of layers
  std::vector<double>                                  relative_grid_resolution_; // Actual grid resolution relative to the wanted grid resolution.
  std::vector<double>                                  dz_rel_;                   // Actual grid resolution relative to simbox 0

};

#endif
