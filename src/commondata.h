/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef COMMONDATA_H
#define COMMONDATA_H

#include <math.h>
#include <string>
#include "src/simbox.h"

class ModelSettings;
class InputFiles;

class CommonData{
public:
  CommonData(ModelSettings    * model_settings,
             InputFiles       * input_files);
  ~ CommonData();

  //SET FUNCTIONS

  void        setEstimationSimbox();
  void        setFullInversionVoulume();

  //GET FUNCTIONS

  Simbox          * getEstimationSimbox()     const { return estimation_simbox_     ;}
  NRLib::Volume   * getFullInversionVolume()  const { return full_inversion_volume_ ;}
  

private:

  void getGeometryFromGridOnFile(const std::string                grid_file,
                                        const TraceHeaderFormat * thf,
                                        SegyGeometry           *& geometry,
                                        std::string             & err_text);

  SegyGeometry * geometryFromCravaFile(const std::string & file_name);

  SegyGeometry * geometryFromStormFile(const std::string    & file_name,
                                       std::string          & err_text,
                                       bool                   scale = false);

  bool createOuterTemporarySimbox(ModelSettings           * model_settings,
                                  InputFiles              * input_files,
                                  std::string             & err_text);

  void writeAreas(const SegyGeometry  * area_params,
                  Simbox              * time_simbox,
                  std::string         & err_text);

  void findSmallestSurfaceGeometry(const double   x0,
                                   const double   y0,
                                   const double   lx,
                                   const double   ly,
                                   const double   rot,
                                   double       & x_min,
                                   double       & y_min,
                                   double       & x_max,
                                   double       & y_max);

  void setSimboxSurfacesMultipleIntervals(Simbox                        *& simbox,
                                                      InputFiles                     * input_files,
                                                      ModelSettings                  * model_settings,
                                                      std::string                    & err_text,
                                                      bool                           & failed);

  void setSimboxSurfacesSingleInterval(Simbox                          *& simbox,
                                                   const std::vector<std::string>   & surf_file,
                                                   ModelSettings                    * model_settings,
                                                   std::string                      & err_text,
                                                   bool                             & failed);

  bool readSeismicData();
  bool readWellData();
  bool blockWellsForEstimation();
  bool setupReflectionMatrixAndTempWavelet();
  bool optimizeWellLocations();
  bool estimateWaveletShape();
  bool estimatePriorCorrelation();
  bool setupEstimationRockPhysics();

  int computeTime(int year, int month, int day) const;

  // Bool variables indicating whether the corresponding data processing
  // succeeded
  bool outer_temp_simbox_;
  bool read_seismic_;
  bool read_wells_;
  bool well_blocks_;
  bool setup_reflection_matrix_;
  bool optimize_well_location_;
  bool wavelet_estimation_shape_;
  bool prior_corr_estimation_;
  bool setup_estimation_rock_physics_;

  Simbox          * estimation_simbox_;
  NRLib::Volume   * full_inversion_volume_;

};

#endif
