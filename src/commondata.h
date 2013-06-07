/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef COMMONDATA_H
#define COMMONDATA_H

#include <math.h>
#include <string>
#include "src/simbox.h"
#include "src/welldata.h"
#include "nrlib/well/well.hpp"
#include "nrlib/segy/segy.hpp"

#include "src/seismicstorage.h"

class ModelSettings;
class ModelGeneral;
class InputFiles;
//class FFTGrid;
//class Surface;

class CommonData{
public:
  CommonData(ModelSettings    * model_settings,
             InputFiles       * input_files);

  ~ CommonData();

  //SET FUNCTIONS

  void        setEstimationSimbox();
  void        setFullInversionVoulume();

  //GET FUNCTIONS

  const Simbox          & getEstimationSimbox()     const { return estimation_simbox_     ;}
  const NRLib::Volume   & getFullInversionVolume()  const { return full_inversion_volume_ ;}



private:

  void getGeometryFromGridOnFile(const std::string         grid_file,
                                 const TraceHeaderFormat * thf,
                                 SegyGeometry           *& geometry,
                                 std::string             & err_text);

  SegyGeometry * geometryFromCravaFile(const std::string & file_name);

  SegyGeometry * geometryFromStormFile(const std::string    & file_name,
                                       std::string          & err_text,
                                       bool                   scale = false);

  bool createOuterTemporarySimbox(ModelSettings           * model_settings,
                                  InputFiles              * input_files,
                                  Simbox                  & estimation_simbox,
                                  NRLib::Volume           & full_inversion_interval,
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

  void setSurfacesMultipleIntervals(Simbox                         & estimation_simbox,
                                    NRLib::Volume                  & full_inversion_volume,
                                    const InputFiles               * input_files,
                                    const ModelSettings            * model_settings,
                                    std::string                    & err_text,
                                    bool                           & failed);

  void setSurfacesSingleInterval(Simbox                           & estimation_simbox,
                                 NRLib::Volume                    & full_inversion_volume,
                                 const std::vector<std::string>   & surf_file,
                                 ModelSettings                    * model_settings,
                                 std::string                      & err_text,
                                 bool                             & failed);

  bool readSeismicData(ModelSettings  * model_settings,
                       InputFiles     * input_files);

  bool       checkThatDataCoverGrid(const SegY   * segy,
                                    float         offset,
                                    const Simbox * timeCutSimbox,
                                    float         guard_zone);

  bool readWellData(ModelSettings  * model_settings,
                    InputFiles     * input_files);

  void        readNorsarWell(const std::string              & wellFileName,
                             NRLib::Well                    & new_well,
                             const std::vector<std::string> & logNames,
                             const std::vector<bool>        & inverseVelocity,
                             bool                             faciesLogGiven,
                             std::string                    & error);

  void          readRMSWell(const std::string              & wellFileName,
                            NRLib::Well                    & new_well,
                            const std::vector<std::string> & logNames,
                            const std::vector<bool>        & inverseVelocity,
                            bool                             faciesLogGiven,
                            std::string                    & error);

  bool blockWellsForEstimation();
  bool setupReflectionMatrixAndTempWavelet(ModelSettings * model_settings,
                                           InputFiles *    input_files);

  float  ** readMatrix(const std::string & fileName,
                              int                 n1,
                              int                 n2,
                              const std::string & readReason,
                              std::string       & errText);

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

  Simbox          estimation_simbox_;
  NRLib::Volume   full_inversion_volume_;

  //std::vector<SeismicStorage> seismic_data_;
  std::map<int, std::vector<SeismicStorage> > seismic_data_;
  std::vector<NRLib::Well>     wells_;

  //Well variables not contained in NRlib::Well
  std::map<std::string, int>                        timemissing_;
  std::map<std::string, double>                     xpos0_;
  std::map<std::string, double>                     ypos0_;
  std::map<std::string, std::vector<std::string> >  faciesnames_;
  std::map<std::string, std::vector<int> >          faciesNr_;
  std::map<std::string, int>                        faciesok_; //Bool?
  std::map<std::string, int>                        nFacies_;

  float                  ** reflectionMatrix_;


  //Simbox                     * estimation_simbox_;
  //NRLib::Volume              * full_inversion_volume_;

};

#endif
