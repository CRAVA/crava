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
  CommonData(ModelSettings    * modelSettings,
             InputFiles       * inputFiles);

  ~ CommonData();

  //SET FUNCTIONS

  void        setEstimationSimbox();
  void        setFullInversionVoulume();

  //GET FUNCTIONS

  Simbox          * getEstimationSimbox()     const { return estimation_simbox_     ;}
  NRLib::Volume   * getFullInversionVolume()  const { return full_inversion_volume_ ;}


private:

  bool createOuterTemporarySimbox();
  bool readSeismicData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles);
  bool readWellData(ModelSettings  * modelSettings,
                    InputFiles     * inputFiles);
  bool blockWellsForEstimation();
  bool setupReflectionMatrixAndTempWavelet();
  bool optimizeWellLocations();
  bool estimateWaveletShape();
  bool estimatePriorCorrelation();
  bool setupEstimationRockPhysics();

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

  //int seismic_type_; ///< Enum seismic types

  //std::vector<NRLib::SegY *>   segy_files_;
  //std::vector<StormContGrid *> stormgrids_;

  std::vector<SeismicStorage> seismic_data_;

  //FFTGrid                   ** seisCube_;              ///< Seismic data cubes
  std::vector<WellData *>      welldata_;                      ///< Well data
  std::vector<NRLib::Well>     wells_;
  Simbox                     * estimation_simbox_;
  NRLib::Volume              * full_inversion_volume_;

};

#endif
