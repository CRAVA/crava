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

  bool createOuterTemporarySimbox(ModelSettings  * modelSettings,
                                  InputFiles     * inputFiles);
  bool readSeismicData();
  bool readWellData();
  bool blockWellsForEstimation();
  bool setupReflectionMatrixAndTempWavelet();
  bool optimizeWellLocations();
  bool estimateWaveletShape();
  bool estimatePriorCorrelation();
  bool setupEstimationRockPhysics();
  bool readSeismicData();

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
