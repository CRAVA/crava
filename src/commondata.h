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
  CommonData(ModelSettings  * modelSettings,
             InputFiles     * inputFiles);
  ~CommonData();

private:

  bool createOuterTemporarySimbox();
  bool readSeismicData();
  bool readWellData();
  bool blockWellsForEstimation();
  bool setupReflectionMatrixAndTempWavelet();
  bool optimizeWellLocations();

  //Bool variables for whether intial data processing went well.
  bool outer_temp_simbox_;
  bool read_seismic_;
  bool read_wells_;
  bool well_blocks_;
  bool setup_reflection_matrix_;
  bool optimize_well_location_;
  bool wavelet_estimation_shape_;
  bool prior_corr_estimation_;
  bool setup_estimation_rock_physics_;

};

#endif
