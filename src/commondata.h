/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef VARIO_H
#define VARIO_H

#include <math.h>
#include <string>

class CommonData{
public:
  CommonData(void);
  virtual ~CommonData();

private:

  void createOuterTemporarySimbox();
  bool readSeismicData();

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
