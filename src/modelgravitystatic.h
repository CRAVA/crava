#ifndef MODELGRAVITYSTATIC_H
#define MODELGRAVITYSTATIC_H

#include <stdio.h>

#include "src/modelsettings.h"
#include "src/inputfiles.h"


class ModelGravityStatic
{
public:
  ModelGravityStatic(ModelSettings        *& modelSettings,
                     ModelGeneral         *& modelGeneral,
                     const InputFiles      * inputFiles);
  ~ModelGravityStatic();

  bool                          GetFailed()                const { return failed_                 ;}
  std::vector<bool>             GetFailedDetails()         const { return failed_details_         ;}

  std::vector<float>            GetGravityResponse()       const { return gravity_response_       ;}
  std::vector<float>            GetGravityStdDev()         const { return gravity_std_dev_        ;}
  Simbox *                      GetUpscaledSimbox()        const { return upscaled_time_simbox_   ;}
  FFTGrid *                     GetUpscalingKernel()       const { return upscaling_kernel_       ;}
  std::vector<std::vector<std::vector<int> > > GetLagIndex() const { return lag_index_            ;}

  int                           GetNxp_upscaled()           const { return nxp_upscaled_          ;}
  int                           GetNyp_upscaled()           const { return nyp_upscaled_          ;}
  int                           GetNzp_upscaled()           const { return nzp_upscaled_          ;}

  // To be used by ModelGravityDynamic as well
  static void ReadGravityDataFile(const std::string   & fileName,
                                  const std::string   & readReason,
                                  int                   nObs,
                                  int                   nColumns,
                                  std::vector <float> & obs_loc_utmx,
                                  std::vector <float> & obs_loc_utmy,
                                  std::vector <float> & obs_loc_depth,
                                  std::vector <float> & gravity_response,
                                  std::vector <float> & gravity_std_dev,
                                  bool                  failed,
                                  std::string         & errText);


  static int  SetPaddingSize(int original_nxp, int upscaling_factor);
  static std::vector<int> findClosestFactorableNumber(int leastint);

private:
  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.

  bool before_injection_start_;   ///< Incicator marking if this first vintage is before start of injection

  std::vector<float> observation_location_utmx_;    ///< Vectors to store observation location coordinates
  std::vector<float> observation_location_utmy_;
  std::vector<float> observation_location_depth_;
  std::vector<float> gravity_response_;             ///< Vector to store base line gravity response
  std::vector<float> gravity_std_dev_;              ///< Vector to store base line gravity standard deviation

  int nx_upscaled_;               ///< Number of cells in each dimension for upscaled FFTGrids.
  int ny_upscaled_;
  int nz_upscaled_;

  int nxp_upscaled_;              ///< Number of cells -including padded region - in each dimension for upscaled FFTGrids
  int nyp_upscaled_;
  int nzp_upscaled_;

  int       x_upscaling_factor_;  ///< Upscaling factors. User input initially, later changed to true value.
  int       y_upscaling_factor_;
  int       z_upscaling_factor_;

  Simbox  * upscaled_time_simbox_;
  FFTGrid * upscaling_kernel_;
  std::vector<std::vector<std::vector<int> > > lag_index_;  // eller ha som klassevar i gravimetric inversion-klassen, hører vel mer hjemme der.

  ModelGeneral * modelGeneral_;

  void MakeUpscaledTimeSimbox(ModelSettings * modelSettings,
                              Simbox * fullTimeSimbox,
                              int nx_upscaled,
                              int ny_upscaled,
                              int nz_upscaled);

  void MakeUpscalingKernel(ModelSettings * modelSettings,
                           Simbox        * fullTimeSimbox);

  void MakeLagIndex(int nx_upscaled, int ny_upscaled, int nz_upscaled);

  void SetUpscaledPaddingSize(ModelSettings * modelSettings);

};

#endif
