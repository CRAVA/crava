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

  ModelGravityStatic(ModelSettings      *& modelSettings,
                     Simbox              * simbox);

  ~ModelGravityStatic();

  bool                          GetFailed()                const { return failed_                 ;}
  std::vector<bool>             GetFailedDetails()         const { return failed_details_         ;}

  std::vector<float>            GetGravityResponse()       const { return gravity_response_       ;}
  std::vector<float>            GetGravityStdDev()         const { return gravity_std_dev_        ;}

  FFTGrid *                     GetUpscalingKernel()       const { return upscaling_kernel_       ;}
  std::vector<std::vector<std::vector<int> > > GetLagIndex() const { return lag_index_            ;}

  int                           GetNx_upscaled()            const { return nx_upscaled_           ;}
  int                           GetNy_upscaled()            const { return ny_upscaled_           ;}
  int                           GetNz_upscaled()            const { return nz_upscaled_           ;}

  int                           GetNxp_upscaled()           const { return nxp_upscaled_          ;}
  int                           GetNyp_upscaled()           const { return nyp_upscaled_          ;}
  int                           GetNzp_upscaled()           const { return nzp_upscaled_          ;}

  double                        GetDx_upscaled()            const { return dx_upscaled_           ;}
  double                        GetDy_upscaled()            const { return dy_upscaled_           ;}
  double                        GetDz_upscaled()            const { return dz_upscaled_           ;}

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


  static int              SetPaddingSize(int original_nxp, int upscaling_factor);
  static std::vector<int> findClosestFactorableNumber(int leastint);

private:
  bool                      debug_;

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.

  bool before_injection_start_;   ///< Incicator marking if this first vintage is before start of injection

  std::vector<float> observation_location_utmx_;    ///< Vectors to store observation location coordinates
  std::vector<float> observation_location_utmy_;
  std::vector<float> observation_location_depth_;
  std::vector<float> gravity_response_;             ///< Vector to store base line gravity response
  std::vector<float> gravity_std_dev_;              ///< Vector to store base line gravity standard deviation

  // This class holds upscaled grid quantities
  int       nx_upscaled_;             ///< Number of cells in each dimension for upscaled FFTGrids. Set equal to nxp_upscaled.
  int       ny_upscaled_;
  int       nz_upscaled_;

  int       nxp_upscaled_;            ///< Number of cells -including padded region - in each dimension for upscaled FFTGrids
  int       nyp_upscaled_;
  int       nzp_upscaled_;

  double    dx_upscaled_;             ///< Coarse grid cell dimensions.
  double    dy_upscaled_;
  double    dz_upscaled_;

  int       x_upscaling_factor_;      ///< Upscaling factors. User input initially, later changed to true value.
  int       y_upscaling_factor_;
  int       z_upscaling_factor_;


  FFTGrid * upscaling_kernel_;
  std::vector<std::vector<std::vector<int> > > lag_index_;  // eller ha som klassevar i gravimetric inversion-klassen, hører vel mer hjemme der.

  ModelGeneral * modelGeneral_;

  void MakeUpscalingKernel(ModelSettings * modelSettings,
                           Simbox        * fullTimeSimbox);

  void MakeLagIndex(int nx_upscaled, int ny_upscaled, int nz_upscaled);

  void SetUpscaledPaddingSize(ModelSettings * modelSettings);

};

#endif
