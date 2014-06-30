/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELAVODYNAMIC_H
#define MODELAVODYNAMIC_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/vario.h"
#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/commondata.h"

struct irapgrid;
class Wavelet;
class Wavelet1D;
class Simbox;
class FFTGrid;
class RandomGen;
class GridMapping;
class InputFiles;
//class ModelAVOStatic;
class ModelGeneral;
class SeismicParametersHolder;

class ModelAVODynamic
{
public:
  ModelAVODynamic(ModelSettings          *& model_settings,
                  //ModelAVOStatic          * model_avo_static,
                  ModelGeneral            * model_general,
                  CommonData              * common_data,
                  SeismicParametersHolder & seismic_parameters,
                  const Simbox            * simbox,
                  int                       t,
                  int                       i_interval);

  ~ModelAVODynamic();

  std::vector<FFTGrid *>        GetSeisCubes()             const { return seis_cubes_                     ;}
  std::vector<Wavelet *>        GetWavelets()              const { return wavelets_                       ;}

  const NRLib::Matrix         & GetAMatrix()               const { return reflection_matrix_              ;}
  Grid2D                      * GetLocalNoiseScale(int i)  const { return local_noise_scale_[i]           ;}
  const std::vector<Grid2D *> & GetLocalNoiseScales()      const { return local_noise_scale_              ;}

  bool                          GetFailed()                const { return failed_                         ;}

  const std::vector<std::vector<float> > & GetAngularCorr() const { return angular_corr_                  ;}

  float                         GetSNRatioAngle(int i)     const { return sn_ratio_[i]                    ;}
  std::vector<float>            GetSNRatio()               const { return sn_ratio_                       ;}
  bool                          GetUseLocalNoise()         const { return use_local_noise_                ;}
  float                         GetAngle(int i)            const { return angle_[i]                       ;}
  //bool                          getEstimateWavelet(int i)  const { return estimateWavelet_[i]             ;}
  //bool                          getMatchEnergies(int i)    const { return matchEnergies_[i]               ;}
  int                           GetNumberOfAngles()        const { return static_cast<int>(angle_.size()) ;}

  float                       * GetThetaDeg()              const { return theta_deg_                      ;}
  float                       * GetDataVariance()          const { return data_variance_                  ;}
  float                       * GetErrorVariance()         const { return error_variance_                 ;}
  float                       * GetModelVariance()         const { return model_variance_                 ;}
  float                       * GetSignalVariance()        const { return signal_variance_                ;}
  float                       * GetTheoSNRatio()           const { return theo_sn_ratio_                  ;}
  double                     ** GetErrThetaCov()           const { return err_theta_cov_                  ;}

  void                          ReleaseGrids();                        // Cuts connection to SeisCube_

private:

  bool             FindTimeGradientSurface(const std::string     & refTimeFile,
                                           const Simbox          * simbox,
                                           NRLib::Grid2D<float>  & refTimeGradX,
                                           NRLib::Grid2D<float>  & refTimeGradY);

  void              ComputeDataVariance(std::vector<FFTGrid *> & seisData,
                                        float                  * dataVariance,
                                        int                      nx,
                                        int                      ny,
                                        int                      nz,
                                        int                      nxp,
                                        int                      nyp,
                                        int                      nzp);

  void              SetupErrorCorrelation(const std::vector<Grid2D *>             & noise_scale,
                                          const float                             * data_variance,
                                          const std::vector<float>                & sn_ratio,
                                          const std::vector<std::vector<float > > & angular_corr,
                                          float                                   * error_variance,
                                          double                                 ** err_theta_cov);

  float             ComputeWDCorrMVar(Wavelet1D* WD,
                                      fftw_real* corrT,
                                      int        nzp);

  void              AddSeismicLogs(std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                   const std::vector<SeismicStorage>          & seismic_data,
                                   const Simbox                               & simbox,
                                   int                                          n_angles);


  int                               number_of_angles_;

  std::vector<Wavelet *>            wavelets_;               ///< Wavelet for angle
  std::vector<FFTGrid *>            seis_cubes_;             ///< Seismic data cubes. Deleted in the commondata destructor

  NRLib::Matrix                     reflection_matrix_;      ///< May specify own Zoeppritz-approximation. Default NULL,
                                                             ///< indicating that standard approximation will be used.

  GridMapping                     * time_depth_mapping_;     ///< Contains both simbox and mapping used for depth conversion

  std::vector<Grid2D *>             local_noise_scale_;      ///< Scale factors for local noise

  bool                              failed_;                 ///< Indicates whether errors occured during construction.

  std::vector<std::vector<float > > angular_corr_;           ///< correlations between angle i and j.

  std::vector<float>                sn_ratio_;
  std::vector<float>                angle_;
  bool                              use_local_noise_;
  int                               this_timelapse_;

  float                           * theta_deg_;
  float                           * data_variance_;
  float                           * error_variance_;
  float                           * model_variance_;
  float                           * signal_variance_;
  float                           * theo_sn_ratio_;      // signal noise ratio from model

  double                         ** err_theta_cov_;

  //GridMapping             * timeCutMapping_;        ///< Simbox and mapping for timeCut*/

  //std::vector<bool>         matchEnergies_;
  //std::vector<bool>         estimateWavelet_;

};

#endif
