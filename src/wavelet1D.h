/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef WAVELET1D_H
#define WAVELET1D_H

#include "fftw.h"
#include "lib/utils.h"
#include "src/wavelet.h"

#include "src/seismicstorage.h"
#include "src/blockedlogscommon.h"

class CovGrid2D;

class Wavelet1D : public Wavelet {
public:
//Constructors and destructor
  Wavelet1D();

  Wavelet1D(const Simbox                                     * simbox,
            SeismicStorage                                   * seismic_data,
            std::map<std::string, BlockedLogsCommon *>       & mapped_blocked_logs,
            const std::vector<Surface *>                     & estimInterval,
            const ModelSettings                              * modelSettings,
            const NRLib::Matrix                              & reflCoef,
            int                                                iAngle,
            std::vector<Wavelet1D *>                         & well_wavelet,
            int                                              & errCode,
            std::string                                      & errTxt,
            bool                                               writing = true);

  Wavelet1D(const std::string & fileName,
            int                 fileFormat,
            const ModelSettings     * modelSettings,
            const float             * reflCoef,
            float               theta,
            int               & errCode,
            std::string       & errText);

  Wavelet1D(Wavelet * wavelet);

  Wavelet1D(std::vector<float>   vec,
            int                  nzp);

  Wavelet1D(fftw_real * vec,
            int         nz,
            int         nzp,
            double      dz,
            bool        fft_order);

  Wavelet1D(Wavelet          * wavelet,
          int                 difftype);

  Wavelet1D(int                 difftype,
          int                 nz,
          int                 nzp);

  Wavelet1D(const ModelSettings * modelSettings,
            const float         * reflCoef,
            float                 theta,
            float                 peakFrequency,
            int                 & errCode);

  void   shiftAndScale(float shift, float scale);
  void   adjustForAmplitudeEffect(double multiplyer, double Halpha);

  virtual ~Wavelet1D();

// Methods that are virtual in Wavelet

  Wavelet1D   * createWavelet1DForErrorNorm(void);
  Wavelet1D   * createLocalWavelet1D(int i,
                                     int j);
  Wavelet1D   * getGlobalWavelet(){ return this;}

  float         findGlobalScaleForGivenWavelet(const ModelSettings                              * modelSettings,
                                               const Simbox                                     * simbox,
                                               SeismicStorage                                   * seismic_data,
                                               const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                                               std::string                                      & err_txt);

  float         calculateSNRatioAndLocalWavelet(const Simbox                                     * inversion_simbox,
                                                const Simbox                                     * estimation_simbox,
                                                const std::vector<std::vector<double> >          & seis_logs,
                                                const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                                                const ModelSettings                              * modelSettings,
                                                std::string                                      & errText,
                                                int                                              & error,
                                                int                                                number,
                                                Grid2D                                          *& noiseScaled,
                                                Grid2D                                          *& shift,
                                                Grid2D                                          *& gain,
                                                float                                              SNRatio,
                                                float                                              waveletScale,
                                                bool                                               estimateSNRatio,
                                                bool                                               estimateGlobalScale,
                                                bool                                               estimateLocalNoise,
                                                bool                                               estimateLocalShift,
                                                bool                                               estimateLocalScale,
                                                bool                                               estimateWavelet);

private:
  float         findOptimalWaveletScale(fftw_real               ** synt_seis_r,
                                        fftw_real               ** seis_r,
                                        int                        nWells,
                                        int                        nzp,
                                        const std::vector<float> & wellWeight,
                                        float                    & err,
                                        std::vector<float>       & errWell,
                                        std::vector<float>       & scaleOptWell,
                                        std::vector<float>       & errWellOptScale)   const;

  float         findOptimalWaveletScale(fftw_real               ** synt_seis_r,
                                        fftw_real               ** seis_r,
                                        int                        nWells,
                                        int                        nzp,
                                        const std::vector<float> & wellWeight)   const;

  void          findLocalNoiseWithGainGiven(fftw_real                                       ** synt_r,
                                            fftw_real                                       ** seis_r,
                                            int                                                nWells,
                                            int                                                nzp,
                                            const std::vector<float>                         & wellWeight,
                                            float                                            & err,
                                            std::vector<float>                               & errWell,
                                            std::vector<float>                               & errWellOptScale,
                                            std::vector<float>                               & scaleOptWell,
                                            Grid2D                                           * gain,
                                            const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                                            const Simbox                                     * simbox)       const;

  void          estimateLocalGain(const CovGrid2D                                  & cov,
                                  Grid2D                                          *& gain,
                                  const std::vector<float>                         & scaleOptWell,
                                  float                                              globalScale,
                                  const std::vector<int>                           & nActiveData,
                                  const Simbox                                     * simbox,
                                  const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs);

  void          estimateLocalShift(const CovGrid2D                                  & cov,
                                   Grid2D                                          *& shift,
                                   const std::vector<float>                         & shiftWell,
                                   const std::vector<int>                           & nActiveData,
                                   const Simbox                                     * simbox,
                                   const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs);

  void          estimateLocalNoise(const CovGrid2D                                  & cov,
                                   Grid2D                                          *& noiseScaled,
                                   float                                              globalNoise,
                                   const std::vector<float>                         & errWellOptScale,
                                   const std::vector<int>                           & nActiveData,
                                   const Simbox                                     * simbox,
                                   const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs);

  float         shiftOptimal(fftw_real                ** ccor_seis_cpp_r,
                             const std::vector<float>  & wellWeight,
                             const std::vector<float>  & dz,
                             int                         nWells,
                             int                         nzp,
                             std::vector<float>        & shiftWell,
                             float                       maxShift);

  void          multiplyPapolouis(fftw_real                ** vec,
                                  const std::vector<float>  & dz,
                                  int                         nWells,
                                  int                         nzp,
                                  float                       waveletLength,
                                  const std::vector<float>  & wellWeight)    const;

  void          adjustLowfrequency(fftw_real * vec_r,
                                   float       dz,
                                   int         nzp,
                                   float       waveletLength) const;

  void          findWavelet(fftw_real               ** ccor_seis_cpp_r,
                           fftw_real                ** cor_cpp_r,
                           fftw_real                ** wavelet_r,
                           const std::vector<float>  & wellWeight,
                           int                         nWells,
                           int                         nt);

  void          CenterWellWavelet(std::vector<fftw_real> & wavelet) const;

  void          writeDebugInfo(fftw_real ** seis_r,
                               fftw_real ** cor_cpp_r,
                               fftw_real ** ccor_seis_cpp_r,
                               fftw_real ** cpp_r,
                               int          nWells) const;
};

#endif
