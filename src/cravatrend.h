#ifndef CRAVATREND_H
#define CRAVATREND_H

#include <stdio.h>

#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/commondata.h"

class ModelSettings;
class InputFiles;
class ModelGeneral;

class CravaTrend
{
public:

  CravaTrend();

  CravaTrend(Simbox                       * timeSimbox,
             Simbox                       * timeCutSimbox,
             ModelSettings                * modelSettings,
             bool                         & failed,
             std::string                  & errTxt,
             const InputFiles             * inputFiles);

  CravaTrend(const Simbox                           * interval_simbox,
             const ModelSettings                    * model_settings,
             const InputFiles                       * input_files,
             const std::string                      & interval_name,
             const std::vector<int>                 & trend_cube_type,
             const std::vector<std::string>         & trend_cube_parameters,
             std::vector<NRLib::Grid<double> >      & trend_cubes,
             //FFTGrid                                * trend_cube,
             //bool                                   & failed,
             std::string                            & err_txt);

  ~CravaTrend();

  std::vector<int>                             GetSizeTrendCubes() const;

  std::vector<double>                          GetTrendPosition(const int & i,
                                                                const int & j,
                                                                const int & k) const;

  const std::vector<std::vector<double> >    & GetTrendCubeSampling()          const   { return trend_cube_sampling_;}

  NRLib::Grid<double>                       &  GetTrendCube(int i)  {return trend_cubes_[i]; }

  void SetTrendCubeSampling(int i, std::vector<double> sampling) {trend_cube_sampling_[i] = sampling;}

private:
  void                                         writeToFile(const Simbox        * timeSimbox,
                                                           FFTGrid             * grid,
                                                           const std::string   & fileName,
                                                           const std::string   & sgriLabel);

  std::vector<NRLib::Grid<double> >    trend_cubes_;               ///< Trend cubes used in rock phyiscs prior model
  std::vector<std::vector<double> >    trend_cube_sampling_;       ///< Common sampling for all trends, corresponding to min/max-values of the trend cubes
  int                                  n_samples_;                 ///< Use default value 1000 for number of samples for trends
  int                                  n_trend_cubes_;

};

#endif
