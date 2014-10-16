#ifndef CRAVATREND_H
#define CRAVATREND_H

#include <stdio.h>

#include "src/simbox.h"
#include "src/fftgrid.h"

class ModelSettings;
class InputFiles;

class CravaTrend
{
public:

  CravaTrend();

  CravaTrend(const Simbox                      * interval_simbox,
             const std::vector<int>            & trend_cube_type,
             const std::vector<std::string>    & trend_cube_parameters,
             std::vector<NRLib::Grid<float> *> & trend_cubes,
             std::string                       & err_txt);

  ~CravaTrend();

  std::vector<int>                             GetSizeTrendCubes() const;

  std::vector<double>                          GetTrendPosition(const int & i,
                                                                const int & j,
                                                                const int & k) const;

  const std::vector<std::vector<double> >    & GetTrendCubeSampling()   const   { return trend_cube_sampling_ ;}

  NRLib::Grid<float>                         * GetTrendCube(int i)      const   { return trend_cubes_[i]      ;}

  void SetTrendCubeSampling(int i, std::vector<double> sampling) { trend_cube_sampling_[i] = sampling ;}

private:

  std::vector<NRLib::Grid<float> *>    trend_cubes_;               ///< Trend cubes used in rock phyiscs prior model
  std::vector<std::vector<double> >    trend_cube_sampling_;       ///< Common sampling for all trends, corresponding to min/max-values of the trend cubes
  int                                  n_samples_;                 ///< Use default value 1000 for number of samples for trends
  int                                  n_trend_cubes_;

};

#endif
