#include "src/cravatrend.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/modelgeneral.h"
#include "src/inputfiles.h"
#include "src/fftgrid.h"

#include <string.h>

CravaTrend::CravaTrend(Simbox                       * timeSimbox,
                       Simbox                       * timeCutSimbox,
                       ModelSettings                * modelSettings,
                       bool                         & failed,
                       std::string                  & errTxt,
                       const InputFiles             * inputFiles)
{

  n_samples_ = 1000;

  const std::vector<std::string> trend_cube_parameters = modelSettings->getTrendCubeParameters();
  n_trend_cubes_                                       = static_cast<int>(trend_cube_parameters.size());

  std::vector<std::string> trendCubeNames(n_trend_cubes_);

  if(n_trend_cubes_ > 0){
    const SegyGeometry      * dummy1 = NULL;
    const TraceHeaderFormat * dummy2 = NULL;
    const float               offset = modelSettings->getSegyOffset(0); //Facies estimation only allowed for one time lapse

    for(int i=0; i<n_trend_cubes_; i++){

      trendCubeNames[i] = inputFiles->getTrendCube(i);

      FFTGrid           * trend_cube = NULL;
      std::string         errorText  = "";
      const std::string   log_name   = "trend cube '"+trend_cube_parameters[i]+"'";

      ModelGeneral::readGridFromFile(trendCubeNames[i],
                                     log_name,
                                     offset,
                                     trend_cube,
                                     dummy1,
                                     dummy2,
                                     FFTGrid::PARAMETER,
                                     timeSimbox,
                                     timeCutSimbox,
                                     modelSettings,
                                     errorText,
                                     true);

      if(errorText != ""){
        errorText += "Reading of file \'"+trendCubeNames[i]+"\' failed\n";
        errTxt    += errorText;
        failed     = true;
      }

      const int rnxp = trend_cube->getRNxp();
      const int nyp  = trend_cube->getNyp();
      const int nzp  = trend_cube->getNzp();
      const int nx   = trend_cube->getNx();
      const int ny   = trend_cube->getNy();
      const int nz   = trend_cube->getNz();

      NRLib::Grid<double> grid_cube(nx, ny, nz);

      for(int k=0; k<nzp; k++) {
        for(int j=0; j<nyp; j++) {
          for(int i=0; i<rnxp; i++) {
            if (i < nx && j < ny && k < nz)
              grid_cube(i,j,k) = trend_cube->getRealValue(i,j,k);
          }
        }
      }

      trend_cubes_.push_back(grid_cube);

      // Calculate trend_cube_sampling_
      // Sample all trends from min to max of the trend cube, using increment_ in the sampling

      trend_cube->calculateStatistics();

      const float  max       = trend_cube->getMaxReal();
      const float  min       = trend_cube->getMinReal();
      const double increment = (max-min)/(n_samples_-1);

      std::vector<double> sampling(n_samples_);

      for(int j=0; j<n_samples_; j++)
        sampling[j] = min + j*increment;

      trend_cube_sampling_.push_back(sampling);

      delete trend_cube;

    }
  }
}

CravaTrend::CravaTrend()
{
}

CravaTrend::~CravaTrend()
{
}

std::vector<double>  CravaTrend::GetTrendPosition(const int & i,
                                                   const int & j,
                                                   const int & k) const
{
  // Find position of the trend corresponding to the value in (i,j,k) of the trend cubes.
  // As all trends are sampled from min to max of the trend cube, using increment_ in the sampling,
  // the position is obtained by subtracting the minimum value of the sampling from the value of the trend cube in (i,j,k)

  std::vector<double> trend_cube_values(2, RMISSING);

  for(int m=0; m<n_trend_cubes_; m++)
    trend_cube_values[m] = trend_cubes_[m](i,j,k) - trend_cube_sampling_[m][0];

  return trend_cube_values;
}
