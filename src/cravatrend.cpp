#include "src/cravatrend.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/modelgeneral.h"
#include "src/inputfiles.h"
#include "src/fftgrid.h"

#include <string.h>

CravaTrend::CravaTrend()
{
}

CravaTrend::CravaTrend(Simbox                       * timeSimbox,
                       Simbox                       * timeCutSimbox,
                       ModelSettings                * modelSettings,
                       bool                         & failed,
                       std::string                  & errTxt,
                       const InputFiles             * inputFiles)
{

  n_samples_ = 1000;

  const std::vector<std::string> trend_cube_parameters = modelSettings->getTrendCubeParameters();
  const std::vector<int>         trend_cube_type       = modelSettings->getTrendCubeType();
  n_trend_cubes_                                       = static_cast<int>(trend_cube_parameters.size());

  std::vector<std::string> trendCubeNames(n_trend_cubes_);

  if(n_trend_cubes_ > 0) {

    std::string errorText  = "";

    const int nx   = timeSimbox->getnx();
    const int ny   = timeSimbox->getny();
    const int nz   = timeSimbox->getnz();
    const int nxp  = nx;
    const int nyp  = ny;
    const int nzp  = nz;
    const int rnxp = 2*(nxp/2 + 1);

    for(int grid_number=0; grid_number<n_trend_cubes_; grid_number++) {

      FFTGrid * trend_cube = NULL;

      const std::string   log_name   = "trend cube '"+trend_cube_parameters[grid_number]+"'";

      if(trend_cube_type[grid_number] == ModelSettings::CUBE_FROM_FILE) {

        trendCubeNames[grid_number] = inputFiles->getTrendCube(grid_number);

        const SegyGeometry      * dummy1     = NULL;
        const TraceHeaderFormat * dummy2     = NULL;
        const float               offset     = modelSettings->getSegyOffset(0); //Facies estimation only allowed for one time lapse

        ModelGeneral::readGridFromFile(trendCubeNames[grid_number],
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

        if(errorText != "") {
          errorText += "Reading of file \'"+trendCubeNames[grid_number]+"\' failed\n";
          errTxt    += errorText;
          failed     = true;
        }
      }

      else if(trend_cube_type[grid_number] == ModelSettings::STRATIGRAPHIC_DEPTH) {

        LogKit::LogFormatted(LogKit::Low,"\nGenerating trend grid \'"+trend_cube_parameters[grid_number]+"\'\n");

        trend_cube = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
        trend_cube->createRealGrid();
        trend_cube->setAccessMode(FFTGrid::WRITE);

        for(int k=0; k<nzp; k++) {
          for(int j=0; j<nyp; j++) {
            for(int i=0; i<rnxp; i++) {
              if(i < nx)
                trend_cube->setRealValue(i, j, k, static_cast<float>(k));
              else
                trend_cube->setRealValue(i, j, k, 0);
            }
          }
        }

        trend_cube->endAccess();
      }

      else if(trend_cube_type[grid_number] == ModelSettings::TWT) {

        LogKit::LogFormatted(LogKit::Low,"\nGenerating trend grid \'"+trend_cube_parameters[grid_number]+"\'\n");

        trend_cube = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
        trend_cube->createRealGrid();
        trend_cube->setAccessMode(FFTGrid::WRITE);

        for(int k=0; k<nzp; k++) {
          for(int j=0; j<nyp; j++) {
            for(int i=0; i<rnxp; i++) {
              if(i < nx) {
                float value = static_cast<float>(timeSimbox->getTop(i,j) + timeSimbox->getdz(i,j)*k);
                trend_cube->setRealValue(i, j, k, value);
              }
              else
                trend_cube->setRealValue(i, j, k, 0);
            }
          }
        }
        trend_cube->endAccess();
      }

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

      for(int j=0; j<n_samples_-1; j++)
        sampling[j] = min + j*increment;

      sampling[n_samples_-1] = max;

      trend_cube_sampling_.push_back(sampling);


      if(modelSettings->getOutputGridsOther() && IO::TREND_CUBES > 0) {
        std::string fileName = IO::PrefixTrendCubes() + trend_cube_parameters[grid_number];
        writeToFile(timeSimbox, trend_cube, fileName, "trend cube");
      }

      delete trend_cube;

    }
  }
}

CravaTrend::~CravaTrend()
{
}

std::vector<double>
CravaTrend::GetTrendPosition(const int & i,
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

std::vector<int>
CravaTrend::GetSizeTrendCubes() const
{
  std::vector<int> gridSize(3);

  int nI = static_cast<int>(trend_cubes_[0].GetNI());
  int nJ = static_cast<int>(trend_cubes_[0].GetNJ());
  int nK = static_cast<int>(trend_cubes_[0].GetNK());
  gridSize[0] = nI;
  gridSize[1] = nJ;
  gridSize[2] = nK;

  return gridSize;
}

void
CravaTrend::writeToFile(const Simbox        * timeSimbox,
                        FFTGrid             * grid,
                        const std::string   & fileName,
                        const std::string   & sgriLabel)
{


  grid ->setAccessMode(FFTGrid::RANDOMACCESS);

  grid->writeFile(fileName,
                  IO::PathToInversionResults(),
                  timeSimbox,
                  sgriLabel);

  grid->endAccess();

}
