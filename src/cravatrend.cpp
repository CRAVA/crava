/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/cravatrend.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
//#include "src/commondata.h"

CravaTrend::CravaTrend()
{
}

CravaTrend::CravaTrend(const Simbox                      * interval_simbox,
                       const std::vector<int>            & trend_cube_type,
                       const std::vector<std::string>    & trend_cube_parameters,
                       std::vector<NRLib::Grid<float> *> & trend_cubes,
                       std::string                       & err_txt)
{
  // Class variables
  (void) err_txt;
  n_samples_      = 1000;
  n_trend_cubes_  = static_cast<int>(trend_cube_parameters.size());

  std::vector<std::string> trend_cube_names(n_trend_cubes_);

  if (n_trend_cubes_ > 0) {

    const int nx   = interval_simbox->getnx();
    const int ny   = interval_simbox->getny();
    const int nz   = interval_simbox->getnz();

    for(int grid_number=0; grid_number<n_trend_cubes_; grid_number++) {

      NRLib::Grid<float> * trend_cube = new NRLib::Grid<float>(nx, ny, nz, 0.0);

      // 1 Trend cube from file ------------------------------------------------
      if(trend_cube_type[grid_number] == ModelSettings::CUBE_FROM_FILE) {
        // Do nothing
        //trend_cube = trend_cubes[grid_number];
        trend_cube = new NRLib::Grid<float>(*trend_cubes[grid_number]);
      }

      // 2 Trend cube from stratigraphic depth  -------------------------------
      else if(trend_cube_type[grid_number] == ModelSettings::STRATIGRAPHIC_DEPTH) {

        LogKit::LogFormatted(LogKit::Low,"\nGenerating trend grid \'"+trend_cube_parameters[grid_number]+"\'\n");

        for(int k=0; k<nz; k++) {
          for(int j=0; j<ny; j++) {
            for(int i=0; i<nx; i++) {
              trend_cube->SetValue(i, j, k, static_cast<float>(k));
            }
          }
        }
      }

      // 3 Trend cube from depth  ------------------------------------------------
      else if(trend_cube_type[grid_number] == ModelSettings::TWT) {

        LogKit::LogFormatted(LogKit::Low,"\nGenerating trend grid \'"+trend_cube_parameters[grid_number]+"\'\n");

        for(int k=0; k<nz; k++) { //nzp
          for(int j=0; j<ny; j++) { //nyp
            for(int i=0; i<nx; i++) { //rnxp
              // value is set to depth from simbox
              float value = static_cast<float>(interval_simbox->getTop(i,j) + interval_simbox->getdz(i,j)*k);
              trend_cube->SetValue(i, j, k, value);
            }
          }
        }
      }

      trend_cubes_.push_back(trend_cube);

      // Calculate trend_cube_sampling_
      // Sample all trends from min to max of the trend cube, using increment_ in the sampling

      float avg = 0.0f;
      float min = 0.0f;
      float max = 0.0f;
      trend_cube->GetAvgMinMax(avg, min, max);
      const double increment = static_cast<double>((max-min)/(n_samples_-1));

      std::vector<double> sampling(n_samples_);

      for(int j=0; j<n_samples_-1; j++)
        sampling[j] = static_cast<double>(min) + j*increment;

      sampling[n_samples_-1] = max;

      trend_cube_sampling_.push_back(sampling);
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
    trend_cube_values[m] = trend_cubes_[m]->GetValue(i,j,k) - trend_cube_sampling_[m][0];

  return trend_cube_values;
}

std::vector<int>
CravaTrend::GetSizeTrendCubes() const
{
  std::vector<int> gridSize(3);

  int nI = static_cast<int>(trend_cubes_[0]->GetNI());
  int nJ = static_cast<int>(trend_cubes_[0]->GetNJ());
  int nK = static_cast<int>(trend_cubes_[0]->GetNK());
  gridSize[0] = nI;
  gridSize[1] = nJ;
  gridSize[2] = nK;

  return gridSize;
}

//void
//CravaTrend::writeToFile(const Simbox        * timeSimbox,
//                        FFTGrid             * grid,
//                        const std::string   & fileName,
//                        const std::string   & sgriLabel)
//{
//
//
//  grid ->setAccessMode(FFTGrid::RANDOMACCESS);
//
//  grid->writeFile(fileName,
//                  IO::PathToInversionResults(),
//                  timeSimbox,
//                  sgriLabel);
//
//  grid->endAccess();
//
//}
