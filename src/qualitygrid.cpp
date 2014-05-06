/***************************************************************************
* Copyright (C) 2008 by Norwegian Computing Center and Statoil *
***************************************************************************/

#include <string.h>

#include "src/parameteroutput.h"
#include "src/krigingdata2d.h"
#include "src/modelsettings.h"
#include "src/qualitygrid.h"
#include "src/kriging2d.h"
#include "src/covgrid2d.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/modelgeneral.h"
#include "src/vario.h"
#include "src/blockedlogscommon.h"
#include "src/seismicparametersholder.h"

QualityGrid::QualityGrid(const std::vector<double>                    pValue,
                         std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                         const Simbox                               * simbox,
                         const ModelSettings                        * modelSettings,
                         //ModelGeneral                             * modelGeneral,
                         SeismicParametersHolder                    & seismicParameters)
: wellValue_(0)
{

  //const int nWells = modelSettings->getNumberOfWells();
  const int nWells = blocked_wells.size();

  wellValue_.resize(nWells);

  for (size_t i = 0; i < blocked_wells.size(); i++) {
    if (pValue[i] >= 0.05)
      wellValue_[i] = 1;
    else if (pValue[i] < 0.05 && pValue[i] >= 0.005)
      wellValue_[i] = static_cast<float>((pValue[i] - 0.005) / 0.045);
    else
      wellValue_[i] = 0;

  }

  FFTGrid * grid;

  generateProbField(grid, blocked_wells, simbox, modelSettings);

  std::string fileName = "Seismic_Quality_Grid";

  //H-Writing
  //ParameterOutput::writeToFile(simbox, modelGeneral, modelSettings, grid, fileName, "");
  seismicParameters.SetQualityGrid(grid);

  delete grid;
}

void QualityGrid::generateProbField(FFTGrid                                   *& grid,
                                    std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                    const Simbox                               * simbox,
                                    const ModelSettings                        * model_settings) const
{

  const int nz = simbox->getnz();
  //const int nWells = modelSettings->getNumberOfWells();

  std::vector<KrigingData2D> krigingData(nz);

  setupKrigingData2D(krigingData, blocked_wells, simbox);

  Vario * vario = model_settings->getBackgroundVario();
  const CovGrid2D & cov = MakeCovGrid2D(simbox, vario);

  const bool isFile = model_settings->getFileGrid();

  makeKrigedProbField(krigingData, grid, simbox, cov, isFile);

  delete &cov;

}

void QualityGrid::setupKrigingData2D(std::vector<KrigingData2D>                 & krigingData,
                                     std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                     const Simbox                               * simbox) const
{
  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    const int nBlocks = blocked_log->GetNumberOfBlocks();
    const std::vector<int> & ipos = blocked_log->GetIposVector();
    const std::vector<int> & jpos = blocked_log->GetJposVector();
    const std::vector<int> & kpos = blocked_log->GetKposVector();

    for (int m=0; m<nBlocks; m++)
    {
      int i = ipos[m];
      int j = jpos[m];
      int k = kpos[m];

      krigingData[k].addData(i, j, wellValue_[w]);
    }
    w++;
  }

  const int nz = simbox->getnz();

  for (int k=0; k<nz; k++)
    krigingData[k].findMeanValues();
}

void QualityGrid::makeKrigedProbField(std::vector<KrigingData2D> & krigingData,
                                      FFTGrid                   *& grid,
                                      const Simbox               * simbox,
                                      const CovGrid2D            & cov,
                                      const bool                   isFile) const
{
  std::string text = "\nBuilding seismic quality grid:";
  LogKit::LogFormatted(LogKit::Low, text);

  const int nx = simbox->getnx();
  const int ny = simbox->getny();
  const int nz = simbox->getnz();

  const int nxp = nx;
  const int nyp = ny;
  const int nzp = nz;
  const int rnxp = 2*(nxp/2 + 1);

  const double x0 = simbox->getx0();
  const double y0 = simbox->gety0();
  const double lx = simbox->getlx();
  const double ly = simbox->getly();

  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);

  const float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n 0% 20% 40% 60% 80% 100%"
    << "\n | | | | | | | | | | | "
    << "\n ^";

  grid = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  grid->createRealGrid();
  grid->setType(FFTGrid::PARAMETER);
  grid->setAccessMode(FFTGrid::WRITE);

  for (int k=0; k<nzp; k++)
  {
    // Set constant value for layer
    surface.Assign(1.0f);

    // Kriging of layer
    Kriging2D::krigSurface(surface, krigingData[k], cov);

    // Set layer in probability field from surface
    for (int j=0; j<nyp; j++){
      for (int i=0; i<rnxp; i++){
        if (i<nxp)
          grid->setNextReal(float(surface(i,j)));
        else
          grid->setNextReal(0); // Dummy in padding
      }
    }

    // Log process
    if (k+1 >= static_cast<int>(nextMonitor))
    {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }

  grid->endAccess();
}

CovGrid2D & QualityGrid::MakeCovGrid2D(const Simbox * simbox,
                                       Vario * vario) const
{
  const int nx = simbox->getnx();
  const int ny = simbox->getny();

  const float dx = static_cast<float>(simbox->getdx());
  const float dy = static_cast<float>(simbox->getdy());

  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);

  return(*cov);
}
