/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>

#include "src/parameteroutput.h"
#include "src/krigingdata2d.h"
#include "src/modelsettings.h"
#include "src/qualitygrid.h"
#include "src/kriging2d.h"
#include "src/covgrid2d.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/modelgeneral.h"
#include "src/vario.h"

QualityGrid::QualityGrid(const std::vector<double>   pValue,
                         std::vector<WellData *>     wells,
                         const Simbox              * simbox,
                         const ModelSettings       * modelSettings,
                         ModelGeneral              * modelGeneral)
: wellValue_(0)
{

  const int nWells = modelSettings->getNumberOfWells();

  wellValue_.resize(nWells);

  for (int i=0; i<nWells; i++){
    if (pValue[i] >= 0.05)
      wellValue_[i] = 1;
    else if (pValue[i] < 0.05 && pValue[i] >= 0.005)
      wellValue_[i] = static_cast<float>((pValue[i] - 0.005) / 0.045);
    else
      wellValue_[i] = 0;
  }

  FFTGrid * grid;

  if(modelSettings->getSeismicQualityGridValue() != RMISSING)
    value_ = static_cast<float>(modelSettings->getSeismicQualityGridValue());
  else {
    double tmpsum = 0.0;
    for (int i=0; i< nWells; i++)
      tmpsum += wellValue_[i];
    value_ = static_cast<float>(tmpsum) / static_cast<float>(nWells);

    LogKit::LogFormatted(LogKit::Low, "The value used for weighting under <seismic-quality-grid> is " + NRLib::ToString(value_) + ", this is an average from all fit values from wells.\n");
 }

  generateProbField(grid, wells, simbox, modelSettings);

  std::string fileName = "Seismic_Quality_Grid";
  ParameterOutput::writeToFile(simbox, modelGeneral, modelSettings, grid, fileName, "");

  delete grid;
}

void QualityGrid::generateProbField(FFTGrid              *& grid,
                                    std::vector<WellData *> wells,
                                    const Simbox          * simbox,
                                    const ModelSettings   * modelSettings) const
{

  const int nz     = simbox->getnz();
  const int nWells = modelSettings->getNumberOfWells();

  Vario * vario = modelSettings->getBackgroundVario();
  CovGrid2D & cov = makeCovGrid2D(simbox, vario);

  if(modelSettings->getSeismicQualityGridRange() != RMISSING) {
    float range1 = static_cast<float>(modelSettings->getSeismicQualityGridRange());
    float range2 = range1;
    float angle = vario->getAngle();
    std::string type = vario->getType();
    Vario * vario_new;

    if(type == "Spherical") {
      vario_new = new SphericalVario(range1, range2, angle);
    }
    else
    {
      GenExpVario * vario_tmp = dynamic_cast<GenExpVario *>(vario);
      float power = vario_tmp->getPower();
      vario_new = new GenExpVario(power, range1, range2, angle);
    }
    cov = makeCovGrid2D(simbox, vario_new);
  }

  std::vector<KrigingData2D> krigingData(nz);
  setupKrigingData2D(krigingData, wells, simbox, nWells);

  const bool isFile = modelSettings->getFileGrid();

  makeKrigedProbField(krigingData, grid, simbox, cov, isFile);

  delete &cov;
}

void QualityGrid::setupKrigingData2D(std::vector<KrigingData2D> & krigingData,
                                     std::vector<WellData *>      wells,
                                     const Simbox               * simbox,
                                     const int                    nWells) const
{
  for (int w=0; w<nWells; w++)
  {
    const BlockedLogs * bl = wells[w]->getBlockedLogsConstThick();

    const int   nBlocks = bl->getNumberOfBlocks();
    const int * ipos    = bl->getIpos();
    const int * jpos    = bl->getJpos();
    const int * kpos    = bl->getKpos();

    for (int m=0; m<nBlocks; m++)
    {
      int i = ipos[m];
      int j = jpos[m];
      int k = kpos[m];

      krigingData[k].addData(i, j, wellValue_[w]);
    }
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

  const int    nxp  = nx;
  const int    nyp  = ny;
  const int    nzp  = nz;
  const int    rnxp = 2*(nxp/2 + 1);

  const double x0   = simbox->getx0();
  const double y0   = simbox->gety0();
  const double lx   = simbox->getlx();
  const double ly   = simbox->getly();

  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);

  const float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  grid = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  grid->createRealGrid();
  grid->setType(FFTGrid::PARAMETER);
  grid->setAccessMode(FFTGrid::WRITE);

  for (int k=0; k<nzp; k++)
  {
    // Set constant value for layer
    surface.Assign(value_);

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

CovGrid2D & QualityGrid::makeCovGrid2D(const Simbox * simbox,
                                       Vario  * vario) const
{
  const int nx = simbox->getnx();
  const int ny = simbox->getny();

  const float dx = static_cast<float>(simbox->getdx());
  const float dy = static_cast<float>(simbox->getdy());

  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);

  return(*cov);
}
