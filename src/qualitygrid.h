/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef QUALITYGRID_H
#define QUALITYGRID_H

#include <vector>

class ModelSettings;
class KrigingData2D;
class CovGrid2D;
class Simbox;
class FFTGrid;
class Vario;
class ModelGeneral;
class BlockedLogsCommon;
class SeismicParametersHolder;

class QualityGrid
{

public:

  QualityGrid(const std::vector<double>                    pValue,
              std::map<std::string, BlockedLogsCommon *> & blocked_wells,
              const Simbox                               * simbox,
              const ModelSettings                        * modelSettings,
              SeismicParametersHolder                    & seismicParameters);

private:

  void generateProbField(FFTGrid                                 *& grid,
                         std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                         const Simbox                             * simbox,
                         const ModelSettings                      * modelSettings) const;

  void setupKrigingData2D(std::vector<KrigingData2D>                 & krigingData,
                          std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                          const Simbox                               * simbox) const;

  void makeKrigedProbField(std::vector<KrigingData2D> & krigingData,
                           FFTGrid                   *& grid,
                           const Simbox               * simbox,
                           const CovGrid2D            & cov,
                           const bool                   isFile) const;

  CovGrid2D & MakeCovGrid2D(const Simbox * simbox,
                            Vario  * vario) const;



  std::vector<float>   wellValue_;

  float value_;
};

#endif
