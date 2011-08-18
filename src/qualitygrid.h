#ifndef QUALITYGRID_H
#define QUALITYGRID_H

#include <vector>

class ModelSettings;
class KrigingData2D;
class CovGrid2D;
class WellData;
class Simbox;
class FFTGrid;
class Vario;
class ModelGeneral;

class QualityGrid 
{

public:
  QualityGrid(const std::vector<double>   pValue,
              WellData                 ** wells,
              const Simbox              * simbox,
              const ModelSettings       * modelSettings,
              ModelGeneral              * modelGeneral);
                      
private:

  void generateProbField(FFTGrid             *& grid, 
                         WellData            ** wells,
                         const Simbox        *  simbox, 
                         const ModelSettings *  modelSettings) const;

  void setupKrigingData2D(std::vector<KrigingData2D> & krigingData, 
                          WellData                  ** wells,
                          const Simbox               * simbox, 
                          const int                    nWells) const; 

  void makeKrigedProbField(std::vector<KrigingData2D> & krigingData, 
                           FFTGrid                   *& grid,
                           const Simbox               * simbox, 
                           const CovGrid2D            & cov, 
                           const bool                   isFile) const;


  CovGrid2D & makeCovGrid2D(const Simbox * simbox,
                            Vario  * vario) const;

  std::vector<float>   wellValue_;
};

#endif
