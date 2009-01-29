#ifndef COVGRID2D_H
#define COVGRID2D_H

#include <vector>
#include "src/vario.h"

class CovGrid2D
{
public:
  CovGrid2D(Vario * vario, 
            int     nx, 
            int     ny, 
            double  dx, 
            double  dy);
  const std::vector<float> & getCov(void);
  float                      getCov(int deltai, int deltaj);
  void                       writeToFile(void);

private:
  std::vector<float>         cov_;
  int                        nx_;
  int                        ny_;
  double                     dx_;
  double                     dy_;
};

#endif
