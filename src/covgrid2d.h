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
  float              getCov(int deltai, int deltaj) const;
  void               writeToFile(const std::string & name) const;

private:
  std::vector<float> cov_;
  int                nx_;
  int                ny_;
  double             dx_;
  double             dy_;
};

#endif
