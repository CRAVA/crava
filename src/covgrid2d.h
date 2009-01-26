#ifndef COVGRID2D_H
#define COVGRID2D_H

#include <vector>
#include "src/vario.h"

class CovGrid2D
{
public:
  CovGrid2D(Vario *vario, int nx, int ny, float dx, float dy);
  float getCov(int deltai, int deltaj);
  void writeToFile();

private:
  std::vector<float> cov_;
  int nx_, ny_;
  float dx_, dy_;

};

#endif
