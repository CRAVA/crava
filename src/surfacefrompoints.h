/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SURFACEFROMPOINTS_H
#define SURFACEFROMPOINTS_H

#include "src/definitions.h"

class SurfaceFromPoints
{

public:

  SurfaceFromPoints(const std::string & filename);
  ~SurfaceFromPoints(void);

  Surface * GetGriddedSurface(void) { return surface_ ;}

private:

  void ReadFile(const std::string   & filename,
                std::vector<double> & x,
                std::vector<double> & y,
                std::vector<double> & twt,
                std::vector<int>    & il,
                std::vector<int>    & xl);

  void CreateGriddedSurface(Surface                   *& surface,
                            const std::vector<double>  & x,
                            const std::vector<double>  & y,
                            const std::vector<double>  & twt,
                            const std::vector<int>     & il,
                            const std::vector<int>     & xl);

  Surface * surface_;
};

#endif
