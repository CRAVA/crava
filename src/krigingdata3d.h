/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef KRIGINGDATA3D_H
#define KRIGINGDATA3D_H

class CBWellPt;
class BlockedLogsCommon;

#include <map>

class KrigingData3D
{
public:
  KrigingData3D(int ntot);
  //KrigingData3D(std::vector<WellData *> wells, int nWells, int type);
  KrigingData3D(std::map<std::string, BlockedLogsCommon *> blocked_wells, int type);
  ~KrigingData3D(void);

  CBWellPt  ** getData(void)         const { return data_ ;}
  const int    getNumberOfData(void) const { return nd_   ;}
  //void         addData(const float * blVp,
  //                     const float * blVs,
  //                     const float * blRho,
  //                     const int   * ipos,
  //                     const int   * jpos,
  //                     const int   * kpos,
  //                     const int     nd);
  void         addData(const std::vector<double> & vp,
                       const std::vector<double> & beta,
                       const std::vector<double> & rho,
                       const std::vector<int>    & ipos,
                       const std::vector<int>    & jpos,
                       const std::vector<int>    & kpos,
                       const int                   nd);
  void         divide(void);
  void         writeToFile(const std::string type);

private:
  int          gotBlock(int i, int j, int k) const;

  CBWellPt  ** data_;   // Array of log blocks
  int          nd_;     // Number of log blocks
};

#endif
