#ifndef SYNTWELLDATA_H
#define SYNTWELLDATA_H

#include <stdlib.h>

#include <vector>
#include <string>

class SyntWellData{

public:
  SyntWellData(double                      const  & trend1,
               double                      const  & trend2,
               int                                  i,
               int                                  j,
               std::vector<float>          const  & alpha,
               std::vector<float>          const  & beta,
               std::vector<float>          const  & rho,
               std::vector<int>            const  & faciesLog,
               std::vector<std::string>    const  & faciesNames);

  ~SyntWellData();

  int                         getWellLength()       const  {return nBins_       ;}
  double                      getTrend1Val()        const  {return trend1Value_ ;}
  double                      getTrend2Val()        const  {return trend2Value_ ;}
  int                         getTrend1Ind()        const  {return trend1Index_ ;}
  int                         getTrend2Ind()        const  {return trend2Index_ ;}
  int                         getFacies(int   i)    const  {return faciesLog_[i];}
  float                       getAlpha(int    i)    const  {return alpha_[i]    ;}
  float                       getBeta (int    i)    const  {return beta_[i]     ;}
  float                       getRho  (int    i)    const  {return rho_[i]      ;}
  const std::vector<int>    & getFaciesLog()        const  {return faciesLog_   ;}
  const std::vector<float>  & getAlpha()            const  {return alpha_       ;}
  const std::vector<float>  & getBeta()             const  {return beta_        ;}
  const std::vector<float>  & getRho()              const  {return rho_         ;}
  const int                 * getIpos()             const  {return ipos_        ;}
  const int                 * getJpos()             const  {return jpos_        ;}
  const int                 * getKpos()             const  {return kpos_        ;}

private:

  int                                    nBins_;
  int                                    nFacies_;
  double                                 trend1Value_;
  double                                 trend2Value_;
  int                                    trend1Index_;
  int                                    trend2Index_;
  std::vector<float>                     alpha_;
  std::vector<float>                     beta_;
  std::vector<float>                     rho_;
  std::vector<int>                       faciesLog_;
  std::vector<std::string>               faciesNames_;

  int                                  * ipos_;                     ///<
  int                                  * jpos_;                     ///< IJK value
  int                                  * kpos_;                     ///<

};

#endif
