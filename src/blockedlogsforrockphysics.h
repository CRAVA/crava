/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BLOCKED_LOGS_FOR_ROCK_PHYSICS_H
#define BLOCKED_LOGS_FOR_ROCK_PHYSICS_H


class WellData;
class Simbox;
class CravaTrend;

class BlockedLogsForRockPhysics
{
public:
  BlockedLogsForRockPhysics(WellData         * well,
                            Simbox           * simbox,
                            const CravaTrend & trend_cubes);

  ~BlockedLogsForRockPhysics(void);

  std::vector<float>          getAlphaForFacies(const std::string & facies_name);
  std::vector<float>          getBetaForFacies(const std::string & facies_name);
  std::vector<float>          getRhoForFacies(const std::string & facies_name);
  std::vector<float>          getBulkForFacies(const std::string & facies_name);
  std::vector<float>          getShearForFacies(const std::string & facies_name);
  std::vector<float>          getPorosityForFacies(const std::string & facies_name);

  const std::vector<double> & getS1(void) {return s1_;}
  const std::vector<double> & getS2(void) {return s2_;}


private:

  void                     calculateBulkShear(const int & nBlocks,
                                              const int & nFacies);

  void                     findTrendPositions(const int        * ipos,
                                              const int        * jpos,
                                              const int        * kpos,
                                              const int        & nBlocks,
                                              const CravaTrend & trend_cubes);

  void                     assignToFacies(const float                      * wellLog,
                                          const int                        * faciesLog,
                                          const int                        * faciesNumbers,
                                          std::vector<std::vector<float> > & blockedLog) const;

  std::vector<std::vector<float> >  alpha_;                    ///<
  std::vector<std::vector<float> >  beta_;                     ///< Raw logs (log-domain)
  std::vector<std::vector<float> >  rho_;                      ///<
  std::vector<std::vector<float> >  porosity_;                 ///<

  std::vector<std::vector<float> >  bulk_modulus_;             ///<
  std::vector<std::vector<float> >  shear_modulus_;            ///< Logs calculated from alpha_, beta_ and rho_

  std::vector<double>               s1_;                       ///< Trend positions corresponding to the first trend cube, same for all facies
  std::vector<double>               s2_;                       ///< Trend positions corresponding to the second trend cube, same for all facies

  std::vector<std::string>          facies_names_;             ///< Names of the facies in the blocked log
};

#endif
