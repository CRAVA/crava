/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BLOCKED_LOGS_FOR_ROCK_PHYSICS_H
#define BLOCKED_LOGS_FOR_ROCK_PHYSICS_H


class WellData;
class Simbox;

class BlockedLogsForRockPhysics
{
public:
  BlockedLogsForRockPhysics(WellData  * well,
                            Simbox    * simbox);

  ~BlockedLogsForRockPhysics(void);

  std::vector<float>       getAlphaForFacies(const std::string & facies_name);
  std::vector<float>       getBetaForFacies(const std::string & facies_name);
  std::vector<float>       getRhoForFacies(const std::string & facies_name);
  std::vector<float>       getBulkForFacies(const std::string & facies_name);
  std::vector<float>       getShearForFacies(const std::string & facies_name);

private:
  void                     assignToFacies(const float                      * wellLog,
                                          const int                        * faciesLog,
                                          const int                        * faciesNumbers,
                                          std::vector<std::vector<float> > & blockedLog);

  std::vector<std::vector<float> >  alpha_;                    ///<
  std::vector<std::vector<float> >  beta_;                     ///< Raw logs (log-domain)
  std::vector<std::vector<float> >  rho_;                      ///<

  std::vector<std::vector<float> > bulk_modulus_;              ///<
  std::vector<std::vector<float> > shear_modulus_;             ///< Logs calculated from alpha_, beta_ and rho_

  std::vector<std::string>          facies_names_;             ///< Names of the facies in the blocked log
};

#endif
