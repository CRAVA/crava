/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef BLOCKEDLOGSCOMMON_H
#define BLOCKEDLOGSCOMMON_H



class BlockedLogsCommon
{

public:
  BlockedLogsCommon(const NRLib::Well                * const well_data,
                    const std::vector<std::string>   & cont_logs_to_be_blocked,
                    const std::vector<std::string>   & disc_logs_to_be_blocked,
                    const Simbox                     * const estimation_simbox,
                    bool                               interpolate, 
                    bool                             & failed,
                    std::string                      & err_text);

  ~ BlockedLogsCommon(void);


  //GET FUNCTIONS --------------------------------

  int                                    GetNumberOfBlocks()   const   { return n_blocks_                                                             ;}
  const std::vector<double>            & GetXpos(void)         const   { return x_pos_blocked_                                                        ;}
  const std::vector<double>            & GetYpos(void)         const   { return y_pos_blocked_                                                        ;}
  const std::vector<double>            & GetZpos(void)         const   { return z_pos_blocked_                                                        ;}
  const std::vector<double>            & GetTVD(void)          const   { return z_pos_blocked_                                                        ;}
  const std::vector<double>            & GetTWT(void)          const   { return twt_blocked_                                                          ;}
  bool                                   HasContLog(std::string s)     { return (continuous_logs_blocked_.find(s) != continuous_logs_blocked_.end())  ;}
  bool                                   HasDiscLog(std::string s)     { return (discrete_logs_blocked_.find(s) != discrete_logs_blocked_.end())      ;}

private:

  // FUNCTIONS------------------------------------

  void    RemoveMissingLogValues(const NRLib::Well                            * const well_data,
                                 std::vector<double>                          & x_pos_unblocked,
                                 std::vector<double>                          & y_pos_unblocked,
                                 std::vector<double>                          & z_pos_unblocked,
                                 std::vector<double>                          & twt_unblocked,
                                 std::map<std::string, std::vector<double> >  & continuous_logs_unblocked,
                                 std::map<std::string, std::vector<int> >     & discrete_logs_unblocked,
                                 const std::vector<std::string>               & cont_logs_to_be_blocked,
                                 const std::vector<std::string>               & disc_logs_to_be_blocked,
                                 unsigned int                                 & n_data,
                                 bool                                         & failed,
                                 std::string                                  & err_text);

  void    BlockWell(const Simbox                                        * const estimation_simbox,
                    const std::map<std::string, std::vector<double> >   & continuous_logs_unblocked,
                    const std::map<std::string, std::vector<int> >      & discrete_logs_unblocked,
                    std::map<std::string, std::vector<double> >         & continuous_logs_blocked,
                    std::map<std::string, std::vector<int> >            & discrete_logs_blocked,
                    unsigned int                                          n_data,
                    bool                                                  interpolate,
                    bool                                                & failed,
                    std::string                                         & err_text);

  void    FindSizeAndBlockPointers(const Simbox         * const estimation_simbox,
                                   std::vector<int>     & b_ind);

  void    FindBlockIJK(const Simbox                     * const estimation_simbox,
                       const std::vector<int>           & b_ind);

  void    BlockCoordinateLog(const std::vector<int>     &  b_ind,
                             const std::vector<double>  &  coord,
                             std::vector<double>        &  blocked_coord);

  void    BlockContinuousLog(const std::vector<int>     &  b_ind,
                             const std::vector<double>  &  well_log,
                             std::vector<double>        &  blocked_log);

  void    InterpolateContinuousLog(std::vector<double>   & blocked_log,
                                   int                     start,
                                   int                     end,
                                   int                     index,
                                   float                   rel);

  void  FindOptimalWellLocation(std::map<int, std::vector<SeismicStorage> >   seismic_data,
                                Simbox                     * timeSimbox,
                                float                     ** reflCoef,
                                int                          nAngles,
                                const std::vector<float>   & angleWeight,
                                float                        maxShift,
                                int                          iMaxOffset,
                                int                          jMaxOffset,
                                const std::vector<Surface *> limits,
                                int                        & iMove,
                                int                        & jMove,
                                float                      & kMove);

  // CLASS VARIABLES -----------------------------

  std::string        well_name_;        // Name of well
  unsigned int       n_data_;           // Number of legal data values in well log
  unsigned int       n_blocks_;         // number of blocks

  // Blocked values

  std::vector<double> x_pos_blocked_;  // Blocked x position
  std::vector<double> y_pos_blocked_;  // Blocked y position
  std::vector<double> z_pos_blocked_;  // Blocked z position
  std::vector<double> twt_blocked_;    // Blocked twt value

  std::vector<int> i_pos_;    // Simbox i position for block
  std::vector<int> j_pos_;    // Simbox j position for block
  std::vector<int> k_pos_;    // Simbox k position for block

  int n_continuous_logs_;     // Number of continuous logs
  int n_discrete_logs_;       // Number of discrete logs

  std::map<std::string, std::vector<double>> continuous_logs_blocked_;  // Map between variable name and blocked continuous log
  std::map<std::string, std::vector<int>> discrete_logs_blocked_;       // Map between variable name and blocked discrete log

  // Unblocked values

  std::vector<double> x_pos_unblocked_;
  std::vector<double> y_pos_unblocked_;
  std::vector<double> z_pos_unblocked_;
  std::vector<double> twt_unblocked_;

  std::map<std::string, std::vector<double> > continuous_logs_unblocked_;  // Map between variable name and unblocked continuous log
  std::map<std::string, std::vector<int> > discrete_logs_unblocked_;       // Map between variable name and unblocked discrete log

  int                       n_layers_;                 ///< Number of layers in estimation_simbox
  float                     dz_;                       ///< Simbox dz value for block

  int                       first_M_;                   ///< First well log entry contributing to blocked well
  int                       last_M_;                    ///< Last well log entry contributing to blocked well

  int                       first_B_;                   ///< First block with contribution from well log
  int                       last_B_;                    ///< Last block with contribution from well log

};

#endif
