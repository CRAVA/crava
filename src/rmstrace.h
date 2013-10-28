/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef RMS_TRACE_H
#define RMS_TRACE_H

#include <string>
#include <vector>

class RMSTrace
{
public:
  RMSTrace(int                 & IL,
           int                 & XL,
           int                 & i_index,
           int                 & j_index,
           double              & utmx,
           double              & utmy,
           std::vector<double> & time,
           std::vector<double> & velocity);

  RMSTrace();

  ~RMSTrace();

  const int                       getIIndex()             const { return i_index_            ;}
  const int                       getJIndex()             const { return j_index_            ;}

  const std::vector<double>       getTime()               const { return time_               ;}
  const std::vector<double>       getVelocity()           const { return velocity_           ;}

private:

  int                 IL_;
  int                 XL_;
  int                 i_index_;                   // i-index in the timeSimbox
  int                 j_index_;                   // j-index in the timeSimbox
  double              utmx_;                      // m
  double              utmy_;                      // m
  std::vector<double> time_;                      // ms
  std::vector<double> velocity_;                  // m/s

};

#endif
