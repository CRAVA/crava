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
           double              & utmx,
           double              & utmy,
           std::vector<double> & time,
           std::vector<double> & velocity);

  ~RMSTrace();

  const double                    getUtmx()               const { return utmx_             ; }
  const double                    getUtmy()               const { return utmy_             ; }

  const std::vector<double>       getTime()               const { return time_             ; }

private:

  int                 IL_;
  int                 XL_;
  double              utmx_;                      // m
  double              utmy_;                      // m
  std::vector<double> time_;                      // ms
  std::vector<double> velocity_;                  // m/s

};

#endif
