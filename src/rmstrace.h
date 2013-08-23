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

private:

  int                 IL_;
  int                 XL_;
  double              utmx_;
  double              utmy_;
  std::vector<double> time_;
  std::vector<double> velocity_;

};

#endif
