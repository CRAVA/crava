/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/rmstrace.h"

//---------------------------------------------------------------------------
RMSTrace::RMSTrace(int                 & IL,
                   int                 & XL,
                   double              & utmx,
                   double              & utmy,
                   std::vector<double> & time,
                   std::vector<double> & velocity)
: IL_(IL),
  XL_(XL),
  utmx_(utmx),
  utmy_(utmy),
  time_(time),
  velocity_(velocity)
{
}

RMSTrace::~RMSTrace()
{
}
