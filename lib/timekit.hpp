#ifndef TIMEKIT_H
#define TIMEKIT_H

#include <time.h>

class TimeKit
{
public:
  static void     getTime(double& wall, double& cpu);
  static void     markTime();
  static double   getPassedTime();

private:
  static clock_t timeMark_;
};

#endif
