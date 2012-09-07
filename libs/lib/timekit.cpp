#include <iostream>

#include "lib/timekit.hpp"

void
TimeKit::getTime(double& wall, double& cpu)
{
  time_t  tmpwall = time(NULL);
  clock_t tmpcpu  = clock();

  if (tmpwall == static_cast<time_t>(-1))  std::cerr << "Unable to get time()\n";
  if (tmpcpu  == static_cast<clock_t>(-1)) std::cerr << "Unable to get clock()\n";

  if (wall==0 && cpu==0)
  {
    wall = static_cast<double>(tmpwall);
    cpu  = static_cast<double>(tmpcpu);
  }
  else
  {
    wall = static_cast<double>(tmpwall) - wall;
    cpu  = (static_cast<double>(tmpcpu)  - cpu)/CLOCKS_PER_SEC;
  }
}


void
TimeKit::markTime()
{
  timeMark_  = clock();
}


double
TimeKit::getPassedTime()
{
  return(static_cast<double>(clock() - timeMark_)/CLOCKS_PER_SEC);
}

clock_t TimeKit::timeMark_   = 0;

