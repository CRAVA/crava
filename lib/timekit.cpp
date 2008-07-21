#include <iostream>

#include "lib/timekit.hpp"

void
TimeKit::getTime(double& wall, double& cpu)
{
  time_t  tmpwall = time(NULL);
  clock_t tmpcpu  = clock();

  if (tmpwall == (time_t)(-1))  std::cerr << "Unable to get time()\n";
  if (tmpcpu  == (clock_t)(-1)) std::cerr << "Unable to get clock()\n";

  if (wall==0 && cpu==0) 
  {
    wall = (double) tmpwall; 
    cpu  = (double) tmpcpu; 
  }
  else 
  {
    wall = (double) tmpwall - wall;
    cpu  = ((double) tmpcpu  - cpu)/CLOCKS_PER_SEC;
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
  return(double(clock() - timeMark_)/CLOCKS_PER_SEC);
}

clock_t TimeKit::timeMark_   = 0;

