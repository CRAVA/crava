#ifndef TIMINGS_H
#define TIMINGS_H

#include "src/definitions.h"
#include "nrlib/iotools/logkit.hpp"

class Timings 
{
public: 
  static void    reportAll(LogKit::MessageLevels logLevel);

  static void    setTimeTotal(double& wall, double& cpu);
  static void    setTimeSeismic(double& wall, double& cpu);
  static void    setTimeResamplingSeismic(double& wall, double& cpu);
  static void    setTimeWells(double& wall, double& cpu);
  static void    setTimeWavelets(double& wall, double& cpu);
  static void    setTimePriorExpectation(double& wall, double& cpu);
  static void    setTimePriorCorrelation(double& wall, double& cpu);
  static void    setTimeStochasticModel(double& wall, double& cpu);
  static void    setTimeInversion(double& wall, double& cpu);
  static void    setTimeSimulation(double& wall, double& cpu);
  static void    setTimeFiltering(double& wall, double& cpu);
  static void    setTimeFaciesProb(double& wall, double& cpu);
  static void    setTimeKriging(double& wall, double& cpu);

private:
  static void    reportOne(const std::string & text, double cpuThis, double wallThis, 
                           double cpuTot, double wallTot, LogKit::MessageLevels logLevel);
  static void    calculateRest(void);

  static double  w_total_;
  static double  c_total_;

  static double  w_rest_;
  static double  c_rest_;

  static double  w_seismic_;
  static double  c_seismic_;

  static double  w_resamplingSeismic_;
  static double  c_resamplingSeismic_;

  static double  w_wells_;
  static double  c_wells_;

  static double  w_wavelets_;
  static double  c_wavelets_;

  static double  w_priorExpectation_;
  static double  c_priorExpectation_;

  static double  w_priorCorrelation_;
  static double  c_priorCorrelation_;

  static double  w_stochasticModel_;
  static double  c_stochasticModel_;

  static double  w_inversion_;
  static double  c_inversion_;

  static double  w_simulation_;
  static double  c_simulation_;

  static double  w_filtering_;
  static double  c_filtering_;

  static double  w_facies_;
  static double  c_facies_;

  static double  w_kriging_;
  static double  c_kriging_;
};

#endif
