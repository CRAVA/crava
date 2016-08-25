/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef TIMINGS_H
#define TIMINGS_H

#include "src/definitions.h"
#include "nrlib/iotools/logkit.hpp"

class Timings
{
public:
  static void    reportAll(LogKit::MessageLevels logLevel);
  static void    reportTotal();

  static void    setTimeTotal(double& wall, double& cpu);
  static void    setTimeOuterModellingGrid(double& wall, double& cpu);
  static void    setTimeReadSeismic(double& wall, double& cpu);
  static void    addTimeResamplingSeismic(double& wall, double& cpu);
  static void    setTimeWells(double& wall, double& cpu);
  static void    setTimeWavelets(double& wall, double& cpu);
  static void    setTimePriorExpectation(double& wall, double& cpu);
  static void    setTimePriorCorrelation(double& wall, double& cpu);
  static void    addTimeStochasticModel(double& wall, double& cpu);
  static void    addTimeInversion(double& wall, double& cpu);
  static void    addTimeSimulation(double& wall, double& cpu);
  static void    setTimeFiltering(double& wall, double& cpu);
  static void    setTimeFaciesProb(double& wall, double& cpu);
  static void    setTimeKrigingPred(double& wall, double& cpu);
  static void    setTimeCombineResults(double& wall, double& cpu);
  static void    setTimeWriteResults(double& wall, double& cpu);
  static void    setTimeDummy(double& wall, double& cpu);
  static void    addToTimeKrigingSim(double& wall, double& cpu);

private:
  static void    reportOne(const std::string & text, double cpuThis, double wallThis,
                           double cpuTot, double wallTot, LogKit::MessageLevels logLevel);
  static void    calculateRest(void);

  static double  w_total_;
  static double  c_total_;

  static double  w_rest_;
  static double  c_rest_;

  static double  w_outerModellingGrid_;
  static double  c_outerModellingGrid_;

  static double  w_readseismic_;
  static double  c_readseismic_;

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

  static double  w_kriging_pred_;
  static double  c_kriging_pred_;

  static double  w_kriging_sim_;
  static double  c_kriging_sim_;

  static double  w_combine_results_;
  static double  c_combine_results_;

  static double  w_write_results_;
  static double  c_write_results_;

  static double  w_dummy_;
  static double  c_dummy_;

};

#endif
