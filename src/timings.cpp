/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <sstream>
#include <iomanip>

#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "lib/timekit.hpp"

#include "src/definitions.h"
#include "src/timings.h"

void
Timings::reportAll(LogKit::MessageLevels logLevel)
{
  LogKit::WriteHeader("Timings summary", logLevel);

  double c_prediction     = c_inversion_ - c_kriging_pred_;
  double w_prediction     = w_inversion_ - w_kriging_pred_;

  double c_simulation     = c_simulation_ - c_kriging_sim_;
  double w_simulation     = w_simulation_ - w_kriging_sim_;

  double c_kriging_tot    = c_kriging_pred_ + c_kriging_sim_;
  double w_kriging_tot    = w_kriging_pred_ + w_kriging_sim_;

  calculateRest();

  LogKit::LogFormatted(logLevel,"\nSection                                    CPU time               Real time");
  LogKit::LogFormatted(logLevel,"\n-----------------------------------------------------------------------------\n");
  reportOne("Setting up outer modelling grid  ", c_outerModellingGrid_, w_outerModellingGrid_, c_total_, w_total_,logLevel);
  reportOne("Loading seismic data             ", c_readseismic_       , w_readseismic_       , c_total_, w_total_,logLevel);
  reportOne("Resampling seismic data          ", c_resamplingSeismic_ , w_resamplingSeismic_ , c_total_, w_total_,logLevel);
  reportOne("Wells                            ", c_wells_             , w_wells_             , c_total_, w_total_,logLevel);
  reportOne("Wavelets                         ", c_wavelets_          , w_wavelets_          , c_total_, w_total_,logLevel);
  reportOne("Prior expectation                ", c_priorExpectation_  , w_priorExpectation_  , c_total_, w_total_,logLevel);
  reportOne("Prior correlation                ", c_priorCorrelation_  , w_priorCorrelation_  , c_total_, w_total_,logLevel);
  reportOne("Building stochastic model        ", c_stochasticModel_   , w_stochasticModel_   , c_total_, w_total_,logLevel);
  reportOne("Inversion                        ", c_prediction         , w_prediction         , c_total_, w_total_,logLevel);
  reportOne("Simulation                       ", c_simulation         , w_simulation         , c_total_, w_total_,logLevel);
  reportOne("Parameter filter                 ", c_filtering_         , w_filtering_         , c_total_, w_total_,logLevel);
  reportOne("Facies probabilities             ", c_facies_            , w_facies_            , c_total_, w_total_,logLevel);
  reportOne("Kriging                          ", c_kriging_tot        , w_kriging_tot        , c_total_, w_total_,logLevel);
  reportOne("Combining results                ", c_combine_results_   , w_combine_results_   , c_total_, w_total_,logLevel);
  reportOne("Writing results                  ", c_write_results_     , w_write_results_     , c_total_, w_total_,logLevel);
  reportOne("Dummy                            ", c_dummy_             , w_dummy_             , c_total_, w_total_,logLevel);
  reportOne("Miscellaneous                    ", c_rest_              , w_rest_              , c_total_, w_total_,logLevel);
  LogKit::LogFormatted(logLevel,  "-----------------------------------------------------------------------------\n");
  reportOne("Total                            ", c_total_             , w_total_             , c_total_, w_total_,logLevel);
}

void
Timings::reportTotal()
{
  LogKit::LogFormatted(LogKit::Low,"\nTotal CPU  time used in CRAVA: %6d seconds",   static_cast<int>(c_total_));
  LogKit::LogFormatted(LogKit::Low,"\nTotal real time used in CRAVA: %6d seconds\n", static_cast<int>(w_total_));
}

void
Timings::reportOne(const std::string & text, double cpuThis, double wallThis,
                   double cpuTot, double wallTot, LogKit::MessageLevels logLevel)
{
  if (wallThis < 0.00001) // To omit stupit zero-treatment in ToString()
    wallThis = 0.00001;

  double percentCPU  = 100.0*cpuThis/cpuTot;
  double percentWall = 100.0*wallThis/wallTot;

  if (cpuThis > 0.01 && percentCPU > 0.01) {
    LogKit::LogFormatted(logLevel,"%s %9.2f  %6.2f ",text.c_str(),cpuThis,percentCPU);
    LogKit::LogMessage(logLevel,"%   ");
    LogKit::LogFormatted(logLevel,"  %9.2f  %6.2f ",wallThis,percentWall);
    LogKit::LogMessage(logLevel,"%\n");
  }
}

void
Timings::calculateRest(void)
{
  //
  // Note that kriging times are included in c_inversion_ and c_simulation_
  //
  c_rest_ = c_total_ - (c_outerModellingGrid_
                        + c_readseismic_
                        + c_wells_
                        + c_wavelets_
                        + c_priorExpectation_
                        + c_priorCorrelation_
                        + c_stochasticModel_
                        + c_inversion_
                        + c_simulation_
                        + c_filtering_
                        + c_facies_
                        + c_combine_results_
                        + c_write_results_
                        + c_dummy_);
  w_rest_ = w_total_ - (w_outerModellingGrid_
                        + w_readseismic_
                        + w_wells_
                        + w_wavelets_
                        + w_priorExpectation_
                        + w_priorCorrelation_
                        + w_stochasticModel_
                        + w_inversion_
                        + w_simulation_
                        + w_filtering_
                        + w_facies_
                        + w_combine_results_
                        + w_write_results_
                        + w_dummy_);
}

void
Timings::setTimeTotal(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_total_ = wall;
  c_total_ = cpu;
}

void
Timings::setTimeOuterModellingGrid(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_outerModellingGrid_ = wall;
  c_outerModellingGrid_ = cpu;
}

void
Timings::setTimeReadSeismic(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_readseismic_ = wall;
  c_readseismic_ = cpu;
}

void
Timings::addTimeResamplingSeismic(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_resamplingSeismic_ += wall; // Sum times used to resample each cube
  c_resamplingSeismic_ += cpu;
}

void
Timings::setTimeWells(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_wells_ = wall;
  c_wells_ = cpu;
}

void
Timings::setTimeWavelets(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_wavelets_ = wall;
  c_wavelets_ = cpu;
}

void
Timings::setTimePriorExpectation(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_priorExpectation_ = wall;
  c_priorExpectation_ = cpu;
}

void
Timings::setTimePriorCorrelation(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_priorCorrelation_ = wall;
  c_priorCorrelation_ = cpu;
}

void
Timings::addTimeStochasticModel(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_stochasticModel_ += wall;
  c_stochasticModel_ += cpu;
}

void
Timings::addTimeInversion(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_inversion_ += wall;
  c_inversion_ += cpu;
}

void
Timings::addTimeSimulation(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_simulation_ += wall;
  c_simulation_ += cpu;
}

void
Timings::setTimeFiltering(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_filtering_ = wall;
  c_filtering_ = cpu;
}

void
Timings::setTimeFaciesProb(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_facies_ = wall;
  c_facies_ = cpu;
}

void
Timings::setTimeKrigingPred(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_kriging_pred_ = wall;
  c_kriging_pred_ = cpu;
}

void
Timings::addToTimeKrigingSim(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_kriging_sim_ += wall;
  c_kriging_sim_ += cpu;
}

void
Timings::setTimeCombineResults(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_combine_results_ = wall;
  c_combine_results_ = cpu;
}

void
Timings::setTimeWriteResults(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_write_results_ = wall;
  c_write_results_ = cpu;
}

void
Timings::setTimeDummy(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_dummy_ = wall;
  c_dummy_ = cpu;
}


double Timings::w_total_              = 0.0;
double Timings::c_total_              = 0.0;

double Timings::w_rest_               = 0.0;
double Timings::c_rest_               = 0.0;

double Timings::w_outerModellingGrid_ = 0.0;
double Timings::c_outerModellingGrid_ = 0.0;

double Timings::w_readseismic_        = 0.0;
double Timings::c_readseismic_        = 0.0;

double Timings::w_resamplingSeismic_  = 0.0;
double Timings::c_resamplingSeismic_  = 0.0;

double  Timings::w_wavelets_          = 0.0;
double  Timings::c_wavelets_          = 0.0;

double Timings::w_wells_              = 0.0;
double Timings::c_wells_              = 0.0;

double Timings::w_priorExpectation_   = 0.0;
double Timings::c_priorExpectation_   = 0.0;

double Timings::w_priorCorrelation_   = 0.0;
double Timings::c_priorCorrelation_   = 0.0;

double  Timings::w_stochasticModel_   = 0.0;
double  Timings::c_stochasticModel_   = 0.0;

double Timings::w_inversion_          = 0.0;
double Timings::c_inversion_          = 0.0;

double Timings::w_simulation_         = 0.0;
double Timings::c_simulation_         = 0.0;

double Timings::w_filtering_          = 0.0;
double Timings::c_filtering_          = 0.0;

double Timings::w_facies_             = 0.0;
double Timings::c_facies_             = 0.0;

double Timings::w_kriging_pred_       = 0.0;
double Timings::c_kriging_pred_       = 0.0;

double Timings::w_kriging_sim_        = 0.0;
double Timings::c_kriging_sim_        = 0.0;

double Timings::w_combine_results_    = 0.0;
double Timings::c_combine_results_    = 0.0;

double Timings::w_write_results_      = 0.0;
double Timings::c_write_results_      = 0.0;

double Timings::w_dummy_              = 0.0;
double Timings::c_dummy_              = 0.0;
