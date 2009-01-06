#include <iostream>
#include <sstream>
#include <iomanip>

#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "lib/timekit.hpp"

#include "src/definitions.h"
#include "src/timings.h"

void 
Timings::reportAll(void)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"\n***********************************************************************");
  LogKit::LogFormatted(LogKit::MEDIUM,"\n***                           Timings                               ***"); 
  LogKit::LogFormatted(LogKit::MEDIUM,"\n***********************************************************************\n");

  calculateRest();

  LogKit::LogFormatted(LogKit::MEDIUM,"\nSection                              CPU time               Wall time");
  LogKit::LogFormatted(LogKit::MEDIUM,"\n---------------------------------------------------------------------\n");
  reportOne("Seismic data             ", c_seismic_         , w_seismic_          , c_total_, w_total_);
  reportOne("Wells                    ", c_wells_           , w_wells_            , c_total_, w_total_);
  reportOne("Wavelets                 ", c_wavelets_        , w_wavelets_         , c_total_, w_total_);
  reportOne("Prior expection          ", c_priorExpectation_, w_priorExpectation_ , c_total_, w_total_);
  reportOne("Prior correlation        ", c_priorCorrelation_, w_priorCorrelation_ , c_total_, w_total_);
  reportOne("Inversion                ", c_inversion_       , w_inversion_        , c_total_, w_total_);
  reportOne("Simulation               ", c_simulation_      , w_simulation_       , c_total_, w_total_);
  reportOne("Parameter filter         ", c_filtering_       , w_filtering_        , c_total_, w_total_);
  reportOne("Facies probabilities     ", c_facies_          , w_facies_           , c_total_, w_total_);
  reportOne("Kriging                  ", c_kriging_         , w_kriging_          , c_total_, w_total_);
  reportOne("Rest                     ", c_rest_            , w_rest_             , c_total_, w_total_);
  LogKit::LogFormatted(LogKit::MEDIUM,  "---------------------------------------------------------------------\n");
  reportOne("Total                    ", c_total_           , w_total_            , c_total_, w_total_);

  LogKit::LogFormatted(LogKit::LOW,"\nTotal CPU  time used in CRAVA: %6d seconds",   static_cast<int>(c_total_));
  LogKit::LogFormatted(LogKit::LOW,"\nTotal Wall time used in CRAVA: %6d seconds\n", static_cast<int>(w_total_));
}

void 
Timings::reportOne(const std::string & text, double cpuThis, double wallThis, double cpuTot, double wallTot)
{
  if (wallThis < 0.00001) // To omit stupit zero-treatment in ToString()
    wallThis = 0.00001;

  double percentCPU  = 100.0*cpuThis/cpuTot;
  double percentWall = 100.0*wallThis/wallTot;

  if (cpuThis > 0.01 && percentCPU > 0.01) {
    LogKit::LogFormatted(LogKit::MEDIUM,"%s %9.2f  %6.2f ",text.c_str(),cpuThis,percentCPU);
    LogKit::LogMessage(LogKit::MEDIUM,"%   ");
    LogKit::LogFormatted(LogKit::MEDIUM,"  %9.2f  %6.2f ",wallThis,percentWall);
    LogKit::LogMessage(LogKit::MEDIUM,"%\n");
  }
}

void 
Timings::calculateRest(void)
{
  w_rest_ = w_total_ - (w_seismic_ + w_wells_ + w_wavelets_ + w_priorExpectation_ + w_priorCorrelation_ 
                        + w_inversion_ + w_simulation_ + w_filtering_ + w_facies_ + w_kriging_);
  c_rest_ = c_total_ - (c_seismic_ + c_wells_ + c_wavelets_ + c_priorExpectation_ + c_priorCorrelation_ 
                        + c_inversion_ + c_simulation_ + w_filtering_ + w_facies_ + c_kriging_);
}

/*
      IF (TIMCPU.GE.D1     .OR.
     &   (TIMEC.GE.DP01 .AND. TIMCPU.GE.DP01 ) .OR.
     &   (TIMEC.GE.DP01 .AND. (REST.OR.TOTAL))      ) THEN
         ICPUHR = INT(TIMCPU/D3600)
         ICPUMN = INT((TIMCPU - D3600*ICPUHR)/D60)
         ICPUSE = INT(TIMCPU - D3600*ICPUHR - D60*ICPUMN)
      ELSE
         ICPUHR = D0
         ICPUMN = D0
         ICPUSE = D0
         TRESTC = TRESTC + TIMCPU
      END IF
C
      IF (TIMWAL.GE.D1     .OR.
     &    (TIMEW.GE.DP01 .AND. TIMWAL.GE.DP01 ) .OR.
     &    (TIMEW.GE.DP01 .AND. (REST.OR.TOTAL))      ) THEN
         IWALHR = INT(TIMWAL/D3600)
         IWALMN = INT((TIMWAL - D3600*IWALHR)/D60)
         IWALSE = INT(TIMWAL - D3600*IWALHR - D60*IWALMN)
      ELSE
         IWALHR = D0
         IWALMN = D0
         IWALSE = D0
         TRESTW = TRESTW + TIMWAL
      END IF
C
      IF (TIMCPU.GT.D1 .OR. 
     &    TIMWAL.GT.D1 .OR. 
     &    TIMEC.GT.DP01 .AND. TIMCPU.GT.DP01 .OR.
     &    TIMEW.GT.DP01 .AND. TIMWAL.GT.DP01 .OR.
     &    (REST.OR.TOTAL) .AND. (TIMEC.GT.DP01 .OR. TIMEW.GT.DP01)) THEN
         WRITE(LUPRI,'(1X,A7,1X,2(5X,I4.4,2(A,I2.2),4X,F6.2,A))')
     &        TEXT,ICPUHR,':',ICPUMN,':',ICPUSE,TIMEC,' %',
     &             IWALHR,':',IWALMN,':',IWALSE,TIMEW,' %'
      END IF
*/


void
Timings::setTimeTotal(double& wall, double& cpu)
{
    TimeKit::getTime(wall,cpu);
  w_total_ = wall;
  c_total_ = cpu;
}

void    
Timings::setTimeSeismic(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_seismic_ = wall;
  c_seismic_ = cpu;
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
Timings::setTimeInversion(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_inversion_ = wall;
  c_inversion_ = cpu;
}

void    
Timings::setTimeSimulation(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_simulation_ = wall;
  c_simulation_ = cpu;
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
Timings::setTimeKriging(double& wall, double& cpu)
{
  TimeKit::getTime(wall,cpu);
  w_kriging_ = wall;
  c_kriging_ = cpu;
}


double Timings::w_total_            = 0.0;
double Timings::c_total_            = 0.0;

double Timings::w_rest_             = 0.0;
double Timings::c_rest_             = 0.0;

double Timings::w_seismic_          = 0.0;
double Timings::c_seismic_          = 0.0;

double  Timings::w_wavelets_        = 0.0;
double  Timings::c_wavelets_        = 0.0;

double Timings::w_wells_            = 0.0;
double Timings::c_wells_            = 0.0;

double Timings::w_priorExpectation_ = 0.0;
double Timings::c_priorExpectation_ = 0.0;

double Timings::w_priorCorrelation_ = 0.0;
double Timings::c_priorCorrelation_ = 0.0;

double Timings::w_inversion_        = 0.0;
double Timings::c_inversion_        = 0.0;

double Timings::w_simulation_       = 0.0;
double Timings::c_simulation_       = 0.0;

double Timings::w_filtering_        = 0.0;
double Timings::c_filtering_        = 0.0;
                
double Timings::w_facies_           = 0.0;
double Timings::c_facies_           = 0.0;

double Timings::w_kriging_          = 0.0;
double Timings::c_kriging_          = 0.0;
