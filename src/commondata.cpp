/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/commondata.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/simbox.h"
#include "src/timeline.h"
#include "nrlib/segy/segy.hpp"

CommonData::CommonData(ModelSettings  * modelSettings,
                       InputFiles     * inputFiles){
  
  createOuterTemporarySimbox(modelSettings, inputFiles);

}

CommonData::~CommonData(){
}

bool CommonData::createOuterTemporarySimbox(ModelSettings   * modelSettings,
                                            InputFiles      * inputFiles){

  // Forward modelling
  if (modelSettings->getForwardModeling() == true){
    int eventIndex = 0;
  // Inversion
  }else{
    TimeLine  * timeLine = new TimeLine();


        bool firstGravimetricEvent = true;
        for(int i=0;i<modelSettings->getNumberOfVintages();i++) {
          //Vintages may have both travel time and AVO
          int time = computeTime(modelSettings->getVintageYear(i),
                                 modelSettings->getVintageMonth(i),
                                 modelSettings->getVintageDay(i));
          // Do gravity first
          if(modelSettings->getGravityTimeLapse(i)){
            if(firstGravimetricEvent){
              // Do not save first gravity event in timeline
              firstGravimetricEvent = false;
            }
            else{
              timeLine->AddEvent(time, TimeLine::GRAVITY, i);
            }
          }
          if(modelSettings->getNumberOfAngles(i) > 0) 
            timeLine->AddEvent(time, TimeLine::AVO, i);
        }

      int  eventType;
      int  eventIndex;
      double oldTime;
      timeLine->ReSet();
      timeLine->GetNextEvent(eventType, eventIndex, oldTime);

      // pick first AVO vintage
      while (eventType!=TimeLine::AVO)
        timeLine->GetNextEvent(eventType, eventIndex, oldTime);

     std::vector<float> angles = modelSettings->getAngle(eventIndex);

     SegyGeometry * geometry = NULL;
  }
  

  return true;
}

bool CommonData::readSeismicData(){
  return true;

}

bool CommonData::readWellData(){
  return true;

}

bool CommonData::blockWellsForEstimation(){
  return true;

}

bool CommonData::setupReflectionMatrixAndTempWavelet(){
  return true;

}

bool CommonData::optimizeWellLocations(){
  return true;
}

bool CommonData::estimateWaveletShape(){

}
  
bool CommonData::estimatePriorCorrelation(){

}
  
bool CommonData::setupEstimationRockPhysics(){

}

int
CommonData::computeTime(int year, int month, int day) const
{
  if(year == IMISSING)
    return(0);

  int deltaYear = year-1900; //Ok baseyear.
  int time = 365*deltaYear+deltaYear/4; //Leap years.
  if(month == IMISSING)
    time += 182;
  else {
    std::vector<int> accDays(12,0);
    accDays[1]  = accDays[0]  + 31;
    accDays[2]  = accDays[1]  + 28;
    accDays[3]  = accDays[2]  + 31;
    accDays[4]  = accDays[3]  + 30;
    accDays[5]  = accDays[4]  + 31;
    accDays[6]  = accDays[5]  + 30;
    accDays[7]  = accDays[6]  + 31;
    accDays[8]  = accDays[7]  + 31;
    accDays[9]  = accDays[8]  + 30;
    accDays[10] = accDays[9]  + 31;
    accDays[11] = accDays[10] + 30;

    time += accDays[month];

    if(day == IMISSING)
      time += 15;
    else
      time += day;
  }
  return(time);
}