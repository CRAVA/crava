#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "lib/global_def.h"
#include "lib/lib_misc.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/modelsettings.h"
#include "src/modelfile.h"
#include "src/model.h"
#include "src/vario.h"
#include "src/definitions.h"

ModelFile::ModelFile(char * fileName)
{
  modelSettings_         = new ModelSettings();

  seedFile_              = NULL;
  backFile_              = NULL;
  wellFile_              = NULL;
  headerList_            = NULL;
  timeSurfFile_          = NULL;
  depthSurfFile_         = NULL;
  seismicFile_           = NULL;
  waveletFile_           = NULL;
  reflMatrFile_          = NULL;
  paramCorrFile_         = NULL;
  seisType_              = NULL;  
  waveScale_             = NULL;
  angle_                 = NULL; // Must be kept as class variable due to reuse of readCommandSeismic()
  noiseEnergy_           = NULL; // Must be kept as class variable due to reuse of readCommandSeismic()
  constBack_             = NULL;
  areaParams_            = NULL;

  nSeisData_             =  0;
  nWaveletTransfArgs_    = -1;
  seed_                  =  0;

  time_dTop_             = 0.0; // doubles
  time_lz_               = 0.0;
  time_dz_               = 0.0;
  time_nz_               = 0;

  depth_dTop_            = 0.0; // doubles
  depth_lz_              = 0.0;
  depth_dz_              = 0.0;
  depth_nz_              = 0;

  faciesLogGiven_        = false;
  hasSignalToNoiseRatio_ = false;
  hasDepthSurfaces_      = false;
  parallelTimeSurfaces_  = false;
  parallelDepthSurfaces_ = false;
  generateBackground_    = true;
  failed_                = false;

  FILE* file = fopen(fileName,"r");
  if(file == 0)
  {
    LogKit::LogFormatted(LogKit::LOW,"Error: Could not open file %s for reading.\n", fileName);
    exit(1);
  }
  int nParam = 0;
  char tmpStr[MAX_STRING];
  while(fscanf(file,"%s",tmpStr) != EOF)
  {
    if(tmpStr[0] != '!')
      nParam++;
    else
      readToEOL(file);
  }
  fclose(file);

  char ** params = new char *[nParam+1];
  file = fopen(fileName,"r");
  for(int i=0;i<nParam;i++)
  {  
    fscanf(file,"%s",tmpStr);
    while(tmpStr[0] == '!')
    {
      readToEOL(file);
      fscanf(file,"%s",tmpStr);
    }
    params[i] = new char[strlen(tmpStr)+1];
    strcpy(params[i],tmpStr);
  }
  fclose(file);

  int nCommands = 31;
  bool * genNeed = new bool[nCommands]; //If run in generate mode, this list holds the necessity.

  //The following variables controls the alternative command system.
  //This is used if command A or B or C is necessary.
  //It allows the commands to be both exclusive and non-exclusive.
  //One of the commands is given as needed, the rest as optional.
  //If you intend to use this, and feel uncertain, talk to Ragnar.

  int *  alternativeUsed = new int[nCommands]; 

  //alternativeUsed[i] = k means that command k was used instead of i.
  //alternativeUsed[i] = -k means k used, i can also be used.

  int ** alternative     = new int * [nCommands]; 

  //alternative[i] == NULL means no alternative to command i, otherwise
  //alternative[i][0] is the number of alternatives, with negative sign if non-exclusive.
  //alternative[i][1] to alternative[i][abs(alternative[i][0])+1] is the number for the
  //  alternative commands.

  char * command;
  char ** commandList = new char*[nCommands];
  for(int i=0;i<nCommands;i++)
  {
    commandList[i] = new char[30];
    genNeed[i] = false;
    alternativeUsed[i] = 0;
    alternative[i] = NULL;
  }

  int neededCommands =  3; // The neededCommands first in the list are normally necessary.
  int extraCommands  = 28; // 

  strcpy(commandList[0],"WELLS");
  strcpy(commandList[1],"DEPTH");                                       // ==> TIME_SURFACES eller TIME
  genNeed[1] = true;
  strcpy(commandList[2],"SEISMIC");                                     // ==> PP_SEISMIC
  genNeed[2] = true;
  alternative[2] = new int[2];
  alternative[2][0] = -1;                  //One alternative, non-exclusive.
  alternative[2][1] = neededCommands+extraCommands;   //Number of alternative command.
  //The following commands are optional.
  strcpy(commandList[neededCommands   ],"ANGULARCORRELATION");          // ==> ANGULAR_CORRELATION
  strcpy(commandList[neededCommands+ 1],"SEED");
  strcpy(commandList[neededCommands+ 2],"LATERALCORRELATION");          // ==> LATERAL_CORRELATION
  strcpy(commandList[neededCommands+ 3],"NSIMULATIONS");
  strcpy(commandList[neededCommands+ 4],"PREDICTION");             
  strcpy(commandList[neededCommands+ 5],"PADDING");                     // ==> FFTGRID_PADDING
  strcpy(commandList[neededCommands+ 6],"PREFIX");
  strcpy(commandList[neededCommands+ 7],"AREA");
  genNeed[neededCommands+7] = true;
  strcpy(commandList[neededCommands+ 8],"WHITENOISE");
  strcpy(commandList[neededCommands+ 9],"OUTPUT");
  strcpy(commandList[neededCommands+10],"SEGYOFFSET");
  strcpy(commandList[neededCommands+11],"FORCEFILE");
  strcpy(commandList[neededCommands+12],"DEBUG");
  strcpy(commandList[neededCommands+13],"KRIGING");
  strcpy(commandList[neededCommands+14],"LOCALWAVELET");
  strcpy(commandList[neededCommands+15],"ENERGYTRESHOLD");
  strcpy(commandList[neededCommands+16],"PARAMETERCORRELATION");        // ==> PARAMETER_CORRELATION
  strcpy(commandList[neededCommands+17],"REFLECTIONMATRIX");            // ==> REFLECTION_MATRIX
  strcpy(commandList[neededCommands+18],"FREQUENCYBAND");
  strcpy(commandList[neededCommands+19],"BACKGROUND");
  genNeed[neededCommands+19] = true;
  strcpy(commandList[neededCommands+20],"MAX_DEVIATION_ANGLE");
  strcpy(commandList[neededCommands+21],"GIVESIGNALTONOISERATIO");
  strcpy(commandList[neededCommands+22],"SEISMICRESOLUTION");           // ==> SEISMIC_RESOLUTION
  strcpy(commandList[neededCommands+23],"WAVELETLENGTH");
  strcpy(commandList[neededCommands+24],"DEPTHSURFACES");
  strcpy(commandList[neededCommands+25],"PSSEISMIC");
  strcpy(commandList[neededCommands+26],"PUNDEF");
  strcpy(commandList[neededCommands+27],"ALLOWED_PARAMETER_VALUES");
  //
  // NBNB-PAL: Here Ragnars 32-bit brain goes full... (no room for additional commands)
  //
  //strcpy(commandList[neededCommands+28],"ALLOWED_PARAMETER_VARIANCES");

  char errText[MAX_STRING];
  char ** errorList = new char*[nCommands];
  int error;
  int nErrors = 0;
  if(params[nParam-1][0] != ';')
  {
    params[nParam] = new char[2];
    strcpy(params[nParam],";");
    nParam++;
    strcpy(errText,"Model file does not end with ';'.\n");
    errorList[nErrors] = new char[strlen(errText)+1];
    strcpy(errorList[nErrors], errText);
    nErrors++;
  }
  int  commandsUsed = 0;
  int  comFlag = 0;
  int  curPos = 0;
  bool wrongCommand = false;

  while(curPos < nParam)
  {
    command = uppercase(params[curPos]);
    curPos++;
    int i;
    for(i=0;i<nCommands;i++)
      if(strcmp(commandList[i],command) == 0)
        break;

    if(i < nCommands) 
    {
      comFlag = int(pow(2.0f,i));
      if((comFlag & commandsUsed) > 0) //Bitwise and
      {
        curPos += getParNum(params, curPos, error, errText, " ", 0, -1) + 1;
        error = 1;
        sprintf(errText, "Command %s specified more than once.\n",commandList[i]);
      }
      else if(alternativeUsed[i] > 0) 
      {
        curPos += getParNum(params, curPos, error, errText, " ", 0, -1) + 1;
        error = 1;
        sprintf(errText, "Can not have both command %s and %s.\n", 
                commandList[i], commandList[alternativeUsed[i]]);
      }
      else
      {
        commandsUsed = (commandsUsed | comFlag);
        if(alternative[i] != NULL) {
          int altSign = (alternative[i][0] > 0) ? 1 : -1;
          int j;
          for(j=0;j<altSign*alternative[i][0];j++)
            alternativeUsed[alternative[i][j+1]] = altSign*i;
        }
        switch(i)
        {
        case 0:
          error = readCommandWells(params, curPos, errText);
          break;
        case 1:
          error = readCommandTimeSurfaces(params, curPos, errText);
          break;
        case 2:
          error = readCommandSeismic(params, curPos, errText);
          break;
        case 3:
          error = readCommandAngularCorr(params, curPos, errText);
          break;
        case 4:
          error = readCommandSeed(params, curPos, errText);
          break;
        case 5:
          error = readCommandLateralCorr(params, curPos, errText);
          break;
        case 6:
          error = readCommandSimulate(params, curPos, errText);
          break;
        case 7:
          curPos += getParNum(params, curPos, error, errText, "PREDICTION", 0)+1;
          modelSettings_->setOutputFlag(1);
          break;
        case 8:
          error = readCommandPadding(params, curPos, errText);
          break;
        case 9:
          error = readCommandPrefix(params, curPos, errText);
          break;
        case 10:
          error = readCommandArea(params, curPos, errText);
          break;
        case 11:
          error = readCommandWhiteNoise(params, curPos, errText);
          break;
        case 12:
          error = readCommandOutput(params, curPos, errText);
          break;
        case 13:
          error = readCommandSegYOffset(params, curPos, errText);
          break;
        case 14:
          error = readCommandForceFile(params, curPos, errText);
          break;
        case 15:
          error = readCommandDebug(params, curPos, errText);
          break;
        case 16:
          error = readCommandKriging(params, curPos, errText);
          break;
        case 17:
          error = readCommandLocalWavelet(params, curPos, errText);
          break;
        case 18:
          error = readCommandEnergyTreshold(params, curPos, errText);
          break;
        case 19:
          error = readCommandParameterCorr(params, curPos, errText);
          break;
        case 20:
          error = readCommandReflectionMatrix(params, curPos, errText);
          break;
        case 21:
          error = readCommandFrequencyBand(params, curPos, errText);
          break;
        case 22:
          error = readCommandBackground(params, curPos, errText);
          break;
        case 23:
          error = readCommandMaxDeviationAngle(params,curPos,errText);
          break;
        case 24:
          error = readCommandGiveSignalToNoiseRatios(params, curPos, errText);
          break;
        case 25:
          error = readCommandSeismicResolution(params, curPos, errText);
          break;
        case 26:
          error = readCommandWaveletTaperingL(params, curPos, errText);
          break;
        case 27:
          error = readCommandDepthSurfaces(params, curPos, errText);
          break;
        case 28:
          error = readCommandSeismic(params, curPos, errText, ModelSettings::PSSEIS);
          break;
        case 29:
          error = readCommandPUndef(params,curPos,errText);
          break;
        case 30:
          error = readCommandAllowedParameterValues(params,curPos,errText);
          break;
        case 31:
          error = readCommandAllowedParameterVariances(params,curPos,errText);
          break;
        default:
          sprintf(errText, "Unknown command: %s.\n",command);
          wrongCommand = true;
          curPos += getParNum(params, curPos, error, errText, "-", 0, -1)+1;
          error = 1;
          break;
        }
      }
    }
    else 
    {
      sprintf(errText, "Unknown command: %s.\n",command);
      wrongCommand = true;
      curPos += getParNum(params, curPos, error, errText, "-", 0, -1)+1;
      error = 1;
    }

    if(error > 0)
    {
      errorList[nErrors] = new char[strlen(errText)+1];
      strcpy(errorList[nErrors], errText);
      nErrors++;
    }
    if(strlen(params[curPos-1]) > 1)
    {
      sprintf(errText,"Found %s instead of ';'\n", params[curPos-1]);
      errorList[nErrors] = new char[strlen(errText)+1];
      strcpy(errorList[nErrors], errText);
      nErrors++;
    }
  }
 
  //
  // Check model file validity 
  //
  if(!modelSettings_->getGenerateSeismic()) 
  {
    for(int i=0;i<neededCommands;i++)
    {
      if((commandsUsed % 2) == 0 && alternativeUsed[i] == 0)
      {
        if(alternative[i] == NULL)
          sprintf(errText,"Necessary command %s not specified\n", commandList[i]);
        else {
          sprintf(errText,"Command %s", commandList[i]);
          int j;
          for(j=1;j<abs(alternative[i][0]);j++)
            sprintf(errText,"%s, %s", errText, commandList[alternative[i][j]]);
          sprintf(errText,"%s or %s must be specified.\n", errText, commandList[alternative[i][j]]);
        }
        errorList[nErrors] = new char[strlen(errText)+1];
        strcpy(errorList[nErrors], errText);
        nErrors++;
      }
      commandsUsed /= 2;
    }
  }
  else for(int i=0;i<nCommands;i++)
  {
    if(genNeed[i] == true && ((commandsUsed % 2) == 0))
    {
      sprintf(errText,"Necessary command %s not specified\n", commandList[i]);
      errorList[nErrors] = new char[strlen(errText)+1];
      strcpy(errorList[nErrors], errText);
      nErrors++;
    }
    commandsUsed /= 2;
  }
  delete [] genNeed;

  if((modelSettings_->getOutputFlag() & (ModelSettings::FACIESPROB + ModelSettings::FACIESPROBRELATIVE)) > 0 
     && faciesLogGiven_==false)
  {
    strcpy(errText,"Facies probabilities can not be generated. Please specify a facies log under WELLS.\n");
    errorList[nErrors] = new char[strlen(errText)+1];
    strcpy(errorList[nErrors], errText);
    nErrors++;
  }
  
  if(modelSettings_->getGenerateSeismic() && generateBackground_)
  {
    strcpy(errText,"Background model and seismic cannot both be generated.\n");
    errorList[nErrors] = new char[strlen(errText)+1];
    strcpy(errorList[nErrors], errText);
    nErrors++;
  }
  
  if(generateBackground_ && modelSettings_->getMaxHzBackground() < modelSettings_->getLowCut())
  {
    sprintf(errText,"The frequency high cut for the background (%.1f) must be larger than the frequency low cut for the inversion (%.1f).\n",
            modelSettings_->getMaxHzBackground(),modelSettings_->getLowCut());
    errorList[nErrors] = new char[strlen(errText)+1];
    strcpy(errorList[nErrors], errText);
    nErrors++;
  }
  
  if(nWaveletTransfArgs_ > 0)
  {
    if(nWaveletTransfArgs_ != nSeisData_)
    {
      sprintf(errText,"The number of arguments in the wavelet transformations differ from the number of seismic grids (%d vs %d).\n",
        nWaveletTransfArgs_, nSeisData_);
      errorList[nErrors] = new char[strlen(errText)+1];
      strcpy(errorList[nErrors], errText);
      nErrors++;
    }
  }

  if(nErrors > 0)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nThe following errors were found when parsing %s:\n", fileName);
    for(int i=0;i<nErrors;i++)
    {
      LogKit::LogFormatted(LogKit::LOW,"%s", errorList[i]);
      delete [] errorList[i];
    }
    failed_ = true;    
  }
 
  if (wrongCommand)
  {
    LogKit::LogFormatted(LogKit::LOW,"\nValid commands are:\n\n");
    LogKit::LogFormatted(LogKit::LOW,"General settings:\n");
    LogKit::LogFormatted(LogKit::LOW,"  AREA\n");
    LogKit::LogFormatted(LogKit::LOW,"  DEBUG\n");
    LogKit::LogFormatted(LogKit::LOW,"  FORCEFILE\n");
    LogKit::LogFormatted(LogKit::LOW,"  FREQUENCYBAND\n");
    LogKit::LogFormatted(LogKit::LOW,"  OUTPUT\n");
    LogKit::LogFormatted(LogKit::LOW,"  PADDING\n");
    LogKit::LogFormatted(LogKit::LOW,"  PREFIX\n");
    LogKit::LogFormatted(LogKit::LOW,"  PUNDEF\n");
    LogKit::LogFormatted(LogKit::LOW,"  REFLECTIONMATRIX\n");
    LogKit::LogFormatted(LogKit::LOW,"  SEED\n");
    LogKit::LogFormatted(LogKit::LOW,"Modelling mode:\n");
    LogKit::LogFormatted(LogKit::LOW,"  KRIGING\n");
    LogKit::LogFormatted(LogKit::LOW,"  NSIMULATIONS\n");
    LogKit::LogFormatted(LogKit::LOW,"  PREDICTION\n");
    LogKit::LogFormatted(LogKit::LOW,"Well data:\n");
    LogKit::LogFormatted(LogKit::LOW,"  ALLOWED_PARAMETER_VALUES\n");
    LogKit::LogFormatted(LogKit::LOW,"  ALLOWED_PARAMETER_VARIANCES\n");
    LogKit::LogFormatted(LogKit::LOW,"  MAX_DEVIATION_ANGLE\n");
    LogKit::LogFormatted(LogKit::LOW,"  SEISMICRESOLUTION\n");
    LogKit::LogFormatted(LogKit::LOW,"  WELLS\n");
    LogKit::LogFormatted(LogKit::LOW,"Surface data:\n");
    LogKit::LogFormatted(LogKit::LOW,"  DEPTH\n");
    LogKit::LogFormatted(LogKit::LOW,"  DEPTHSURFACES\n");
    LogKit::LogFormatted(LogKit::LOW,"Seismic data:\n");
    LogKit::LogFormatted(LogKit::LOW,"  ENERGYTRESHOLD\n");
    LogKit::LogFormatted(LogKit::LOW,"  LOCALWAVELET\n");
    LogKit::LogFormatted(LogKit::LOW,"  PSSEISMIC\n");
    LogKit::LogFormatted(LogKit::LOW,"  SEGYOFFSET\n");
    LogKit::LogFormatted(LogKit::LOW,"  SEISMIC\n");
    LogKit::LogFormatted(LogKit::LOW,"  WAVELETLENGTH\n");
    LogKit::LogFormatted(LogKit::LOW,"Prior model:\n");
    LogKit::LogFormatted(LogKit::LOW,"  BACKGROUND\n");
    LogKit::LogFormatted(LogKit::LOW,"  LATERALCORRELATION\n");
    LogKit::LogFormatted(LogKit::LOW,"  PARAMETERCORRELATION\n");
    LogKit::LogFormatted(LogKit::LOW,"Error model:\n");
    LogKit::LogFormatted(LogKit::LOW,"  ANGULARCORRELATION\n");
    LogKit::LogFormatted(LogKit::LOW,"  GIVESIGNALTONOISERATIO\n");
    LogKit::LogFormatted(LogKit::LOW,"  WHITENOISE\n");
  }

  for(int i=0;i<nParam;i++)
    delete [] params[i];
  delete [] params;

  for(int i=0;i<nCommands;i++)
    delete [] commandList[i];
  delete [] commandList;

  for(int i=0;i<nCommands;i++)
    delete [] alternative[i];
  delete [] alternative;
  
  delete [] alternativeUsed;
  delete [] errorList;
}

ModelFile::~ModelFile()
{
  if(wellFile_ != NULL) 
  {
    for(int i=0 ; i < modelSettings_->getNumberOfWells() ; i++)
      delete [] wellFile_[i];
    delete [] wellFile_;
  }

  if(headerList_ != NULL)
  {
    for(int i=0 ; i < 5 ; i++)
      if (headerList_[i] != NULL)
        delete [] headerList_[i];
    delete [] headerList_;
  }

  if(timeSurfFile_ != NULL) 
  {
    delete [] timeSurfFile_[0];   // top surface
    if (timeSurfFile_[1] != NULL)
      delete [] timeSurfFile_[1]; // base surface
    delete [] timeSurfFile_;
  }
  
  if(depthSurfFile_ != NULL) 
  {
    delete [] depthSurfFile_[0];   // top surface
    if (depthSurfFile_[1] != NULL)
      delete [] depthSurfFile_[1]; // base surface
    delete [] depthSurfFile_;
  }

  if (seismicFile_ != NULL)
  {
    for(int i=0;i<nSeisData_;i++)
      if(seismicFile_[i] != NULL)
        delete [] seismicFile_[i];
    delete [] seismicFile_;
  }

  if (waveletFile_ != NULL)
  {
    for(int i=0;i<nSeisData_;i++)
      if(waveletFile_[i] != NULL)
        delete [] waveletFile_[i];
    delete [] waveletFile_;  
  }

  if (areaParams_ != NULL)
    delete [] areaParams_;

  if(seedFile_!=NULL)
    delete seedFile_;

  if (constBack_ != NULL)
    delete [] constBack_;

  if(backFile_!=NULL)
  {
    for(int i=0;i<3;i++)
      if(backFile_[i]!=NULL)
        delete [] backFile_[i];
    delete [] backFile_; 
  }

  if(paramCorrFile_!=NULL)
    delete [] paramCorrFile_;
  delete [] angle_;
  delete [] noiseEnergy_;
  delete [] waveScale_;
  delete [] seisType_;
}

int 
ModelFile::readCommandWells(char ** params, int & pos, char * errText)
{
  sprintf(errText,"%c",'\0');
  int error=0;
  char *commandName = params[pos-1];
  char *command = new char[MAX_STRING];
  strcpy(command,params[pos]);
  if(strcmp(uppercase(command),"HEADERS")==0)
  {
    int nHeader = getParNum(params, pos+1, error, errText, "HEADERS in WELLS", 4,5);
    if(nHeader==5)// facies log present
      faciesLogGiven_ = true;
    // Read headers from model file
    if(error==0)
    { 
      headerList_ = new char*[5]; // changed from 4 to allow for facies 

      for(int i=0 ; i<nHeader ; i++)
      {
        headerList_[i] = new char[MAX_STRING];
        sprintf(headerList_[i],uppercase(params[pos+1+i])); 
      } 
      if(nHeader<5)// no facies log, dummy name
      {
        headerList_[4] = new char[MAX_STRING]; 
        strcpy(headerList_[4],"FACIES");
      }
    }
    pos+= nHeader+2;
  }
  delete [] command;

  int tmperror = 0;
  char tmpErrText[MAX_STRING];

  int nWells = getParNum(params, pos, tmperror, tmpErrText, commandName, 1, -1);

  if(tmperror > 0)
  {
    sprintf(errText,"%s%s",errText,tmpErrText);
    pos += 1;
    return(1);
  }

  wellFile_ = new char*[nWells];

  for(int i=0 ; i<nWells ; i++)
  {
    wellFile_[i] = new char[strlen(params[pos+i])+1];
    strcpy(wellFile_[i],params[pos+i]);
  }
  modelSettings_->setNumberOfWells(nWells);

  pos += nWells+1;
  return(error);
}

int 
ModelFile::readCommandBackground(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 2, -1);
  if(error == 1)
  {
    pos += nPar+1;
    return(1);
  }
  //
  // If the number of parameters are exactly 3 we assume the background model
  // to be given either on file or as constants.
  //
  if (nPar == 3) 
  {
    generateBackground_ = false;

    backFile_  = new char *[3];
    constBack_ = new float[3];

    for(int i=0 ; i<3 ; i++)
    {
      if(isNumber(params[pos+i])) // constant background
      { 
        constBack_[i] = (float)(atof(params[pos+i]));  
      }
      else // not constant, read from file.
      {
        backFile_[i] = new char[strlen(params[pos+i])+1];
        strcpy(backFile_[i], params[pos+i]);
        if(checkFileOpen(&(backFile_[i]), 1, params[pos-1], errText) > 0) 
        {
          error = 1;
          delete [] backFile_[i];
          backFile_[i]  = NULL; //Indicates problem with reading file.
          constBack_[i] = -1;  //Indicates model should be read from file.
        }
        else
          constBack_[i] = float(checkFileType(backFile_[i]));
      }
    }
  }
  else
  {
    //
    // Currently we have to generate all or none parameters. The reason
    // is that the kriging algorithm handles all parameters simulataneously.
    //
    generateBackground_ = true;

    int nSubCommands = 2;
    char ** subCommand = new char * [nSubCommands];
    for(int i=0 ; i<nSubCommands ; i++)
      subCommand[i] = new char[30];
    strcpy(subCommand[0],"VARIOGRAM");
    strcpy(subCommand[1],"FREQUENCY");
    int subCom, comPar, curPar = 0;
    while(curPar < nPar && error == 0) {
      comPar = 1;
      while(isNumber(params[pos+curPar+comPar]))
        comPar++;
      subCom = 0;
      while(subCom < nSubCommands && strcmp(params[pos+curPar], subCommand[subCom]) != 0)
        subCom++;
      switch(subCom) {
      case 0:
        if(comPar > 1) {
          error = 1;
          sprintf(errText,"Subcommand VARIOGRAM in command %s needs variogram type.\n",params[pos-1]);
        }
        else {
          comPar = 2;
          while(isNumber(params[pos+curPar+comPar]))
            comPar++;
          if(strcmp(params[pos+curPar],"VARIOGRAM") == 0) {

            Vario * vario = createVario(&(params[pos+curPar+1]), comPar-1, params[pos-1], errText);
            if(vario == NULL)
            {
              error = 1;
            }
            else
            {
              modelSettings_->setBackgroundVario(vario);
              //
              // NBNB-PAL: Temporary setting until 2D-kriging of BG allows other variograms to be used 
              //
              if(strcmp(params[pos+curPar+1],"GENEXP") != 0)
              {
                sprintf(errText,"Subcommand VARIOGRAM in command %s only allows variogram type \'genexp\'.",
                        params[pos-1]);
                error = 1;
              }
            }
          }
        }
        break;
      case 1:
        if(comPar != 2) {
          error = 1;
          sprintf(errText,"Subcommand FREQUENCY in command %s takes 1 parameter, not %d as was given.\n",
                  params[pos-1], comPar-1);
        }
        else 
          modelSettings_->setMaxHzBackground(float(atof(params[pos+curPar+1])));
        break;
      default: 
        error = 1;
        sprintf(errText,"Unknown subcommand %s found in command %s.\n", params[pos+curPar],params[pos-1]);
        break;
      }
      curPar += comPar;
    }
    for(int i=0;i<nSubCommands;i++)
      delete [] subCommand[i];
    delete [] subCommand;
  }
  pos += nPar+1;
  return(error);
}

int 
ModelFile::readCommandBackgroundControl(char ** params, int & pos, char * errText)
{
  //printf("\nDont used readCommandBackgroundControl....");
  //printf("\nDont used readCommandBackgroundControl....");
  //printf("\nDont used readCommandBackgroundControl....\n");
  //exit(1);

  generateBackground_ = true;

  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 2, -1);
  if(error == 0)
  {
    int i, nSubCommands = 2;
    char ** subCommand = new char * [nSubCommands];
    for(i=0;i<nSubCommands;i++)
      subCommand[i] = new char[30];
    strcpy(subCommand[0],"VARIOGRAM");
    strcpy(subCommand[1],"FREQUENCY");
    int subCom, comPar, curPar = 0;
    while(curPar < nPar && error == 0) {
      comPar = 1;
      while(isNumber(params[pos+curPar+comPar]))
        comPar++;
      subCom = 0;
      while(subCom < nSubCommands && strcmp(params[pos+curPar], subCommand[subCom]) != 0)
        subCom++;
      switch(subCom) {
        case 0:
          if(comPar > 1) {
            error = 1;
            sprintf(errText,"Subcommand VARIOGRAM in command %s needs variogram type.\n",params[pos-1]);
          }
          else {
            comPar = 2;
            while(isNumber(params[pos+curPar+comPar]))
              comPar++;
            if(strcmp(params[pos+curPar],"VARIOGRAM") == 0) {
              Vario * vario = createVario(&(params[pos+curPar+1]), comPar-1, params[pos-1], errText);
              if(vario == NULL)
                error = 1;
              else
                modelSettings_->setBackgroundVario(vario);
            }
          }
          break;
        case 1:
          if(comPar != 2) {
            error = 1;
            sprintf(errText,"Subcommand FREQUENCY in command %s takes 1 parameter, not %d as was given.\n",
              params[pos-1], comPar-1);
          }
          else 
            modelSettings_->setMaxHzBackground(float(atof(params[pos+curPar+1])));
          break;
        default: 
          error = 1;
          sprintf(errText,"Unknown subcommand %s found in command %s.\n", params[pos+curPar],params[pos-1]);
          break;
      }
      curPar += comPar;
    }
    for(i=0;i<nSubCommands;i++)
      delete [] subCommand[i];
    delete [] subCommand;
  }
  pos += nPar+1;

  return(error);
}

int 
ModelFile::readCommandArea(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 7);

  if(error==1)
  {
    pos += nPar+1;
    return(1);
  }
  areaParams_ = new double[7];
  areaParams_[0] = static_cast<double>(atof(params[pos]));            // x0 
  areaParams_[1] = static_cast<double>(atof(params[pos+1]));          // y0 
  areaParams_[2] = static_cast<double>(atof(params[pos+2]));          // lx 
  areaParams_[3] = static_cast<double>(atof(params[pos+3]));          // ly 
  areaParams_[4] = static_cast<double>(PI*atof(params[pos+4])/180.0); // rot
  areaParams_[5] = static_cast<double>(atof(params[pos+5]));          // dx 
  areaParams_[6] = static_cast<double>(atof(params[pos+6]));          // dy 

  pos += nPar+1;
  return(error);
}


int 
ModelFile::readCommandTimeSurfaces(char ** params, int & pos, char * errText)
{
  int error = 0;
  int nPar  = getParNum(params, pos, error, errText, params[pos-1], 3,4);

  timeSurfFile_    = new char*[2]; // top and base surface
  timeSurfFile_[0] = new char[strlen(params[pos])+1];
  strcpy(timeSurfFile_[0],params[pos]);

  if(isNumber(params[pos+1])) //Only one reference surface
  {
    timeSurfFile_[1] = NULL;
    if(nPar == 4)
    {
      parallelTimeSurfaces_  = true;
      time_dTop_ = static_cast<double>(atof(params[pos+1]));
      time_lz_   = static_cast<double>(atof(params[pos+2]));
      time_dz_   = static_cast<double>(atof(params[pos+3]));
    }
    else
    {
      sprintf(errText,"Command %s with only one grid takes %d arguments; %d given in file.\n", 
              params[pos-1], 4, nPar);
      error = 1;
    }
  }
  else
  {
    if(nPar == 3)
    {
      timeSurfFile_[1] = new char[strlen(params[pos+1])+1];
      strcpy(timeSurfFile_[1],params[pos+1]);
      time_nz_ = atoi(params[pos+2]);
    }
    else
    {
      sprintf(errText,"Command %s with two grids takes %d arguments; %d given in file.\n", 
              params[pos-1], 3, nPar);
      error = 1;
    }
  }

  pos += nPar+1;
  return(error);
}

int 
ModelFile::readCommandDepthSurfaces(char ** params, int & pos, char * errText)
{
  hasDepthSurfaces_ = true;

  int error = 0;
  int nPar= getParNum(params, pos, error, errText, params[pos-1], 3,4);

  depthSurfFile_    = new char*[2]; // top and base surface
  depthSurfFile_[0] = new char[strlen(params[pos])+1];
  strcpy(depthSurfFile_[0],params[pos]);

  if(isNumber(params[pos+1])) //Only a reference surface
  {
    if(nPar == 4)
    {
      parallelDepthSurfaces_ = true;
      depth_dTop_ = static_cast<double>(atof(params[pos+1]));
      depth_lz_   = static_cast<double>(atof(params[pos+2]));
      depth_dz_   = static_cast<double>(atof(params[pos+3]));
    }
    else
    {
      sprintf(errText,"Command %s with only one grid takes %d arguments; %d given in file.\n", 
              params[pos-1], 4, nPar);
      error = 1;
    }
  }
  else
  {
    if(nPar == 3)
    {
      depthSurfFile_[1] = new char[strlen(params[pos+1])+1];
      strcpy(depthSurfFile_[1],params[pos+1]);
      depth_nz_ = atoi(params[pos+2]);
    }
    else
    {
      sprintf(errText,"Command %s with two grids takes %d arguments; %d given in file.\n", 
              params[pos-1], 3, nPar);
      error = 1;
    }
  }
  pos += nPar+1;
  return(error);
}

int 
ModelFile::readCommandSeismic(char ** params, int & pos, char * errText, int seisType)
{
  int nCol = 5;
  int i, error;
  int nSeisData = getParNum(params, pos, error, errText, params[pos-1], 2, -1);

  if(error == 0)
  {
    if((nSeisData % nCol) != 0)
    {
      sprintf(errText, "Number of parameters in command %s must be multiple of %d (%d were found).\n", 
        params[pos-1], nCol, nSeisData);
      error = 1;
    }
    else
      nSeisData = nSeisData/nCol;
  }
  if(error == 1)
  {
    pos += nSeisData+1;
    return(1);
  }

  char  ** seisFile        = new char*[2*nSeisData];
  char  ** waveFile        = new char*[nSeisData];  
  float *  angle           = new float[nSeisData];
  float *  noiseEnergy     = new float[nSeisData];
  float *  waveScale       = new float[nSeisData];

  for(i=0;i<nSeisData;i++)
  {
    // 1. file containing seismic
    seisFile[i] = new char[strlen(params[pos+nCol*i])+1];
    strcpy(seisFile[i], params[pos+nCol*i]);
    if(seisFile[i][0] == '?')
      modelSettings_->setGenerateSeismic(true);

    // 2. angle
    angle[i] = float(atof(params[pos+nCol*i+1])*PI/180.0);

    // 3. noiseEnergy (possibly signalToNoiseRatio)
    if (params[pos+nCol*i+2][0] == '*') 
      noiseEnergy[i] = RMISSING;
    else
      noiseEnergy[i] = float(atof(params[pos+nCol*i+2]));
    
    // 4. file containing wavelet
    waveFile[i] = new char[strlen(params[pos+nCol*i+3])+1];   
    strcpy(waveFile[i], params[pos+nCol*i+3]);

    // 5. wavelet scale
    if (params[pos+nCol*i+4][0] == '*') 
      waveScale[i] = RMISSING;
    else
      waveScale[i] = float(atof(params[pos+nCol*i+4]));
  }

  if(modelSettings_->getGenerateSeismic() == false) 
    error = checkFileOpen(seisFile, nSeisData, "SEISMIC", errText);

  if(modelSettings_->getGenerateSeismic() || error == 0)
    for(i=0;i<nSeisData;i++)
      if (waveFile[i][0] != '*' || waveScale[i] == RMISSING)
        sprintf(errText, "For seismic data generation both wavelet and wavelet scale must be given\n");

  if (error == 0)
  {
    int flag = int(pow(2.0f,nSeisData));
    for(i=0;i<nSeisData;i++)
    {
      if((error & flag) == 0 && waveFile[i][0] != '*')
        error = checkFileOpen(&(waveFile[i]), 1, "WAVELET", errText);
      flag *= 2;
    }
  }

  if(nSeisData_ == 0) 
  {
    seismicFile_     = seisFile;
    waveletFile_     = waveFile;
    noiseEnergy_     = noiseEnergy;
    waveScale_       = waveScale;
    angle_           = angle;
    nSeisData_       = nSeisData;
    seisType_        = new int[nSeisData];
    for(i=0;i<nSeisData;i++)
      seisType_[i] = seisType;
  }
  else 
  {
    char    ** seisFile2        = seismicFile_;
    char    ** waveFile2        = waveletFile_;
    float   *  noiseEnergy2     = noiseEnergy_;
    float   *  waveScale2       = waveScale_;
    float   *  angle2           = angle_;
    int     *  seisType2        = seisType_;

    int nTot = nSeisData + nSeisData_;
    seismicFile_     = new char *[nTot];
    waveletFile_     = new char *[nTot];
    noiseEnergy_     = new float[nTot];
    waveScale_       = new float[nTot];
    angle_           = new float[nTot];
    seisType_        = new int[nTot];

    //Copy old values.
    for(i=0;i<nSeisData_;i++)  // Note nSeisData_ not updated yet.
    {
      seismicFile_[i] = seisFile2[i];
      waveletFile_[i] = waveFile2[i];
      noiseEnergy_[i] = noiseEnergy2[i];
      waveScale_[i]   = waveScale2[i];
      angle_[i]       = angle2[i];
      seisType_[i]    = seisType2[i];
    }
    delete [] seisFile2;
    delete [] waveFile2;
    delete [] angle2;
    delete [] noiseEnergy2;
    delete [] waveScale2;
    delete [] seisType2;

    //Insert new values:
    for(i=0;i<nSeisData;i++) 
    {
      seismicFile_[i+nSeisData_] = seisFile[i];
      waveletFile_[i+nSeisData_] = waveFile[i];
      noiseEnergy_[i+nSeisData_] = noiseEnergy[i];
      waveScale_[i+nSeisData_]   = waveScale[i];
      angle_[i+nSeisData_]       = angle[i];
      seisType_[i+nSeisData_]    = seisType;
    }
    delete [] seisFile;
    delete [] waveFile;
    delete [] angle;
    delete [] noiseEnergy;
    delete [] waveScale;

    nSeisData_ += nSeisData;
  }

  modelSettings_->setNumberOfAngles(nSeisData_);
  modelSettings_->setAngle(angle_,nSeisData_);
  modelSettings_->setNoiseEnergy(noiseEnergy_,nSeisData_);
  modelSettings_->setMatchEnergies(waveScale_,nSeisData_);

  pos += nCol*nSeisData+1;
  return(error);
}

int 
ModelFile::readCommandAngularCorr(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 2, 3);
  if(error == 0)
  {
    Vario * vario = createVario(&(params[pos]), nPar, params[pos-1], errText);
    if(vario != NULL)
    {
      vario->convertRangesFromDegToRad();
      modelSettings_->setAngularCorr(vario);
    }
    else
      error = 1;
  }
  pos += nPar+1;

  return(error);
}

int 
ModelFile::readCommandLateralCorr(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 4,5);
  if(error == 0)
  {
    Vario * vario = createVario(&(params[pos]), nPar, params[pos-1], errText); 
    if(vario != NULL)
      modelSettings_->setLateralCorr(vario);
    else
      error = 1;
  }
  pos += nPar+1;
  return(error);
}


int 
ModelFile::readCommandSimulate(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  int nSimulations = atoi(params[pos]);
  modelSettings_->setNumberOfSimulations(nSimulations);
  pos += nPar+1;
  return(error);
}


int 
ModelFile::readCommandPadding(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 3);
  if(error == 0)
  {
    sprintf(errText,"%c",'\0');
    float xPad = float(atof(params[pos]));
    if(xPad > 1 || xPad < 0 || isNumber(params[pos]) == 0)
    {
      error = 1;
      sprintf(errText,"Padding in x-direction must be a number between 0 and 1 (found %s)\n",params[pos]);
    }
    else
    {
      modelSettings_->setXpad(xPad);
    }
    float yPad = float(atof(params[pos+1]));
    if(yPad > 1 || yPad <= 0 || isNumber(params[pos+1]) == 0)
    {
      error = 1;
      sprintf(errText,"%sPadding in y-direction must be a number between 0 and 1 (found %s)\n",
        errText,params[pos+1]);
    }
    else
    {
      modelSettings_->setYpad(yPad);
    }
    float zPad = float(atof(params[pos+2]));
    if(zPad > 1 || zPad <= 0 || isNumber(params[pos+2]) == 0)
    {
      error = 1;
      sprintf(errText,"%sPadding in z-direction must be a number between 0 and 1 (found %s)\n",
        errText,params[pos+2]);
    }
    else
    {
      modelSettings_->setZpad(zPad);
    }
  } 
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandSeed(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    if(isNumber(params[pos]))
      seed_ = atoi(params[pos]);
    else
    {
      error = checkFileOpen(&(params[pos]), 1, params[pos-1], errText);
      if(error == 0)
      {
        seedFile_ = new char[MAX_STRING];
        strcpy(seedFile_,params[pos]);
      }
    }
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandPrefix(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
    ModelSettings::setFilePrefix(params[pos]);
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandOutput(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1, -1);
  if(error == 0)
  {
    int i, key, nKeys = 23;
    char ** keywords;
    keywords = new char*[nKeys];
    for(i=0;i<nKeys;i++)
      keywords[i] = new char[40];
    i = 0;
    strcpy(keywords[i++],"STORM");
    strcpy(keywords[i++],"SEGY");
    strcpy(keywords[i++],"STORMASCII");
    strcpy(keywords[i++],"CORRELATION");
    strcpy(keywords[i++],"RESIDUALS");
    strcpy(keywords[i++],"VP");
    strcpy(keywords[i++],"VS");
    strcpy(keywords[i++],"RHO");
    strcpy(keywords[i++],"LAMELAMBDA");
    strcpy(keywords[i++],"LAMEMU");
    strcpy(keywords[i++],"POISSONRATIO");
    strcpy(keywords[i++],"AI");
    strcpy(keywords[i++],"SI");
    strcpy(keywords[i++],"VPVSRATIO");
    strcpy(keywords[i++],"MURHO");
    strcpy(keywords[i++],"LAMBDARHO");
    strcpy(keywords[i++],"PRIORCORRELATIONS");
    strcpy(keywords[i++],"BACKGROUND");
    strcpy(keywords[i++],"WELLS");
    strcpy(keywords[i++],"WAVELETS");
    strcpy(keywords[i++],"NOTIME");
    strcpy(keywords[i++],"FACIESPROB");
    strcpy(keywords[i++],"FACIESPROBRELATIVE");

    if (i != nKeys)
    { 
      printf("ERROR: i != nKeys\n");
      exit(1);
    }

    int outputFlag = 0;
    int formatFlag = 0;

    char * flag;
    sprintf(errText,"%c",'\0');

    for(i=0; i<nPar; i++)
    {
      flag = uppercase(params[pos+i]);
      for(key = 0;key < nKeys;key++)
        if(strcmp(flag,keywords[key]) == 0)
          break;
      if(key < 3)   
        formatFlag = (formatFlag | int(pow(2.0f,key)));
      else if(key < nKeys)
        outputFlag = (outputFlag | int(pow(2.0f,(key-2))));
      else
      {
        error = 1;
        sprintf(errText,"%sUnknown option '%s' in command %s. Valid options are\n", errText, flag,params[pos-1]);
        for (i = 0 ; i < nKeys ; i++)
          sprintf(errText,"%s  %s\n", errText, keywords[i]);
      }
    }
    for(i=0;i<nKeys;i++)
      delete [] keywords[i];
    delete [] keywords;

    if((outputFlag & ModelSettings::FACIESPROB) >0 && (outputFlag & ModelSettings::FACIESPROBRELATIVE)>0)
    {
      outputFlag -=1048576;
      LogKit::LogFormatted(LogKit::LOW,"Warning: Both FACIESPROB and FACIESPROBRELATIVE arr wanted as output. Only FACIESPROB is given.\n");
    }

    modelSettings_->setFormatFlag(formatFlag);
    modelSettings_->setOutputFlag(outputFlag);
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandWhiteNoise(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    float wnc = float(atof(params[pos]));
    if(wnc < 0.00005 || wnc > 0.1000001)
    {
      error = 1;
      sprintf(errText,"Error: Value given in command %s must be between 0.005 and 0.1\n",
        params[pos-1]);
    }
    else
      modelSettings_->setWNC(wnc);
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandSegYOffset(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
    modelSettings_->setSegyOffset(float(atof(params[pos])));
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandForceFile(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    int fileGrid = atoi(params[pos]);
    if(fileGrid < -1 || fileGrid > 1)
    {
      error = 1;
      sprintf(errText,"Error: Value given in command %s must be between -1 and 1\n",
        params[pos-1]);
    }
    modelSettings_->setFileGrid(fileGrid);
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandKriging(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1, -1);
  if(error == 0)
  {
    int i;
    float * krigingParams = new float[nPar];
    for(i=0;i<nPar;i++)
      krigingParams[i] = float(atof(params[pos+i]));
    modelSettings_->setKrigingParameters(krigingParams,nPar);
    delete [] krigingParams;
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandLocalWavelet(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 2, -1);
  if(error == 0)
  {
    int i, key, nKeys = 2;
    char ** keywords = new char*[nKeys];
    int * keyIndex = new int[nKeys];
    int * nArgs = new int[nKeys];
    for(i=0;i<nKeys;i++)
    {
      keywords[i] = new char[40];
      keyIndex[i] = -1;
    }
    i = 0;
    strcpy(keywords[i++],"SHIFT");
    strcpy(keywords[i++],"GAIN");

    char * param;
    sprintf(errText,"%c",'\0');
    int prevKey = -1;
    for(i=0; i<nPar; i++)
    {
      param = uppercase(params[pos+i]);
      for(key = 0;key < nKeys;key++)
        if(strcmp(param,keywords[key]) == 0)
          break;
      if(key < nKeys)
      {
        if(keyIndex[key] != -1)
        {
          error = 1;
          sprintf(errText,"%sError: Keyword %s given more than once in command %s.\n",
            errText, keywords[key], params[pos-1]);
        }
        keyIndex[key] = i;
        nArgs[key] = nPar-i-1;
        if(prevKey != -1)
          nArgs[prevKey] = keyIndex[key] - keyIndex[prevKey] - 1;
        prevKey = key;
      }
      else if(i == 0)
      {
        error = 1;
        sprintf(errText,"%sError: Unrecognized keyword %s given in command %s.\n",
          errText, keywords[key], params[pos-1]);
        break;
      }
    }
    if(error == 0) //Read keywords. Check number of arguments.
    {
      int j;
      for(i=0;i<nKeys;i++)
      {
        if(nArgs[i] == 0)
        {
          error = 1;
          sprintf(errText,"%sError: Keyword %s has 0 arguments in command %s.\n",
            errText, keywords[i], params[pos-1]);
        }
        for(j=i+1;j<nKeys;j++)
          if(nArgs[i] > 0 && nArgs[j] > 0 && nArgs[i] != nArgs[j])
          {
            error = 1;
            sprintf(errText,"%sError: Keywords %s and %s has different number of arguments (%d vs %d) in command %s.\n",
              errText, keywords[i], keywords[j], nArgs[i], nArgs[j], params[pos-1]);
          }
      }
    }
    if(error == 0) //Check file reading
    {
      for(i=0;i<nKeys;i++)
        if(keyIndex[i] > -1)
        {
          nWaveletTransfArgs_ = nArgs[i];
          int openError = checkFileOpen(&(params[pos+keyIndex[i]+1]), nArgs[i], params[pos-1], errText);
          if(openError != 0)
            error = 1;

          LogKit::LogFormatted(LogKit::LOW,"ERROR: Keyword LOCALWAVELET has temporarily been deactivated..\n");
          exit(1);

          /* 
          int j, tmpErr, flag = 1;
          irapgrid ** grids = new irapgrid * [nArgs[i]];
          for(j=0;j<nArgs[i];j++)
          {
            if((openError & flag) == 0 && params[pos+keyIndex[i]+1+j][0] !='*')
            {
              grids[j] = irapgridRead(params[pos+keyIndex[i]+1+j], &tmpErr);
              if(tmpErr != 0)
              {
                sprintf(errText, "Reading of %s grid from file %s failed in command %s.\n", 
                  keywords[i], params[pos+keyIndex[i]+1+j], params[pos-1]);
                error = 1;
              }
            }
            else
              grids[j] = NULL;
            flag *= 2;
          }
          switch(i) {
          case 0:
            shiftGrids_ = grids;
            break;
          case 1:
            gainGrids_ = grids;
          }
          */

        }
    }
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandEnergyTreshold(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    float energyTreshold = float(atof(params[pos]));
    if(energyTreshold < 0.0f || energyTreshold > 0.1000001)
    {
      error = 1;
      sprintf(errText,"Error: Value given in command %s must be between 0 and 0.1\n",
        params[pos-1]);
    }
    else
    {
      modelSettings_->setEnergyThreshold(energyTreshold);
    }
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandParameterCorr(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  paramCorrFile_ = new char[MAX_STRING]; 

  if(error == 0)
  {
    error = checkFileOpen(&(params[pos]), 1, params[pos-1], errText);
    if(error == 0)
      strcpy(paramCorrFile_,params[pos]);
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandReflectionMatrix(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
    error = checkFileOpen(&(params[pos]), 1, params[pos-1], errText);
  if(error == 0)
    strcpy(reflMatrFile_,params[pos]);
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandFrequencyBand(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 2, 2);
  if(error == 0)
  {
    float lowCut  = float(atof(params[pos]));
    float highCut = float(atof(params[pos+1]));
    if(lowCut >= highCut)
    {
      error = 1;
      sprintf(errText,"Error: High cut must be lager than low cut in command %s\n",params[pos-1]);
    }
    else
    {
      modelSettings_->setLowCut(lowCut);
      modelSettings_->setHighCut(highCut);
    }
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandDebug(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 0, 1);
  if(error == 0)
  {
    int debugFlag = 1;
    if(nPar > 0)
      debugFlag = int(atoi(params[pos]));
    modelSettings_->setDebugFlag(debugFlag);
  }
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandGiveSignalToNoiseRatios(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 0);
  if(error == 0)
  {
    hasSignalToNoiseRatio_ = true;
  }
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandSeismicResolution(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    float frequency = float(atof(params[pos]));
    if (frequency < 30.0f || frequency > 60.0f)
    {
      error = 1;
      sprintf(errText,"Error: The seismic resolution must be in the range 30Hz - 60Hz\n");
    }
    else
      modelSettings_->setMaxHzSeismic(frequency);
  }
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandWaveletTaperingL(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    modelSettings_->setWaveletTaperingL(float(atof(params[pos])));
  }
  pos += nPar+1;
  return(error);
}

int 
ModelFile::readCommandPUndef(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    modelSettings_->setPundef(float(atof(params[pos])));
  }
  pos += nPar+1;
  return(error);
}


int
ModelFile::readCommandMaxDeviationAngle(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 1);
  if(error == 0)
  {
    float angle = float(atof(params[pos]));
    if (angle <= 0.0f || angle > 90.0f)
    {
      error = 1;
      sprintf(errText,"Error: The deviation angle must be in the rang 0-90 degrees.\n");
    }
    else
      modelSettings_->setMaxDevAngle(angle);
  }
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandAllowedParameterValues(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 6);
  if(error == 0)
  {
    modelSettings_->setAlphaMin(float(atof(params[pos+0])));
    modelSettings_->setAlphaMax(float(atof(params[pos+1])));
    modelSettings_->setBetaMin (float(atof(params[pos+2])));
    modelSettings_->setBetaMax (float(atof(params[pos+3])));
    modelSettings_->setRhoMin  (float(atof(params[pos+4])));
    modelSettings_->setRhoMax  (float(atof(params[pos+5])));
  }
  pos += nPar+1;
  return(error);
}

int
ModelFile::readCommandAllowedParameterVariances(char ** params, int & pos, char * errText)
{
  int error;
  int nPar = getParNum(params, pos, error, errText, params[pos-1], 6);
  if(error == 0)
  {
    modelSettings_->setVarAlphaMin(float(atof(params[pos+0])));
    modelSettings_->setVarAlphaMax(float(atof(params[pos+1])));
    modelSettings_->setVarBetaMin (float(atof(params[pos+2])));
    modelSettings_->setVarBetaMax (float(atof(params[pos+3])));
    modelSettings_->setVarRhoMin  (float(atof(params[pos+4])));
    modelSettings_->setVarRhoMax  (float(atof(params[pos+5])));
  }
  pos += nPar+1;
  return(error);
}

int 
ModelFile::getParNum(char ** params, int pos, int & error, char * errText,
                 const char * command, int min, int max)
{
  int i=0;
  error = 0;
  while(params[pos+i][0] != ';')
    i++;
  if(max == 0) //min is exact limit)
  {
    if(i != min)
    {
      error = 1;
      sprintf(errText,"Command %s takes %d arguments; %d given in file.\n",
        command, min, i);
    }
  }
  else if(i < min)//May have no upper limit; does not matter.
  {
    error = 1;
    sprintf(errText,
      "Command %s takes at least %d arguments; %d given in file.\n",
      command, min, i);
  }
  else if(max > 0 && i > max)
  {
    error = 1;
    sprintf(errText,
      "Command %s takes at most %d arguments; %d given in file.\n",
      command, max, i);
  }
  return i;
}


// checkFileOpen: checks if the files in fNames table can be opened. Error message is
//                given in errText if not, and error is binary coded numbering of failures.  
int 
ModelFile::checkFileOpen(char ** fNames, int nFiles, const char * command, char * errText, int start,
                     bool details)
{
  if(details == true && nFiles >= 32)
  {
    details = false;
    LogKit::LogFormatted(LogKit::LOW,"\nINFO: Log details for checkFileOpen() has been turned off since nFiles>=32\n");
  }
  int i, error = 0;
  int flag = 1;
  int nErr = 0;
  char errFiles[MAX_STRING];
  strcpy(errFiles," ");
  for(i=start;i<nFiles+start;i++)
  {
    if(fNames[i][0] != '*' && fNames[i][0] != '?')
    {
      FILE * file = fopen(fNames[i],"r");
      if(file == 0)
      {
        error+= flag;
        nErr++;
        sprintf(errFiles,"%s%s ", errFiles, fNames[i]);
      }
      else
        fclose(file);
      if(details == true)
        flag *= 2;
    }
  }
  if(error > 0)
  {
    if(nErr == 1)
      sprintf(errText,"Could not open file %s (command %s).\n",errFiles, command);
    else
      sprintf(errText,"Could not open these files:%s (command %s).\n",errFiles, command);
  }
  return(error);
}

//NBNB The following routine is stupid, and assumes SEGY if file does not
//start with 'storm_petro_binary'

int
ModelFile::checkFileType(char * fileName)
{
  FILE * file = fopen(fileName,"r");
  char header[19];
  fread(header, 1, 18, file);
  fclose(file);
  header[18] = 0;
  if(strcmp(header,"storm_petro_binary") == 0)
    return(STORMFILE);
  else
    return(SEGYFILE);
}

Vario *
ModelFile::createVario(char ** param, int nPar, const char * command, char * errText)
{
  int i, nVario = 2; //New vario: Increase this.
  char * vTypeIn = uppercase(param[0]);
  char ** vTypes = new char *[nVario];
  int * varPar = new int[nVario];
  for(i=0;i<nVario;i++)
    vTypes[i] = new char[20];
  strcpy(vTypes[0],"SPHERICAL"); //Name of vario
  varPar[0] = 3;                 //Number of parameters for vario
  strcpy(vTypes[1],"GENEXP"); 
  varPar[1] = 4; //New vario: Add below.

  for(i=0;i<nVario;i++)
    if(strcmp(vTypeIn,vTypes[i])==0)
      break;

  Vario * result = NULL;
  if(i == nVario || nPar-1 == varPar[i] || nPar-1 == varPar[i] - 2)
    {
      switch(i)
        {
        case 0:
          if(nPar-1==varPar[i])
            result = new SphericalVario(float(atof(param[1])), float(atof(param[2])), 
                                        float(atof(param[3])));
          else
            result = new SphericalVario(float(atof(param[1])));
          break;
        case 1:
          if(nPar-1==varPar[i])
            result = new GenExpVario(float(atof(param[1])), float(atof(param[2])), 
                                     float(atof(param[3])), float(atof(param[4])));
          else
            result = new GenExpVario(float(atof(param[1])), float(atof(param[2])));
          break;
        default: //New vario: Add here.
          sprintf(errText, "Unknown variogram type %s in command %s.\n",
                  param[0], command);
          break;
        }
    }
  else
    sprintf(errText, 
            "Variogram %s in command %s needs %d parameters (%d given).\n",
            param[0], command, varPar[i], nPar-1);
  
  for(i=0;i<nVario;i++)
    delete [] vTypes[i];
  delete [] vTypes;
  delete [] varPar;
  return(result);
}

