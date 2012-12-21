#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "nrlib/exception/exception.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/flens/nrlib_flens.hpp"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/normal.hpp"

#include "src/modelsettings.h"
#include "src/xmlmodelfile.h"
#include "src/definitions.h"
#include "src/inputfiles.h"
#include "src/wavelet.h"
#include "src/fftgrid.h"
#include "src/vario.h"
#include "tasklist.h"
#include "src/io.h"

#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"

XmlModelFile::XmlModelFile(const std::string & fileName)
{
  modelSettings_ = new ModelSettings();
  inputFiles_    = new InputFiles();
  failed_        = false;
  surveyFailed_  = false;

  std::ifstream file;
  try {
    NRLib::OpenRead(file,fileName);
  }
  catch (NRLib::IOError e) {
    LogKit::LogMessage(LogKit::Error,"\nERROR: "+std::string(e.what()));
    exit(1);
  }

  //Remove all comments, since this convention is outside xml.
  std::string line;
  std::string active;
  std::string clean;
  while (std::getline(file,line))
  {
    active = line.substr(0, line.find_first_of("#"));
    clean = clean+active+"\n";
  }
  file.close();

  TiXmlDocument doc;
  doc.Parse(clean.c_str());

  if (doc.Error() == true) {
    LogKit::WriteHeader("Invalid XML file");
    LogKit::LogFormatted(LogKit::Error,"\n%s is not a valid XML file. %s In line %d, column %d.",
                         fileName.c_str(), doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
    if (doc.ErrorId() == 9) { // Not very robust check, but a start
      LogKit::LogFormatted(LogKit::Error,"\nPossible cause: Mis-spelled or forgotten end tag.");
    }
    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
    failed_ = true;
  }
  else {
    std::string errTxt = "";
    if(parseCrava(&doc, errTxt) == false)
      errTxt = "'"+std::string(fileName)+"' is not a crava model file (lacks the <crava> keyword.)\n";
    std::vector<std::string> legalCommands(1);
    checkForJunk(&doc, errTxt, legalCommands);
    checkConsistency(errTxt);

    if (errTxt == "") {
      setDerivedParameters(errTxt);
    }
    if (errTxt != "") {
      LogKit::WriteHeader("Invalid model file");
      LogKit::LogFormatted(LogKit::Error,"\n%s is not a valid model file:\n",fileName.c_str());
      LogKit::LogMessage(LogKit::Error, errTxt);
      LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
      failed_ = true;
    }
  }
}

XmlModelFile::~XmlModelFile()
{
  //Both modelSettings and inputFiles are taken out and deleted outside.
}

bool
XmlModelFile::parseCrava(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("crava");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("actions");
  legalCommands.push_back("project-settings");
  legalCommands.push_back("survey");
  legalCommands.push_back("well-data");
  legalCommands.push_back("prior-model");

  parseActions(root, errTxt);
  parseWellData(root, errTxt);

  size_t checkFailed = errTxt.size();
  while(parseSurvey(root, errTxt) == true);
  if(errTxt.size() > checkFailed)
    surveyFailed_ = true;

  parseProjectSettings(root, errTxt);
  parsePriorModel(root, errTxt);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseActions(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("actions");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("mode");
  legalCommands.push_back("inversion-settings");
  legalCommands.push_back("estimation-settings");

  std::string mode;
  if(parseValue(root, "mode", mode, errTxt) == false)
    errTxt += "Command <mode> must be given under command <"+
      root->ValueStr()+">"+lineColumnText(root)+".\n";
  else {
    if(mode == "forward")
      modelSettings_->setForwardModeling(true);
    else if(mode == "estimation") {
      modelSettings_->setEstimationMode(true);
      modelSettings_->setEstimateCorrelations(true);
      modelSettings_->setEstimateWaveletNoise(true);
      modelSettings_->setEstimateBackground(true);
      if (modelSettings_->getOutputGridsDefaultInd()) {
        modelSettings_->setOutputGridsElastic(IO::BACKGROUND);
        modelSettings_->setOutputGridsDefaultInd(false);
      }
      else {
        int paramFlag = modelSettings_->getOutputGridsElastic() & IO::BACKGROUND;
        modelSettings_->setOutputGridsElastic(paramFlag);
      }
    }
    else if(mode != "inversion")
      errTxt += "String after <mode> must be either <inversion>, <estimation> or <forward>, found <"+
        mode+"> under command <"+root->ValueStr()+">"+lineColumnText(root)+".\n";
  }

  if (parseInversionSettings(root, errTxt) && (mode == "forward" || mode == "estimation"))
    errTxt += "Inversion settings can only be given with the mode 'inversion'.\n";

  if (parseEstimationSettings(root, errTxt) && (mode == "forward" || mode == "inversion"))
    errTxt += "Estimation settings can only be given with the mode 'estimation'.\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseInversionSettings(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("inversion-settings");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("prediction");
  legalCommands.push_back("simulation");
  legalCommands.push_back("kriging-to-wells");
  legalCommands.push_back("facies-probabilities");

  bool value;
  if(parseBool(root, "prediction", value, errTxt) == true)
    modelSettings_->setWritePrediction(value);

  parseSimulation(root, errTxt);

  if(parseBool(root, "kriging-to-wells", value, errTxt) == true) {
    if(value == true) {
      if(modelSettings_->getKrigingParameter() == 0)
        modelSettings_->setKrigingParameter(250); //Default value.
    }
    else
      modelSettings_->setKrigingParameter(-1);
  }

  if(parseBool(root, "facies-probabilities", value, errTxt) == true && value == true) {
    modelSettings_->setFaciesProbRelative(true);
    modelSettings_->setEstimateFaciesProb(true);
    if (modelSettings_->getOutputGridsDefaultInd()) {
      modelSettings_->setOutputGridsOther(IO::FACIESPROB);
      modelSettings_->setOutputGridsDefaultInd(false);
      modelSettings_->setOutputGridsElastic(0);
    }
    else {
      int paramFlag = modelSettings_->getOutputGridsElastic() & IO::FACIESPROB;
      modelSettings_->setOutputGridsElastic(paramFlag);
    }
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseSimulation(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("simulation");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("seed");
  legalCommands.push_back("seed-file");
  legalCommands.push_back("number-of-simulations");

  int seed;
  bool seedGiven = parseValue(root, "seed", seed, errTxt);
  if(seedGiven == true)
    modelSettings_->setSeed(seed);

  std::string filename;
  if(parseFileName(root, "seed-file", filename, errTxt) == true) {
    if(seedGiven == true)
      errTxt += "Both seed and seed file given in command <"+
        root->ValueStr()+"> "+lineColumnText(root)+".\n";
    inputFiles_->setSeedFile(filename);
  }

  int value;
  bool numberGiven = parseValue(root, "number-of-simulations", value, errTxt);
  if(numberGiven == true)
    modelSettings_->setNumberOfSimulations(value);
  else
    modelSettings_->setNumberOfSimulations(1);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseEstimationSettings(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("estimation-settings");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("estimate-background");
  legalCommands.push_back("estimate-correlations");
  legalCommands.push_back("estimate-wavelet-or-noise");

  bool estimateBG = true;
  parseBool(root,"estimate-background", estimateBG, errTxt);
  modelSettings_->setEstimateBackground(estimateBG);

  bool estimateCorr = true;
  parseBool(root,"estimate-correlations", estimateCorr, errTxt);
  modelSettings_->setEstimateCorrelations(estimateCorr);

  bool estimateWN = true;
  parseBool(root,"estimate-wavelet-or-noise", estimateWN, errTxt);
  modelSettings_->setEstimateWaveletNoise(estimateWN);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWellData(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("well-data");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("log-names");
  legalCommands.push_back("well");
  legalCommands.push_back("high-cut-seismic-resolution");
  legalCommands.push_back("allowed-parameter-values");
  legalCommands.push_back("maximum-deviation-angle");
  legalCommands.push_back("maximum-rank-correlation");
  legalCommands.push_back("maximum-merge-distance");
  legalCommands.push_back("maximum-offset");
  legalCommands.push_back("maximum-shift");
  legalCommands.push_back("well-move-data-interval");

  parseLogNames(root, errTxt);

  int nWells = 0;
  while(parseWell(root, errTxt) == true)
    nWells++;
  modelSettings_->setNumberOfWells(nWells);

  float value;
  if(parseValue(root, "high-cut-seismic-resolution", value, errTxt) == true)
    modelSettings_->setHighCut(value);

  parseAllowedParameterValues(root, errTxt);

  if(parseValue(root, "maximum-deviation-angle", value, errTxt) == true)
    modelSettings_->setMaxDevAngle(value);

  if(parseValue(root, "maximum-rank-correlation", value, errTxt) == true)
    modelSettings_->setMaxRankCorr(value);

  if(parseValue(root, "maximum-merge-distance", value, errTxt) == true)
    modelSettings_->setMaxMergeDist(value);

  float offset;
  if(parseValue(root, "maximum-offset", offset, errTxt) == true)
    modelSettings_->setMaxWellOffset(offset);

  if(parseValue(root, "maximum-shift", value, errTxt) == true)
    modelSettings_->setMaxWellShift(value);

  parseWellMoveDataInterval(root, errTxt);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseLogNames(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("log-names");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("time");
  legalCommands.push_back("vp");
  legalCommands.push_back("dt");
  legalCommands.push_back("vs");
  legalCommands.push_back("dts");
  legalCommands.push_back("density");
  legalCommands.push_back("facies");

  std::string value;
  if(parseValue(root, "time", value, errTxt) == true)
    modelSettings_->setLogName(0, value);

  bool vp = parseValue(root, "vp", value, errTxt);
  if(vp == true) {
    modelSettings_->setLogName(1, value);
    modelSettings_->setInverseVelocity(0, false);
  }
  if(parseValue(root, "dt", value, errTxt) == true) {
    if(vp == true)
      errTxt += "Both vp and dt given as logs in command <"
        +root->ValueStr()+"> "+lineColumnText(root)+".\n";
    else {
      modelSettings_->setLogName(1, value);
      modelSettings_->setInverseVelocity(0, true);
    }
  }

  bool vs = parseValue(root, "vs", value, errTxt);
  if(vp == true) {
    modelSettings_->setLogName(3, value);
    modelSettings_->setInverseVelocity(1, false);
  }
  if(parseValue(root, "dts", value, errTxt) == true) {
    if(vs == true)
      errTxt += "Both vs and dts given as logs in command <"
        +root->ValueStr()+"> "+lineColumnText(root)+".\n";
    else {
      modelSettings_->setLogName(3, value);
      modelSettings_->setInverseVelocity(1, true);
    }
  }

  if(parseValue(root, "density", value, errTxt) == true)
    modelSettings_->setLogName(2, value);
  if(parseValue(root, "facies", value, errTxt) == true) {
    modelSettings_->setLogName(4, value);
    modelSettings_->setFaciesLogGiven(true);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWell(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("well");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("use-for-wavelet-estimation");
  legalCommands.push_back("use-for-background-trend");
  legalCommands.push_back("use-for-facies-probabilities");
  legalCommands.push_back("synthetic-vs-log");
  legalCommands.push_back("optimize-position");
  legalCommands.push_back("filter-elastic-logs");

  std::string tmpErr = "";
  std::string value;
  if(parseValue(root, "file-name", value, tmpErr) == true) {
    inputFiles_->addWellFile(value);
  }
  else {
    inputFiles_->addWellFile(""); //Dummy to keep tables balanced.
    tmpErr += "The file name for the well is not given.\n";
  }

  int useForBackgroundTrend   = ModelSettings::NOTSET;
  int useForWaveletEstimation = ModelSettings::NOTSET;
  int useForFaciesProbability = ModelSettings::NOTSET;
  int filterElasticLogs       = ModelSettings::NOTSET;
  int hasRealVs               = ModelSettings::NOTSET;

  bool use;
  if(parseBool(root, "use-for-background-trend", use, tmpErr)) {
    if(use)
      useForBackgroundTrend = ModelSettings::YES;
    else
      useForBackgroundTrend = ModelSettings::NO;
  }

  if(parseBool(root, "use-for-wavelet-estimation", use, tmpErr)) {
    if(use)
      useForWaveletEstimation = ModelSettings::YES;
    else
      useForWaveletEstimation = ModelSettings::NO;
  }

  if(parseBool(root, "use-for-facies-probabilities", use, tmpErr)) {
    if(use)
      useForFaciesProbability = ModelSettings::YES;
    else
      useForFaciesProbability = ModelSettings::NO;
  }

  bool synth = false;
  if(parseBool(root, "synthetic-vs-log", synth, tmpErr)) {
    if(synth == false)
      hasRealVs = ModelSettings::YES;   //Note the inversion. The keyword is most logical this way,
    else                                //but the default should be that the log is real.
      hasRealVs = ModelSettings::NO;
  }

  if(parseBool(root, "filter-elastic-logs", use, tmpErr)) {
    if(use)
      filterElasticLogs = ModelSettings::YES;
    else
      filterElasticLogs = ModelSettings::NO;
  }

  modelSettings_->addIndicatorBGTrend(useForBackgroundTrend);
  modelSettings_->addIndicatorWavelet(useForWaveletEstimation);
  modelSettings_->addIndicatorFacies(useForFaciesProbability);
  modelSettings_->addIndicatorRealVs(hasRealVs);
  modelSettings_->addIndicatorFilter(filterElasticLogs);

  modelSettings_->clearMoveWell();
  while(parseOptimizeLocation(root, tmpErr) == true);
  modelSettings_->addMoveWell();

  checkForJunk(root, errTxt, legalCommands, true); //Allow duplicates

  errTxt += tmpErr;
  return(true);
}

bool
XmlModelFile::parseOptimizeLocation(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("optimize-position");
  if(root == 0)
    return(false);

  modelSettings_->setOptimizeWellLocation(true);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("angle");
  legalCommands.push_back("weight");

  float value;
  if(parseValue(root, "angle", value, errTxt) == true)
    modelSettings_->addMoveAngle(value*float(NRLib::Pi/180));
  else
    modelSettings_->addMoveAngle(RMISSING);

  if(parseValue(root, "weight", value, errTxt) == true)
    modelSettings_->addMoveWeight(value);
  else
    modelSettings_->addMoveWeight(1); // Default

  checkForJunk(root, errTxt, legalCommands, true);

  return(true);
}


bool
XmlModelFile::parseAllowedParameterValues(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("allowed-parameter-values");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("minimum-vp");
  legalCommands.push_back("maximum-vp");
  legalCommands.push_back("minimum-vs");
  legalCommands.push_back("maximum-vs");
  legalCommands.push_back("minimum-density");
  legalCommands.push_back("maximum-density");
  legalCommands.push_back("minimum-variance-vp");
  legalCommands.push_back("maximum-variance-vp");
  legalCommands.push_back("minimum-variance-vs");
  legalCommands.push_back("maximum-variance-vs");
  legalCommands.push_back("minimum-variance-density");
  legalCommands.push_back("maximum-variance-density");
  legalCommands.push_back("minimum-vp-vs-ratio");
  legalCommands.push_back("maximum-vp-vs-ratio");

  float value;
  if(parseValue(root, "minimum-vp", value, errTxt) == true)
    modelSettings_->setAlphaMin(value);
  if(parseValue(root, "maximum-vp", value, errTxt) == true)
    modelSettings_->setAlphaMax(value);
  if(parseValue(root, "minimum-vs", value, errTxt) == true)
    modelSettings_->setBetaMin(value);
  if(parseValue(root, "maximum-vs", value, errTxt) == true)
    modelSettings_->setBetaMax(value);
  if(parseValue(root, "minimum-density", value, errTxt) == true)
    modelSettings_->setRhoMin(value);
  if(parseValue(root, "maximum-density", value, errTxt) == true)
    modelSettings_->setRhoMax(value);

  if(parseValue(root, "minimum-variance-vp", value, errTxt) == true)
    modelSettings_->setVarAlphaMin(value);
  if(parseValue(root, "maximum-variance-vp", value, errTxt) == true)
    modelSettings_->setVarAlphaMax(value);
  if(parseValue(root, "minimum-variance-vs", value, errTxt) == true)
    modelSettings_->setVarBetaMin(value);
  if(parseValue(root, "maximum-variance-vs", value, errTxt) == true)
    modelSettings_->setVarBetaMax(value);
  if(parseValue(root, "minimum-variance-density", value, errTxt) == true)
    modelSettings_->setVarRhoMin(value);
  if(parseValue(root, "maximum-variance-density", value, errTxt) == true)
    modelSettings_->setVarRhoMax(value);

  if(parseValue(root, "minimum-vp-vs-ratio", value, errTxt) == true)
    modelSettings_->setVpVsRatioMin(value);
  if(parseValue(root, "maximum-vp-vs-ratio", value, errTxt) == true)
    modelSettings_->setVpVsRatioMax(value);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseWellMoveDataInterval(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("well-move-data-interval");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-surface-file");
  legalCommands.push_back("base-surface-file");

  std::string filename;
  if(parseFileName(root, "top-surface-file", filename, errTxt) == true)
    inputFiles_->setWellMoveIntFile(0, filename);
  if(parseFileName(root, "base-surface-file", filename, errTxt) == true)
    inputFiles_->setWellMoveIntFile(1, filename);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseSurvey(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("survey");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("angular-correlation");
  legalCommands.push_back("segy-start-time");
  legalCommands.push_back("angle-gather");
  legalCommands.push_back("wavelet-estimation-interval");
  legalCommands.push_back("time-gradient-settings");
  legalCommands.push_back("vintage");

  Vario * vario = NULL;
  modelSettings_->addDefaultAngularCorr();
  if(parseVariogram(root, "angular-correlation", vario, errTxt) == true) {
    if (vario != NULL) {
      vario->convertRangesFromDegToRad();
      modelSettings_->setLastAngularCorr(vario);
    }
  }

  float value;
  if(parseValue(root, "segy-start-time", value, errTxt) == true)
    modelSettings_->addSegyOffset(value);
  else
    modelSettings_->addDefaultSegyOffset();

  modelSettings_->clearTimeLapse();
  inputFiles_->clearTimeLapse();
  while(parseAngleGather(root, errTxt) == true);
  modelSettings_->addTimeLapse();
  inputFiles_->addTimeLapse();
  if(modelSettings_->getNumberOfTimeLapses() > 1)
    modelSettings_->setDo4DInversion(true);

  if(parseWaveletEstimationInterval(root, errTxt) == false){
    inputFiles_->addDefaultWaveletEstIntFileTop();
    inputFiles_->addDefaultWaveletEstIntFileBase();
  }

  if(parseTimeGradientSettings(root,errTxt) == false)
    modelSettings_->addDefaultTimeGradientSettings();

  if(parseVintage(root, errTxt) == false)
    modelSettings_->addDefaultVintage();

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}


bool
XmlModelFile::parseAngleGather(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("angle-gather");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("offset-angle");
  legalCommands.push_back("seismic-data");
  legalCommands.push_back("wavelet");
  legalCommands.push_back("wavelet-3d");
  legalCommands.push_back("match-energies");
  legalCommands.push_back("signal-to-noise-ratio");
  legalCommands.push_back("local-noise-scaled");
  legalCommands.push_back("estimate-local-noise");

  float value;
  if(parseValue(root, "offset-angle", value, errTxt) == true)
    modelSettings_->addAngle(value*float(NRLib::Pi/180));
  else
    errTxt += "Need offset angle for gather"+lineColumnText(root)+".\n";

  if(parseSeismicData(root, errTxt) == false) {
    //Go for defaults. Assume forward model (check later)
    inputFiles_->addSeismicFile("");
    modelSettings_->addSeismicType(ModelSettings::STANDARDSEIS);
  }

  bool oneDwavelet(false);
  if(parseWavelet(root, errTxt) == false) {
    modelSettings_->addWaveletScale(1.0); // NBNB OK  why RMISSING here????
    modelSettings_->addEstimateGlobalWaveletScale(false);
    modelSettings_->addUseRickerWavelet(false);
    modelSettings_->addRickerPeakFrequency(RMISSING);
    inputFiles_->addShiftFile("");
    inputFiles_->addScaleFile("");
    modelSettings_->addEstimateLocalShift(false);
    modelSettings_->addEstimateLocalScale(false);
  }
  else //1D-wavelet is given
    oneDwavelet = true;

  if(parseWavelet3D(root, errTxt) == false) {
    inputFiles_->addWaveletFilterFile("");
    modelSettings_->addWaveletDim(Wavelet::ONE_D);
    if (oneDwavelet == false) { //Neither 1D-wavelet nor 3D-wavelet are given
      inputFiles_->addWaveletFile("");
      modelSettings_->addEstimateWavelet(true);
      modelSettings_->addStretchFactor(1.0);
    }
  }
  else //3D-wavelet is given
    if (oneDwavelet) //1D-wavelet is also given. Arrays might be unbalanced
      errTxt += "Both wavelet and wavelet-3d given in command <"+
        root->ValueStr()+"> "+lineColumnText(root)+".\n";

  bool bVal = false;
  if(parseBool(root, "match-energies", bVal, errTxt) == true)
    modelSettings_->addMatchEnergies(bVal);
  else
    modelSettings_->addMatchEnergies(false);

  if(parseValue(root, "signal-to-noise-ratio", value, errTxt) == true) {
    modelSettings_->addEstimateSNRatio(false);
    modelSettings_->addSNRatio(value);
  }
  else {
    modelSettings_->addEstimateSNRatio(true);
    modelSettings_->addSNRatio(RMISSING);
  }
  std::string fileName;
  bool localNoiseGiven = parseFileName(root, "local-noise-scaled", fileName, errTxt);
  modelSettings_->setDefaultUseLocalNoise();
  if(localNoiseGiven) {
    inputFiles_->addNoiseFile(fileName);
    modelSettings_->setUseLocalNoise(true);
  }
  else {
    inputFiles_->addNoiseFile("");
  }
  bool estimate;
  std::string tmpErr;
  if(parseBool(root, "estimate-local-noise",estimate,tmpErr) == false || tmpErr != "") {
    modelSettings_->addEstimateLocalNoise(false);
    errTxt += tmpErr;
  }
  else {
    if(estimate == true && localNoiseGiven == true) {
      errTxt = errTxt +"Error: Can not give both file and ask for estimation of local noise"+
        lineColumnText(node);
      modelSettings_->addEstimateLocalNoise(false);
    }
    else if (estimate==false && localNoiseGiven==false)
      modelSettings_->addEstimateLocalNoise(false);
    else {
      modelSettings_->setUseLocalNoise(true);
      modelSettings_->addEstimateLocalNoise(estimate);
    }
  }
  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}


bool
XmlModelFile::parseSeismicData(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("seismic-data");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("start-time");
  legalCommands.push_back("segy-format");
  legalCommands.push_back("type");


  std::string value;
  if(parseFileName(root, "file-name", value, errTxt) == true)
    inputFiles_->addSeismicFile(value);
  else
  {
    inputFiles_->addSeismicFile(""); //Keeping tables balanced.
    errTxt += "The file name for the seismic data is not given.\n";
  }

  if(parseValue(root, "type", value, errTxt) == true) {
    if(value == "pp")
      modelSettings_->addSeismicType(ModelSettings::STANDARDSEIS);
    else if(value == "ps")
      modelSettings_->addSeismicType(ModelSettings::PSSEIS);
    else
      errTxt += "Only 'pp' and 'ps' are valid seismic types, found '"+value+"' on line "
        +NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
  }
  else
    modelSettings_->addSeismicType(ModelSettings::STANDARDSEIS);

  float fVal;
  if(parseValue(root, "start-time", fVal, errTxt) == true)
    modelSettings_->addLocalSegyOffset(fVal);
  else
    modelSettings_->addLocalSegyOffset(-1);

  TraceHeaderFormat * thf = NULL;
  if(parseTraceHeaderFormat(root, "segy-format", thf, errTxt) == true) {
    modelSettings_->addTraceHeaderFormat(thf);
  }
  else
    modelSettings_->addTraceHeaderFormat(NULL);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWavelet(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("wavelet");
  if(root == 0)
  {
    return(false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("ricker");
  legalCommands.push_back("scale");
  legalCommands.push_back("estimate-scale");
  legalCommands.push_back("local-wavelet");

  modelSettings_->addStretchFactor(1.0);
  std::string value;
  bool estimate = false;
  float peakFrequency;
  if(parseFileName(root, "file-name", value, errTxt) == true) {
    inputFiles_->addWaveletFile(value);
    modelSettings_->addEstimateWavelet(false);
    modelSettings_->addUseRickerWavelet(false);
    modelSettings_->addRickerPeakFrequency(RMISSING);
  }
  else if(parseValue(root, "ricker", peakFrequency, errTxt) == true){
    inputFiles_->addWaveletFile(""); //Keeping tables balanced.
    modelSettings_->addRickerPeakFrequency(peakFrequency);
    modelSettings_->addEstimateWavelet(false);
    modelSettings_->addUseRickerWavelet(true);
  }
  else {
    inputFiles_->addWaveletFile(""); //Keeping tables balanced.
    modelSettings_->addEstimateWavelet(true);
    modelSettings_->addUseRickerWavelet(false);
    modelSettings_->addRickerPeakFrequency(RMISSING);
    estimate = true;
  }

  float scale;
  bool scaleGiven = false;
  if(parseValue(root, "scale", scale, errTxt) == true)
  {
    scaleGiven = true;
    modelSettings_->addWaveletScale(scale);
    modelSettings_->addEstimateGlobalWaveletScale(false);
    if (estimate)
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: The wavelet scale specified in the model file ("
           +NRLib::ToString(scale,2)
           +") has no effect when the wavelet\n         is estimated and not read from file\n\n");
  }
  std::string tmpErr;

  if(parseBool(root, "estimate-scale",estimate,tmpErr) == false || tmpErr != "") {
   if(scaleGiven==false) // no commands given
   {
    modelSettings_->addWaveletScale(1);
    modelSettings_->addEstimateGlobalWaveletScale(false);
   }
    errTxt += tmpErr;
  }
  else
    if(estimate == true && scaleGiven == true) {
      errTxt += "Error: Can not give both value and ask for estimation of global wavelet scale"+
        lineColumnText(node);
     // modelSettings_->addEstimateGlobalWaveletScale(false);
    }
    else if(estimate==true && scaleGiven==false)// estimate==true eller false, scaleGiven==false
    {
      modelSettings_->addEstimateGlobalWaveletScale(true);
      modelSettings_->addWaveletScale(1);
    }
    else if(estimate==false && scaleGiven==false)// estimate given as no, no scale command given
    {
      modelSettings_->addWaveletScale(1);
      modelSettings_->addEstimateGlobalWaveletScale(false);
    }


  if(parseLocalWavelet(root, errTxt)==false)
  {
    inputFiles_->addShiftFile("");
    inputFiles_->addScaleFile("");
    modelSettings_->addEstimateLocalShift(false);
    modelSettings_->addEstimateLocalScale(false);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseLocalWavelet(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("local-wavelet");
  if(root == 0)
  {
    return(false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("shift-file");
  legalCommands.push_back("scale-file");
  legalCommands.push_back("estimate-shift");
  legalCommands.push_back("estimate-scale");

  modelSettings_->setUseLocalWavelet(false);
  std::string filename;
  bool shiftGiven = parseFileName(root, "shift-file", filename, errTxt);
  if(shiftGiven)
  {
    inputFiles_->addShiftFile(filename);
    modelSettings_->setUseLocalWavelet(true);
  }
  else
    inputFiles_->addShiftFile("");

  bool scaleGiven = parseFileName(root, "scale-file", filename, errTxt);
  if(scaleGiven)
  {
    inputFiles_->addScaleFile(filename);
    modelSettings_->setUseLocalWavelet(true);
  }
  else
    inputFiles_->addScaleFile("");

  bool estimate;
  std::string tmpErr;
  if(parseBool(root, "estimate-shift",estimate,tmpErr) == false || tmpErr != "") {
    modelSettings_->addEstimateLocalShift(false);
    errTxt += tmpErr;
  }
  else
    if(estimate == true && shiftGiven == true) {
      errTxt += "Can not give both file and ask for estimation of local wavelet shift"+
        lineColumnText(node);
      modelSettings_->addEstimateLocalShift(false);
    }
    else if (estimate == true)
    {
      modelSettings_->addEstimateLocalShift(true);
      modelSettings_->setUseLocalWavelet(true);
    }
    else
      modelSettings_->addEstimateLocalShift(false);

  tmpErr = "";
  if(parseBool(root, "estimate-scale",estimate,tmpErr) == false || tmpErr != "") {
    modelSettings_->addEstimateLocalScale(false);
    errTxt += tmpErr;
  }
  else
    if(estimate == true && scaleGiven == true) {
      errTxt += "Can not give both file and ask for estimation of local wavelet scale"+
        lineColumnText(node);
      modelSettings_->addEstimateLocalScale(false);
    }
    else if (estimate == true)
    {
      modelSettings_->addEstimateLocalScale(true);
      modelSettings_->setUseLocalWavelet(true);
    }
    else
      modelSettings_->addEstimateLocalScale(false);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWaveletEstimationInterval(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("wavelet-estimation-interval");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-surface-file");
  legalCommands.push_back("base-surface-file");

  std::string filenameTop  = "";
  std::string filenameBase = "";
  parseFileName(root, "top-surface-file", filenameTop, errTxt);
  parseFileName(root, "base-surface-file", filenameBase, errTxt);

  inputFiles_->addWaveletEstIntFileTop(filenameTop);
  inputFiles_->addWaveletEstIntFileBase(filenameBase);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWavelet3D(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("wavelet-3d");
  if(root == 0)
  {
    return(false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("processing-factor-file-name");
  legalCommands.push_back("propagation-factor-file-name");
  legalCommands.push_back("stretch-factor");
  legalCommands.push_back("estimation-range-x-direction");
  legalCommands.push_back("estimation-range-y-direction");

  std::string value;
  if(parseFileName(root, "file-name", value, errTxt) == true) {
    inputFiles_->addWaveletFile(value);
    modelSettings_->addEstimateWavelet(false);
  }
  else {
    inputFiles_->addWaveletFile(""); //Keeping tables balanced.
    modelSettings_->addEstimateWavelet(true);
    modelSettings_->setEstimate3DWavelet(true);
  }

  if(parseFileName(root, "processing-factor-file-name", value, errTxt) == true)
    inputFiles_->addWaveletFilterFile(value);
  else {
    inputFiles_->addWaveletFilterFile(""); //Keeping tables balanced.
    errTxt += "No processing factor file for 3D-wavelet is given in<"
              +root->ValueStr()+"> "+ lineColumnText(root)+".\n";
  }

  if(parseFileName(root, "propagation-factor-file-name", value, errTxt) == true)
    inputFiles_->addWaveletCorrFile(value);
  else
    inputFiles_->addWaveletCorrFile(""); //Keeping tables balanced.

  float stretch;
  if(parseValue(root, "stretch-factor", stretch, errTxt) == true)
    modelSettings_->addStretchFactor(stretch);
  else
    modelSettings_->addStretchFactor(1.0);

  float range;
  if(parseValue(root, "estimation-range-x-direction", range, errTxt) == true)
    modelSettings_->addEstRangeX(range);
  else
    modelSettings_->addEstRangeX(0.0);

  if(parseValue(root, "estimation-range-y-direction", range, errTxt) == true)
    modelSettings_->addEstRangeY(range);
  else
    modelSettings_->addEstRangeY(0.0);

  modelSettings_->addWaveletDim(Wavelet::THREE_D);
  modelSettings_->setUse3DWavelet(true);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}
bool
XmlModelFile::parseTimeGradientSettings(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("time-gradient-settings");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("distance");
  legalCommands.push_back("sigma");

  float distance, sigma_m, value;
  if(parseValue(root,"distance", value, errTxt) == true)
    distance = value;
  else
    distance = 100.0;
  if(parseValue(root,"sigma", value, errTxt) == true)
    sigma_m = value;
  else
    sigma_m = 1.0;

  modelSettings_->addTimeGradientSettings(distance,sigma_m);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseVintage(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("vintage");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("year");
  legalCommands.push_back("month");
  legalCommands.push_back("day-of-month");

  int value;
  int year  = IMISSING;
  int month = IMISSING;
  int day   = IMISSING;

  if(parseValue(root, "year", value, errTxt) == true)
    year = value;

  if(parseValue(root, "month", value, errTxt) == true){
    if(year == IMISSING)
      errTxt += "<year> needs to be specified if <month> is specified in <vintage>.\n";
    else
      month = value;
  }

  if(parseValue(root, "day-of-month", value, errTxt) == true){
    if(year == IMISSING || month == IMISSING)
      errTxt += "Both <year> and <month> need to be specified if <day-of-month> is specifeid in <vintage>.\n";
    else
      day = value;
  }

  modelSettings_->addVintage(year, month, day);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parsePriorModel(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("prior-model");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("background");
  legalCommands.push_back("lateral-correlation");
  legalCommands.push_back("temporal-correlation");
  legalCommands.push_back("parameter-correlation");
  legalCommands.push_back("correlation-direction");
  legalCommands.push_back("facies-probabilities");
  legalCommands.push_back("earth-model");
  legalCommands.push_back("local-wavelet");
  legalCommands.push_back("rock-physics");

  parseBackground(root, errTxt);
  parseEarthModel(root,errTxt);
  parsePriorLocalWavelet(root, errTxt);

  Vario* vario = NULL;
  if(parseVariogram(root, "lateral-correlation", vario, errTxt) == true) {
    if (vario != NULL) {
      modelSettings_->setLateralCorr(vario);
    }
  }

  std::string filename;
  if(parseFileName(root, "temporal-correlation", filename, errTxt) == true)
    inputFiles_->setTempCorrFile(filename);

  if(parseFileName(root, "parameter-correlation", filename, errTxt) == true)
    inputFiles_->setParamCorrFile(filename);

  if(parseFileName(root, "correlation-direction", filename, errTxt) == true)
    inputFiles_->setCorrDirFile(filename);

  parseFaciesProbabilities(root, errTxt);
  parseRockPhysics(root, errTxt);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parsePriorLocalWavelet(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("local-wavelet");
  if(root == 0)
    return(false);
  std::vector<std::string> legalCommands;
  legalCommands.push_back("lateral-correlation");
  Vario* vario = NULL;
  if(parseVariogram(root, "lateral-correlation", vario, errTxt) == true) {
    if (vario != NULL) {
      modelSettings_->setLocalWaveletVario(vario);
    }
  }
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseEarthModel(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("earth-model");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("vp-file");
  legalCommands.push_back("vs-file");
  legalCommands.push_back("density-file");
  legalCommands.push_back("ai-file");
  legalCommands.push_back("si-file");
  legalCommands.push_back("vp-vs-ratio-file");

  if (modelSettings_->getBackgroundType() != "")
    errTxt += "Both background and earth-model can not be given. Under forward mode, the earth-model must be used.\n";
  modelSettings_->setBackgroundType("earth model");

  std::string filename;
  bool vp = parseFileName(root, "vp-file", filename, errTxt);
  if(vp == true) {
    inputFiles_->setBackFile(0, filename);
    modelSettings_->setConstBackValue(0, -1);
  }

  bool vs = parseFileName(root, "vs-file", filename, errTxt);
  if(vs == true) {
    inputFiles_->setBackFile(1, filename);
    modelSettings_->setConstBackValue(1, -1);
  }

  bool rho = parseFileName(root, "density-file", filename, errTxt);
  if(rho == true) {
    inputFiles_->setBackFile(2, filename);
    modelSettings_->setConstBackValue(2, -1);
  }

  bool ai = parseFileName(root, "ai-file", filename, errTxt);
  if(ai == true) {
    inputFiles_->setBackFile(0, filename);     // Store AI earth model in Vp slot
    modelSettings_->setConstBackValue(0, -1);  //
    modelSettings_->setUseAIBackground(true);
  }

  bool si = parseFileName(root, "si-file", filename, errTxt);
  if(si == true) {
    inputFiles_->setBackFile(1, filename);     // Store SI earth model in Vs slot
    modelSettings_->setConstBackValue(1, -1);  //
    modelSettings_->setUseSIBackground(true);
  }

  bool vpvs = parseFileName(root, "vp-vs-ratio-file", filename, errTxt);
  if(vpvs == true) {
    inputFiles_->setBackFile(1, filename);     // Store Vp/Vs earth model in Vs slot
    modelSettings_->setConstBackValue(1, -1);  //
    modelSettings_->setUseVpVsBackground(true);
  }

  if (vp && ai) {
      errTxt += "Earth model for both AI and Vp has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  if (vs && si) {
      errTxt += "Earth model for both SI and Vs has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  if (vs && vpvs) {
      errTxt += "Earth model for both Vs and Vp/Vs has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  if (si && vpvs) {
      errTxt += "Earth model for both SI and Vp/Vs has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  bool earthModelGiven =
    (vp & vs   & rho)  ||
    (vp & si   & rho)  ||
    (vp & vpvs & rho)  ||
    (ai & vs   & rho)  ||
    (ai & si   & rho)  ||
    (ai & vpvs & rho);

  if (!earthModelGiven) {
      errTxt += "Earth model has not been completely specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please specify {Vp,Vs,Rho}, {AI,SI,Rho}, or {AI,Vp/Vs,Rho}.\n";
  }

  modelSettings_->setGenerateBackground(false);
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseBackground(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("background");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("vp-file");
  legalCommands.push_back("vs-file");
  legalCommands.push_back("density-file");
  legalCommands.push_back("ai-file");
  legalCommands.push_back("si-file");
  legalCommands.push_back("vp-vs-ratio-file");
  legalCommands.push_back("vp-constant");
  legalCommands.push_back("vs-constant");
  legalCommands.push_back("density-constant");
  legalCommands.push_back("velocity-field");
  legalCommands.push_back("lateral-correlation");
  legalCommands.push_back("high-cut-background-modelling");
  legalCommands.push_back("multizone-model");

  modelSettings_->setBackgroundType("background");

  std::string filename;
  bool vp = parseFileName(root, "vp-file", filename, errTxt);
  if(vp == true) {
    inputFiles_->setBackFile(0, filename);
    modelSettings_->setConstBackValue(0, -1);
  }

  bool vs = parseFileName(root, "vs-file", filename, errTxt);
  if(vs == true) {
    inputFiles_->setBackFile(1, filename);
    modelSettings_->setConstBackValue(1, -1);
  }

  bool rho = parseFileName(root, "density-file", filename, errTxt);
  if(rho == true) {
    inputFiles_->setBackFile(2, filename);
    modelSettings_->setConstBackValue(2, -1);
  }

  bool ai = parseFileName(root, "ai-file", filename, errTxt);
  if(ai == true) {
    inputFiles_->setBackFile(0, filename);     // Store AI background in Vp slot
    modelSettings_->setConstBackValue(0, -1);  //
    modelSettings_->setUseAIBackground(true);
  }

  bool si = parseFileName(root, "si-file", filename, errTxt);
  if(si == true) {
    inputFiles_->setBackFile(1, filename);     // Store SI background in Vs slot
    modelSettings_->setConstBackValue(1, -1);  //
    modelSettings_->setUseSIBackground(true);
  }

  bool vpvs = parseFileName(root, "vp-vs-ratio-file", filename, errTxt);
  if(vpvs == true) {
    inputFiles_->setBackFile(1, filename);     // Store Vp/Vs background in Vs slot
    modelSettings_->setConstBackValue(1, -1);  //
    modelSettings_->setUseVpVsBackground(true);
  }

  float value;
  if(parseValue(root, "vp-constant", value, errTxt) == true) {
    if(vp == true)
      errTxt += "Both file and constant given for vp in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+".\n";
    else {
      vp = true;
      modelSettings_->setConstBackValue(0, value);
    }
  }

  if(parseValue(root, "vs-constant", value, errTxt) == true) {
    if(vs == true)
      errTxt += "Both file and constant given for vs in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+".\n";
    else {
      vs = true;
      modelSettings_->setConstBackValue(1, value);
    }
  }

  if(parseValue(root, "density-constant", value, errTxt) == true) {
    if(rho == true)
      errTxt += "Both file and constant given for density in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+".\n";
    else {
      rho = true;
      modelSettings_->setConstBackValue(2, value);
    }
  }

  if (vp && ai) {
      errTxt += "Both AI background and Vp background has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  if (vs && si) {
      errTxt += "Both SI background and Vs background has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  if (vs && vpvs) {
      errTxt += "Both Vs background and Vp/Vs background has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  if (si && vpvs) {
      errTxt += "Both SI background and Vp/Vs background has been specified in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+". Please give only one.\n";
  }

  bool bgGiven =
    (vp & vs   & rho)  ||
    (vp & si   & rho)  ||
    (vp & vpvs & rho)  ||
    (ai & vs   & rho)  ||
    (ai & si   & rho)  ||
    (ai & vpvs & rho);

  bool estimate = !(vp | vs | rho);
  modelSettings_->setGenerateBackground(estimate);
  if((bgGiven | estimate) == false) {
    errTxt +=  "Either all or no background parameters must be given in command <"+root->ValueStr()+"> "+
        lineColumnText(root)+".\n";
  }
  if(parseFileName(root, "velocity-field", filename, errTxt) == true)
    inputFiles_->setBackVelFile(filename);

  Vario * vario = NULL;
  if(parseVariogram(root, "lateral-correlation", vario, errTxt) == true) {
    if (vario != NULL) {
      modelSettings_->setBackgroundVario(vario);
    }
  }

  if(parseValue(root, "high-cut-background-modelling", value, errTxt) == true)
    modelSettings_->setMaxHzBackground(value);

  parseMultizoneModel(root,errTxt);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseMultizoneModel(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("multizone-model");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-surface-file");
  legalCommands.push_back("top-surface-erosion-priority");
  legalCommands.push_back("zone");

  std::string file_name;
  if(parseFileName(root, "top-surface-file", file_name, errTxt) == true)
    inputFiles_->addMultizoneSurfaceFile(file_name);
  else
    errTxt += "<top-surface-file> needs to be given in <background><multizone-model>\n";

  int erosion_priority;
  if(parseValue(root, "top-surface-erosion-priority", erosion_priority, errTxt) == true)
    modelSettings_->addErosionPriority(erosion_priority);
  else
    errTxt += "<top-surface-erosion-priority> needs to be given in <background><multizone-model>\n";

  modelSettings_->addCorrelationStructure(ModelSettings::TOP); //Dummy variable used for balance in vector

  while(parseZone(root,errTxt) == true);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseZone(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("zone");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("base-surface-file");
  legalCommands.push_back("erosion-priority");
  legalCommands.push_back("correlation-structure");

  std::string file_name;
  if(parseFileName(root, "base-surface-file", file_name, errTxt) == true)
    inputFiles_->addMultizoneSurfaceFile(file_name);
  else
    errTxt += "<base-surface-file> needs to be given for all zones in <background><multizone-model><zone>\n";

  int erosion_priority;
  if(parseValue(root, "erosion-priority", erosion_priority, errTxt) == true)
    modelSettings_->addErosionPriority(erosion_priority);
  else
    errTxt += "<erosion-priority> needs to be given for all zones in <background><multizone-model><zone>\n";

  std::string correlation_structure;
  if(parseValue(root, "correlation-structure", correlation_structure, errTxt) == true) {
    if(correlation_structure == "top")
      modelSettings_->addCorrelationStructure(ModelSettings::TOP);
    else if(correlation_structure == "base")
      modelSettings_->addCorrelationStructure(ModelSettings::BASE);
    else if(correlation_structure == "compaction")
      modelSettings_->addCorrelationStructure(ModelSettings::COMPACTION);
    else
      errTxt += "<correlation-structure> in <background><multizone-model><zone> needs to be either 'top', 'base' or 'compaction'\n";
  }
  else
    errTxt += "<correlation-structure> needs to be given for all zones in <background><multizone-model><zone>\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseFaciesProbabilities(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("facies-probabilities");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("estimation-interval");
  legalCommands.push_back("prior-probabilities");
  legalCommands.push_back("uncertainty-level");
  legalCommands.push_back("use-vs");
  legalCommands.push_back("use-prediction");
  legalCommands.push_back("use-absolute-elastic-parameters");

  parseFaciesEstimationInterval(root, errTxt);

  parsePriorFaciesProbabilities(root, errTxt);

  float value;
  if(parseValue(root, "uncertainty-level", value, errTxt) == true)
    modelSettings_->setPundef(value);

  bool useVs = true;
  if(parseBool(root, "use-vs", useVs, errTxt) == true && useVs == false)
    modelSettings_->setNoVsFaciesProb(true); //Note: What we store is inverse of what we ask for.

  bool useInversion = false;
  if(parseBool(root, "use-prediction", useInversion, errTxt) == true && useInversion == true)
    modelSettings_->setUseFilterForFaciesProb(false); //Note: What we store is inverse of what we ask for.

  bool absolute;
  if(parseBool(root, "use-absolute-elastic-parameters", absolute, errTxt) == true && absolute == true)
    modelSettings_->setFaciesProbRelative(false);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseFaciesEstimationInterval(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("estimation-interval");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-surface-file");
  legalCommands.push_back("base-surface-file");

  std::string filename;
  if(parseFileName(root, "top-surface-file", filename, errTxt) == true)
    inputFiles_->setFaciesEstIntFile(0, filename);
  else
    errTxt += "Must specify <top-surface-file> in command <"+root->ValueStr()+"> "+
      lineColumnText(root)+".\n";
  if(parseFileName(root, "base-surface-file", filename, errTxt) == true)
    inputFiles_->setFaciesEstIntFile(1, filename);
  else
    errTxt += "Must specify <base-surface-file> in command <"+root->ValueStr()+">"+
      lineColumnText(root)+".\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parsePriorFaciesProbabilities(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("prior-probabilities");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("facies");

  float sum       = 0.0;
  int   status    = 0;
  int   oldStatus = 0;

  while(parseFacies(root,errTxt)==true)
  {
    status = modelSettings_->getIsPriorFaciesProbGiven();
    if(oldStatus!=0 &&oldStatus!=status)
      errTxt+= "Prior facies probability must be given in the same way for all facies.\n";

    oldStatus = status;

  }
  if(status==1)
  {
    typedef std::map<std::string,float> mapType;
    mapType myMap = modelSettings_->getPriorFaciesProb();
    for(mapType::const_iterator it=myMap.begin();it!=myMap.end();++it)
    {
      if (it->second < 0 || it->second > 1)
        errTxt += "Prior facies probabilities must be between 0 and 1.\n";
      sum+=(*it).second;
    }
    if(sum!=1.0)
    {
      errTxt+="Prior facies probabilities must sum to 1.0. They sum to "+ NRLib::ToString(sum) +".\n";
    }
  }
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseFacies(TiXmlNode * node, std::string & errTxt)
{
TiXmlNode * root = node->FirstChildElement("facies");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("name");
  legalCommands.push_back("probability");
  legalCommands.push_back("probability-cube");

  std::string faciesname;
  std::string filename;
  float value;
  parseValue(root, "name", faciesname,errTxt, true);

  if(parseValue(root,"probability",value,errTxt,true)==true)
  {
    modelSettings_->setPriorFaciesProbGiven(ModelSettings::FACIES_FROM_MODEL_FILE);
    modelSettings_->addPriorFaciesProb(faciesname,value);
  }
  else if(parseValue(root,"probability-cube",filename,errTxt,true)==true)
  {
    modelSettings_->setPriorFaciesProbGiven(ModelSettings::FACIES_FROM_CUBES);
    inputFiles_->setPriorFaciesProb(faciesname,filename);
  }

 checkForJunk(root, errTxt, legalCommands, true); //allow duplicates
  return(true);


}

bool
XmlModelFile::parseRockPhysics(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("rock-physics");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("reservoir");
  legalCommands.push_back("evolve");
  legalCommands.push_back("predefinitions");
  legalCommands.push_back("rock");
  legalCommands.push_back("trend-cube");

  // The order of the parse-commands should not be changed
  parseReservoir(root, errTxt);
  while(parseEvolve(root, errTxt));
  parsePredefinitions(root, errTxt);
  std::string dummy;
  while(parseRock(root, dummy, errTxt));

  int count = 0;
  while(parseTrendCube(root, errTxt) == true)
    count++;
  if(count > 2)
    errTxt += "The maximum allowed number of trend cubes in <rock-physics><trend-cube> is two.\n";

  modelSettings_->setFaciesProbFromRockPhysics(true);

  modelSettings_->setFaciesProbRelative(false);

  modelSettings_->setBackgroundFromRockPhysics(true);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parsePredefinitions(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("predefinitions");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("fluid");
  legalCommands.push_back("solid");
  legalCommands.push_back("dry-rock");

  std::string label;
  while(parseFluid(root, label, errTxt));
  while(parseSolid(root, label, errTxt));
  while(parseDryRock(root, label, errTxt));

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}
bool
XmlModelFile::parseRock(TiXmlNode * node, std::string & label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("rock");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("label");
  legalCommands.push_back("use");
  legalCommands.push_back("tabulated");
  legalCommands.push_back("reuss");
  legalCommands.push_back("voigt");
  legalCommands.push_back("hill");
  legalCommands.push_back("dem");
  legalCommands.push_back("gassmann");
  legalCommands.push_back("bounding");

  label = "";
  parseValue(root, "label", label, errTxt);

  int given = 0;

  std::string use = "";
  parseValue(root, "use", use, errTxt);

  if(use != "" && label != "")
    errTxt += "Both <label> and <use> can not be given in <rock>\n";
  else if(use != "") {
    label = use;
    given++;
  }
  int constituent_type = ModelSettings::ROCK;

  std::vector<DistributionWithTrendStorage *> dummy_porosity(1,NULL);
  std::string                                 dummy_moduli   = "";

  if(parseTabulated(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseReuss(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseVoigt(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseHill(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseDEM(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseGassmann(root, constituent_type, label, errTxt) == true)
    given++;
  if(parseBounding(root, constituent_type, label, errTxt) == true)
    given++;

  if(given == 0)
    errTxt += "A theory or <use> needs to follow after <label> in <rock>\n";
  else if(given > 1)
    errTxt += "Only one theory or <use> can be given after <label> in <rock>\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseSolid(TiXmlNode * node, std::string & label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("solid");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("label");
  legalCommands.push_back("use");
  legalCommands.push_back("tabulated");
  legalCommands.push_back("reuss");
  legalCommands.push_back("voigt");
  legalCommands.push_back("hill");
  legalCommands.push_back("dem");
  //legalCommands.push_back("hashin-shtrikman");
  //legalCommands.push_back("walton");

  label = "";
  parseValue(root, "label", label, errTxt);

  int given = 0;

  std::string use = "";
  parseValue(root, "use", use, errTxt);

  if(use != "" && label != "")
    errTxt += "Both <label> and <use> can not be given in <solid>\n";
  else if(use != "") {
    label = use;
    given++;
  }

  int constituent_type = ModelSettings::SOLID;

  std::vector<DistributionWithTrendStorage *> dummy_porosity(1,NULL);
  std::string                                 dummy_moduli = "";

  if(parseTabulated(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseReuss(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseVoigt(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseHill(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseDEM(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;

  if(given == 0)
    errTxt += "A theory or <use> needs to follow after <label> in <solid>\n";
  else if(given > 1)
    errTxt += "Only one theory or <use> can be given after <label> in <solid>\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseDryRock(TiXmlNode * node, std::string & label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("dry-rock");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("label");
  legalCommands.push_back("use");
  legalCommands.push_back("tabulated");
  legalCommands.push_back("reuss");
  legalCommands.push_back("voigt");
  legalCommands.push_back("hill");
  legalCommands.push_back("dem");
  //legalCommands.push_back("hashin-shtrikman");
  //legalCommands.push_back("walton");
  legalCommands.push_back("total-porosity");
  legalCommands.push_back("mineral-moduli");

  label = "";
  parseValue(root, "label", label, errTxt);

  int given = 0;

  std::string use = "";
  parseValue(root, "use", use, errTxt);

  if(use != "" && label != "")
    errTxt += "Both <label> and <use> can not be given in <dry-rock>\n";
  else if(use != "") {
    label = use;
    given++;
  }

  std::string                                 dummy = "";
  std::vector<DistributionWithTrendStorage *> total_porosity;

  if(parseDistributionWithTrend(root, "total-porosity", total_porosity, dummy, false, errTxt) == false)
      errTxt += "The total porosity must be given for the dry-rock\n";

  std::string moduli;
  if(parseMineralModuli(root, moduli, errTxt) == false)
    errTxt += "The mineral moduli must be given for the dry-rock\n";

  int constituent_type = ModelSettings::DRY_ROCK;

  if(parseTabulated(root, constituent_type, label, total_porosity, moduli, errTxt) == true)
    given++;
  if(parseReuss(root, constituent_type, label, total_porosity, moduli, errTxt) == true)
    given++;
  if(parseVoigt(root, constituent_type, label, total_porosity, moduli, errTxt) == true)
    given++;
  if(parseHill(root, constituent_type, label, total_porosity, moduli, errTxt) == true)
    given++;
  if(parseDEM(root, constituent_type, label, total_porosity, moduli, errTxt) == true)
    given++;

 if(given == 0)
    errTxt += "A theory or <use> needs to follow after <label> in <dry-rock>\n";
  else if(given > 1)
    errTxt += "Only one theory or <use> can be given after <label> in <dry-rock>\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseFluid(TiXmlNode * node, std::string & label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("fluid");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("label");
  legalCommands.push_back("use");
  legalCommands.push_back("tabulated");
  legalCommands.push_back("reuss");
  legalCommands.push_back("voigt");
  legalCommands.push_back("hill");
  legalCommands.push_back("batzle-wang-brine");
  //legalCommands.push_back("hashin-shtrikman");
  //legalCommands.push_back("walton");

  label = "";
  parseValue(root, "label", label, errTxt);

  int given = 0;

  std::string use = "";
  parseValue(root, "use", use, errTxt);

  if(use != "" && label != "")
    errTxt += "Both <label> and <use> can not be given in <fluid>\n";
  else if(use != "") {
    label = use;
    given++;
  }

  int constituent_type = ModelSettings::FLUID;

  std::vector<DistributionWithTrendStorage *> dummy_porosity;
  std::string                                 dummy_moduli = "";

  if(parseTabulatedFluid(root, constituent_type, label, errTxt) == true)
    given++;
  if(parseReuss(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseVoigt(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseHill(root, constituent_type, label, dummy_porosity, dummy_moduli, errTxt) == true)
    given++;
  if(parseBatzleWangBrine(root, constituent_type, label, errTxt) == true)
    given++;

  if(given == 0)
    errTxt += "A theory or <use> needs to follow after <label> in <fluid>\n";
  else if(given > 1)
    errTxt += "Only one theory or <use> can be given after <label> in <fluid>\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseReuss(TiXmlNode                                   * node,
                         int                                           constituent,
                         std::string                                   label,
                         std::vector<DistributionWithTrendStorage *>   total_porosity,
                         std::string                                   moduli,
                         std::string                                 & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("reuss");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("constituent");

  std::vector<std::string>                                  constituent_label;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_fraction;

  std::string                                 this_label;
  std::vector<DistributionWithTrendStorage *> volume_fraction;

  while(parseConstituent(root, this_label, volume_fraction, errTxt) == true) {
    constituent_label.push_back(this_label);
    constituent_fraction.push_back(volume_fraction);
    volume_fraction.clear();
  }

  if(constituent_label.size() < 2)
    errTxt += "At least two constituents must be given in the Reuss rock physics model for "+label+"\n";

  if(constituent == ModelSettings::FLUID) {
    DistributionsFluidStorage * fluid = new ReussFluidStorage(constituent_label, constituent_fraction);
    modelSettings_->addFluid(label, fluid);
  }
  else if(constituent == ModelSettings::SOLID) {
    DistributionsSolidStorage * solid = new ReussSolidStorage(constituent_label, constituent_fraction);
    modelSettings_->addSolid(label, solid);
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    DistributionsDryRockStorage * dry_rock = new ReussDryRockStorage(constituent_label, constituent_fraction, total_porosity, moduli);
    modelSettings_->addDryRock(label, dry_rock);
  }
  else if(constituent == ModelSettings::ROCK) {
    DistributionsRockStorage * rock = new ReussRockStorage(constituent_label, constituent_fraction);
    modelSettings_->addRock(label, rock);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}
bool
XmlModelFile::parseVoigt(TiXmlNode                                   * node,
                         int                                           constituent,
                         std::string                                   label,
                         std::vector<DistributionWithTrendStorage *>   total_porosity,
                         std::string                                   moduli,
                         std::string                                 & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("voigt");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("constituent");

  std::vector<std::string>                                  constituent_label;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_fraction;

  std::string                                 this_label;
  std::vector<DistributionWithTrendStorage *> volume_fraction;

  while(parseConstituent(root, this_label, volume_fraction, errTxt) == true) {
    constituent_label.push_back(this_label);
    constituent_fraction.push_back(volume_fraction);
    volume_fraction.clear();
  }

  if(constituent_label.size() < 2)
    errTxt += "At least two constituents must be given in the Voigt rock physics model for "+label+"\n";

  if(constituent == ModelSettings::FLUID) {
    DistributionsFluidStorage * fluid = new VoigtFluidStorage(constituent_label, constituent_fraction);
    modelSettings_->addFluid(label, fluid);
  }
  else if(constituent == ModelSettings::SOLID) {
    DistributionsSolidStorage * solid = new VoigtSolidStorage(constituent_label, constituent_fraction);
    modelSettings_->addSolid(label, solid);
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    DistributionsDryRockStorage * dry_rock = new VoigtDryRockStorage(constituent_label, constituent_fraction, total_porosity, moduli);
    modelSettings_->addDryRock(label, dry_rock);
  }
  else if(constituent == ModelSettings::ROCK) {
    DistributionsRockStorage * rock = new VoigtRockStorage(constituent_label, constituent_fraction);
    modelSettings_->addRock(label, rock);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}
bool
XmlModelFile::parseHill(TiXmlNode                                   * node,
                        int                                           constituent,
                        std::string                                   label,
                        std::vector<DistributionWithTrendStorage *>   total_porosity,
                        std::string                                   moduli,
                        std::string                                 & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("hill");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("constituent");

  std::vector<std::string>                                  constituent_label;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_fraction;

  std::string                                 this_label;
  std::vector<DistributionWithTrendStorage *> volume_fraction;

  while(parseConstituent(root, this_label, volume_fraction, errTxt) == true) {
    constituent_label.push_back(this_label);
    constituent_fraction.push_back(volume_fraction);
    volume_fraction.clear();
  }

  if(constituent_label.size() < 2)
    errTxt += "At least two constituents must be given in the Hill rock physics model for "+label+"\n";

  if(constituent == ModelSettings::FLUID) {
    DistributionsFluidStorage * fluid = new HillFluidStorage(constituent_label, constituent_fraction);
    modelSettings_->addFluid(label, fluid);
  }
  else if(constituent == ModelSettings::SOLID) {
    DistributionsSolidStorage * solid = new HillSolidStorage(constituent_label, constituent_fraction);
    modelSettings_->addSolid(label, solid);
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    DistributionsDryRockStorage * dry_rock = new HillDryRockStorage(constituent_label, constituent_fraction, total_porosity, moduli);
    modelSettings_->addDryRock(label, dry_rock);
  }
  else if(constituent == ModelSettings::ROCK) {
    DistributionsRockStorage * rock = new HillRockStorage(constituent_label, constituent_fraction);
    modelSettings_->addRock(label, rock);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseConstituent(TiXmlNode                                   * node,
                               std::string                                 & constituent_label,
                               std::vector<DistributionWithTrendStorage *> & volume_fraction,
                               std::string                                 & errTxt)
{

  TiXmlNode * root = node->FirstChildElement("constituent");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("solid");
  legalCommands.push_back("fluid");
  legalCommands.push_back("dry-rock");
  legalCommands.push_back("volume-fraction");

  int constituent_given = 0;
  while(parseSolid(root, constituent_label, errTxt) == true)
    constituent_given++;
  while(parseFluid(root, constituent_label, errTxt) == true)
    constituent_given++;
  while(parseDryRock(root, constituent_label, errTxt) == true)
    constituent_given++;

  if(constituent_given == 0)
    errTxt += "A solid, fluid or dry-rock needs to be given in each constituent of the rock physics model\n";
  else if(constituent_given > 1)
    errTxt += "There can only be one element in each constituent of the rock physics model\n";

  std::string volume_label = "";
  if(parseDistributionWithTrend(root, "volume-fraction", volume_fraction, volume_label, false, errTxt) == false)
    volume_fraction.push_back(NULL);

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseBatzleWangBrine(TiXmlNode * node, int constituent, std::string label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("batzle-wang-brine");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("pore-pressure");
  legalCommands.push_back("temperature");
  legalCommands.push_back("salinity");

  std::string dummy;

  std::vector<DistributionWithTrendStorage *> pore_pressure;
  if(parseDistributionWithTrend(root, "pore-pressure", pore_pressure, dummy, false, errTxt) == false)
    errTxt += "The pore pressure must be given in the Batzle-Wang brine model\n";

  std::vector<DistributionWithTrendStorage *> temperature;
  if(parseDistributionWithTrend(root, "temperature", temperature, dummy, false, errTxt) == false)
    errTxt += "The temperature must be given in the Batzle-Wang brine model\n";

  std::vector<DistributionWithTrendStorage *> salinity;
  if(parseDistributionWithTrend(root, "salinity", salinity, dummy, false, errTxt) == false)
    errTxt += "The salinity must be given in the Batzle-Wang brine model\n";

  if(constituent == ModelSettings::FLUID) {
    DistributionsFluidStorage * fluid = new BatzleWangFluidStorage(pore_pressure, temperature, salinity);
    modelSettings_->addFluid(label, fluid);
  }
  else if(constituent == ModelSettings::SOLID) {
    errTxt += "Implementation error: The Batzle-Wang brine model can not be used to make a solid\n";
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    errTxt += "Implementation error: The Batzle-Wang brine model can not be used to make a dry-rock\n";
  }
  else if(constituent == ModelSettings::ROCK) {
    errTxt += "Implementation error: The Batzle-Wang brine model can not be used to make a rock\n";
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseDEM(TiXmlNode                                  * node,
                       int                                          constituent,
                       std::string                                  label,
                       std::vector<DistributionWithTrendStorage *>  total_porosity,
                       std::string                                  moduli,
                       std::string                                & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("dem");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("host");
  legalCommands.push_back("inclusion");

  std::vector<DistributionWithTrendStorage *> host_volume;
  std::string                                 host_label  = "";

  bool missing_vol_frac_host = false;
  if(parseDEMHost(root, host_label, host_volume, errTxt, missing_vol_frac_host) == false)
    errTxt += "The host must be given in the DEM rock physics model for "+label+"\n";

  if (missing_vol_frac_host == true)
    host_volume = std::vector<DistributionWithTrendStorage* >(1, NULL);

  std::vector<std::string>                                  inclusion_label;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume;
  std::vector<std::vector<DistributionWithTrendStorage *> > aspect_ratio;

  std::vector<DistributionWithTrendStorage *> volume;
  std::vector<DistributionWithTrendStorage *> aspect;
  std::string                                 this_label = "";

  bool missing_vol_frac = false;
  int  counter_missing_vol_frac = 0;

  while(parseDEMInclusion(root, this_label, aspect, volume, errTxt, missing_vol_frac) == true) {
    inclusion_label.push_back(this_label);
    if (missing_vol_frac == true) {
      counter_missing_vol_frac++;
      inclusion_volume.push_back(std::vector<DistributionWithTrendStorage* >(1, NULL));
    }
    else
      inclusion_volume.push_back(volume);

    aspect_ratio.push_back(aspect);
  }

  //check missing volume fractions
  if (missing_vol_frac_host == false) {
    if (counter_missing_vol_frac != 1)
      errTxt += "All the volume fractions except one must be given for the inclusions of the DEM model when the volume fraction for host is given.\n";
  }
  else {
    if (counter_missing_vol_frac != 0)
      errTxt += "All the volume fractions must be given for the inclusions of the DEM model when the volume fraction for host is not given.\n";
  }

  if(inclusion_label.size() < 1)
    errTxt += "At least one inclusion must be given in the DEM rock physics model for "+label+"\n";

  if(constituent == ModelSettings::FLUID) {
    errTxt += "Implementation error: The DEM model can not be used to mix a fluid\n";
  }
  else if(constituent == ModelSettings::SOLID) {
    DistributionsSolidStorage * solid = new DEMSolidStorage(host_label, host_volume, inclusion_label, inclusion_volume, aspect_ratio);
    modelSettings_->addSolid(label, solid);
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    DistributionsDryRockStorage * dry_rock = new DEMDryRockStorage(host_label, host_volume, inclusion_label, inclusion_volume, aspect_ratio, total_porosity, moduli);
    modelSettings_->addDryRock(label, dry_rock);
  }
  else if(constituent == ModelSettings::ROCK) {
    DistributionsRockStorage * rock = new DEMRockStorage(host_label, host_volume, inclusion_label, inclusion_volume, aspect_ratio);
    modelSettings_->addRock(label, rock);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseDEMHost(TiXmlNode                                   * node,
                           std::string                                 & label,
                           std::vector<DistributionWithTrendStorage *> & volume_fraction,
                           std::string                                 & errTxt,
                           bool                                        & missing_vol_frac)
{

  TiXmlNode * root = node->FirstChildElement("host");
  if(root == 0)
    return(false);

  missing_vol_frac = false;

  std::vector<std::string> legalCommands;
  legalCommands.push_back("solid");
  legalCommands.push_back("fluid");
  legalCommands.push_back("dry-rock");
  legalCommands.push_back("volume-fraction");

  int host_given = 0;
  while(parseSolid(root, label, errTxt) == true)
    host_given++;
  while(parseFluid(root, label, errTxt) == true)
    host_given++;
  while(parseDryRock(root, label, errTxt) == true)
    host_given++;

  if(host_given > 1)
    errTxt += "There can only be one constituent in the host of the DEM model\n";

  std::string dummy;
  if(parseDistributionWithTrend(root, "volume-fraction", volume_fraction, dummy, false, errTxt) == false) {
    //errTxt += "The volume fraction must be given for the host of the DEM model\n";
    missing_vol_frac = true;
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseDEMInclusion(TiXmlNode                                   * node,
                                std::string                                 & label,
                                std::vector<DistributionWithTrendStorage *> & aspect_ratio,
                                std::vector<DistributionWithTrendStorage *> & volume_fraction,
                                std::string                                 & errTxt,
                                bool                                        & missing_vol_frac)
{
  missing_vol_frac = false;
  TiXmlNode * root = node->FirstChildElement("inclusion");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("solid");
  legalCommands.push_back("fluid");
  legalCommands.push_back("dry-rock");
  legalCommands.push_back("volume-fraction");
  legalCommands.push_back("aspect-ratio");

  int inclusion_given = 0;
  while(parseSolid(root, label, errTxt) == true)
    inclusion_given++;
  while(parseFluid(root, label, errTxt) == true)
    inclusion_given++;
  while(parseDryRock(root, label, errTxt) == true)
    inclusion_given++;

  if(inclusion_given > 1)
    errTxt += "There can only be one constituent for each of the inclusions of the DEM model\n";

  std::string dummy;
  if(parseDistributionWithTrend(root, "volume-fraction", volume_fraction, dummy, false, errTxt) == false) {
    //errTxt += "The volume fraction must be given for the inclusions of the DEM model\n";
    missing_vol_frac = true;
  }

  if(parseDistributionWithTrend(root, "aspect-ratio", aspect_ratio, dummy, false, errTxt) == false)
    errTxt += "The aspect ratio must be given for the inclusions of the DEM model\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseGassmann(TiXmlNode * node, int constituent, std::string label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("gassmann");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("dry-rock");
  legalCommands.push_back("fluid");

  int dry_rock_given = 0;
  std::string dry_rock;
  while(parseDryRock(root, dry_rock, errTxt) == true)
    dry_rock_given++;

  if(dry_rock_given == 0)
    errTxt += "A dry-rock must be given in the Gassmann model\n";
  else if(dry_rock_given > 1)
    errTxt += "Only one dry-rock can be given in the Gassmann model\n";

  int fluid_given = 0;
  std::string fluid;
  while(parseFluid(root, fluid, errTxt) == true)
    fluid_given++;

  if(fluid_given == 0)
    errTxt += "A fluid must be given in the Gassmann model\n";
  else if(fluid_given > 1)
    errTxt += "Only one fluid can be given in the Gassmann model\n";

  if(constituent == ModelSettings::FLUID) {
    errTxt += "Implementation error: The Gassmann model can not be used to make a fluid\n";
  }
  else if(constituent == ModelSettings::SOLID) {
    errTxt += "Implementation error: The Gassmann model can not be used to make a solid\n";
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    errTxt += "Implementation error: The Gassmann model can not be used to make a dry-rock\n";
  }
  else if(constituent == ModelSettings::ROCK) {
    DistributionsRockStorage * rock = new GassmannRockStorage(dry_rock, fluid);
    modelSettings_->addRock(label, rock);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseBounding(TiXmlNode * node, int constituent, std::string label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("bounding");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("upper-bound");
  legalCommands.push_back("lower-bound");
  legalCommands.push_back("bulk-modulus-weight");
  legalCommands.push_back("shear-modulus-weight");
  legalCommands.push_back("correlation-weights");
  legalCommands.push_back("porosity");

  std::string upper_bound;
  if(parseUpperBound(root, upper_bound, errTxt) == false)
    errTxt += "<upper-bound> needs to be given in <rock><bounding>\n";

  std::string lower_bound;
  if(parseLowerBound(root, lower_bound, errTxt) == false)
    errTxt += "<lower-bound> needs to be given in <rock><bounding>\n";

  std::string dummy;
  std::vector<DistributionWithTrendStorage *> bulk_weight;
  if(parseDistributionWithTrend(root, "bulk-modulus-weight", bulk_weight, dummy, false, errTxt, true) == false)
    errTxt += "<bulk-modulus-weight> needs to be given in <rock><bounding>\n";

  std::vector<DistributionWithTrendStorage *> shear_weight;
  if(parseDistributionWithTrend(root, "shear-modulus-weight", shear_weight, dummy, false, errTxt, true) == false)
    errTxt += "<shear-weight> needs to be given in <rock><bounding>\n";

  double correlation;
  if(parseValue(root, "correlation-weights", correlation, errTxt) == false)
    correlation = 0;

  std::vector<DistributionWithTrendStorage *> porosity;
  if(parseDistributionWithTrend(root, "porosity", porosity, dummy, false, errTxt, true) == false)
    errTxt += "<porosity> needs to be given in <rock><bounding>\n";

  if(constituent == ModelSettings::FLUID) {
    errTxt += "Implementation error: The Bounding model can not be used to make a fluid\n";
  }
  else if(constituent == ModelSettings::SOLID) {
    errTxt += "Implementation error: The Bounding model can not be used to make a solid\n";
  }
  else if(constituent == ModelSettings::DRY_ROCK) {
    errTxt += "Implementation error: The Boundning model can not be used to make a dry-rock\n";
  }
  else if(constituent == ModelSettings::ROCK) {
    DistributionsRockStorage * rock = new BoundingRockStorage(upper_bound, lower_bound, porosity, bulk_weight, shear_weight, correlation);
    modelSettings_->addRock(label, rock);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseUpperBound(TiXmlNode * node, std::string & label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("upper-bound");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("rock");

  int rock_given = 0;
  while(parseRock(root, label, errTxt) == true)
    rock_given++;

  if(rock_given == 0)
    errTxt += "A rock must be given for the upper bound in the Bounding model\n";
  else if(rock_given > 1)
    errTxt += "Only one rock can be given for the upper bound in the Bounding model\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseLowerBound(TiXmlNode * node, std::string & label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("lower-bound");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("rock");

  int rock_given = 0;
  while(parseRock(root, label, errTxt) == true)
    rock_given++;

  if(rock_given == 0)
    errTxt += "A rock must be given for the lower bound in the Bounding model\n";
  else if(rock_given > 1)
    errTxt += "Only one rock can be given for the lower bound in the Bounding model\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseMineralModuli(TiXmlNode * node, std::string & moduli, std::string & errTxt)
{

  TiXmlNode * root = node->FirstChildElement("mineral-moduli");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("use");

  moduli = "";
  if(parseValue(root, "use", moduli, errTxt) == false)
    errTxt += "The mineral moduli in the dry rock must be specified with a variable using <use>\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseTabulated(TiXmlNode                                   * node,
                             int                                           constituent,
                             std::string                                   label,
                             std::vector<DistributionWithTrendStorage *>   total_porosity,
                             std::string                                   moduli,
                             std::string                                 & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("tabulated");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("density");
  legalCommands.push_back("vp");
  legalCommands.push_back("vs");
  legalCommands.push_back("correlation-vp-vs");
  legalCommands.push_back("correlation-vp-density");
  legalCommands.push_back("correlation-vs-density");
  legalCommands.push_back("bulk-modulus");
  legalCommands.push_back("shear-modulus");
  legalCommands.push_back("correlation-bulk-shear");
  legalCommands.push_back("correlation-bulk-density");
  legalCommands.push_back("correlation-shear-density");

  std::string dummy;

  bool use_vp      = false;
  bool use_modulus = false;

  std::vector<DistributionWithTrendStorage *> vp;
  if(parseDistributionWithTrend(root, "vp", vp, dummy, false, errTxt, true) == true)
    use_vp = true;

  std::vector<DistributionWithTrendStorage *> vs;
  if(parseDistributionWithTrend(root, "vs", vs, dummy, false, errTxt, true) == false && use_vp == true)
    errTxt +="Both <vp> and <vs> need to be given in <solid><tabulated>\n";

  std::vector<DistributionWithTrendStorage *> bulk_modulus;
  if(parseDistributionWithTrend(root, "bulk-modulus", bulk_modulus, dummy, false, errTxt, true) == true)
    use_modulus = true;

  std::vector<DistributionWithTrendStorage *> shear_modulus;
  if(parseDistributionWithTrend(root, "shear-modulus", shear_modulus, dummy, false, errTxt, true) == false && use_modulus == true)
    errTxt +="Both <bulk-modulus> and <shear-modulus> need to be given in <solid><tabulated>\n";

  if(use_vp == true && use_modulus == true)
    errTxt += "Both <vp> and <bulk-modulus> can not be used in <solid><tabulated>\n";
  else if(use_vp == false && use_modulus == false)
    errTxt += "One of <vp> or <bulk-modulus> must be used in <solid><tabulated>\n";

  std::vector<DistributionWithTrendStorage *> density;
  if(parseDistributionWithTrend(root, "density", density, dummy, false, errTxt, true) == false)
    errTxt += "<density> needs to be specified in <solid><tabulated>\n";

  std::vector<DistributionWithTrendStorage *> correlation_vp_vs_with_trend;
  std::vector<double>                         correlation_vp_vs;
  if(parseDistributionWithTrend(root, "correlation-vp-vs", correlation_vp_vs_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_vp_vs_with_trend, "correlation", correlation_vp_vs, errTxt);
    for(size_t i=0; i<correlation_vp_vs.size(); i++) {
      if(correlation_vp_vs[i] > 1 || correlation_vp_vs[i] < -1)
        errTxt += "<correlation-vp-vs> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_vp_vs_with_trend[0]->GetIsShared() == false)
      delete correlation_vp_vs_with_trend[0];
  }
  else
    correlation_vp_vs.push_back(modelSettings_->getDefaultCorrelationVpVs());

  std::vector<DistributionWithTrendStorage *> correlation_vp_density_with_trend;
  std::vector<double>                         correlation_vp_density;
  if(parseDistributionWithTrend(root, "correlation-vp-density", correlation_vp_density_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_vp_density_with_trend, "correlation", correlation_vp_density, errTxt);
    for(size_t i=0; i<correlation_vp_density.size(); i++) {
      if(correlation_vp_density[i] > 1 || correlation_vp_density[i] < -1)
        errTxt += "<correlation-vp-density> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_vp_density_with_trend[0]->GetIsShared() == false)
      delete correlation_vp_density_with_trend[0];
  }
  else
    correlation_vp_density.push_back(0.0);

  std::vector<DistributionWithTrendStorage *> correlation_vs_density_with_trend;
  std::vector<double>                         correlation_vs_density;
  if(parseDistributionWithTrend(root, "correlation-vs-density", correlation_vs_density_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_vs_density_with_trend, "correlation", correlation_vs_density, errTxt);
    for(size_t i=0; i<correlation_vs_density.size(); i++) {
      if(correlation_vs_density[i] > 1 || correlation_vs_density[i] < -1)
        errTxt += "<correlation-vs-density> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_vs_density_with_trend[0]->GetIsShared() == false)
      delete correlation_vs_density_with_trend[0];
  }
  else
    correlation_vs_density.push_back(0.0);

  std::vector<DistributionWithTrendStorage *> correlation_bulk_shear_with_trend;
  std::vector<double>                         correlation_bulk_shear;
  if(parseDistributionWithTrend(root, "correlation-bulk-shear", correlation_bulk_shear_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_bulk_shear_with_trend, "correlation", correlation_bulk_shear, errTxt);
    for(size_t i=0; i<correlation_bulk_shear.size(); i++) {
      if(correlation_bulk_shear[i] > 1 || correlation_bulk_shear[i] < -1)
        errTxt += "<correlation-bulk-shear> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_bulk_shear_with_trend[0]->GetIsShared() == false)
      delete correlation_bulk_shear_with_trend[0];
  }
  else
    correlation_bulk_shear.push_back(modelSettings_->getDefaultCorrelationVpVs());

  std::vector<DistributionWithTrendStorage *> correlation_bulk_density_with_trend;
  std::vector<double>                         correlation_bulk_density;
  if(parseDistributionWithTrend(root, "correlation-bulk-density", correlation_bulk_density_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_bulk_density_with_trend, "correlation", correlation_bulk_density, errTxt);
    for(size_t i=0; i<correlation_bulk_density.size(); i++) {
      if(correlation_bulk_density[i] > 1 || correlation_bulk_density[i] < -1)
        errTxt += "<correlation-bulk-density> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_bulk_density_with_trend[0]->GetIsShared() == false)
      delete correlation_bulk_density_with_trend[0];
  }
  else
    correlation_bulk_density.push_back(0.0);


  std::vector<DistributionWithTrendStorage *> correlation_shear_density_with_trend;
  std::vector<double>                         correlation_shear_density;
  if(parseDistributionWithTrend(root, "correlation-shear-density", correlation_shear_density_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_shear_density_with_trend, "correlation", correlation_shear_density, errTxt);
    for(size_t i=0; i<correlation_shear_density.size(); i++) {
      if(correlation_shear_density[i] > 1 || correlation_shear_density[i] < -1)
        errTxt += "<correlation-shear-density> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_shear_density_with_trend[0]->GetIsShared() == false)
      delete correlation_shear_density_with_trend[0];
  }
  else
    correlation_shear_density.push_back(0.0);


  if(use_vp) {
    if(constituent == ModelSettings::FLUID) {
      errTxt += "Implementation error: parseTabulated() can not be used for fluids. Use parseTabulatedFluid()\n";
    }
    else if(constituent == ModelSettings::SOLID) {
      DistributionsSolidStorage * solid = new TabulatedVelocitySolidStorage(vp, vs, density, correlation_vp_vs, correlation_vp_density, correlation_vs_density);
      modelSettings_->addSolid(label, solid);
    }
    else if(constituent == ModelSettings::DRY_ROCK) {
      DistributionsDryRockStorage * dry_rock = new TabulatedVelocityDryRockStorage(vp, vs, density, correlation_vp_vs, correlation_vp_density, correlation_vs_density, total_porosity, moduli);
      modelSettings_->addDryRock(label, dry_rock);
  }
    else if(constituent == ModelSettings::ROCK) {
      DistributionsRockStorage * rock = new TabulatedVelocityRockStorage(vp, vs, density, correlation_vp_vs, correlation_vp_density, correlation_vs_density);
      modelSettings_->addRock(label, rock);
    }
  }
  else {
    if(constituent == ModelSettings::FLUID) {
      errTxt += "Implementation error: parseTabulated() can not be used for fluids. Use parseTabulatedFluid()\n";
    }
    else if(constituent == ModelSettings::SOLID) {
      DistributionsSolidStorage * solid = new TabulatedModulusSolidStorage(bulk_modulus, shear_modulus, density, correlation_bulk_shear, correlation_bulk_density, correlation_shear_density);
      modelSettings_->addSolid(label, solid);
    }
    else if(constituent == ModelSettings::DRY_ROCK) {
      DistributionsDryRockStorage * dry_rock = new TabulatedModulusDryRockStorage(bulk_modulus, shear_modulus, density, correlation_bulk_shear, correlation_bulk_density, correlation_shear_density, total_porosity, moduli);
      modelSettings_->addDryRock(label, dry_rock);
  }
    else if(constituent == ModelSettings::ROCK) {
      DistributionsRockStorage * rock = new TabulatedModulusRockStorage(bulk_modulus, shear_modulus, density, correlation_bulk_shear, correlation_bulk_density, correlation_shear_density);
      modelSettings_->addRock(label, rock);
    }
  }
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseTabulatedFluid(TiXmlNode * node, int constituent, std::string label, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("tabulated");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("density");
  legalCommands.push_back("vp");
  legalCommands.push_back("correlation-vp-density");
  legalCommands.push_back("bulk-modulus");
  legalCommands.push_back("correlation-bulk-density");

  std::string dummy;

  bool use_vp      = false;
  bool use_modulus = false;

  std::vector<DistributionWithTrendStorage *> vp;
  if(parseDistributionWithTrend(root, "vp", vp, dummy, false, errTxt, true) == true)
    use_vp = true;

  std::vector<DistributionWithTrendStorage *> bulk_modulus;
  if(parseDistributionWithTrend(root, "bulk-modulus", bulk_modulus, dummy, false, errTxt, true) == true)
    use_modulus = true;

  if(use_vp == true && use_modulus == true)
    errTxt += "Both <vp> and <bulk-modulus> can not be used in <fluid><tabulated>\n";
  else if(use_vp == false && use_modulus == false)
    errTxt += "One of <vp> or <bulk-modulus> must be used in <fluid><tabulated>\n";

  std::vector<DistributionWithTrendStorage *> density;
  if(parseDistributionWithTrend(root, "density", density, dummy, false, errTxt, true) == false)
    errTxt += "<density> needs to be specified in <fluid><tabulated>\n";

  std::vector<DistributionWithTrendStorage *> correlation_vp_density_with_trend;
  std::vector<double>                         correlation_vp_density;
  if(parseDistributionWithTrend(root, "correlation-vp-density", correlation_vp_density_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_vp_density_with_trend, "correlation", correlation_vp_density, errTxt);
    for(size_t i=0; i<correlation_vp_density.size(); i++) {
      if(correlation_vp_density[i] > 1 || correlation_vp_density[i] < -1)
        errTxt += "<correlation-vp-density> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_vp_density_with_trend[0]->GetIsShared() == false)
      delete correlation_vp_density_with_trend[0];
  }
  else
    correlation_vp_density.push_back(0.0);


  std::vector<DistributionWithTrendStorage *> correlation_bulk_density_with_trend;
  std::vector<double>                         correlation_bulk_density;
  if(parseDistributionWithTrend(root, "correlation-bulk-density", correlation_bulk_density_with_trend, dummy, false, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(correlation_bulk_density_with_trend, "correlation", correlation_bulk_density, errTxt);
    for(size_t i=0; i<correlation_bulk_density.size(); i++) {
      if(correlation_bulk_density[i] > 1 || correlation_bulk_density[i] < -1)
        errTxt += "<correlation-bulk-density> should be in the interval [-1,1] in the tabulated model\n";
    }
    if(correlation_bulk_density_with_trend[0]->GetIsShared() == false)
      delete correlation_bulk_density_with_trend[0];
  }
  else
    correlation_bulk_density.push_back(0.0);


  if(use_vp) {
    if(constituent == ModelSettings::FLUID) {
      DistributionsFluidStorage * fluid = new TabulatedVelocityFluidStorage(vp, density, correlation_vp_density);
      modelSettings_->addFluid(label, fluid);
    }
    else {
      errTxt += "Implementation error: parseTabulatedFluid() can only be used for fluids. Use parseTabulated()\n";
    }
  }
  else {
    if(constituent == ModelSettings::FLUID) {
      DistributionsFluidStorage * fluid = new TabulatedModulusFluidStorage(bulk_modulus, density, correlation_bulk_density);
      modelSettings_->addFluid(label, fluid);
    }
    else {
      errTxt += "Implementation error: parseTabulatedFluid() can only be used for fluids. Use parseTabulated()\n";
    }
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

void
XmlModelFile::FindDoubleValueFromDistributionWithTrend(const std::vector<DistributionWithTrendStorage *> & dist_with_trend,
                                                       std::string                                         type,
                                                       std::vector<double>                               & value,
                                                       std::string                                       & errTxt) const
{
  int n_vintages = static_cast<int>(dist_with_trend.size());

  value.resize(n_vintages);

  for(int i=0; i<n_vintages; i++) {
    if(typeid(*(dist_with_trend[i])) == typeid(DeltaDistributionWithTrendStorage)) {
      const DeltaDistributionWithTrendStorage * d1 = dynamic_cast<const DeltaDistributionWithTrendStorage *>(dist_with_trend[i]);
      if(typeid((*d1->GetMean())) == typeid(NRLib::TrendConstantStorage)) {
        const NRLib::TrendConstantStorage * t1 = dynamic_cast<const NRLib::TrendConstantStorage *>(d1->GetMean());
        value[i] = t1->GetMean();
      }
    }
    else {
      if(n_vintages == 1)
        errTxt += "All "+type+" variables need to be double values, not trends or distributions\n";
      else
        errTxt += "All "+type+" variables need to be double values in vintage "+NRLib::ToString(i+1)+", not trends or distributions\n";
    }
  }
}

bool
XmlModelFile::parseReservoir(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("reservoir");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("variable");

  std::string label;
  std::vector<DistributionWithTrendStorage *> distributionWithTrend; //Deleted in ~Modelsettings

  while(parseDistributionWithTrend(root, "variable", distributionWithTrend, label, true, errTxt) == true) {
    if(label == "")
      errTxt += "All reservoir variables need to be defined using <label>\n";

    modelSettings_->addReservoirVariable(label, distributionWithTrend);

    for(size_t i=0; i<distributionWithTrend.size(); i++)
      distributionWithTrend[i] = NULL;
    distributionWithTrend.clear();
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseEvolve(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("evolve");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("reservoir-variable");
  legalCommands.push_back("one-year-correlation");
  legalCommands.push_back("vintage");

  std::string variable;
  if(parseValue(root, "reservoir-variable", variable, errTxt) == false)
    errTxt += "The keyword <reservoir-variable> telling which variable to be evolved needs to be given in evlove.\n";

  typedef std::map<std::string, std::vector<DistributionWithTrendStorage *> > my_map;
  const my_map& reservoir_variable = modelSettings_->getReservoirVariable();
  my_map::const_iterator it = reservoir_variable.find(variable);

  std::vector<DistributionWithTrendStorage *> evolving_variable;

  if(it != reservoir_variable.end())
    evolving_variable = it->second;
  else {
    errTxt += "The variable "+variable+" used in <evolve> is not defined in <reservoir>.\n";
    errTxt += "  Note that <reservoir> needs to be given before <evolve> in <rock-physics>\n";
  }

  while(parseEvolveVintage(root, evolving_variable, errTxt) == true);

  double correlation;
  if(parseValue(root, "one-year-correlation", correlation, errTxt) == true) {
    if(correlation <= -1 || correlation >= 1)
      errTxt += "The <one-year-correlation> of the <reservoir-variable> should be in the interval (-1,1) in <evolve>\n";
  }
  else
    correlation = 1.0;


  for(size_t i=0; i<evolving_variable.size(); i++)
    evolving_variable[i]->SetOneYearCorrelation(correlation);

  modelSettings_->addReservoirVariable(variable, evolving_variable); //Replace the variable


  // Check consistency
  int evolve_size = static_cast<int>(evolving_variable.size());

  std::vector<int> vintage_year(evolve_size);
  for(int i=0; i<evolve_size; i++)
    vintage_year[i] = evolving_variable[i]->GetVintageYear();

  if(evolve_size > 1) {
    if(vintage_year[0] < 1)
      errTxt += "The vintage years need to be larger than zero in <evolve><vintage><vintage-year> in the rock physics model\n";

    int compare = vintage_year[0];
    for(int i=1; i<evolve_size; i++) {
      if(vintage_year[i] <= compare) {
        errTxt += "The vintage years need to be given in ascending order in <evolve><vintage><vintage-year> in the rock physics model\n";
        break;
      }
      else
        compare = vintage_year[i];
    }
  }

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}


bool
XmlModelFile::parseEvolveVintage(TiXmlNode * node, std::vector<DistributionWithTrendStorage *> & reservoir_variable, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("vintage");

  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("distribution");
  legalCommands.push_back("vintage-year");

  std::string dummy_label;
  if(parseDistributionWithTrend(root, "distribution", reservoir_variable, dummy_label, true, errTxt) == false)
    errTxt += "<distribution> needs to be given in <evolve><vintage>.\n";

  size_t variable_size = reservoir_variable.size();

  int vintage_year;
  if(parseValue(root, "vintage-year", vintage_year, errTxt) == false)
    errTxt += "A unique integer vintage year needs to be defined for each <evolve><vintage>\n";

  reservoir_variable[variable_size-1]->SetVintageYear(vintage_year);


  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}


bool
XmlModelFile::parseDistributionWithTrend(TiXmlNode                                   * node,
                                         const std::string                           & keyword,
                                         std::vector<DistributionWithTrendStorage *> & storage,
                                         std::string                                 & label,
                                         bool                                          is_shared,       //True for the variables in reservoir
                                         std::string                                 & errTxt,
                                         bool                                          allowDistribution)
{
  TiXmlNode * root = node->FirstChildElement(keyword);
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("label");
  legalCommands.push_back("reservoir-variable");
  legalCommands.push_back("value");
  legalCommands.push_back("trend-1d");
  legalCommands.push_back("trend-2d");
  if(allowDistribution == true) {
    legalCommands.push_back("gaussian");
    legalCommands.push_back("uniform");
    legalCommands.push_back("beta");
  }
  legalCommands.push_back("divide");

  label = "";
  parseValue(root, "label", label, errTxt);

  DistributionWithTrendStorage * distWithTrend = NULL;

  int trendGiven = 0;

  double value;
  if(root->FirstChildElement() == NULL) { //We have an explicit value
    parseValue(node, keyword, value, errTxt);
    distWithTrend = new DeltaDistributionWithTrendStorage(value, is_shared);
    storage.push_back(distWithTrend);
    trendGiven++;
    return(true);
  }

  if(parseValue(root, "value", value, errTxt) == true) {
    distWithTrend = new DeltaDistributionWithTrendStorage(value, is_shared);
    storage.push_back(distWithTrend);
    trendGiven++;
  }

  NRLib::TrendStorage * trend;
  if(parse1DTrend(root, "trend-1d", trend, errTxt) == true) {
    distWithTrend = new DeltaDistributionWithTrendStorage(trend, is_shared);
    storage.push_back(distWithTrend);
    trendGiven++;
  }

  if(parse2DTrend(root, "trend-2d", trend, errTxt) == true) {
    distWithTrend = new DeltaDistributionWithTrendStorage(trend, is_shared);
    storage.push_back(distWithTrend);
    trendGiven++;
  }

  if(allowDistribution == true) {
    if(parseGaussianWithTrend(root, storage, is_shared, errTxt) == true)
      trendGiven++;

    std::string beta;
    if(parseBetaWithTrend(root, storage, is_shared, errTxt) == true)
      trendGiven++;

    std::string uniform;
    if(parseValue(root, "uniform", uniform, errTxt) == true)
      errTxt += "The uniform distribution has not been implemented\n";
  }

  std::string variable;
  if(parseValue(root, "reservoir-variable", variable, errTxt) == true) {
    typedef std::map<std::string, std::vector<DistributionWithTrendStorage *> > my_map;
    const my_map& reservoir_variable = modelSettings_->getReservoirVariable();
    my_map::const_iterator it = reservoir_variable.find(variable);
    if(it != reservoir_variable.end()) {
      std::vector<DistributionWithTrendStorage *> store = it->second;
      for(size_t i=0; i<store.size(); i++)
        storage.push_back(store[i]);
      label = variable;
      trendGiven++;
    }
    else {
      errTxt += "The variable "+variable+" is not defined in <reservoir>.\n";
      errTxt += "  Note that <reservoir> needs to be given before <rock> and <predefinitions> in <rock-physics>\n";
    }
  }

  double divide;
  if(parseValue(root, "divide", divide, errTxt) == true)
    errTxt += "The <divide> functionality in 'DistributionWithTrend' has not been implemented by NR\n";

  if(trendGiven == 0)
    errTxt += "Need at least one definition for the variable in <"+keyword+">\n";
  else if(trendGiven > 1)
    errTxt += "There can only be one definition for the variable in <"+keyword+">\n";

  checkForJunk(root, errTxt, legalCommands, true);
  return(true);
}

bool
XmlModelFile::parseGaussianWithTrend(TiXmlNode                                   * node,
                                     std::vector<DistributionWithTrendStorage *> & storage,
                                     bool                                          is_shared,
                                     std::string                                 & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("gaussian");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("mean");
  legalCommands.push_back("variance");

  std::vector<DistributionWithTrendStorage *> mean_storage;
  std::string label;
  bool ok = true;
  if(parseDistributionWithTrend(root, "mean", mean_storage, label, is_shared, errTxt, false) == false) {
    errTxt += "Keyword <mean> is lacking"+lineColumnText(root)+" for Gaussian distribution.\n";
    ok = false;
  }

  std::vector<DistributionWithTrendStorage *> variance_storage;
  if(parseDistributionWithTrend(root, "variance", variance_storage, label, is_shared, errTxt, false) == false) {
    errTxt += "Keyword <variance> is lacking "+lineColumnText(root)+" for Gaussian distribution.\n";
    ok = false;
  }

  if (ok) {
    const NRLib::TrendStorage * mean     = mean_storage[0]->CloneMean();
    const NRLib::TrendStorage * variance = variance_storage[0]->CloneMean();

    delete mean_storage[0];
    delete variance_storage[0];

    DistributionWithTrendStorage * dist = new NormalDistributionWithTrendStorage(mean, variance, is_shared);

    delete mean;
    delete variance;

    storage.push_back(dist);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseBetaWithTrend(TiXmlNode                                   * node,
                                 std::vector<DistributionWithTrendStorage *> & storage,
                                 bool                                          is_shared,
                                 std::string                                 & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("beta");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("mean");
  legalCommands.push_back("variance");
  legalCommands.push_back("lower-limit");
  legalCommands.push_back("upper-limit");

  std::vector<DistributionWithTrendStorage *> mean_storage;
  std::string label;
  bool ok = true;

  if(parseDistributionWithTrend(root, "mean", mean_storage, label, is_shared, errTxt, false) == false) {
    errTxt += "Keyword <mean> is lacking"+lineColumnText(root)+" for beta distribution.\n";
    ok = false;
  }

  std::vector<DistributionWithTrendStorage *> variance_storage;
  if(parseDistributionWithTrend(root, "variance", variance_storage, label, is_shared, errTxt, false) == false) {
    errTxt += "Keyword <variance> is lacking "+lineColumnText(root)+" for Gaussian distribution.\n";
    ok = false;
  }

  std::vector<DistributionWithTrendStorage *> lower_limit_storage;
  std::vector<double>                         lower_limit;
  if(parseDistributionWithTrend(root, "lower-limit", lower_limit_storage, label, is_shared, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(lower_limit_storage, "lower-limit", lower_limit, errTxt);
    delete lower_limit_storage[0];
  }
  else
    lower_limit.push_back(0.0);

  std::vector<DistributionWithTrendStorage *> upper_limit_storage;
  std::vector<double>                         upper_limit;
  if(parseDistributionWithTrend(root, "upper-limit", upper_limit_storage, label, is_shared, errTxt, false) == true) {
    FindDoubleValueFromDistributionWithTrend(upper_limit_storage, "upper-limit", upper_limit, errTxt);
    delete upper_limit_storage[0];
  }
  else
    upper_limit.push_back(1.0);

  if (ok) {
    const NRLib::TrendStorage * mean     = mean_storage[0]    ->CloneMean();
    const NRLib::TrendStorage * variance = variance_storage[0]->CloneMean();

    delete mean_storage[0];
    delete variance_storage[0];

    DistributionWithTrendStorage * dist = new BetaDistributionWithTrendStorage(mean, variance, lower_limit[0], upper_limit[0], is_shared);

    delete mean;
    delete variance;

    storage.push_back(dist);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parse1DTrend(TiXmlNode * node, const std::string & keyword, NRLib::TrendStorage *& trend, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("trend-1d");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("reference-parameter");

  std::string reference_parameter;
  if(parseValue(root, "reference-parameter", reference_parameter, errTxt) == false)
    errTxt += "The <reference-parameter> for <"+keyword+"> needs to be given in <trend-1d>\n";

  std::string file_name;
  if(parseValue(root, "file-name", file_name, errTxt) == false)
    errTxt += "The <file-name> for <"+keyword+"> needs to be given in <trend-1d>\n";

  trend = new NRLib::Trend1DStorage(file_name,reference_parameter);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parse2DTrend(TiXmlNode * node, const std::string & keyword, NRLib::TrendStorage *& trend, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("trend-2d");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("reference-parameter-first-axis");
  legalCommands.push_back("reference-parameter-second-axis");

  std::string first_reference;
  if(parseValue(root, "reference-parameter-first-axis", first_reference, errTxt) == false)
    errTxt += "The <reference-parameter-first-axis> for <"+keyword+"> needs to be given in <trend-2d>\n";

  std::string second_reference;
  if(parseValue(root, "reference-parameter-second-axis", second_reference, errTxt) == false)
    errTxt += "The <reference-parameter-second-axis> for <"+keyword+"> needs to be given in <trend-2d>n";

  std::string file_name;
  if(parseValue(root, "file-name", file_name, errTxt) == false)
    errTxt += "The <file-name> for <"+keyword+"> needs to be given in <trend-2d>\n";

  trend = new NRLib::Trend2DStorage(file_name, first_reference, second_reference);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseTrendCube(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("trend-cube");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("parameter-name");
  legalCommands.push_back("file-name");

  std::string name;
  if(parseValue(root, "parameter-name", name, errTxt) == true)
    modelSettings_->addTrendCubeParameter(name);
  else
    errTxt += "<parameter-name> needs to be specified in <trend-cube> when <rock-physics> is used.\n";

  if(parseValue(root, "file-name", name, errTxt) == true)
    inputFiles_->addTrendCubes(name);
  else
    errTxt += "<file-name> needs to be specified in <trend-cube> when <rock-physics> is used.\n";

  checkForJunk(root, errTxt, legalCommands, true); //allow duplicates
  return(true);
}

bool
XmlModelFile::parseProjectSettings(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("project-settings");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("output-volume");
  legalCommands.push_back("time-to-depth-mapping-for-3d-wavelet");
  legalCommands.push_back("io-settings");
  legalCommands.push_back("advanced-settings");

  if(parseOutputVolume(root, errTxt) == false)
    errTxt += "Command <output-volume> is needed in command <"+
      root->ValueStr()+">"+lineColumnText(root)+".\n";

  if (parseTime3DMapping(root, errTxt))
    modelSettings_->setHasTime3DMapping(true);
  parseIOSettings(root, errTxt);
  parseAdvancedSettings(root, errTxt);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseOutputVolume(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("output-volume");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("interval-two-surfaces");
  legalCommands.push_back("interval-one-surface");
  legalCommands.push_back("area-from-surface");
  legalCommands.push_back("utm-coordinates");
  legalCommands.push_back("inline-crossline-numbers");

  bool interval = parseIntervalTwoSurfaces(root, errTxt);

  if(parseIntervalOneSurface(root, errTxt) == interval) { //Either both or none given
    if(interval == true)
      errTxt += "Time interval specified in more than one way in command <"+root->ValueStr()+"> "
        +lineColumnText(root)+".\n";
    else
      errTxt += "No time interval specified in command <"+root->ValueStr()+"> "
        +lineColumnText(root)+".\n";
  }

  bool area1 = parseAreaFromSurface(root, errTxt);
  bool area2 = parseILXLArea(root, errTxt);
  if(area1 && area2)
    errTxt+= "Area can not be given both by a surface and by inline and crossline numbers\n";
  bool area3 = parseUTMArea(root, errTxt);
  if(area1 && area3)
    errTxt+= "Area can no be given both by a surface and by coordinates\n";
  if(area2 && area3)
    errTxt+= "Area can not be given both by inline crossline numbers and UTM coordinates\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseIntervalTwoSurfaces(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("interval-two-surfaces");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-surface");
  legalCommands.push_back("base-surface");
  legalCommands.push_back("number-of-layers");
  legalCommands.push_back("velocity-field");
  legalCommands.push_back("velocity-field-from-inversion");

  modelSettings_->setParallelTimeSurfaces(false);

  if(parseTopSurface(root, errTxt) == false)
    errTxt += "Top surface not specified in command <"+root->ValueStr()+"> "
      +lineColumnText(root)+".\n";
  bool topDepthGiven;
  if(inputFiles_->getDepthSurfFile(0) == "")
    topDepthGiven = false;
  else
    topDepthGiven = true;

  if(parseBaseSurface(root, errTxt) == false)
    errTxt += "Base surface not specified in command <"+root->ValueStr()+"> "
      +lineColumnText(root)+".\n";
  bool baseDepthGiven;
  if(inputFiles_->getDepthSurfFile(1) == "")
    baseDepthGiven = false;
  else
    baseDepthGiven = true;

  int value = 0;
  if(parseValue(root, "number-of-layers", value, errTxt) == false)
    errTxt += "Number of layers not specified in command <"+root->ValueStr()+">"
      +lineColumnText(root)+".\n";
  else
  {
    if (value>0)
      modelSettings_->setTimeNz(value);
    else
      errTxt += "The number of layers needs to be larger than 0.\n";
  }

  std::string filename;
  bool externalField = parseFileName(root, "velocity-field", filename, errTxt);
  if(externalField == true)
    inputFiles_->setVelocityField(filename);

  bool inversionField = false;
  parseBool(root, "velocity-field-from-inversion", inversionField, errTxt);


  if(inversionField == true) {
    modelSettings_->setVelocityFromInversion(true);
    if(externalField == true)
      errTxt += "Both <velocity-field> and <velocity-field-from-inversion> given in command <"
        +root->ValueStr()+">"+lineColumnText(root)+".\n";
  }

  if(externalField == false && inversionField == false) {
    if(topDepthGiven != baseDepthGiven)
      errTxt += "Only one depth surface given, and no velocity field in command <"
        +root->ValueStr()+"> "+lineColumnText(root)+".\n";
    else if(topDepthGiven == true) // Both top and bottom given.
      modelSettings_->setDepthDataOk(true);
  }
  else {
    if(topDepthGiven == false && baseDepthGiven == false) {
      errTxt += "Velocity field, but no depth surface given in command <"
        +root->ValueStr()+"> "+lineColumnText(root)+".\n";
    }
    else
      modelSettings_->setDepthDataOk(true); //One or two surfaces, and velocity field.
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseTopSurface(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("top-surface");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("time-file");
  legalCommands.push_back("time-value");
  legalCommands.push_back("depth-file");

  std::string filename;
  bool timeFile = parseFileName(root,"time-file", filename, errTxt);
  if(timeFile == true)
    inputFiles_->addTimeSurfFile(filename);

  float value;
  bool timeValue = parseValue(root,"time-value", value, errTxt);
  if(timeValue == true) {
    if(timeFile == false)
      inputFiles_->addTimeSurfFile(NRLib::ToString(value));
    else
      errTxt += "Both file and value given for top time in command <"
        +root->ValueStr()+"> "+lineColumnText(root)+".\n";
  }
  else if(timeFile == false) {
    inputFiles_->addTimeSurfFile("");

    errTxt += "No time surface given in command <"+root->ValueStr()+"> "
      +lineColumnText(root)+".\n";
  }

  if(parseFileName(root,"depth-file", filename, errTxt) == true)
    inputFiles_->setDepthSurfFile(0, filename);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseBaseSurface(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("base-surface");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("time-file");
  legalCommands.push_back("time-value");
  legalCommands.push_back("depth-file");

  std::string filename;
  bool timeFile = parseFileName(root,"time-file", filename, errTxt);
  if(timeFile == true)
    inputFiles_->addTimeSurfFile(filename);

  float value;
  bool timeValue = parseValue(root,"time-value", value, errTxt);
  if(timeValue == true) {
    if(timeFile == false)
      inputFiles_->addTimeSurfFile(NRLib::ToString(value));
    else
      errTxt += "Both file and value given for base time in command <"+
        root->ValueStr()+"> "+lineColumnText(root)+".\n";
  }
  else if(timeFile == false) {
    inputFiles_->addTimeSurfFile("");
    errTxt += "No time surface given in command <"+root->ValueStr()+"> "
      +lineColumnText(root)+".\n";
  }

  if(parseFileName(root,"depth-file", filename, errTxt) == true)
    inputFiles_->setDepthSurfFile(1, filename);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseIntervalOneSurface(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("interval-one-surface");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("reference-surface");
  legalCommands.push_back("shift-to-interval-top");
  legalCommands.push_back("thickness");
  legalCommands.push_back("sample-density");

  modelSettings_->setParallelTimeSurfaces(true);

  std::string filename;
  if(parseFileName(root, "reference-surface", filename, errTxt) == true)
    inputFiles_->addTimeSurfFile(filename);
  else
    errTxt += "No reference surface given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double value;
  if(parseValue(root, "shift-to-interval-top", value, errTxt) == true)
    modelSettings_->setTimeDTop(value);

  if(parseValue(root, "thickness", value, errTxt) == true)
    modelSettings_->setTimeLz(value);
  else
    errTxt += "No thickness given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  if(parseValue(root, "sample-density", value, errTxt) == true)
    modelSettings_->setTimeDz(value);
  else
    errTxt += "No sample density given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseAreaFromSurface(TiXmlNode * node, std::string & errTxt)
{
 TiXmlNode * root = node->FirstChildElement("area-from-surface");
  if(root == 0)
    return(false);
  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("snap-to-seismic-data");

  std::string filename;
  if(parseFileName(root, "file-name", filename, errTxt) == true)
  {
    inputFiles_->setAreaSurfaceFile(filename);
  }
  bool snapToSeismicData = false;
  if (parseBool(root, "snap-to-seismic-data", snapToSeismicData, errTxt) == true) {
    modelSettings_->setSnapGridToSeismicData(true);
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseILXLArea(TiXmlNode * node, std::string & errTxt)
{
TiXmlNode * root = node->FirstChildElement("inline-crossline-numbers");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("il-start");
  legalCommands.push_back("il-end");
  legalCommands.push_back("xl-start");
  legalCommands.push_back("xl-end");
  legalCommands.push_back("il-step");
  legalCommands.push_back("xl-step");

  std::vector<int> ilxlnumbers(6);
  int value;
  if(parseValue(root, "il-start", value, errTxt) == true)
    ilxlnumbers[0] = value;
  else
    ilxlnumbers[0] = IMISSING;
  if(parseValue(root, "il-end", value, errTxt) == true)
  {
    if(value<ilxlnumbers[0])
    {
      errTxt+="il-end must be larger than il-start.";
    }
    ilxlnumbers[1] = value;
  }
  else
    ilxlnumbers[1] = IMISSING;
  if(parseValue(root, "xl-start", value, errTxt) == true)
    ilxlnumbers[2] = value;
  else
    ilxlnumbers[2] = IMISSING;
  if(parseValue(root, "xl-end", value, errTxt) == true)
  {
    if(value<ilxlnumbers[2])
    {
      errTxt+="xl-end must be larger than xl-start.";
    }
    ilxlnumbers[3] = value;
  }
  else
    ilxlnumbers[3] = IMISSING;
  if(parseValue(root, "il-step", value, errTxt) == true)
    ilxlnumbers[4] = value;
  else
    ilxlnumbers[4] = IMISSING;
  if(parseValue(root, "xl-step", value, errTxt) == true)
    ilxlnumbers[5] = value;
  else
    ilxlnumbers[5] = IMISSING;

  modelSettings_->setAreaILXLParameters(ilxlnumbers);
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseUTMArea(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("utm-coordinates");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("reference-point-x");
  legalCommands.push_back("reference-point-y");
  legalCommands.push_back("length-x");
  legalCommands.push_back("length-y");
  legalCommands.push_back("sample-density-x");
  legalCommands.push_back("sample-density-y");
  legalCommands.push_back("angle");
  legalCommands.push_back("snap-to-seismic-data");

  double x0 = 0;
  if(parseValue(root, "reference-point-x", x0, errTxt) == false)
    errTxt += "Reference x-coordinate must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double y0 = 0;
  if(parseValue(root, "reference-point-y", y0, errTxt) == false)
    errTxt += "Reference y-coordinate must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double lx = 0;
  if(parseValue(root, "length-x", lx, errTxt) == false)
    errTxt += "X-length must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double ly = 0;
  if(parseValue(root, "length-y", ly, errTxt) == false)
    errTxt += "Y-length must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  bool snapToSeismicData = parseBool(root, "snap-to-seismic-data", snapToSeismicData, errTxt);

  double dx = RMISSING;
  double dy = RMISSING;

  bool densX = parseValue(root, "sample-density-x", dx, errTxt);
  bool densY = parseValue(root, "sample-density-y", dy, errTxt);

  if (snapToSeismicData) {
    modelSettings_->setSnapGridToSeismicData(true);
    dx = lx;
    dy = ly;
    if(densX)
      TaskList::addTask("Keyword <sample-density-x> has no effect when <snap-to-seismic-data> has been specified.");
    if(densY)
      TaskList::addTask("Keyword <sample-density-y> has no effect when <snap-to-seismic-data> has been specified.");
  }
  else {
    if (!densX)
      errTxt += "Sample density for x must be given in command <"+ root->ValueStr()+"> "+lineColumnText(root)+".\n";
    if (!densY)
      errTxt += "Sample density for y must be given in command <"+ root->ValueStr()+"> "+lineColumnText(root)+".\n";
  }

  double angle = 0;
  if(parseValue(root, "angle", angle, errTxt) == false)
    errTxt += "Rotation angle must be given in command <"+root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double rot = (-1)*angle*(NRLib::Pi/180.0);
  int nx = static_cast<int>(lx/dx);
  int ny = static_cast<int>(ly/dy);
  SegyGeometry * geometry = new SegyGeometry(x0, y0, dx, dy, nx, ny, rot);
  modelSettings_->setAreaParameters(geometry);
  delete geometry;

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseTime3DMapping(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("time-to-depth-mapping-for-3d-wavelet");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands(3);
  legalCommands[0]="reference-depth";
  legalCommands[1]="average-velocity";
  legalCommands[2]="reference-time-surface";

  float value;
  if (parseValue(root, "reference-depth", value, errTxt) == false)
    errTxt += "Reference depth must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";
  else
    modelSettings_->setRefDepth(value);

  if (parseValue(root, "average-velocity", value, errTxt) == false)
    errTxt += "Average velocity must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";
  else
    modelSettings_->setAverageVelocity(value);

  std::string filename;

  if (parseFileName(root, "reference-time-surface", filename, errTxt) == false)
    errTxt += "Reference time surface must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";
  else
    inputFiles_->setRefSurfaceFile(filename);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseIOSettings(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("io-settings");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-directory");
  legalCommands.push_back("input-directory");
  legalCommands.push_back("output-directory");
  legalCommands.push_back("file-output-prefix");
  legalCommands.push_back("log-level");
  legalCommands.push_back("grid-output");
  legalCommands.push_back("well-output");
  legalCommands.push_back("wavelet-output");
  legalCommands.push_back("other-output");


  std::string topDir = IO::TopDirectory();
  parseValue(root, "top-directory", topDir, errTxt);
  ensureTrailingSlash(topDir);

  std::string inputDir = IO::InputDirectory();
  parseValue(root, "input-directory", inputDir, errTxt);
  inputDir = topDir+inputDir;
  ensureTrailingSlash(inputDir);
  inputFiles_->setInputDirectory(inputDir);

  std::string outputDir = IO::OutputDirectory();
  parseValue(root, "output-directory", outputDir, errTxt);
  outputDir = topDir+outputDir;
  ensureTrailingSlash(outputDir);
  IO::setOutputPath(outputDir);

  parseGridOutput(root, errTxt);
  parseWellOutput(root, errTxt);
  parseWaveletOutput(root, errTxt);
  parseOtherOutput(root, errTxt);
  std::string value;
  if(parseValue(root, "file-output-prefix", value, errTxt) == true) {
    IO::setFilePrefix(value);
  }

  std::string level;
  if(parseValue(root, "log-level", level, errTxt) == true) {
    int logLevel = LogKit::Error;
    if(level=="error")
      logLevel = LogKit::L_Error;
    else if(level=="warning")
      logLevel = LogKit::L_Warning;
    else if(level=="low")
      logLevel = LogKit::L_Low;
    else if(level=="medium")
      logLevel = LogKit::L_Medium;
    else if(level=="high")
      logLevel = LogKit::L_High;
    else if(level=="debuglow")
      logLevel = LogKit::L_DebugLow;
    else if(level=="debughigh")
      logLevel = LogKit::L_DebugHigh;
    else {
      errTxt += "Unknown log level " + level + " in command <log-level>. ";
      errTxt += "Choose from: error, warning, low, medium, and high\n";
      checkForJunk(root, errTxt, legalCommands);
      return(false);
    }
    modelSettings_->setLogLevel(logLevel);
  }

  //Test-open log file, to check valid path.
  std::string logFileName = IO::makeFullFileName("",IO::FileLog()+IO::SuffixTextFiles());
  std::ofstream file;
  try {
    NRLib::OpenWrite(file, logFileName);
  }
  catch(NRLib::Exception & e) {
    errTxt += std::string(e.what()) + "\n";
  }
  file.close();

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


void
XmlModelFile::ensureTrailingSlash(std::string & directory)
{
  std::string slash("/");

  if (directory.find_last_of(slash) != directory.length() - 1) // There are no slashes present
    directory += slash;
}



bool
XmlModelFile::parseGridOutput(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("grid-output");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("domain");
  legalCommands.push_back("format");
  legalCommands.push_back("elastic-parameters");
  legalCommands.push_back("seismic-data");
  legalCommands.push_back("other-parameters");

  parseGridFormats(root, errTxt);
  parseGridDomains(root, errTxt);
  parseGridElasticParameters(root, errTxt);
  parseGridSeismicData(root, errTxt);
  parseGridOtherParameters(root, errTxt);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseGridDomains(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("domain");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("depth");
  legalCommands.push_back("time");

  bool useDomain = false;
  int domainFlag = modelSettings_->getOutputGridDomain();
  if(parseBool(root, "time", useDomain, errTxt) == true) {
    if(useDomain == true)
      domainFlag = (domainFlag | IO::TIMEDOMAIN);
    else if((domainFlag & IO::TIMEDOMAIN) > 0)
      domainFlag -= IO::TIMEDOMAIN;
  }
  if(parseBool(root, "depth", useDomain, errTxt) == true) {
    if(useDomain == true) {
      domainFlag = (domainFlag | IO::DEPTHDOMAIN);
    }
    else if((domainFlag & IO::DEPTHDOMAIN) > 0)
      domainFlag -= IO::DEPTHDOMAIN;
  }
  modelSettings_->setOutputGridDomain(domainFlag);

  if(domainFlag == 0)
    errTxt += "Both time and depth domain output turned off after command <"
      +root->ValueStr()+"> "+lineColumnText(root)+".\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseGridFormats(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("format");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("segy-format");
  legalCommands.push_back("segy");
  legalCommands.push_back("storm");
  legalCommands.push_back("ascii");
  legalCommands.push_back("sgri");
  legalCommands.push_back("crava");
  TraceHeaderFormat *thf = NULL;
  bool segyFormat = parseTraceHeaderFormat(root, "segy-format",thf, errTxt);
  if(segyFormat==true)
    modelSettings_->setTraceHeaderFormatOutput(thf);
  bool useFormat = false;
  int formatFlag = 0;
  bool stormSpecified = false;  //Default format, error if turned off and no others turned on.
  if((parseBool(root, "segy", useFormat, errTxt) == true && useFormat == true) || segyFormat==true)
    formatFlag += IO::SEGY;
  if(parseBool(root, "storm", useFormat, errTxt) == true) {
    stormSpecified = true;
    if(useFormat == true)
      formatFlag += IO::STORM;
  }
  if(parseBool(root, "ascii", useFormat, errTxt) == true && useFormat == true)
    formatFlag += IO::ASCII;
  if(parseBool(root, "sgri", useFormat, errTxt) == true && useFormat == true)
    formatFlag += IO::SGRI;
  if(parseBool(root, "crava", useFormat, errTxt) == true && useFormat == true)
    formatFlag += IO::CRAVA;

  if(formatFlag > 0 || stormSpecified == true)
    modelSettings_->setOutputGridFormat(formatFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseGridElasticParameters(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("elastic-parameters");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("vp");
  legalCommands.push_back("vs");
  legalCommands.push_back("density");
  legalCommands.push_back("lame-lambda");
  legalCommands.push_back("lame-mu");
  legalCommands.push_back("poisson-ratio");
  legalCommands.push_back("ai");
  legalCommands.push_back("si");
  legalCommands.push_back("vp-vs-ratio");
  legalCommands.push_back("murho");
  legalCommands.push_back("lambdarho");
  legalCommands.push_back("background");
  legalCommands.push_back("background-trend");

  bool value     = false;
  int  paramFlag = 0;

  if(parseBool(root, "vp", value, errTxt) == true && value == true)
    paramFlag += IO::VP;
  if(parseBool(root, "vs", value, errTxt) == true && value == true)
    paramFlag += IO::VS;
  if(parseBool(root, "density", value, errTxt) == true && value == true)
    paramFlag += IO::RHO;
  if(parseBool(root, "lame-lambda", value, errTxt) == true && value == true)
    paramFlag += IO::LAMELAMBDA;
  if(parseBool(root, "lame-mu", value, errTxt) == true && value == true)
    paramFlag += IO::LAMEMU;
  if(parseBool(root, "poisson-ratio", value, errTxt) == true && value == true)
    paramFlag += IO::POISSONRATIO;
  if(parseBool(root, "ai", value, errTxt) == true && value == true)
    paramFlag += IO::AI;
  if(parseBool(root, "si", value, errTxt) == true && value == true)
    paramFlag += IO::SI;
  if(parseBool(root, "vp-vs-ratio", value, errTxt) == true && value == true)
    paramFlag += IO::VPVSRATIO;
  if(parseBool(root, "murho", value, errTxt) == true && value == true)
    paramFlag += IO::MURHO;
  if(parseBool(root, "lambdarho", value, errTxt) == true && value == true)
    paramFlag += IO::LAMBDARHO;
  if(parseBool(root, "background", value, errTxt) == true && value == true)
    paramFlag += IO::BACKGROUND;
  if(parseBool(root, "background-trend", value, errTxt) == true && value == true)
    paramFlag += IO::BACKGROUND_TREND;

  if (modelSettings_->getOutputGridsDefaultInd()) {
    if (paramFlag > 0) {
      modelSettings_->setOutputGridsElastic(paramFlag);
      modelSettings_->setOutputGridsDefaultInd(false);
    }
  }
  else {
    paramFlag = paramFlag | modelSettings_->getOutputGridsElastic();
    modelSettings_->setOutputGridsElastic(paramFlag);
  }
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseGridSeismicData(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("seismic-data");
  if(root == 0)
    return(false);
  std::vector<std::string> legalCommands;
  legalCommands.push_back("original");
  legalCommands.push_back("synthetic");
  legalCommands.push_back("residuals");
  legalCommands.push_back("synthetic-residuals");

  bool value    = false;
  int paramFlag = 0;

  if(parseBool(root, "residuals", value, errTxt) == true && value == true)
    paramFlag += IO::RESIDUAL;

  if(parseBool(root, "original", value, errTxt) == true && value == true)
    paramFlag += IO::ORIGINAL_SEISMIC_DATA;

  if(parseBool(root, "synthetic", value, errTxt) == true && value == true)
  {
    paramFlag += IO::SYNTHETIC_SEISMIC_DATA;
    modelSettings_->setGenerateSeismicAfterInv(true);
  }

  if(parseBool(root, "synthetic-residuals", value, errTxt) == true && value == true)
  {
    paramFlag += IO::SYNTHETIC_RESIDUAL;
    modelSettings_->setGenerateSeismicAfterInv(true);
  }

  if (modelSettings_->getOutputGridsDefaultInd() && paramFlag > 0){
    modelSettings_->setOutputGridsDefaultInd(false);
    modelSettings_->setOutputGridsElastic(0);
  }
  modelSettings_->setOutputGridsSeismic(paramFlag);
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseGridOtherParameters(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("other-parameters");
  if(root == 0)
    return(false);
  std::vector<std::string> legalCommands;
  legalCommands.push_back("facies-probabilities");
  legalCommands.push_back("facies-probabilities-with-undef");
  legalCommands.push_back("facies-likelihood");
  legalCommands.push_back("correlations");
  legalCommands.push_back("time-to-depth-velocity");
  legalCommands.push_back("extra-grids");
  legalCommands.push_back("seismic-quality-grid");

  bool facies           = true;
  bool faciesUndef      = false;
  bool faciesLH         = false;
  bool faciesValue      = false;
  bool faciesValueUndef = false;
  int  paramFlag        = 0;

  facies      = parseBool(root, "facies-probabilities", faciesValue, errTxt);
  faciesUndef = parseBool(root, "facies-probabilities-with-undef", faciesValueUndef, errTxt);
  faciesLH    = parseBool(root, "facies-likelihood", faciesValueUndef, errTxt);

  if (modelSettings_->getEstimateFaciesProb()){
    int tmpFlag = 0;
    if(faciesLH == true)
      tmpFlag += IO::FACIES_LIKELIHOOD;
    if(faciesUndef == true)
      tmpFlag += IO::FACIESPROB_WITH_UNDEF;
    if(facies == true || tmpFlag == 0)
      tmpFlag += IO::FACIESPROB;
    paramFlag += tmpFlag;
  }
  else{
    if (facies || faciesUndef || faciesLH){
      if (faciesValue || faciesValueUndef)
        errTxt += "Facies probability related cubes can not be specified under <other-output> when facies probabilities are not estimated.\n";
    }
  }

  bool value = false;
  if(parseBool(root, "correlations", value, errTxt) == true && value == true)
    paramFlag += IO::CORRELATION;
  if(parseBool(root, "extra-grids", value, errTxt) == true && value == true)
    paramFlag += IO::EXTRA_GRIDS;
  if(parseBool(root, "time-to-depth-velocity", value, errTxt) == true && value == true)
    paramFlag += IO::TIME_TO_DEPTH_VELOCITY;
  if(parseBool(root, "seismic-quality-grid", value, errTxt) == true && value == true)
    paramFlag += IO::SEISMIC_QUALITY_GRID;


  if (modelSettings_->getOutputGridsDefaultInd() && paramFlag > 0){
    modelSettings_->setOutputGridsDefaultInd(false);
    modelSettings_->setOutputGridsElastic(0);
  }

  modelSettings_->setOutputGridsOther(paramFlag);
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWellOutput(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("well-output");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("format");
  legalCommands.push_back("wells");
  legalCommands.push_back("blocked-wells");
  legalCommands.push_back("blocked-logs");

  parseWellFormats(root, errTxt);

  bool value;
  int wellFlag = 0;
  if(parseBool(root, "wells", value, errTxt) == true && value == true)
    wellFlag += IO::WELLS;
  if(parseBool(root, "blocked-wells", value, errTxt) == true && value == true)
    wellFlag += IO::BLOCKED_WELLS;
  if(parseBool(root, "blocked-logs", value, errTxt) == true && value == true)
    wellFlag += IO::BLOCKED_LOGS;

  modelSettings_->setWellOutputFlag(wellFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseWellFormats(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("format");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("rms");
  legalCommands.push_back("norsar");

  bool useFormat = false;
  int formatFlag = 0;
  bool rmsSpecified = false;  //Default format, check if turned off.
  if(parseBool(root, "rms", useFormat, errTxt) == true) {
    rmsSpecified = true;
    if(useFormat == true)
      formatFlag += IO::RMSWELL;
  }
  if(parseBool(root, "norsar", useFormat, errTxt) == true && useFormat == true)
    formatFlag += IO::NORSARWELL;

  if(formatFlag > 0 || rmsSpecified == true)
    modelSettings_->setWellFormatFlag(formatFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseWaveletOutput(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("wavelet-output");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("format");
  legalCommands.push_back("well-wavelets");
  legalCommands.push_back("global-wavelets");
  legalCommands.push_back("local-wavelets");

  parseWaveletFormats(root, errTxt);

  bool value;
  int waveletFlag = 0;
  if(parseBool(root, "well-wavelets", value, errTxt) == true && value == true)
    waveletFlag += IO::WELL_WAVELETS;
  if(parseBool(root, "global-wavelets", value, errTxt) == true && value == true)
    waveletFlag += IO::GLOBAL_WAVELETS;
  if(parseBool(root, "local-wavelets", value, errTxt) == true && value == true)
    waveletFlag += IO::LOCAL_WAVELETS;

  modelSettings_->setWaveletOutputFlag(waveletFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseWaveletFormats(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("format");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("jason");
  legalCommands.push_back("norsar");

  modelSettings_->setWaveletFormatManual(true);
  bool useFormat = false;
  int formatFlag = 0;
  bool jasonSpecified = false;  //Default format, check if turned off.
  if(parseBool(root, "jason", useFormat, errTxt) == true) {
    jasonSpecified = true;
    if(useFormat == true)
      formatFlag += IO::JASONWAVELET;
  }
  if(parseBool(root, "norsar", useFormat, errTxt) == true && useFormat == true)
    formatFlag += IO::NORSARWAVELET;

  if(formatFlag > 0 || jasonSpecified == true)
    modelSettings_->setWaveletFormatFlag(formatFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseOtherOutput(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("other-output");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("extra-surfaces");
  legalCommands.push_back("prior-correlations");
  legalCommands.push_back("background-trend-1d");
  legalCommands.push_back("local-noise");
  legalCommands.push_back("rock-physics-distributions");
  legalCommands.push_back("error-file");
  legalCommands.push_back("task-file");

  bool value;
  int otherFlag = 0;
  if(parseBool(root, "extra-surfaces", value, errTxt) == true && value == true)
    otherFlag += IO::EXTRA_SURFACES;
  if(parseBool(root, "prior-correlations", value, errTxt) == true && value == true)
    otherFlag += IO::PRIORCORRELATIONS;
  if(parseBool(root, "background-trend-1d", value, errTxt) == true && value == true)
    otherFlag += IO::BACKGROUND_TREND_1D;
  if(parseBool(root, "local-noise", value, errTxt) == true && value == true)
    otherFlag += IO::LOCAL_NOISE;
  if(parseBool(root, "rock-physics-distributions", value, errTxt) == true && value == true)
    otherFlag += IO::ROCK_PHYSICS;
  if(parseBool(root, "error-file", value, errTxt) == true && value == true)
    otherFlag += IO::ERROR_FILE;
  if(parseBool(root, "task-file", value, errTxt) == true && value == true)
    otherFlag += IO::TASK_FILE;

  modelSettings_->setOtherOutputFlag(otherFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseAdvancedSettings(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("advanced-settings");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("fft-grid-padding");
  legalCommands.push_back("vp-vs-ratio");
  legalCommands.push_back("vp-vs-ratio-from-wells");
  legalCommands.push_back("use-intermediate-disk-storage");
  legalCommands.push_back("maximum-relative-thickness-difference");
  legalCommands.push_back("frequency-band");
  legalCommands.push_back("energy-threshold");
  legalCommands.push_back("wavelet-tapering-length");
  legalCommands.push_back("minimum-relative-wavelet-amplitude");
  legalCommands.push_back("maximum-wavelet-shift");
  legalCommands.push_back("minimum-sampling-density");
  legalCommands.push_back("minimum-horizontal-resolution");
  legalCommands.push_back("white-noise-component");
  legalCommands.push_back("reflection-matrix");
  legalCommands.push_back("kriging-data-limit");
  legalCommands.push_back("debug-level");
  legalCommands.push_back("smooth-kriged-parameters");
  legalCommands.push_back("rms-panel-mode");
  legalCommands.push_back("guard-zone");
  legalCommands.push_back("3d-wavelet-tuning-factor");
  legalCommands.push_back("gradient-smoothing-range");
  legalCommands.push_back("estimate-well-gradient-from-seismic");

  parseFFTGridPadding(root, errTxt);

  float ratio = RMISSING;
  if(parseValue(root,"vp-vs-ratio", ratio, errTxt) == true)
    modelSettings_->setVpVsRatio(ratio);

  bool ratio_from_wells = false;
  if(parseBool(root,"vp-vs-ratio-from-wells", ratio_from_wells, errTxt) == true)
    modelSettings_->setVpVsRatioFromWells(ratio_from_wells);

  if (ratio_from_wells && ratio != RMISSING) {
    errTxt += "You cannot both specify a Vp/Vs ratio (" + NRLib::ToString(ratio,2)
              + ") and ask the ratio to be estimated from well data.\n";
  }

  bool fileGrid;
  if(parseBool(root, "use-intermediate-disk-storage", fileGrid, errTxt) == true)
    modelSettings_->setFileGrid(fileGrid);

  double limit;
  if(parseValue(root,"maximum-relative-thickness-difference", limit, errTxt) == true)
    modelSettings_->setLzLimit(limit);

  parseFrequencyBand(root, errTxt);

  float value = 0.0f;
  if(parseValue(root, "energy-threshold", value, errTxt) == true)
    modelSettings_->setEnergyThreshold(value);
  if(parseValue(root, "wavelet-tapering-length", value, errTxt) == true)
    modelSettings_->setWaveletTaperingL(value);
  if(parseValue(root, "minimum-relative-wavelet-amplitude", value, errTxt) == true)
    modelSettings_->setMinRelWaveletAmp(value);
  if(parseValue(root, "maximum-wavelet-shift", value, errTxt) == true)
    modelSettings_->setMaxWaveletShift(value);
  if(parseValue(root, "minimum-sampling-density", value, errTxt) == true)
    modelSettings_->setMinSamplingDensity(value);
  if(parseValue(root, "minimum-horizontal-resolution", value, errTxt) == true)
    modelSettings_->setMinHorizontalRes(value);
  if(parseValue(root, "white-noise-component", value, errTxt) == true)
    modelSettings_->setWNC(value);
  std::string filename;
  if(parseFileName(root, "reflection-matrix", filename, errTxt) == true)
    inputFiles_->setReflMatrFile(filename);
  int kLimit = 0;
  if(parseValue(root, "kriging-data-limit", kLimit, errTxt) == true) {
    if(modelSettings_->getKrigingParameter() >= 0)
      modelSettings_->setKrigingParameter(kLimit);
    else
      errTxt += "The number of data in neighbourhood when doing kriging must be larger than or equal to zero\n";
  }
  int level = 0;
  if(parseValue(root, "debug-level", level, errTxt) == true)
    modelSettings_->setDebugFlag(level);

  bool smooth = false;
  if(parseBool(root, "smooth-kriged-parameters", smooth, errTxt) == true)
    modelSettings_->setDoSmoothKriging(smooth);

  bool panel = false;
  if(parseBool(root, "rms-panel-mode", panel, errTxt) == true)
    modelSettings_->setRunFromPanel(panel);

  if(parseValue(root, "guard-zone", value, errTxt) == true) {
    float smooth_length = modelSettings_->getSmoothLength();
    smooth_length = std::min(smooth_length, value);
    modelSettings_->setGuardZone(value);
    modelSettings_->setSmoothLength(smooth_length);
  }
  if(parseValue(root, "3d-wavelet-tuning-factor", value, errTxt) == true) {
    modelSettings_->setWavelet3DTuningFactor(value);
    if (value < 1.0f || value > 500.0f)
      errTxt += "The 3D wavelet tuning factor must be in range [1.0, 500.0]\n";
  }
  if(parseValue(root, "gradient-smoothing-range", value, errTxt) == true) {
    modelSettings_->setGradientSmoothingRange(value);
    if (value < 25.0f || value > 500.0f)
      errTxt += "The gradient smoothing range for 3D wavelet estimation/inversion must be in range [25.0, 500.0]\n";
  }
  bool estimate = false;
  if(parseBool(root, "estimate-well-gradient-from-seismic", estimate, errTxt) == true)
    modelSettings_->setEstimateWellGradientFromSeismic(estimate);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseFFTGridPadding(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("fft-grid-padding");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("x-fraction");
  legalCommands.push_back("y-fraction");
  legalCommands.push_back("z-fraction");

  double padding;
  bool estXpad = true;
  bool estYpad = true;
  if(parseValue(root, "x-fraction", padding, errTxt) == true) {
    modelSettings_->setXPadFac(padding);
    estXpad = false;
  }
  if(parseValue(root, "y-fraction", padding, errTxt) == true) {
    modelSettings_->setYPadFac(padding);
    estYpad = false;
  }
  if(parseValue(root, "z-fraction", padding, errTxt) == true) {
    modelSettings_->setZPadFac(padding);
    modelSettings_->setEstimateZPadding(false);
  }

  if ((estXpad && !estYpad) || (!estXpad && estYpad))
    errTxt += "The lateral padding factors must either both be estimated or both given.";
  if (!estXpad && !estYpad)
    modelSettings_->setEstimateXYPadding(false);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseFrequencyBand(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("frequency-band");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("low-cut");
  legalCommands.push_back("high-cut");

  float value;
  if(parseValue(root, "low-cut", value, errTxt) == true)
    modelSettings_->setLowCut(value);
  if(parseValue(root, "high-cut", value, errTxt) == true)
    modelSettings_->setHighCut(value);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseTraceHeaderFormat(TiXmlNode * node, const std::string & keyword, TraceHeaderFormat *& thf, std::string & errTxt)
{
  thf = NULL;
  TiXmlNode * root = node->FirstChildElement(keyword);
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("standard-format");
  legalCommands.push_back("location-x");
  legalCommands.push_back("location-y");
  legalCommands.push_back("location-il");
  legalCommands.push_back("location-xl");
  legalCommands.push_back("location-scaling-coefficient");
  legalCommands.push_back("bypass-coordinate-scaling");

  std::string stdFormat;
  if(parseValue(root, "standard-format", stdFormat, errTxt) == true) {
    if(NRLib::Uppercase(stdFormat) == NRLib::Uppercase("seisworks"))
      thf = new TraceHeaderFormat(TraceHeaderFormat::SEISWORKS);
    else if(NRLib::Uppercase(stdFormat) == NRLib::Uppercase("iesx"))
      thf = new TraceHeaderFormat(TraceHeaderFormat::IESX);
    else if(NRLib::Uppercase(stdFormat) == NRLib::Uppercase("SIP"))
      thf = new TraceHeaderFormat(TraceHeaderFormat::SIP);
    else if(NRLib::Uppercase(stdFormat) == NRLib::Uppercase("charisma"))
      thf = new TraceHeaderFormat(TraceHeaderFormat::CHARISMA);
    else {
      errTxt += "Unknown segy-format '"+stdFormat+"' found on line"+
        NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n Known formats are 'seisworks', 'iesx', 'SIP', 'charisma'.\n";
      thf = new TraceHeaderFormat(TraceHeaderFormat::SEISWORKS);
    }
  }
  else
    thf = new TraceHeaderFormat(TraceHeaderFormat::SEISWORKS);

  int value;
  if(parseValue(root,"location-x",value, errTxt) == true)
    thf->SetUtmxLoc(value);
  if(parseValue(root,"location-y",value, errTxt) == true)
    thf->SetUtmyLoc(value);
  if(parseValue(root,"location-il",value, errTxt) == true)
    thf->SetInlineLoc(value);
  if(parseValue(root,"location-xl",value, errTxt) == true)
    thf->SetCrosslineLoc(value);
  if(parseValue(root,"location-scaling-coefficient",value, errTxt) == true)
    thf->SetScaleCoLoc(value);

  bool bypass;
  if(parseBool(root,"bypass-coordinate-scaling", bypass, errTxt) == true)
    if(bypass == true)
      thf->SetScaleCoLoc(-1);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseVariogram(TiXmlNode * node, const std::string & keyword, Vario * & vario, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement(keyword);
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("range");
  legalCommands.push_back("subrange");
  legalCommands.push_back("angle");
  legalCommands.push_back("power");
  legalCommands.push_back("variogram-type");

  float value, range = 1;
  float subrange=1;
  float angle = 0;
  float expo = 1;
  if(parseValue(root,"range", value, errTxt) == true)
    range = value;
  else
    errTxt += "Keyword <range> is lacking for variogram under command <"+keyword+
      ">, line "+NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
  if(parseValue(root,"subrange", value, errTxt) == true)
    subrange = value;
  if(parseValue(root,"angle", value, errTxt) == true)
    angle = 90.0f-value; //From geological to mathematical.
  bool power = parseValue(root,"power", value, errTxt);
  if(power == true)
    expo = value;

  if (range==0.0) {
    errTxt += "The value of <range> given for command <"+keyword+"> must be greater than zero.\n";
  }
  if (subrange==0.0) {
    errTxt += "The value of <subrange> given for command <"+keyword+"> must be greater than zero.\n";
  }

  std::string vType;
  vario = NULL;
  if(parseValue(root,"variogram-type", vType, errTxt) == false)
    errTxt += "Keyword <variogram-type> is lacking for variogram under command <"+keyword+
      ">, line "+NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
  else {
    if(vType == "genexp") {
      if(power == false) {
        errTxt += "Keyword <power> is lacking for gen. exp. variogram under command <"+keyword+
          ">, line "+NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
      }
      vario = new GenExpVario(expo, range, subrange, angle);
    }
    else if(vType == "spherical") {
      if(power == true) {
        errTxt += "Keyword <power> is given for spherical variogram under command <"+keyword+
          ">, line "+NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
      }
      vario = new SphericalVario(range, subrange, angle);
    }
    else
      errTxt += "Variogram type "+vType+" is not a valid variogram type. Please choose between genexp and spherical.\n";
  }

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseBool(TiXmlNode * node, const std::string & keyword, bool & value, std::string & errTxt, bool allowDuplicates)
{
  std::string tmpVal;
  std::string tmpErr = "";
  if(parseValue(node, keyword, tmpVal, tmpErr, allowDuplicates) == false)
    return(false);

  //Keyword is found.
  if(tmpErr == "") {
    if(tmpVal == "yes")
      value = true;
    else if(tmpVal == "no")
      value = false;
    else {
      tmpErr = "Found '"+tmpVal+"' under keyword <"+keyword+">, expected 'yes' or 'no'. This happened in command <"+
        node->ValueStr()+"> on line "+NRLib::ToString(node->Row())+", column "+NRLib::ToString(node->Column())+".\n";
    }
  }

  //No junk-clearing call, done in parseValue.
  errTxt += tmpErr;
  return(true);
}


bool
XmlModelFile::parseFileName(TiXmlNode * node, const std::string & keyword, std::string & filename, std::string & errTxt, bool allowDuplicates)
{
  filename = "";
  std::string value;
  std::string tmpErr;
  if(parseValue(node, keyword, value, tmpErr, allowDuplicates) == false)
    return(false);

  filename = value;

  //No junk-clearing, done in parseValue
  errTxt += tmpErr;
  return(true);
}


void
XmlModelFile::checkForJunk(TiXmlNode * root, std::string & errTxt, const std::vector<std::string> & legalCommands,
                    bool allowDuplicates)
{
  TiXmlNode * child = root->FirstChild();
  unsigned int startLength = static_cast<unsigned int>(errTxt.size());
  while(child != NULL) {
    switch(child->Type()) {
      case TiXmlNode::COMMENT :
      case TiXmlNode::DECLARATION :
        break;
      case TiXmlNode::TEXT :
        errTxt = errTxt + "Unexpected value '"+child->Value()+"' is not part of command <"+root->Value()+
          "> on line "+NRLib::ToString(child->Row())+", column "+NRLib::ToString(child->Column())+".\n";
        break;
      case TiXmlNode::ELEMENT :
        errTxt = errTxt + "Unexpected command <"+child->Value()+"> is not part of command <"+root->Value()+
          "> on line "+NRLib::ToString(child->Row())+", column "+NRLib::ToString(child->Column())+".\n";
        break;
      default :
        errTxt = errTxt + "Unexpected text '"+child->Value()+"' is not part of command <"+root->Value()+
          "> on line "+NRLib::ToString(child->Row())+", column "+NRLib::ToString(child->Column())+".\n";
        break;
    }
    root->RemoveChild(child);
    child = root->FirstChild();
  }

  if(startLength<errTxt.size()){
    errTxt = errTxt + "Legal commands are:\n";
    for(unsigned int i=0;i<legalCommands.size();i++)
      errTxt = errTxt +"  <"+legalCommands[i]+">\n";
  }

  TiXmlNode * parent = root->Parent();

  if(parent != NULL) {
    std::string cmd = root->ValueStr();
    parent->RemoveChild(root);

    if(allowDuplicates == false) {
      root = parent->FirstChildElement(cmd);
      if(root != NULL) {
        int n = 0;
        while(root != NULL) {
          n++;
          parent->RemoveChild(root);
          root = parent->FirstChildElement(cmd);
        }
        errTxt += "Found "+NRLib::ToString(n)+" extra occurences of command <"+cmd+"> under command <"+parent->Value()+
          "> on line "+NRLib::ToString(parent->Row())+", column "+NRLib::ToString(parent->Column())+".\n";
      }
    }
  }
}


std::string
XmlModelFile::lineColumnText(TiXmlNode * node)
{
  std::string result = " on line "+NRLib::ToString(node->Row())+", column "+NRLib::ToString(node->Column());
  return(result);
}


void
XmlModelFile::setDerivedParameters(std::string & errTxt)
{
  int areaSpecification;
  if(modelSettings_->getAreaParameters() != NULL) {
    areaSpecification = ModelSettings::AREA_FROM_UTM;
    for(int i=0; i<modelSettings_->getNumberOfTimeLapses(); i++){
      if (modelSettings_->getNoSeismicNeeded() && inputFiles_->getNumberOfSeismicFiles(i)>0)
        errTxt += "Seismic data should not be given when estimating background or correlations. \nExceptions are for optimization of well locations or if the area is taken from the seismic data.\n";
    }
    if (modelSettings_->getSnapGridToSeismicData()) {
      areaSpecification = ModelSettings::AREA_FROM_GRID_DATA_AND_UTM;
    }
  }
  else if(inputFiles_->getAreaSurfaceFile() != "") {
    areaSpecification = ModelSettings::AREA_FROM_SURFACE;
    for(int i=0; i<modelSettings_->getNumberOfTimeLapses(); i++){
      if (modelSettings_->getNoSeismicNeeded() && inputFiles_->getNumberOfSeismicFiles(i)>0)
        errTxt += "Seismic data should not be given when estimating background or correlations. \nExceptions are for optimization of well locations or if the area is taken from the seismic data.";
    }
    if (modelSettings_->getSnapGridToSeismicData()) {
      areaSpecification = ModelSettings::AREA_FROM_GRID_DATA_AND_SURFACE;
    }
  }
  else {
    areaSpecification = ModelSettings::AREA_FROM_GRID_DATA; // inversion:seismic data, forward modelling: Vp
    for(int i=0; i<modelSettings_->getNumberOfTimeLapses(); i++){
      if (modelSettings_->getNoSeismicNeeded() && inputFiles_->getNumberOfSeismicFiles(i)==0)
        errTxt += "The area needs to be defined from seismic data, a surface or UTM-coordinates.\n";
    }
  }
  modelSettings_->setAreaSpecification(areaSpecification);

  if (modelSettings_->getEstimateFaciesProb()) {
    // Cannot be placed under parseWellData() since we do not there if getUseFilterForFaciesProb() has been set.
    bool useFilter = modelSettings_->getUseFilterForFaciesProb();
    for (int i=0 ; i < modelSettings_->getNumberOfWells() ; i++) {
      int filterElasticLogs       = modelSettings_->getIndicatorFilter(i);
      int useForFaciesProbability = modelSettings_->getIndicatorFacies(i);
      if (useFilter && useForFaciesProbability != ModelSettings::NO && filterElasticLogs == ModelSettings::NO) {
        modelSettings_->setIndicatorFilter(i, ModelSettings::YES);
      }
    }
  }
}


void
XmlModelFile::checkConsistency(std::string & errTxt) {
  if(modelSettings_->getForwardModeling() == true)
    checkForwardConsistency(errTxt);
  else {
    if (modelSettings_->getEstimationMode())
      checkEstimationConsistency(errTxt);
    else
      checkInversionConsistency(errTxt);
  }
  if(modelSettings_->getLocalWaveletVario()==NULL)
      modelSettings_->copyBackgroundVarioToLocalWaveletVario();
  if(modelSettings_->getOptimizeWellLocation()==true)
    checkAngleConsistency(errTxt);
  checkIOConsistency(errTxt);
  if(modelSettings_->getDo4DInversion() && surveyFailed_ == false)
    checkTimeLapseConsistency(errTxt);

  if (inputFiles_->getReflMatrFile() != "") {
    if (modelSettings_->getVpVsRatio() != RMISSING) {
      errTxt += "You cannot specify a Vp/Vs ratio when a reflection matrix is read from file";
    }
    else if (modelSettings_->getVpVsRatioFromWells()) {
      errTxt += "You cannot ask the Vp/Vs ratio to be calculated from well data when a reflection matrix is read from file";
    }
  }
  if (modelSettings_->getVpVsRatio() != RMISSING) {
    double vpvs    = modelSettings_->getVpVsRatio();
    double vpvsMin = modelSettings_->getVpVsRatioMin();
    double vpvsMax = modelSettings_->getVpVsRatioMax();
    if (vpvs < vpvsMin) {
      errTxt += "Specified Vp/Vs of "+NRLib::ToString(vpvs,2)
                +" is less than minimum allowed value of "+NRLib::ToString(vpvsMin,2);
    }
    if (vpvs > vpvsMax) {
      errTxt += "Specified Vp/Vs of "+NRLib::ToString(vpvs,2)
                +" is larger than maximum allowed value of "+NRLib::ToString(vpvsMax,2);
    }
  }
}

void
XmlModelFile::checkForwardConsistency(std::string & errTxt)
{

  if (modelSettings_->getDo4DInversion())
    errTxt += "Forward modeling can not be done in 4D.\n";

  if (modelSettings_->getNumberOfAngles(0) == 0) {
    errTxt+="At least one wavelet must be specified for forward modelling\n";
  }

  if (modelSettings_->getNumberOfTraceHeaderFormats(0) == 0)
    modelSettings_->addTraceHeaderFormat(NULL);

  std::vector<float>  angle = modelSettings_->getAngle(0);
  std::vector<bool> wavelet = modelSettings_->getEstimateWavelet(0);

  for(int i=0; i<modelSettings_->getNumberOfAngles(0); i++){
    modelSettings_->setEstimateSNRatio(0,i,false);
    if(wavelet[i] == true)
      errTxt+="Wavelet must be given when doing forward modeling. Wavelet is not given for angle "+NRLib::ToString(angle[i]*(180/NRLib::Pi),1)+".\n";

    if(inputFiles_->getSeismicFile(0,i)!="")
      errTxt+="Seismic data should not be given when doing forward modeling.\n";
  }

  if(modelSettings_->getUseLocalNoise(0)==true)
    errTxt+="Local noise can not be used in forward modeling.\n";

  if(modelSettings_->getUseLocalWavelet()==true)
    errTxt+="Local wavelet can not be used in forward modeling.\n";

  if (modelSettings_->getNumberOfWells() > 0)
    errTxt +="Wells should not be given when doing forward modeling.\n";

  if (modelSettings_->getBackgroundType() == "background")
    errTxt += "An earth model needs to be given when doing forward modeling. The background model should not be given.\n";
  else if (modelSettings_->getBackgroundType() == "")
    errTxt += "An earth model needs to be given when doing forward modeling.\n";

  if (modelSettings_->getEstimate3DWavelet() && !modelSettings_->getHasTime3DMapping())
    errTxt += "Time 3D mapping must be given when 3D wavelet is to be estimated.\n";
}

void
XmlModelFile::checkEstimationConsistency(std::string & errTxt) {
  if (modelSettings_->getNumberOfWells()==0) {
    if (modelSettings_->getEstimateBackground())
      errTxt += "Wells are needed for estimation of background.\n";
    if (modelSettings_->getEstimateCorrelations())
      errTxt += "Wells are needed for estimation of correlations.\n";
    if (modelSettings_->getEstimateWaveletNoise())
      errTxt += "Wells are needed for estimation of wavelet and noise.\n";
  }
  if (modelSettings_->getEstimateWaveletNoise()==false && modelSettings_->getOptimizeWellLocation()==false)
    modelSettings_->setNoSeismicNeeded(true);
  else {
    if (inputFiles_->getNumberOfSeismicFiles(0)==0) {
      if (modelSettings_->getOptimizeWellLocation())
        errTxt += "Seismic data are needed for optimizing well locations.\n";
      if (modelSettings_->getEstimateWaveletNoise())
        errTxt += "Seismic data are needed for estimation of wavelets and noise.\n";
    }
  }
  if (modelSettings_->getEstimate3DWavelet() && !modelSettings_->getHasTime3DMapping())
    errTxt += "Time 3D mapping must be given when 3D wavelet is to be estimated.\n";
}

void
XmlModelFile::checkInversionConsistency(std::string & errTxt) {
  if (inputFiles_->getNumberOfSeismicFiles(0)==0)
    errTxt += "Seismic data are needed for inversion.\n";

  if (modelSettings_->getNumberOfWells() == 0) {
    std::vector<bool> useRicker = modelSettings_->getUseRickerWavelet(0);
    if ((inputFiles_->getBackFile(0)!="" || modelSettings_->getGenerateBackground() == false) &&
        (inputFiles_->getWaveletFile(0,0)!="" || useRicker[0] == true) &&
        inputFiles_->getTempCorrFile()!="" &&
        inputFiles_->getParamCorrFile()!="" ||
        modelSettings_->getGenerateBackgroundFromRockPhysics())
      modelSettings_->setNoWellNeeded(true);
    else {
      errTxt += "Wells are needed for the inversion. ";
      errTxt += "Alternatively, all of background, wavelet, \ntemporal correlation and parameter correlation must be given.\n";
    }
    if (modelSettings_->getKrigingParameter()>0)
      errTxt += "Wells are needed for kriging.\n";
    if (modelSettings_->getEstimateFaciesProb() && !modelSettings_->getFaciesProbFromRockPhysics())
      errTxt += "Wells are needed for facies probabilities.\n";
  }

  if (modelSettings_->getOutputGridsElastic() == 0 &&
      modelSettings_->getOutputGridsSeismic() == 0 &&
      modelSettings_->getOutputGridsOther()   == 0 &&
      modelSettings_->getWellOutputFlag()     == 0 &&
      modelSettings_->getWaveletOutputFlag()  == 0 &&
      modelSettings_->getOtherOutputFlag()    == 0)
    errTxt += "No output is specified for the inversion model.\n";

  if (modelSettings_->getEstimateFaciesProb()                                 &&
     (modelSettings_->getOutputGridsOther() & IO::FACIESPROB)            == 0 &&
     (modelSettings_->getOutputGridsOther() & IO::FACIESPROB_WITH_UNDEF) == 0 &&
     (modelSettings_->getWellOutputFlag()   & IO::BLOCKED_WELLS)         == 0)
  {
    errTxt += "Grid output for facies probabilities or facies probabilities with undefined value,\n";
    errTxt += "or blocked wells needs to be specified when doing facies estimation.\n";
  }
  if (modelSettings_->getEstimateFaciesProb() == false && modelSettings_->getFaciesProbRelative() == false)
    errTxt += "Absolute facies probabilities can not be requested without requesting facies probabilities under inversion settings.\n";
  if (modelSettings_->getEstimateFaciesProb() == false && (modelSettings_->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID))
    errTxt += "Seismic quality grid can not be estimated without requesting facies probabilities under inversion settings.\n";
  if(modelSettings_->getFaciesProbFromRockPhysics() == true  && (modelSettings_->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID))
    errTxt += "Seismic quality grid can not be estimated when facies probabilities are calculated using rock physics models\n";

  //Rock physics consistency
  if(modelSettings_->getFaciesProbFromRockPhysics()) {

    if(modelSettings_->getEstimateFaciesProb() == false && modelSettings_->getNumberOfVintages() == 0)
      errTxt += "Rocks in the rock physics prior model should not be given \nwithout requesting facies probabilities under inversion settings.\n";

    const std::map<std::string, DistributionsRockStorage *>& rock_storage = modelSettings_->getRockStorage();

    // Compare names in rock physics model with names given in .xml-file
    if(modelSettings_->getIsPriorFaciesProbGiven() == ModelSettings::FACIES_FROM_MODEL_FILE) {
      std::map<std::string, float> facies_probabilities = modelSettings_->getPriorFaciesProb();

      for(std::map<std::string, float>::const_iterator it = facies_probabilities.begin(); it != facies_probabilities.end(); it++) {
        std::map<std::string, DistributionsRockStorage *>::const_iterator iter = rock_storage.find(it->first);

        if(iter == rock_storage.end())
          errTxt += "Problem with rock physics prior model. Facies '"+it->first+"' is not one of the rocks given in the rock physics model.\n";
      }
    }

    // Compare names in rock physics model with names given as input in proability cubes
    else if(modelSettings_->getIsPriorFaciesProbGiven() == ModelSettings::FACIES_FROM_CUBES) {
      const std::map<std::string,std::string>& facies_probabilities = inputFiles_->getPriorFaciesProbFile();

      for(std::map<std::string, std::string>::const_iterator it = facies_probabilities.begin(); it != facies_probabilities.end(); it++) {
        std::map<std::string, DistributionsRockStorage *>::const_iterator iter = rock_storage.find(it->first);
        if (iter == rock_storage.end())
          errTxt += "Problem with rock physics prior model. Facies "+it->first+" is not one of the rocks given in the rock physics model.\n";
      }
    }
    else
      errTxt += "Prior facies probabilities must be given in the prior model when rock physics models are used\n";
  }

  if(modelSettings_->getGenerateBackgroundFromRockPhysics()) {
    // In case of setting background/prior model from rock physics, the file inputFiles_->getParamCorrFile() should not be set.
    // This assumption is relevant for the function ModelGeneral::processPriorCorrelations.
    if(inputFiles_->getParamCorrFile()!=""){
      errTxt += "Parameter correlation should not be specified in file when background/prior model is build from rock physics.";
    }

    if(modelSettings_->getGenerateBackground() == false) {
      errTxt += "The background model can not be given when rock physics models are used.\n";
      errTxt += "The background model is estimated from the rock physics models\n";
    }

    std::vector<bool> wavelet = modelSettings_->getEstimateWavelet(0);
    std::vector<float>  angle = modelSettings_->getAngle(0);
    for(int i=0; i<modelSettings_->getNumberOfAngles(0); i++) {
      if(wavelet[i] == true) {
        errTxt += "Wavelet must be given when rock physics models are used.\n";
        errTxt += "  Wavelet is not given for angle "+NRLib::ToString(angle[i]*(180/NRLib::Pi),1)+".\n";
      }
    }
    if((modelSettings_->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0)
      errTxt += "The backround trend can not be written to file when rock physics models are used\n";

  }
}

void
XmlModelFile::checkAngleConsistency(std::string & errTxt) { //Wells can not be moved for time lapse data

  float angle;
  int   i,j,w;
  int   nMoveAngles;
  int   nSeismicAngles = modelSettings_->getNumberOfAngles(0);
  int   nWells         = modelSettings_->getNumberOfWells();
  std::vector<float> seismicAngle = modelSettings_->getAngle(0);

  for(w=0; w<nWells; w++){
    nMoveAngles = modelSettings_->getNumberOfWellAngles(w);
    std::vector<bool> compare(nMoveAngles);

    for( i=0; i<nMoveAngles; i++ ){
      compare[i] = false;
      angle   = modelSettings_->getWellMoveAngle(w,i);

      for( j=0; j<nSeismicAngles; j++ ){
        if( angle == seismicAngle[j])
        {
          if (inputFiles_->getSeismicFile(0,j)=="")
            errTxt += "Seismic data are needed for angle "+NRLib::ToString(angle/float(NRLib::Pi/180))+" to optimize the well locations.\n";
          compare[i] = true;
          break;
        }
      }
      if( compare[i]==false ){
        errTxt += "Unexpected angle "+NRLib::ToString(angle/float(NRLib::Pi/180))+" in <optimize-position> is not equal to any of seismic offset-angles"+".\n";
      }
    }
  }
}

void
XmlModelFile::checkTimeLapseConsistency(std::string & errTxt)
{
  int nTimeLapse = modelSettings_->getNumberOfTimeLapses();

  if(modelSettings_->getEstimateFaciesProb())
    errTxt += "Facies estimation is not allowed for time lapse data.\n";

  for(int i=0; i<nTimeLapse-1; i++){
    if(modelSettings_->getSegyOffset(i) != modelSettings_->getSegyOffset(i+1))
      errTxt += "<segy-start-time> in <survey> need to be the same for all time lapses.\n";

    if(inputFiles_->getWaveletEstIntFileTop(i) != inputFiles_->getWaveletEstIntFileTop(i+1))
      errTxt += "When <top-surface-file> in <wavelet-estimation-interval> is given, it needs to be the same for all time lapses.\n";
    if(inputFiles_->getWaveletEstIntFileBase(i) != inputFiles_->getWaveletEstIntFileBase(i+1))
      errTxt += "When <base-surface-file> in <wavelet-estimation-interval> is given, it needs to be the same for all time lapses.\n";

    if(modelSettings_->getUseLocalNoise(i) != modelSettings_->getUseLocalNoise(i+1))
      errTxt += "If local noise is used for one time lapse, it needs to be used for all time lapses.\n";
  }
  for(int i=0; i<modelSettings_->getNumberOfTimeLapses(); i++){
    if(modelSettings_->getNumberOfAngles(i) == 0)
      errTxt += "Need at least one <angle-gather> in <survey>.\n";
    if(modelSettings_->getVintageYear(i) == RMISSING)
      errTxt += "<year> in <vintage> needs to be specified when using time lapse data.\n";
  }

  if(modelSettings_->getUse3DWavelet()==true)
    errTxt += "3D wavelets can not be used for time lapse data.\n";

  if (modelSettings_->getOptimizeWellLocation())
    errTxt += "The well locations can not be optimized with time lapse data.\n";

  TraceHeaderFormat * thf1;
  TraceHeaderFormat * thf2;
  for(int i=0; i<nTimeLapse-1; i++){
    int nThf1 = modelSettings_->getNumberOfAngles(i);
    int nThf2 = modelSettings_->getNumberOfAngles(i+1);
    int j = 0;
    while(j<nThf1 && j<nThf2){
      thf1 = modelSettings_->getTraceHeaderFormat(i,j);
      thf2 = modelSettings_->getTraceHeaderFormat(i+1,j);
      if(thf1 != NULL && thf2 != NULL){
        if(thf1->GetUtmxLoc() != thf2->GetUtmxLoc())
          errTxt +="When <location-x> in <segy-format> is given, it needs to be the same for all time lapses.\n";
        if(thf1->GetUtmyLoc() != thf2->GetUtmyLoc())
          errTxt += "When <location-y> in <segy-format> is given, it needs to be the same for all time lapses.\n";
        if(thf1->GetInlineLoc() != thf2->GetInlineLoc())
          errTxt += "When <location-il> in <segy-format> is given, it needs to be the same for all time lapses.\n";
        if(thf1->GetCrosslineLoc() != thf2->GetCrosslineLoc())
          errTxt += "When <location-xl> in <segy-format> is given, it needs to be the same for all time lapses.\n";
        if(thf1->GetScalCoLoc() != thf2->GetScalCoLoc())
          errTxt += "When <bypass-coordinate-scaling> in <segy-format> is given, it needs to be the same for all time lapses.\n";
        if(thf1->GetFormatName() != thf2->GetFormatName())
          errTxt += "When <standard-format> in <segy-format> is given, it needs to be the same for all time lapses.\n";
      }
      else if (thf1 != NULL && thf2 == NULL)
        errTxt += "When <segy-format> in <seismic-data> is given, it needs to be identical for all time lapses.\n";
      else if (thf1 == NULL && thf2 != NULL)
        errTxt += "When <segy-format> in <seismic-data> is given, it needs to be identical for all time lapses.\n";
      j++;
    }
  }

  // Check vintages in <rock-physics><evolve>
  // Compare vintages in <reservoir> with vintages in <survey>
  int n_surveys = modelSettings_->getNumberOfVintages();

  std::vector<int> survey_year(n_surveys);
  for(int i=0; i<n_surveys; i++)
    survey_year[i] = modelSettings_->getVintageYear(i);

  std::sort(survey_year.begin(), survey_year.end());

  typedef std::map<std::string, std::vector<DistributionWithTrendStorage *> > my_map;
  const my_map& reservoir_variable = modelSettings_->getReservoirVariable();
  for(my_map::const_iterator it = reservoir_variable.begin(); it != reservoir_variable.end(); it++) {
    if(n_surveys != static_cast<int>(it->second.size()))
      errTxt += "The number of vintages for reservoir variable '"+it->first+"' must be the same as the number of survey vintages\n";
    else {
      std::vector<DistributionWithTrendStorage *> reservoir_evolve = it->second;
      for(int i=1; i<n_surveys; i++) {
        if(reservoir_evolve[i]->GetVintageYear() != survey_year[i]) {
          errTxt += "The vintage years for reservoir variable '"+it->first+"' must be the same as the <vintage><year> given in <survey>\n";
          break;
        }
      }
    }
  }
}

void
XmlModelFile::checkIOConsistency(std::string & errTxt)
{
  if ((modelSettings_->getOtherOutputFlag() & IO::LOCAL_NOISE)>0 && modelSettings_->getUseLocalNoise(0)==false) //When local noise is used, it is used for all time lapses
  {
    LogKit::LogFormatted(LogKit::Warning, "\nWarning: Local noise can not be written to file when <local-noise-scaled> or <estimate-local-noise> is not requested.");
    TaskList::addTask("Remove <local-noise> from <other-output> in the model file if local noise is not used.");
  }
  if ((modelSettings_->getWaveletOutputFlag() & IO::LOCAL_WAVELETS)>0 && modelSettings_->getUseLocalWavelet()==false)
  {
   LogKit::LogFormatted(LogKit::Warning, "\nWarning: Local wavelets can not be written to file when <local-wavelet> is not requested for the angle gathers.");
   TaskList::addTask("Remove <local-wavelets> from <wavelet-output> in the model file if local wavelets are not used.");
  }
  if (((modelSettings_->getWaveletFormatFlag() & IO::NORSARWAVELET)  > 0   ||
       (modelSettings_->getWaveletFormatFlag() & IO::JASONWAVELET)   > 0 ) &&
       (modelSettings_->getWaveletOutputFlag() & IO::LOCAL_WAVELETS)  == 0 &&
       (modelSettings_->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) == 0 &&
       (modelSettings_->getWaveletOutputFlag() & IO::WELL_WAVELETS)   == 0 &&
        modelSettings_->getWaveletFormatManual() == true                   &&
        modelSettings_->getEstimationMode() == false)
  {
    errTxt += "A format is requested in wavelet-output without specifying any of the wavelet\n";
    errTxt += " outputs <well-wavelets>, <global-wavelets> nor <local-wavelets>.\n";
  }
}
