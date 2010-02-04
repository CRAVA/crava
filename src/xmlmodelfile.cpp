#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "lib/global_def.h"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/segy.hpp"

#include "src/modelsettings.h"
#include "src/xmlmodelfile.h"
#include "src/definitions.h"
#include "src/inputfiles.h"
#include "src/wavelet.h"
#include "src/fftgrid.h"
#include "src/vario.h"
#include "lib/utils.h"
#include "src/io.h"


XmlModelFile::XmlModelFile(const std::string & fileName)
{
  modelSettings_ = new ModelSettings();
  inputFiles_    = new InputFiles();
  failed_        = false;

  std::ifstream file;
  NRLib::OpenRead(file,fileName);

  if (!file) {
    LogKit::LogFormatted(LogKit::ERROR,"\nERROR: Could not open file %s for reading.\n\n", fileName.c_str());
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
    Utils::writeHeader("Invalid XML file");
    LogKit::LogFormatted(LogKit::ERROR,"\n%s is not a valid XML file. %s In line %d, column %d.", 
                         fileName.c_str(), doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
    if (doc.ErrorId() == 9) { // Not very robust check, but a start 
      LogKit::LogFormatted(LogKit::ERROR,"\nPossible cause: Mis-spelled or forgotten end tag.");
    }
    LogKit::LogFormatted(LogKit::ERROR,"\nAborting\n");
    failed_ = true;
  }
  else {
    std::string errTxt = "";
    if(parseCrava(&doc, errTxt) == false)
      errTxt = "'"+std::string(fileName)+"' is not a crava model file (lacks the <crava> keyword.)\n";

    std::vector<std::string> legalCommands(1);
    checkForJunk(&doc, errTxt, legalCommands);

    checkConsistency(errTxt);

    if(errTxt != "") {
      Utils::writeHeader("Invalid model file");
      LogKit::LogFormatted(LogKit::ERROR,"\n%s is not a valid model file:\n",fileName.c_str());
      LogKit::LogMessage(LogKit::ERROR, errTxt);
      LogKit::LogFormatted(LogKit::ERROR,"\nAborting\n");
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

  parseWellData(root, errTxt);
  parseSurvey(root, errTxt);
  parseProjectSettings(root, errTxt);
  parsePriorModel(root, errTxt);
  parseActions(root, errTxt);

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
    else if(mode == "estimation")
      modelSettings_->setEstimationMode(true);
    else if(mode != "inversion")
      errTxt += "String after <mode> must be either <inversion>, <estimation> or <forward>, found <"+
        mode+"> under command <"+root->ValueStr()+">"+lineColumnText(root)+".\n";
  }

  if (parseInversionSettings(root, errTxt) && mode == "forward")
    errTxt += "Inversion can not be done with the mode 'forward'.\n";

  if (parseEstimationSettings(root, errTxt) && mode == "forward")
    errTxt += "Estimation can not be done with the mode 'forward'.\n";
  
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

  std::string facprob;
  if(parseValue(root, "facies-probabilities", facprob, errTxt) == true) {
    int flag = 0;
    if(modelSettings_->getDefaultGridOutputInd() == false)
      flag = modelSettings_->getGridOutputFlag();
    if(facprob == "absolute") {
      flag = (flag | IO::FACIESPROB);
      modelSettings_->setFaciesProbRelative(false);
    }
    else if(facprob == "relative") {
      flag = (flag | IO::FACIESPROBRELATIVE);
      modelSettings_->setFaciesProbRelative(true);
    }
    else
      errTxt += "Illegal type \'"+facprob+"\' specified for keyword <facies-probabilities>."
        +" Choose from \'absolute\' and \'relative\'\n";

    if(flag != 0) {
      modelSettings_->setGridOutputFlag(flag);
      modelSettings_->setDefaultGridOutputInd(false);
      modelSettings_->setEstimateFaciesProb(true);
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
  legalCommands.push_back("optimize-location-to");

  std::string tmpErr = "";
  std::string value;
  if(parseValue(root, "file-name", value, tmpErr) == true) {
    inputFiles_->addWellFile(value);
  }
  else
    inputFiles_->addWellFile(""); //Dummy to keep tables balanced.

  bool use;
  if(parseBool(root, "use-for-wavelet-estimation", use, tmpErr) == true) {
    if(use == false)
      modelSettings_->addIndicatorWavelet(ModelSettings::NO);
    else
      modelSettings_->addIndicatorWavelet(ModelSettings::YES);
  }
  else
    modelSettings_->addIndicatorWavelet(ModelSettings::NOTSET);


  if(parseBool(root, "use-for-background-trend", use, tmpErr) == true) {
    if(use == false)
      modelSettings_->addIndicatorBGTrend(ModelSettings::NO);
    else
      modelSettings_->addIndicatorBGTrend(ModelSettings::YES);
  }
  else
    modelSettings_->addIndicatorBGTrend(ModelSettings::NOTSET);

  if(parseBool(root, "use-for-facies-probabilities", use, tmpErr) == true) {
    if(use == false)
      modelSettings_->addIndicatorFacies(ModelSettings::NO);
    else
      modelSettings_->addIndicatorFacies(ModelSettings::YES);
  }
  else
    modelSettings_->addIndicatorFacies(ModelSettings::NOTSET);

  bool synth = false;
  if(parseBool(root, "synthetic-vs-log", synth, tmpErr) == true) {
    if(synth == false)
      modelSettings_->addIndicatorRealVs(ModelSettings::YES); //Note the inversion. The keyword is most logical this way,
    else                                                      //but the default should be that the log is real.
      modelSettings_->addIndicatorRealVs(ModelSettings::NO);
  }
  else
    modelSettings_->addIndicatorRealVs(ModelSettings::NOTSET);

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
  TiXmlNode * root = node->FirstChildElement("optimize-location-to");
  if(root == 0)
    return(false);

  modelSettings_->setOptimizeWellLocation(true);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("angle");
  legalCommands.push_back("weight");

  float value;
  if(parseValue(root, "angle", value, errTxt) == true)
    modelSettings_->addMoveAngle(value*float(M_PI/180));
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

  Vario * vario = NULL;
  if(parseVariogram(root, "angular-correlation", vario, errTxt) == true) {
    if (vario != NULL) {
      vario->convertRangesFromDegToRad();
      modelSettings_->setAngularCorr(vario);
    }
  }

  float value;
  if(parseValue(root, "segy-start-time", value, errTxt) == true)
    modelSettings_->setSegyOffset(value);

  while(parseAngleGather(root, errTxt) == true);

  if(modelSettings_->getNumberOfAngles() == 0)
    errTxt +=  "Need at least one angle gather in command <"+root->ValueStr()+">"+
    lineColumnText(root)+".\n";

  parseWaveletEstimationInterval(root, errTxt);

  checkForJunk(root, errTxt, legalCommands);
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
    modelSettings_->addAngle(value*float(M_PI/180));
  else
    errTxt += "Need offset angle for gather"+lineColumnText(root)+".\n";

  if(parseSeismicData(root, errTxt) == false) {
    //Go for defaults. Assume forward model (check later)
    inputFiles_->addSeismicFile("");
    modelSettings_->addSeismicType(ModelSettings::STANDARDSEIS);
  }

  bool oneDwavelet(false);
  if(parseWavelet(root, errTxt) == false) {
    modelSettings_->addWaveletScale(RMISSING);
    modelSettings_->addEstimateGlobalWaveletScale(false); 
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
    inputFiles_->addSeismicFile(""); //Keeping tables balanced.

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
  legalCommands.push_back("scale");
  legalCommands.push_back("estimate-scale");
  legalCommands.push_back("local-wavelet");

  modelSettings_->addStretchFactor(1.0);
  std::string value;
  if(parseFileName(root, "file-name", value, errTxt) == true) {
    inputFiles_->addWaveletFile(value);
    modelSettings_->addEstimateWavelet(false);
  }
  else {
    inputFiles_->addWaveletFile(""); //Keeping tables balanced.
    modelSettings_->addEstimateWavelet(true);
  }
  
  float scale;
  bool scaleGiven = false;
  if(parseValue(root, "scale", scale, errTxt) == true)
  {
    scaleGiven = true;
    modelSettings_->addWaveletScale(scale);
    modelSettings_->addEstimateGlobalWaveletScale(false);
  }
  bool estimate;
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
      errTxt = errTxt +"Error: Can not give both value and ask for estimation of global wavelet scale"+
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

  modelSettings_->setUseLocalWavelet(true);
  std::string filename;
  bool shiftGiven = parseFileName(root, "shift-file", filename, errTxt);
  if(shiftGiven)
    inputFiles_->addShiftFile(filename);
  else
    inputFiles_->addShiftFile("");

  bool scaleGiven = parseFileName(root, "scale-file", filename, errTxt);
  if(scaleGiven)
    inputFiles_->addScaleFile(filename);
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
    else
      modelSettings_->addEstimateLocalShift(estimate);

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
    else
      modelSettings_->addEstimateLocalScale(estimate);

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

  std::string filename;
  if(parseFileName(root, "top-surface-file", filename, errTxt) == true)
    inputFiles_->setWaveletEstIntFile(0, filename);
  if(parseFileName(root, "base-surface-file", filename, errTxt) == true)
    inputFiles_->setWaveletEstIntFile(1, filename);

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
  bool estimate = false;
  if(parseFileName(root, "file-name", value, errTxt) == true) {
    inputFiles_->addWaveletFile(value);
    modelSettings_->addEstimateWavelet(false);
  }
  else {
    inputFiles_->addWaveletFile(""); //Keeping tables balanced.
    modelSettings_->addEstimateWavelet(true);
    estimate = true;
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
  
  if (modelSettings_->getBackgroundType() != "")
    errTxt += "Both background and earth-model can not be given. Under forward mode, the earth-model must be used.\n";
  modelSettings_->setBackgroundType("earth model");

  std::string filename;
  bool vp = parseFileName(root, "vp-file", filename, errTxt);
  if(vp == true) {
    inputFiles_->setBackFile(0, filename);
    modelSettings_->setConstBackValue(0, -1);
  }
  else
    errTxt += "Vp is not given in command earth-model.\n";

  bool vs = parseFileName(root, "vs-file", filename, errTxt);
  if(vs == true) {
    inputFiles_->setBackFile(1, filename);
    modelSettings_->setConstBackValue(1, -1);
  }
  else
    errTxt += "Vs is not given in command earth-model.";

  bool rho = parseFileName(root, "density-file", filename, errTxt);
  if(rho == true) {
    inputFiles_->setBackFile(2, filename);
    modelSettings_->setConstBackValue(2, -1);
  }
  else
    errTxt += "Density is not given in command earth-model.";

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
  legalCommands.push_back("vp-constant");
  legalCommands.push_back("vs-constant");
  legalCommands.push_back("density-constant");
  legalCommands.push_back("velocity-field");
  legalCommands.push_back("lateral-correlation");
  legalCommands.push_back("high-cut-background-modelling");

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

  bool bgGiven = (vp & vs & rho);
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

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseFaciesProbabilities(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("facies-probabilities");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("facies-estimation-interval");
  legalCommands.push_back("prior-facies-probabilities");
  legalCommands.push_back("facies-probability-undefined-value");
  legalCommands.push_back("use-vs-for-facies-probabilities");

  parseFaciesEstimationInterval(root, errTxt);

  parsePriorFaciesProbabilities(root, errTxt);

  float value;
  if(parseValue(root, "facies-probability-undefined-value", value, errTxt) == true)
    modelSettings_->setPundef(value);

  bool useVs = true;
  if(parseBool(root, "use-vs-for-facies-probabilities", useVs, errTxt) == true && useVs == false)
    modelSettings_->setNoVsFaciesProb(true);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}


bool
XmlModelFile::parseFaciesEstimationInterval(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("facies-estimation-interval");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-file-name");
  legalCommands.push_back("base-file-name");

  std::string filename;
  if(parseFileName(root, "top-file-name", filename, errTxt) == true)
    inputFiles_->setFaciesEstIntFile(0, filename);
  else
    errTxt += "Must specify <top-file-name> in command <"+root->ValueStr()+"> "+
      lineColumnText(root)+".\n";
  if(parseFileName(root, "base-file-name", filename, errTxt) == true)
    inputFiles_->setFaciesEstIntFile(1, filename);
  else
    errTxt += "Must specify <base-file-name> in command <"+root->ValueStr()+">"+
      lineColumnText(root)+".\n";

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parsePriorFaciesProbabilities(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("prior-facies-probabilities");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("facies");
  
  float sum;
  int status = 0;
  sum = 0.0;
  int oldStatus = 0;
 
 // std::map<std::string,float> pp;
 // std::map<std::string, std::string> pp2;
  modelSettings_->setPriorFaciesProbGiven(0);
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
     modelSettings_->setPriorFaciesProbGiven(1);
     modelSettings_->addPriorFaciesProb(faciesname,value);
     // status = 1;
    }
    else if(parseValue(root,"probability-cube",filename,errTxt,true)==true)
    {
     modelSettings_->setPriorFaciesProbGiven(2);
     inputFiles_->setPriorFaciesProb(faciesname,filename);
     // status = 2;
    }
 
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
  legalCommands.push_back("time-3D-mapping");
  legalCommands.push_back("io-settings");
  legalCommands.push_back("advanced-settings");

  if(parseOutputVolume(root, errTxt) == false)
    errTxt += "Command <output-volume> is needed in command <"+
      root->ValueStr()+">"+lineColumnText(root)+".\n";

  parseTime3DMapping(root, errTxt);
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
    modelSettings_->setTimeNz(value);

  std::string filename;
  bool externalField = parseFileName(root, "velocity-field", filename, errTxt);
  if(externalField == true)
    inputFiles_->setVelocityField(filename);

  bool inversionField = false;
  parseBool(root, "velocity-field-from-inversion", inversionField, errTxt);
  

  if(inversionField == true) {
    modelSettings_->setVelocityFromInversion(true);
    if(externalField == true)
      errTxt += "Both <velcoity-field> and <velocity-field-from-inversion> given in command <"
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
  std::string filename;
  if(parseFileName(root, "file-name", filename, errTxt) == true)
  {
    inputFiles_->setAreaSurfaceFile(filename);
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

  double dx = 0;
  if(parseValue(root, "sample-density-x", dx, errTxt) == false)
    errTxt += "Sample density for x must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double dy = 0;
  if(parseValue(root, "sample-density-y", dy, errTxt) == false)
    errTxt += "Sample density for y must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double angle = 0;
  if(parseValue(root, "angle", angle, errTxt) == false)
    errTxt += "Rotation angle must be given in command <"+
      root->ValueStr()+"> "+lineColumnText(root)+".\n";

  double rot = (-1)*angle*(M_PI/180.0);
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
  TiXmlNode * root = node->FirstChildElement("time-3D-mapping");
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
  parseOtherOutput(root, errTxt);
  std::string value;
  if(parseValue(root, "file-output-prefix", value, errTxt) == true) {
    IO::setFilePrefix(value);
  }

  std::string level;
  if(parseValue(root, "log-level", level, errTxt) == true) {
    int logLevel = LogKit::ERROR;
    if(level=="error")
      logLevel = LogKit::L_ERROR;
    else if(level=="warning")
      logLevel = LogKit::L_WARNING;
    else if(level=="low")
      logLevel = LogKit::L_LOW;
    else if(level=="medium")
      logLevel = LogKit::L_MEDIUM;
    else if(level=="high")
      logLevel = LogKit::L_HIGH;
    else if(level=="debuglow")
      logLevel = LogKit::L_DEBUGLOW;
    else if(level=="debughigh")
      logLevel = LogKit::L_DEBUGHIGH;
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
    errTxt +=  "Cannot open file '" + logFileName +"' : " + e.what()+"\n";
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
  int paramFlag = 0;
  parseGridElasticParameters(root, paramFlag, errTxt);
  parseGridSeismicData(root, paramFlag, errTxt);
  parseGridOtherParameters(root, paramFlag, errTxt);
  if(paramFlag!=0)
    modelSettings_->setGridOutputFlag(paramFlag);

  //Set output for all FFTGrids.
  FFTGrid::setOutputFlags(modelSettings_->getGridOutputFormat(), 
                          modelSettings_->getGridOutputDomain());

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
  int domainFlag = modelSettings_->getGridOutputDomain();
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
  modelSettings_->setGridOutputDomain(domainFlag);

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
    modelSettings_->setGridOutputFormat(formatFlag);

  checkForJunk(root, errTxt, legalCommands);
  return(true);
}



bool
XmlModelFile::parseGridElasticParameters(TiXmlNode * node, int & paramFlag, std::string & errTxt)
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
  
  bool value = false;
  //int paramFlag = 0;
  if(modelSettings_->getDefaultGridOutputInd() == false) //May have set faciesprobs.
    paramFlag = modelSettings_->getGridOutputFlag();

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
 
  modelSettings_->setDefaultGridOutputInd(false);
 
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseGridSeismicData(TiXmlNode * node, int & paramFlag, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("seismic-data");
  if(root == 0)
    return(false);
  std::vector<std::string> legalCommands;
  legalCommands.push_back("original");
  legalCommands.push_back("synthetic");
  legalCommands.push_back("residuals");
  bool value = false;
  if(parseBool(root, "residuals", value, errTxt) == true && value == true)
    paramFlag += IO::RESIDUAL;
  if(parseBool(root, "original", value, errTxt) == true && value == true)
    paramFlag += IO::ORIGINAL_SEISMIC_DATA;  

  if(parseBool(root, "synthetic", value, errTxt) == true && value == true)
  {
    paramFlag += IO::SYNTHETIC_SEISMIC_DATA;
    modelSettings_->setGenerateSeismicAfterInv(true);
  }
  modelSettings_->setDefaultGridOutputInd(false);
  checkForJunk(root, errTxt, legalCommands);
  return(true);
}

bool
XmlModelFile::parseGridOtherParameters(TiXmlNode * node, int & paramFlag, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("other-parameters");
  if(root == 0)
    return(false);
  std::vector<std::string> legalCommands;
  legalCommands.push_back("correlations");
  legalCommands.push_back("time-to-depth-velocity");
  legalCommands.push_back("extra-grids");

  bool value = false;
  if(parseBool(root, "correlations", value, errTxt) == true && value == true)
    paramFlag += IO::CORRELATION;
  if(parseBool(root, "extra-grids", value, errTxt) == true && value == true)
    paramFlag += IO::EXTRA_GRIDS;
  if(parseBool(root, "time-to-depth-velocity", value, errTxt) == true && value == true)
    paramFlag += IO::TIME_TO_DEPTH_VELOCITY;

  modelSettings_->setDefaultGridOutputInd(false);
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
XmlModelFile::parseOtherOutput(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("other-output");
  if(root == 0)
    return(false);

  std::vector<std::string> legalCommands;
  legalCommands.push_back("wavelets");
  legalCommands.push_back("extra-surfaces");
  legalCommands.push_back("prior-correlations");
  legalCommands.push_back("background-trend-1d");

  bool value;
  int otherFlag = 0;
  if(parseBool(root, "wavelets", value, errTxt) == true && value == true)
    otherFlag += IO::WAVELETS;
  if(parseBool(root, "extra-surfaces", value, errTxt) == true && value == true)
    otherFlag += IO::EXTRA_SURFACES;
  if(parseBool(root, "prior-correlations", value, errTxt) == true && value == true)
    otherFlag += IO::PRIORCORRELATIONS;
  if(parseBool(root, "background-trend-1d", value, errTxt) == true && value == true)
    otherFlag += IO::BACKGROUND_TREND_1D;

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
  legalCommands.push_back("use-intermediate-disk-storage");
  legalCommands.push_back("maximum-relative-thickness-difference");
  legalCommands.push_back("frequency-band");
  legalCommands.push_back("energy-threshold");
  legalCommands.push_back("wavelet-tapering-length");
  legalCommands.push_back("minimum-relative-wavelet-amplitude");
  legalCommands.push_back("maximum-wavelet-shift");
  legalCommands.push_back("white-noise-component");
  legalCommands.push_back("reflection-matrix");
  legalCommands.push_back("kriging-data-limit");
  legalCommands.push_back("debug-level");

  parseFFTGridPadding(root, errTxt);

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
  if(parseValue(root, "white-noise-component", value, errTxt) == true)
    modelSettings_->setWNC(value);
  std::string filename;
  if(parseFileName(root, "reflection-matrix", filename, errTxt) == true)
    inputFiles_->setReflMatrFile(filename);
  int kLimit = 0;
  if(parseValue(root, "kriging-data-limit", kLimit, errTxt) == true) {
    if(modelSettings_->getKrigingParameter() >= 0)
      modelSettings_->setKrigingParameter(kLimit);
  }
  int level = 0;
  if(parseValue(root, "debug-level", level, errTxt) == true)
    modelSettings_->setDebugFlag(level);

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
  else
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


/*
bool
XmlModelFile::parse(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement();
  if(root == 0)
    return(false);

  checkForJunk(root, errTxt);
  return(true);
}
*/


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
    if(stdFormat == "seisworks")
      thf = new TraceHeaderFormat(TraceHeaderFormat::SEISWORKS);
    else if(stdFormat == "iesx")
      thf = new TraceHeaderFormat(TraceHeaderFormat::IESX);
    else if(stdFormat == "SIP")
      thf = new TraceHeaderFormat(TraceHeaderFormat::SIP);
    else {
      errTxt += "Unknown segy-format '"+stdFormat+"' found on line"+
        NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
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
  unsigned int startLength = errTxt.size();
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
XmlModelFile::checkConsistency(std::string & errTxt) {
  if(modelSettings_->getForwardModeling() == true)
    checkForwardConsistency(errTxt);
  else
    checkEstimationInversionConsistency(errTxt);
  if(modelSettings_->getLocalWaveletVario()==NULL) 
      modelSettings_->copyBackgroundVarioToLocalWaveletVario();
  if(modelSettings_->getOptimizeWellLocation()==true)
    checkAngleConsistency(errTxt);
}


void
XmlModelFile::checkForwardConsistency(std::string & errTxt) {
  if(modelSettings_->getForwardModeling() == true) {
    int i;
    for(i=0;i<modelSettings_->getNumberOfAngles();i++)
    {
      modelSettings_->setEstimateSNRatio(i,false);
      if(modelSettings_->getEstimateWavelet(i)==true)
        errTxt+="Wavelet must be given when doing forward modeling. Wavelet is not given for angle "+NRLib::ToString(modelSettings_->getAngle(i)*(180/M_PI),1)+".\n";

    if(inputFiles_->getSeismicFile(i)!="")
      errTxt+="Seismic data should not be given when doing forward modeling.\n";
    }
   if(modelSettings_->getUseLocalNoise()==true)
      errTxt+="Local noise can not be used in forward modeling.\n";

    if(modelSettings_->getUseLocalWavelet()==true)
      errTxt+="Local wavelet can not be used in forward modeling.\n";

    if (modelSettings_->getNumberOfWells() > 0)
      errTxt +="Wells should not be given when doing forward modeling.\n";
    
    if (modelSettings_->getBackgroundType() == "background")
      errTxt += "An earth model needs to be given when doing forward modeling. The background model should not be given.\n";
    else if (modelSettings_->getBackgroundType() == "")
      errTxt += "An earth model needs to be given when doing forward modeling.\n";
  }
}

void
XmlModelFile::checkEstimationInversionConsistency(std::string & /*errTxt*/) {
  //Quite a few things to check here.
}

void
XmlModelFile::checkAngleConsistency(std::string & errTxt) {

  float angle;
  int   i,j,w;
  int   nMoveAngles;
  int   nSeismicAngles = modelSettings_->getNumberOfAngles();
  int   nWells         = modelSettings_->getNumberOfWells();

  for(w=0; w<nWells; w++){
    nMoveAngles = modelSettings_->getNumberOfWellAngles(w); 
    std::vector<bool> compare(nMoveAngles);

    for( i=0; i<nMoveAngles; i++ ){
      compare[i] = false;
      angle   = modelSettings_->getWellMoveAngle(w,i);

      for( j=0; j<nSeismicAngles; j++ ){
        if( angle==modelSettings_->getAngle(j)){
          compare[i] = true;
          break;
        }
      }

      if( compare[i]==false ){
        errTxt += "Unexpected angle "+NRLib::ToString(angle/float(M_PI/180))+" in <optimize-location-to> is not equal to any of seismic offset-angles"+".\n";
      }
    }
  }
}
