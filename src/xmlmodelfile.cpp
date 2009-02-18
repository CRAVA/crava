#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "lib/global_def.h"
#include "lib/lib_misc.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/segy.hpp"

#include "src/inputfiles.h"
#include "src/modelsettings.h"
#include "src/xmlmodelfile.h"
#include "src/vario.h"
#include "src/definitions.h"
#include "src/fftgrid.h"

XmlModelFile::XmlModelFile(const char * fileName)
{
  modelSettings_         = new ModelSettings();
  inputFiles_            = new InputFiles();
  failed_                = false;

  TiXmlDocument doc;
	bool loadOkay = doc.LoadFile(fileName);

  if (loadOkay == false) {
    LogKit::LogFormatted(LogKit::ERROR,"Error: Reading of xml-file %s failed. %s Line %d, column %d\n", 
      fileName, doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
    exit(1);
  }

  std::string errTxt = "";
  if(parseCrava(&doc, errTxt) == false)
    errTxt = "Error: '"+std::string(fileName)+"' is not a crava model file (lacks the 'crava' keyword.\n";

  checkForJunk(&doc, errTxt);

  LogKit::LogMessage(LogKit::ERROR, errTxt);
  if(errTxt != "")
    failed_ = true;
  
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

  std::string cmdErr = "";
  parseWellData(root, cmdErr);

  checkForJunk(root, errTxt);
  errTxt += cmdErr;
  return(true);
}

  
bool
XmlModelFile::parseWellData(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("well-data");
  if(root == 0)
    return(false);
  

  parseLogNames(root, errTxt);

  while(parseWell(root, errTxt) == true);

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

  checkForJunk(root, errTxt);
  return(true);
}

bool
XmlModelFile::parseLogNames(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("log-names");
  if(root == 0)
    return(false);

  std::string value;
  if(parseValue(root, "time-log-name", value, errTxt) == true)
    inputFiles_->setLogName(0, value);
  if(parseValue(root, "vp-log-name", value, errTxt) == true)
    inputFiles_->setLogName(1, value);
  if(parseValue(root, "vs-log-name", value, errTxt) == true)
    inputFiles_->setLogName(2, value);
  if(parseValue(root, "density-log-name", value, errTxt) == true)
    inputFiles_->setLogName(3, value);
  if(parseValue(root, "facies-log-name", value, errTxt) == true)
    inputFiles_->setLogName(4, value);
  
  checkForJunk(root, errTxt);
  return(true);
}


bool
XmlModelFile::parseWell(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("well");
  if(root == 0)
    return(false);

  std::string tmpErr = "";
  std::string value;
  if(parseValue(root, "filename", value, tmpErr) == true) {
    inputFiles_->addWellFile(value);
    if(tmpErr == "")
      checkFileOpen(value, root, tmpErr);
  }
  else
    inputFiles_->addWellFile(""); //Dummy to keep tables balanced.

  bool use;
  if(parseBool(root, "use-in-estimate-of-wavelet", use, tmpErr) == true && use == false)
    modelSettings_->addIndicatorWavelet(0);
  else
    modelSettings_->addIndicatorWavelet(1);

  if(parseBool(root, "use-in-estimate-of-background-trend", use, tmpErr) == true && use == false)
    modelSettings_->addIndicatorBGTrend(0);
  else
    modelSettings_->addIndicatorBGTrend(1);

  if(parseBool(root, "use-in-estimate-of-facies-probabilities", use, tmpErr) == true && use == false)
    modelSettings_->addIndicatorFacies(0);
  else
    modelSettings_->addIndicatorFacies(1);

  checkForJunk(root, errTxt, true); //Allow duplicates

  errTxt += tmpErr;
  return(true);
}

bool
XmlModelFile::parseAllowedParameterValues(TiXmlNode * node, std::string & errTxt)
{
  TiXmlNode * root = node->FirstChildElement("allowed-parameter-values");
  if(root == 0)
    return(false);

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

  checkForJunk(root, errTxt);
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
      tmpErr = "Error: Found '"+tmpVal+"' under keyword '"+keyword+"', expected 'yes' or 'no'. This happened in command '"+
        node->ValueStr()+"' on line "+NRLib2::ToString(node->Row())+", column "+NRLib2::ToString(node->Column())+".\n";
    }
  }

  //No junk-clearing call, done in parseValue.
  errTxt += tmpErr;
  return(true);
}

void 
XmlModelFile::checkForJunk(TiXmlNode * root, std::string & errTxt, bool allowDuplicates)
{
  TiXmlNode * child = root->FirstChild();
  while(child != NULL) {
    switch(child->Type()) {
      case TiXmlNode::COMMENT :
      case TiXmlNode::DECLARATION :
        break;
      case TiXmlNode::TEXT :
        errTxt = errTxt + "Error: Unexpected value '"+child->Value()+"' is not part of command '"+root->Value()+
          "' on line "+NRLib2::ToString(child->Row())+", column "+NRLib2::ToString(child->Column())+".\n";
        break;
      case TiXmlNode::ELEMENT :
        errTxt = errTxt + "Error: Unexpected command '"+child->Value()+"' is not part of command '"+root->Value()+
          "' on line "+NRLib2::ToString(child->Row())+", column "+NRLib2::ToString(child->Column())+".\n";
        break;
      default :
        errTxt = errTxt + "Error: Unexpected text '"+child->Value()+"' is not part of command '"+root->Value()+
          "' on line "+NRLib2::ToString(child->Row())+", column "+NRLib2::ToString(child->Column())+".\n";
        break;
    }
    root->RemoveChild(child);
    child = root->FirstChild();
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
          root = parent->FirstChildElement(root->Value());
        }
        errTxt = errTxt +"Error: Found "+NRLib2::ToString(n)+" extra occurences of command '"+root->Value()+"' under command '"+parent->Value()+
          " on line "+NRLib2::ToString(parent->Row())+", column "+NRLib2::ToString(parent->Column())+".\n";
      }
    }
  }
}

bool 
XmlModelFile::checkFileOpen(const std::string & fName, TiXmlNode * node, std::string & errText)
{
  bool error = false;
  if(fName != "*" && fName != "?")
  {
    FILE * file = fopen(fName.c_str(),"r");
    if(file == 0)
    {
      error = true;
      errText = errText +"Error: Failed to open file '"+fName+"' in command '"+node->ValueStr()+
        "' on line "+NRLib2::ToString(node->Row())+", column "+NRLib2::ToString(node->Column())+".\n";
    }
    else
      fclose(file);
  }
  return(error);
}

    