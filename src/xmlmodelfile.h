#ifndef XMLMODELFILE_H
#define XMLMODELFILE_H

#include <stdio.h>

#include "nrlib/tinyxml/tinyxml.h"

class Vario;
class ModelSettings;
class InputFiles;

class XmlModelFile
{
public:
  XmlModelFile(const char * fileName);
    ~XmlModelFile(void);

  ModelSettings  * getModelSettings(void)  const { return modelSettings_ ;}
  InputFiles     * getInputFiles(void)     const { return inputFiles_    ;}
  bool             getParsingFailed(void)  const { return failed_        ;}

private:  
  bool parseCrava(TiXmlNode * node, std::string & errTxt);

  bool parseActions(TiXmlNode * node, std::string & errTxt);
  bool   parseInversionSettings(TiXmlNode * node, std::string & errTxt);
  bool     parseSimulation(TiXmlNode * node, std::string & errTxt);
  bool   parseEstimationSettings(TiXmlNode * node, std::string & errTxt);

  bool parseWellData(TiXmlNode * node, std::string & errTxt);
  bool   parseLogNames(TiXmlNode * node, std::string & errTxt);
  bool   parseWell(TiXmlNode * node, std::string & errTxt);
  bool   parseOptimizeLocation(TiXmlNode * node, std::string & errTxt);

  bool   parseAllowedParameterValues(TiXmlNode * node, std::string & errTxt);

  bool parseSurvey(TiXmlNode * node, std::string & errTxt);
  bool   parseAngleGather(TiXmlNode * node, std::string & errTxt);
  bool     parseSeismicData(TiXmlNode * node, std::string & errTxt);
  bool     parseWavelet(TiXmlNode * node, std::string & errTxt);
  bool       parseLocalWavelet(TiXmlNode * node, std::string & errTxt);
  bool       parseWaveletEstimationInterval(TiXmlNode * node, std::string & errTxt);
  bool     parseWavelet3D(TiXmlNode * node, std::string & errTxt);

  bool parsePriorModel(TiXmlNode * node, std::string & errTxt);
  bool   parseBackground(TiXmlNode * node, std::string & errTxt);
  bool   parseFaciesProbabilities(TiXmlNode * node, std::string & errTxt);
  bool   parsePriorFaciesProbabilities(TiXmlNode * node, std::string & errTxt);
  bool   parseFaciesEstimationInterval(TiXmlNode * node, std::string & errTxt);  

  bool parseProjectSettings(TiXmlNode * node, std::string & errTxt);
  bool   parseOutputVolume(TiXmlNode * node, std::string & errTxt);
  bool     parseIntervalTwoSurfaces(TiXmlNode * node, std::string & errTxt);
  bool       parseTopSurface(TiXmlNode * node, std::string & errTxt);
  bool       parseBaseSurface(TiXmlNode * node, std::string & errTxt);
  bool     parseIntervalOneSurface(TiXmlNode * node, std::string & errTxt);
  bool     parseArea(TiXmlNode * node, std::string & errTxt);
  bool   parseTime3DMapping(TiXmlNode * node, std::string & errTxt);
  bool   parseIOSettings(TiXmlNode * node, std::string & errTxt);
  bool     parseOutputTypes(TiXmlNode * node, std::string & errTxt);
  bool       parseGridOutput(TiXmlNode * node, std::string & errTxt);
  bool         parseGridDomains(TiXmlNode * node, std::string & errTxt);
  bool         parseGridFormats(TiXmlNode * node, std::string & errTxt);
  bool         parseGridParameters(TiXmlNode * node, std::string & errTxt);
  bool       parseWellOutput(TiXmlNode * node, std::string & errTxt);
  bool         parseWellFormats(TiXmlNode * node, std::string & errTxt);
  bool       parseOtherOutput(TiXmlNode * node, std::string & errTxt);
  bool   parseAdvancedSettings(TiXmlNode * node, std::string & errTxt);
  bool     parseFFTGridPadding(TiXmlNode * node, std::string & errTxt);
  bool     parseFrequencyBand(TiXmlNode * node, std::string & errTxt);
  bool     parseFacies(TiXmlNode * node, std::string & errTxt);
  template <typename T>
  bool parseValue(TiXmlNode * node, const std::string & keyword, T & value, std::string & errTxt, bool allowDuplicates = false);
  bool parseBool(TiXmlNode * node, const std::string & keyword, bool & value, std::string & errTxt, bool allowDuplicates = false);
  bool parseVariogram(TiXmlNode * node, const std::string & keyword, Vario * & vario, std::string & errTxt);
  bool parseTraceHeaderFormat(TiXmlNode * node, const std::string & keyword, TraceHeaderFormat *& thf, std::string & errTxt);
  bool parseFileName(TiXmlNode * node, const std::string & keyword, std::string & filename, std::string & errTxt, bool allowDuplicates = false);
  void checkAngleConsistency(std::string & errTxt);

  void ensureTrailingSlash(std::string & directory);
  void checkForJunk(TiXmlNode * root, std::string & errTxt, const std::vector<std::string> & legalCommands,
                    bool allowDuplicates = false);
  std::string lineColumnText(TiXmlNode * node);

  void checkConsistency(std::string & errTxt);
  void checkForwardConsistency(std::string & errTxt);
  void checkEstimationInversionConsistency(std::string & errTxt);

  void setMissing(int & value)         { value = IMISSING ;}
  void setMissing(float & value)       { value = RMISSING ;}
  void setMissing(double & value)      { value = RMISSING ;}
  void setMissing(std::string & value) { value = ""       ;}

private:

  ModelSettings  * modelSettings_;
  InputFiles     * inputFiles_;
  
  bool             failed_;                // Indicates whether errors ocuured during construction. 
};

template <typename T>
bool
XmlModelFile::parseValue(TiXmlNode * node, const std::string & keyword, T & value, std::string & errTxt, bool allowDuplicates)
{
  setMissing(value);
  TiXmlNode * root = node->FirstChild(keyword);
  if(root == NULL)
    return(false);

  std::vector<std::string> legalCommands(1);

  TiXmlNode * valNode = root->FirstChild();
  while(valNode != NULL && valNode->Type() != TiXmlNode::TEXT)
    valNode = valNode->NextSibling();

  std::string tmpErr = "";
  if(valNode == NULL)
    tmpErr = "Error: No value found under keyword '"+keyword+
      " on line "+NRLib::ToString(root->Row())+", column "+NRLib::ToString(root->Column())+".\n";
  else {
    try {
      value = NRLib::ParseType<T>(valNode->ValueStr());
    }
    catch(NRLib::Exception & e) {
      setMissing(value);
      tmpErr = "Error: "+std::string(e.what())+" on line "+NRLib::ToString(valNode->Row())+", column "+NRLib::ToString(valNode->Column())+".\n";
    }
    root->RemoveChild(valNode);
  }

  checkForJunk(root, tmpErr, legalCommands, allowDuplicates);

  errTxt += tmpErr;
  return(true);
}
#endif
