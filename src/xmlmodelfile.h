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

  ModelSettings  * getModelSettings(void)         const { return modelSettings_ ;}
  InputFiles     * getInputFiles(void)            const { return inputFiles_    ;}

  bool             getParsingFailed(void)         const { return failed_        ;}

private:  
  bool parseCrava(TiXmlNode * node, std::string & errTxt);
  bool parseWellData(TiXmlNode * node, std::string & errTxt);
  bool parseLogNames(TiXmlNode * node, std::string & errTxt);
  bool parseWell(TiXmlNode * node, std::string & errTxt);
  bool parseAllowedParameterValues(TiXmlNode * node, std::string & errTxt);


  template <typename T>
  bool parseValue(TiXmlNode * node, const std::string & keyword, T & value, std::string & errTxt, bool allowDuplicates = false);
  bool parseBool(TiXmlNode * node, const std::string & keyword, bool & value, std::string & errTxt, bool allowDuplicates = false);
  void checkForJunk(TiXmlNode * root, std::string & errTxt, bool allowDuplicates = false);
  bool checkFileOpen(const std::string & fName, TiXmlNode * node, std::string & errText);
 
  void setMissing(int & value) {value = IMISSING;}
  void setMissing(float & value) {value = RMISSING;}
  void setMissing(double & value) {value = RMISSING;}
  void setMissing(std::string & value) {value = "";}



//  Vario          * createVario(TiXmlNode * root, const std::string & command, std::string & errText);

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

  TiXmlNode * valNode = root->FirstChild();
  while(valNode != NULL && valNode->Type() != TiXmlNode::TEXT)
    valNode = valNode->NextSibling();

  std::string tmpErr = "";
  if(valNode == NULL)
    tmpErr = "Error: No value found under keyword '"+keyword+
      " on line "+NRLib2::ToString(root->Row())+", column "+NRLib2::ToString(root->Column())+".\n";
  else {
    try {
      value = NRLib2::ParseType<T>(valNode->ValueStr());
    }
    catch(NRLib2::Exception & e) {
      setMissing(value);
      tmpErr = "Error: "+std::string(e.what())+" on line "+NRLib2::ToString(valNode->Row())+", column "+NRLib2::ToString(valNode->Column())+".\n";
    }
    root->RemoveChild(valNode);
  }

  checkForJunk(root, tmpErr, allowDuplicates);

  errTxt += tmpErr;
  return(true);
}


#endif
