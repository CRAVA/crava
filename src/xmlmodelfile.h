/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef XMLMODELFILE_H
#define XMLMODELFILE_H

#include <stdio.h>

#include "nrlib/segy/traceheader.hpp"
#include "nrlib/tinyxml/tinyxml.h"
#include "src/definitions.h"
#include "nrlib/trend/trendstorage.hpp"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionssolidstorage.h"

class Vario;
class ModelSettings;
class InputFiles;

class XmlModelFile
{
public:
  XmlModelFile(const std::string & fileName);
    ~XmlModelFile(void);

  ModelSettings * getModelSettings(void)  const { return modelSettings_ ;}
  InputFiles    * getInputFiles(void)     const { return inputFiles_    ;}
  bool            getParsingFailed(void)  const { return failed_        ;}

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
  bool   parseWellMoveDataInterval(TiXmlNode * node, std::string & errTxt);

  bool parseSurvey(TiXmlNode * node, std::string & errTxt);
  bool parseAngleGather(TiXmlNode * node, std::string & errTxt);
  bool parseSeismicData(TiXmlNode * node, std::string & errTxt);
  bool parseWavelet(TiXmlNode * node, std::string & errTxt);
  bool parseLocalWavelet(TiXmlNode * node, std::string & errTxt);
  bool parseWaveletEstimationInterval(TiXmlNode * node, std::string & errTxt);
  bool parseWavelet3D(TiXmlNode * node, std::string & errTxt);
  bool parseTravelTime(TiXmlNode * node, std::string & errTxt);
  bool parseRMSVelocities(TiXmlNode * node, std::string & errTxt);
  bool parseTimeGradientSettings(TiXmlNode * node, std::string & errTxt);
  bool parseGravimetry(TiXmlNode * node, std::string & errTxt);
  bool parseVintage(TiXmlNode * node, std::string & errTxt);

  bool parsePriorModel(TiXmlNode * node, std::string & errTxt);
    bool parseCorrelationDirection(TiXmlNode * node, std::string & errTxt);
      bool parseIntervalCorrelationDirection(TiXmlNode * node, std::string & errTxt);
  bool parseEarthModel(TiXmlNode * node, std::string & errTxt);
  bool parsePriorLocalWavelet(TiXmlNode * node, std::string & errTxt);
  bool   parseBackground(TiXmlNode * node, std::string & errTxt);
  bool parseMultizoneModel(TiXmlNode * node, std::string & errTxt);
  bool parseZone(TiXmlNode * node, std::string & errTxt);
  bool   parseFaciesProbabilities(TiXmlNode * node, std::string & errTxt);
  bool   parsePriorFaciesProbabilities(TiXmlNode * node, std::string & errTxt);
  bool    parseFaciesInterval(TiXmlNode * node, std::string & errTxt);
  bool      parseFaciesPerInterval(TiXmlNode * node, std::map<std::string, float> & facies_map, float prob, std::string & errTxt);
  bool   parseVolumeFractions(TiXmlNode * node, std::string & errTxt);
  bool   parseVolumeFractionsInterval(TiXmlNode * node, std::string & errTxt);
  bool      parseVolumeFractionsPerInterval(TiXmlNode * node, std::map<std::string, float> & fraction_map, float prob, std::string & errTxt);
  bool    parseFaciesVolumeFractions(TiXmlNode * node, std::string & errTxt);
  bool   parseFaciesEstimationInterval(TiXmlNode * node, std::string & errTxt);
  bool parseRockPhysics(TiXmlNode * node, std::string & errTxt);
  bool parseRock(TiXmlNode * node, std::string & label, std::string & errTxt);
  bool parsePredefinitions(TiXmlNode * node, std::string & errTxt);
  bool parseSolid(TiXmlNode * node, std::string & label, std::string & errTxt);
  bool parseDryRock(TiXmlNode * node, std::string & label, std::string & errTxt);
  bool parseFluid(TiXmlNode * node, std::string & label, std::string & errTxt);
  bool parseTabulated(TiXmlNode * node, int constituent, std::string label, std::vector<DistributionWithTrendStorage *> total_porosity, std::vector<DistributionWithTrendStorage *>   mineral_k, std::string & errTxt);
  bool parseTabulatedDryRock(TiXmlNode * node, int constituent, std::string label, std::vector<DistributionWithTrendStorage *> total_porosity, std::vector<DistributionWithTrendStorage *>   mineral_k, std::string & errTxt);
  bool parseTabulatedFluid(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseReuss(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseVoigt(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseHill(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseConstituent(TiXmlNode * node, std::string & constituent_label, std::vector<DistributionWithTrendStorage *> & volume_fraction, std::string & errTxt);
  bool parseBatzleWangBrine(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseWalton(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseDEM(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseDEMHost(TiXmlNode * node, std::string & label, std::vector<DistributionWithTrendStorage *> & volume_fraction, std::string & errTxt, bool & missing_vol_frac);
  bool parseDEMInclusion(TiXmlNode * node, std::string & label, std::vector<DistributionWithTrendStorage *> & aspect_ratio, std::vector<DistributionWithTrendStorage *> & volume_fraction, std::string & errTxt, bool & missing_vol_frac);
  bool parseGassmann(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseBounding(TiXmlNode * node, int constituent, std::string label, std::string & errTxt);
  bool parseUpperBound(TiXmlNode * node, std::string & label, std::string & errTxt);
  bool parseLowerBound(TiXmlNode * node, std::string & label, std::string & errTxt);
  bool parseInclusion(TiXmlNode * node, std::string & errTxt);
  bool parseReservoir(TiXmlNode * node, std::string & errTxt);
  bool parseEvolve(TiXmlNode * node, std::string & errTxt);
  bool parseEvolveVintage(TiXmlNode * node, std::vector<DistributionWithTrendStorage *> & reservoir_variable, std::string & errTxt);
  bool parseGaussianWithTrend(TiXmlNode * node, std::vector<DistributionWithTrendStorage *> & storage, bool is_shared, std::string & errTxt);
  bool parseBetaWithTrend(TiXmlNode * node, std::vector<DistributionWithTrendStorage *> & storage, bool is_shared, std::string & errTxt);
  bool parseBetaEndMassWithTrend(TiXmlNode * node, std::vector<DistributionWithTrendStorage *> & storage, bool is_shared, std::string & errTxt);
  bool parse1DTrend(TiXmlNode * node, const std::string & keyword, NRLib::TrendStorage *& trend, std::string & errTxt);
  bool parse2DTrend(TiXmlNode * node, const std::string & keyword, NRLib::TrendStorage *& trend, std::string & errTxt);
  bool parseTrendCube(TiXmlNode * node, std::string & errTxt);

  bool parseProjectSettings(TiXmlNode * node, std::string & errTxt);
  bool   parseOutputVolume(TiXmlNode * node, std::string & errTxt);
  bool     parseIntervalTwoSurfaces(TiXmlNode * node, std::string & errTxt);
  bool       parseTopSurface(TiXmlNode * node, std::string & errTxt);
  bool       parseBaseSurface(TiXmlNode * node, std::string & errTxt);
  bool     parseIntervalOneSurface(TiXmlNode * node, std::string & errTxt);
  bool     parseMultipleIntervals(TiXmlNode * node, std::string & errTxt);
  bool      parseInterval(TiXmlNode * node, std::string & errTxt);
  bool        parseIntervalBaseSurface(TiXmlNode * node, std::string & interval_name, std::string & errTxt);
  bool     parseAreaFromSurface(TiXmlNode * node, std::string & errTxt);
  bool     parseUTMArea(TiXmlNode * node, std::string & errTxt);
  bool     parseILXLArea(TiXmlNode * node, std::string & errTxt);
  bool   parseTime3DMapping(TiXmlNode * node, std::string & errTxt);
  bool   parseIOSettings(TiXmlNode * node, std::string & errTxt);
  bool       parseGridOutput(TiXmlNode * node, std::string & errTxt);
  bool         parseGridDomains(TiXmlNode * node, std::string & errTxt);
  bool         parseGridFormats(TiXmlNode * node, std::string & errTxt);
  bool         parseGridElasticParameters(TiXmlNode * node, std::string & errTxt);
  bool         parseGridSeismicData(TiXmlNode * node, std::string & errTxt);
  bool         parseGridOtherParameters(TiXmlNode * node, std::string & errTxt);
  bool         parseWellOutput(TiXmlNode * node, std::string & errTxt);
  bool         parseWellFormats(TiXmlNode * node, std::string & errTxt);
  bool         parseWaveletOutput(TiXmlNode * node, std::string & errTxt);
  bool         parseWaveletFormats(TiXmlNode * node, std::string & errTxt);
  bool       parseOtherOutput(TiXmlNode * node, std::string & errTxt);
  bool   parseAdvancedSettings(TiXmlNode * node, std::string & errTxt);
  bool     parseFFTGridPadding(TiXmlNode * node, std::string & errTxt);
  bool     parseVpVsRatio(TiXmlNode * node, std::string & errTxt);
  bool       parseIntervalVpVs(TiXmlNode * node, std::string interval_name, std::string & errTxt);
  bool     parseFrequencyBand(TiXmlNode * node, std::string & errTxt);
  bool     parseFacies(TiXmlNode * node, std::string & errTxt);
  template <typename T>
  bool parseValue(TiXmlNode * node, const std::string & keyword, T & value, std::string & errTxt, bool allowDuplicates = false);
  bool parseDistributionWithTrend(TiXmlNode                                   * node,
                                  const std::string                           & keyword,
                                  std::vector<DistributionWithTrendStorage *> & storage,
                                  std::string                                 & label,
                                  bool                                          is_shared,
                                  std::string                                 & errTxt,
                                  bool                                          allowDistribution = true);
  bool parseBool(TiXmlNode * node, const std::string & keyword, bool & value, std::string & errTxt, bool allowDuplicates = false);
  bool parseVariogram(TiXmlNode * node, const std::string & keyword, Vario * & vario, std::string & errTxt);
  bool parseTraceHeaderFormat(TiXmlNode * node, const std::string & keyword, TraceHeaderFormat *& thf, std::string & errTxt);
  bool parseFileName(TiXmlNode * node, const std::string & keyword, std::string & filename, std::string & errTxt, bool allowDuplicates = false);

  void FindDoubleValueFromDistributionWithTrend(const std::vector<DistributionWithTrendStorage *> & dist_with_trend,
                                                std::string                                         type,
                                                std::vector<double>                               & value,
                                                std::string                                       & errTxt) const;

  void checkAngleConsistency(std::string & errTxt);

  void ensureTrailingSlash(std::string & directory);
  void checkForJunk(TiXmlNode * root, std::string & errTxt, const std::vector<std::string> & legalCommands,
                    bool allowDuplicates = false);
  std::string lineColumnText(TiXmlNode * node);

  void setDerivedParameters(std::string & errTxt);
  void checkConsistency(std::string & errTxt);
  void checkForwardConsistency(std::string & errTxt);
  void checkEstimationConsistency(std::string & errTxt);
  void checkInversionConsistency(std::string & errTxt);
  void checkTimeLapseConsistency(std::string & errTxt);
  void checkRockPhysicsConsistency(std::string & errTxt);
  void checkIOConsistency(std::string & errTxt);
  void checkMultizoneBackgroundConsistency(std::string & errTxt);

  void setMissing(int & value)         { value = IMISSING ;}
  void setMissing(float & value)       { value = RMISSING ;}
  void setMissing(double & value)      { value = RMISSING ;}
  void setMissing(std::string & value) { value = ""       ;}

private:

  ModelSettings  * modelSettings_;
  InputFiles     * inputFiles_;

  bool             failed_;                // Indicates whether errors ocuured during construction.
  bool             surveyFailed_;
};

template <typename T>
bool
XmlModelFile::parseValue(TiXmlNode         * node,
                         const std::string & keyword,
                         T                 & value,
                         std::string       & errTxt,
                         bool                allowDuplicates)
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
