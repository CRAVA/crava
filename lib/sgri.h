#ifndef SGRI_H
#define SGRI_H

#include <string>
#include <vector>
#include "nrlib/grid/grid.hpp"

class Sgri : public NRLib::Grid<float> 
{
public:
  Sgri(const std::string & fileName, 
       char              * errText, 
       int               & errCode); //For reading
  ~Sgri();
  
  float         getWaveletValue(float x, float y, float z) const;
  int           sizeOk(float xLim, float yLim, float zLim);
  
private:
  bool          readHeaderFile(const std::string & fileName, char *errText, int &errCode);
  void          readBinaryFile(int n, char *errText, int &errCode);
  void          write1DWL();

  int           dim_;
  int           nGrid_;

  float         x0_;
  float         y0_;
  float         z0_; 

  float         dX_;
  float         dY_;
  float         dZ_;

  float         scaleX_;
  float         scaleY_;
  float         scaleZ_;

  float         rotAngle_;
  float         dipAngle_;

  char          binFileName_[MAX_STRING];

  int           hasComplex_;

  std::string * axisLabels_;
  std::string * gridLabels_;
  std::string   gridValLabel_;

  float         unDef_;
};

#endif
