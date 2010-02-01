#ifndef VARIO_H
#define VARIO_H

#include <math.h>
#include <string>

class Vario{
public:
  Vario(float range1, float range2 = 1.0, float angle = 0.0);
  virtual ~Vario();

  virtual float  corr(float deltaX, float deltaY) const = 0;
  std::string    getType(void)        const { return type_   ;}
  float          getRange(void)       const { return range1_ ;}
  float          getSubRange(void)    const { return range2_ ;}
  float          getAngle(void)       const { return angle_  ;}
  bool           getAnisotropic(void) const { return (fabs(range1_ - range2_) > 1.0f) ;}
  void           getParams(float& range1, float& range2, float& angle) const { range1 = range1_; range2 = range2_; angle = angle_;}
  void           convertRangesFromDegToRad();
  void           rotateCounterClockwise(float rotAngle);

protected:
  float findDist(float deltaX, float deltaY) const;

protected:
  std::string type_;

private:
  float  range1_;
  float  range2_;
  float  angle_;
};


class SphericalVario : public Vario
{
public:
  SphericalVario(float range1, float range2 = 1.0, float angle =0.0);
  ~SphericalVario();

  float  corr(float deltaX, float deltaY) const; 
};



class GenExpVario : public Vario
{
public:
  GenExpVario(float pot, float range1, float range2=1.0, float angle=0.0);
  ~GenExpVario();

  float  corr(float deltaX, float deltaY) const; 
  float  getPower() const { return pot_; }

private:
  float pot_;

};


#endif
