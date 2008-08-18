// $Id: regularsurface.hpp 69 2008-05-23 13:56:27Z perroe $

#ifndef NRLIB_REGULARSURFACE_HPP
#define NRLIB_REGULARSURFACE_HPP

#include "surface.hpp"
#include "../grid/grid2d.hpp"
#include "../iotools/stringtools.hpp"

#include <cmath>
#include <algorithm>

namespace NRLib2 {

/// Surface represented as a regular grid, where each of the grid cells
/// are modelled as bilinear surfaces.

template <class A>
class RegularSurface : public Grid2D<A>, public Surface {
public:
  RegularSurface();

  RegularSurface(double x0, double y0, double lx, double ly, int nx, int ny, 
                 const A& value = A());
  
  Surface* Clone() const 
  { return new RegularSurface<A>(*this); }

  ///Return z. If outside, throws exception when failOutside is true, returns missing otherwise.
  double GetZ(double x, double y, bool failOutside = true) const;
  
  bool EnclosesRectangle(double x_min, double x_max, 
                         double y_min, double y_max) const;

  void Add(double c) {
    typename Grid2D<A>::iterator i;
    for(i=this->begin();i<this->end();i++)
      *i += c;
  }

  void Multiply(double c) {
    typename Grid2D<A>::iterator i;
    for(i=this->begin();i<this->end();i++)
      *i *= c;
  }

  bool Add(const Surface * s2);
  bool Subtract(const Surface * s2);


  double Min() const {
    return(double(*(std::min_element(this->begin(),this->end()))));
  }

  double Max() const {
    return(double(*(std::max_element(this->begin(),this->end()))));
  }


  /// Returns vector with the four corners around the point (x,y).
  std::vector<A> GetCorners(double x, double y, bool failOutside = true) const;

  /// Gets the index for the nearest point up, and to the left of
  /// (x,y), the corner point with the lowest i and j values.
  void GetIndex(double x, double y, int& i, int& j, bool failOutside = true) const;
  void GetXY(int i, int j, double & x, double & y) const {
    x = x_min_+i*dx_;
    y = y_min_+j*dy_;
  }
  void GetXY(int index, double & x, double & y) const {
    int j = index/this->GetNI();
    int i = index-j*this->GetNI();
    GetXY(i,j,x,y);
  }

  double GetXMin() const     
  { return x_min_; }
  double GetYMin() const
  { return y_min_; }
  double GetXMax() const
  { return x_min_ + lx_; } 
  double GetYMax() const
  { return y_min_ + ly_; } 
  double GetDX() const
  { return dx_; }
  double GetDY() const
  { return dy_; }
  double GetLengthX() const
  { return lx_; }
  double GetLengthY() const
  { return ly_; }

  void SetDimensions(double x_min, double y_min, 
                     double lx, double ly);

  /// Resize grid. Overrides Grid2D's resize.
  void Resize(int ni, int nj, const A& val = A());

  /// Check if grid value is missing
  bool IsMissing(A val) const {
    return val == missing_val_;
  }

private:
  double CellRelX(double x) const;
  double CellRelY(double y) const;

  double x_min_;
  double y_min_;
  double lx_;
  double ly_;
  double dx_;
  double dy_;
  
  /// Missing value
  A missing_val_;
};

template <class A>
RegularSurface<A>::RegularSurface()
  : Grid2D<A>(0, 0),
    x_min_(0),
    y_min_(0),
    lx_(0),
    ly_(0),
    dx_(0),
    dy_(0),
    missing_val_(-999.0)
{}

template <class A>
RegularSurface<A>::RegularSurface(double x_min, double y_min, 
                                  double lx, double ly, 
                                  int nx, int ny, 
                                  const A& value)
  : Grid2D<A>(nx, ny, value),
    x_min_(x_min),
    y_min_(y_min),
    lx_(lx),
    ly_(ly),
    missing_val_(-999.0)
{
  dx_ = lx / (nx-1);
  dy_ = ly / (ny-1);
}


template <class A>
double RegularSurface<A>::GetZ(double x, double y, bool failOutside) const
{
  // Hack to prevent problems when we are at the edge of the surface.
  if(x == GetXMax())
    x = x - 0.00001*dx_;
  if(y == GetYMax())
    y = y - 0.00001*dy_;

  std::vector<A> corner = GetCorners(x, y);

  int n_missing = 0;
  A sum = 0;

  // double sum = Statistics::Sum(corners, Missing);
  // int n_missing = Statistics::NMissing(corners, Missing);
	for (int pix = 0; pix < 4; pix++)
	{
    if (IsMissing(corner[pix]))
			n_missing++;
    else
      sum += corner[pix];
	}

  if (n_missing == 4) {
    // All surrounding values are missing
    if(failOutside == true)
      throw Exception("Grid cell containing (" + ToString(x) + ", " 
        + ToString(y) + ") is missing.");
    else
      return(missing_val_);
  }
  else if (n_missing == 0) {
    double x1 = CellRelX(x);
    double y1 = CellRelY(y);
    return double(GetBilinearZ(corner, x1, y1));
  }
  else {
    return double(sum/(4-n_missing));
  }
}

template <class A>
bool RegularSurface<A>::Add(const Surface * s2)
{
  double x,y,value;
  int i;
  for(i=0;i<this->GetN();i++) {
    GetXY(i, x, y);
    value = s2->GetZ(x,y,false);
    if(s2->IsMissing(value) == false)
      (*this)(i) += value;
    else
      (*this)(i) = missing_val_;
  }
  return(true);
}


template <class A>
bool RegularSurface<A>::Subtract(const Surface * s2)
{
  double x,y,value;
  int i;
  for(i=0;i<this->GetN();i++) {
    GetXY(i, x, y);
    value = s2->GetZ(x,y,false);
    if(s2->IsMissing(value) == false)
      (*this)(i) -= value;
    else
      (*this)(i) = missing_val_;
  }
  return(true);
}


template <class A>
bool RegularSurface<A>::EnclosesRectangle(double x_min, double x_max, 
                                          double y_min, double y_max) const
{
  if (x_min < GetXMin() || x_max > GetXMax() || 
      y_min < GetYMin() || y_max > GetYMax()) {
    return false;
  }

  return true;
}

template <class A>
void RegularSurface<A>::GetIndex(double x, double y, int& i, int& j,
                                 bool failOutside) const
{
  if(failOutside == true && 
     (x < GetXMin() || x > GetXMax() || y < GetYMin() || y > GetYMax())) {
    throw Exception("Point (" + ToString(x) + ", " + ToString(y) + ") is outside the grid.");  
  }

  i = int((x-GetXMin())/GetDX());
  j = int((y-GetYMin())/GetDY());
}


template <class A>
std::vector<A> RegularSurface<A>::GetCorners(double x, double y, bool failOutside) const
{
  int i, j;
  GetIndex(x, y, i , j, failOutside);
  std::vector<A> corners(4);
  
  if(failOutside == true) {
    corners[0] = (*this)(i, j);
    corners[1] = (*this)(i+1, j);
    corners[2] = (*this)(i, j+1);
    corners[3] = (*this)(i+1, j+1);
  }
  else {
    if(i+1 >= 0 && i < this->GetNI() && j+1 >= 0 && j < this->GetNJ()) {
      if(i >=0 && j>= 0)
        corners[0] = (*this)(i, j);
      else
        corners[0] = missing_val_;
      if(i+1 < this->GetNI() && j>= 0)
        corners[1] = (*this)(i+1, j);
      else
        corners[1] = missing_val_;
      if(i >=0 && j+1 < this->GetNJ())
        corners[2] = (*this)(i, j+1);
      else
        corners[2] = missing_val_;
      if(i+1 < this->GetNI() && j+1 < this->GetNJ())
        corners[3] = (*this)(i+1, j+1);
      else
        corners[3] = missing_val_;
    }
    else {
      corners[0] = missing_val_;
      corners[1] = missing_val_;
      corners[2] = missing_val_;
      corners[3] = missing_val_;
    }
  }
  return corners;
}


/// \todo Move to a more suitable location.
template <class A>
A GetBilinearZ(std::vector<A> corners, double x1, double y1) {
  if (corners.size() != 4) {
    throw Exception("Bug: wrong number of corners.");
  }
  
  A a = corners[0];
  A b = corners[1] - a;
  A c = corners[2] - a;
  A d = corners[3] - a - b - c;
  
  return static_cast<A>(a + b*x1 + c*y1 + d*x1*y1);
}


template <class A>
void RegularSurface<A>::SetDimensions(double x_min, double y_min, 
                                      double lx, double ly)
{
  x_min_ = x_min;
  y_min_ = y_min;
  lx_ = lx;
  ly_ = ly;
  dx_ = lx / Grid2D<A>::GetNI() - 1;
  dy_ = ly / Grid2D<A>::GetNJ() - 1;
}


template <class A>
void RegularSurface<A>::Resize(int nx, int ny, const A& value)
{
  Grid2D<A>::Resize(nx, ny, value);
  dx_ = lx_ / nx;
  dy_ = ly_ / ny;
}


template <class A>
double RegularSurface<A>::CellRelX(double x) const
{
  double i;
  return std::modf((x-GetXMin())/GetDX(), &i);
}


template <class A>
double RegularSurface<A>::CellRelY(double y) const
{
  double j;
  return std::modf((y-GetYMin())/GetDY(), &j);
}

} // namespace NRLib2

#endif // NRLIB_REGULARSURFACE_HPP
