#ifndef SGRI_H
#define SGRI_H

#include <string>
#include <vector>
#include "lib/grid.hpp"

class Sgri : public NRLib2::Grid<float> {
public:
	Sgri(char * fileName, char *errText, int &errCode); //For reading
	~Sgri();

  int                 GetDim() const {return (dim_);}
  double              GetX0() const {return (x0_);}
  double              GetY0() const {return (y0_);}
  double              GetZ0() const {return (z0_);}
  double              GetDX() const {return (dX_);}
  double              GetDY() const {return (dY_);}
  double              GetDZ() const {return (dZ_);}
  
private:
	bool				        readHeaderFile(char * fileName, char *errText, int &errCode);
	void				        readBinaryFile(int n, char *errText, int &errCode);
	int					        dim_;
	int					        nGrid_;
	double				      x0_, y0_, z0_; 
	double				      dX_, dY_, dZ_;
	double				      rotAngle_, dipAngle_;
	double				      scaleX_, scaleY_, scaleZ_;
	char				        binFileName_[MAX_STRING];
	int 				        hasComplex_;
	std::string*	      axisLabels_;
	std::string*	      gridLabels_;
	std::string		      gridValLabel_;
	float				        unDef_;
};

#endif
