#ifndef PDF3DGAUSSIAN_H
#define PDF3DGAUSSIAN_H

#include <rplib/pdf3d.h>
#include <rplib/trinormalwith2dtrend.h>

// Class for holding a 3D Gaussian pdf.

class Pdf3DGaussian : public Pdf3D {
public:

  Pdf3DGaussian(TriNormalWith2DTrend * tri_normal_distr);

  virtual ~Pdf3DGaussian();

  virtual double density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const;

private:

  TriNormalWith2DTrend * tri_normal_distr_; //deleted in DistributionsRockTriNormal

};

#endif
