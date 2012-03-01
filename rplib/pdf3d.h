#ifndef PDF3D_H
#define PDF3D_H

// Abstract class for holding a 3D pdf, either empirical or based on a distribution,
// as decided by the derived class. Provides the density for the given set of values.
// Returns missing if values are outside definition area.

class Pdf3D {
public:

  Pdf3D() {}

  virtual ~Pdf3D() {}

  virtual double density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const = 0;

};

#endif
