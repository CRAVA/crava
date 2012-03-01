#ifndef PDF3DEMPIRICAL_H
#define PDF3DEMPIRICAL_H

// Class for holding a 3D empirical pdf.

class Pdf3DEmpirical : Public Pdf3D {
public:

  Pdf3DEmpirical(const std::vector<double> & d1,
                 const std::vector<double> & d2,
                 const std::vector<double> & d3,
                 int n1,
                 int n2,
                 int n3,
                 smooth_var1,
                 smooth_var2,
                 smooth_var3,
                 smooth_corr12,
                 smooth_corr13,
                 smooth_corr23);
  virtual ~Pdf3DEmpirical()

  // Returns missing if values are outside definition area.
  virtual double density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const;

private:

  FFTGrid * histogram_;

  int n1_;       //Resolution for each varaiable.
  int n2_;
  int n3_;

  double v1_min; //Limits for variable 1
  double v1_max;
  double dv1;    //dv1 = (v1_max-v1_min)/n1. To find index, take floor((v-v_min)/dv)).

  double v2_min;
  double v2_max;
  double dv2;

  double v3_min;
  double v3_max;
  double dv3;
};

#endif
