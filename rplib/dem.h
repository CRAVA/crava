#ifndef RPLIB_DEM_H
#define RPLIB_DEM_H

#include <vector>

class DEM {
public:
  DEM(const std::vector<double>&       bulk_modulus,
      const std::vector<double>&       shear_modulus,
      const std::vector<double>&       aspect_ratio,
      std::vector<double>&             concentration,
      double                           bulk_modulus_bg,
      double                           shear_modulus_bg);

  virtual ~DEM();

   /*
   function [KEff, MuEff] = geqDEM(K, MU, Aspect, Concentration)
   calculate bulk and shear modulus using DEM

   KEff            - effective bulk modulus
   MuEff           - effective shear modulus
   K               - bulk modulus; [host incl1 incl2 ... inclN]
   MU              - shear modulus; [host incl1 incl2 ... inclN]
   Aspect          - aspect ratios; [incl1 incl2 ... inclN]
   Concentration   - asp.r.conc.; sum([incl1 incl2 ... inclN]) == 1-host

   If multiple inclusions are used, then aspect ratio and concentration_ of
   each must be specified.
   If multiple inclusions with various aspect ratios are used for the same
   inclusion, then moduli of the inclusion must be repeated the same
   number of times as aspect ratios.

   Original code by T. Mukerji, 1997.
   Extended by Erling H. Jensen (2009) to allow for multiple aspect ratios
   and multiple inclusions.
   */
  // first version is possible a slow version which mirrors the original matlab implementation
  void CalcEffectiveModulus(double&                    effective_bulk_modulus,
                            double&                    effective_shear_modulus);

  std::vector<double> GEQDEMYPrime(std::vector<double>&       y,
                                   double                     t);

private:
  double                           bulk_modulus_bg_;
  double                           shear_modulus_bg_;
  const std::vector<double>&       bulk_modulus_;
  const std::vector<double>&       shear_modulus_;
  const std::vector<double>&       aspect_ratio_;
  std::vector<double>&             concentration_;
};


#endif
