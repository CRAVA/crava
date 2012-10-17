#ifndef RPLIB_ROCK_BOUNDING_H
#define RPLIB_ROCK_BOUNDING_H

#include "rplib/rock.h"

#include <vector>

class RockBounding : public Rock {
public:

  RockBounding(const Rock          * upper_rock,
               const Rock          * lower_rock,
               double                porosity,
               double                K_weight,
               double                M_weight,
               std::vector<double>   u);

  RockBounding();

  RockBounding(const RockBounding & old_rock);

  virtual ~RockBounding();

  virtual Rock                        * Clone()                                                 const;

  virtual Rock                        * Evolve(const std::vector<int>         & delta_time,
                                               const std::vector< Rock * >    & rock)           const;

  virtual void                          SetPorosity(double porosity);

  Rock                                * GetUpperRock() const { return upper_rock_; }
  Rock                                * GetLowerRock() const { return lower_rock_; }

private:

  void                                  ComputeSeismicVariables();

  Rock * upper_rock_;                 // Voigt
  Rock * lower_rock_;                 // Reuss
  double porosity_;
  double K_weight_;
  double M_weight_;

};

#endif
