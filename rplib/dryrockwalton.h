#ifndef RPLIB_DRYROCK_WALTON_H
#define RPLIB_DRYROCK_WALTON_H


#include "rplib/dryrock.h"

class Solid;

class DryRockWalton : public DryRock {
public:

  DryRockWalton(const Solid                       * solid,
                const double                      & friction_weight,
                const double                      & pressure,
                const double                      & porosity,
                const double                      & coord_number, //if coord_number is negative->interpolatevalue from dataset from Murphy
                const std::vector<double>         & u);

  DryRockWalton();

  virtual ~DryRockWalton();

  // Assignment operator.
  DryRockWalton                         & operator=(const DryRockWalton& rhs);

  virtual DryRock                       * Clone() const;

  const Solid                           * GetSolid() const { return solid_;}

  virtual void                            SetTotalPorosity(double porosity);

private:
  //Copy constructor for getting base class variables , used by Clone:
  DryRockWalton(const DryRockWalton & rhs) : DryRock(rhs) {}

  // Calculate elastic and seismic parameters, to be
  // used whenever new information is sent to class.
  void                                    ComputeElasticParams();
  double                                  CalcCoordNumber() const;


  Solid                                 * solid_;
  double                                  friction_weight_; //0 = no friction, 1 = high friction
  double                                  pressure_;
  double                                  coord_number_;
  bool                                    calc_coord_number_;


};

#endif
