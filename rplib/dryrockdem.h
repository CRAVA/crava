#ifndef RPLIB_DRYROCK_DEM_H
#define RPLIB_DRYROCK_DEM_H


#include "rplib/dryrock.h"

#include "rplib/demmodelling.h"

#include <vector>
#include <cstring>

class DryRockDEM : public DryRock {
public:

  DryRockDEM(const DryRock                       * dryrock,
             const std::vector< DryRock* >       & dryrock_inc,
             const std::vector<double>           & inclusion_spectrum,
             const std::vector<double>           & inclusion_concentration, // the first element is the concentration of the host
             const std::vector<double>           & u);

  DryRockDEM();

  virtual ~DryRockDEM();

  // Assignment operator.
  DryRockDEM                            & operator=(const DryRockDEM& rhs);

  virtual DryRock                       * Clone() const;

  const DryRock                         * GetDryRockHost()              const { return dryrock_        ;}

  const std::vector<DryRock*>           & GetDryRockInclusion()         const { return dryrock_inc_    ;}
  const DryRock                         * GetDryRockInclusion(size_t i) const { return dryrock_inc_[i] ;} // no error checking on valid index range

private:
  //Copy constructor for getting base class variables , used by Clone:
  DryRockDEM(const DryRockDEM & rhs) : DryRock(rhs) {}

  // Calculate elastic and seismic parameters, to be
  // used whenever new information is sent to class.
  void                                    ComputeElasticParams();
  void                                    Clone(const std::vector< DryRock* > & dryrock_in);
  void                                    DeleteInclusion();


  DryRock                               * dryrock_;                 // Owned and deleted by this class.
  std::vector< DryRock* >                 dryrock_inc_;             // Owned and deleted by this class.
  std::vector<double>                     inclusion_spectrum_;
  std::vector<double>                     inclusion_concentration_; // the first element is the concentration of the host
};

#endif
