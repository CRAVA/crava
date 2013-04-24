#ifndef RPLIB_DRYROCK_MIX_H
#define RPLIB_DRYROCK_MIX_H

#include <vector>
#include <cstring>

#include "rplib/dryrock.h"
#include "rplib/demmodelling.h"

class Solid;

//This file contains two classes DryRockMix and DryRockMixOfDryRockAndSolid.

//-------------------------------------- DryRockMix ---------------------------------------------------------

class DryRockMix : public DryRock {
public:

  DryRockMix(const std::vector<DryRock*>    & dryrock,
             const std::vector<double>      & volume_fraction,
             const std::vector<double>      & u,
             DEMTools::MixMethod              mix_method);


  virtual ~DryRockMix();

  // Assignment operator.
  DryRockMix                        & operator=(const DryRockMix& rhs);

  virtual DryRock                   * Clone() const;

  DryRock                           * GetSubDryRock(size_t i) const { return dryrock_[i]; }

private:
  //Copy constructor for getting base class variables , used by Clone:
  DryRockMix(const DryRockMix & rhs) : DryRock(rhs) {}

  std::vector<DryRock*>               dryrock_;           // Owned and deleted by this class.
  std::vector<double>                 volume_fraction_;
  DEMTools::MixMethod                 mix_method_;
};

//-------------------------------------- DryRockMixOfDryRockAndSolid ---------------------------------------------------------

class DryRockMixOfDryRockAndSolid : public DryRock {
public:

  DryRockMixOfDryRockAndSolid(const std::vector<DryRock*> &  dryrock,
                              const std::vector<Solid*>   &  solid,
                              const std::vector<double>   &  volume_fraction_dryrock,
                              const std::vector<double>   &  volume_fraction_solid,
                              const std::vector<double>   &  u,
                              DEMTools::MixMethod            mix_method);


  virtual ~DryRockMixOfDryRockAndSolid();

  // Assignment operator.
  DryRockMixOfDryRockAndSolid       & operator=(const DryRockMixOfDryRockAndSolid& rhs);

  virtual DryRock                   * Clone() const;

  DryRock                           * GetSubDryRock(size_t i) const { return dryrock_[i]; }
  Solid                             * GetSubSolid(size_t i) const { return solid_[i]; }

private:
  //Copy constructor for getting base class variables , used by Clone:
  DryRockMixOfDryRockAndSolid(const DryRockMixOfDryRockAndSolid & rhs) : DryRock(rhs) {}

  void                                ComputeSeismicVariables();

  std::vector< DryRock* >             dryrock_;               // Owned and deleted by this class.
  std::vector< Solid* >               solid_;                 // Owned and deleted by this class.
  std::vector<double>                 volume_fraction_dryrock_;
  std::vector<double>                 volume_fraction_solid_;
  DEMTools::MixMethod                 mix_method_;
};

#endif
