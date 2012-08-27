
#include "src/definitions.h"

#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionwithtrendstorage.h"

void CheckVolumeConsistency(const std::vector<DistributionWithTrendStorage *> & volume_fraction,
                            std::string                                       & errTxt)
{
  int n_constituents = static_cast<int>(volume_fraction.size());

  if(n_constituents > 2) {
    for(int i = 0; i<n_constituents; i++) {
      if(volume_fraction[i]->GetIsDistribution() == true)
        errTxt += "The volume fractions can not be defined by a distribution when more than two constituents are used in a rock physics model\n";
    }
  }

  int n_missing = 0;

  for(int i=0; i<n_constituents; i++) {
    if(volume_fraction[i] == NULL)
      n_missing++;
  }

  if(n_missing == 0)
    errTxt += "One of the volume frations must be unspecified in the rock physics models where elements with corresponding volume frations are given\n";
  else if(n_missing > 1)
    errTxt += "All but one of the volume frations must be defined in the rock physics models where elements with corresponding volume frations are given\n";
}
