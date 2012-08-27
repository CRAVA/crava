#ifndef RPLIB_DISTRIBUTIONS_STORAGE_KIT_H
#define RPLIB_DISTRIBUTIONS_STORAGE_KIT_H

#include <string>
#include <vector>

#include "rplib/distributionwithtrendstorage.h"

void CheckVolumeConsistency(const std::vector<DistributionWithTrendStorage *> & volume_fraction,
                            std::string                                       & errTxt);

#endif
