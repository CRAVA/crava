#ifndef RPLIB_DISTRIBUTIONS_STORAGE_KIT_H
#define RPLIB_DISTRIBUTIONS_STORAGE_KIT_H

#include <string>
#include <vector>

#include "rplib/distributionwithtrendstorage.h"

std::vector<double>
getVolume(const std::vector<DistributionWithTrendStorage *> & volume_fraction,
          const std::string                                 & path,
          const std::vector<std::string>                    & trend_cube_parameters,
          const std::vector<std::vector<double> >           & trend_cube_sampling,
          std::string                                       & errTxt);


#endif
