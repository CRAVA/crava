
#include "src/definitions.h"

#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionwithtrendstorage.h"

std::vector<double>
getVolume(const std::vector<DistributionWithTrendStorage *> & volume_fraction,
          const std::string                                 & path,
          const std::vector<std::string>                    & trend_cube_parameters,
          const std::vector<std::vector<double> >           & trend_cube_sampling,
          std::string                                       & errTxt)
{
  int n_constituents = static_cast<int>(volume_fraction.size());

  std::vector<double> volume(n_constituents);

  for(int i=0; i<n_constituents; i++) {
    if(volume_fraction[i]->GetIsDouble() == true) {
      NRLib::TrendStorage * volume_storage = volume_fraction[i]->CloneMean();
      NRLib::Trend        * volume_trend   = volume_storage    ->GenerateTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
      volume[i]                            = volume_trend      ->GetValue(0,0,0);
      delete volume_storage;
      delete volume_trend;
    }
    else {
      errTxt += "The volume fractions in the constituents of the rock physics model\n"
      "can not be trends nor distributions, and need to be given by a value\n";
      volume[i] = 0;
    }
  }

  double sum = 0;
  for(int i=0; i<n_constituents-1; i++) {
    if(volume[i] != RMISSING)
      sum += volume[i];
    else
      errTxt += "The volume fraction needs to be given for all but the last constituent of the rock physics model\n";
  }

  if(volume[n_constituents-1] == RMISSING)
    volume[n_constituents-1] = 1-sum;
  else {
    sum += volume[n_constituents-1];
    if(sum != 1)
      errTxt += "The volume fractions in the constituents of the rock physics model need to sum to one\n";
  }

  return(volume);
}
