#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/iotools/stringtools.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsdryrocktabulated.h"
#include "rplib/distributionsdryrockdem.h"
#include "rplib/distributionsdryrockmix.h"
#include "rplib/distributionsdryrockwalton.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"

#include "src/modelsettings.h"



DistributionsDryRockStorage::DistributionsDryRockStorage()
{
}

DistributionsDryRockStorage::~DistributionsDryRockStorage()
{
}

std::vector<DistributionsDryRock *>
DistributionsDryRockStorage::CreateDistributionsDryRockMix(const int                                                       & n_vintages,
                                                           const std::string                                               & path,
                                                           const std::vector<std::string>                                  & trend_cube_parameters,
                                                           const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                           const std::map<std::string, DistributionsDryRockStorage *>      & model_dryrock_storage,
                                                           const std::map<std::string, DistributionsSolidStorage *>        & model_solid_storage,
                                                           const std::vector<std::string>                                  & constituent_label,
                                                           const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                           DEMTools::MixMethod                                               mix_method,
                                                           std::string                                                     & errTxt) const
{
  int n_constituents = static_cast<int>(constituent_label.size());

  bool mix_dryrock  = false;
  bool mix_solid    = false;

  std::vector<int> constituent_type(n_constituents);

  FindMixTypesForDryRock(constituent_label,
                         n_constituents,
                         model_solid_storage,
                         model_dryrock_storage,
                         mix_dryrock,
                         mix_solid,
                         constituent_type,
                         errTxt);

  if (mix_dryrock == true && mix_solid == false) { //mix of only dryrocks
    std::vector<int> n_vintages_constit(n_constituents);
    for(int i=0; i<n_constituents; i++)
      n_vintages_constit[i] = static_cast<int>(constituent_volume_fraction[i].size());

    std::vector<double> alpha(n_constituents);
    for(int i=0; i<n_constituents; i++) {
      if(constituent_volume_fraction[i][0] != NULL)
        alpha[i] = constituent_volume_fraction[i][0]->GetOneYearCorrelation();
      else
        alpha[i] = 1;
    }

    std::vector<std::vector<DistributionsDryRock *> > distr_dryrock(n_vintages);
    for(int i=0; i<n_vintages; i++)
      distr_dryrock[i].resize(n_constituents, NULL);

    for (int s = 0; s < n_constituents; s++) {
      std::vector<DistributionsDryRock *> distr_dryrock_all_vintages = ReadDryRock(n_vintages,
                                                                                   constituent_label[s],
                                                                                   path,
                                                                                   trend_cube_parameters,
                                                                                   trend_cube_sampling,
                                                                                   model_dryrock_storage,
                                                                                   model_solid_storage,
                                                                                   errTxt);

      for(int i=0; i<n_vintages; i++) {
        if(i < static_cast<int>(distr_dryrock_all_vintages.size()))
          distr_dryrock[i][s] = distr_dryrock_all_vintages[i];
        else
          distr_dryrock[i][s] = distr_dryrock[i-1][s]->Clone();
      }
    }

    const std::vector<std::vector<double> > dummy_blocked_logs;
    const std::vector<std::vector<double> > dummy_s1;
    const std::vector<std::vector<double> > dummy_s2;
    const int                               dummy_output_other = -999;

    for(int i=0; i<n_constituents; i++)
      CheckValuesInZeroOne(constituent_volume_fraction[i], "volume-fraction", path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, errTxt);

    std::vector<DistributionsDryRock *>                final_dist_dryrock(n_vintages, NULL);
    std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);

    for(int i=0; i<n_vintages; i++)
      all_volume_fractions[i].resize(n_constituents, NULL);

    for(int i=0; i<n_vintages; i++) {

      for (int s=0; s<n_constituents; s++) {

        if(i < n_vintages_constit[s]) {
          if(constituent_volume_fraction[s][i] != NULL)
            all_volume_fractions[i][s] = constituent_volume_fraction[s][i]->GenerateDistributionWithTrend(path,
                                                                                                          trend_cube_parameters,
                                                                                                          trend_cube_sampling,
                                                                                                          dummy_blocked_logs,
                                                                                                          dummy_s1,
                                                                                                          dummy_s2,
                                                                                                          dummy_output_other,
                                                                                                          "dummy",
                                                                                                          errTxt);
        }
        else {
          if(all_volume_fractions[i-1][s] != NULL)
            all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
        }
      }

      CheckVolumeConsistency(all_volume_fractions[i], errTxt);
    }

    if (errTxt == "") {
      for(int i=0; i<n_vintages; i++)
        final_dist_dryrock[i] = new DistributionsDryRockMix(distr_dryrock[i], all_volume_fractions[i], mix_method, alpha);

      for(int i=0; i<n_vintages; i++) {
        for(size_t s=0; s<distr_dryrock.size(); s++)
          delete distr_dryrock[i][s];

        for(size_t s=0; s<all_volume_fractions.size(); s++) {
          if(all_volume_fractions[i][s] != NULL) {
            if(all_volume_fractions[i][s]->GetIsShared() == false)
              delete all_volume_fractions[i][s];
          }
        }
      }
    }

    return(final_dist_dryrock);
  }
  else if (mix_dryrock == true && mix_solid == true) { //mix of dryrocks and solids

    std::vector<std::vector<DistributionsDryRock *> >    distr_dryrock(n_vintages);
    std::vector<std::vector<DistributionsSolid *> >      distr_solid(n_vintages);
    std::vector<std::vector<DistributionWithTrend *> >   dryrock_volume_fractions(n_vintages);
    std::vector<std::vector<DistributionWithTrend *> >   solid_volume_fractions(n_vintages);
    std::vector<double>                                  solid_alpha;
    std::vector<double>                                  dryrock_alpha;

    std::vector<int> n_vintages_constit(n_constituents);
    for(int i=0; i<n_constituents; i++)
      n_vintages_constit[i] = static_cast<int>(constituent_volume_fraction[i].size());

    std::vector<double> alpha(n_constituents);
    for(int i=0; i<n_constituents; i++) {
      if(constituent_volume_fraction[i][0] != NULL)
        alpha[i] = constituent_volume_fraction[i][0]->GetOneYearCorrelation();
      else
        alpha[i] = 1;
    }

    const std::vector<std::vector<double> > dummy_blocked_logs;
    const std::vector<std::vector<double> > dummy_s1;
    const std::vector<std::vector<double> > dummy_s2;
    const int                               dummy_output_other = -999;

    std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);

    for(int i=0; i<n_vintages; i++)
      all_volume_fractions[i].resize(n_constituents, NULL);

    for(int i=0; i<n_vintages; i++) {

      for (int s=0; s<n_constituents; s++) {

        if(i < n_vintages_constit[s]) {
          if(constituent_volume_fraction[s][i] != NULL)
            all_volume_fractions[i][s] = constituent_volume_fraction[s][i]->GenerateDistributionWithTrend(path,
                                                                                                          trend_cube_parameters,
                                                                                                          trend_cube_sampling,
                                                                                                          dummy_blocked_logs,
                                                                                                          dummy_s1,
                                                                                                          dummy_s2,
                                                                                                          dummy_output_other,
                                                                                                          "dummy",
                                                                                                          errTxt);
        }
        else {
          if(all_volume_fractions[i-1][s] != NULL)
            all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
        }
      }

      CheckVolumeConsistency(all_volume_fractions[i], errTxt);
    }

    for(int i=0; i<n_vintages; i++)
      distr_dryrock[i].resize(n_constituents, NULL);

    int n_dryrocks = 0;
    int n_solids = 0;

    for (int s = 0; s < n_constituents; s++) {
      if (constituent_type[s] == ModelSettings::DRY_ROCK) {
        std::vector<DistributionsDryRock *> distr_dryrock_all_vintages = ReadDryRock(n_vintages,
                                                                                     constituent_label[s],
                                                                                     path,
                                                                                     trend_cube_parameters,
                                                                                     trend_cube_sampling,
                                                                                     model_dryrock_storage,
                                                                                     model_solid_storage,
                                                                                     errTxt);
        int n_vintages_dryrock = static_cast<int>(distr_dryrock_all_vintages.size());

        for(int i=0; i<n_vintages; i++)
          distr_dryrock[i].resize(n_dryrocks+1);

        for(int i=0; i<n_vintages; i++) {
          if(i < n_vintages_dryrock)
            distr_dryrock[i][n_dryrocks] = distr_dryrock_all_vintages[i];
          else
            distr_dryrock[i][n_dryrocks] = distr_dryrock[i-1][n_dryrocks]->Clone();
        }

        for(int i=0; i<n_vintages; i++)
          dryrock_volume_fractions[i].push_back(all_volume_fractions[i][s]);

        dryrock_alpha.push_back(alpha[s]);

        n_dryrocks++;

      }
      else if (constituent_type[s] == ModelSettings::SOLID) {
        std::vector<DistributionsSolid *> constit_solid_all_vintages = ReadSolid(n_vintages,
                                                                                 constituent_label[s],
                                                                                 path,
                                                                                 trend_cube_parameters,
                                                                                 trend_cube_sampling,
                                                                                 model_solid_storage,
                                                                                 errTxt);

        int n_vintages_solid = static_cast<int>(constit_solid_all_vintages.size());

        for(int i=0; i<n_vintages_solid; i++)
          distr_solid[i].resize(n_solids+1);

        for(int i=0; i<n_vintages; i++) {
          if(i < n_vintages_solid)
            distr_solid[i][n_solids] = constit_solid_all_vintages[i];
          else
            distr_solid[i][n_solids] = distr_solid[i-1][n_solids]->Clone();
        }

        for(int i=0; i<n_vintages; i++)
          solid_volume_fractions[i].push_back(all_volume_fractions[i][s]);

        solid_alpha.push_back(alpha[s]);
        n_solids++;
      }
    }
    // Now we have read the subsolids and subdryrocks.

    for(int i=0; i<n_constituents; i++)
      CheckValuesInZeroOne(constituent_volume_fraction[i], "volume-fraction", path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, errTxt);

    std::vector<DistributionsDryRock *>                  final_dist_dryrock(n_vintages, NULL);

    if (errTxt == "") {
      // Reuse alpha, but now it is sorted
      for(int i=0; i<n_solids; i++)
        alpha[i] = solid_alpha[i];
      for(int i=0; i<n_dryrocks; i++)
        alpha[i+n_solids] = dryrock_alpha[i];

      for(int i=0; i<n_vintages; i++)
        final_dist_dryrock[i] = new DistributionsDryRockMixOfDryRockAndSolid(distr_dryrock[i],
                                                                             distr_solid[i],
                                                                             dryrock_volume_fractions[i],
                                                                             solid_volume_fractions[i],
                                                                             mix_method,
                                                                             alpha);


      for(int i=0; i<n_vintages; i++) {
        for(size_t s=0; s<distr_solid[i].size(); s++)
          delete distr_solid[i][s];

        for(size_t s=0; s<distr_dryrock[i].size(); s++)
          delete distr_dryrock[i][s];

        for(size_t s=0; s<solid_volume_fractions[i].size(); s++) {
          if(solid_volume_fractions[i][s] != NULL) {
            if(solid_volume_fractions[i][s]->GetIsShared() == false)
              delete solid_volume_fractions[i][s];
          }
        }

        for(size_t s=0; s<dryrock_volume_fractions[i].size(); s++) {
          if(dryrock_volume_fractions[i][s] != NULL) {
            if(dryrock_volume_fractions[i][s]->GetIsShared() == false)
              delete dryrock_volume_fractions[i][s];
          }
        }
      }

    }

    return(final_dist_dryrock);
  } // end if mix_dryrock == true && mix_solid == true
  else {
    errTxt += "Not valid mix of constituents for creating a dry rock\n";
    return std::vector<DistributionsDryRock *>(n_vintages, NULL);
  }
}
//----------------------------------------------------------------------------------//

TabulatedVelocityDryRockStorage::TabulatedVelocityDryRockStorage(std::vector<DistributionWithTrendStorage *> vp,
                                                                 std::vector<DistributionWithTrendStorage *> vs,
                                                                 std::vector<DistributionWithTrendStorage *> density,
                                                                 std::vector<DistributionWithTrendStorage *> correlation_vp_vs,
                                                                 std::vector<DistributionWithTrendStorage *> correlation_vp_density,
                                                                 std::vector<DistributionWithTrendStorage *> correlation_vs_density,
                                                                 std::vector<DistributionWithTrendStorage *> total_porosity,
                                                                 std::vector<DistributionWithTrendStorage *> mineral_k)
: vp_(vp),
  vs_(vs),
  density_(density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density),
  total_porosity_(total_porosity),
  mineral_k_(mineral_k)
{
}

TabulatedVelocityDryRockStorage::~TabulatedVelocityDryRockStorage()
{
  if(vp_[0]->GetIsShared() == false)
    delete vp_[0];

  if(vs_[0]->GetIsShared() == false)
    delete vs_[0];

  if(density_[0]->GetIsShared() == false)
    delete density_[0];

  if(correlation_vp_vs_[0]->GetIsShared() == false)
    delete correlation_vp_vs_[0];

  if(correlation_vp_density_[0]->GetIsShared() == false)
    delete correlation_vp_density_[0];

  if(correlation_vs_density_[0]->GetIsShared() == false)
    delete correlation_vs_density_[0];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];

  if(mineral_k_[0]->GetIsShared() == false)
    delete mineral_k_[0];
}

std::vector<DistributionsDryRock *>
TabulatedVelocityDryRockStorage::GenerateDistributionsDryRock(const int                                                 & n_vintages,
                                                              const std::string                                         & path,
                                                              const std::vector<std::string>                            & trend_cube_parameters,
                                                              const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                              const std::map<std::string, DistributionsDryRockStorage *>& /*model_dryrock_storage*/,
                                                              const std::map<std::string, DistributionsSolidStorage *>  & /*model_solid_storage*/,
                                                              std::string                                               & errTxt) const
{
  std::vector<double> alpha(3);
  alpha[0] = vp_[0]     ->GetOneYearCorrelation();
  alpha[1] = vs_[0]     ->GetOneYearCorrelation();
  alpha[2] = density_[0]->GetOneYearCorrelation();

  int n_vintages_vp         = static_cast<int>(vp_.size());
  int n_vintages_vs         = static_cast<int>(vs_.size());
  int n_vintages_density    = static_cast<int>(density_.size());
  int n_vintages_porosity   = static_cast<int>(total_porosity_.size());
  int n_vintages_mineral_k  = static_cast<int>(mineral_k_.size());
  int n_vintages_vp_vs      = static_cast<int>(correlation_vp_vs_.size());
  int n_vintages_vp_density = static_cast<int>(correlation_vs_density_.size());
  int n_vintages_vs_density = static_cast<int>(correlation_vs_density_.size());

  std::string tmpErrTxt = "";

  for(int i=0; i<n_vintages_vp; i++) {
    if(vp_[i]->GetEstimate() == true)
      tmpErrTxt += "Vp can not be estimated from wells for a <dry-rock>\n";
  }
  for(int i=0; i<n_vintages_vs; i++) {
    if(vs_[i]->GetEstimate() == true)
      tmpErrTxt += "Vs can not be estimated from wells for a <dry-rock>\n";
  }
  for(int i=0; i<n_vintages_density; i++) {
    if(density_[i]->GetEstimate() == true)
      tmpErrTxt += "Density can not be estimated from wells for a <dry-rock>\n";
  }

  std::vector<double> corr_vp_vs(n_vintages, 0);
  std::vector<double> corr_vp_density(n_vintages, 0);
  std::vector<double> corr_vs_density(n_vintages, 0);

  for(int i=0; i<n_vintages_vp_vs; i++) {
    if(correlation_vp_vs_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-vp-vs> can not be estimated from wells\n";

    else
      FindDoubleValueFromDistributionWithTrend(correlation_vp_vs_[i], "correlation", corr_vp_vs[i], errTxt);

    if(corr_vp_vs[i] > 1 || corr_vp_vs[i] < -1)
        errTxt += "<correlation-vp-vs> should be in the interval [-1,1] in the tabulated model\n";
  }

  for(int i=0; i<n_vintages_vp_density; i++) {
    if(correlation_vp_density_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-vp-density> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_vp_density_[i], "correlation", corr_vp_density[i], errTxt);

    if(corr_vp_density[i] > 1 || corr_vp_density[i] < -1)
        errTxt += "<correlation-vp-density> should be in the interval [-1,1] in the tabulated model\n";
  }

  for(int i=0; i<n_vintages_vs_density; i++) {
    if(correlation_vs_density_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-vs-density> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_vs_density_[i], "correlation", corr_vs_density[i], errTxt);

    if(corr_vs_density[i] > 1 || corr_vs_density[i] < -1)
        errTxt += "<correlation-vs-density> should be in the interval [-1,1] in the tabulated model\n";
  }

  const std::vector<std::vector<double> > dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;
  const int                               dummy_output_other = -999;

  std::vector<DistributionsDryRock *>  dist_dryrock(n_vintages, NULL);
  std::vector<DistributionWithTrend *> vp_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> vs_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> porosity_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> mineral_k_dist_with_trend(n_vintages, NULL);

  if(tmpErrTxt == "") {

    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_vp)
        vp_dist_with_trend[i] = vp_[i]->GenerateDistributionWithTrend(path,
                                                                      trend_cube_parameters,
                                                                      trend_cube_sampling,
                                                                      dummy_blocked_logs,
                                                                      dummy_s1,
                                                                      dummy_s2,
                                                                      dummy_output_other,
                                                                      "dummy",
                                                                      tmpErrTxt);
      else
        vp_dist_with_trend[i] = vp_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_vs)
        vs_dist_with_trend[i] = vs_[i]->GenerateDistributionWithTrend(path,
                                                                      trend_cube_parameters,
                                                                      trend_cube_sampling,
                                                                      dummy_blocked_logs,
                                                                      dummy_s1,
                                                                      dummy_s2,
                                                                      dummy_output_other,
                                                                      "dummy",
                                                                      tmpErrTxt);
      else
        vs_dist_with_trend[i] = vs_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_density)
        density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path,
                                                                                trend_cube_parameters,
                                                                                trend_cube_sampling,
                                                                                dummy_blocked_logs,
                                                                                dummy_s1,
                                                                                dummy_s2,
                                                                                dummy_output_other,
                                                                                "dummy",
                                                                                tmpErrTxt);
      else
        density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_mineral_k)
        mineral_k_dist_with_trend[i] = mineral_k_[i]->GenerateDistributionWithTrend(path,
                                                                                    trend_cube_parameters,
                                                                                    trend_cube_sampling,
                                                                                    dummy_blocked_logs,
                                                                                    dummy_s1,
                                                                                    dummy_s2,
                                                                                    dummy_output_other,
                                                                                    "dummy",
                                                                                    tmpErrTxt);
      else
        mineral_k_dist_with_trend[i] = mineral_k_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_porosity)
        porosity_dist_with_trend[i] = total_porosity_[i]->GenerateDistributionWithTrend(path,
                                                                                        trend_cube_parameters,
                                                                                        trend_cube_sampling,
                                                                                        dummy_blocked_logs,
                                                                                        dummy_s1,
                                                                                        dummy_s2,
                                                                                        dummy_output_other,
                                                                                        "dummy",
                                                                                        tmpErrTxt);
      else
        porosity_dist_with_trend[i] = porosity_dist_with_trend[i-1]->Clone();

      if(i >= n_vintages_vp_vs)
        corr_vp_vs[i] = corr_vp_vs[i-1];

      if(i >= n_vintages_vp_density)
        corr_vp_density[i] = corr_vp_density[i-1];

      if(i >= n_vintages_vs_density)
        corr_vs_density[i] = corr_vs_density[i-1];
    }

    for(int i=0; i<n_vintages; i++) {
      std::string corrErrTxt = "";
      CheckPositiveDefiniteCorrMatrix(corr_vp_vs[i], corr_vp_density[i], corr_vs_density[i], corrErrTxt);
      if(corrErrTxt != "") {
        if(n_vintages > 1)
          tmpErrTxt += "Vintage "+NRLib::ToString(i+1)+":";
        tmpErrTxt += corrErrTxt;
      }
    }
  }

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      DistributionsDryRock * dryrock = new DistributionsDryRockTabulated(vp_dist_with_trend[i],
                                                                         vs_dist_with_trend[i],
                                                                         density_dist_with_trend[i],
                                                                         mineral_k_dist_with_trend[i],
                                                                         porosity_dist_with_trend[i],
                                                                         corr_vp_vs[i],
                                                                         corr_vp_density[i],
                                                                         corr_vs_density[i],
                                                                         DEMTools::Velocity,
                                                                         alpha);

      dist_dryrock[i] = dryrock;
    }

    for(int i=0; i<n_vintages; i++) {
      if(vp_dist_with_trend[i]->GetIsShared() == false)
        delete vp_dist_with_trend[i];
      if(vs_dist_with_trend[i]->GetIsShared() == false)
        delete vs_dist_with_trend[i];
      if(density_dist_with_trend[i]->GetIsShared() == false)
        delete density_dist_with_trend[i];
      if(porosity_dist_with_trend[i]->GetIsShared() == false)
        delete porosity_dist_with_trend[i];
      if(mineral_k_dist_with_trend[i]->GetIsShared() == false)
        delete mineral_k_dist_with_trend[i];
    }
  }

  else {
    errTxt += "\nProblems with the Tabulated rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(dist_dryrock);
}

//----------------------------------------------------------------------------------//
TabulatedModulusDryRockStorage::TabulatedModulusDryRockStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                                                               std::vector<DistributionWithTrendStorage *> shear_modulus,
                                                               std::vector<DistributionWithTrendStorage *> density,
                                                               std::vector<DistributionWithTrendStorage *> correlation_bulk_shear,
                                                               std::vector<DistributionWithTrendStorage *> correlation_bulk_density,
                                                               std::vector<DistributionWithTrendStorage *> correlation_shear_density,
                                                               std::vector<DistributionWithTrendStorage *> total_porosity,
                                                               std::vector<DistributionWithTrendStorage *> mineral_k)
: bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  density_(density),
  correlation_bulk_shear_(correlation_bulk_shear),
  correlation_bulk_density_(correlation_bulk_density),
  correlation_shear_density_(correlation_shear_density),
  total_porosity_(total_porosity),
  mineral_k_(mineral_k)
{
}

TabulatedModulusDryRockStorage::~TabulatedModulusDryRockStorage()
{

 if(bulk_modulus_[0]->GetIsShared() == false)
    delete bulk_modulus_[0];

  if(shear_modulus_[0]->GetIsShared() == false)
    delete shear_modulus_[0];

  if(density_[0]->GetIsShared() == false)
    delete density_[0];

  if(correlation_bulk_shear_[0]->GetIsShared() == false)
    delete correlation_bulk_shear_[0];

  if(correlation_bulk_density_[0]->GetIsShared() == false)
    delete correlation_bulk_density_[0];

  if(correlation_shear_density_[0]->GetIsShared() == false)
    delete correlation_shear_density_[0];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];

  if(mineral_k_[0]->GetIsShared() == false)
    delete mineral_k_[0];
}

std::vector<DistributionsDryRock *>
TabulatedModulusDryRockStorage::GenerateDistributionsDryRock(const int                                                 & n_vintages,
                                                             const std::string                                         & path,
                                                             const std::vector<std::string>                            & trend_cube_parameters,
                                                             const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                             const std::map<std::string, DistributionsDryRockStorage *>& /*model_dryrock_storage*/,
                                                             const std::map<std::string, DistributionsSolidStorage *>  & /*model_solid_storage*/,
                                                             std::string                                               & errTxt) const
{
  std::vector<double> alpha(3);
  alpha[0] = bulk_modulus_[0] ->GetOneYearCorrelation();
  alpha[1] = shear_modulus_[0]->GetOneYearCorrelation();
  alpha[2] = density_[0]      ->GetOneYearCorrelation();

  int n_vintages_bulk          = static_cast<int>(bulk_modulus_.size());
  int n_vintages_shear         = static_cast<int>(shear_modulus_.size());
  int n_vintages_density       = static_cast<int>(density_.size());
  int n_vintages_porosity      = static_cast<int>(total_porosity_.size());
  int n_vintages_mineral_k     = static_cast<int>(mineral_k_.size());
  int n_vintages_bulk_shear    = static_cast<int>(correlation_bulk_shear_.size());
  int n_vintages_bulk_density  = static_cast<int>(correlation_bulk_density_.size());
  int n_vintages_shear_density = static_cast<int>(correlation_shear_density_.size());

  std::string tmpErrTxt = "";

  for(int i=0; i<n_vintages_bulk; i++) {
    if(bulk_modulus_[i]->GetEstimate() == true)
      tmpErrTxt += "Bulk modulus can not be estimated from wells for a <dry-rock>\n";
  }
  for(int i=0; i<n_vintages_shear; i++) {
    if(shear_modulus_[i]->GetEstimate() == true)
      tmpErrTxt += "Shear modulus can not be estimated from wells for a <dry-rock>\n";
  }
  for(int i=0; i<n_vintages_density; i++) {
    if(density_[i]->GetEstimate() == true)
      tmpErrTxt += "Density can not be estimated from wells for a <dry-rock>\n";
  }

  std::vector<double> corr_bulk_shear(n_vintages, 0);
  std::vector<double> corr_bulk_density(n_vintages, 0);
  std::vector<double> corr_shear_density(n_vintages, 0);

  for(int i=0; i<n_vintages_bulk_shear; i++) {
    if(correlation_bulk_shear_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-bulk-shear> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_bulk_shear_[i], "correlation", corr_bulk_shear[i], errTxt);

    if(corr_bulk_shear[i] > 1 || corr_bulk_shear[i] < -1)
        errTxt += "<correlation-bulk-shear> should be in the interval [-1,1] in the tabulated model\n";
  }

  for(int i=0; i<n_vintages_bulk_density; i++) {
    if(correlation_bulk_density_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-bulk-density> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_bulk_density_[i], "correlation", corr_bulk_density[i], errTxt);

    if(corr_bulk_density[i] > 1 || corr_bulk_density[i] < -1)
        errTxt += "<correlation-bulk-density> should be in the interval [-1,1] in the tabulated model\n";
  }

  for(int i=0; i<n_vintages_shear_density; i++) {
    if(correlation_shear_density_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-shear-density> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_shear_density_[i], "correlation", corr_shear_density[i], errTxt);

    if(corr_shear_density[i] > 1 || corr_shear_density[i] < -1)
        errTxt += "<correlation-shear-density> should be in the interval [-1,1] in the tabulated model\n";
  }

  const std::vector<std::vector<double> > dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;
  const int                               dummy_output_other = -999;

  std::vector<DistributionsDryRock *>  dist_dryrock(n_vintages, NULL);
  std::vector<DistributionWithTrend *> bulk_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> shear_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> porosity_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> mineral_k_dist_with_trend(n_vintages, NULL);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_bulk)
        bulk_dist_with_trend[i] = bulk_modulus_[i]->GenerateDistributionWithTrend(path,
                                                                                  trend_cube_parameters,
                                                                                  trend_cube_sampling,
                                                                                  dummy_blocked_logs,
                                                                                  dummy_s1,
                                                                                  dummy_s2,
                                                                                  dummy_output_other,
                                                                                  "dummy",
                                                                                  tmpErrTxt);
      else
        bulk_dist_with_trend[i] = bulk_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_shear)
        shear_dist_with_trend[i] = shear_modulus_[i]->GenerateDistributionWithTrend(path,
                                                                                    trend_cube_parameters,
                                                                                    trend_cube_sampling,
                                                                                    dummy_blocked_logs,
                                                                                    dummy_s1,
                                                                                    dummy_s2,
                                                                                    dummy_output_other,
                                                                                    "dummy",
                                                                                    tmpErrTxt);
      else
        shear_dist_with_trend[i] = shear_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_density)
        density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path,
                                                                                trend_cube_parameters,
                                                                                trend_cube_sampling,
                                                                                dummy_blocked_logs,
                                                                                dummy_s1,
                                                                                dummy_s2,
                                                                                dummy_output_other,
                                                                                "dummy",
                                                                                tmpErrTxt);
      else
        density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_porosity)
        porosity_dist_with_trend[i] = total_porosity_[i]->GenerateDistributionWithTrend(path,
                                                                                        trend_cube_parameters,
                                                                                        trend_cube_sampling,
                                                                                        dummy_blocked_logs,
                                                                                        dummy_s1,
                                                                                        dummy_s2,
                                                                                        dummy_output_other,
                                                                                        "dummy",
                                                                                        tmpErrTxt);
      else
        porosity_dist_with_trend[i] = porosity_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_mineral_k)
        mineral_k_dist_with_trend[i] = mineral_k_[i]->GenerateDistributionWithTrend(path,
                                                                                    trend_cube_parameters,
                                                                                    trend_cube_sampling,
                                                                                    dummy_blocked_logs,
                                                                                    dummy_s1,
                                                                                    dummy_s2,
                                                                                    dummy_output_other,
                                                                                    "dummy",
                                                                                    tmpErrTxt);
      else
        mineral_k_dist_with_trend[i] = mineral_k_dist_with_trend[i-1]->Clone();


      double lower_mega = 1.0e+5; //Finn grenser fra modelsettings
      double upper_mega = 1.0e+8;
      double test_bulk  = bulk_dist_with_trend[0]->ReSample(0,0);
      double test_shear = shear_dist_with_trend[0]->ReSample(0,0);
      if(test_bulk < lower_mega || test_bulk > upper_mega)
        tmpErrTxt += "Bulk modulus need to be given in kPa\n";
      if(test_shear < lower_mega || test_shear > upper_mega)
        tmpErrTxt += "Shear modulus need to be given in kPa\n";

      if(i >= n_vintages_bulk_shear)
        corr_bulk_shear[i] = corr_bulk_shear[i-1];

      if(i >= n_vintages_bulk_density)
        corr_bulk_density[i] = corr_bulk_density[i-1];

      if(i >= n_vintages_shear_density)
        corr_shear_density[i] = corr_shear_density[i-1];
    }


    for(int i=0; i<n_vintages; i++) {
      std::string corrErrTxt = "";
      CheckPositiveDefiniteCorrMatrix(corr_bulk_shear[i], corr_bulk_density[i], corr_shear_density[i], corrErrTxt);
      if(corrErrTxt != "") {
        if(n_vintages > 1)
          tmpErrTxt += "Vintage "+NRLib::ToString(i+1)+":";
        tmpErrTxt += corrErrTxt;
      }
    }
  }

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      DistributionsDryRock * dryrock = new DistributionsDryRockTabulated(bulk_dist_with_trend[i],
                                                                         shear_dist_with_trend[i],
                                                                         density_dist_with_trend[i],
                                                                         mineral_k_dist_with_trend[i],
                                                                         porosity_dist_with_trend[i],
                                                                         corr_bulk_shear[i],
                                                                         corr_bulk_density[i],
                                                                         corr_shear_density[i],
                                                                         DEMTools::Modulus,
                                                                         alpha);

      dist_dryrock[i] = dryrock;
    }

    for(int i=0; i<n_vintages; i++) {
      if(bulk_dist_with_trend[i]->GetIsShared() == false)
        delete bulk_dist_with_trend[i];
      if(shear_dist_with_trend[i]->GetIsShared() == false)
        delete shear_dist_with_trend[i];
      if(density_dist_with_trend[i]->GetIsShared() == false)
        delete density_dist_with_trend[i];
      if(porosity_dist_with_trend[i]->GetIsShared() == false)
        delete porosity_dist_with_trend[i];
      if(mineral_k_dist_with_trend[i]->GetIsShared() == false)
        delete mineral_k_dist_with_trend[i];
    }
  }
  else {
    errTxt += "\nProblems with the Tabulated rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(dist_dryrock);
}

//----------------------------------------------------------------------------------//
ReussDryRockStorage::ReussDryRockStorage(std::vector<std::string >                                 constituent_label,
                                         std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussDryRockStorage::~ReussDryRockStorage()
{
  for (int i = 0; i < static_cast<int>(constituent_volume_fraction_[0].size()); i++) {
    if (constituent_volume_fraction_[0][i] != NULL) {
      if (constituent_volume_fraction_[0][i]->GetIsShared() == false)
        delete constituent_volume_fraction_[0][i];
    }
  }

}

std::vector<DistributionsDryRock *>
ReussDryRockStorage::GenerateDistributionsDryRock(const int                                                 & n_vintages,
                                                  const std::string                                         & path,
                                                  const std::vector<std::string>                            & trend_cube_parameters,
                                                  const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                  const std::map<std::string, DistributionsDryRockStorage *>& model_dryrock_storage,
                                                  const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                  std::string                                               & errTxt) const
{
  std::vector<DistributionsDryRock *> dryrock;

  std::string tmpErrTxt = "";

  for(size_t i=0; i<constituent_volume_fraction_.size(); i++) {
    for(size_t j=0; j<constituent_volume_fraction_[i].size(); j++) {
      if(constituent_volume_fraction_[i][j] != NULL && constituent_volume_fraction_[i][j]->GetEstimate() == true) {
        tmpErrTxt += "The volume fractions can not be estimated from wells\n";
        break;
      }
    }
    if(tmpErrTxt != "")
      break;
  }

  if(tmpErrTxt == "")
    dryrock = CreateDistributionsDryRockMix(n_vintages,
                                            path,
                                            trend_cube_parameters,
                                            trend_cube_sampling,
                                            model_dryrock_storage,
                                            model_solid_storage,
                                            constituent_label_,
                                            constituent_volume_fraction_,
                                            DEMTools::Reuss,
                                            tmpErrTxt);

  else {
    errTxt += "\nProblems with the Reuss rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(dryrock);
}

//----------------------------------------------------------------------------------//
VoigtDryRockStorage::VoigtDryRockStorage(std::vector<std::string>                                  constituent_label,
                                         std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtDryRockStorage::~VoigtDryRockStorage()
{
  for (int i = 0; i < static_cast<int>(constituent_volume_fraction_[0].size()); i++) {
    if (constituent_volume_fraction_[0][i] != NULL) {
      if (constituent_volume_fraction_[0][i]->GetIsShared() == false)
        delete constituent_volume_fraction_[0][i];
    }
  }
}

std::vector<DistributionsDryRock *>
VoigtDryRockStorage::GenerateDistributionsDryRock(const int                                                 & n_vintages,
                                                  const std::string                                         & path,
                                                  const std::vector<std::string>                            & trend_cube_parameters,
                                                  const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                  const std::map<std::string, DistributionsDryRockStorage *>& model_dryrock_storage,
                                                  const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                  std::string                                               & errTxt) const
{
  std::vector<DistributionsDryRock *> dryrock;

  std::string tmpErrTxt = "";

  for(size_t i=0; i<constituent_volume_fraction_.size(); i++) {
    for(size_t j=0; j<constituent_volume_fraction_[i].size(); j++) {
      if(constituent_volume_fraction_[i][j] != NULL && constituent_volume_fraction_[i][j]->GetEstimate() == true) {
        tmpErrTxt += "The volume fractions can not be estimated from wells\n";
        break;
      }
    }
    if(tmpErrTxt != "")
      break;
  }

  if(tmpErrTxt == "")
    dryrock = CreateDistributionsDryRockMix(n_vintages,
                                            path,
                                            trend_cube_parameters,
                                            trend_cube_sampling,
                                            model_dryrock_storage,
                                            model_solid_storage,
                                            constituent_label_,
                                            constituent_volume_fraction_,
                                            DEMTools::Voigt,
                                            tmpErrTxt);

  else {
    errTxt += "\nProblems with the Voigt rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(dryrock);
}

//----------------------------------------------------------------------------------//
HillDryRockStorage::HillDryRockStorage(std::vector<std::string>                                  constituent_label,
                                       std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillDryRockStorage::~HillDryRockStorage()
{
  for (int i = 0; i < static_cast<int>(constituent_volume_fraction_[0].size()); i++) {
    if (constituent_volume_fraction_[0][i] != NULL) {
      if (constituent_volume_fraction_[0][i]->GetIsShared() == false)
        delete constituent_volume_fraction_[0][i];
    }
  }
}

std::vector<DistributionsDryRock *>
HillDryRockStorage::GenerateDistributionsDryRock(const int                                                 & n_vintages,
                                                 const std::string                                         & path,
                                                 const std::vector<std::string>                            & trend_cube_parameters,
                                                 const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                 const std::map<std::string, DistributionsDryRockStorage *>& model_dryrock_storage,
                                                 const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                 std::string                                               & errTxt) const
{
  std::vector<DistributionsDryRock *> dryrock;

  std::string tmpErrTxt = "";

  for(size_t i=0; i<constituent_volume_fraction_.size(); i++) {
    for(size_t j=0; j<constituent_volume_fraction_[i].size(); j++) {
      if(constituent_volume_fraction_[i][j] != NULL && constituent_volume_fraction_[i][j]->GetEstimate() == true) {
        tmpErrTxt += "The volume fractions can not be estimated from wells\n";
        break;
      }
    }
    if(tmpErrTxt != "")
      break;
  }

  if(tmpErrTxt == "")
    dryrock = CreateDistributionsDryRockMix(n_vintages,
                                            path,
                                            trend_cube_parameters,
                                            trend_cube_sampling,
                                            model_dryrock_storage,
                                            model_solid_storage,
                                            constituent_label_,
                                            constituent_volume_fraction_,
                                            DEMTools::Hill,
                                            tmpErrTxt);

  else {
    errTxt += "\nProblems with the Hill rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(dryrock);
}

//----------------------------------------------------------------------------------//

DEMDryRockStorage::DEMDryRockStorage(std::string                                               host_label,
                                     std::vector<DistributionWithTrendStorage *>               host_volume_fraction,
                                     std::vector<std::string>                                  inclusion_label,
                                     std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction,
                                     std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio)
{
}

DEMDryRockStorage::~DEMDryRockStorage()
{
  if(host_volume_fraction_[0] != NULL) {
    if(host_volume_fraction_[0]->GetIsShared() == false)
      delete host_volume_fraction_[0];
  }

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_[0].size()); i++) {
    if(inclusion_volume_fraction_[0][i] != NULL) {
      if(inclusion_volume_fraction_[0][i]->GetIsShared() == false)
        delete inclusion_volume_fraction_[0][i];
    }
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_[0].size()); i++) {
    if(inclusion_aspect_ratio_[0][i]->GetIsShared() == false)
      delete inclusion_aspect_ratio_[0][i];
  }
}

std::vector<DistributionsDryRock *>
DEMDryRockStorage::GenerateDistributionsDryRock(const int                                                 & n_vintages,
                                                const std::string                                         & path,
                                                const std::vector<std::string>                            & trend_cube_parameters,
                                                const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                const std::map<std::string, DistributionsDryRockStorage *>& model_dryrock_storage,
                                                const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                std::string                                               & errTxt) const
{
  std::string tmpErrTxt = "";

  // Remember: Host info is included first in constituent vectors
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());
  int n_constituents = n_inclusions + 1;

  std::vector<std::vector<DistributionWithTrendStorage *> > volume_fractions(n_constituents);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=0; i<n_inclusions; i++)
    volume_fractions[i+1] = inclusion_volume_fraction_[i];

  std::vector<int> n_vintages_aspect(n_constituents);
  for(int i=0; i<n_inclusions; i++)
    n_vintages_aspect[i] = static_cast<int>(inclusion_aspect_ratio_[i].size());

  std::vector<int> n_vintages_volume(n_constituents);
  for(int i=0; i<n_constituents; i++)
    n_vintages_volume[i] = static_cast<int>(volume_fractions[i].size());

  // Order in alpha: aspect_ratios, host_volume_fraction, inclusion_volume_fractions
  std::vector<double> alpha(n_inclusions + n_constituents);

  for(int i=0; i<n_inclusions; i++)
    alpha[i] = inclusion_aspect_ratio_[i][0]->GetOneYearCorrelation();

  for(int i=0; i<n_constituents; i++) {
    if(volume_fractions[0][i] != NULL)
      alpha[i+n_inclusions] = volume_fractions[i][0]->GetOneYearCorrelation();
    else
      alpha[i+n_inclusions] = 1;
  }

  //Read host label
  std::vector<DistributionsDryRock *> final_distr_dryrock (n_vintages);
  std::vector<DistributionsDryRock *> distr_dryrock;

  distr_dryrock = ReadDryRock(n_vintages,
                              host_label_,
                              path,
                              trend_cube_parameters,
                              trend_cube_sampling,
                              model_dryrock_storage,
                              model_solid_storage,
                              tmpErrTxt);

  for(int i=0; i<n_vintages; i++) {
    if(i < static_cast<int>(distr_dryrock.size()))
      final_distr_dryrock[i] = distr_dryrock[i];
    else
      final_distr_dryrock[i] = final_distr_dryrock[i-1]->Clone();
  }

  //Read inclusion label
  std::vector<std::vector<DistributionsDryRock *> > final_distr_dryrock_inc(n_vintages);
  for(int i=0; i<n_vintages; i++)
    final_distr_dryrock_inc[i].resize(n_inclusions, NULL);

  for (int s = 0; s < n_inclusions; s++) {
    std::vector<DistributionsDryRock *> distr_dryrock_all_vintages = ReadDryRock(n_vintages,
                                                                           inclusion_label_[s],
                                                                           path,
                                                                           trend_cube_parameters,
                                                                           trend_cube_sampling,
                                                                           model_dryrock_storage,
                                                                           model_solid_storage,
                                                                           tmpErrTxt);

    for(int i=0; i<n_vintages; i++) {
      if(i < static_cast<int>(distr_dryrock_all_vintages.size()))
        final_distr_dryrock_inc[i][s] = distr_dryrock_all_vintages[i];
      else
        final_distr_dryrock_inc[i][s] = final_distr_dryrock_inc[i-1][s]->Clone();
    }
  }

  const std::vector<std::vector<double> > dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;
  const int                               dummy_output_other = -999;

  std::vector<DistributionsDryRock *>                  final_dist_dryrock(n_vintages, NULL);
  std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);
  std::vector<std::vector<DistributionWithTrend *> > all_aspect_ratios(n_vintages);

  for(size_t i=0; i<volume_fractions.size(); i++) {
    for(size_t j=0; j<volume_fractions[i].size(); j++) {
      if(volume_fractions[i][j] != NULL && volume_fractions[i][j]->GetEstimate() == true)
        tmpErrTxt += "Volume fractions can not be estimated from wells\n";
    }
  }
  for(size_t i=0; i<inclusion_aspect_ratio_.size(); i++) {
    for(size_t j=0; j<inclusion_aspect_ratio_[i].size(); j++) {
      if(inclusion_aspect_ratio_[i][j]->GetEstimate() == true)
        tmpErrTxt += "Aspect ratios can not be estimated from wells\n";
    }
  }

  for(int i=0; i<n_vintages; i++) {
    all_volume_fractions[i].resize(n_constituents, NULL);
    all_aspect_ratios[i].resize(n_inclusions, NULL);
  }

  for(int i=0; i<n_constituents; i++)
    CheckValuesInZeroOne(volume_fractions[i], "volume-fraction", path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, tmpErrTxt);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {

      for (int s = 0; s < n_inclusions; s++) {

        if(i < n_vintages_aspect[s]) {
          if(inclusion_aspect_ratio_[s][i] != NULL)
            all_aspect_ratios[i][s] = inclusion_aspect_ratio_[s][i]->GenerateDistributionWithTrend(path,
                                                                                                   trend_cube_parameters,
                                                                                                   trend_cube_sampling,
                                                                                                   dummy_blocked_logs,
                                                                                                   dummy_s1,
                                                                                                   dummy_s2,
                                                                                                   dummy_output_other,
                                                                                                   "dummy",
                                                                                                   tmpErrTxt);
        }
        else
          all_aspect_ratios[i][s] = all_aspect_ratios[i-1][s]->Clone();
      }

      for (int s = 0; s < n_constituents; s++) {

        if(i < n_vintages_volume[s]) {
          if(volume_fractions[s][i] != NULL)
            all_volume_fractions[i][s] = volume_fractions[s][i]->GenerateDistributionWithTrend(path,
                                                                                               trend_cube_parameters,
                                                                                               trend_cube_sampling,
                                                                                               dummy_blocked_logs,
                                                                                               dummy_s1,
                                                                                               dummy_s2,
                                                                                               dummy_output_other,
                                                                                               "dummy",
                                                                                               tmpErrTxt);
        }
        else
          all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
      }

      CheckVolumeConsistency(all_volume_fractions[i], tmpErrTxt);

    }
  }

  if (tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++)
      final_dist_dryrock[i] = new DistributionsDryRockDEM(final_distr_dryrock[i],
                                                      final_distr_dryrock_inc[i],
                                                      all_aspect_ratios[i],
                                                      all_volume_fractions[i],
                                                      alpha);

    for(int i=0; i<n_vintages; i++) {
      delete final_distr_dryrock[i];

      for(size_t s=0; s<final_distr_dryrock_inc[i].size(); s++)
        delete final_distr_dryrock_inc[i][s];

      for(size_t s=0; s<all_aspect_ratios[i].size(); s++) {
        if(all_aspect_ratios[i][s]->GetIsShared() == false)
          delete all_aspect_ratios[i][s];
      }

      for(size_t s=0; s<all_volume_fractions[i].size(); s++) {
        if(all_volume_fractions[i][s] != NULL) {
          if(all_volume_fractions[i][s]->GetIsShared() == false)
            delete all_volume_fractions[i][s];
        }
      }
    }
  }
  else {
    errTxt += "\nProblems with the DEM rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(final_dist_dryrock);
}

//----------------------------------------------------------------------------------------------------------------------------------//

WaltonDryRockStorage::WaltonDryRockStorage(const std::string                           & solid_label,
                                           std::vector<DistributionWithTrendStorage *> & distr_friction_weight,
                                           std::vector<DistributionWithTrendStorage *> & distr_pressure,
                                           std::vector<DistributionWithTrendStorage *> & distr_porosity,
                                           std::vector<DistributionWithTrendStorage *> & distr_coord_number)
: solid_label_(solid_label),
  distr_friction_weight_(distr_friction_weight),
  distr_pressure_(distr_pressure),
  distr_porosity_(distr_porosity),
  distr_coord_number_(distr_coord_number)
{
}

WaltonDryRockStorage::~WaltonDryRockStorage()
{
  for(size_t i = 0; i < distr_coord_number_.size(); i++) {
    if(distr_coord_number_[i]->GetIsShared() == false)
      delete distr_coord_number_[i];
  }

  for(size_t i = 0; i < distr_porosity_.size(); i++) {
    if(distr_porosity_[i]->GetIsShared() == false)
      delete distr_porosity_[i];
  }

  for(size_t i = 0; i < distr_pressure_.size(); i++) {
    if(distr_pressure_[i]->GetIsShared() == false)
      delete distr_pressure_[i];
  }

  for(size_t i = 0; i < distr_friction_weight_.size(); i++) {
    if(distr_friction_weight_[i]->GetIsShared() == false)
      delete distr_friction_weight_[i];
  }

}

std::vector<DistributionsDryRock *>
WaltonDryRockStorage::
GenerateDistributionsDryRock(const int                                                 & n_vintages,
                              const std::string                                         & path,
                              const std::vector<std::string>                            & trend_cube_parameters,
                              const std::vector<std::vector<double> >                   & trend_cube_sampling,
                              const std::map<std::string, DistributionsDryRockStorage *>& /*model_dryrock_storage*/,
                              const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                              std::string                                               & errTxt) const
{
  std::string tmpErrTxt = "";

  std::vector<double> alpha(4);
  // order: friction_weight, pressure, porosity, coord_number
  alpha[0] = distr_friction_weight_[0]->GetOneYearCorrelation();
  alpha[1] = distr_pressure_[0]->GetOneYearCorrelation();
  alpha[2] = distr_porosity_[0]->GetOneYearCorrelation();
  alpha[3] = distr_coord_number_[0]->GetOneYearCorrelation();

  for(size_t i=0; i<distr_friction_weight_.size(); i++) {
    if(distr_friction_weight_[i]->GetEstimate() == true)
      tmpErrTxt += "No slip can not be estimated from wells for a <dry-rock>\n";
  }
  for(size_t i=0; i<distr_pressure_.size(); i++) {
    if(distr_pressure_[i]->GetEstimate() == true)
      tmpErrTxt += "Pressure can not be estimated from wells for a <dry-rock>\n";
  }
  for(size_t i=0; i<distr_porosity_.size(); i++) {
    if(distr_porosity_[i]->GetEstimate() == true)
      tmpErrTxt += "Porosity can not be estimatedfrom wells for a <dry-rock>\n";
  }
  for(size_t i=0; i<distr_coord_number_.size(); i++) {
    if(distr_coord_number_[i]->GetEstimate() == true)
      tmpErrTxt += "Coord nr can not be estimatedfrom wells for a <dry-rock>\n";
  }


  int n_vintages_friction_weight          = static_cast<int>(distr_friction_weight_.size());
  int n_vintages_pressure                 = static_cast<int>(distr_pressure_.size());
  int n_vintages_porosity                 = static_cast<int>(distr_porosity_.size());
  int n_vintages_coord_number             = static_cast<int>(distr_coord_number_.size());

  const std::vector<std::vector<double> > dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;
  const int                               dummy_output_other = -999;

  std::vector<DistributionsDryRock *>    dist_solid(n_vintages, NULL);
  std::vector<DistributionWithTrend *>   dwt_friction_weight(n_vintages, NULL);
  std::vector<DistributionWithTrend *>   dwt_pressure(n_vintages, NULL);
  std::vector<DistributionWithTrend *>   dwt_porosity(n_vintages, NULL);
  std::vector<DistributionWithTrend *>   dwt_coord_number(n_vintages, NULL);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_friction_weight)
        dwt_friction_weight[i] = distr_friction_weight_[i]->GenerateDistributionWithTrend(path,
                                                                                          trend_cube_parameters,
                                                                                          trend_cube_sampling,
                                                                                          dummy_blocked_logs,
                                                                                          dummy_s1,
                                                                                          dummy_s2,
                                                                                          dummy_output_other,
                                                                                          "dummy",
                                                                                          tmpErrTxt);
      else
        dwt_friction_weight[i] = dwt_friction_weight[i-1]->Clone();

      if(i < n_vintages_pressure)
        dwt_pressure[i] = distr_pressure_[i]->GenerateDistributionWithTrend(path,
                                                                            trend_cube_parameters,
                                                                            trend_cube_sampling,
                                                                            dummy_blocked_logs,
                                                                            dummy_s1,
                                                                            dummy_s2,
                                                                            dummy_output_other,
                                                                            "dummy",
                                                                            tmpErrTxt);
      else
        dwt_pressure[i] = dwt_pressure[i-1]->Clone();

      if(i < n_vintages_porosity)
        dwt_porosity[i] = distr_porosity_[i]->GenerateDistributionWithTrend(path,
                                                                            trend_cube_parameters,
                                                                            trend_cube_sampling,
                                                                            dummy_blocked_logs,
                                                                            dummy_s1,
                                                                            dummy_s2,
                                                                            dummy_output_other,
                                                                            "dummy",
                                                                            tmpErrTxt);
      else
        dwt_porosity[i] = dwt_porosity[i-1]->Clone();

      if(i < n_vintages_coord_number)
        dwt_coord_number[i] = distr_coord_number_[i]->GenerateDistributionWithTrend(path,
                                                                                    trend_cube_parameters,
                                                                                    trend_cube_sampling,
                                                                                    dummy_blocked_logs,
                                                                                    dummy_s1,
                                                                                    dummy_s2,
                                                                                    dummy_output_other,
                                                                                    "dummy",
                                                                                    tmpErrTxt);
      else
        dwt_coord_number[i] = dwt_coord_number[i-1]->Clone();
    }
  }

  std::vector<DistributionsSolid *> final_distr_solid(n_vintages);

  std::vector<DistributionsSolid *> distr_solid = ReadSolid(n_vintages,
                                                            solid_label_,
                                                            path,
                                                            trend_cube_parameters,
                                                            trend_cube_sampling,
                                                            model_solid_storage,
                                                            tmpErrTxt);

  int n_vintages_solid = static_cast<int>(distr_solid.size());

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_solid)
      final_distr_solid[i] = distr_solid[i];
    else
      final_distr_solid[i] = final_distr_solid[i-1]->Clone();
  }

  std::vector<DistributionsDryRock*> distr_dryrock(n_vintages, NULL);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      DistributionsDryRock * dryrock = new DistributionsDryRockWalton(distr_solid[i],
                                                                      dwt_friction_weight[i],
                                                                      dwt_pressure[i],
                                                                      dwt_porosity[i],
                                                                      dwt_coord_number[i],
                                                                      alpha);

      distr_dryrock[i] = dryrock;
    }

    //clean up
    for(int i=0; i<n_vintages; i++) {
      delete distr_solid[i];

      if (dwt_friction_weight[i]->GetIsShared() == false)
        delete dwt_friction_weight[i];

      if (dwt_pressure[i]->GetIsShared() == false)
        delete dwt_pressure[i];

      if (dwt_porosity[i]->GetIsShared() == false)
        delete dwt_porosity[i];

      if (dwt_coord_number[i]->GetIsShared() == false)
        delete dwt_coord_number[i];

    }
  }

  else {
    errTxt += "\nProblems with the Walton rock physics model for <dry-rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(distr_dryrock);
}
