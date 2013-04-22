#include "src/correlatedrocksamples.h"

#include "rplib/rock.h"
#include <nrlib/flens/nrlib_flens.hpp>


CorrelatedRockSamples::CorrelatedRockSamples()
{
}

CorrelatedRockSamples::~CorrelatedRockSamples()
{
}

std::vector< std::vector< std::vector<double> > >
CorrelatedRockSamples::CreateSamples(int                                     i_max,
                                     TimeLine                              & time_line,
                                     const std::vector<DistributionsRock*> & dist_rock)
{
  std::list<int> time;
  time_line.GetAllTimes(time);
  int k_max = static_cast<int>( time.size() );

  // Set up time vectors, the same for all sets of correlated samples.
  time_line.ReSet();
  int et_dummy, edi_dummy;
  double dt;
  std::vector< std::vector<double> > delta_time(k_max);
  for (int k = 0; k < k_max; ++k){
    time_line.GetNextEvent(et_dummy, edi_dummy, dt);// dt in days
    for (int l = k; l < k_max; ++l)
      delta_time[l].push_back(dt); // dt in years
  }

  // Set up data structures.
  // The order of indices is chosen to make extraction of all samples for a given time instance easy.
  std::vector< std::vector< std::vector<double> > > m(k_max);
  std::vector< std::vector< Rock * > > rock(k_max);
  for (int k = 0; k < k_max; ++k){
    m[k].resize(i_max);
    rock[k].resize(i_max);
    for (int i = 0; i < i_max; ++i)
      m[k][i].resize(3);
  }

  const std::vector<double> trend_params_dummy(2,0);
  // Finding the sets of correlated samples.
  // Each set of samples for a specific i are correlated in time.
  dist_rock[0]->SetResamplingLevel(DistributionWithTrend::Full);
  for (int i = 0; i < i_max; ++i){
    rock[0][i] = dist_rock[0]->GenerateSample(trend_params_dummy);
    rock[0][i]->GetSeismicParams(m[0][i][0], m[0][i][1], m[0][i][2]);
    m[0][i][0]=std::log(m[0][i][0]);
    m[0][i][1]=std::log(m[0][i][1]);
    m[0][i][2]=std::log(m[0][i][2]);
    std::vector< Rock * > rock_seen(1, rock[0][i]);
    for (int k = 1; k < k_max; ++k){
      dist_rock[k]->SetResamplingLevel(DistributionWithTrend::Full);
      rock[k][i] = dist_rock[k]->EvolveSample(delta_time[k][k], *(rock[k-1][i])); // delta_time info also for the rock to be found.
      rock[k][i]->GetSeismicParams(m[k][i][0], m[k][i][1], m[k][i][2]);
      m[k][i][0]=std::log(m[k][i][0]);
      m[k][i][1]=std::log(m[k][i][1]);
      m[k][i][2]=std::log(m[k][i][2]);
      rock_seen.push_back(rock[k][i]);
    }
  }

  // Must delete memory allocated by classes DistributionsRock and Rock.
  for (int k = 0; k < k_max; ++k){
    for (int i = 0; i < i_max; ++i){
      if (rock[k][i] != NULL)
        delete rock[k][i];
    }
  }
  return m;
}




std::vector< std::vector< std::vector<double> > >
CorrelatedRockSamples::CreateSamplesExtended(int                                     i_max,
                                             TimeLine                              & time_line,
                                             const std::vector<DistributionsRock*> & dist_rock)
{
  std::list<int> time;
  time_line.GetAllTimes(time);
  int k_max = static_cast<int>( time.size() );

  // Set up time vectors, the same for all sets of correlated samples.
  time_line.ReSet();
  int et_dummy, edi_dummy;
  double dt_year;
  std::vector< std::vector<double> > delta_time(k_max);
  for (int k = 0; k < k_max; ++k){
    time_line.GetNextEvent(et_dummy, edi_dummy, dt_year);
    for (int l = k; l < k_max; ++l)
      delta_time[l].push_back(dt_year);
  }

  int nReservoirVariables =  dist_rock[0]->GetNumberOfReservoirVariables();

  // Set up data structures.
  // The order of indices is chosen to make extraction of all samples for a given time instance easy.
  std::vector< std::vector< std::vector<double> > > m(k_max);
  std::vector< std::vector< Rock * > > rock(k_max);
  for (int k = 0; k < k_max; ++k){
    m[k].resize(i_max);
    rock[k].resize(i_max);
    for (int i = 0; i < i_max; ++i)
      m[k][i].resize(3+nReservoirVariables);
  }

  const std::vector<double> trend_params_dummy(2,0);
  std::vector<double>reservoirVariables(nReservoirVariables,0);
  // Finding the sets of correlated samples.
  // Each set of samples for a specific i are correlated in time.
  dist_rock[0]->SetResamplingLevel(DistributionWithTrend::Full);
  for (int i = 0; i < i_max; ++i){
    rock[0][i] = dist_rock[0]->GenerateSampleAndReservoirVariables(trend_params_dummy,reservoirVariables);
    rock[0][i]->GetSeismicParams(m[0][i][0], m[0][i][1], m[0][i][2]);
    m[0][i][0]=std::log(m[0][i][0]);
    m[0][i][1]=std::log(m[0][i][1]);
    m[0][i][2]=std::log(m[0][i][2]);
    for(int l =0;l< nReservoirVariables;l++)
        m[0][i][l+3]=reservoirVariables[l];

    std::vector< Rock * > rock_seen(1, rock[0][i]);
    for (int k = 1; k < k_max; ++k){
      dist_rock[k]->SetResamplingLevel(DistributionWithTrend::Full);
      rock[k][i] = dist_rock[k]->EvolveSampleAndReservoirVaribles(delta_time[k][k], *(rock[k-1][i]),reservoirVariables); // delta_time info also for the rock to be found.
      rock[k][i]->GetSeismicParams(m[k][i][0], m[k][i][1], m[k][i][2]);
      m[k][i][0]=std::log(m[k][i][0]);
      m[k][i][1]=std::log(m[k][i][1]);
      m[k][i][2]=std::log(m[k][i][2]);
      for(int l =0;l< nReservoirVariables;l++)
        m[k][i][l+3]=reservoirVariables[l];
      rock_seen.push_back(rock[k][i]);
    }
  }

  // Must delete memory allocated by classes DistributionsRock and Rock.
  for (int k = 0; k < k_max; ++k){
    for (int i = 0; i < i_max; ++i){
      if (rock[k][i] != NULL)
        delete rock[k][i];
    }
  }

  return m;
}


