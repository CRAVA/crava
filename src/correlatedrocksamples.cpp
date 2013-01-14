#include "src/correlatedrocksamples.h"

#include "rplib/rock.h"

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
  int et_dummy, edi_dummy, dt;
  std::vector< std::vector<int> > delta_time(k_max);
  for (int k = 0; k < k_max; ++k){
    time_line.GetNextEvent(et_dummy, edi_dummy, dt);
    for (int l = k; l < k_max; ++l)
      delta_time[l].push_back(dt);
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
  for (int i = 0; i < i_max; ++i){
    rock[0][i] = dist_rock[0]->GenerateSample(trend_params_dummy);
    rock[0][i]->GetSeismicParams(m[0][i][0], m[0][i][1], m[0][i][2]);
    std::vector< Rock * > rock_seen(1, rock[0][i]);
    for (int k = 1; k < k_max; ++k){
      rock[k][i] = dist_rock[k]->EvolveSample(delta_time[k][k], *(rock[k-1][i])); // delta_time info also for the rock to be found.
      rock[k][i]->GetSeismicParams(m[k][i][0], m[k][i][1], m[k][i][2]);
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

