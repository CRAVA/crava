#ifndef CORRELATEDROCKSAMPLES_H
#define CORRELATEDROCKSAMPLES_H

#include "src/timeline.h"
#include "rplib/distributionsrock.h"

// Class that creates K X I samples of seismic parameters.
// K = number of time steps.
// I = number of samples per time step.
// Each sample is a 3-dim vector [vp, vs, rho].
// Each set of samples for a specific i in [0:I-1] are correlated in time.

class CorrelatedRockSamples {
public:

  CorrelatedRockSamples();

  ~CorrelatedRockSamples();

  std::vector< std::vector< std::vector<double> > > CreateSamples(int                                             i_max,
                                                                  TimeLine                                      & time_line,
                                                                  const std::vector<DistributionsRock*>         & dist_rock);

  std::vector< std::vector< std::vector<double> > > CreateSamplesExtended(int                                     i_max,
                                                                          TimeLine                              & time_line,
                                                                          const std::vector<DistributionsRock*> & dist_rock);
};

#endif
