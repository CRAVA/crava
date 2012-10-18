#include "rplib/distributionsfluid.h"
#include "rplib/fluid.h"

Fluid * DistributionsFluid::EvolveSample(double        time,
                                         const Fluid & fluid) const
{
    const std::vector<double> trend(2);
    return UpdateSample(time, true, trend, &fluid);
}
