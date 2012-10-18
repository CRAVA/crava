#include "rplib/distributionssolid.h"
#include "rplib/solid.h"

Solid * DistributionsSolid::EvolveSample(double        time,
                                         const Solid & solid) const
{
    const std::vector<double> trend(2);
    return UpdateSample(time, true, trend, &solid);
}
