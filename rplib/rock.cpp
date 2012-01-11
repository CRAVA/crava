#include "rplib/rock.h"

std::vector<Fluid *> Rock::fluid_(0);

Rock::Rock(const std::vector<double> & param, const std::vector<double> & saturation)
: param_(param), saturation_(saturation)
{
  assert(fluid_.size() == saturation_.size());
}


Rock::~Rock()
{
}
