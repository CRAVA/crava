#include "rplib/rock.h"


Rock::Rock(const std::vector<double> & param, const std::vector<double> & saturation) 
: param_(param), saturation_(saturation)
{
  assert(fluid_.size() == saturation_.size());
}


Rock::~Rock()
{
}
