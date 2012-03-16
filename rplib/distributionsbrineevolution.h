#ifndef DISTRIBUTIONSBRINEEVOLUTION_H
#define DISTRIBUTIONSBRINEEVOLUTION_H


// Parallel classes are Brine and DistributionsBrineT0.
class DistributionsBrineEvolution : public DistributionsFluidEvolution {
public:

  DistributionsBrineEvolution() : DistributionsFluidEvolution(){}

  virtual ~DistributionsBrineEvolution(){}

  virtual void GetParameters(std::vector<double> & param_fluid_evolve) const {
    param_fluid_evolve.resize(1);  //FAKE!
    param_fluid_evolve[0] = 1.0;   //FAKE!
  }

};

#endif
