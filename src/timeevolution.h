#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>

class SeismicParametersHolder;
class TimeLine;
class DistributionsRockT0;
class DistributionsSaturation;
class DistributionsGeochemical;

class TimeEvolution
{
public:
  TimeEvolution(int number_of_timesteps,
                int i_max,
                TimeLine & time_line,
                const DistributionsRockT0 * dist_rock,
                const DistributionsSaturation * dist_sat,
                const DistributionsGeochemical * dist_geochem);
  void Evolve(int time_step, SeismicParametersHolder &m);

private:
  int number_of_timesteps_;

  std::vector< NRLib::Matrix> evolution_matrix_;
  std::vector< NRLib::Vector> delta_mu_tilde_;
  std::vector< NRLib::Matrix> delta_tilde_;

  void Split  (int fourier_point);
  void Predict(int fourier_point);
  void Merge  (int fourier_point);
  void SetUpEvolutionMatrix(std::vector< NRLib::Matrix> &evolution_matrix,
                            int i_max, TimeLine & time_line,
                            const DistributionsRockT0 * dist_rock,
                            const DistributionsSaturation * dist_sat,
                            const DistributionsGeochemical * dist_geochem);
  void SetUpDeltaMuTilde   (std::vector< NRLib::Vector> &delta_mu_tilde);
  void SetUpDeltaTilde     (std::vector< NRLib::Matrix> &delta_tilde);
};

#endif