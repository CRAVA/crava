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
  std::vector< NRLib::Vector> mean_correction_term_;
  std::vector< NRLib::Matrix> cov_correction_term_;

  void Split(SeismicParametersHolder &m, SeismicParametersHolder &m_dynamic, SeismicParametersHolder &m_static);
  void Merge(SeismicParametersHolder &m, SeismicParametersHolder &m_dynamic, SeismicParametersHolder &m_static);

  void SetUpEvolutionMatrices(std::vector< NRLib::Matrix> &evolution_matrix,
                             std::vector< NRLib::Matrix> &cov_correction_term,
                             std::vector< NRLib::Vector> &mean_correction_term,
                             int i_max,
                             TimeLine & time_line,
                             const DistributionsRockT0 * dist_rock,
                             const DistributionsSaturation * dist_sat,
                             const DistributionsGeochemical * dist_geochem);
};

#endif
