#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>

class SeismicParametersHolder;

class TimeEvolution
{
public:
  TimeEvolution(int number_of_timesteps);
  void Evolve(int time_step, SeismicParametersHolder &m);

private:
  int number_of_timesteps_;

  std::vector< NRLib::Matrix> evolution_matrix_;
  std::vector< NRLib::Vector> delta_mu_tilde_;
  std::vector< NRLib::Matrix> delta_tilde_;

  void Split  (int fourier_point);
  void Predict(int fourier_point);
  void Merge  (int fourier_point);
  void SetUpEvolutionMatrix(std::vector< NRLib::Matrix> &evolution_matrix);
  void SetUpDeltaMuTilde   (std::vector< NRLib::Vector> &delta_mu_tilde);
  void SetUpDeltaTilde     (std::vector< NRLib::Matrix> &delta_tilde);
};

#endif
