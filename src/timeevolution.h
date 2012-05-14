#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>

class SeismicParametersHolder;
class State4D;
class TimeLine;
class DistributionsRock;
class FFTGrid;


class TimeEvolution
{
public:
  TimeEvolution(int i_max,
                TimeLine & time_line,
                const DistributionsRock       * dist_rock);
  //void Split(const SeismicParametersHolder &m_combined, State4D & state4D);
  void Evolve(int time_step, State4D & state4D);
  //void Merge(const State4D & state4D, SeismicParametersHolder &m_combined);

private:
  int number_of_timesteps_;

  std::vector< NRLib::Matrix> evolution_matrix_;
  std::vector< NRLib::Vector> mean_correction_term_;
  std::vector< NRLib::Matrix> cov_correction_term_;

  //Expand every 3-vector into a 6-vector of 3 static and 3 dynamic elements.
  void SplitSamplesStaticDynamic(std::vector<std::vector<std::vector<double> > > & m_ik);


  // Estimate time evolution matrices and correction term mean and covariance:
  void SetUpEvolutionMatrices(std::vector< NRLib::Matrix>   & evolution_matrix,
                             std::vector< NRLib::Matrix>    & cov_correction_term,
                             std::vector< NRLib::Vector>    & mean_correction_term,
                             int                              i_max,
                             TimeLine                       & time_line,
                             const DistributionsRock        * dist_rock);

  // Adjusting diagonal of matrix to be inverted if necessary,
  // and adjusting the other matrices similarly to ensure block form of evolution
  // matrices and correction term covariance:
  void DoRobustInversion(NRLib::Matrix & SigmaInv,
                         NRLib::Matrix & Cov_mk_mk,
                         NRLib::Matrix & Cov_mk_mkm1,
                         NRLib::Matrix & Cov_mkm1_mkm1,
                         int             dim,
                         double          adjustment_factor);
  bool DiagonalOK(const NRLib::Matrix & A, int dim);
  void FixMatrices(NRLib::Matrix & Cov_mk_mk,
                   NRLib::Matrix & Cov_mk_mkm1,
                   NRLib::Matrix & Cov_mkm1_mkm1,
                   int             dim,
                   double          adjustment_factor);

  // Removing numerical noise from known blocks:
  void AdjustMatrixAForm(NRLib::Matrix & A_k, int dim);
  void AdjustMatrixDeltaForm(NRLib::Matrix & delta_k, int dim);

  void PrintToScreen(NRLib::Matrix, int dim1, int dim2, std::string name);
  void PrintToScreen(NRLib::Vector, int dim, std::string name);
};

#endif
