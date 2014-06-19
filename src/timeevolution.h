#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>
#include <string>

class SeismicParametersHolder;
class State4D;
class TimeLine;
class DistributionsRock;
class FFTGrid;


class TimeEvolution
{
public:
  TimeEvolution() {}
  TimeEvolution(int                                     i_max,
                TimeLine                              & time_line,
                const std::vector<DistributionsRock*> & dist_rock);
  //void Split(const SeismicParametersHolder &m_combined, State4D & state4D);
  //void Evolve(int time_step, State4D & state4D);
  //void Merge(const State4D & state4D, SeismicParametersHolder &m_combined);

  NRLib::Matrix getEvolutionMatrix(int time_step)          const { return evolution_matrix_[time_step]; }
  NRLib::Vector getMeanCorrectionTerm(int time_step)       const { return mean_correction_term_[time_step];}
  NRLib::Matrix getCovarianceCorrectionTerm(int time_step) const { return cov_correction_term_[time_step]; }
  std::vector<std::vector<std::vector<double> > > returnCorrelatedSample(int                                      i_max,
                                                                         TimeLine                               & time_line,
                                                                         const std::vector<DistributionsRock*>  & dist_rock);
  NRLib::Vector computePriorMeanStaticAndDynamicLastTimeStep();
  NRLib::Matrix computePriorCovStaticAndDynamicLastTimeStep();
  void SetInitialMean(NRLib::Vector initialMean){initial_mean_=initialMean;}
  void SetInitialCov(NRLib::Matrix initialCov){initial_cov_=initialCov;}
  int  GetNTimSteps(){return number_of_timesteps_;}

private:
  int number_of_timesteps_;

  NRLib::Matrix initial_cov_;
  NRLib::Vector initial_mean_;
  std::vector< NRLib::Matrix> evolution_matrix_;
  std::vector< NRLib::Vector> mean_correction_term_;
  std::vector< NRLib::Matrix> cov_correction_term_;

  //Expand every 3-vector into a 6-vector of 3 static and 3 dynamic elements.
  void SplitSamplesStaticDynamic(std::vector<std::vector<std::vector<double> > > & m_ik) const;

  // Estimate time evolution matrices and correction term mean and covariance:
  void SetUpEvolutionMatrices(std::vector< NRLib::Matrix>          & evolution_matrix,
                             std::vector< NRLib::Matrix>           & cov_correction_term,
                             std::vector< NRLib::Vector>           & mean_correction_term,
                             int                                     i_max,
                             TimeLine                              & time_line,
                             const std::vector<DistributionsRock*> & dist_rock);

  // Adjusting diagonal of matrix to be inverted if necessary,
  // and adjusting the other matrices similarly to ensure block form of evolution
  // matrices and correction term covariance:
  void DoRobustInversion(NRLib::Matrix & SigmaInv,
                       //  NRLib::Matrix & Cov_mk_mk,
                       //  NRLib::Matrix & Cov_mk_mkm1,
                         const NRLib::Matrix  Cov_mkm1_mkm1,
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
