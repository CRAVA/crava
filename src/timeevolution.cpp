#include <src/timeevolution.h>
#include <assert.h>
#include <float.h>

#include "nrlib/statistics/statistics.hpp"

#include "src/seismicparametersholder.h"
#include "src/state4d.h"
#include "src/fftgrid.h"
#include "src/timeline.h"
#include "src/correlatedrocksamples.h"
#include "src/tasklist.h"

#include "rplib/distributionsrock.h"


TimeEvolution::TimeEvolution(int i_max,
                             TimeLine & time_line,
                             const DistributionsRock        * dist_rock)
{
  std::list<int> time;
  time_line.GetAllTimes(time);
  number_of_timesteps_ = static_cast<int>( time.size() );
  SetUpEvolutionMatrices(evolution_matrix_, cov_correction_term_, mean_correction_term_, i_max, time_line, dist_rock);
}


void TimeEvolution::SetUpEvolutionMatrices(std::vector< NRLib::Matrix>    & evolution_matrix,
                                           std::vector< NRLib::Matrix>    & cov_correction_term,
                                           std::vector< NRLib::Vector>    & mean_correction_term,
                                           int                              i_max,
                                           TimeLine                       & time_line,
                                           const DistributionsRock        * dist_rock)
{
  // A note on variable names in this function:
  // _mk in the variable names refers to the seismic parameter m with subscript k, time index.
  // _mkm1 in the variable names refers to the seismic parameter m with subscript k-1.

  // A_k:        Denotes the time evolution matrix in the expression m_{k} = A_{k}*m_{k-1} + \delta m_{k}
  // delta_k:    Denotes the covariance of the correction term \delta m_{k}, denoted with symbol \delta_k
  // delta_mu_k: Denotes the mean of the correction term \delta m_{k}, denoted with symbol \delta \mu_{k}

  // E_mk:         Denotes the expected value of m_{k},   E(m_{k})
  // E_mkm1:       Denotes the expected value of m_{k-1}, E(m_{k-1})
  // Cov_mk_mk:     Denotes the covariance of m_{k} and m_{k},     Cov(m_{k}, m_{k})
  // Cov_mk_mkm1:   Denotes the covariance of m_{k} and m_{k-1},   Cov(m_{k}, m_{k-1})
  // Cov_mkm1_mkm1: Denotes the covariance of m_{k-1} and m_{k-1}, Cov(m_{k-1}, m_{k-1})

  CorrelatedRockSamples correlated_rock_samples;
  std::vector<std::vector<std::vector<double> > > m_ik = correlated_rock_samples.CreateSamples(i_max, time_line, dist_rock);
  // The dimension of m_ik[k][i] is expected to be equal to 3, in other words we do not expect to receive samples splitted into dynamic and static parts.
  // Now do the separation, which expands m_ik[k][i] to double size:
  SplitSamplesStaticDynamic(m_ik);

  int dim = static_cast<int>(m_ik[0][0].size());
  int K = number_of_timesteps_;

  // Data structures for evolution matrix and correction term
  NRLib::Matrix A_k(dim, dim);
  NRLib::Matrix delta_k(dim, dim);
  NRLib::Vector delta_mu_k(dim);

  NRLib::Vector E_mk(dim);
  NRLib::Vector E_mkm1(dim);
  NRLib::Matrix Cov_mk_mk(dim,dim);
  NRLib::Matrix Cov_mk_mkm1(dim,dim);
  NRLib::Matrix Cov_mkm1_mkm1(dim,dim);

  // For conversion of m_ik to NRLib::Vector
  std::vector<NRLib::Vector> m_k(dim);
  std::vector<NRLib::Vector> m_km1(dim);

  // For all time steps, find A_k, delta_k, delta_mu_k
  for (int k = 1; k < K; ++k)
  {
    // Data spooling: Get m_{k} and m_{k-1} from samples matrix
    NRLib::Vector temp_m_k(i_max);
    NRLib::Vector temp_m_km1(i_max);
    for (int d = 0; d < dim; d++)
    {
      for (int i = 0; i < i_max; ++i) {
        temp_m_k(i)   = m_ik[k][i][d];
        temp_m_km1(i) = m_ik[k-1][i][d];
      }
      m_k[d]   = temp_m_k;
      m_km1[d] = temp_m_km1;
    }

    // Estimation of mean and covariance elements:
    for (int d1 = 0; d1<dim; d1++)
    {
      E_mk(d1)   = NRLib::Mean(m_k[d1]);
      E_mkm1(d1) = NRLib::Mean(m_km1[d1]);

      for (int d2=0; d2<dim; d2++)
      {
        Cov_mk_mk(d1,d2)     = NRLib::Cov(m_k[d1], m_k[d2]);
        Cov_mk_mkm1(d1,d2)   = NRLib::Cov(m_k[d1], m_km1[d2]);
        Cov_mkm1_mkm1(d1,d2) = NRLib::Cov(m_km1[d1], m_km1[d2]);
      }
    }

    NRLib::Matrix SigmaInv;
    double adjustment_factor = 0.001;         // To be changed if DoRobustInversion() produce warnings.
    DoRobustInversion(SigmaInv, Cov_mk_mk, Cov_mk_mkm1, Cov_mkm1_mkm1, dim, adjustment_factor);
    A_k = Cov_mk_mkm1*SigmaInv;               // A_k = Cov(m_{k}, m_{k-1})\Sigma_{k-1}^{-1}
    AdjustMatrixAForm(A_k, dim);              // Should be done before the matrix is used further.

    NRLib::Matrix A_kT = flens::transpose(A_k);
    delta_k = Cov_mk_mk - Cov_mk_mkm1*A_kT;   // \Sigma_{k} - Cov(m_{k}m m_{k-1})A_{k}^T
    delta_mu_k = E_mk - A_k*E_mkm1;           // \mu_k - A_k\mu_{k-1}
    AdjustMatrixDeltaForm(delta_k, dim);

    // Push back of computed matrices and vectors into return parameters.
    evolution_matrix.push_back(A_k);
    cov_correction_term.push_back(delta_k);
    mean_correction_term.push_back(delta_mu_k);
  }
}

void TimeEvolution::DoRobustInversion(NRLib::Matrix & SigmaInv,
                                      NRLib::Matrix & Cov_mk_mk,
                                      NRLib::Matrix & Cov_mk_mkm1,
                                      NRLib::Matrix & Cov_mkm1_mkm1,
                                      int             dim,
                                      double          adjustment_factor)
{
  // The inversion may need to be stabilized.
  // The adjustment made for this is repeated also for the other covariance
  // matrices to ensure block form of evolution matrix and correction term covariance.

  int max_counter     = 10;
  int counter         = 0;
  bool all_ok         = false;
  NRLib::Matrix Sigma = Cov_mkm1_mkm1;
  std::string message("");

  while (!all_ok && counter < max_counter) {
    counter += 1;
    flens::DenseVector<flens::Array<int> > p(Sigma.numRows());
    try {
      NRLib::TriangleFactorize(Sigma, p);            // Upper triangle
      if (DiagonalOK(Sigma, dim)){
        NRLib::TriangleInvert(Sigma, p);
        SigmaInv = Sigma;
        all_ok = true;
      }
      else {
        FixMatrices(Cov_mk_mk, Cov_mk_mkm1, Cov_mkm1_mkm1, dim, adjustment_factor);
        Sigma = Cov_mkm1_mkm1;
      }
    }
    catch (NRLib::Exception & e) {
      FixMatrices(Cov_mk_mk, Cov_mk_mkm1, Cov_mkm1_mkm1, dim, adjustment_factor);
      Sigma = Cov_mkm1_mkm1;
      message += std::string(e.what())+"\n";
    }
  }

  if (!all_ok) {
    std::string text("");
    text += "\nWARNING: The inverse covariance matrix needed for time development could not be found.\n";
    text += "Tried " + NRLib::ToString(counter) + " stabilizations, each with an adjustment factor of " + NRLib::ToString(adjustment_factor) + ".\n";
    text += "Accumulated catch messages from NRLib are: " + message + "\n";
    LogKit::LogFormatted(LogKit::High,text);
    text = "";
    text += "Consider using another adjustment factor for matrix inversion in TimeEvolution::DoRobustInversion(). \n";
    TaskList::addTask(text);
  }

  /*NRLib::Matrix IdentL = SigmaInv*Cov_mkm1_mkm1;
  NRLib::Matrix IdentR = Cov_mkm1_mkm1*SigmaInv;
  PrintToScreen(Cov_mkm1_mkm1, dim, dim, "Cov_mkm1_mkm1: ");
  PrintToScreen(SigmaInv, dim, dim, "SigmaInv: ");
  PrintToScreen(IdentL, dim, dim, "IdentL: ");
  PrintToScreen(IdentR, dim, dim, "IdentR: ");*/
}

void TimeEvolution::FixMatrices(NRLib::Matrix & Cov_mk_mk,
                                NRLib::Matrix & Cov_mk_mkm1,
                                NRLib::Matrix & Cov_mkm1_mkm1,
                                int             dim,
                                double          adjustment_factor)
{
  assert (dim % 2 == 0);
  int dim2 = dim/2;
  NRLib::Matrix tmp_ss = Cov_mk_mk;  // We will only use the static-static block of this covariance matrix.
  for (int d = 0; d < dim2; ++d) {
    double x                          = tmp_ss(d,d) * adjustment_factor;
    if (x == 0.0)
      x = adjustment_factor;
    Cov_mkm1_mkm1(d,d)               += x;
    Cov_mk_mkm1(d,d)                 += x;
    Cov_mk_mk(d,d)                   += x;
    Cov_mkm1_mkm1(d + dim2,d + dim2) += x;
    Cov_mk_mkm1(d + dim2,d + dim2)   += x;
    Cov_mk_mk(d + dim2,d + dim2)     += x;
  }
}

bool TimeEvolution::DiagonalOK(const NRLib::Matrix & A, int dim)
{
  bool d_ok = true;
  double low_limit = FLT_EPSILON; // The value for FLT_EPSILON is defined in float.h and is 1.192092896e-07F in Visual Studio 8, Windows 7.
  for (int d = 0; d < dim; ++d) {
    if (std::abs(A(d, d)) < low_limit)
      d_ok = false;
  }
  return d_ok;
}

void TimeEvolution::AdjustMatrixAForm(NRLib::Matrix & A_k, int dim)
{
  assert(dim % 2 == 0);
  int dim2 = dim/2;
  bool all_ok = true;
  double sensitivity_limit = FLT_EPSILON;

  // 1st block is supposed to be identity matrix:
  for (int d1 = 0; d1 < dim2; d1++) {
    for (int d2 = 0; d2 < dim2; d2++) {
      if (d2 != d1) {
        if (A_k(d1, d2) - 0.0 > sensitivity_limit || A_k(d1, d2) - 0.0 < -sensitivity_limit) // Detecting if form is far off.
          all_ok = false;
        A_k(d1, d2) = 0.0;
      }
      else {
        if (A_k(d1, d2) - 1.0 > sensitivity_limit || A_k(d1, d2) - 1.0 < -sensitivity_limit) // Detecting if form is far off.
          all_ok = false;
        A_k(d1, d2) = 1.0;
      }
    }
  }
  // 2nd block is supposed to be zero matrix:
  for (int d1 = 0; d1 < dim2; d1++) {
    for (int d2 = dim2; d2 < dim; d2++) {
      if (A_k(d1, d2) - 0.0 > sensitivity_limit || A_k(d1, d2) - 0.0 < -sensitivity_limit) // Detecting if form is far off.
          all_ok = false;
      A_k(d1, d2) = 0.0;
    }
  }

  if (!all_ok) {
    std::string text("");
    text += "\nWARNING: Evolution matrix off by too much before adjustment.\n";
    LogKit::LogFormatted(LogKit::High,text);
    text = "";
    text += "Check whether sensitivity is set too restrictive or there is something seriously wrong with matrix in TimeEvolution::AdjustMatrixAForm(). \n";
    TaskList::addTask(text);
  }
}

void TimeEvolution::AdjustMatrixDeltaForm(NRLib::Matrix & delta_k, int dim)
{
  // 1st and 3rd blocks are guaranteed to be zero matrices already, given that AdjustMatrixAForm() has been done.
  assert(dim % 2 == 0);
  int dim2 = dim/2;
  bool all_ok = true;
  double sensitivity_limit = FLT_EPSILON;

  // 2nd block is to be a zero matrix:
  for (int d1 = 0; d1 < dim2; d1++) {
    for (int d2 = dim2; d2 < dim; d2++) {
      if (delta_k(d1, d2) - 0.0 > sensitivity_limit || delta_k(d1, d2) - 0.0 < -sensitivity_limit) // Detecting if form is far off.
          all_ok = false;
      delta_k(d1, d2) = 0.0;
    }
  }

  if (!all_ok) {
    std::string text("");
    text += "\nWARNING: Covariance correction matrix off by too much before adjustment.\n";
    LogKit::LogFormatted(LogKit::High,text);
    text = "";
    text += "Check whether sensitivity is set too restrictive or there is something seriously wrong with matrix in TimeEvolution::AdjustMatrixDeltaForm(). \n";
    TaskList::addTask(text);
  }
}




void TimeEvolution::SplitSamplesStaticDynamic(std::vector<std::vector<std::vector<double> > > & m_ik)
{
  // Splits according to the convention that for a given set of time correlated seismic parameters.
  // the first time instance is considered static and all the rest deviates from this static part by a dynamic component.
  // That is: m_ik[0][i] is the static component, m_ik[k][i] - m_ik[0][i] is the dynamic component for time instance k.

  // Order of indices: m_ik[k][i][d]
  int k_max = int(m_ik.size());
  int i_max = int(m_ik[0].size());
  int dim   = int(m_ik[0][0].size());

  std::vector<double> m_static(dim), m_dynamic(dim);
  for (int i = 0; i < i_max; ++i) {
    m_static = m_ik[0][i];
    //for (int d = 0; d < dim; ++d)                                                         //Random addition just for testing matrix forms!!!
    //  m_static[d] += 100.0 * double(std::rand())/double(RAND_MAX);
    for (int k = 0; k < k_max; ++k) {
      for (int d = 0; d < dim; ++d)
        m_dynamic[d] = m_ik[k][i][d] - m_static[d] /*+ 100.0 * double(std::rand())/double(RAND_MAX)*/; //Random addition just for testing matrix forms!!!
      m_ik[k][i].resize(2*dim);                // Expanding the vector.
      for (int d = 0; d < dim; ++d) {
        m_ik[k][i][d] = m_static[d];
        m_ik[k][i][d + dim] = m_dynamic[d];
      }
    }
  }
  return;
}

void TimeEvolution::PrintToScreen(NRLib::Matrix a, int dim1, int dim2, std::string name)
{
  std::cout << name << std::endl;
  for (int d1 = 0; d1 < dim1; ++d1) {
    std::ostringstream strs;
    for (int d2 = 0; d2 < dim2; ++d2) {
      strs << a(d1,d2);
      strs << " ";
    }
    std::cout << strs.str() << std::endl;
    std::cout << std::endl;
  }
}

void TimeEvolution::PrintToScreen(NRLib::Vector v, int dim, std::string name)
{
  std::cout << name << std::endl;
  std::ostringstream strs;
  for (int d = 0; d < dim; ++d) {
    strs << v(d);
    strs << " ";
  }
  std::cout << strs.str() << std::endl;
  std::cout << std::endl;
}
