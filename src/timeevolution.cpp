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


TimeEvolution::TimeEvolution(int                                     i_max,
                             TimeLine                              & time_line,
                             const std::vector<DistributionsRock*> & dist_rock)
{
  LogKit::WriteHeader("Setting up matrices for time evolution");
  std::list<int> time;
  time_line.GetAllTimes(time);
  number_of_timesteps_ = static_cast<int>(dist_rock.size());
  SetUpEvolutionMatrices(evolution_matrix_, cov_correction_term_, mean_correction_term_, i_max, time_line, dist_rock); //NBNB to be replaced
}


std::vector<std::vector<std::vector<double> > >

  TimeEvolution::returnCorrelatedSample(int                                   i_max,
                                      TimeLine                              & time_line,
                                      const std::vector<DistributionsRock*> & dist_rock)
{
  CorrelatedRockSamples correlated_rock_samples;
  std::vector<std::vector<std::vector<double> > > sample= correlated_rock_samples.CreateSamplesExtended(i_max, time_line, dist_rock);
  // The dimension of m_ik[k][i] is expected to be equal to 3, in other words we do not return samples splitted into dynamic and static parts.
  return sample;
}


void TimeEvolution::SetUpEvolutionMatrices(std::vector< NRLib::Matrix>           & evolution_matrix,
                                           std::vector< NRLib::Matrix>           & cov_correction_term,
                                           std::vector< NRLib::Vector>           & mean_correction_term,
                                           int                                     i_max,
                                           TimeLine                              & time_line,
                                           const std::vector<DistributionsRock*> & dist_rock)
{
  // Covariance matrices are made robust before the
  // Equations from NR Note SAND/22/11 are used
  // A note on variable names in this function:

  // A_k:        Denotes the time evolution matrix in the expression m_{k} = A_{k}*m_{k-1} + \delta m_{k}
  // delta_k:    Denotes the covariance of the correction term \delta m_{k}, denoted with symbol \delta_k
  // delta_mu_k: Denotes the mean of the correction term \delta m_{k}, denoted with symbol \delta \mu_{k}

  // E_mk:         Denotes the expected value of m_{k},   E(m_{k})
  // E_mkm1:       Denotes the expected value of m_{k-1}, E(m_{k-1})
  // Cov_mk_mk:     Denotes the covariance of m_{k} and m_{k},     Cov(m_{k}, m_{k})
  // Cov_mk_mkm1:   Denotes the covariance of m_{k} and m_{k-1},   Cov(m_{k}, m_{k-1})
  // Cov_mkm1_mkm1: Denotes the covariance of m_{k-1} and m_{k-1}, Cov(m_{k-1}, m_{k-1})
  double adjustment_factor=1e-6;

  CorrelatedRockSamples correlated_rock_samples;
  std::vector<std::vector<std::vector<double> > > m_ik = correlated_rock_samples.CreateSamples(i_max, time_line, dist_rock);

  //write seismic parameters to check ok
   if(true)
  {
    NRLib::Matrix rockSamples(i_max,3*number_of_timesteps_);

    for(int i = 0;i<i_max;i++)
      for(int j=0;j < number_of_timesteps_;j++ )
        for(int k=0;k<3;k++)
        {
          int ind=k+j*3;
          rockSamples(i,ind)= m_ik[j][i][k];
        }

    NRLib::WriteMatrixToFile("SeisParSampleEvolution.dat", rockSamples);
  }


  // The dimension of m_ik[k][i] is expected to be equal to 3 (dim),
  // Now do the separation, which expands m_ik[k][i] to size dim*number_of_timesteps_

  int dim = static_cast<int>(m_ik[0][0].size());
  int K = number_of_timesteps_;

  std::vector<std::vector<double> > vectorSample = SplitSamplesStaticDynamic(m_ik);

  // computes covariance off all variables
  NRLib::Matrix Cov_all(dim*K,dim*K);
 // NRLib::Matrix Cov_robust1(dim*K,dim*K);
  NRLib::Matrix Cov_robust(dim*K,dim*K);
  NRLib::Vector E_all(dim*K);

  for(int k=0;k < K*dim; k ++)
  {
    E_all(k)       = NRLib::Mean(vectorSample[k]);
    for(int l=0;l<K*dim;l++)
    {
      Cov_all(k,l)     = NRLib::Cov(vectorSample[k], vectorSample[l]);
    }
  }

  NRLib::WriteMatrixToFile("Cov_all_Not_adjust.dat",Cov_all);
  NRLib::WriteVectorToFile("E_all.dat",E_all);
  // makes the covariance robust.
  for(int k=1; k<K; k++)
    for(int d=0;d<dim;d++)
       Cov_all(k*dim+d,k*dim+d) += Cov_all(k*dim+d,k*dim+d)*0.01+Cov_all(d,d)*0.001;// makes dynamic part robust by adding to  diagonal
  for(int d=0;d<dim;d++)
    Cov_all(d,d)*=1.01;    // makes static part robust by adding to  diagonal

  NRLib::WriteMatrixToFile("Cov_all.dat",Cov_all);

  Cov_robust = makeCovRobust( E_all ,Cov_all, vectorSample, dim, K*dim-1 ,adjustment_factor ); // makes dynamic part robust by adjusting for sample
  //Cov_robust = makeCovRobust( E_all ,Cov_robust1, vectorSample, 0, dim-1 ,adjustment_factor ); // makes dynamic part robust by adjusting for sample
  NRLib::WriteMatrixToFile("Cov_robust.dat",Cov_robust);


  //std::cout << "Warning: Adjustment under construction in TimeEvolution::SetUpEvolutionMatrices " << std::endl;

  // Data structures for evolution matrix and correction term
  NRLib::Matrix A_k(2*dim, 2*dim);
  NRLib::Matrix delta_k(2*dim, 2*dim);
  NRLib::Vector delta_mu_k(2*dim);

  NRLib::Vector E_mk(2*dim);
  NRLib::Vector E_mkm1(2*dim);
  NRLib::Matrix Cov_mk_mk(2*dim,2*dim);
  NRLib::Matrix Cov_mk_mkm1(2*dim,2*dim);
  NRLib::Matrix Cov_mkm1_mkm1(2*dim,2*dim);

  // For all time steps, find A_k, delta_k, delta_mu_k
  for (int k = 1; k < K; ++k)
  {

    // Setting mean and covariance elements:
    for (int d1 = 0; d1<dim; d1++)
    {
      E_mkm1(d1)       = E_all(d1);
      if(k>1)
        E_mkm1(d1+dim) = E_all(d1+(k-1)*dim);
      else
        E_mkm1(d1+dim) = 0;

      E_mk(d1)     =  E_all(d1);
      E_mk(d1+dim) =  E_all(d1+k*dim);
      for (int d2=0; d2<dim; d2++)
      {
        Cov_mk_mk(d1,d2)         = Cov_robust(d1,d2);
        Cov_mk_mk(d1+dim,d2)     = Cov_robust(d1+k*dim,d2);
        Cov_mk_mk(d1,d2+dim)     = Cov_robust(d1,d2+k*dim);
        Cov_mk_mk(d1+dim,d2+dim) = Cov_robust(d1+k*dim,d2+k*dim);


        Cov_mk_mkm1(d1,d2)         = Cov_robust(d1,d2);
        Cov_mkm1_mkm1(d1,d2)       = Cov_robust(d1,d2);
        Cov_mk_mkm1(d1+dim,d2)     = Cov_robust(d1+k*dim,d2);

        if(k>1){
          Cov_mk_mkm1(d1,d2+dim)     = Cov_robust(d1,d2+(k-1)*dim);
          Cov_mk_mkm1(d1+dim,d2+dim) = Cov_robust(d1+k*dim,d2+(k-1)*dim);

          Cov_mkm1_mkm1(d1+dim,d2)     = Cov_robust(d1+(k-1)*dim,d2);
          Cov_mkm1_mkm1(d1,d2+dim)     = Cov_robust(d1,d2+(k-1)*dim);
          Cov_mkm1_mkm1(d1+dim,d2+dim) = Cov_robust(d1+(k-1)*dim,d2+(k-1)*dim);
        }else{
          Cov_mk_mkm1(d1,d2+dim)     = 0;
          Cov_mk_mkm1(d1+dim,d2+dim) = 0;

          Cov_mkm1_mkm1(d1+dim,d2)     = 0;
          Cov_mkm1_mkm1(d1,d2+dim)     = 0;
          Cov_mkm1_mkm1(d1+dim,d2+dim) = 0;
        }
      }
    }
    if(k==1){
      initial_cov_  = Cov_mkm1_mkm1;
      initial_mean_ = E_mkm1;
    }
    NRLib::Matrix SigmaInv;

    NRLib::WriteMatrixToFile("Cov_mk_mk.dat",Cov_mk_mk);
    NRLib::WriteMatrixToFile("Cov_mk_mkm1.dat",Cov_mk_mkm1);
    NRLib::WriteMatrixToFile("Cov_mkm1_mkm1.dat",Cov_mkm1_mkm1);
    DoRobustInversion(SigmaInv, Cov_mkm1_mkm1, 2*dim,adjustment_factor);
    A_k = Cov_mk_mkm1*SigmaInv;               // A_k = Cov(m_{k}, m_{k-1})\Sigma_{k-1}^{-1}
    AdjustMatrixAForm(A_k, 2*dim);              // Should be done before the matrix is used further.

    NRLib::Matrix A_kT = flens::transpose(A_k);
    delta_k = Cov_mk_mk - Cov_mk_mkm1*A_kT;   // \Sigma_{k} - A_{k}Cov(m_{k}m m_{k-1})A_{k}^T
    delta_mu_k = E_mk - A_k*E_mkm1;           // \mu_k - A_k\mu_{k-1}
    AdjustMatrixDeltaForm(delta_k, 2*dim);
    NRLib::WriteMatrixToFile("A_k.dat",A_k);
    NRLib::WriteMatrixToFile("delta_k.dat",delta_k);
    NRLib::WriteVectorToFile("E_mk.dat",E_mk);
    NRLib::WriteVectorToFile("E_mkm1.dat",E_mkm1);
    NRLib::WriteVectorToFile("delta_mu_k.dat",delta_mu_k);

    // Push back of computed matrices and vectors into return parame
    evolution_matrix.push_back(A_k);
    cov_correction_term.push_back(delta_k);
    mean_correction_term.push_back(delta_mu_k);
  }
}

void TimeEvolution::SetUpEvolutionMatrices2(std::vector< NRLib::Matrix>           & evolution_matrix,
                                           std::vector< NRLib::Matrix>           & cov_correction_term,
                                           std::vector< NRLib::Vector>           & mean_correction_term,
                                           int                                     i_max,
                                           TimeLine                              & time_line,
                                           const std::vector<DistributionsRock*> & dist_rock)
{
  // Old version

  // Equations from NR Note SAND/22/11
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

  //write seismic parameters to check ok
   if(true)
  {
    NRLib::Matrix rockSamples(i_max,3*number_of_timesteps_);

    for(int i = 0;i<i_max;i++)
      for(int j=0;j < number_of_timesteps_;j++ )
        for(int k=0;k<3;k++)
        {
          int ind=k+j*3;
          rockSamples(i,ind)= m_ik[j][i][k];
        }

    NRLib::WriteMatrixToFile("SeisParSampleEvolution.dat", rockSamples);
  }


  // The dimension of m_ik[k][i] is expected to be equal to 3, in other words we do not expect to receive samples splitted into dynamic and static parts.
  // Now do the separation, which expands m_ik[k][i] to double size:


  SplitSamplesStaticDynamic2(m_ik);

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

  initial_cov_.resize(dim,dim);
  initial_mean_.resize(dim);

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
      if(k==1)
        initial_mean_(d1)=E_mkm1(d1);

      for (int d2=0; d2<dim; d2++)
      {
        Cov_mk_mk(d1,d2)     = NRLib::Cov(m_k[d1], m_k[d2]);
        Cov_mk_mkm1(d1,d2)   = NRLib::Cov(m_k[d1], m_km1[d2]);
        Cov_mkm1_mkm1(d1,d2) = NRLib::Cov(m_km1[d1], m_km1[d2]);

        if(d1>3 && d2>3){  // makes predictions with less bias (we think coupled with exp below)
          Cov_mk_mk(d1,d2)     += E_mk(d1) *E_mk(d2) ;
          Cov_mk_mkm1(d1,d2)   += E_mk(d1) *E_mkm1(d2) ;
          Cov_mkm1_mkm1(d1,d2) += E_mkm1(d1) *E_mkm1(d2) ;
        }
      }
      if(d1>3){  // makes predictions more robust
        E_mk(d1) =0.0;
        E_mkm1(d1)=0.0;
      }

    }
    if(k==1)
    {
      initial_cov_=Cov_mkm1_mkm1;
      initial_mean_=E_mkm1;
    }

    NRLib::Matrix SigmaInv;
    double adjustment_factor = 1e-6;         // To be changed if DoRobustInversion() produce warnings.
   // NRLib::WriteMatrixToFile("Cov_mk_mk.dat",Cov_mk_mk);
    //NRLib::WriteMatrixToFile("Cov_mk_mkm1.dat",Cov_mk_mkm1);
    //NRLib::WriteMatrixToFile("Cov_mkm1_mkm1.dat",Cov_mkm1_mkm1);
    DoRobustInversion(SigmaInv, Cov_mkm1_mkm1, dim,adjustment_factor);
    A_k = Cov_mk_mkm1*SigmaInv;               // A_k = Cov(m_{k}, m_{k-1})\Sigma_{k-1}^{-1}
    AdjustMatrixAForm(A_k, dim);              // Should be done before the matrix is used further.

    NRLib::Matrix A_kT = flens::transpose(A_k);
    delta_k = Cov_mk_mk - Cov_mk_mkm1*A_kT;   // \Sigma_{k} - A_{k}Cov(m_{k}m m_{k-1})A_{k}^T
    delta_mu_k = E_mk - A_k*E_mkm1;           // \mu_k - A_k\mu_{k-1}
    AdjustMatrixDeltaForm(delta_k, dim);
    //NRLib::WriteMatrixToFile("A_k.dat",A_k);
    //NRLib::WriteMatrixToFile("delta_k.dat",delta_k);
    //NRLib::WriteVectorToFile("E_mk.dat",E_mk);
    //NRLib::WriteVectorToFile("E_mkm1.dat",E_mkm1);

    // Push back of computed matrices and vectors into return parameters.
    evolution_matrix.push_back(A_k);
    cov_correction_term.push_back(delta_k);
    mean_correction_term.push_back(delta_mu_k);
  }
}

void TimeEvolution::DoRobustInversion(NRLib::Matrix & SigmaInv,
                                   //   NRLib::Matrix & Cov_mk_mk,
                                   //   NRLib::Matrix & Cov_mk_mkm1,
                                      const NRLib::Matrix  Cov_mkm1_mkm1,
                                      int             dim,
                                      double          adjustment_factor)
{
  // The inversion may need to be stabilized.
  // The adjustment made for this is repeated also for the other covariance
  // matrices to ensure block form of evolution matrix and correction term covariance.
  SigmaInv.resize(dim,dim);
  for(int i=0;i<dim;i++)
      for(int j=0;j<dim;j++)
        SigmaInv(i,j)=0.0;

  NRLib::Vector eVals(dim);
  NRLib::Matrix eVec(dim,dim);
  NRLib::Matrix tmp(dim,dim);
  tmp=Cov_mkm1_mkm1;
  NRLib::ComputeEigenVectors(tmp,eVals,eVec);
  double maxEval = eVals(0);
  for(int k=1;k<dim;k++)
    if(maxEval < eVals(k))
      maxEval = eVals(k);

  for(int k=0;k<dim;k++)
    for(int i=0;i<dim;i++)
      for(int j=0;j<dim;j++)
      {
        double multiplyer = eVals(k) > maxEval*adjustment_factor ? 1.0/eVals(k) : 0.0;
        SigmaInv(i,j)+= multiplyer*eVec(i,k)*eVec(j,k);
      }
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
  // adding on diagonal of all theree matrixes corresponds to an assumption that the additional part
  // has correlation one between timesteps
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

  for (int d1 = 0; d1 < dim; d1++) { //NBNB OK makes pos def
    delta_k(d1, d1) = delta_k(d1, d1) *1.1;
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

std::vector<std::vector<double> >
TimeEvolution::SplitSamplesStaticDynamic(std::vector<std::vector<std::vector<double> > >  m_ik) const
{
  // Alligns Data by removing one dimension in the data,
  // Keeps the static in first three collumns, keeps the dynamic part in the remaining.
  int k_max = int(m_ik.size());// number of surveys
  int i_max = int(m_ik[0].size());// number of samples
  int dim   = int(m_ik[0][0].size());// 3 when vp, vs,rho

  std::vector<std::vector<double> > vectorData(k_max*dim);

  for (int k = 0; k < k_max*dim; ++k){
    vectorData[k].resize(i_max);
  }

  for (int i = 0; i < i_max; ++i) {
    for(int d=0;d<dim;d++)
      vectorData[d][i] = m_ik[0][i][d]; // static part in first three collumns

    for(int k=1;k<k_max;k++)
      for(int d=0;d<dim;d++)
        vectorData[k*3+d][i] =m_ik[k][i][d] -m_ik[0][i][d]; // dynamic part in remaining collumns
  }
  return vectorData;
}

void TimeEvolution::SplitSamplesStaticDynamic2(std::vector<std::vector<std::vector<double> > > & m_ik) const
{
  // Splits according to the convention that for a given set of time correlated seismic parameters.
  // the first time instance is considered static and all the rest deviates from this static part by a dynamic component.
  // That is: m_ik[0][i] is the static component, m_ik[k][i] - m_ik[0][i] is the dynamic component for time instance k.

  // Order of indices: m_ik[k][i][d]
  int k_max = int(m_ik.size());// number of surveys
  int i_max = int(m_ik[0].size());// number of samples
  int dim   = int(m_ik[0][0].size());// 3 when vp, vs,rho

  std::vector<double> m_static(dim), m_dynamic(dim);
  for (int i = 0; i < i_max; ++i) {
    m_static = m_ik[0][i];
    for (int k = 0; k < k_max; ++k) {
      for (int d = 0; d < dim; ++d)
        m_dynamic[d] = m_ik[k][i][d] - m_static[d];
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

NRLib::Vector
  TimeEvolution::computePriorMeanStaticAndDynamicLastTimeStep()
{
 NRLib::Vector priorMean=initial_mean_;
 NRLib::Vector tmp;

 for(int time_step=0;time_step<number_of_timesteps_-1;time_step++)
 {
    NRLib::Matrix evolution_matrix     = getEvolutionMatrix(time_step);
    NRLib::Vector mean_correction_term = getMeanCorrectionTerm(time_step);
    tmp  = evolution_matrix*priorMean;
    priorMean = tmp+ mean_correction_term; // Note mean_correction_term adjust for the fact that the
 }                                         // static part influence the dynamic, hence it is not zero everywhere
 return priorMean;
}

NRLib::Matrix
  TimeEvolution::computePriorCovStaticAndDynamicLastTimeStep()
{
  NRLib::Matrix priorCov=initial_cov_;
  NRLib::Matrix tmp;

  for(int time_step=0;time_step<number_of_timesteps_-1;time_step++)
  {
    NRLib::Matrix evolution_matrix     = getEvolutionMatrix(time_step);
    NRLib::Matrix cov_correction_term  = getCovarianceCorrectionTerm(time_step);
    tmp=evolution_matrix*priorCov;
    priorCov=tmp*transpose(evolution_matrix);
    priorCov+= cov_correction_term;
  }


  return priorCov;
}

NRLib::Matrix
TimeEvolution::makeCovRobust(NRLib::Vector  mu, NRLib::Matrix sigma, std::vector<std::vector<double>  > &sample,
                                  int start_component,int last_component, double  adjustment_factor)
{
  // make covariance matrix robust towards sample.
  // by adjusting the eigen values of the covariance matrix of dynamic values


  double eps =adjustment_factor;

  int nd    = last_component-start_component+1;
  int n     = sample[0].size();
  int nFull = sample.size();

  NRLib::Matrix AdjustedSigmaFull(nFull,nFull);
  for(int i=0;i<nFull;i++)
    for(int j=0;j<nFull;j++)
    {
      AdjustedSigmaFull(i,j)=sigma(i,j);
    }


  NRLib::Matrix AdjustedSigma(nd,nd);
  for(int i=start_component;i<=last_component;i++)
    for(int j=start_component;j<=last_component;j++)
    {
      AdjustedSigma(i-start_component,j-start_component)=sigma(i,j);
    }

  // Computes and robustify  Eigen Vectors and values
  NRLib::Vector eVals(nd);
  NRLib::Matrix eVec(nd,nd);
  NRLib::Matrix tmp(nd,nd);
  NRLib::Matrix dummySigma=AdjustedSigma; // to avoid overwriting AdjustSigma
  NRLib::ComputeEigenVectors(dummySigma,eVals,eVec);

  double maxVal = 0;
  for(int i=0;i<nd;i++)
    maxVal=std::max(maxVal,eVals(i));
  for(int i=0;i<nd;i++)
    eVals(i) =std::max(eVals(i),maxVal*eps);

  NRLib::Vector delta(nd);
  NRLib::Vector isNotOK(n);
  double maxLim = nd+2.0*sqrt(2.0*nd);

  double maxFactor =  maxLim / nd;

  for(int i=0;i<n;i++){

      for(int j=0;j<nd;j++){ // get current deviation from mean
        delta(j)=sample[start_component+j][i]-mu(start_component+j);
      }

      NRLib::Vector factor(nd);
      for(int j=0;j<nd;j++){ //
        factor(j)=0;
        for(int k=0;k<nd;k++){
           factor(j) +=  delta(k)*eVec(k,j);
        }
      }

      double b=0.0;
      for(int j=0;j<nd;j++){
        b += factor(j)*factor(j)/eVals(j);
      }

      if( b > maxLim  ){// if   mahalanobis distance is too large we must update covariance.
       // we identify those factors which are too large
       // scale them to be such that they get little weigth if the deviation is small

       NRLib::Vector proposedUChange(nd);
       for(int j=0;j<nd;j++){
         double contributionj = factor(j)*factor(j)/eVals(j);
         if(contributionj > maxFactor){
           proposedUChange(j)=factor(j)*sqrt((contributionj -maxFactor)/contributionj );
         }else{
          proposedUChange(j)=0.0;
         }
       }

       NRLib::Vector v(nd);
       for(int j=0;j<nd;j++){ // this gives the change in "original"
         v(j)=0.0;
         for(int k=0;k<nd;k++){
           v(j) +=  proposedUChange(k)*eVec(j,k);
         }
       }

       double c=0.0;
       double d=0.0;
       for(int j=0;j<nd;j++){
         c += proposedUChange(j)*proposedUChange(j)/eVals(j);
         d += factor(j)*proposedUChange(j)/eVals(j);
       }

       double a = 1/(d*d/(b- maxLim)-c); // this balances such thatthe norm is maxLim

       // Updates covariance and eigen values/vectors
       for(int j=0;j<nd;j++){
         for(int k=0;k<nd;k++){
           AdjustedSigma(j,k) +=  a*v(j)*v(k);
         }
       }
       dummySigma=AdjustedSigma;
       NRLib::ComputeEigenVectors(dummySigma,eVals,eVec);

       maxVal = 0;
       for(int i=0;i<nd;i++)
         maxVal=std::max(maxVal,eVals(i));
       for(int i=0;i<nd;i++)
         eVals(i) =std::max(eVals(i),maxVal*eps);
      }
  }



  for(int i=start_component;i<=last_component;i++){
    for(int j=start_component;j<=last_component;j++){
      AdjustedSigmaFull(i,j)=AdjustedSigma(i-start_component,j-start_component);
    }
  }

  return AdjustedSigmaFull;
}
