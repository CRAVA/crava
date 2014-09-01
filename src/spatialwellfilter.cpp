/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "lib/timekit.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "lib/lib_matr.h"

#include "src/spatialwellfilter.h"
#include "src/definitions.h"
#include "src/timings.h"
#include "src/modelsettings.h"
#include "src/avoinversion.h"
#include "src/seismicparametersholder.h"
#include "src/blockedlogscommon.h"


SpatialWellFilter::SpatialWellFilter()
{
}

SpatialWellFilter::SpatialWellFilter(int    nwells,
                                     bool   cov_estimated)
{
  nWells_ = nwells;

  // If the prior cov has been estimated, we do not use a common corr_t
  if (cov_estimated){
    prior_cov_vp_.resize(nWells_);
    prior_cov_vs_.resize(nWells_);
    prior_cov_rho_.resize(nWells_);
    prior_cov_vpvs_.resize(nWells_);
    prior_cov_vprho_.resize(nWells_);
    prior_cov_vsrho_.resize(nWells_);
    priorSpatialCorr_ = NULL;
  }
  else{
    priorSpatialCorr_ = new double **[nWells_];
    for(int i = 0 ; i < nWells_ ; i++)
      priorSpatialCorr_[i] = NULL;
  }

  n_                  = new int[nWells_];

}

SpatialWellFilter::~SpatialWellFilter()
{

}



//------------------------------------------------------------------
void SpatialWellFilter::updateSigmaE(NRLib::Matrix       & sigmae,
                                     const NRLib::Matrix & Filter,
                                     const NRLib::Matrix & PostCov,
                                     int                   n)
//------------------------------------------------------------------
{
  NRLib::Matrix sigmaeW = Filter * PostCov;

  for(int i=0 ; i < n ; i++)
  {
    sigmae(0,0) += sigmaeW(i      , i      );
    sigmae(1,0) += sigmaeW(i +   n, i      );
    sigmae(2,0) += sigmaeW(i + 2*n, i      );
    sigmae(1,1) += sigmaeW(i +   n, i +   n);
    sigmae(2,1) += sigmaeW(i + 2*n, i +   n);
    sigmae(2,2) += sigmaeW(i + 2*n, i + 2*n);
  }
  // sigmae Is normalized (1/n) in completeSigmaE, Here well by well is added.
}

//-------------------------------------------------------------------------------
void SpatialWellFilter::completeSigmaE(std::vector<NRLib::Matrix>  & sigmae,
                                       int                           lastn,
                                       const AVOInversion          * avoInversionResult,
                                       const std::vector<Grid2D *> & noiseScale)
//-------------------------------------------------------------------------------
{
  // finds the scale at  default inversion (all minimum noise in case of local noise)
  sigmae[0](0,0) /= lastn;
  sigmae[0](1,0) /= lastn;
  sigmae[0](1,1) /= lastn;
  sigmae[0](2,0) /= lastn;
  sigmae[0](2,1) /= lastn;
  sigmae[0](2,2) /= lastn;
  sigmae[0](0,1)  = sigmae[0](1,0);
  sigmae[0](0,2)  = sigmae[0](2,0);
  sigmae[0](1,2)  = sigmae[0](2,1);

  adjustDiagSigma(sigmae[0]);

  if (sigmae.size() > 1) { // We have local noise
    int nAng = static_cast<int>(noiseScale.size());

    NRLib::Vector maxScale(nAng);
    for (int angle=0 ; angle < nAng ; angle++) {
      double maxS = noiseScale[angle]->FindMax(RMISSING);
      double minS = noiseScale[angle]->FindMin(RMISSING);
      maxScale(angle) = maxS/minS;
    }

    NRLib::Vector scale(nAng);
    for (int angle=0 ; angle<nAng ; angle++)
      scale(angle) = 1.0;

    NRLib::Matrix G(nAng, 3);
    avoInversionResult->computeG(G);

    NRLib::Matrix dummy(3,3);
    NRLib::Matrix H(3,3);

    avoInversionResult->newPosteriorCovPointwise(dummy,
                                                 G,
                                                 scale,
                                                 H);

    NRLib::Matrix postCovAdj = H * H;

    NRLib::SymmetricMatrix symPostCovAdj(3);
    for(int i = 0 ; i < 3 ; i++)
      for(int j = 0 ; j <= i ; j++)
        symPostCovAdj(j,i) = postCovAdj(j,i);

    NRLib::SymmetricMatrix SigmaPri0  = avoInversionResult->getSymmetricPriorVar0();
    NRLib::Matrix          filter     = avoInversionResult->computeFilter(SigmaPri0, symPostCovAdj);

    NRLib::Matrix          sigmaE0    = filter * postCovAdj;

    // Idealy sigmaE0 = sigmae[0] but this is not the case due to
    // spatial effects.  We therfor adjust the sigmae[0]
    // by   sqrt(sigmaETmp*sigmaE0^-1)*sigmae[0]* sqrt(sigmaE0^-1*sigmaETmp)

    // Interpret conf counter bitwise - 0 means min value noise for that component, 1 means max.

    for(unsigned int conf = 1; conf < sigmae.size() ; conf++) {

      // Compute pointwise filter
      int factor = 1;
      for (int angle=0 ; angle<nAng ; angle++) {
        if ((conf & factor) > 0)
          scale(angle) = maxScale(angle);
        else
          scale(angle) = 1.0;
        factor *= 2;
      }

      avoInversionResult->newPosteriorCovPointwise(dummy,
                                                   G,
                                                   scale,
                                                   H);    // H is the square root of postCovAdj

      postCovAdj = H * H;

      for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j <= i ; j++)
          symPostCovAdj(j,i) = postCovAdj(j,i);

      filter = avoInversionResult->computeFilter(SigmaPri0, symPostCovAdj);

      NRLib::Matrix sigmaETmp = filter * postCovAdj;

      computeSigmaEAdjusted(sigmae[0],
                            sigmaE0,
                            sigmaETmp,
                            sigmae[conf]);

      adjustDiagSigma(sigmae[conf]);
    }
  }
}



void
SpatialWellFilter::computeSigmaEAdjusted(NRLib::Matrix & sigmae,
                                         NRLib::Matrix & sigmaE0,
                                         NRLib::Matrix & sigmaETmp,
                                         NRLib::Matrix & sigmaEAdj)
{
  // sigmaEAdj = sqrt(sigmaETmp*sigmaE0^-1) * sigmae * sqrt(sigmaE0^-1*sigmaETmp)

  int n = sigmae.numRows();

  NRLib::Vector Eval0(n);
  NRLib::Matrix Evec0(n,n);

  NRLib::ComputeEigenVectors(sigmaE0, Eval0, Evec0);

  NRLib::Vector EvalTmp(n);
  NRLib::Matrix EvecTmp(n,n);

  NRLib::ComputeEigenVectors(sigmaETmp, EvalTmp, EvecTmp);

  double max0   = NRLib::FindLargestElement(Eval0);
  double maxTmp = NRLib::FindLargestElement(EvalTmp);
  for (int i=0 ; i<n ; i++) {
    Eval0(i)   = std::max(Eval0(i),max0/1000);
    EvalTmp(i) = std::max(EvalTmp(i),maxTmp/1000);
  }

  // Computes: = sqrt(sigmaETmp*sigmaE0^-1) = UTmp*LambdaTmp^(1/2)*UTmpT*U0*Lambda0^(-1/2)*U0T

  NRLib::Matrix T1(n,n);
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++) {
      T1(i,j) = EvecTmp(i,j) * sqrt(EvalTmp(j));   // T1 = UTmp*LambdaTmp^(1/2)
    }
  }

  NRLib::Matrix T2 = NRLib::transpose(EvecTmp);
  NRLib::Matrix T3 = T1 * T2;                      // T3 = UTmp*LambdaTmp^(1/2)*UTmpT

  T1 = T3 * Evec0;                                 // T1 = UTmp*LambdaTmp^(1/2)*UTmpT*U0

  for(int i=0 ; i<n ; i++) {
    for(int j=0 ; j<n ; j++) {
      T1(i,j) = T1(i,j) * sqrt(1.0/Eval0(j));      // T1 = UTmp*LambdaTmp^(1/2)*UTmpT*U0*Lambda0^(-1/2)
    }
  }

  T2 = NRLib::transpose(Evec0);
  T3 = T1 * T2;                                    // T3 = UTmp*LambdaTmp^(1/2)*UTmpT*U0*Lambda0^(-1/2)*U0T
  T1 = T3 * sigmae;                                // T1 = sqrt(sigmaETmp*sigmaE0^-1)*sigmae

  T2 = NRLib::transpose(T3);

  sigmaEAdj = T1 * T2;                             // sigmaEAdj = sqrt(sigmaETmp*sigmaE0^-1)*sigmae*sqrt(sigmaE0^-1*sigmaETmp)
}

void SpatialWellFilter::doVpRhoFiltering(std::vector<NRLib::Matrix> &  sigmaeVpRho,
                                         double                     ** sigmapri,
                                         double                     ** sigmapost,
                                         const int                     n,
                                         BlockedLogsCommon          *  blockedLogs)
//---------------------------------------------------------------------------------
{
  int m = 2*n;

  NRLib::Matrix tmp1(m,m);
  NRLib::Matrix tmp2(m,m);

  for (int j=0 ; j<n ; j++) {
    for (int i=0 ; i<n ; i++) {
      tmp1(i,   j  ) = sigmapri [i  ][j  ];
      tmp1(i+n, j  ) = sigmapri [i+m][j  ];
      tmp1(i  , j+n) = sigmapri [i  ][j+m];
      tmp1(i+n, j+n) = sigmapri [i+m][j+m];

      tmp2(i,   j  ) = sigmapost[i  ][j  ];
      tmp2(i+n, j  ) = sigmapost[i+m][j  ];
      tmp2(i,   j+n) = sigmapost[i  ][j+m];
      tmp2(i+n, j+n) = sigmapost[i+m][j+m];
    }
  }

  NRLib::SymmetricMatrix Sprior2(m);
  NRLib::SymmetricMatrix Spost2(m);

  for(int j = 0 ; j < m ; j++)
    for(int i = 0 ; i <= j ; i++)
      Sprior2(i,j) = tmp1(i,j);

  for(int j = 0 ; j < m ; j++)
    for(int i = 0 ; i <= j ; i++)
      Spost2(i,j) = tmp2(i,j);

  NRLib::Matrix Aw;
  NRLib::Matrix I = NRLib::IdentityMatrix(m);
  NRLib::CholeskySolve(Sprior2, I);
  Aw = Spost2 * I;
  Aw = Aw * (-1);
  for(int i=0 ; i<m ; i++) {
    Aw(i,i) += 1.0;
  }

  calculateFilteredLogs(Aw, blockedLogs, n, false);

  updateSigmaEVpRho(sigmaeVpRho,
                    Aw,
                    Spost2,
                    static_cast<int>(sigmae_.size()),
                    n);
}

//---------------------------------------------------------------------------------
void SpatialWellFilter::updateSigmaEVpRho(std::vector<NRLib::Matrix> & sigmaeVpRho,
                                          const NRLib::Matrix        & Aw,
                                          const NRLib::Matrix        & Spost,
                                          int                          nDim,
                                          int                          n)
//---------------------------------------------------------------------------------
{
  if (sigmaeVpRho.size() == 0) { // then first time alocate memory
    sigmaeVpRho.resize(nDim);

    NRLib::Matrix sigmae(2, 2);
    NRLib::InitializeMatrix(sigmae, 0.0);

    for (int k=0 ; k < nDim ; k++) {
      sigmaeVpRho[k] = sigmae;
    }
  }

  NRLib::Matrix sigma = Aw * Spost;

  //
  // NBNB-PAL: Bug? førsteindeksen på sigmaeVpRho[0][0][0] står
  // stille hele tiden. Det er ingen n-avhengighet.
  //
  for(int i=0 ; i < n ; i++) {
    sigmaeVpRho[0](0,0) += sigma(i    , i    );
    sigmaeVpRho[0](1,0) += sigma(i + n, i    );
    sigmaeVpRho[0](1,1) += sigma(i + n, i + n);
  }
}

//------------------------------------------------------------------------------------
void SpatialWellFilter::completeSigmaEVpRho(std::vector<NRLib::Matrix>  & sigmaeVpRho,
                                            int                           lastn,
                                            const AVOInversion          * avoInversionResult,
                                            const std::vector<Grid2D *> & noiseScale)
//------------------------------------------------------------------------------------
{
  sigmaeVpRho[0](0,0) /= lastn;
  sigmaeVpRho[0](1,0) /= lastn;
  sigmaeVpRho[0](1,1) /= lastn;
  sigmaeVpRho[0](0,1)  = sigmae_[0](1,0);

  adjustDiagSigma(sigmaeVpRho[0]);

  if(sigmaeVpRho.size() > 1) { // We have local noise
    // initialization
    int nAng = static_cast<int>(noiseScale.size());
    std::vector<double> maxScale(nAng);
    double maxS, minS;
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      maxScale[angle] = maxS/minS;
    }

    NRLib::Matrix sigmaPri0 = avoInversionResult->getPriorVar0();

    NRLib::SymmetricMatrix priCovVpRho(2);
    priCovVpRho(0,0) = sigmaPri0(0,0);
    priCovVpRho(0,1) = sigmaPri0(0,2);
    priCovVpRho(1,0) = sigmaPri0(2,0);
    priCovVpRho(1,1) = sigmaPri0(2,2);

    NRLib::Matrix dummy(3,3);
    NRLib::Matrix help(nAng,3);

    NRLib::Vector scale(nAng);
    for(int angle=0 ; angle<nAng ; angle++)
      scale(angle) = 1.0;

    NRLib::Matrix G(nAng, 3);

    avoInversionResult->computeG(G);
    avoInversionResult->newPosteriorCovPointwise(dummy,
                                                 G,
                                                 scale,
                                                 help);

    NRLib::Matrix postCovAdj = help * help;

    NRLib::SymmetricMatrix postCovVpRho(2);
    postCovVpRho(0,0) = postCovAdj(0,0);
    postCovVpRho(1,0) = postCovAdj(2,0);
    postCovVpRho(1,1) = postCovAdj(2,2);

    NRLib::Matrix filter = avoInversionResult->computeFilter(priCovVpRho,
                                                             postCovVpRho);

    NRLib::Matrix sigmaE0 = filter * postCovVpRho;

    for(unsigned int conf = 1; conf < sigmaeVpRho.size();conf++) {

      //Compute pointwise filter
      int factor = 1;
      for(int angle=0;angle<nAng;angle++) {
        if((conf & factor) > 0)
          scale(angle) = maxScale[angle];
        else
          scale(angle) = 1.0;
        factor *= 2;
      }

      avoInversionResult->newPosteriorCovPointwise(dummy, G, scale,  help);  // help is the square root of postCovAdj

      postCovAdj = help * help;

      postCovVpRho(0,0) = postCovAdj(0,0);
      postCovVpRho(0,1) = postCovAdj(0,2);
      postCovVpRho(1,0) = postCovAdj(2,0);
      postCovVpRho(1,1) = postCovAdj(2,2);

      filter = avoInversionResult->computeFilter(priCovVpRho,
                                                 postCovVpRho);

      NRLib::Matrix sigmaETmp = filter * postCovVpRho;

      computeSigmaEAdjusted(sigmaeVpRho[0],
                            sigmaE0,
                            sigmaETmp,
                            sigmaeVpRho[conf]);
      adjustDiagSigma(sigmaeVpRho[conf]);
    }
  }
}

void SpatialWellFilter::calculateFilteredLogs(const NRLib::Matrix & Aw,
                                              BlockedLogsCommon   * blockedlogs,
                                              int                   n,
                                              bool                  useVs)
//------------------------------------------------------------------------------
{
  int nLogs = 2;
  if(useVs == true)
    nLogs++;

  NRLib::Vector residuals(nLogs*n);

  int currentEnd = 0;
  const std::vector<double> & vp    = blockedlogs->GetVpBlocked();
  const std::vector<double> & bg_vp = blockedlogs->GetVpHighCutBackground();
  MakeInterpolatedResiduals(vp, bg_vp, n, currentEnd, residuals);
  currentEnd += n;

  const std::vector<double> & vs    = blockedlogs->GetVsBlocked();
  const std::vector<double> & bg_vs = blockedlogs->GetVsHighCutBackground();
  if(useVs == true) {
    MakeInterpolatedResiduals(vs, bg_vs, n, currentEnd, residuals);
    currentEnd += n;
  }
  const std::vector<double> & rho    = blockedlogs->GetRhoBlocked();
  const std::vector<double> & bg_rho = blockedlogs->GetRhoHighCutBackground();
  MakeInterpolatedResiduals(rho, bg_rho, n, currentEnd, residuals);

  NRLib::Vector filteredVal = Aw * residuals;

  std::vector<double> vpFiltered(n);
  std::vector<double> vsFiltered(n);
  std::vector<double> rhoFiltered(n);

  for(int i=0;i<n;i++)
  {
    int offset = 0;
    if(vp[i] == RMISSING)
      vpFiltered[i] = 0.0;
    else
      vpFiltered[i] = static_cast<double>(filteredVal(i+offset));
    offset += n;

    if(useVs == true) {
      if(vs[i] == RMISSING)
        vsFiltered[i] = 0.0;
      else
        vsFiltered[i] = static_cast<double>(filteredVal(i+offset));
      offset += n;
    }

    if(rho[i] == RMISSING)
      rhoFiltered[i] = 0.0;
    else
      rhoFiltered[i] = static_cast<double>(filteredVal(i+offset));
  }

  if(useVs == true) {
    blockedlogs->SetSpatialFilteredLogs(vpFiltered,  n, "VP_SEISMIC_RESOLUTION", bg_vp);
    blockedlogs->SetSpatialFilteredLogs(vsFiltered,  n, "VS_SEISMIC_RESOLUTION",  bg_vs);
    blockedlogs->SetSpatialFilteredLogs(rhoFiltered, n, "RHO_SEISMIC_RESOLUTION",   bg_rho);
  }
  else {
    blockedlogs->SetSpatialFilteredLogs(vpFiltered,  n, "VP_FOR_FACIES", bg_vp);
    blockedlogs->SetSpatialFilteredLogs(rhoFiltered, n, "RHO_FOR_FACIES"  , bg_rho);
  }

}

void SpatialWellFilter::MakeInterpolatedResiduals(const std::vector<double> & bwLog,
                                                  const std::vector<double> & bwLogBG,
                                                  const int                   n,
                                                  const int                   offset,
                                                  NRLib::Vector &             residuals)
//--------------------------------------------------------------------------
{
  //
  // When the log starts with a missing value
  //
  int first_nonmissing = 0;

  if (bwLog[0] == RMISSING)
  {
    int i = 1;
    while (bwLog[i] == RMISSING)
      i++;

    first_nonmissing = i;
    double first_residual = static_cast<double>(bwLog[i] - bwLogBG[i]);

    for (i = 0 ; i < first_nonmissing ; i++)
      residuals(offset + i) = first_residual;
  }

  //
  // The general case (also handles logs ending with missing values)
  //
  int nmiss = 0;
  for(int i=first_nonmissing ; i<n ; i++)
  {
    if(bwLog[i] != RMISSING)
    {
      double res_i = double(bwLog[i] - bwLogBG[i]);
      residuals(offset + i) = res_i;

      if(nmiss>0)
      {
        for(int j=1 ; j<=nmiss ; j++)
        {
          double w = static_cast<double>(j)/static_cast<double>(nmiss + 1);
          residuals(offset + i - j) *= w;
          residuals(offset + i - j) += (1.0 - w)*res_i;
        }
      }
      nmiss = 0;
    }
    else
    {
      nmiss++;
      residuals(offset + i) = residuals(offset + i - 1);
    }
  }
}



// The variances used for smootihng in faciesprob might be very small.
// Therefore eigenvalues are adjusted in order to be able to invert matrix.

void SpatialWellFilter::adjustDiagSigma(NRLib::Matrix & sigmae)
{
  int    n   = sigmae.numRows();
  double eps = 0.0001;

  NRLib::Vector Eval(n);
  NRLib::Matrix Evec(n,n);

  NRLib::ComputeEigenVectors(sigmae, Eval, Evec);

  for (int i=1 ; i < n ; i++) {
    if (Eval(i)/Eval(0) < eps) {
      Eval(i) = eps*Eval(0);
    }
  }

  NRLib::Matrix EvalMat(n,n);
  NRLib::InitializeMatrix(EvalMat, 0.0);

  for (int i=0 ; i < n ; i++) {
    EvalMat(i,i) = Eval(i);
  }

  NRLib::Matrix H;
  NRLib::Matrix EvecT;

  H      = Evec * EvalMat;
  EvecT  = NRLib::transpose(Evec);
  sigmae = H * EvecT;
}
