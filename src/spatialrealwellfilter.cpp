/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/spatialwellfilter.h"
#include "src/spatialrealwellfilter.h"
#include "avoinversion.h"


SpatialRealWellFilter::SpatialRealWellFilter()
{
}

SpatialRealWellFilter::SpatialRealWellFilter(int    nwells,
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

SpatialRealWellFilter::~SpatialRealWellFilter()
{
  int i,j;
  for(i=0;i<nWells_;i++)
  {
    for(j=0;j<n_[i];j++)
      delete [] priorSpatialCorr_[i][j];
    delete [] priorSpatialCorr_[i];
  }
  delete [] priorSpatialCorr_;

  delete [] n_;
}


void SpatialRealWellFilter::SetPriorSpatialCovariance(const BlockedLogsCommon   * blocked_log,
                                                      int                         wellnr,
                                                      const FFTGrid             * cov_vp,
                                                      const FFTGrid             * cov_vs,
                                                      const FFTGrid             * cov_rho,
                                                      const FFTGrid             * cov_vpvs,
                                                      const FFTGrid             * cov_vprho,
                                                      const FFTGrid             * cov_vsrho)
{
  int n_blocks = blocked_log->GetNumberOfBlocks();

  prior_cov_vp_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vs_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_rho_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vpvs_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vprho_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vsrho_[wellnr].Resize(n_blocks, n_blocks, 0);

  int i1,j1,k1, i2,j2,k2;
  const std::vector<int> ipos = blocked_log->GetIposVector();
  const std::vector<int> jpos = blocked_log->GetJposVector();
  const std::vector<int> kpos = blocked_log->GetKposVector();
  for(int l1 = 0; l1 < n_blocks; l1++)
  {
    i1 = ipos[l1];
    j1 = jpos[l1];
    k1 = kpos[l1];
    for(int l2 = 0; l2 <= l1; l2++)
    {
      i2 = ipos[l2];
      j2 = jpos[l2];
      k2 = kpos[l2];
      prior_cov_vp_[wellnr](l1,l2) = cov_vp->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vp_[wellnr](l2,l1) = prior_cov_vp_[wellnr](l1,l2);
      prior_cov_vs_[wellnr](l1,l2) = cov_vs->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vs_[wellnr](l2,l1) = prior_cov_vs_[wellnr](l1,l2);
      prior_cov_rho_[wellnr](l1,l2) = cov_rho->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_rho_[wellnr](l2,l1) = prior_cov_rho_[wellnr](l1,l2);
      
      prior_cov_vpvs_[wellnr](l1,l2) = cov_vpvs->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vpvs_[wellnr](l2,l1) = prior_cov_vpvs_[wellnr](l1,l2);
      prior_cov_vprho_[wellnr](l1,l2) = cov_vprho->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vprho_[wellnr](l2,l1) = prior_cov_vprho_[wellnr](l1,l2);
      prior_cov_vsrho_[wellnr](l1,l2) = cov_vsrho->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vsrho_[wellnr](l2,l1) = prior_cov_vsrho_[wellnr](l1,l2);
    }
  }
}

void
SpatialRealWellFilter::fillValuesInSigmapost(double    ** sigmapost,
                                             const int *  ipos,
                                             const int *  jpos,
                                             const int *  kpos,
                                             FFTGrid   *  covgrid,
                                             int          n,
                                             int          ni,
                                             int          nj)
{
  covgrid->setAccessMode(FFTGrid::RANDOMACCESS);
  for (int l1=0 ; l1<n ; l1++) {
    int i1 = ipos[l1];
    int j1 = jpos[l1];
    int k1 = kpos[l1];
    for (int l2=0 ; l2<n ; l2++) {
      int i2 = ipos[l2];
      int j2 = jpos[l2];
      int k2 = kpos[l2];
      sigmapost[l1+ni][l2+nj] = covgrid->getRealValueCyclic(i1-i2, j1-j2, k1-k2);
    }
  }
  covgrid->endAccess();
}

void SpatialRealWellFilter::setPriorSpatialCorr(FFTGrid             * parSpatialCorr,
                                                BlockedLogsCommon   * blocked_log,
                                                int                   wellnr)
{
  float constant = parSpatialCorr->getRealValue(0, 0, 0);

  int n = blocked_log->GetNumberOfBlocks();
  priorSpatialCorr_[wellnr] = new double *[n];
  n_[wellnr] = n;
  for(int i=0;i<n;i++)
    priorSpatialCorr_[wellnr][i] = new double[n];

  int i1,j1,k1,l1, i2,j2,k2,l2;
  const std::vector<int> & ipos = blocked_log->GetIposVector();
  const std::vector<int> & jpos = blocked_log->GetJposVector();
  const std::vector<int> & kpos = blocked_log->GetKposVector();
  for(l1=0;l1<n;l1++)
  {
    i1 = ipos[l1];
    j1 = jpos[l1];
    k1 = kpos[l1];
    for(l2=0;l2<=l1;l2++)
    {
      i2 = ipos[l2];
      j2 = jpos[l2];
      k2 = kpos[l2];
      priorSpatialCorr_[wellnr][l1][l2] = parSpatialCorr->getRealValueCyclic(i1-i2,j1-j2,k1-k2) / constant;
      priorSpatialCorr_[wellnr][l2][l1] = priorSpatialCorr_[wellnr][l1][l2];
    }
  }
}

//-------------------------------------------------------------------------------
void SpatialRealWellFilter::doFiltering(std::map<std::string, BlockedLogsCommon *> blocked_logs,
                                        bool                                       useVpRhoFilter,
                                        int                                        nAngles,
                                        const AVOInversion                       * avoInversionResult,
                                        const std::vector<Grid2D *>              & noiseScale,
                                        SeismicParametersHolder                  & seismicParameters)
{
  LogKit::WriteHeader("Creating spatial multi-parameter filter");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  std::vector<NRLib::Matrix> sigmaeVpRho;

  double ** sigmapost;
  double ** sigmapri;

  int lastn = 0;
  int n = 0;
  int nDim = 1;
  for(int i=0;i<nAngles;i++)
    nDim *= 2;

  if(sigmae_.size() == 0) {
    sigmae_.resize(nDim);

    NRLib::Matrix sigmae(3,3);
    NRLib::InitializeMatrix(sigmae, 0.0);

    for(int k=0 ; k < nDim ; k++) {
      sigmae_[k] = sigmae;
    }
  }

  NRLib::Matrix priorCov0 = avoInversionResult->getPriorVar0();

  bool no_wells_filtered = true;
  int w1 = 0;

  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    n = blocked_log->GetNumberOfBlocks();

    if (blocked_log->GetUseForFiltering() == true)
    {
      LogKit::LogFormatted(LogKit::Low,"\nFiltering well "+blocked_log->GetWellName());
      no_wells_filtered = false;

      sigmapost = new double * [3*n];
      for(int i=0;i<3*n;i++)
        sigmapost[i] = new double[3*n];

      sigmapri = new double * [3*n];
      for(int i=0;i<3*n;i++)
        sigmapri[i] = new double[3*n];

      const std::vector<int> & ipos = blocked_log->GetIposVector();
      const std::vector<int> & jpos = blocked_log->GetJposVector();
      const std::vector<int> & kpos = blocked_log->GetKposVector();

      float regularization = Definitions::SpatialFilterRegularisationValue();

      fillValuesInSigmapost(sigmapost, &ipos[0], &jpos[0], &kpos[0], seismicParameters.GetCovVp()     , n, 0  , 0   );
      fillValuesInSigmapost(sigmapost, &ipos[0], &jpos[0], &kpos[0], seismicParameters.GetCovVs()     , n, n  , n   );
      fillValuesInSigmapost(sigmapost, &ipos[0], &jpos[0], &kpos[0], seismicParameters.GetCovRho()    , n, 2*n, 2*n );
      fillValuesInSigmapost(sigmapost, &ipos[0], &jpos[0], &kpos[0], seismicParameters.GetCrCovVpVs() , n, 0  , n   );
      fillValuesInSigmapost(sigmapost, &ipos[0], &jpos[0], &kpos[0], seismicParameters.GetCrCovVpRho(), n, 0  , 2*n );
      fillValuesInSigmapost(sigmapost, &ipos[0], &jpos[0], &kpos[0], seismicParameters.GetCrCovVsRho(), n, 2*n, n   );

      for(int l1=0 ; l1 < n ; l1++) {
        for(int l2=0 ; l2 < n ; l2++) {
          sigmapost[l2 + n  ][l1      ] = sigmapost[l1      ][l2 + n  ];
          sigmapost[l2 + 2*n][l1      ] = sigmapost[l1      ][l2 + 2*n];
          sigmapost[l2 + n  ][l1 + 2*n] = sigmapost[l1 + 2*n][l2 + n  ];
          sigmapri [l1      ][l2      ] = prior_cov_vp_[w1](l1,l2);//priorCov0(0,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri [l1 + n  ][l2 + n  ] = prior_cov_vs_[w1](l1,l2);//priorCov0(1,1)*priorSpatialCorr_[w1][l1][l2];
          sigmapri [l1 + 2*n][l2 + 2*n] = prior_cov_rho_[w1](l1,l2);//priorCov0(2,2)*priorSpatialCorr_[w1][l1][l2];
          if(l1==l2)
          {
            sigmapost[l1      ][l2      ] += regularization*sigmapost[l1      ][l2      ]/sigmapri[l1      ][l2      ];
            sigmapost[l1 + n  ][l2 + n  ] += regularization*sigmapost[l1 + n  ][l2 + n  ]/sigmapri[l1 + n  ][l2 + n  ];
            sigmapost[l1 + 2*n][l2 + 2*n] += regularization*sigmapost[l1 + 2*n][l2 + 2*n]/sigmapri[l1 + 2*n][l2 + 2*n];
            sigmapri [l1      ][l2      ] += regularization;
            sigmapri [l1 + n  ][l2 + n  ] += regularization;
            sigmapri [l1 + 2*n][l2 + 2*n] += regularization;
          }
          sigmapri[l1 + n  ][l2      ] = prior_cov_vpvs_[w1](l1,l2);//priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l2      ][l1 + n  ] = sigmapri[l1 + n  ][l2      ];//priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l1 + 2*n][l2      ] = prior_cov_vprho_[w1](l1,l2);//priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l2      ][l1 + 2*n] = sigmapri[l1 + 2*n][l2      ];//priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l1 + n  ][l2 + 2*n] = prior_cov_vsrho_[w1](l1,l2);//priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l2 + 2*n][l1 + n  ] = sigmapri[l1 + n  ][l2 + 2*n];//priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
        }
      }

      NRLib::SymmetricMatrix Sprior(3*n);
      NRLib::SymmetricMatrix Spost(3*n);

      for(int i = 0 ; i < 3*n ; i++)
        for(int j = 0 ; j <= i ; j++)
          Sprior(j,i) = sigmapri[j][i];

      for(int i = 0 ; i < 3*n ; i++)
        for(int j = 0 ; j <= i ; j++)
          Spost(j,i) = sigmapost[j][i];

      if(useVpRhoFilter == true) //Only additional
        doVpRhoFiltering(sigmaeVpRho,
                         sigmapri,
                         sigmapost,
                         n,
                         blocked_log); //Must do before Cholesky of sigmapri.


      NRLib::WriteMatrixToFile("Sprior_2014.dat", Sprior);
      NRLib::WriteMatrixToFile("Spost_2014.dat", Spost);

      //
      // Filter = I - Sigma_post * inv(Sigma_prior)
      //
      NRLib::Matrix I = NRLib::IdentityMatrix(3*n);
      NRLib::CholeskySolve(Sprior, I);

      NRLib::Matrix Aw = Spost * I;
      //NRLib::WriteMatrixToFile("SpostxSpriinv_2014.dat", Aw);
      Aw = Aw * (-1);
      for(int i=0 ; i<3*n ; i++) {
        Aw(i,i) += 1.0;
      }

      if(useVpRhoFilter == false) { //Save time, since below is not needed then.
        updateSigmaE(sigmae_[0],
                     Aw,
                     Spost,
                     n);
      }

      //NRLib::WriteMatrixToFile("Aw_2014.dat", Aw);

      //NRLib::Vector eigenvalues(3*n);
      //NRLib::Matrix eigenvectors(3*n, 3*n);
      //NRLib::ComputeEigenVectors(Aw, eigenvalues, eigenvectors);

      calculateFilteredLogs(Aw,
                            blocked_log,
                            n,
                            true);

      lastn += n;

      for(int i=0;i<3*n;i++)
      {
        delete [] sigmapost[i];
        delete [] sigmapri[i];
      }
      delete [] sigmapri;
      delete [] sigmapost;
    }
    w1++;
  }

  if(no_wells_filtered == false)
    completeSigmaE(sigmae_,
                   lastn,
                   avoInversionResult,
                   noiseScale);

  if(useVpRhoFilter == true)
    completeSigmaEVpRho(sigmaeVpRho,
                        lastn,
                        avoInversionResult,
                        noiseScale);

  if (no_wells_filtered) {
    LogKit::LogFormatted(LogKit::Low,"\nNo wells have been filtered.\n");
  }

  Timings::setTimeFiltering(wall,cpu);
}
