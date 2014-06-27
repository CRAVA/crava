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


SpatialWellFilter::SpatialWellFilter(int nwells)
{
  nWells_ = nwells;
  priorSpatialCorr_ = new double **[nWells_];
  n_                = new int[nWells_];
  for(int i = 0 ; i < nWells_ ; i++)
    priorSpatialCorr_[i] = NULL;
}

SpatialWellFilter::~SpatialWellFilter()
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

void SpatialWellFilter::setPriorSpatialCorrSyntWell(FFTGrid             * parSpatialCorr,
                                                    SyntWellData        * well,
                                                    int                   wellnr)
{
  parSpatialCorr->setAccessMode(FFTGrid::RANDOMACCESS);

  double minValue = std::pow(10.0,-9);
  int nz = parSpatialCorr->getNz();
  double factorNorm = 1/sqrt(2*NRLib::Pi);
  double endValue = parSpatialCorr->getRealValueCyclic(0,0,nz-2);
  double div = 1/(nz*.1);
  // tapering limit at 90%
  double smoothLimit = nz*.9;

  int n = well->getWellLength();
  priorSpatialCorr_[wellnr] = new double *[n];
  n_[wellnr] = n;
  for(int i=0;i<n;i++)
    priorSpatialCorr_[wellnr][i] = new double[n];

  int i1,j1,k1,l1, i2,j2,k2,l2;
  const int *ipos = well->getIpos();
  const int *jpos = well->getJpos();
  const int *kpos = well->getKpos();

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
      if(abs(k2-k1)>smoothLimit){
        priorSpatialCorr_[wellnr][l1][l2] = std::max(endValue*factorNorm*exp(-(div*(abs(k2-k1)-smoothLimit))*(div*(abs(k2-k1)-smoothLimit))*0.5), minValue);
        priorSpatialCorr_[wellnr][l2][l1] = priorSpatialCorr_[wellnr][l1][l2];
      }
      else{
        priorSpatialCorr_[wellnr][l1][l2] = parSpatialCorr->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        priorSpatialCorr_[wellnr][l2][l1] = priorSpatialCorr_[wellnr][l1][l2];
      }
    }
  }
  parSpatialCorr->endAccess();
}

void SpatialWellFilter::setPriorSpatialCorr(FFTGrid *parSpatialCorr, BlockedLogsCommon * blocked_log, int wellnr)
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

void SpatialWellFilter::doFilteringSyntWells(std::vector<SyntWellData *>              & syntWellData,
                                             const std::vector<std::vector<double> >  & v,
                                             SeismicParametersHolder                  & seismicParameters,
                                             int                                        nWells,
                                             const NRLib::Matrix                      & priorVar0)
{
  (void) v;

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  double ** sigmapost;
  double ** sigmapri;
  double ** imat;
  double ** Aw;

  int lastn = 0;

  // nDim is always 1 for synthetic wells
  int nDim = 1;

  //number of dimensions of transformed elastic variables (2 or 3)
  //int nDimElasticVar = static_cast<int>(V.size());

  if(sigmaeSynt_.size() == 0) {
    sigmaeSynt_.resize(nDim);
    double ** Se = new double *[3];
    for(int i=0; i<3; i++) {
      Se[i] = new double[3] ;
      for(int j=0; j<3; j++)
        Se[i][j] = 0;
    }
    sigmaeSynt_[0] = Se;
  }

  NRLib::Matrix priorCov0 = priorVar0;

  bool no_wells_filtered = true;

  for(int w1=0;w1<nWells;w1++){
    //LogKit::LogFormatted(LogKit::Low,"\nFiltering synthetic well number " + NRLib::ToString(w1+1,1) + "...");
    no_wells_filtered = false;

    int n = syntWellData[w1]->getWellLength();

    sigmapost = new double * [3*n];
    for(int i=0;i<3*n;i++) {
      sigmapost[i] = new double[3*n];
      for(int j=0; j<3*n; j++)
        sigmapost[i][j] = RMISSING;
    }
    sigmapri = new double * [3*n];
    for(int i=0;i<3*n;i++)
      sigmapri[i] = new double[3*n];
    imat = new double * [3*n];
    for(int i=0;i<3*n;i++)
      imat[i] = new double[3*n];
    Aw = new double * [3*n];
    for(int i=0;i<3*n;i++)
      Aw[i] = new double[3*n];

    for(int i=0;i<3*n;i++)
      for(int j=0;j<3*n;j++)
        if(i==j)
          imat[i][j] = 1.0;
        else
          imat[i][j] = 0.0;
    int i1,j1,k1,l1, i2,j2,k2,l2;
    const int *ipos = syntWellData[w1]->getIpos();
    const int *jpos = syntWellData[w1]->getJpos();
    const int *kpos = syntWellData[w1]->getKpos();
    float regularization = Definitions::SpatialFilterRegularisationValue();

    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCovVp(),      n, 0,   0);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCovVs(),      n, n,   n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCovRho(),     n, 2*n, 2*n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCrCovVpVs(),  n, 0,   n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCrCovVpRho(), n, 0,   2*n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCrCovVsRho(), n, 2*n, n);

    // In case the synthetic well is longer than the vertical size of covgrid,
    // set correlation for the relevant grid points to 0
    for(l1=0;l1<3*n;l1++){
      for(l2=0;l2<3*n;l2++){
        if(sigmapost[l1][l2] == RMISSING)
          sigmapost[l1][l2] = 0.0;
      }
    }

    for(l1=0;l1<n;l1++){
      i1 = ipos[l1];
      j1 = jpos[l1];
      k1 = kpos[l1];
      for(l2=0;l2<n;l2++){
        i2 = ipos[l2];
        j2 = jpos[l2];
        k2 = kpos[l2];

        sigmapost[l2 + n  ][l1      ] = sigmapost[l1][n+l2];
        sigmapost[l2 + 2*n][l1      ] = sigmapost[l1][2*n+l2];
        sigmapost[l2 + n  ][l1 + 2*n] = sigmapost[2*n+l1][n+l2];
        sigmapri [l1      ][l2      ] = priorCov0(0,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + n  ][l2 + n  ] = priorCov0(1,1)*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + 2*n][l2 + 2*n] = priorCov0(2,2)*priorSpatialCorr_[w1][l1][l2];
        if(l1==l2){
          sigmapost[l1      ][l2      ] += regularization*sigmapost[l1][l2]/sigmapri[l1][l2];
          sigmapost[l1 + n  ][l2 + n  ] += regularization*sigmapost[n+l1][n+l2]/sigmapri[n+l1][n+l2];
          sigmapost[l1 + 2*n][l2 + 2*n] += regularization*sigmapost[2*n+l1][2*n+l2]/sigmapri[2*n+l1][2*n+l2];
          sigmapri [l1      ][l2      ] += regularization;
          sigmapri [l1 + n  ][l2 + n  ] += regularization;
          sigmapri [l1 + 2*n][l2 + 2*n] += regularization;
        }
        sigmapri[l1 + n  ][l2      ] = priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2      ][l1 + n  ] = priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l1 + 2*n][l2      ] = priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2      ][l1 + 2*n] = priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l1 + n  ][l2 + 2*n] = priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2 + 2*n][l1 + n  ] = priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
      }
    }


    //LogKit::LogFormatted(LogKit::Low,"\n  Cholesky decomposition ...");
    lib_matrCholR(3*n, sigmapri);
    //LogKit::LogFormatted(LogKit::Low,"\n  Equation solving ...");
    lib_matrAXeqBMatR(3*n, sigmapri, imat, 3*n);
    //LogKit::LogFormatted(LogKit::Low,"\n  Matrix multiplication ...\n");
    lib_matr_prod(sigmapost,imat,3*n,3*n,3*n,Aw);

    for(int i=0;i<3*n;i++) {
      for(int j=0;j<3*n;j++){
        Aw[i][j] *=-1.0;
        if(i==j)
          Aw[i][j]+=1.0;
      }
    }

    updateSigmaeSynt(Aw, sigmapost, n);

    lastn += n;

    for(int i=0;i<3*n;i++){
      delete [] Aw[i];
      delete [] sigmapri[i];
      delete [] sigmapost[i];
      delete [] imat[i];
    }
    delete [] Aw;
    delete [] sigmapri;
    delete [] sigmapost;
    delete [] imat;
  }

  if(no_wells_filtered == false){
    // finds the scale at default inversion (all minimum noise in case of local noise)
    NRLib::Matrix Se(3,3);// = sigmaeSynt_[0]; Marit
    Se(0,0) /= lastn;
    Se(1,0) /= lastn;
    Se(1,1) /= lastn;
    Se(2,0) /= lastn;
    Se(2,1) /= lastn;
    Se(2,2) /= lastn;
    Se(0,1)  = Se(1,0);
    Se(0,2)  = Se(2,0);
    Se(1,2)  = Se(2,1);
    adjustDiagSigma(Se);
  }

  if (no_wells_filtered) {
    LogKit::LogFormatted(LogKit::Low,"\nNo synthetic wells have been filtered.\n");
  }

  Timings::setTimeFiltering(wall,cpu);
}

void SpatialWellFilter::doFiltering(std::map<std::string, BlockedLogsCommon *> blocked_logs,
                                    bool                                       useVpRhoFilter,
                                    int                                        nAngles,
                                    const AVOInversion                       * avoInversionResult,
                                    const std::vector<Grid2D *>              & noiseScale,
                                    SeismicParametersHolder                  & seismicParameters)
//-------------------------------------------------------------------------------
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
          sigmapri [l1      ][l2      ] = priorCov0(0,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri [l1 + n  ][l2 + n  ] = priorCov0(1,1)*priorSpatialCorr_[w1][l1][l2];
          sigmapri [l1 + 2*n][l2 + 2*n] = priorCov0(2,2)*priorSpatialCorr_[w1][l1][l2];
          if(l1==l2)
          {
            sigmapost[l1      ][l2      ] += regularization*sigmapost[l1      ][l2      ]/sigmapri[l1      ][l2      ];
            sigmapost[l1 + n  ][l2 + n  ] += regularization*sigmapost[l1 + n  ][l2 + n  ]/sigmapri[l1 + n  ][l2 + n  ];
            sigmapost[l1 + 2*n][l2 + 2*n] += regularization*sigmapost[l1 + 2*n][l2 + 2*n]/sigmapri[l1 + 2*n][l2 + 2*n];
            sigmapri [l1      ][l2      ] += regularization;
            sigmapri [l1 + n  ][l2 + n  ] += regularization;
            sigmapri [l1 + 2*n][l2 + 2*n] += regularization;
          }
          sigmapri[l1 + n  ][l2      ] = priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l2      ][l1 + n  ] = priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l1 + 2*n][l2      ] = priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l2      ][l1 + 2*n] = priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l1 + n  ][l2 + 2*n] = priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
          sigmapri[l2 + 2*n][l1 + n  ] = priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
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


      //
      // Filter = I - Sigma_post * inv(Sigma_prior)
      //
      NRLib::Matrix Aw;
      NRLib::Matrix I = NRLib::IdentityMatrix(3*n);
      NRLib::CholeskySolve(Sprior, I);

      Aw = Spost * I;
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


void
SpatialWellFilter::fillValuesInSigmapostSyntWell(double **sigmapost, const int *ipos, const int *jpos, const int *kpos, FFTGrid *covgrid, int n, int ni, int nj)
{
  double minValue = std::pow(10.0,-9);
  int nz = covgrid->getNz();
  double factorNorm = 1/sqrt(2*NRLib::Pi);
  double endValue = covgrid->getRealValueCyclic(0,0,nz-1);
  double div = 1/(nz*.1);
  // tapering limit at 90%
  double smoothLimit = nz*.9;

  covgrid->setAccessMode(FFTGrid::RANDOMACCESS);
  int i1, j1, k1, l1, i2, j2, k2, l2;
  for(l1=0;l1<n;l1++)
  {
    i1 = ipos[l1];
    j1 = jpos[l1];
    k1 = kpos[l1];
    for(l2=0;l2<n;l2++)
    {
      i2 = ipos[l2];
      j2 = jpos[l2];
      k2 = kpos[l2];
      if(abs(k2-k1) > smoothLimit){
        sigmapost[l1+ni][l2+nj] = std::max(endValue*factorNorm*exp(-(div*(abs(k2-k1)-smoothLimit))*(div*(abs(k2-k1)-smoothLimit))*0.5), minValue);
      }
      else{
        sigmapost[l1+ni][l2+nj] = covgrid->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      }
    }
  }
  covgrid->endAccess();


}

void
SpatialWellFilter::fillValuesInSigmapost(double    ** sigmapost,
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

void
SpatialWellFilter::updateSigmaeSynt(double ** filter, double ** postCov,  int n)
{
  double **sigmaeW;
  sigmaeW = new double * [3*n];
  for(int i=0;i<3*n;i++)
  {
    sigmaeW[i] = new double[3*n];
  }

  lib_matr_prod(filter,postCov,3*n,3*n,3*n,sigmaeW);

  for(int i=0; i<n; i++) {
    NRLib::Matrix Se(3,3); //= sigmaeSynt_[0]; Marit
    Se(0,0) += sigmaeW[i      ][i     ];
    Se(1,0) += sigmaeW[i +   n][i     ];
    Se(2,0) += sigmaeW[i + 2*n][i     ];
    Se(1,1) += sigmaeW[i +   n][i +  n];
    Se(2,1) += sigmaeW[i + 2*n][i +  n];
    Se(2,2) += sigmaeW[i + 2*n][i + 2*n];
  }
  // sigmaeSynt_ Is normalized (1/n) in completeSigmaE, Here well by well is added.

  for(int i=0;i<3*n;i++)
    delete [] sigmaeW[i];

  delete [] sigmaeW;
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
    blockedlogs->SetSpatialFilteredLogs(vpFiltered,  n, "ALPHA_SEISMIC_RESOLUTION", bg_vp);
    blockedlogs->SetSpatialFilteredLogs(vsFiltered,  n, "BETA_SEISMIC_RESOLUTION",  bg_vs);
    blockedlogs->SetSpatialFilteredLogs(rhoFiltered, n, "RHO_SEISMIC_RESOLUTION",   bg_rho);
  }
  else {
    blockedlogs->SetSpatialFilteredLogs(vpFiltered,  n, "ALPHA_FOR_FACIES", bg_vp);
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
