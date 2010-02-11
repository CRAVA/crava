#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "lib/lib_matr.h"
#include "lib/timekit.hpp"

#include "nrlib/iotools/logkit.hpp"

#include "src/spatialwellfilter.h"
#include "src/welldata.h"
#include "src/timings.h"
#include "src/model.h"
#include "src/corr.h"
#include "src/model.h"
#include "src/crava.h"

SpatialWellFilter::SpatialWellFilter(int nwells)
{
  nWells_ = nwells;
  priorSpatialCorr_ = new double **[nWells_];
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

  for(j=0;j<static_cast<int>(sigmae_.size());j++) {
    for(i=0;i<3;i++)
      delete [] sigmae_[j][i];
    delete [] sigmae_[j];
  }

  for(j=0;j<static_cast<int>(sigmaeVpRho_.size());j++) {
    for(i=0;i<2;i++)
      delete [] sigmaeVpRho_[j][i];
    delete [] sigmaeVpRho_[j];
  }

  delete [] n_;
}


void SpatialWellFilter::setPriorSpatialCorr(FFTGrid *parSpatialCorr, WellData *well, int wellnr)
{
  int n = well->getBlockedLogsOrigThick()->getNumberOfBlocks();
  priorSpatialCorr_[wellnr] = new double *[n];
  for(int i=0;i<n;i++)
    priorSpatialCorr_[wellnr][i] = new double[n];
  
  int i1,j1,k1,l1, i2,j2,k2,l2;
  const int *ipos = well->getBlockedLogsOrigThick()->getIpos();
  const int *jpos = well->getBlockedLogsOrigThick()->getJpos();
  const int *kpos = well->getBlockedLogsOrigThick()->getKpos();
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
      priorSpatialCorr_[wellnr][l1][l2] = parSpatialCorr->getRealValueCyclic(i1-i2,j1-j2,k1-k2); 
      priorSpatialCorr_[wellnr][l2][l1] = priorSpatialCorr_[wellnr][l1][l2]; 
    }
  }
}

void SpatialWellFilter::doFiltering(Corr *corr, WellData **wells, int nWells, bool useVpRhoFilter, int nAngles,
                                    const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale)
{
 Utils::writeHeader("Creating spatial multi-parameter filter");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  n_ = new int[nWells];
  double ** sigmapost;
  double ** sigmapri;
  double ** imat;
  double ** Aw;
  nData_ = 0;
  for(int w1=0;w1<nWells;w1++)
    nData_ += wells[w1]->getBlockedLogsOrigThick()->getNumberOfBlocks();

  int lastn = 0;
  int n = 0;
  int nDim = 1;
  for(int i=0;i<nAngles;i++)
    nDim *= 2;
  if(sigmae_.size() == 0) {
    sigmae_.resize(nDim);
    double ** sigmae = new double * [3];
    for(int i=0;i<3;i++) {
      sigmae[i] = new double[3];
      for(int j=0;j<3;j++)
        sigmae[i][j] = 0;
    }
    sigmae_[0] = sigmae;
    
    for(int k=1;k<nDim;k++) {
      sigmae = new double * [3];
      for(int i=0;i<3;i++) {
        sigmae[i] = new double[3];
        for(int j=0;j<3;j++)
          sigmae[i][j] = 0;
      }
      sigmae_[k] = sigmae;
    }
  }

  for(int w1=0;w1<nWells;w1++)
  {   
    LogKit::LogFormatted(LogKit::LOW,"\nFiltering well "+wells[w1]->getWellname());

    n = wells[w1]->getBlockedLogsOrigThick()->getNumberOfBlocks();

    sigmapost = new double * [3*n];
    for(int i=0;i<3*n;i++)
      sigmapost[i] = new double[3*n];
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
    const int *ipos = wells[w1]->getBlockedLogsOrigThick()->getIpos();
    const int *jpos = wells[w1]->getBlockedLogsOrigThick()->getJpos();
    const int *kpos = wells[w1]->getBlockedLogsOrigThick()->getKpos();
    float regularization = Definitions::SpatialFilterRegularisationValue();
    fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCovAlpha(), n, 0, 0);
    fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCovBeta(), n, n, n);
    fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCovRho(), n, 2*n, 2*n);
    fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCrCovAlphaBeta(), n, 0, n);
    fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCrCovAlphaRho(), n, 0, 2*n);
    fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCrCovBetaRho(), n, 2*n, n);

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
        // sigmapost[l1      ][l2      ] = corr->getPostCovAlpha()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        // sigmapost[l1 + n  ][l2 + n  ] = corr->getPostCovBeta()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        // sigmapost[l1 + 2*n][l2 + 2*n] = corr->getPostCovRho()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        // sigmapost[l1      ][l2 + n  ] = corr->getPostCrCovAlphaBeta()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        // sigmapost[l1      ][l2 + 2*n] = corr->getPostCrCovAlphaRho()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);       
        // sigmapost[l1 + 2*n][l2 + n  ] = corr->getPostCrCovBetaRho()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        sigmapost[l2 + n  ][l1      ] = sigmapost[l1][n+l2];
        sigmapost[l2 + 2*n][l1      ] = sigmapost[l1][2*n+l2];
        sigmapost[l2 + n  ][l1 + 2*n] = sigmapost[2*n+l1][n+l2];
        sigmapri [l1      ][l2      ] = corr->getPriorVar0()[0][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + n  ][l2 + n  ] = corr->getPriorVar0()[1][1]*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + 2*n][l2 + 2*n] = corr->getPriorVar0()[2][2]*priorSpatialCorr_[w1][l1][l2];
        if(l1==l2)
        {
          sigmapost[l1      ][l2      ] += regularization*sigmapost[l1][l2]/sigmapri[l1][l2];
          sigmapost[l1 + n  ][l2 + n  ] += regularization*sigmapost[n+l1][n+l2]/sigmapri[n+l1][n+l2];
          sigmapost[l1 + 2*n][l2 + 2*n] += regularization*sigmapost[2*n+l1][2*n+l2]/sigmapri[2*n+l1][2*n+l2];
          sigmapri [l1      ][l2      ] += regularization;
          sigmapri [l1 + n  ][l2 + n  ] += regularization;
          sigmapri [l1 + 2*n][l2 + 2*n] += regularization;
        }
        sigmapri[l1 + n  ][l2      ] = corr->getPriorVar0()[1][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2      ][l1 + n  ] = corr->getPriorVar0()[1][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l1 + 2*n][l2      ] = corr->getPriorVar0()[2][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2      ][l1 + 2*n] = corr->getPriorVar0()[2][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l1 + n  ][l2 + 2*n] = corr->getPriorVar0()[2][1]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2 + 2*n][l1 + n  ] = corr->getPriorVar0()[2][1]*priorSpatialCorr_[w1][l1][l2];
      }
    }

    if(useVpRhoFilter == true) //Only additional
      doVpRhoFiltering(const_cast<const double **>(sigmapri), 
                       const_cast<const double **>(sigmapost), 
                       n, 
                       wells[w1]->getBlockedLogsOrigThick(),
                       cravaResult,
                       noiseScale); //Must do before Cholesky of sigmapri.

    double ** sigmapriKeep = NULL;
    if(nDim > 1) {
      sigmapriKeep = new double * [3*n];
      for(int i=0;i<3*n;i++)
      {
        sigmapriKeep[i] = new double[3*n];
        for(int j=0;j<3*n;j++)
          sigmapriKeep[i][j] = sigmapri[i][j];
      }
    }

    LogKit::LogFormatted(LogKit::LOW,"\n  Performing a Cholesky decomposition ...");
    lib_matrCholR(3*n, sigmapri);
    LogKit::LogFormatted(LogKit::LOW,"\n  Solving an equation ...");
    lib_matrAXeqBMatR(3*n, sigmapri, imat, 3*n);
    LogKit::LogFormatted(LogKit::LOW,"\n  Make a matrix multiplication ...\n");
    lib_matr_prod(sigmapost,imat,3*n,3*n,3*n,Aw);

    for(int i=0;i<3*n;i++)
      for(int j=0;j<3*n;j++)
        {
        Aw[i][j] *=-1.0;
        if(i==j)
          Aw[i][j]+=1.0;
      }

    if(useVpRhoFilter == false) //Save time, since below is not needed then.
      updateSigmaE(Aw, sigmapriKeep, sigmapost, n, cravaResult, noiseScale);

    calculateFilteredLogs(Aw, wells[w1]->getBlockedLogsOrigThick(), n, true);

    n_[w1] = n;
    lastn += n;

    for(int i=0;i<3*n;i++)
    {
      delete [] Aw[i];
      delete [] sigmapri[i];
      if(sigmapriKeep != NULL)
        delete [] sigmapriKeep[i];
      delete [] sigmapost[i];
      delete [] imat[i];
    }
    delete [] Aw;
    delete [] sigmapriKeep;
    delete [] sigmapri;
    delete [] sigmapost;
    delete [] imat;  
  }

  completeSigmaE(lastn);

  if(useVpRhoFilter == true) 
    completeSigmaEVpRho(lastn);

  Timings::setTimeFiltering(wall,cpu);
}

void
SpatialWellFilter::fillValuesInSigmapost(double **sigmapost, const int *ipos, const int *jpos, const int *kpos, FFTGrid *covgrid, int n, int ni, int nj)
{
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
        sigmapost[l1+ni][l2+nj] = covgrid->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      }
  }
  covgrid->endAccess();


}

void
SpatialWellFilter::updateSigmaE(double ** filter, double ** priCov, double ** postCov, 
                                int n, const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale)
{
  double **sigmaeW;
  sigmaeW = new double * [3*n];
  for(int i=0;i<3*n;i++)
  {
    sigmaeW[i] = new double[3*n];
  }

  lib_matr_prod(filter,postCov,3*n,3*n,3*n,sigmaeW);

  for(int i=0;i<n;i++)
  {
    sigmae_[0][0][0] += sigmaeW[i      ][i     ];
    sigmae_[0][1][0] += sigmaeW[i +   n][i     ];
    sigmae_[0][2][0] += sigmaeW[i + 2*n][i     ];
    sigmae_[0][1][1] += sigmaeW[i +   n][i +  n];
    sigmae_[0][2][1] += sigmaeW[i + 2*n][i +  n];
    sigmae_[0][2][2] += sigmaeW[i + 2*n][i + 2*n];
  }

  if(sigmae_.size() > 1) {
    int nAng = noiseScale.size();
    std::vector<double> maxScale(nAng);
    double maxS, minS;
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      maxScale[angle] = maxS/minS;
    }

    double **G = new double*[nAng];
    for(int i=0;i<nAng;i++)
      G[i] = new double[3];

    cravaResult->computeG(G);

    double ** H    = new double * [3];
    double ** junk = new double * [3];
    for(int i=0;i<3;i++) {
      H[i]    = new double[3];
      junk[i] = new double[3];
    }
    double ** fullH = new double * [3*n];
    double ** HA    = new double * [3*n];
    double ** Sdiff = new double * [3*n];
    double ** SmdU  = new double * [3*n];

    for(int i=0;i<3*n;i++) {
      fullH[i] = new double[3*n];
      HA[i]    = new double[3*n];
      Sdiff[i] = new double[3*n];
      SmdU[i]  = new double[3*n];
    }

    //Interpret conf counter bitwise - 0 means min value noise for that component, 1 means max.
    for(unsigned int conf = 1; conf < sigmae_.size();conf++) { 
      //Compute H matrix
      std::vector<double> scale(nAng);
      int factor = 1;
      for(int angle=0;angle<nAng;angle++) {
        if((conf & factor) > 0)
          scale[angle] = maxScale[angle];
        else
          scale[angle] = 1.0;
        factor *= 2;
      }

      cravaResult->newPosteriorCovPointwise(H, G, scale, junk);

      for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++)
          if(i==j) {
            for(int k=0;k<3;k++)
              for(int l=0;l<3;l++)
                fullH[i+n*k][j+n*l] = H[k][l];
          }
          else {
            for(int k=0;k<3;k++)
              for(int l=0;l<3;l++)
                fullH[i+n*k][j+n*l] = 0.0;
          }
      }

      //H is shrinkage. Updated filter is HA,
      //updated posterior covariance is Sm+H(Sm|d-Sm).
      //Then Se = ASm|d becomes
      //Se = HA(Sm + H(Sm|d-Sm))
      lib_matrCopy(postCov, 3*n, 3*n, Sdiff);
      lib_matrSubtMat(priCov, 3*n, 3*n, Sdiff);
      lib_matr_prod(fullH, Sdiff, 3*n, 3*n, 3*n, SmdU);
      lib_matrAddMat(priCov, 3*n, 3*n, SmdU); //SmdU is now updated posterior cov.
      lib_matr_prod(fullH, filter, 3*n, 3*n, 3*n, HA);
      lib_matr_prod(HA, SmdU, 3*n, 3*n, 3*n, sigmaeW);

      for(int i=0;i<n;i++)
      {
        sigmae_[conf][0][0] += sigmaeW[i      ][i     ];
        sigmae_[conf][1][0] += sigmaeW[i +   n][i     ];
        sigmae_[conf][2][0] += sigmaeW[i + 2*n][i     ];
        sigmae_[conf][1][1] += sigmaeW[i +   n][i +  n];
        sigmae_[conf][2][1] += sigmaeW[i + 2*n][i +  n];
        sigmae_[conf][2][2] += sigmaeW[i + 2*n][i + 2*n];
      }
    }

    for(int i=0;i<3;i++) {
      delete [] junk[i];
      delete [] H[i];
    }
    delete [] junk;
    delete [] H;

    for(int i=0;i<nAng;i++)
      delete [] G[i];
    delete [] G;

    for(int i=0;i<3*n;i++) {
      delete [] fullH[i];
      delete [] HA[i];
      delete [] Sdiff[i];
      delete [] SmdU[i];
    }
    delete [] fullH;
    delete [] HA;
    delete [] Sdiff;
    delete [] SmdU;
  }

  for(int i=0;i<3*n;i++) 
    delete [] sigmaeW[i];

  delete [] sigmaeW;
}

void
SpatialWellFilter::completeSigmaE(int lastn)
{
  for(unsigned int i=0;i<sigmae_.size();i++) { 
    sigmae_[i][0][0] /= lastn;
    sigmae_[i][1][0] /= lastn;
    sigmae_[i][1][1] /= lastn;
    sigmae_[i][2][0] /= lastn;
    sigmae_[i][2][1] /= lastn;
    sigmae_[i][2][2] /= lastn;
    sigmae_[i][0][1]  = sigmae_[i][1][0];
    sigmae_[i][0][2]  = sigmae_[i][2][0];
    sigmae_[i][1][2]  = sigmae_[i][2][1];
   
    adjustDiagSigma(sigmae_[i], 3);
  }
}


void
SpatialWellFilter::doVpRhoFiltering(const double ** sigmapri, const double ** sigmapost, int n,
                                    BlockedLogs * blockedLogs, const Crava * cravaResult,
                                    const std::vector<Grid2D *> & noiseScale)
{
  double ** sigmapri2  = new double *[2*n];
  double ** sigmapost2 = new double *[2*n];
  double ** imat       = new double *[2*n];
  double ** sigma      = new double *[2*n];
  double ** Aw         = new double *[2*n];

  for(int i=0;i<2*n;i++) {
    sigmapri2[i]  = new double[2*n];
    sigmapost2[i] = new double[2*n];
    imat[i]       = new double[2*n];
    sigma[i]      = new double[2*n];
    Aw[i]         = new double[2*n];
  }

  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++) {
      sigmapri2[i][j]      = sigmapri[i][j];
      sigmapost2[i][j]     = sigmapost[i][j];

      sigmapri2[i+n][j]    = sigmapri[i+2*n][j];
      sigmapost2[i+n][j]   = sigmapost[i+2*n][j];

      sigmapri2[i][j+n]    = sigmapri[i][j+2*n];
      sigmapost2[i][j+n]   = sigmapost[i][j+2*n];

      sigmapri2[i+n][j+n]  = sigmapri[i+2*n][j+2*n];
      sigmapost2[i+n][j+n] = sigmapost[i+2*n][j+2*n];

      imat[i][j]     = 0.0;
      imat[i+n][j]   = 0.0;
      imat[i][j+n]   = 0.0;
      imat[i+n][j+n] = 0.0;
    }
    imat[i][i]     = 1.0;
    imat[i+n][i+n] = 1.0;
  }

  double ** sigmapri2Keep = NULL;
  if(sigmae_.size() > 1) {//Interpolate for local noise
    sigmapri2Keep = new double * [2*n];
    for(int j=0;j<2*n;j++) {
      sigmapri2Keep[j] = new double[2*n];
      for(int i=0;i<2*n;i++)
        sigmapri2Keep[j][i] = sigmapri2[j][i];
    }
  }
  lib_matrCholR(2*n, sigmapri2);
  lib_matrAXeqBMatR(2*n, sigmapri2, imat, 2*n);
  lib_matr_prod(sigmapost2,imat,2*n,2*n,2*n,Aw);
  for(int i=0;i<2*n;i++)
    for(int j=0;j<2*n;j++)
    {
      Aw[i][j] *=-1.0;
      if(i==j)
        Aw[i][j]+=1.0;
    }

  calculateFilteredLogs(Aw, blockedLogs, n, false);

  updateSigmaEVpRho(Aw, sigmapri2Keep, sigmapost2, sigmae_.size(), n, cravaResult, noiseScale);
  
  for(int i=0;i<2*n;i++)
  {
    delete [] Aw[i];
    delete [] sigmapri2[i];
    if(sigmapri2Keep != NULL)
      delete [] sigmapri2Keep[i];
    delete [] sigmapost2[i];
    delete [] imat[i];
    delete [] sigma[i];
  }
  delete [] Aw;
  delete [] sigmapri2;
  delete [] sigmapri2Keep;
  delete [] sigmapost2;
  delete [] imat;  
  delete [] sigma;
}

void
SpatialWellFilter::updateSigmaEVpRho(double ** filter, 
                                     double ** priCov, 
                                     double ** postCov, 
                                     int nDim, 
                                     int n, 
                                     const Crava * cravaResult,
                                     const std::vector<Grid2D *> & noiseScale)
{
  if(sigmaeVpRho_.size() == 0) {
    sigmaeVpRho_.resize(nDim);
    double ** sigmae = new double * [2];
    for(int i=0;i<2;i++) {
      sigmae[i] = new double[2];
      for(int j=0;j<2;j++)
        sigmae[i][j] = 0;
    }
    sigmaeVpRho_[0] = sigmae;
    
    for(int k=1;k<nDim;k++) {
      sigmae = new double * [2];
      for(int i=0;i<2;i++) {
        sigmae[i] = new double[2];
        for(int j=0;j<2;j++)
          sigmae[i][j] = 0;
      }
      sigmaeVpRho_[k] = sigmae;
    }
  }

  double **sigma;
  sigma = new double * [2*n];
  for(int i=0;i<2*n;i++)
  {
    sigma[i] = new double[2*n];
  }

  lib_matr_prod(filter,postCov,2*n,2*n,2*n,sigma);

  for(int i=0;i<n;i++)
  {
    sigmaeVpRho_[0][0][0] += sigma[i      ][i     ];
    sigmaeVpRho_[0][1][0] += sigma[i +   n][i     ];
    sigmaeVpRho_[0][1][1] += sigma[i +   n][i +  n];
  }

  if(sigmaeVpRho_.size() > 1) {
    int nAng = noiseScale.size();
    std::vector<double> maxScale(nAng);
    double maxS, minS;
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      maxScale[angle] = maxS/minS;
    }

    double **G = new double*[nAng];
    for(int i=0;i<nAng;i++)
      G[i] = new double[3];

    cravaResult->computeG(G);

    double ** H    = new double * [3];
    double ** junk = new double * [3];
    for(int i=0;i<nAng;i++) {
      H[i]    = new double[3];
      junk[i] = new double[3];
    }
    double ** fullH = new double * [2*n];
    double ** HA    = new double * [2*n];
    double ** Sdiff = new double * [2*n];
    double ** SmdU  = new double * [2*n];

    for(int i=0;i<2*n;i++) {
      fullH[i] = new double[2*n];
      HA[i]    = new double[2*n];
      Sdiff[i] = new double[2*n];
      SmdU[i]  = new double[2*n];
    }

    //Interpret conf counter bitwise - 0 means min value noise for that component, 1 means max.
    for(unsigned int conf = 1; conf < sigmae_.size();conf++) { 
      //Compute H matrix
      std::vector<double> scale(nAng);
      int factor = 1;
      for(int angle=0;angle<nAng;angle++) {
        if((conf & factor) > 0)
          scale[angle] = maxScale[angle];
        else
          scale[angle] = 1.0;
        factor *= 2;
      }

      cravaResult->newPosteriorCovPointwise(H, G, scale, junk);

      for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++)
          for(int k=0;k<2;k++) {
            int k3 = 2*k;
            for(int l=0;l<2;l++) {
              int l3 = 2*l;
              fullH[i+n*k][j+n*l] = H[k3][l3];
            }
          }
      }

      //H is shrinkage. Updated filter is HA,
      //updated posterior covariance is Sm+H(Sm|d-Sm).
      //Then Se = ASm|d becomes
      //Se = HA(Sm + H(Sm|d-Sm))
      lib_matrCopy(postCov, 2*n, 2*n, Sdiff);
      lib_matrSubtMat(priCov, 2*n, 2*n, Sdiff);
      lib_matr_prod(fullH, Sdiff, 2*n, 2*n, 2*n, SmdU);
      lib_matrAddMat(priCov, 2*n, 2*n, SmdU); //SmdU is now updated posterior cov.
      lib_matr_prod(fullH, filter, 2*n, 2*n, 2*n, HA);
      lib_matr_prod(HA, SmdU, 2*n, 2*n, 2*n, sigma);

      for(int i=0;i<n;i++)
      {
        sigmaeVpRho_[0][0][0] += sigma[i      ][i     ];
        sigmaeVpRho_[0][1][0] += sigma[i +   n][i     ];
        sigmaeVpRho_[0][1][1] += sigma[i +   n][i +  n];
      }
    }

    for(int i=0;i<nAng;i++)
      delete [] G[i];
    delete [] G;

    for(int i=0;i<3;i++) {
      delete junk[i];
      delete H[i];
    }
    delete [] junk;
    delete [] H;

    for(int i=0;i<2*n;i++) {
      delete [] fullH[i];
      delete [] HA[i];
      delete [] Sdiff[i];
      delete [] SmdU[i];
    }
    delete [] fullH;
    delete [] HA;
    delete [] Sdiff;
    delete [] SmdU;
  }

  for(int i=0;i<2*n;i++)
    delete [] sigma[i];
  delete [] sigma;
}

void
SpatialWellFilter::completeSigmaEVpRho(int lastn)
{
  for(unsigned int i=0;i<sigmaeVpRho_.size();i++) {
    sigmaeVpRho_[i][0][0] /= lastn;
    sigmaeVpRho_[i][1][0] /= lastn;
    sigmaeVpRho_[i][1][1] /= lastn;
    sigmaeVpRho_[i][0][1]  = sigmae_[i][1][0];

    adjustDiagSigma(sigmaeVpRho_[i], 2);
  }
}

void SpatialWellFilter::calculateFilteredLogs(double **Aw, BlockedLogs *blockedlogs, int n, bool useVs)
{
  int nLogs = 2;
  if(useVs == true)
    nLogs++;

  double **filterval;
  filterval = new double *[nLogs*n];
  int i;
  for(i=0;i<nLogs*n;i++)
    filterval[i] = new double[1];
  
  double **residuals;
  residuals = new double * [nLogs*n];
  for(i=0;i<nLogs*n;i++)
    residuals[i] = new double[1];

  int currentEnd = 0;

  const float * alpha    = blockedlogs->getAlpha();
  const float * bgAlpha  = blockedlogs->getAlphaHighCutBackground();
  MakeInterpolatedResiduals(alpha, bgAlpha, n, currentEnd  , residuals);
  currentEnd += n;
  
  const float * beta   = blockedlogs->getBeta();
  const float * bgBeta = blockedlogs->getBetaHighCutBackground();
  if(useVs == true) {
    MakeInterpolatedResiduals(beta, bgBeta, n, currentEnd, residuals);
    currentEnd += n;
  }
  const float * rho      = blockedlogs->getRho();
  const float * bgRho    = blockedlogs->getRhoHighCutBackground();
  MakeInterpolatedResiduals(rho, bgRho, n, currentEnd, residuals);

  lib_matr_prod(Aw, residuals, nLogs*n, nLogs*n, 1, filterval);

  float * alphaFiltered = new float[n];
  float * betaFiltered  = new float[n];
  float * rhoFiltered   = new float[n];

  for(i=0;i<n;i++)
  {
    int offset = 0;
    if(alpha[i] == RMISSING)
      alphaFiltered[i] = 0.0;
    else
      alphaFiltered[i] = float(filterval[i+offset][0]);
    offset += n;

    if(useVs == true) {
      if(beta[i] == RMISSING)
        betaFiltered[i] = 0.0;
      else
        betaFiltered[i] = float(filterval[i + offset][0]);
      offset += n;
    }

    if(rho[i] == RMISSING)
      rhoFiltered[i] = 0.0;
    else
      rhoFiltered[i] = float(filterval[i + offset][0]);
  }

  if(useVs == true) {
    blockedlogs->setSpatialFilteredLogs(alphaFiltered, n, "ALPHA_SEISMIC_RESOLUTION",bgAlpha);
    blockedlogs->setSpatialFilteredLogs(betaFiltered , n, "BETA_SEISMIC_RESOLUTION" ,bgBeta);
    blockedlogs->setSpatialFilteredLogs(rhoFiltered  , n, "RHO_SEISMIC_RESOLUTION"  ,bgRho);
  }
  else {
    blockedlogs->setSpatialFilteredLogs(alphaFiltered, n, "ALPHA_FOR_FACIES",bgAlpha);
    blockedlogs->setSpatialFilteredLogs(rhoFiltered  , n, "RHO_FOR_FACIES"  ,bgRho);
  }

  delete [] alphaFiltered;
  delete [] betaFiltered;
  delete [] rhoFiltered;

  for(i=0;i<nLogs*n;i++)
  {
    delete [] residuals[i];
    delete [] filterval[i];
  }

  delete [] filterval;
  delete [] residuals;
}


void SpatialWellFilter::MakeInterpolatedResiduals(const float * bwLog, 
                                                  const float * bwLogBG,
                                                  const int     n,
                                                  const int     offset,
                                                  double     ** residuals)
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
      residuals[offset + i][0] = first_residual;
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
      residuals[offset + i][0] = res_i;
      
      if(nmiss>0)
      {
        for(int j=1 ; j<=nmiss ; j++)
        {
          double w = static_cast<double>(j)/static_cast<double>(nmiss + 1);
          residuals[offset + i - j][0] *= w;
          residuals[offset + i - j][0] += (1.0 - w)*res_i;
        }
      }
      nmiss = 0;
    }
    else
    {
      nmiss++;
      residuals[offset + i][0] = residuals[offset + i - 1][0];
    }
  }
}


// The variances used for smootihng in faciesprob might be very small. 
// Therefore eigenvalues are adjusted in order to be able to invert matrix.
void SpatialWellFilter::adjustDiagSigma(double ** sigmae, int n)
{
  double    eps    = 0.0001;
  double  * eigval = new double[n];
  int     * error  = new int[1];
  double ** eigvec = new double *[n];
  for(int i=0;i<n;i++)
    eigvec[i] = new double[n];

  lib_matr_eigen(sigmae,n,eigvec,eigval,error);
  delete [] error;

  for(int i=1;i<n;i++)
    if(eigval[i]/eigval[0]<eps)
      eigval[i]=eps*eigval[0];

  double ** help        = new double *[n];
  double ** eigvalmat   = new double *[n];
  double ** eigvectrans = new double *[n];
  for(int i=0;i<n;i++)
  {
    help[i]        = new double[n];
    eigvalmat[i]   = new double[n];
    eigvectrans[i] = new double[n];
  }
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      if(i==j) 
        eigvalmat[i][j] = eigval[i];
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,n,n,n,help);
  lib_matrTranspose(eigvec,n,n,eigvectrans);
  lib_matr_prod(help,eigvectrans,n,n,n,sigmae);
  
  for(int i=0;i<n;i++)
  {
    delete [] eigvec[i];
    delete [] help[i];
    delete [] eigvectrans[i];
    delete [] eigvalmat[i];
  }
  delete [] eigval;
  delete [] eigvec;
  delete [] help;
  delete [] eigvectrans;
  delete [] eigvalmat;
}
