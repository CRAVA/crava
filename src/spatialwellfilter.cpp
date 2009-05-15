#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include "lib/lib_matr.h"
#include "lib/kriging1d.h"
#include "src/welldata.h"
#include "src/spatialwellfilter.h"
#include "src/corr.h"
#include "src/model.h"
#include "nrlib/iotools/logkit.hpp"

SpatialWellFilter::SpatialWellFilter(int nwells)
{
  nWells_ = nwells;
  priorSpatialCorr_ = new double **[nWells_];
  sigmae_ = new double *[3];
  for(int i=0;i<3;i++)
    sigmae_[i] = new double[3];
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      sigmae_[i][j] = 0.0;
}

SpatialWellFilter::~SpatialWellFilter()
{
  delete [] alphaFiltered_;
  delete [] betaFiltered_;
  delete [] rhoFiltered_;
  int i,j;
  for(i=0;i<nWells_;i++)
  {
    for(j=0;j<n_[i];j++)
      delete [] priorSpatialCorr_[i][j];
    delete [] priorSpatialCorr_[i];
  }
  for(i=0;i<3;i++)
    delete [] sigmae_[i];
  delete [] sigmae_;
  delete [] priorSpatialCorr_;
  delete [] n_;


}
void SpatialWellFilter::setPriorSpatialCorr(FFTGrid *parSpatialCorr, WellData *well, int wellnr)
{
  int n = well->getBlockedLogsOrigThick()->getNumberOfBlocks();
  priorSpatialCorr_[wellnr] = new double *[n];
  int i;
  for(i=0;i<n;i++)
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

void SpatialWellFilter::doFiltering(Corr *corr, WellData **wells, int nWells, int relative)
{
  n_ = new int[nWells];
  double **sigmapost;
  double **sigmapri;
  int w1, n;
  double **imat;
  double **Aw;
  nData_ = 0;
  for(w1=0;w1<nWells;w1++)
    nData_ += wells[w1]->getBlockedLogsOrigThick()->getNumberOfBlocks();

  alphaFiltered_ = new float[nData_];
  betaFiltered_ = new float[nData_];
  rhoFiltered_ = new float[nData_];
  
  int lastn = 0;
  n = 0;
  for(w1=0;w1<nWells;w1++)
  {   
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
    float regularization = float(0.0001);
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
        sigmapost[l1][l2] = corr->getPostCovAlpha()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        sigmapost[n+l1][n+l2] = corr->getPostCovBeta()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        sigmapost[2*n+l1][2*n+l2] = corr->getPostCovRho()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        sigmapost[l1][n+l2] = corr->getPostCrCovAlphaBeta()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        sigmapost[n+l2][l1] = sigmapost[l1][n+l2];
        sigmapost[l1][2*n+l2] = corr->getPostCrCovAlphaRho()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);       
        sigmapost[2*n+l2][l1] = sigmapost[l1][2*n+l2];
        sigmapost[2*n+l1][n+l2] = corr->getPostCrCovBetaRho()->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
        sigmapost[n+l2][2*n+l1] = sigmapost[2*n+l1][n+l2];
        sigmapri[l1][l2] = corr->getPriorVar0()[0][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[n+l1][n+l2] = corr->getPriorVar0()[1][1]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[2*n+l1][2*n+l2] = corr->getPriorVar0()[2][2]*priorSpatialCorr_[w1][l1][l2];
        if(l1==l2)
        {
          sigmapost[l1][l2]+=regularization*sigmapost[l1][l2]/sigmapri[l1][l2];
          sigmapost[n+l1][n+l2]+=regularization*sigmapost[n+l1][n+l2]/sigmapri[n+l1][n+l2];
          sigmapost[2*n+l1][2*n+l2]+=regularization*sigmapost[2*n+l1][2*n+l2]/sigmapri[2*n+l1][2*n+l2];
          sigmapri[l1][l2]+=regularization;
          sigmapri[n+l1][n+l2]+=regularization;
          sigmapri[2*n+l1][2*n+l2]+=regularization;
        }
        sigmapri[n+l1][l2] = corr->getPriorVar0()[1][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2][n+l1] = corr->getPriorVar0()[1][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[2*n+l1][l2] = corr->getPriorVar0()[2][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2][2*n+l1] = corr->getPriorVar0()[2][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[n+l1][2*n+l2] = corr->getPriorVar0()[2][1]*priorSpatialCorr_[w1][l1][l2];
        sigmapri[2*n+l2][n+l1] = corr->getPriorVar0()[2][1]*priorSpatialCorr_[w1][l1][l2];
      }
    }


    double **test;
    test = new double * [3*n];
    for(int i=0;i<3*n;i++)
    {
      test[i] = new double[3*n];
    }
   
    lib_matrCholR(3*n, sigmapri);
    lib_matrAXeqBMatR(3*n, sigmapri, imat, 3*n);
    lib_matr_prod(sigmapost,imat,3*n,3*n,3*n,Aw);
    for(int i=0;i<3*n;i++)
      for(int j=0;j<3*n;j++)
      {
        Aw[i][j] *=-1.0;
        if(i==j)
          Aw[i][j]+=1.0;
      }

    lib_matr_prod(Aw,sigmapost,3*n,3*n,3*n,test);
    for(int i=0;i<n;i++)
    {
      sigmae_[0][0] += test[i][i];
      sigmae_[1][0] += test[i+n][i];
      sigmae_[2][0] += test[i+2*n][i];
      sigmae_[1][1] += test[i+n][i+n];
      sigmae_[2][1] += test[i+2*n][i+n];
      sigmae_[2][2] += test[i+2*n][i+2*n];
    }
    
    calculateFilteredLogs(Aw,wells[w1]->getBlockedLogsOrigThick(), n, lastn, relative);
    
    n_[w1] = n;
    lastn = n;
    for(int i=0;i<3*n;i++)
    {
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
  sigmae_[0][0]/=(n+lastn);
  sigmae_[1][0]/=(n+lastn);
  sigmae_[1][1]/=(n+lastn);
  sigmae_[2][0]/=(n+lastn);
  sigmae_[2][1]/=(n+lastn);
  sigmae_[2][2]/=(n+lastn);
  sigmae_[0][1] = sigmae_[1][0];
  sigmae_[0][2] = sigmae_[2][0];
  sigmae_[1][2] = sigmae_[2][1];
    
  adjustDiagSigmae();

}

void SpatialWellFilter::calculateFilteredLogs(double **Aw, BlockedLogs *blockedlogs, int n, int lastn, int relative)
{

  double **filterval;
  filterval = new double *[3*n];
  int i;
  for(i=0;i<3*n;i++)
  filterval[i] = new double[1];
  
  double **log;
  log = new double * [3*n];
  for(i=0;i<3*n;i++)
    log[i] = new double[1];
 

  int nmiss = 0;

  for(i=0;i<n;i++)
  {
    if(blockedlogs->getAlpha()[i]!=RMISSING)
    {
      log[i][0] = double(blockedlogs->getAlpha()[i]-blockedlogs->getAlphaHighCutBackground()[i]);
      if(nmiss>0)
      {
        for(int j=1;j<=nmiss;j++)
        {
          log[i-j][0]*=double(j*1.0/(nmiss+1));
          log[i-j][0]+=double(blockedlogs->getAlpha()[i]-blockedlogs->getAlphaHighCutBackground()[i])*(nmiss+1-j)/(nmiss+1);
        }
      }
      nmiss = 0;
    }
    else if(i-nmiss>=0)
    {
      nmiss++;
      log[i][0] = double(blockedlogs->getAlpha()[i-nmiss]-blockedlogs->getAlphaHighCutBackground()[i-nmiss]);
    }
    else
    {
      nmiss++;
      log[i][0] = 0.0;
    }
  }
  nmiss = 0;
  for(i=0;i<n;i++)
  {
    if(blockedlogs->getBeta()[i]!=RMISSING)
    {
      log[i+n][0] = double(blockedlogs->getBeta()[i]-blockedlogs->getBetaHighCutBackground()[i]);
      if(nmiss>0)
      {
        for(int j=1;j<=nmiss;j++)
        {
          log[i+n-j][0]*=double(j*1.0/(nmiss+1));
          log[i+n-j][0]+=double(blockedlogs->getBeta()[i]-blockedlogs->getBetaHighCutBackground()[i])*(nmiss+1-j)/(nmiss+1);
        }
      }
      nmiss = 0;
    }
    else if(i-nmiss>=0)
    {
      nmiss++;
      log[i+n][0] = double(blockedlogs->getBeta()[i-nmiss]-blockedlogs->getBetaHighCutBackground()[i-nmiss]);
    }
    else
    {
      nmiss++;
      log[i+n][0] = 0.0;
    }
  }
  nmiss = 0;
  for(i=0;i<n;i++)
  {
    if(blockedlogs->getRho()[i]!=RMISSING)
    {
      log[i+2*n][0] = double(blockedlogs->getRho()[i]-blockedlogs->getRhoHighCutBackground()[i]);
      if(nmiss>0)
      {
        for(int j=1;j<=nmiss;j++)
        {
          log[i+2*n-j][0]*=double(j*1.0/(nmiss+1));
          log[i+2*n-j][0]+=double(blockedlogs->getRho()[i]-blockedlogs->getRhoHighCutBackground()[i])*(nmiss+1-j)/(nmiss+1);
        }
      }
      nmiss = 0;
    }
    else if(i-nmiss>=0)
    {
      nmiss++;
      log[i+2*n][0] = double(blockedlogs->getRho()[i-nmiss]-blockedlogs->getRhoHighCutBackground()[i-nmiss]);
    }
    else
    {
      nmiss++;
      log[i+2*n][0] = 0.0;
    }
  }
  lib_matr_prod(Aw,log,3*n,3*n,1,filterval);
  for(i=0;i<n;i++)
  {
    if(blockedlogs->getAlpha()[i]==RMISSING)
      alphaFiltered_[i+lastn] = 0.0;
    else
      alphaFiltered_[i+lastn] = float(filterval[i][0]);
    if(blockedlogs->getBeta()[i]==RMISSING)
      betaFiltered_[i+lastn] = 0.0;
    else
      betaFiltered_[i+lastn] = float(filterval[i+n][0]);
    if(blockedlogs->getRho()[i]==RMISSING)
      rhoFiltered_[i+lastn] = 0.0;
    else
      rhoFiltered_[i+lastn] = float(filterval[i+2*n][0]);
  }
   blockedlogs->setSpatialFilteredLogs(alphaFiltered_,lastn,n+lastn,"ALPHA_SEISMIC_RESOLUTION",blockedlogs->getAlphaHighCutBackground());
   blockedlogs->setSpatialFilteredLogs(betaFiltered_,lastn,n+lastn,"BETA_SEISMIC_RESOLUTION",blockedlogs->getBetaHighCutBackground());
   blockedlogs->setSpatialFilteredLogs(rhoFiltered_,lastn,n+lastn,"RHO_SEISMIC_RESOLUTION",blockedlogs->getRhoHighCutBackground());
  if(relative==0)
    {
      for(i=0;i<n;i++)
      {
        alphaFiltered_[i+lastn]+=blockedlogs->getAlphaHighCutBackground()[i];
        betaFiltered_[i+lastn]+=blockedlogs->getBetaHighCutBackground()[i];
        rhoFiltered_[i+lastn]+=blockedlogs->getRhoHighCutBackground()[i];
      }
    }
 
  for(i=0;i<3*n;i++)
  {
    delete [] log[i];
    delete [] filterval[i];
  }

  delete [] filterval;
  delete [] log;
}

// The variances used for smootihng in faciesprob might be very small. 
// Therefore eigenvalues are adjusted in order to be able to invert matrix.
void SpatialWellFilter::adjustDiagSigmae()
{
  int n = 3;
  double **eigvec = new double *[3];
  int i,j;
  for(i=0;i<3;i++)
    eigvec[i] = new double[3];
  double *eigval = new double[3];
  int *error =new int[1];
  double eps = 0.0001;
  lib_matr_eigen(sigmae_,n,eigvec,eigval,error);
  if(eigval[1]/eigval[0]<eps)
    eigval[1]=eps*eigval[0];
  if(eigval[2]/eigval[0]<eps)
    eigval[2]=eps*eigval[0];
  double **help = new double *[3];
  double  **eigvalmat = new double *[3];
  double **eigvectrans = new double *[3];
  for(i=0;i<3;i++)
  {
    help[i] = new double[3];
    eigvalmat[i] = new double[3];
    eigvectrans[i] = new double [3];
  }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j) 
        eigvalmat[i][j] = eigval[i];
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmae_);

  for(i=0;i<3;i++)
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
