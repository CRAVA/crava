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
#include "src/definitions.h"
#include "src/welldata.h"
#include "src/timings.h"
#include "src/modelsettings.h"
#include "src/corr.h"
#include "src/crava.h"

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

  for(j=0; j<static_cast<int>(sigmaeSynt_.size()); j++) {
    double ** Se = sigmaeSynt_[j];
    for(i=0;i<3;i++)
      delete [] Se[i];
    delete [] Se;
  }

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


void SpatialWellFilter::setPriorSpatialCorrSyntWell(FFTGrid             * parSpatialCorr,
                                                    SyntWellData        * well,
                                                    int                   wellnr)
{

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
}

void SpatialWellFilter::setPriorSpatialCorr(FFTGrid *parSpatialCorr, WellData *well, int wellnr)
{
  int n = well->getBlockedLogsOrigThick()->getNumberOfBlocks();
  priorSpatialCorr_[wellnr] = new double *[n];
  n_[wellnr] = n;
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

void SpatialWellFilter::doFilteringSyntWells(Corr                                     * corr,
                                             std::vector<SyntWellData *>              & syntWellData,
                                             const std::vector<std::vector<double> >  & v,
                                             int                                        nWells)
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
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, corr->getPostCovAlpha(), n, 0, 0);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, corr->getPostCovBeta(), n, n, n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, corr->getPostCovRho(), n, 2*n, 2*n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, corr->getPostCrCovAlphaBeta(), n, 0, n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, corr->getPostCrCovAlphaRho(), n, 0, 2*n);
    fillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, corr->getPostCrCovBetaRho(), n, 2*n, n);

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
        sigmapri [l1      ][l2      ] = corr->getPriorVar0()[0][0]*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + n  ][l2 + n  ] = corr->getPriorVar0()[1][1]*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + 2*n][l2 + 2*n] = corr->getPriorVar0()[2][2]*priorSpatialCorr_[w1][l1][l2];
        if(l1==l2){
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
    // finds the scale at  default inversion (all minimum noise in case of local noise)
    double ** Se = sigmaeSynt_[0];
    Se[0][0] /= lastn;
    Se[1][0] /= lastn;
    Se[1][1] /= lastn;
    Se[2][0] /= lastn;
    Se[2][1] /= lastn;
    Se[2][2] /= lastn;
    Se[0][1]  = Se[1][0];
    Se[0][2]  = Se[2][0];
    Se[1][2]  = Se[2][1];
    adjustDiagSigma(Se, 3);
  }

  if (no_wells_filtered) {
    LogKit::LogFormatted(LogKit::Low,"\nNo synthetic wells have been filtered.\n");
  }

  Timings::setTimeFiltering(wall,cpu);
}

void SpatialWellFilter::doFiltering(Corr                        * corr,
                                    std::vector<WellData  *>      wells,
                                    int                           nWells,
                                    bool                          useVpRhoFilter,
                                    int                           nAngles,
                                    const Crava                 * cravaResult,
                                    const std::vector<Grid2D *> & noiseScale)
{
  LogKit::WriteHeader("Creating spatial multi-parameter filter");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  double ** sigmapost;
  double ** sigmapri;
  double ** imat;
  double ** Aw;

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

  bool no_wells_filtered = true;

  for(int w1=0;w1<nWells;w1++)
  {
    n = wells[w1]->getBlockedLogsOrigThick()->getNumberOfBlocks();

    if (wells[w1]->getUseForFiltering() == true)
    {
      LogKit::LogFormatted(LogKit::Low,"\nFiltering well "+wells[w1]->getWellname());
      no_wells_filtered = false;

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
                         wells[w1]->getBlockedLogsOrigThick()); //Must do before Cholesky of sigmapri.

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

      LogKit::LogFormatted(LogKit::Low,"\n  Cholesky decomposition ...");
      lib_matrCholR(3*n, sigmapri);
      LogKit::LogFormatted(LogKit::Low,"\n  Equation solving ...");
      lib_matrAXeqBMatR(3*n, sigmapri, imat, 3*n);
      LogKit::LogFormatted(LogKit::Low,"\n  Matrix multiplication ...\n");
      lib_matr_prod(sigmapost,imat,3*n,3*n,3*n,Aw);
      //NBNB-PAL: Aktiver etterhvert
      //lib_matr_prod_sym(sigmapost,imat,3*n,3*n,3*n,Aw);

      for(int i=0;i<3*n;i++) {
        for(int j=0;j<3*n;j++)
        {
          Aw[i][j] *=-1.0;
          if(i==j)
            Aw[i][j]+=1.0;
        }
      }

      if(useVpRhoFilter == false) //Save time, since below is not needed then.
        updateSigmaE(Aw, sigmapost, n);

      calculateFilteredLogs(Aw, wells[w1]->getBlockedLogsOrigThick(), n, true);

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
  }

  if(no_wells_filtered == false)
    completeSigmaE(lastn, cravaResult,noiseScale);

  if(useVpRhoFilter == true)
    completeSigmaEVpRho(lastn,cravaResult,noiseScale);

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
    double ** Se = sigmaeSynt_[0];
    Se[0][0] += sigmaeW[i      ][i     ];
    Se[1][0] += sigmaeW[i +   n][i     ];
    Se[2][0] += sigmaeW[i + 2*n][i     ];
    Se[1][1] += sigmaeW[i +   n][i +  n];
    Se[2][1] += sigmaeW[i + 2*n][i +  n];
    Se[2][2] += sigmaeW[i + 2*n][i + 2*n];
  }
  // sigmaeSynt_ Is normalized (1/n) in completeSigmaE, Here well by well is added.

  for(int i=0;i<3*n;i++)
    delete [] sigmaeW[i];

  delete [] sigmaeW;
}

void
SpatialWellFilter::updateSigmaE(double ** filter, double ** postCov,  int n)
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
  // sigmae_ Is normalized (1/n) in completeSigmaE, Here well by well is added.

  for(int i=0;i<3*n;i++)
    delete [] sigmaeW[i];

  delete [] sigmaeW;
}

void
SpatialWellFilter::completeSigmaE(int lastn, const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale)
{
  // finds the scale at  default inversion (all minimum noise in case of local noise)
  sigmae_[0][0][0] /= lastn;
  sigmae_[0][1][0] /= lastn;
  sigmae_[0][1][1] /= lastn;
  sigmae_[0][2][0] /= lastn;
  sigmae_[0][2][1] /= lastn;
  sigmae_[0][2][2] /= lastn;
  sigmae_[0][0][1]  = sigmae_[0][1][0];
  sigmae_[0][0][2]  = sigmae_[0][2][0];
  sigmae_[0][1][2]  = sigmae_[0][2][1];
  adjustDiagSigma(sigmae_[0], 3);


  if(sigmae_.size() > 1) { // then we have local noise
    // initialization
    int nAng = static_cast<int>(noiseScale.size());
    std::vector<double> maxScale(nAng);
    double maxS, minS;
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      maxScale[angle] = maxS/minS;
    }

    // initialize
    float **sigmaPri0 = cravaResult->getPriorVar0();// not to be deleted

    double ** dummy   = new double * [3];
    double ** help   = new double * [3];
    double ** filter  = new double * [3];
    double ** postCovAdj = new double * [3];
    double ** sigmaEAdj = new double * [3];
    double ** sigmaETmp = new double * [3];
    double ** sigmaE0 = new double * [3];

    for(int i=0;i<3;i++) {
      help[i]       = new double[3];
      dummy[i]      = new double[3];
      filter[i]     = new double[3];
      postCovAdj[i] = new double[3];
      sigmaEAdj[i]  = new double[3];
      sigmaE0[i]    = new double[3];
      sigmaETmp[i]  = new double[3];
    }

    double **G = new double*[nAng];
    for(int i=0;i<nAng;i++)
      G[i] = new double[3];

    std::vector<double> scale(nAng);
    for(int angle=0;angle<nAng;angle++)
      scale[angle] = 1.0;

    cravaResult->computeG(G);
    cravaResult->newPosteriorCovPointwise(dummy, G, scale, help);
    lib_matr_prod(help,help, 3, 3, 3, postCovAdj);
    cravaResult->computeFilter(sigmaPri0,postCovAdj,3,filter);

    lib_matr_prod(filter, postCovAdj, 3, 3, 3, sigmaE0);
    // note idealy  sigmaE0 =sigmae_[0] but this is not the case due to
    // spatial effects.  We therfor adjust the sigmae_[0]
    // by   sqrt(sigmaETmp*sigmaE0^-1)*sigmae_[0]* sqrt(sigmaE0^-1*sigmaETmp)

    //Interpret conf counter bitwise - 0 means min value noise for that component, 1 means max.
    for(unsigned int conf = 1; conf < sigmae_.size();conf++) {

      //Compute pointwise filter
      int factor = 1;
      for(int angle=0;angle<nAng;angle++) {
        if((conf & factor) > 0)
          scale[angle] = maxScale[angle];
        else
          scale[angle] = 1.0;
        factor *= 2;
      }

      cravaResult->newPosteriorCovPointwise(dummy, G, scale,  help);  // help is the square root of postCovAdj
      lib_matr_prod(help,help, 3, 3, 3, postCovAdj);
      cravaResult->computeFilter(sigmaPri0,postCovAdj,3,filter);
      lib_matr_prod(filter, postCovAdj, 3, 3, 3, sigmaETmp);

      computeSigmaEAdjusted(sigmae_[0],sigmaE0,sigmaETmp,3,sigmaEAdj);

      sigmae_[conf][0][0] = sigmaEAdj[0][0];
      sigmae_[conf][0][1] = sigmaEAdj[0][1];
      sigmae_[conf][0][2] = sigmaEAdj[0][2];
      sigmae_[conf][1][0] = sigmaEAdj[1][0];
      sigmae_[conf][1][1] = sigmaEAdj[1][1];
      sigmae_[conf][1][2] = sigmaEAdj[1][2];
      sigmae_[conf][2][0] = sigmaEAdj[2][0];
      sigmae_[conf][2][1] = sigmaEAdj[2][1];
      sigmae_[conf][2][2] = sigmaEAdj[2][2];
      adjustDiagSigma(sigmae_[conf], 3);

    }

    for(int i=0;i<3;i++) {
      delete [] dummy[i];
      delete [] help[i];
      delete [] filter[i];
      delete [] postCovAdj[i];
      delete [] sigmaEAdj[i];
      delete [] sigmaE0[i];
      delete [] sigmaETmp[i];
    }
    delete [] dummy;
    delete [] help;
    delete [] filter;
    delete [] postCovAdj;
    delete [] sigmaEAdj;
    delete [] sigmaE0;
    delete [] sigmaETmp;

    for(int i=0;i<nAng;i++)
      delete [] G[i];
    delete [] G;

  }
}


void
SpatialWellFilter::computeSigmaEAdjusted(double** sigmae ,double** sigmaE0,double** sigmaETmp,int n,double** sigmaEAdj)
{
  // sigmaEAdj  =sqrt(sigmaETmp*sigmaE0^-1)*sigmae* sqrt(sigmaE0^-1*sigmaETmp)

  int     * error         = new int[1];
  double  * eigval0       = new double[n];
  double ** eigvec0       = new double * [n];
  double  * eigvalTmp       = new double[n];
  double ** eigvecTmp       = new double * [n];
  double ** tmp1    = new double * [n];
  double ** tmp2    = new double * [n];
  double ** tmp3    = new double * [n];

  for(int i=0;i<n;i++)
  {
    eigvec0[i]   = new double[n];
    eigvecTmp[i] = new double[n];
    tmp1[i]      = new double[n];
    tmp2[i]      = new double[n];
    tmp3[i]      = new double[n];
  }


  lib_matr_eigen(sigmaE0,n,eigvec0,eigval0,error);
  lib_matr_eigen(sigmaETmp,n,eigvecTmp,eigvalTmp,error);

  double max0=0.0;
  double maxTmp=0.0;
  for(int i=0;i<n;i++)
  {
    max0=std::max(max0,eigval0[i]);
    maxTmp=std::max(maxTmp,eigvalTmp[i]);
  }
  for(int i=0;i<n;i++)
  {
    eigval0[i]=std::max(eigval0[i],max0/1000);
    eigvalTmp[i]=std::max(eigvalTmp[i],maxTmp/1000);
  }

  // computes: =sqrt(sigmaETmp*sigmaE0^-1) = UTmp*LambdaTmp^(1/2)*UTmpT*U0*Lambda0^(-1/2)*U0T

  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      tmp1[i][j] = eigvecTmp[i][j]*sqrt(eigvalTmp[j]);
  // tmp1=UTmp*LambdaTmp^(1/2)

  lib_matrTranspose(eigvecTmp, n, n, tmp2);
  lib_matr_prod(tmp1,tmp2,n,n,n,tmp3);
  // tmp3=UTmp*LambdaTmp^(1/2)*UTmpT

  lib_matr_prod(tmp3,eigvec0,n,n,n,tmp1);
  // tmp1 = UTmp*LambdaTmp^(1/2)*UTmpT*U0

  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      tmp1[i][j] = tmp1[i][j]*sqrt(1./eigval0[j]);
  //tmp1=UTmp*LambdaTmp^(1/2)*UTmpT*U0*Lambda0^(-1/2)

  lib_matrTranspose(eigvec0, n, n, tmp2);
  lib_matr_prod(tmp1,tmp2,n,n,n,tmp3);
  //tmp3=UTmp*LambdaTmp^(1/2)*UTmpT*U0*Lambda0^(-1/2)*U0T

  lib_matr_prod(tmp3,sigmae,n,n,n,tmp1);
  // tmp1 = sqrt(sigmaETmp*sigmaE0^-1)*sigmae

  lib_matrTranspose(tmp3, n, n, tmp2);
  lib_matr_prod(tmp1,tmp2,n,n,n,sigmaEAdj);
  // sigmaEAdj  =sqrt(sigmaETmp*sigmaE0^-1)*sigmae* sqrt(sigmaE0^-1*sigmaETmp)

  for(int i=0;i<n;i++)
  {
    delete [] eigvec0[i];
    delete [] eigvecTmp[i];
    delete [] tmp1[i];
    delete [] tmp2[i];
    delete [] tmp3[i];
  }
  delete [] error;
  delete [] eigval0 ;
  delete [] eigvec0;
  delete [] eigvalTmp;
  delete [] eigvecTmp;
  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp3;
}



void
SpatialWellFilter::doVpRhoFiltering(const double ** sigmapri, const double ** sigmapost, int n,
                                    BlockedLogs * blockedLogs)
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
  for(int i=0;i<2*n;i++) {
    for(int j=0;j<2*n;j++)
    {
      Aw[i][j] *=-1.0;
      if(i==j)
        Aw[i][j]+=1.0;
    }
  }
  calculateFilteredLogs(Aw, blockedLogs, n, false);

  updateSigmaEVpRho(Aw, sigmapost2, static_cast<int>(sigmae_.size()), n);

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
                                     double ** postCov,
                                     int nDim,
                                     int n)
{
  if(sigmaeVpRho_.size() == 0) { // then first time alocate memory
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


  for(int i=0;i<2*n;i++)
    delete [] sigma[i];
  delete [] sigma;
}

void
SpatialWellFilter::completeSigmaEVpRho(int lastn, const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale)
{
  sigmaeVpRho_[0][0][0] /= lastn;
  sigmaeVpRho_[0][1][0] /= lastn;
  sigmaeVpRho_[0][1][1] /= lastn;
  sigmaeVpRho_[0][0][1]  = sigmae_[0][1][0];

  adjustDiagSigma(sigmaeVpRho_[0], 2);

  if(sigmaeVpRho_.size() > 1) { // then we have local noise
    // initialization
    int nAng = static_cast<int>(noiseScale.size());
    std::vector<double> maxScale(nAng);
    double maxS, minS;
    for(int angle=0;angle<nAng;angle++) {
      minS = noiseScale[angle]->FindMin(RMISSING);
      maxS = noiseScale[angle]->FindMax(RMISSING);
      maxScale[angle] = maxS/minS;
    }

    double ** dummy   = new double * [3];
    double ** help   = new double * [3];
    double ** postCovAdj = new double * [3];

    double ** sigmaEAdj = new double * [2];
    double ** sigmaETmp = new double * [2];
    double ** sigmaE0 = new double * [2];
    float  ** priCovVpRho= new float * [2];
    double ** postCovVpRho= new double * [2];
    double ** filter  = new double * [2];

    for(int i=0;i<3;i++) {
      help[i]       = new double[3];
      dummy[i]      = new double[3];
      postCovAdj[i] = new double[3];
    }

    for(int i=0;i<2;i++) {
      priCovVpRho[i] = new float[2];
      postCovVpRho[i]= new double[2];
      filter[i]      = new double[2];
      sigmaEAdj[i]   = new double[2];
      sigmaE0[i]     = new double[2];
      sigmaETmp[i]   = new double[2];
    }

    float **sigmaPri0 = cravaResult->getPriorVar0();// not to be deleted
    priCovVpRho[0][0]=sigmaPri0[0][0];
    priCovVpRho[0][1]=sigmaPri0[0][2];
    priCovVpRho[1][0]=sigmaPri0[2][0];
    priCovVpRho[1][1]=sigmaPri0[2][2];


    double **G = new double*[nAng];
    for(int i=0;i<nAng;i++)
      G[i] = new double[3];

    std::vector<double> scale(nAng);
    for(int angle=0;angle<nAng;angle++)
      scale[angle] = 1.0;

    cravaResult->computeG(G);
    cravaResult->newPosteriorCovPointwise(dummy, G, scale, help);
    lib_matr_prod(help,help, 3, 3, 3, postCovAdj);
    postCovVpRho[0][0]=postCovAdj[0][0];
    postCovVpRho[0][1]=postCovAdj[0][2];
    postCovVpRho[1][0]=postCovAdj[2][0];
    postCovVpRho[1][1]=postCovAdj[2][2];

    cravaResult->computeFilter(priCovVpRho,postCovVpRho,2,filter);
    lib_matr_prod(filter, postCovVpRho, 2, 2, 2, sigmaE0);

    for(unsigned int conf = 1; conf < sigmaeVpRho_.size();conf++) {

      //Compute pointwise filter
      int factor = 1;
      for(int angle=0;angle<nAng;angle++) {
        if((conf & factor) > 0)
          scale[angle] = maxScale[angle];
        else
          scale[angle] = 1.0;
        factor *= 2;
      }
      cravaResult->newPosteriorCovPointwise(dummy, G, scale,  help);  // help is the square root of postCovAdj
      lib_matr_prod(help,help, 3, 3, 3, postCovAdj);
      postCovVpRho[0][0]=postCovAdj[0][0];
      postCovVpRho[0][1]=postCovAdj[0][2];
      postCovVpRho[1][0]=postCovAdj[2][0];
      postCovVpRho[1][1]=postCovAdj[2][2];

      cravaResult->computeFilter( priCovVpRho,postCovVpRho,2,filter);
      lib_matr_prod(filter, postCovVpRho, 2, 2, 2, sigmaETmp);

      computeSigmaEAdjusted(sigmaeVpRho_[0],sigmaE0,sigmaETmp,2,sigmaEAdj);

      sigmaeVpRho_[conf][0][0] = sigmaEAdj[0][0];
      sigmaeVpRho_[conf][0][1] = sigmaEAdj[0][1];
      sigmaeVpRho_[conf][1][0] = sigmaEAdj[1][0];
      sigmaeVpRho_[conf][1][1] = sigmaEAdj[1][1];
      adjustDiagSigma(sigmaeVpRho_[conf], 2);
    }
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
