/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <float.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "lib/lib_matr.h"
#include "lib/timekit.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

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

  delete [] n_;
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

//--------------------------------------------------------------------------------
void SpatialWellFilter::doFiltering(Corr                        * corr,
                                    WellData                   ** wells,
                                    int                           nWells,
                                    bool                          useVpRhoFilter,
                                    int                           nAngles,
                                    const Crava                 * cravaResult,
                                    const std::vector<Grid2D *> & noiseScale)
//-------------------------------------------------------------------------------
{
  LogKit::WriteHeader("Creating spatial multi-parameter filter");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

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

  bool no_wells_filtered = true;

  for(int w1=0 ; w1 < nWells ; w1++)
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

      const int *ipos = wells[w1]->getBlockedLogsOrigThick()->getIpos();
      const int *jpos = wells[w1]->getBlockedLogsOrigThick()->getJpos();
      const int *kpos = wells[w1]->getBlockedLogsOrigThick()->getKpos();

      float regularization = Definitions::SpatialFilterRegularisationValue();

      fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCovAlpha()      , n, 0  , 0   );
      fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCovBeta()       , n, n  , n   );
      fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCovRho()        , n, 2*n, 2*n );
      fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCrCovAlphaBeta(), n, 0  , n   );
      fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCrCovAlphaRho() , n, 0  , 2*n );
      fillValuesInSigmapost(sigmapost, ipos, jpos, kpos, corr->getPostCrCovBetaRho()  , n, 2*n, n   );

      for(int l1=0 ; l1 < n ; l1++) {
        for(int l2=0 ; l2 < n ; l2++) {
          sigmapost[l2 + n  ][l1      ] = sigmapost[l1      ][l2 + n  ];
          sigmapost[l2 + 2*n][l1      ] = sigmapost[l1      ][l2 + 2*n];
          sigmapost[l2 + n  ][l1 + 2*n] = sigmapost[l1 + 2*n][l2 + n  ];
          sigmapri [l1      ][l2      ] = corr->getPriorVar0()[0][0]*priorSpatialCorr_[w1][l1][l2];
          sigmapri [l1 + n  ][l2 + n  ] = corr->getPriorVar0()[1][1]*priorSpatialCorr_[w1][l1][l2];
          sigmapri [l1 + 2*n][l2 + 2*n] = corr->getPriorVar0()[2][2]*priorSpatialCorr_[w1][l1][l2];
          if(l1==l2)
          {
            sigmapost[l1      ][l2      ] += regularization*sigmapost[l1      ][l2      ]/sigmapri[l1      ][l2      ];
            sigmapost[l1 + n  ][l2 + n  ] += regularization*sigmapost[l1 + n  ][l2 + n  ]/sigmapri[l1 + n  ][l2 + n  ];
            sigmapost[l1 + 2*n][l2 + 2*n] += regularization*sigmapost[l1 + 2*n][l2 + 2*n]/sigmapri[l1 + 2*n][l2 + 2*n];
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

      NRLib::SymmetricMatrix Sprior(3*n);
      NRLib::SymmetricMatrix Spost(3*n);
      for(int i = 0 ; i < 3*n ; i++)
        for(int j = 0 ; j <= i ; j++)
          Sprior(j,i) = sigmapri[j][i];

      for(int i = 0 ; i < 3*n ; i++)
        for(int j = 0 ; j <= i ; j++)
          Spost(j,i) = sigmapost[j][i];

      if(useVpRhoFilter == true) //Only additional
        doVpRhoFiltering(Sprior,
                         Spost,
                         n,
                         wells[w1]->getBlockedLogsOrigThick()); //Must do before Cholesky of sigmapri.


      NRLib::Matrix I(3*n, 3*n);
      NRLib::InitializeMatrix(I, 0.0);
      for(int i=0 ; i<3*n ; i++)
        I(i,i) = 1.0;

      //
      // Filter = I - Sigma_post * inv(Sigma_prior)
      //
      NRLib::Matrix Aw;
      NRLib::CholeskySolve(Sprior, I);
      Aw = Spost * I;
      Aw = Aw * (-1);
      for(int i=0 ; i<3*n ; i++) {
        Aw(i,i) += 1.0;
      }

      if(useVpRhoFilter == false) { //Save time, since below is not needed then.
        updateSigmaE(Aw,
                     Spost,
                     sigmae_[0],
                     n);
      }

      calculateFilteredLogs(Aw,
                            wells[w1]->getBlockedLogsOrigThick(),
                            n,
                            true);

      lastn += n;

      for(int i=0;i<3*n;i++)
      {
        delete [] sigmapost[i];
      }
      delete [] sigmapri;
      delete [] sigmapost;
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

//------------------------------------------------------------------
void SpatialWellFilter::updateSigmaE(const NRLib::Matrix & Filter,
                                     const NRLib::Matrix & PostCov,
                                     NRLib::Matrix       & sigmae,
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
  // sigmae_ Is normalized (1/n) in completeSigmaE, Here well by well is added.
}

//-------------------------------------------------------------------------------
void SpatialWellFilter::completeSigmaE(int                           lastn,
                                       const Crava                 * cravaResult,
                                       const std::vector<Grid2D *> & noiseScale)
//-------------------------------------------------------------------------------
{
  // finds the scale at  default inversion (all minimum noise in case of local noise)
  sigmae_[0](0,0) /= lastn;
  sigmae_[0](1,0) /= lastn;
  sigmae_[0](1,1) /= lastn;
  sigmae_[0](2,0) /= lastn;
  sigmae_[0](2,1) /= lastn;
  sigmae_[0](2,2) /= lastn;
  sigmae_[0](0,1)  = sigmae_[0](1,0);
  sigmae_[0](0,2)  = sigmae_[0](2,0);
  sigmae_[0](1,2)  = sigmae_[0](2,1);

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

    double ** dummy      = new double * [3];
    double ** help       = new double * [3];
    double ** filter     = new double * [3];
    double ** postCovAdj = new double * [3];
    double ** sigmaEAdj  = new double * [3];
    double ** sigmaETmp  = new double * [3];
    double ** sigmaE0    = new double * [3];

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

      sigmae_[conf](0,0) = sigmaEAdj[0][0];
      sigmae_[conf](0,1) = sigmaEAdj[0][1];
      sigmae_[conf](0,2) = sigmaEAdj[0][2];
      sigmae_[conf](1,0) = sigmaEAdj[1][0];
      sigmae_[conf](1,1) = sigmaEAdj[1][1];
      sigmae_[conf](1,2) = sigmaEAdj[1][2];
      sigmae_[conf](2,0) = sigmaEAdj[2][0];
      sigmae_[conf](2,1) = sigmaEAdj[2][1];
      sigmae_[conf](2,2) = sigmaEAdj[2][2];

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
SpatialWellFilter::computeSigmaEAdjusted(NRLib::Matrix &  sigmaeX,
                                         double        ** sigmaE0,
                                         double        ** sigmaETmp,
                                         int              n,
                                         double        ** sigmaEAdj)
{
  // sigmaEAdj  =sqrt(sigmaETmp*sigmaE0^-1)*sigmae* sqrt(sigmaE0^-1*sigmaETmp)



  double ** sigmae = new double*[n];
  for (int i = 0 ; i < n ; i++) {
    sigmae[i] = new double[n];
  }

 NRLib::Set2DArrayFromMatrix(sigmaeX, sigmae);


  int     * error     = new int[1];
  double  * eigval0   = new double[n];
  double ** eigvec0   = new double * [n];
  double  * eigvalTmp = new double[n];
  double ** eigvecTmp = new double * [n];
  double ** tmp1      = new double * [n];
  double ** tmp2      = new double * [n];
  double ** tmp3      = new double * [n];

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

  for (int i = 0 ; i < n ; i++) {
    delete [] sigmae[i];
  }
  delete [] sigmae;

}

//----------------------------------------------------------------------------------
void SpatialWellFilter::doVpRhoFiltering(const NRLib::SymmetricMatrix & sigmapri,
                                         const NRLib::SymmetricMatrix & sigmapost,
                                         const int                      n,
                                         BlockedLogs                  * blockedLogs)
//----------------------------------------------------------------------------------
{
  int m = 2*n;
  NRLib::SymmetricMatrix Sprior2(m);
  NRLib::SymmetricMatrix Spost2(m);

  for (int j=0 ; j<n ; j++) {
    for (int i=0 ; i<j ; i++) {
      Sprior2(i,   j  ) = sigmapri(i,   j  );
      Sprior2(i+n, j  ) = sigmapri(i+m, j  );
      Sprior2(i  , j+n) = sigmapri(i  , j+m);
      Sprior2(i+n, j+n) = sigmapri(i+m, j+m);

      Spost2(i,   j  )  = sigmapost(i  , j  );
      Spost2(i+n, j  )  = sigmapost(i+m, j  );
      Spost2(i,   j+n)  = sigmapost(i  , j+m);
      Spost2(i+n, j+n)  = sigmapost(i+m, j+m);
    }
  }

  NRLib::Matrix I(m, m);
  NRLib::InitializeMatrix(I, 0);
  for(int i=0 ; i<m ; i++)
    I(i,i) = 1.0;

  NRLib::Matrix Aw;
  NRLib::CholeskySolve(Sprior2, I);
  Aw = Spost2 * I;
  Aw = Aw * (-1);
  for(int i=0 ; i<m ; i++) {
    Aw(i,i) += 1.0;
  }

  calculateFilteredLogs(Aw, blockedLogs, n, false);

  updateSigmaEVpRho(Aw,
                    Spost2,
                    static_cast<int>(sigmae_.size()),
                    n);
}

//---------------------------------------------------------------------
void SpatialWellFilter::updateSigmaEVpRho(const NRLib::Matrix & Aw,
                                          const NRLib::Matrix & Spost,
                                          int                   nDim,
                                          int                   n)
//---------------------------------------------------------------------
{
  if (sigmaeVpRho_.size() == 0) { // then first time alocate memory
    sigmaeVpRho_.resize(nDim);

    NRLib::Matrix sigmae(2, 2);
    NRLib::InitializeMatrix(sigmae, 0.0);

    for (int k=0 ; k < nDim ; k++) {
      sigmaeVpRho_[k] = sigmae;
    }
  }

  NRLib::Matrix sigma = Aw * Spost;

  //
  // NBNB-PAL: Bug? førsteindeksen på sigmaeVpRho_[0][0][0] står
  // stille hele tiden. Det er ingen n-avhengighet.
  //
  for(int i=0 ; i < n ; i++) {
    sigmaeVpRho_[0](0,0) += sigma(i    , i    );
    sigmaeVpRho_[0](1,0) += sigma(i + n, i    );
    sigmaeVpRho_[0](1,1) += sigma(i + n, i + n);
  }
}

//------------------------------------------------------------------------------------
void SpatialWellFilter::completeSigmaEVpRho(int                           lastn,
                                            const Crava                 * cravaResult,
                                            const std::vector<Grid2D *> & noiseScale)
//------------------------------------------------------------------------------------
{
  sigmaeVpRho_[0](0,0) /= lastn;
  sigmaeVpRho_[0](1,0) /= lastn;
  sigmaeVpRho_[0](1,1) /= lastn;
  sigmaeVpRho_[0](0,1)  = sigmae_[0](1,0);

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

    double ** dummy        = new double * [3];
    double ** help         = new double * [3];
    double ** postCovAdj   = new double * [3];

    double ** sigmaEAdj    = new double * [2];
    double ** sigmaETmp    = new double * [2];
    double ** sigmaE0      = new double * [2];
    float  ** priCovVpRho  = new float  * [2];
    double ** postCovVpRho = new double * [2];
    double ** filter       = new double * [2];

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

      sigmaeVpRho_[conf](0,0) = sigmaEAdj[0][0];
      sigmaeVpRho_[conf](0,1) = sigmaEAdj[0][1];
      sigmaeVpRho_[conf](1,0) = sigmaEAdj[1][0];
      sigmaeVpRho_[conf](1,1) = sigmaEAdj[1][1];
      adjustDiagSigma(sigmaeVpRho_[conf], 2);
    }
  }
}

//------------------------------------------------------------------------------
void SpatialWellFilter::calculateFilteredLogs(const NRLib::Matrix & Aw,
                                              BlockedLogs         * blockedlogs,
                                              int                   n,
                                              bool                  useVs)
//------------------------------------------------------------------------------
{
  int nLogs = 2;
  if(useVs == true)
    nLogs++;

  NRLib::Vector residuals(nLogs*n);

  int currentEnd = 0;
  const float * alpha   = blockedlogs->getAlpha();
  const float * bgAlpha = blockedlogs->getAlphaHighCutBackground();
  MakeInterpolatedResiduals(alpha, bgAlpha, n, currentEnd, residuals);
  currentEnd += n;

  const float * beta    = blockedlogs->getBeta();
  const float * bgBeta  = blockedlogs->getBetaHighCutBackground();
  if(useVs == true) {
    MakeInterpolatedResiduals(beta, bgBeta, n, currentEnd, residuals);
    currentEnd += n;
  }
  const float * rho     = blockedlogs->getRho();
  const float * bgRho   = blockedlogs->getRhoHighCutBackground();
  MakeInterpolatedResiduals(rho, bgRho, n, currentEnd, residuals);

  NRLib::Vector filteredVal = Aw * residuals;

  float * alphaFiltered = new float[n];
  float * betaFiltered  = new float[n];
  float * rhoFiltered   = new float[n];

  for(int i=0;i<n;i++)
  {
    int offset = 0;
    if(alpha[i] == RMISSING)
      alphaFiltered[i] = 0.0;
    else
      alphaFiltered[i] = static_cast<float>(filteredVal(i+offset));
    offset += n;

    if(useVs == true) {
      if(beta[i] == RMISSING)
        betaFiltered[i] = 0.0;
      else
        betaFiltered[i] = static_cast<float>(filteredVal(i+offset));
      offset += n;
    }

    if(rho[i] == RMISSING)
      rhoFiltered[i] = 0.0;
    else
      rhoFiltered[i] = static_cast<float>(filteredVal(i+offset));
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
}

//--------------------------------------------------------------------------
void SpatialWellFilter::MakeInterpolatedResiduals(const float   * bwLog,
                                                  const float   * bwLogBG,
                                                  const int       n,
                                                  const int       offset,
                                                  NRLib::Vector & residuals)
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

void SpatialWellFilter::adjustDiagSigma(NRLib::Matrix & sigmae,
                                        int             n)
{
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
