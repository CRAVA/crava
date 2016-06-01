/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/spatialsyntwellfilter.h"
#include "src/spatialwellfilter.h"
#include "nrlib/random/random.hpp"

class FFTGrid;


SpatialSyntWellFilter::SpatialSyntWellFilter()
{
}

SpatialSyntWellFilter::SpatialSyntWellFilter(const std::map<std::string, DistributionsRock *>        & rock_distributions,
                                             const std::vector<std::string>                          & facies_names,
                                             const std::vector<double>                               & trend_min,
                                             const std::vector<double>                               & trend_max,
                                             int                                                       n_synt_wells,
                                             int                                                       n_wells_pr_combination_trend,
                                             double                                                    dz,
                                             int                                                       n_bins_trend,
                                             int                                                       syntWellLength,
                                             bool                                                      cov_estimated)
{
  nWells_                            = n_synt_wells;
  n_bins_trend_                      = n_bins_trend;
  nWellsToBeFiltered_                = n_synt_wells;
  nWellsPerCombinationOfTrendParams_ = n_wells_pr_combination_trend;
  n_                                 = new int[nWells_];

  // If the prior cov has been estimated, we do not use a common corr_t
  if (cov_estimated) {
    prior_cov_vp_.resize(nWells_);
    prior_cov_vs_.resize(nWells_);
    prior_cov_rho_.resize(nWells_);
    prior_cov_vpvs_.resize(nWells_);
    prior_cov_vprho_.resize(nWells_);
    prior_cov_vsrho_.resize(nWells_);
    priorSpatialCorr_ = NULL;
  }
  else {
    priorSpatialCorr_ = new double **[nWells_];
    for(int i = 0 ; i < nWells_ ; i++)
      priorSpatialCorr_[i] = NULL;
  }

  if ((trend_max[0] - trend_min[0]) > 0.0 && (trend_max[1] - trend_min[1]) > 0.0) {

    trend_1_.resize(n_bins_trend_);
    trend_2_.resize(n_bins_trend_);

    trend_1_bin_size_ = (trend_max[0] - trend_min[0])/n_bins_trend_;
    trend_2_bin_size_ = (trend_max[1] - trend_min[1])/n_bins_trend_;

    for (int i = 0; i < n_bins_trend_; i++) {
      trend_1_[i] = trend_min[0] + i*trend_1_bin_size_;
      trend_2_[i] = trend_min[1] + i*trend_2_bin_size_;
    }
    //it is not possible to have a trend2 but no trend1
  }
  else if (trend_max[1] - trend_min[1] == 0.0 && trend_max[0] - trend_min[0] > 0.0) {
    trend_1_.resize(n_bins_trend_);
    trend_2_.resize(1, trend_min[1]);

    trend_1_bin_size_ = (trend_max[0] - trend_min[0]) / n_bins_trend_;

    for (int i = 0; i < n_bins_trend_; i++)
      trend_1_[i] = trend_min[0] + i*trend_1_bin_size_;
    // if trend2max == trend1min and trend1max == trend1min
  }
  else {
    trend_1_.resize(1, trend_min[0]);
    trend_2_.resize(1, trend_min[1]);
  }

  syntWellData_.resize((trend_1_.size() * trend_2_.size() * nWellsPerCombinationOfTrendParams_), NULL);

  GenerateSyntWellData(rock_distributions,
                       facies_names,
                       dz,
                       syntWellLength);
}

SpatialSyntWellFilter::~SpatialSyntWellFilter()
{

  if (priorSpatialCorr_ != NULL) {
    for (int i = 0; i < nWells_; i++) {
      for (int j = 0; j < n_[i]; j++)
        delete [] priorSpatialCorr_[i][j];
      delete [] priorSpatialCorr_[i];
    }
    delete [] priorSpatialCorr_;
  }
  delete [] n_;
}

void  SpatialSyntWellFilter::SetPriorSpatialCovarianceSyntWell(const FFTGrid               * cov_vp,
                                                               const FFTGrid               * cov_vs,
                                                               const FFTGrid               * cov_rho,
                                                               const FFTGrid               * cov_vpvs,
                                                               const FFTGrid               * cov_vprho,
                                                               const FFTGrid               * cov_vsrho,
                                                               int                           wellnr)
{
  int n_blocks = syntWellData_[wellnr]->getWellLength();

  prior_cov_vp_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vs_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_rho_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vpvs_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vprho_[wellnr].Resize(n_blocks, n_blocks, 0);
  prior_cov_vsrho_[wellnr].Resize(n_blocks, n_blocks, 0);

  int i1,j1,k1, i2,j2,k2;
  const int * ipos = syntWellData_[wellnr]->getIpos();
  const int * jpos = syntWellData_[wellnr]->getJpos();
  const int * kpos = syntWellData_[wellnr]->getKpos();
  for(int l1 = 0; l1 < n_blocks; l1++)
  {
    i1 = ipos[l1];
    j1 = jpos[l1];
    k1 = kpos[l1];
    for(int l2 = l1; l2 < n_blocks; l2++)
    {
      i2 = ipos[l2];
      j2 = jpos[l2];
      k2 = kpos[l2];
      // Entries on the diagonal of the prior covariance matrix
      prior_cov_vp_[wellnr](l1,l2) = cov_vp->getRealValueCyclic(i1-i2,j1-j2,k2-k1);
      prior_cov_vp_[wellnr](l2,l1) = prior_cov_vp_[wellnr](l1,l2);
      prior_cov_vs_[wellnr](l1,l2) = cov_vs->getRealValueCyclic(i1-i2,j1-j2,k2-k1);
      prior_cov_vs_[wellnr](l2,l1) = prior_cov_vs_[wellnr](l1,l2);
      prior_cov_rho_[wellnr](l1,l2) = cov_rho->getRealValueCyclic(i1-i2,j1-j2,k2-k1);
      prior_cov_rho_[wellnr](l2,l1) = prior_cov_rho_[wellnr](l1,l2);

      // Cross covariance entries in the prior cov matrix
      prior_cov_vpvs_[wellnr](l1,l2) = cov_vpvs->getRealValueCyclic(i1-i2,j1-j2,k2-k1);
      prior_cov_vpvs_[wellnr](l2,l1) = cov_vpvs->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vprho_[wellnr](l1,l2) = cov_vprho->getRealValueCyclic(i1-i2,j1-j2,k2-k1);
      prior_cov_vprho_[wellnr](l2,l1) = cov_vprho->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
      prior_cov_vsrho_[wellnr](l1,l2) = cov_vsrho->getRealValueCyclic(i1-i2,j1-j2,k2-k1);
      prior_cov_vsrho_[wellnr](l2,l1) = cov_vsrho->getRealValueCyclic(i1-i2,j1-j2,k1-k2);
    }
  }
}


void  SpatialSyntWellFilter::DoFilteringSyntWells(SeismicParametersHolder                  & seismicParameters,
                                                  const NRLib::Matrix                      & priorVar0)
{

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

  if(sigmae_.size() == 0) {
    sigmae_.resize(nDim);
    sigmae_[0].resize(3,3);
  }

  NRLib::Matrix priorCov0 = priorVar0;

  bool no_wells_filtered = true;

  for(int w1=0;w1<nWellsToBeFiltered_;w1++){
    //LogKit::LogFormatted(LogKit::Low,"\nFiltering synthetic well number " + NRLib::ToString(w1+1,1) + "...");
    no_wells_filtered = false;

    int n = syntWellData_[w1]->getWellLength();

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
    const int *ipos = syntWellData_[w1]->getIpos();
    const int *jpos = syntWellData_[w1]->getJpos();
    const int *kpos = syntWellData_[w1]->getKpos();
    float regularization = Definitions::SpatialFilterRegularisationValue();

    FillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCovVp(),      n, 0,   0);
    FillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCovVs(),      n, n,   n);
    FillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCovRho(),     n, 2*n, 2*n);
    FillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCrCovVpVs(),  n, 0,   n);
    FillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCrCovVpRho(), n, 0,   2*n);
    FillValuesInSigmapostSyntWell(sigmapost, ipos, jpos, kpos, seismicParameters.GetCrCovVsRho(), n, 2*n, n);

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

        //sigmapost[l2 + n  ][l1      ] = sigmapost[l1][n+l2];
        //sigmapost[l2 + 2*n][l1      ] = sigmapost[l1][2*n+l2];
        //sigmapost[l2 + n  ][l1 + 2*n] = sigmapost[2*n+l1][n+l2];
        sigmapri [l1      ][l2      ] = prior_cov_vp_[w1](l1,l2);//priorCov0(0,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + n  ][l2 + n  ] = prior_cov_vs_[w1](l1,l2);//priorCov0(1,1)*priorSpatialCorr_[w1][l1][l2];
        sigmapri [l1 + 2*n][l2 + 2*n] = prior_cov_rho_[w1](l1,l2);//priorCov0(2,2)*priorSpatialCorr_[w1][l1][l2];
        if(l1==l2){
          sigmapost[l1      ][l2      ] += regularization*sigmapost[l1][l2]/sigmapri[l1][l2];
          sigmapost[l1 + n  ][l2 + n  ] += regularization*sigmapost[n+l1][n+l2]/sigmapri[n+l1][n+l2];
          sigmapost[l1 + 2*n][l2 + 2*n] += regularization*sigmapost[2*n+l1][2*n+l2]/sigmapri[2*n+l1][2*n+l2];
          sigmapri [l1      ][l2      ] += regularization;
          sigmapri [l1 + n  ][l2 + n  ] += regularization;
          sigmapri [l1 + 2*n][l2 + 2*n] += regularization;
        }
        sigmapri[l1      ][l2 + n  ] = prior_cov_vpvs_[w1](l1,l2);
        sigmapri[l2      ][l1 + n  ] = prior_cov_vpvs_[w1](l2,l1);
        sigmapri[l1][l2+2*n]         = prior_cov_vprho_[w1](l1,l2);
        sigmapri[l2][l1+2*n]         = prior_cov_vprho_[w1](l2,l1);
        sigmapri[l1 + n  ][l2 + 2*n] = prior_cov_vsrho_[w1](l1,l2);
        sigmapri[l2 + n  ][l1 + 2*n] = prior_cov_vsrho_[w1](l2,l1);
        /*
        sigmapri[l1 + n  ][l2      ] = priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2      ][l1 + n  ] = priorCov0(1,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l1 + 2*n][l2      ] = priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2      ][l1 + 2*n] = priorCov0(2,0)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l1 + n  ][l2 + 2*n] = priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
        sigmapri[l2 + 2*n][l1 + n  ] = priorCov0(2,1)*priorSpatialCorr_[w1][l1][l2];
        */
      }
    }

    NRLib::SymmetricMatrix Sprior(3*n);
    NRLib::SymmetricMatrix Spost(3*n);

    for (int i = 0; i < 3*n; i++)
      for (int j = 0; j <= i; j++)
        Sprior(j, i) = sigmapri[j][i];

    for (int i = 0; i < 3*n; i++)
      for (int j = i; j <= i; j++)
        Spost(j, i) = sigmapost[j][i];

    //
    // Filter = I - Sigma_post * inv(Sigma_prior)
    //
    NRLib::Matrix Aw;
    NRLib::Matrix I = NRLib::IdentityMatrix(3*n);
    NRLib::CholeskySolve(Sprior, I);


    Aw = Spost * I;
    //NRLib::WriteMatrixToFile("SpostxSpriinv_2014.dat", Aw);
    Aw = Aw * (-1);
    for(int i=0 ; i<3*n ; i++) {
      Aw(i,i) += 1.0;
    }

    updateSigmaE(sigmae_[0],
                 Aw,
                 Spost,
                 n);


    lastn += n;

    for(int i=0;i<3*n;i++){
      //delete [] Aw[i];
      delete [] sigmapri[i];
      delete [] sigmapost[i];
      delete [] imat[i];
    }
    //delete [] Aw;
    delete [] sigmapri;
    delete [] sigmapost;
    delete [] imat;
  }

  if(no_wells_filtered == false){
    // finds the scale at default inversion (all minimum noise in case of local noise)
    NRLib::Matrix Se(3,3);
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

void SpatialSyntWellFilter::FillValuesInSigmapostSyntWell(double      ** sigmapost,
                                                          const int    * ipos,
                                                          const int    * jpos,
                                                          const int    * kpos,
                                                          FFTGrid      * covgrid,
                                                          int            n,
                                                          int            ni,
                                                          int            nj)
{
  //double minValue = std::pow(10.0,-9);
  //int nz = covgrid->getNz();
  //double factorNorm = 1/sqrt(2*NRLib::Pi);
  //double endValue = covgrid->getRealValueCyclic(0,0,nz-1);
  //double div = 1/(nz*.1);
  // tapering limit at 90%
  //double smoothLimit = nz*.9;

  covgrid->setAccessMode(FFTGrid::RANDOMACCESS);
  int i1, j1, k1, l1, i2, j2, k2, l2;
  for(l1 = 0; l1<n; l1++)
  {
    i1 = ipos[l1];
    j1 = jpos[l1];
    k1 = kpos[l1];
    for (l2 = l1 ; l2 < n ; l2++) {
      i2 = ipos[l2];
      j2 = jpos[l2];
      k2 = kpos[l2];

      sigmapost[l1+ni][l2+nj] = covgrid->getRealValueCyclic(i1-i2, j1-j2, k2-k1);
      sigmapost[l2+ni][l1+nj] = covgrid->getRealValueCyclic(i1-i2, j1-j2, k1-k2);
    }
  }
  covgrid->endAccess();

}

/*
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
    NRLib::Matrix Se(3,3);
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
*/

void    SpatialSyntWellFilter::GenerateSyntWellData (const std::map<std::string, DistributionsRock *>        & rock_distributions,
                                                     const std::vector<std::string>                          & facies_names,
                                                     double                                                    dz,
                                                     int                                                       syntWellLength)
{
  int nFacies = static_cast<int>(facies_names.size());

  // 10.0 is the desired mean length
  double p = 0.0;
  // lambda is set below, initialized to 1.0
  double lambda = 1.0;
  // corr is initially set to 0.5
  double corr = 0.5;
  if (dz < 10.0) { //Erik: dz<log(10)
    // this ensures that the mean of the geometric distribution is 10.0/dz
    p = dz / 10.0; // Erik: exp(dz)/10.0;
    // calculate lambda for the exponential distribution used below
    lambda = - std::log((1.0 - p));
  }
  else
    throw NRLib::Exception("Facies probabilities: dz is too large to generate synthetic well data, need dz < 10");

  int nWell = 0;
  for (int i = 0; i < static_cast<int>(trend_1_.size()); i++) {
    for (int j = 0; j < static_cast<int>(trend_2_.size()); j++) {
      for (int m = 0; m < nWellsPerCombinationOfTrendParams_; m++) {

        std::vector<double> trend_params(2);

        trend_params[0] = trend_1_[i];
        trend_params[1] = trend_2_[j];

        std::vector<float> vp;
        std::vector<float> vs;
        std::vector<float> rho;
        std::vector<int>   facies;

        int k = 0;

        while (k < syntWellLength ) {
          // pick a random (uniform) facies
          int randomFacies = static_cast<int>(floor(nFacies * NRLib::Random::Unif01()));

          // draw the facies length from a geometric distribution
          double uRan            = NRLib::Random::Unif01();
          int randomFaciesLength = static_cast<int>(ceil(-std::log(uRan)/lambda));

          std::vector<double> vp_temp(randomFaciesLength);
          std::vector<double> vs_temp(randomFaciesLength);
          std::vector<double> rho_temp(randomFaciesLength);

          std::map<std::string, DistributionsRock *>::const_iterator iter = rock_distributions.find(facies_names[randomFacies]);

          iter->second->GenerateWellSample(corr, vp_temp, vs_temp, rho_temp, trend_params);

          for (int l = 0; l < randomFaciesLength; l++){
            facies.push_back(randomFacies);

            vp.push_back(static_cast<float>(log(vp_temp[l])));
            vs.push_back(static_cast<float>(log(vs_temp[l])));
            rho.push_back(static_cast<float>(log(rho_temp[l])));
          }

          k += randomFaciesLength;
        }
        syntWellData_[nWell] = new SyntWellData(trend_1_[i], trend_2_[j], i, j, vp, vs, rho, facies, facies_names);
        nWell++;
      }
    }
  }

}
