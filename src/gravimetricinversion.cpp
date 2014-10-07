/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

//#include "rfftw.h"

#include "src/gravimetricinversion.h"
#include "src/modelgeneral.h"
#include "src/modelgravitystatic.h"
#include "src/modelgravitydynamic.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/simbox.h"
#include "src/definitions.h"
#include "src/io.h"
#include "src/parameteroutput.h"

#include "lib/timekit.hpp"
//#include "lib/lib_matr.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>
#include <string>
#include <algorithm>

GravimetricInversion::GravimetricInversion()
{
}

/*
GravimetricInversion::GravimetricInversion(ModelGeneral            *  modelGeneral,
                                           ModelGravityStatic      *  modelGravityStatic,
                                           ModelGravityDynamic     *& modelGravityDynamic,
                                           SeismicParametersHolder &  seismicParameters)
{
  LogKit::WriteHeader("Building Stochastic Model for Gravimetric Inversion");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  State4D * state4d = modelGeneral->getState4D();
  lag_index_              = modelGravityStatic->GetLagIndex();

  int nxp               = state4d->getMuVpStatic()->getNxp();
  int nyp               = state4d->getMuVpStatic()->getNyp();
  int nzp               = state4d->getMuVpStatic()->getNzp();

  int nx_upscaled       = modelGravityStatic->GetNx_upscaled();
  int ny_upscaled       = modelGravityStatic->GetNy_upscaled();
  int nz_upscaled       = modelGravityStatic->GetNz_upscaled();

  int nxp_upscaled      = modelGravityStatic->GetNxp_upscaled();
  int nyp_upscaled      = modelGravityStatic->GetNyp_upscaled();
  int nzp_upscaled      = modelGravityStatic->GetNzp_upscaled();
  int Np_up             = nxp_upscaled*nyp_upscaled*nzp_upscaled;

  bool   include_level_shift = true;
  double shift_parameter    = 0;
  double level_shift        = 0;


  FFTGrid * upscaling_kernel_conj = new FFTGrid(modelGravityStatic->GetUpscalingKernel());
  if(upscaling_kernel_conj->getIsTransformed() == false){
    upscaling_kernel_conj->fftInPlace();
  }
  upscaling_kernel_conj->conjugate();  // Conjugate only in FFT domain.

  FFTGrid * upscaling_kernel_abs = new FFTGrid(modelGravityStatic->GetUpscalingKernel());
  upscaling_kernel_abs->realAbs();

  // Find distribution of rho_current
  FFTGrid * mean_rho_current  = new FFTGrid(seismicParameters.GetMuRho());    // At this stage we are on log scale, \mu_{log rho^c}
  FFTGrid * cov_rho_current   = new FFTGrid(seismicParameters.GetCovRho());   // at this stage: on log scale

  float sigma_squared = GetSigmaForTransformation(cov_rho_current);
  float mean_current  = GetMeanForTransformation (mean_rho_current);

  //Need to be in real domain from transforming from log domain
  if(mean_rho_current->getIsTransformed())
    mean_rho_current->invFFTInPlace();
  if(cov_rho_current->getIsTransformed())
    cov_rho_current->invFFTInPlace();

  MeanExpTransform(mean_rho_current, sigma_squared);
  CovExpTransform (cov_rho_current,  mean_current);

  // Compute joint distributions of [mc, ms] = exp([rho_current_log, rho_static_log])
  FFTGrid * mean_rho_static = new FFTGrid(state4d->getMuRhoStatic());
  FFTGrid * cov_rho_static  = new FFTGrid(state4d->getCovRhoRhoStaticStatic());

  sigma_squared      = GetSigmaForTransformation(cov_rho_static);
  float mean_static  = GetMeanForTransformation(mean_rho_static);

  //Need to be in real domain from transforming from log domain
  if(mean_rho_static->getIsTransformed())
    mean_rho_static->invFFTInPlace();

  if(cov_rho_static->getIsTransformed())
    cov_rho_static->invFFTInPlace();

  MeanExpTransform(mean_rho_static, sigma_squared);
  CovExpTransform (cov_rho_static,  mean_static);

  // Cov of exp(rho_static_current_log) cov current = cov static static + cov static dynamic
  FFTGrid * cov_rhorho_static_current = new FFTGrid(state4d->getCovRhoRhoStaticStatic());
  FFTGrid * cov_rhorho_static_dynamic = new FFTGrid(state4d->getCovRhoRhoStaticDynamic());

  if(cov_rhorho_static_current->getIsTransformed())
    cov_rhorho_static_current->invFFTInPlace();

  if(cov_rhorho_static_dynamic->getIsTransformed())
    cov_rhorho_static_dynamic->invFFTInPlace();

  cov_rhorho_static_current->add(cov_rhorho_static_dynamic);
  CrCovExpTransform(cov_rhorho_static_current, mean_static, mean_current);

  // Mean of new parameter
  if(mean_rho_current->getIsTransformed())
    mean_rho_current->invFFTInPlace();

  if(mean_rho_static->getIsTransformed())
    mean_rho_static->invFFTInPlace();

  FFTGrid * mean_rho_total = new FFTGrid(mean_rho_current);
  mean_rho_total->subtract(mean_rho_static);

  // Covariance of new parameter
  FFTGrid * cov_rho_total = new FFTGrid(cov_rho_static);
  cov_rho_total->add(cov_rho_current);
  cov_rho_total->subtract(cov_rhorho_static_current);

  cov_rho_total            ->setAccessMode(FFTGrid::READANDWRITE);
  cov_rhorho_static_current->setAccessMode(FFTGrid::READ);

  float value, value2;
  for(int i = 0; i < nxp; i++){
    for(int j = 0; j < nyp; j++){
      for(int k = 0; k < nzp; k++){
        value  = cov_rho_total->getRealValue(i,j,k);
        value2 = cov_rhorho_static_current->getRealValueCyclic(-i, -j, -k);
        cov_rho_total->setRealValue(i, j, k, value-value2);
      }
    }
  }
  cov_rho_total           ->endAccess();
  cov_rhorho_static_current->endAccess();

// Now we have a parameter for inversion: meanRhoTotal and covRhoTotal
  if(mean_rho_total->getIsTransformed() == false)
    mean_rho_total->fftInPlace();

  // Upscale Rho: Convolution in FFT domain and pointwise multiplication
  mean_rho_total->multiply(upscaling_kernel_conj);  // Now is expMeanRhoTotal smoothed

  // Subsample in FFTDomain
  FFTGrid * upscaled_mean_rho_total;
  Subsample(upscaled_mean_rho_total, mean_rho_total,
            nx_upscaled,  ny_upscaled,  nz_upscaled,
            nxp_upscaled, nyp_upscaled, nzp_upscaled);

  // Upscale Covariance
  if(cov_rho_total->getIsTransformed() == false)
    cov_rho_total->fftInPlace();

  // Convolution in the FFTdomain;
  cov_rho_total->multiply(upscaling_kernel_abs);  // Now is expCovRhoTotal smoothed
  cov_rho_total->multiply(upscaling_kernel_abs);  // Abs value (or use complex conjugate);

  // Subsample in FFTDomain
  FFTGrid * upscaled_cov_rho_total;
  Subsample(upscaled_cov_rho_total, cov_rho_total,
            nx_upscaled, ny_upscaled, nz_upscaled,
            nxp_upscaled, nyp_upscaled, nzp_upscaled);


  if(upscaled_mean_rho_total->getIsTransformed())
    upscaled_mean_rho_total->invFFTInPlace();

  if(upscaled_cov_rho_total->getIsTransformed())
    upscaled_cov_rho_total->invFFTInPlace();


  LogKit::WriteHeader("Performing Gravimetric Inversion");

  NRLib::Matrix G = modelGravityDynamic->GetGMatrix();
  ExpandMatrixWithZeros(G, Np_up, include_level_shift);

  NRLib::Vector    Rho(Np_up);
  VectorizeFFTGrid(Rho, upscaled_mean_rho_total);

  NRLib::Matrix            Sigma(Np_up, Np_up);
  NRLib::InitializeMatrix (Sigma, 0.0);
  ReshapeCovAccordingToLag(Sigma, upscaled_cov_rho_total);

  NRLib::Vector      gravity_data(30);
  std::vector<float> d = modelGravityDynamic->GetGravityResponse();

  NRLib::Matrix           Sigma_error(30, 30);
  NRLib::InitializeMatrix(Sigma_error, 0.0);
  std::vector<float> std_dev = modelGravityDynamic->GetGravityStdDev();

  //spool over data from std vec to NRLib Vector and summing the squares of the data values for use in level shift
    for(int i = 0; i<gravity_data.length(); i++){
      gravity_data(i)  = d[i];
      shift_parameter += d[i]*d[i];
      Sigma_error(i,i) = std_dev[i];
    }
    shift_parameter = shift_parameter*10;

  if(include_level_shift){
    // Expand prior mean with one element equal to 0
    int l = Rho.length();
    NRLib::Vector RhoNew(l+1);
    for(int i = 0; i<Rho.length(); i++){
      RhoNew(i) = Rho(i);
    }
    RhoNew(l) = 0;   // set last value
    Rho = RhoNew;

    // Expand prior covariance matrix
    ExpandCovMatrixWithLevelShift(Sigma, shift_parameter);
  }

  NRLib::WriteVectorToFile("Rho_prior.txt", Rho);

  NRLib::Vector Rho_posterior  (Np_up);
  NRLib::Matrix Sigma_posterior(Np_up, Np_up);

  NRLib::Matrix GT         = NRLib::transpose(G);
  NRLib::Matrix G_Sigma    = G * Sigma;
  NRLib::Matrix G_Sigma_GT = G_Sigma * GT;
  NRLib::Matrix Sigma_GT   = Sigma * GT;

  NRLib::Matrix inv_G_Sigma_GT_plus_Sigma_error = G_Sigma_GT + Sigma_error;
  NRLib::Invert(inv_G_Sigma_GT_plus_Sigma_error);

  NRLib::Vector temp_1 = gravity_data - G*Rho;
  NRLib::Vector temp_2 = inv_G_Sigma_GT_plus_Sigma_error * temp_1;
  temp_1               = Sigma_GT*temp_2;
  Rho_posterior        = Rho + temp_1;


  NRLib::Matrix temp_3 = inv_G_Sigma_GT_plus_Sigma_error * G_Sigma;
  NRLib::Matrix temp_4 = Sigma_GT*temp_3;
  Sigma_posterior      = Sigma - temp_4;

  // Remove shift parameter
  if(include_level_shift){
    RemoveLevelShiftFromVector(Rho_posterior, level_shift);
    RemoveLevelShiftFromCovMatrix(Sigma_posterior);
  }
  NRLib::WriteVectorToFile("Rho_posterior.txt", Rho_posterior);

  // Reshape back to FFTGrid
  ReshapeVectorToFFTGrid(upscaled_mean_rho_total, Rho_posterior);

  // For backsampling and deconvolution, need to be in FFT-domain
  if(mean_rho_total->getIsTransformed()==false)
    mean_rho_total->fftInPlace();

  if(upscaled_mean_rho_total->getIsTransformed()==false)
    upscaled_mean_rho_total->fftInPlace();

  //Backsample and deconvolve in Fourier domain
  Backsample(upscaled_mean_rho_total, mean_rho_total); //Now meanRhoTotal is posterior!
  Divide(mean_rho_total, upscaling_kernel_conj);

  /// Posterior upscaled covariance
  FFTGrid * posterior_upscaled_cov_rho_total = new FFTGrid(nx_upscaled, ny_upscaled, nz_upscaled,
                                                           nxp_upscaled, nyp_upscaled, nzp_upscaled);
  posterior_upscaled_cov_rho_total->createRealGrid();
  posterior_upscaled_cov_rho_total->setType(FFTGrid::PARAMETER);

  ReshapeCovMatrixToFFTGrid(posterior_upscaled_cov_rho_total, Sigma_posterior);


  // Odds algorithme
  // Pick elements that corresponds to effect of inversion, not from adjusting to pos def matrix, in FFT domain
  FFTGrid * fft_factor = new FFTGrid(nx_upscaled,  ny_upscaled,  nz_upscaled,
                                     nxp_upscaled, nyp_upscaled, nzp_upscaled);
  fft_factor->createRealGrid();
  fft_factor->setType(FFTGrid::PARAMETER);
  fft_factor->fillInConstant(1.0);  // ikke i padded område

  // Need to be in FFT-domain
  if(upscaled_cov_rho_total->getIsTransformed() == false)
    upscaled_cov_rho_total->fftInPlace();
  if(posterior_upscaled_cov_rho_total->getIsTransformed() == false)
    posterior_upscaled_cov_rho_total->fftInPlace();

  fftw_complex reference;
  double nu = 0.05;
  fftw_complex prior = upscaled_cov_rho_total          ->getFirstComplexValue();
  fftw_complex post  = posterior_upscaled_cov_rho_total->getFirstComplexValue();
  reference.re = prior.re - post.re;
  reference.im = prior.im - post.im;
  reference.re = static_cast<fftw_real>(reference.re*nu);
  reference.im = static_cast<fftw_real>(reference.im*nu);

  fft_factor                      ->setAccessMode(FFTGrid::WRITE);
  upscaled_cov_rho_total          ->setAccessMode(FFTGrid::READ);
  posterior_upscaled_cov_rho_total->setAccessMode(FFTGrid::READ);

  for(int k = 0; k < nz_upscaled; k++){
    for(int j = 0; j < ny_upscaled; j++){
      for(int i = 0; i < nx_upscaled; i++){
        prior = upscaled_cov_rho_total          ->getComplexValue(i,j,k);
        post  = posterior_upscaled_cov_rho_total->getComplexValue(i,j,k);

        fftw_complex m;
        m.re = prior.re - post.re;
        m.im = prior.im - post.im;

        if(m.re > reference.re){  // Need only consider real part
          fftw_complex ratio;
          ratio.re = posterior_upscaled_cov_rho_total->getComplexValue(i,j,k).re/upscaled_cov_rho_total->getComplexValue(i,j,k).re;

          float value = static_cast<float>(ratio.re);
          if(value < 1)
            fft_factor->setRealValue(i, j, k, value);
          else
            fft_factor->setRealValue(i, j, k, 1);
        }
        else{
          fft_factor->setRealValue(i, j, k, 1);
        }
      }
    }
  }
  fft_factor                      ->endAccess();
  upscaled_cov_rho_total          ->endAccess();
  posterior_upscaled_cov_rho_total->endAccess();

  posterior_upscaled_cov_rho_total = upscaled_cov_rho_total; // Set posterior equal to prior and introduce effects of inversion through factors in fft_factor-grid

  posterior_upscaled_cov_rho_total->multiply(fft_factor);  // pointwise multiplication in real domain

  // For backsampling and deconvolution need to be in Fourier domain
   if(cov_rho_total->getIsTransformed() == false)
    cov_rho_total->fftInPlace();

   if(posterior_upscaled_cov_rho_total->getIsTransformed() == false)
     posterior_upscaled_cov_rho_total->fftInPlace();

  //Backsample in Fourier domain
  Backsample(posterior_upscaled_cov_rho_total, cov_rho_total);  // covRhoTotal is now posterior!

  // Deconvolution: Do pointwise division in FFTdomain - twice
  Divide(cov_rho_total, upscaling_kernel_abs);
  Divide(cov_rho_total, upscaling_kernel_abs);


   // Only for debugging purposes
  if(posterior_upscaled_cov_rho_total->getIsTransformed() == true)
    posterior_upscaled_cov_rho_total->invFFTInPlace();
  NRLib::Matrix Post_sigma_temp(Np_up, Np_up);
  NRLib::InitializeMatrix (Post_sigma_temp, 0.0);

  ReshapeCovAccordingToLag(Post_sigma_temp, posterior_upscaled_cov_rho_total);
  NRLib::WriteMatrixToFile("Sigma_posterior.txt", Post_sigma_temp);


  // For transforming back to log-domain, need to be in real domain
  if(cov_rho_total->getIsTransformed())
    cov_rho_total->invFFTInPlace();

  if(mean_rho_total->getIsTransformed())
    mean_rho_total->invFFTInPlace();

  // Transform back to log scale
  float mean_total = GetMeanForTransformation(mean_rho_total);
  CovLogTransform(cov_rho_total, mean_total);

  sigma_squared = GetSigmaForTransformation(cov_rho_total);
  MeanLogTransform(mean_rho_total, sigma_squared);

  // Kind of ready to put back into seismicparametersholder and state4d.
  // Need some merge functions or split or something
  // Sorry about uncomplete code.

  ComputeSyntheticGravimetry(mean_rho_total, modelGravityDynamic, level_shift);

  // Delete all FFTGrids create with "new"
  delete upscaling_kernel_conj;
  delete upscaling_kernel_abs;
  delete mean_rho_static;
  delete cov_rho_static;
  delete cov_rhorho_static_current;
  delete cov_rhorho_static_dynamic;
  delete mean_rho_total;
  delete cov_rho_total;

  delete fft_factor;
}

GravimetricInversion::~GravimetricInversion()
{
}

float
  GravimetricInversion::GetSigmaForTransformation(FFTGrid * sigma)
{
  float sigma_squared = 0;
  if(sigma->getIsTransformed() == false){
    sigma_squared = sigma->getFirstRealValue();   // \sigma^2
  }
  else{
  // Loop though grid in complex domain.
    sigma_squared = 0;

    fftw_complex sum;
    sum.re = 0.0;
    sum.im = 0.0;
    int f = 1;

    sigma->setAccessMode(FFTGrid::READ);

    for(int k=0; k<sigma->getNzp(); k++){
      for(int j=0; j<sigma->getNyp(); j++){
        for(int i=0; i<sigma->getCNxp(); i++){
          f=2;
          if(i == 0)
            f = 1;
          if(sigma->getCNxp() % 2 == 1){
            if(i==sigma->getCNxp()-1)
              f=1;
          }

          fftw_complex value = sigma->getNextComplex();
          sum.re += f*value.re;
          sum.im += f*value.im;   // Blir ikke null, fordi summerer ikke konjugerte par i praksis
        }
      }
    }
    sigma->endAccess();

    // Due to summing complex conjugate numbers, then imaginary part is zero
    int N = sigma->getNxp()*sigma->getNyp()*sigma->getNzp();
    sigma_squared = sum.re/N;
  }
  return(sigma_squared);
}

float
  GravimetricInversion::GetMeanForTransformation(FFTGrid * grid)
{
  float mean = 0;

  if(grid->getIsTransformed() == true){
    int nxp = grid->getNxp();
    int nyp = grid->getNyp();
    int nzp = grid->getNzp();

    fftw_complex mean_tmp = grid->getFirstComplexValue();          // Standard way of getting mean value of a FFTGrid (including padded region)
    mean                  = mean_tmp.re/pow(static_cast<float>(nxp*nyp*nzp), 0.5f);   // Note the scaling
  }
  else{
    // Loop through grid in real domain
    grid->setAccessMode(FFTGrid::READ);
    float sum = 0.0;
    for(int k=0;k<grid->getNzp();k++) {
      for(int j=0;j<grid->getNyp();j++) {
        for(int i=0;i<grid->getNxp();i++) {
          sum += grid->getNextReal();
        }
        for(int i = 0; i<grid->getRNxp()-grid->getNxp(); i++){
          grid->getNextReal();
        }
      }
    }
    grid->endAccess();

    float N = static_cast<float>(grid->getNxp()*grid->getNyp()*grid->getNzp());
    mean = sum/N;

  }

  return(mean);
}

void
  GravimetricInversion::MeanExpTransform(FFTGrid * log_mean, float sigma_squared)
{
  assert(log_mean->getIsTransformed() == false);

  log_mean->addScalar(0.5f*sigma_squared);  // \mu_{log rho^c} + 0.5*\sigma^2
  log_mean->expTransf();  //exp{\mu_{log rho^c} + 0.5*\sigma^2}. Finished transformation.
}

void
  GravimetricInversion::CovExpTransform(FFTGrid  * log_cov, float mean)
{
  assert(log_cov->getIsTransformed() == false);

  log_cov->expTransf();
  log_cov->addScalar(-1);
  log_cov->multiplyByScalar(mean*mean);
}

void
  GravimetricInversion::CrCovExpTransform(FFTGrid * log_cov,  float mean_a, float mean_b)
{
  assert(log_cov->getIsTransformed() == false);

  log_cov->expTransf();
  log_cov->addScalar(-1);
  log_cov->multiplyByScalar(mean_a*mean_b);
}

void
  GravimetricInversion::MeanLogTransform(FFTGrid * mean,    float sigma_squared)
{
  assert(mean->getIsTransformed() == false);

  mean->logTransf();
  mean->addScalar(-0.5f*sigma_squared);
}

void
  GravimetricInversion::CovLogTransform(FFTGrid  * cov,     float mean)
{
  assert(cov->getIsTransformed() == false);

  cov->multiplyByScalar(1.0f/(mean*mean));
  cov->addScalar(1);
  cov->logTransf();
}

void
  GravimetricInversion::Subsample(FFTGrid *& upscaled_grid, FFTGrid * original_grid, int nx_up, int ny_up, int nz_up, int nxp_up, int nyp_up, int nzp_up)
{
  assert(original_grid->getIsTransformed());

  upscaled_grid = new FFTGrid(nx_up, ny_up, nz_up, nxp_up, nyp_up, nzp_up);
  upscaled_grid->createComplexGrid();
  upscaled_grid->setType(FFTGrid::PARAMETER);

  upscaled_grid->setAccessMode(FFTGrid::WRITE);
  original_grid->setAccessMode(FFTGrid::READ);

  int nxp_up_half = static_cast<int>(ceil(static_cast<double>(nxp_up)/2));
  int nyp_up_half = static_cast<int>(ceil(static_cast<double>(nyp_up)/2));
  int nzp_up_half = static_cast<int>(ceil(static_cast<double>(nzp_up)/2));
  int nxp = original_grid->getNxp();
  int nyp = original_grid->getNyp();
  int nzp = original_grid->getNzp();

  // Set up index-vectorer
  std::vector<int> x_indices(nxp_up);
  std::vector<int> y_indices(nyp_up);
  std::vector<int> z_indices(nzp_up);

  for(int i = 0; i<nxp_up; i++){
    if(i<nxp_up_half)
      x_indices[i] = i;
    else
      x_indices[i] = nxp - (nxp_up - i);
  }
  for(int i = 0; i<nyp_up; i++){
    if(i<nyp_up_half)
      y_indices[i] = i;
    else
      y_indices[i] = nyp - (nyp_up - i);
  }
  for(int i = 0; i<nzp_up; i++){
    if(i<nzp_up_half)
      z_indices[i] = i;
    else
      z_indices[i] = nzp - (nzp_up - i);
  }

  // Subsample the grid in FFT domain - use also padded region.
  for(int i = 0; i < nxp_up; i++){
    for(int j = 0; j < nyp_up; j++){
      for(int k = 0; k < nzp_up; k++){
        upscaled_grid->setComplexValue(i, j, k, original_grid->getComplexValue(x_indices[i], y_indices[j], z_indices[k], true));
      }
    }
  }
  upscaled_grid->endAccess();
  original_grid->endAccess();
}

void
  GravimetricInversion::Backsample(FFTGrid * upscaled_grid, FFTGrid * new_full_grid)
{
  assert(upscaled_grid->getIsTransformed());
  assert(new_full_grid->getIsTransformed());

  new_full_grid->setAccessMode(FFTGrid::WRITE);
  upscaled_grid->setAccessMode(FFTGrid::READ);

  int nxp_up = upscaled_grid->getNxp();
  int nyp_up = upscaled_grid->getNyp();
  int nzp_up = upscaled_grid->getNzp();

  int nxp = new_full_grid->getNxp();
  int nyp  = new_full_grid->getNyp();
  int nzp  = new_full_grid->getNzp();

  int nxp_up_half = static_cast<int>(ceil(static_cast<double>(nxp_up)/2));
  int nyp_up_half = static_cast<int>(ceil(static_cast<double>(nyp_up)/2));
  int nzp_up_half = static_cast<int>(ceil(static_cast<double>(nzp_up)/2));

  // Set up index-vectorer
  std::vector<int> x_indices(nxp_up);
  std::vector<int> y_indices(nyp_up);
  std::vector<int> z_indices(nzp_up);

  for(int i = 0; i<nxp_up; i++){
    if(i<nxp_up_half)
      x_indices[i] = i;
    else
      x_indices[i] = nxp - (nxp_up - i);
  }
  for(int i = 0; i<nyp_up; i++){
    if(i<nyp_up_half)
      y_indices[i] = i;
    else
      y_indices[i] = nyp - (nyp_up - i);
  }
  for(int i = 0; i<nzp_up; i++){
    if(i<nzp_up_half)
      z_indices[i] = i;
    else
      z_indices[i] = nzp - (nzp_up - i);
  }

  // Subsample the grid in FFT domain - use also padded region
  for(int i = 0; i < nxp_up; i++){
    for(int j = 0; j < nyp_up; j++){
      for(int k = 0; k < nzp_up; k++){
        new_full_grid->setComplexValue(x_indices[i], y_indices[j], z_indices[k], upscaled_grid->getComplexValue(i,j,k, true));
      }
    }
  }

  new_full_grid->endAccess();
  upscaled_grid->endAccess();
  }

void
  GravimetricInversion::VectorizeFFTGrid(NRLib::Vector &vec, FFTGrid * grid, bool withPadding)
{
  assert(grid->getIsTransformed() == false);

  int nx, ny, nz;
  if(withPadding){
  nx = grid->getNxp();
  ny = grid->getNyp();
  nz = grid->getNzp();
  }
  else{
  nx = grid->getNx();
  ny = grid->getNy();
  nz = grid->getNz();
  }

  grid->setAccessMode(FFTGrid::READ);

  int I = 0;
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      for(int k = 1; k <= nz; k++){
        I =  i + (j-1)*nx + (k-1)*nx*ny;
        vec(I-1) = grid->getNextReal();
      }
    }
  }
  grid->endAccess();
}

void
  GravimetricInversion::ReshapeVectorToFFTGrid(FFTGrid * grid, NRLib::Vector vec)
{
  assert(grid->getIsTransformed() == false);

  grid->setAccessMode(FFTGrid::WRITE);

  int nxp = grid->getNxp();
  int nyp = grid->getNyp();
  int nzp = grid->getNzp();
  int I = 0;
  for(int i = 1; i <= nxp; i++){
    for(int j = 1; j <= nyp; j++){
      for(int k = 1; k <= nzp; k++){
        I =  i + (j-1)*nxp + (k-1)*nxp*nyp;
        grid->setNextReal(static_cast<float>(vec(I-1)));
      }
    }
  }
  grid->endAccess();
}

void
  GravimetricInversion::ReshapeCovAccordingToLag(NRLib::Matrix &CovMatrix, FFTGrid * covGrid)
{
  // Reshape according to lag
  assert(covGrid->getIsTransformed() == false);

  int nxp = covGrid->getNxp();
  int nyp = covGrid->getNyp();
  int nzp = covGrid->getNzp();

  covGrid->setAccessMode(FFTGrid::READ);
  int I;
  int J;
  for(int k1 = 1; k1 <= nzp; k1++){
    for(int j1 = 1; j1 <= nyp; j1++){
      for(int i1 = 1; i1 <= nxp; i1++){
        I =  i1 + (j1-1)*nxp + (k1-1)*nxp*nyp;

        for(int k2 = 1; k2 <= nzp; k2++){
          for(int j2 = 1; j2 <= nyp; j2++){
            for(int i2 = 1; i2 <= nxp; i2++){
              J = i2 + (j2-1)*nxp + (k2-1)*nxp*nyp;

              if(lag_index_[I-1][J-1][0]==-1 && lag_index_[I-1][J-1][1]==-1 && lag_index_[I-1][J-1][2]==-1)
                CovMatrix(I-1,J-1) = 0.0;
              else{
                CovMatrix(I-1,J-1) = covGrid->getRealValue(lag_index_[I-1][J-1][0], lag_index_[I-1][J-1][1], lag_index_[I-1][J-1][2], true);
              }
            }
          }
        }

      }
    }
  }
  covGrid->endAccess();
}

void
  GravimetricInversion::ReshapeCovMatrixToFFTGrid(FFTGrid * cov_grid, NRLib::Matrix cov_matrix)
{
  // Reshape back from covariance matrix to 3D cube

  assert(cov_grid->getIsTransformed() == false);

  int nx = cov_grid->getNx();
  int ny = cov_grid->getNy();
  int nz = cov_grid->getNz();

  int nxp = cov_grid->getNxp();
  int nyp = cov_grid->getNyp();
  int nzp = cov_grid->getNzp();

  // initilize arrays
  std::vector<std::vector<std::vector<float> > > sum;
  std::vector<std::vector<std::vector<int> > >   counter;
  counter.resize(nxp);
  sum    .resize(nxp);
  for (int i = 0; i < nxp; ++i) {
    counter[i].resize(nyp);
    sum[i]    .resize(nyp);
    for (int j = 0; j < nyp; ++j){
      counter[i][j].resize(nzp);
      sum[i][j]    .resize(nzp);
    }
  }

  cov_grid->setAccessMode(FFTGrid::WRITE);
  // Intended: One set of for loops over padded region as well, the other set of for loops over nx, ny, nz
  int I, J;
  int i, j, k;
  for(int k1 = 1; k1 <= nzp; k1++){
    for(int j1 = 1; j1 <= nyp; j1++){
      for(int i1 = 1; i1 <= nxp; i1++){
        I =  i1 + (j1-1)*nxp + (k1-1)*nxp*nyp;
        for(int k2 = 1; k2 <= nz; k2++){
          for(int j2 = 1; j2 <= ny; j2++){
            for(int i2 = 1; i2 <= nx; i2++){
              J = i2 + (j2-1)*nx + (k2-1)*nx*ny;
              if(lag_index_[I-1][J-1][0] >= 0 && lag_index_[I-1][J-1][1] >= 0 && lag_index_[I-1][J-1][2] >= 0){
                i = lag_index_[I-1][J-1][0];
                j = lag_index_[I-1][J-1][1];
                k = lag_index_[I-1][J-1][2];
                sum[i][j][k]    += static_cast<float>(cov_matrix(I-1,J-1));
                counter[i][j][k]++;
              }
            }
          }
        }
      }
    }
  }

  for(int k1 = 0; k1 < nzp; k1++){
    for(int j1 = 0; j1 < nyp; j1++){
      for(int i1 = 0; i1 < nxp; i1++){
        if(counter[i1][j1][k1]>0){
          float value = sum[i1][j1][k1]/counter[i1][j1][k1];
          cov_grid->setRealValue(i1, j1, k1, value);
        }
      }
    }
  }
  cov_grid->endAccess();
}

void
  GravimetricInversion::ExpandMatrixWithZeros(NRLib::Matrix &G, int Np, bool include_level_shift)
{
  //Expanding matrix to include padded region

  int r = G.rows().length();
  NRLib::Matrix G_star(r, Np);

  if(include_level_shift)
    G_star.resize(r, Np+1);

  // Initilize to zero
  NRLib::InitializeMatrix(G_star, 0.0);

  // Copy all elements
  for(int i = 0; i<G.rows().length(); i++)
    for(int j = 0; j<G.cols().length(); j++)
    {
      G_star(i,j) = G(i,j);
    }

    if(include_level_shift){
      // Last column with ones
      for(int i = 0; i<G.rows().length(); i++)
        G_star(i,Np) = 1;
    }

    G = G_star;
}

void
  GravimetricInversion::ExpandCovMatrixWithLevelShift(NRLib::Matrix &Sigma, double shift_parameter)
{
    int r = Sigma.rows().length();
    int c = Sigma.cols().length();

    NRLib::Matrix Sigma_star(r+1, c+1);
    NRLib::InitializeMatrix(Sigma_star, 0.0);

    // Copy all values except last row and last column - they are left to be zero.
    for(int i = 0; i<r; i++)
      for(int j = 0; j<c; j++){
        Sigma_star(i,j) = Sigma(i,j);
    }
    // Set last element
    Sigma_star(r,c) = shift_parameter;

    Sigma = Sigma_star;
}

void
  GravimetricInversion::RemoveLevelShiftFromVector(NRLib::Vector &rho, double level_shift)
{
    int r       = rho.length();

    // Level shift is found in the last element of the vector
    level_shift = rho(r-1);

    // Copy all elements except last element
    NRLib::Vector rho_new(r-1);
    for(int i = 0; i<r-1; i++){
      rho_new(i) = rho(i);
    }
    rho = rho_new;
}

void
  GravimetricInversion::RemoveLevelShiftFromCovMatrix(NRLib::Matrix &Sigma)
{
  int r = Sigma.rows().length();
  int c = Sigma.cols().length();

  // Copy all elements except last row and last column
  NRLib::Matrix new_Sigma(r-1, c-1);
  for(int i = 0; i<r-1; i++)
    for(int j = 0; j<c-1; j++){
      new_Sigma(i,j) = Sigma(i,j);
    }

  Sigma = new_Sigma;
}



void
GravimetricInversion::Divide(FFTGrid *& fftGrid_numerator, FFTGrid * fftGrid_denominator)
{
  assert(fftGrid_numerator->getNxp() == fftGrid_denominator->getNxp());

  if(fftGrid_numerator->getIsTransformed()==true && fftGrid_denominator->getIsTransformed()==true)
  {
    // Division of complex numbers:
    // \frac{a + bi}{c + di} =
    // \frac{(a + bi)(c - di)}{(c + di)(c - di)} =  // <- Multiply with conjugate denominator
    // \frac{(ac + bd)}{(c^2 + d^2)} + \frac{(bc - ad)}{(c^2 + d^2)}i

    for(int i=0;i<fftGrid_numerator->getcsize();i++)
    {
      fftw_complex numerator   = fftGrid_numerator->getNextComplex();
      fftw_complex denominator = fftGrid_denominator->getNextComplex();;

      fftw_complex tmp;
      tmp.re = numerator.re*denominator.re + numerator.im*denominator.im;
      tmp.im = numerator.im*denominator.re - numerator.re*denominator.im;

      float denominator_tmp = denominator.re*denominator.re + denominator.im*denominator.im;

      if(denominator_tmp != 0 ){
        fftw_complex answer;
        answer.re = tmp.re / denominator_tmp;
        answer.im = tmp.im / denominator_tmp;
        fftGrid_numerator->setNextComplex(answer);
      }
      else{
        // Do nothing when denominator is zero
        fftGrid_numerator->setNextComplex(numerator);
      }
    }
  }

   if(fftGrid_numerator->getIsTransformed()==false && fftGrid_denominator->getIsTransformed()==false)
   {
    for(int i=0;i < fftGrid_numerator->getrsize();i++)
    {
      float numerator   = fftGrid_numerator  ->getNextReal();
      float denominator = fftGrid_denominator->getNextReal();

      if(denominator != 0 ){
        fftGrid_numerator->setNextReal(numerator/denominator);
      }
      else{
        // Do nothing when denominator is zero
        fftGrid_numerator->setNextReal(numerator);
      }
    }
    }
}

void
  GravimetricInversion::ComputeSyntheticGravimetry(FFTGrid * rho, ModelGravityDynamic *& modelGravityDynamic, double level_shift)
{
  LogKit::WriteHeader("Compute Synthetic Gravimetry (and Residuals?)");

  NRLib::Matrix G_fullsize = modelGravityDynamic->GetGMatrixFullSize();
  int N                    = G_fullsize.cols().length();
  int n_obs                = G_fullsize.rows().length();

  NRLib::Vector rho_vec(N);
  VectorizeFFTGrid(rho_vec, rho, false); // Not use padded region

  NRLib::Vector level_shift_vec(n_obs);
  level_shift_vec.initialize(level_shift);

  NRLib::Vector synthetic_data = G_fullsize*rho_vec + level_shift_vec;
  modelGravityDynamic->SetSyntheticData(synthetic_data);

  // Dump synthetic data and full size matrix
  NRLib::WriteVectorToFile("Synthetic_data.txt",  synthetic_data);
}
*/