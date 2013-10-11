/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef GRAVIMETRICINVERSION_H
#define GRAVIMETRICINVERSION_H

#include "fftw.h"
#include "definitions.h"
#include "libs/nrlib/flens/nrlib_flens.hpp"

class ModelSettings;
class ModelGeneral;
class ModelGravityStatic;
class ModelGravityDynamic;
class FFTGrid;
class Simbox;
class SeismicParametersHolder;
class State4D;

class GravimetricInversion
{
public:
  GravimetricInversion(ModelGeneral            * modelGeneral,
                       ModelGravityStatic      * modelGravityStatic,
                       ModelGravityDynamic     *& modelGravityDynamic,
                       SeismicParametersHolder & seismicParameters,
                       ModelSettings           * modelSettings);  // NB: Using ModelSettings only for debug purposes.

  ~GravimetricInversion();

  // Subsample in Fourier domain
  static void            Subsample(FFTGrid *& upscaled_grid, FFTGrid * original_grid,
                                   int nx_up,  int ny_up,  int nz_up,
                                   int nxp_up, int nyp_up, int nzp_up);

  // Functions for transforming the mean from log-domain
  static float           GetSigmaForTransformation(FFTGrid * sigma);
  static void            MeanExpTransform         (FFTGrid * log_mean, float sigma);   // Updates FFTGrid


private:

  // Functions for transformint covariance from log-domain
  float                  GetMeanForTransformation(FFTGrid * mean);
  void                   CovExpTransform         (FFTGrid * log_cov, float mean);                 // Updates FFTGrid
  void                   CrCovExpTransform       (FFTGrid * log_cov, float mean_a, float mean_b); // Updates FFTGrid

  // Functions for transforming mean and covariance to log doamin 
  void                   MeanLogTransform (FFTGrid * mean,   float sigma);                // Updates FFTGrid
  void                   CovLogTransform  (FFTGrid * cov,    float mean);                 // Updates FFTGrid

  // Functions related to upsampling, backsampling and vectorizing of the system
  void                   Backsample(FFTGrid * upscaled_grid, FFTGrid * new_full_grid);
  void                   VectorizeFFTGrid(NRLib::Vector &vec, FFTGrid * grid, bool with_padding = true);
  void                   ReshapeVectorToFFTGrid(FFTGrid * grid, NRLib::Vector vec);
  void                   ReshapeCovAccordingToLag(NRLib::Matrix &cov_matrix, FFTGrid * cov_grid);
  void                   ReshapeCovMatrixToFFTGrid(FFTGrid * cov_grid, NRLib::Matrix cov_matrix);

  // Functions related to expanding linear system with level_shift unknown
  void                   ExpandMatrixWithZeros(NRLib::Matrix &G, int Np, bool include_level_shift);
  void                   ExpandCovMatrixWithLevelShift(NRLib::Matrix &Sigma, double shift_parameter);
  void                   RemoveLevelShiftFromVector(NRLib::Vector &rho, double level_shift);
  void                   RemoveLevelShiftFromCovMatrix(NRLib::Matrix &Sigma);

  void                   Divide(FFTGrid *& fftGrid_numerator, FFTGrid * fftGrid_denominator);

  void                   ComputeSyntheticGravimetry(FFTGrid * rho, ModelGravityDynamic *& modelGravityDynamic, double level_shift);


   std::vector<std::vector<std::vector<int> > > lag_index_;   // class variable for faster access
};

#endif
