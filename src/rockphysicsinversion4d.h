/***************************************************************************
*      Copyright (C) 2013 by Norwegian Computing Center and Statoil        *
***************************************************************************/
#ifndef ROCKPHYSICSINVERSION4D_H
#define ROCKPHYSICSINVERSION4D_H


#include "nrlib/grid/grid2d.hpp"
#include "nrlib/trend/trend.hpp"
#include "rfftw.h"
#include <nrlib/flens/nrlib_flens.hpp>
#include <vector>

class FFTGrid;
class Simbox;

class RockPhysicsInversion4D
{
public:
  RockPhysicsInversion4D();
  RockPhysicsInversion4D(NRLib::Vector                      priorMean,
                         NRLib::Matrix                      priorCov,
                         NRLib::Matrix                      posteriorCov,
                         std::vector<std::vector<double> >  mSamp);

  ~RockPhysicsInversion4D();
  void     makeNewPredictionTable(std::vector<std::vector<double> >  mSamp,std::vector<double>   rSamp);
  FFTGrid* makePredictions(std::vector<FFTGrid *> mu_static_,
                         std::vector<FFTGrid *> mu_dynamic_ );

  double   getPredictedValue(NRLib::Vector f);

  void     allocatePredictionTables( );
  void     fillInTable( std::vector<std::vector<double> >  mSamp,std::vector<double>   rSamp,int tableInd);
  void     smoothAllDirectionsAndNormalize();
  void     DivideAndSmoothTable(int tableInd,std::vector<std::vector<double> > priorDistribution, std::vector<fftw_complex*> smoothingFilter);
  fftw_complex*        MakeSmoothingFilter(double posteriorVariance,double  df);
  std::vector<double>  MakeGaussKernel(double mean, double variance, double minf, double  df,int nf);
  // another option is to use data reference


  void     SolveGEVProblem( NRLib::Matrix sigma_prior,
                            NRLib::Matrix sigma_posterior,
                            NRLib::Matrix & v);
  double   GetGridValue(int TableNr,int i1,int i2,int i3,int i4);
  void     SetGridValue(int TableNr,int i1,int i2,int i3,int i4,double value);
  void     AddToGridValue(int TableNr, int i1,int i2,int i3,int i4,double value);
  void     writeTableTofile(std::string fileName);

private:

  void GetLowerIndexAndW(double minValue,double maxValue,int nValue,double value,int& index, double& w);
  int GetLowerIndex(double minValue,double maxValue,int nValue,double value);

  void ClearContentInPredictionTable( );
  NRLib::Grid2D<FFTGrid *> meanRockPrediction_;

  std::vector<fftw_complex*> smoothingFilter_; // all have length nfp_/2+1
  std::vector<std::vector<double> > priorDistribution_;
  //Linear transform 1 to 4 [f1 f2 f3 f4  ]=v_*m  (v_1 4x6 m 6x1 ),
  NRLib::Matrix v_;

  // Grid resolution for the transformed variables.
  std::vector<int> nf_;
  int nfp_;

  // Limits for factor  1 to 4
  NRLib::Vector minf_;
  NRLib::Vector maxf_;
  NRLib::Vector meanf_;
  rfftwnd_plan fftplan1_;
  rfftwnd_plan fftplan2_;

};
#endif

