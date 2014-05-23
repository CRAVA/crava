/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/parameteroutput.h"
#include "src/modelsettings.h"
#include "src/simbox.h"
//#include "src/modelgeneral.h"
#include "src/gridmapping.h"
#include "src/io.h"

//void
//ParameterOutput::WriteParameters(const Simbox        * simbox,
//                                 ModelGeneral         * modelGeneral,
//                                 const ModelSettings * model_settings,
//                                 FFTGrid       * vp,
//                                 FFTGrid       * vs,
//                                 FFTGrid       * rho,
//                                 int                   output_flag,
//                                 bool                  file_grid,
//                                 int                   sim_num,
//                                 bool                  kriged)
//{
//  std::string prefix;
//  std::string suffix;
//  std::string file_name;
//
//  if(sim_num >= 0) {
//    prefix = IO::PrefixSimulations();
//    suffix = "_"+NRLib::ToString(sim_num+1);
//  }
//  else {
//    prefix = IO::PrefixPredictions();
//    suffix = "";
//  }
//
//  if(kriged)
//    suffix = "_Kriged"+suffix;
//
//  if((output_flag & IO::MURHO) > 0)
//  {
//    file_name = prefix+"MuRho"+suffix;
//    computeMuRho(simbox, modelGeneral, model_settings, vp, vs, rho, file_grid, file_name);
//  }
//  if((outputFlag & IO::LAMBDARHO) > 0)
//  {
//    fileName = prefix+"LambdaRho"+suffix;
//    computeLambdaRho(simbox, modelGeneral, modelSettings, vp, vs, rho, fileGrid, fileName);
//  }
//  if((outputFlag & IO::LAMELAMBDA) > 0)
//  {
//    fileName = prefix+"LameLambda"+suffix;
//    computeLameLambda(simbox, modelGeneral, modelSettings, vp, vs, rho, fileGrid, fileName);
//  }
//  if((outputFlag & IO::LAMEMU) > 0)
//  {
//    fileName = prefix+"LameMu"+suffix;
//    computeLameMu(simbox, modelGeneral,  modelSettings, vs, rho, fileGrid, fileName);
//  }
//  if((outputFlag & IO::POISSONRATIO) > 0)
//  {
//    fileName = prefix+"PoissonRatio"+suffix;
//    computePoissonRatio(simbox, modelGeneral, modelSettings, vp, vs, fileGrid, fileName);
//  }
//  if((outputFlag & IO::AI) > 0)
//  {
//    fileName = prefix+"AI"+suffix;
//    computeAcousticImpedance(simbox, modelGeneral, modelSettings, vp, rho, fileGrid, fileName);
//  }
//  if((outputFlag & IO::SI) > 0)
//  {
//    fileName = prefix+"SI"+suffix;
//    computeShearImpedance(simbox, modelGeneral, modelSettings, vs, rho, fileGrid, fileName);
//  }
//  if((outputFlag & IO::VPVSRATIO) > 0)
//  {
//    fileName = prefix+"VpVsRatio"+suffix;
//    computeVpVsRatio(simbox, modelGeneral, modelSettings, vp, vs, fileGrid, fileName);
//  }
//  if((outputFlag & IO::VP) > 0)
//  {
//    fileName = prefix+"Vp"+suffix;
//    vp->setAccessMode(FFTGrid::RANDOMACCESS);
//    vp->expTransf();
//    writeToFile(simbox, modelGeneral, modelSettings, vp, fileName, "Inverted Vp");
//    if(simNum<0) //prediction, need grid unharmed.
//      vp->logTransf();
//    vp->endAccess();
//  }
//  if((outputFlag & IO::VS) > 0)
//  {
//    fileName = prefix+"Vs"+suffix;
//    vs->setAccessMode(FFTGrid::RANDOMACCESS);
//    vs->expTransf();
//    writeToFile(simbox, modelGeneral, modelSettings, vs, fileName, "Inverted Vs");
//    if(simNum<0) //prediction, need grid unharmed.
//      vs->logTransf();
//    vs->endAccess();
//  }
//  if((outputFlag & IO::RHO) > 0)
//  {
//    fileName = prefix+"Rho"+suffix;
//    rho->setAccessMode(FFTGrid::RANDOMACCESS);
//    rho->expTransf();
//    writeToFile(simbox, modelGeneral, modelSettings, rho, fileName, "Inverted density");
//    if(simNum<0) //prediction, need grid unharmed.
//      rho->logTransf();
//    rho->endAccess();
//  }
//}

void
ParameterOutput::WriteParameters(const Simbox        * simbox,
                                 GridMapping         * time_depth_mapping,
                                 const ModelSettings * model_settings,
                                 StormContGrid       * vp,
                                 StormContGrid       * vs,
                                 StormContGrid       * rho,
                                 int                   output_flag,
                                 //bool                  file_grid,
                                 int                   sim_num,
                                 bool                  kriged)
{
  std::string prefix;
  std::string suffix;
  std::string file_name;

  if(sim_num >= 0) {
    prefix = IO::PrefixSimulations();
    suffix = "_"+NRLib::ToString(sim_num+1);
  }
  else {
    prefix = IO::PrefixPredictions();
    suffix = "";
  }

  if(kriged)
    suffix = "_Kriged"+suffix;

  if((output_flag & IO::MURHO) > 0) {
    file_name = prefix+"MuRho"+suffix;
    ComputeMuRho(simbox, time_depth_mapping, model_settings, vp, vs, rho, file_name);
  }
  if((output_flag & IO::LAMBDARHO) > 0) {
    file_name = prefix+"LambdaRho"+suffix;
    ComputeLambdaRho(simbox, time_depth_mapping, model_settings, vp, vs, rho, file_name);
  }
  if((output_flag & IO::LAMELAMBDA) > 0) {
    file_name = prefix+"LameLambda"+suffix;
    ComputeLameLambda(simbox, time_depth_mapping, model_settings, vp, vs, rho, file_name);
  }
  if((output_flag & IO::LAMEMU) > 0) {
    file_name = prefix+"LameMu"+suffix;
    ComputeLameMu(simbox, time_depth_mapping,  model_settings, vs, rho, file_name);
  }
  if((output_flag & IO::POISSONRATIO) > 0) {
    file_name = prefix+"PoissonRatio"+suffix;
    ComputePoissonRatio(simbox, time_depth_mapping, model_settings, vp, vs, file_name);
  }
  if((output_flag & IO::AI) > 0) {
    file_name = prefix+"AI"+suffix;
    ComputeAcousticImpedance(simbox, time_depth_mapping, model_settings, vp, rho, file_name);
  }
  if((output_flag & IO::SI) > 0) {
    file_name = prefix+"SI"+suffix;
    ComputeShearImpedance(simbox, time_depth_mapping, model_settings, vs, rho, file_name);
  }
  if((output_flag & IO::VPVSRATIO) > 0) {
    file_name = prefix+"VpVsRatio"+suffix;
    ComputeVpVsRatio(simbox, time_depth_mapping, model_settings, vp, vs, file_name);
  }
  if((output_flag & IO::VP) > 0) {
    file_name = prefix+"Vp"+suffix;

    ExpTransf(vp);

    WriteToFile(simbox, time_depth_mapping, model_settings, vp, file_name, "Inverted Vp");
    //if (sim_num < 0) //prediction, need grid unharmed.
    //  vp->logTransf();

  }
  if((output_flag & IO::VS) > 0) {
    file_name = prefix+"Vs"+suffix;

    ExpTransf(vs);

    WriteToFile(simbox, time_depth_mapping, model_settings, vs, file_name, "Inverted Vs");
    //if (sim_num < 0) //prediction, need grid unharmed.
    //  vs->logTransf();

  }
  if((output_flag & IO::RHO) > 0) {
    file_name = prefix+"Rho"+suffix;

    ExpTransf(rho);

    WriteToFile(simbox, time_depth_mapping, model_settings, rho, file_name, "Inverted density");
    //if (sim_num < 0) //prediction, need grid unharmed.
    //  rho->logTransf();

  }
}

//void
//ParameterOutput::computeAcousticImpedance(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                                          FFTGrid * Vp, FFTGrid * Rho ,
//                                          bool file_grid, const std::string & file_name)
//{
//  if(Vp->getIsTransformed()) Vp->invFFTInPlace();
//  if(Rho->getIsTransformed()) Rho->invFFTInPlace();
//
//  Vp->setAccessMode(FFTGrid::READ);
//  Rho->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* prImpedance;
//  prImpedance = createFFTGrid(Vp, file_grid);
//  prImpedance->setType(FFTGrid::PARAMETER);
//  prImpedance->createRealGrid();
//  prImpedance->setAccessMode(FFTGrid::WRITE);
//
//  int i;
//  int rSize =  prImpedance->getrsize();
//  double ijkA, ijkR, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkA = Vp->getNextReal();
//    ijkR = Rho->getNextReal();
//    compVal = exp(ijkA + ijkR);
//    prImpedance->setNextReal(float( compVal));
//  }
//
//  Vp->endAccess();
//  Rho->endAccess();
//
//  prImpedance->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, prImpedance, file_name, "Acoustic Impedance");
//  delete prImpedance;
//}

void
ParameterOutput::ComputeAcousticImpedance(const Simbox        * simbox,
                                          GridMapping         * time_depth_mapping,
                                          const ModelSettings * model_settings,
                                          StormContGrid       * vp,
                                          StormContGrid       * rho,
                                          const std::string   & file_name)
{
  StormContGrid * pr_impedance = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = exp(ijk_a + ijk_r);

        pr_impedance->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, pr_impedance, file_name, "Acoustic Impedance");

  delete pr_impedance;
}


//void
//ParameterOutput::computeShearImpedance(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                                       FFTGrid * Vs, FFTGrid * Rho,
//                                       bool file_grid, const std::string & file_name)
//{
//
//  if(Vs->getIsTransformed()) Vs->invFFTInPlace();
//  if(Rho->getIsTransformed()) Rho->invFFTInPlace();
//
//  Vs->setAccessMode(FFTGrid::READ);
//  Rho->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* shImpedance;
//  shImpedance  = createFFTGrid(Vs, file_grid);
//  shImpedance->setType(FFTGrid::PARAMETER);
//  shImpedance->createRealGrid();
//  shImpedance->setAccessMode(FFTGrid::WRITE);
//  int i;
//  int rSize =  shImpedance->getrsize();
//  double ijkB, ijkR, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkB = Vs->getNextReal();
//    ijkR = Rho->getNextReal();
//    compVal = exp(ijkB + ijkR);
//    shImpedance->setNextReal(float( compVal));
//  }
//
//  Vs->endAccess();
//  Rho->endAccess();
//
//  shImpedance->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, shImpedance, file_name, "Shear impedance");
//  delete shImpedance;
//}

void
ParameterOutput::ComputeShearImpedance(const Simbox        * simbox,
                                       GridMapping         * time_depth_mapping,
                                       const ModelSettings * model_settings,
                                       StormContGrid       * vs,
                                       StormContGrid       * rho,
                                       const std::string   & file_name)
{
  StormContGrid * sh_impedance = new StormContGrid(*vs);

  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vs->GetNI(); i++) {
    for (size_t j = 0; j < vs->GetNJ(); j++) {
      for (size_t k = 0; k < vs->GetNK(); k++) {

        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = exp(ijk_b + ijk_r);

        sh_impedance->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, sh_impedance, file_name, "Shear impedance");

  delete sh_impedance;
}


//void
//ParameterOutput::computeVpVsRatio(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                                  FFTGrid * Vp, FFTGrid * Vs,
//                                  bool file_grid, const std::string & file_name)
//{
//  if(Vp->getIsTransformed()) Vp->invFFTInPlace();
//  if(Vs->getIsTransformed())  Vs->invFFTInPlace();
//
//  Vp->setAccessMode(FFTGrid::READ);
//  Vs->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* ratioVpVs;
//  ratioVpVs = createFFTGrid(Vp, file_grid);
//  ratioVpVs->setType(FFTGrid::PARAMETER);
//  ratioVpVs->createRealGrid();
//  ratioVpVs->setAccessMode(FFTGrid::WRITE);
//  int i;
//  int rSize =  ratioVpVs->getrsize();
//  double ijkA, ijkB, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkA = Vp->getNextReal();
//    ijkB = Vs->getNextReal();
//    compVal = exp(ijkA - ijkB);
//    ratioVpVs->setNextReal(float( compVal));
//  }
//
//  Vp->endAccess();
//  Vs->endAccess();
//
//  ratioVpVs->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, ratioVpVs, file_name, "Vp-Vs ratio");
//  delete ratioVpVs;
//}

void
ParameterOutput::ComputeVpVsRatio(const Simbox        * simbox,
                                  GridMapping         * time_depth_mapping,
                                  const ModelSettings * model_settings,
                                  StormContGrid       * vp,
                                  StormContGrid       * vs,
                                  const std::string   & file_name)
{
  StormContGrid * ratio_vp_vs = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_b    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_b    = vs->GetValue(i, j, k);
        comp_val = exp(ijk_a - ijk_b);

        ratio_vp_vs->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, ratio_vp_vs, file_name, "Vp-Vs ratio");

  delete ratio_vp_vs;
}



//void
//ParameterOutput::computePoissonRatio(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                                     FFTGrid * Vp, FFTGrid * Vs,
//                                     bool file_grid, const std::string & file_name)
//{
//  if(Vp->getIsTransformed()) Vp->invFFTInPlace();
//  if(Vs->getIsTransformed()) Vs->invFFTInPlace();
//
//  Vp->setAccessMode(FFTGrid::READ);
//  Vs->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* poiRat;
//  poiRat  = createFFTGrid(Vp, file_grid);
//  poiRat->setType(FFTGrid::PARAMETER);
//  poiRat->createRealGrid();
//  poiRat->setAccessMode(FFTGrid::WRITE);
//
//  int i;
//  int rSize =  poiRat->getrsize();
//  double ijkA, ijkB, compVal, vRatioSq;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkA      = Vp->getNextReal();
//    ijkB      = Vs->getNextReal();
//    vRatioSq  = exp(2*(ijkA-ijkB));
//    compVal   = 0.5*(vRatioSq - 2)/(vRatioSq - 1);
//    poiRat->setNextReal(float( compVal));
//  }
//
//  Vp->endAccess();
//  Vs->endAccess();
//
//  poiRat->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, poiRat, file_name, "Poisson ratio");
//  delete poiRat;
//}

void
ParameterOutput::ComputePoissonRatio(const Simbox        * simbox,
                                     GridMapping         * time_depth_mapping,
                                     const ModelSettings * model_settings,
                                     StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     const std::string   & file_name)
{
  StormContGrid * poi_rat = new StormContGrid(*vp);

  float ijk_a      = 0.0f;
  float ijk_b      = 0.0f;
  float comp_val   = 0.0f;
  float v_ratio_sq = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a      = vp->GetValue(i, j, k);
        ijk_b      = vs->GetValue(i, j, k);
        v_ratio_sq = exp(2*(ijk_a-ijk_b));
        comp_val   = static_cast<float>(0.5*(v_ratio_sq - 2)/(v_ratio_sq - 1));

        poi_rat->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, poi_rat, file_name, "Poisson ratio");

  delete poi_rat;
}



//void
//ParameterOutput::computeLameMu(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                               FFTGrid * Vs, FFTGrid * Rho,
//                               bool file_grid, const std::string & file_name )
//{
//  if(Vs->getIsTransformed()) Vs->invFFTInPlace();
//  if(Rho->getIsTransformed()) Rho->invFFTInPlace();
//
//  Vs->setAccessMode(FFTGrid::READ);
//  Rho->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* mu;
//  mu  = createFFTGrid(Vs, file_grid);
//  mu->setType(FFTGrid::PARAMETER);
//  mu->createRealGrid();
//  mu->setAccessMode(FFTGrid::WRITE);
//
//  int i;
//  int rSize =  mu->getrsize();
//  double ijkB, ijkR, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkB = Vs->getNextReal();
//    ijkR = Rho->getNextReal();
//    compVal = exp(ijkR+2*ijkB-13.81551); // -13.81551 in the exponent divides by 1 000 000
//    mu->setNextReal(float( compVal));
//  }
//
//  Vs->endAccess();
//  Rho->endAccess();
//  mu->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, mu, file_name, "Lame mu");
//
//  delete mu;
//}

void
ParameterOutput::ComputeLameMu(const Simbox        * simbox,
                               GridMapping         * time_depth_mapping,
                               const ModelSettings * model_settings,
                               StormContGrid       * vs,
                               StormContGrid       * rho,
                               const std::string   & file_name)
{
  StormContGrid * mu = new StormContGrid(*mu);

  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vs->GetNI(); i++) {
    for (size_t j = 0; j < vs->GetNJ(); j++) {
      for (size_t k = 0; k < vs->GetNK(); k++) {

        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = static_cast<float>(exp(ijk_r+2*ijk_b-13.81551)); // -13.81551 in the exponent divides by 1 000 000

        mu->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, mu, file_name, "Lame mu");

  delete mu;
}

//void
//ParameterOutput::computeLameLambda(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                                   FFTGrid * Vp, FFTGrid * Vs, FFTGrid * Rho,
//                                   bool file_grid, const std::string & file_name)
//{
//  if(Vp->getIsTransformed()) Vp->invFFTInPlace();
//  if(Vs->getIsTransformed()) Vs->invFFTInPlace();
//  if(Rho->getIsTransformed()) Rho->invFFTInPlace();
//
//  Vp->setAccessMode(FFTGrid::READ);
//  Vs->setAccessMode(FFTGrid::READ);
//  Rho->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* lambda;
//  lambda  = createFFTGrid(Vp, file_grid);
//  lambda->setType(FFTGrid::PARAMETER);
//  lambda->createRealGrid();
//  lambda->setAccessMode(FFTGrid::WRITE);
//
//  int i;
//  int rSize =  lambda->getrsize();
//  double ijkA, ijkB, ijkR, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkA = Vp->getNextReal();
//    ijkB = Vs->getNextReal();
//    ijkR = Rho->getNextReal();
//    compVal = exp(ijkR)*(exp(2*ijkA-13.81551)-2*exp(2*ijkB-13.81551)); // -13.81551 in the exponent divides by 1 000 000
//    lambda->setNextReal(float( compVal));
//  }
//
//  Vp->endAccess();
//  Vs->endAccess();
//  Rho->endAccess();
//
//  lambda->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, lambda, file_name, "Lame lambda");
//
//  delete lambda;
//}

void
ParameterOutput::ComputeLameLambda(const Simbox        * simbox,
                                   GridMapping         * time_depth_mapping,
                                   const ModelSettings * model_settings,
                                   StormContGrid       * vp,
                                   StormContGrid       * vs,
                                   StormContGrid       * rho,
                                   const std::string   & file_name)
{
  StormContGrid * lambda = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = static_cast<float>(exp(ijk_r)*(exp(2*ijk_a-13.81551)-2*exp(2*ijk_b-13.81551))); // -13.81551 in the exponent divides by 1 000 000

        lambda->SetValue(i, j, k, comp_val);
      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, lambda, file_name, "Lame lambda");

  delete lambda;
}

//void
//ParameterOutput::computeLambdaRho(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                                  FFTGrid * Vp, FFTGrid * Vs, FFTGrid * Rho,
//                                  bool file_grid, const std::string & file_name)
//{
//  if(Vp->getIsTransformed()) Vp->invFFTInPlace();
//  if(Vs->getIsTransformed()) Vs->invFFTInPlace();
//  if(Rho->getIsTransformed()) Rho->invFFTInPlace();
//
//  Vp->setAccessMode(FFTGrid::READ);
//  Vs->setAccessMode(FFTGrid::READ);
//  Rho->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* lambdaRho;
//  lambdaRho  = createFFTGrid(Vp, file_grid);
//  lambdaRho->setType(FFTGrid::PARAMETER);
//  lambdaRho->createRealGrid();
//  lambdaRho->setAccessMode(FFTGrid::WRITE);
//
//  int i;
//  int rSize =  lambdaRho->getrsize();
//  double ijkA, ijkB, ijkR, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkA = Vp->getNextReal();
//    ijkB = Vs->getNextReal();
//    ijkR = Rho->getNextReal();
//    compVal = exp(2.0*(ijkA +ijkR)-13.81551)-2.0*exp(2.0*(ijkB +ijkR)-13.81551); // -13.81551 in the exponent divides by 1e6=(1 000 000)
//    lambdaRho->setNextReal(float( compVal));
//  }
//
//  Vp->endAccess();
//  Vs->endAccess();
//  Rho->endAccess();
//
//  lambdaRho->endAccess();
//
//  writeToFile(simbox, modelGeneral, model_settings, lambdaRho, file_name, "Lambda rho");
//  delete lambdaRho;
//}

void
ParameterOutput::ComputeLambdaRho(const Simbox        * simbox,
                                  GridMapping         * time_depth_mapping,
                                  const ModelSettings * model_settings,
                                  StormContGrid       * vp,
                                  StormContGrid       * vs,
                                  StormContGrid       * rho,
                                  const std::string   & file_name)
{

  StormContGrid * lambda_rho = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = static_cast<float>(exp(2.0*(ijk_a +ijk_r)-13.81551)-2.0*exp(2.0*(ijk_b +ijk_r)-13.81551)); // -13.81551 in the exponent divides by 1e6=(1 000 000)

        lambda_rho->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, lambda_rho, file_name, "Lambda rho");

  delete lambda_rho;
}

//void
//ParameterOutput::computeMuRho(const Simbox * simbox, ModelGeneral * modelGeneral, const ModelSettings * model_settings,
//                              FFTGrid * Vp, FFTGrid * Vs, FFTGrid * Rho,
//                              bool file_grid, const std::string & file_name)
//{
//  if(Vs->getIsTransformed()) Vs->invFFTInPlace();
//  if(Rho->getIsTransformed()) Rho->invFFTInPlace();
//
//  Vs->setAccessMode(FFTGrid::READ);
//  Rho->setAccessMode(FFTGrid::READ);
//
//  FFTGrid* muRho;
//  muRho = createFFTGrid(Vp, file_grid);
//  muRho->setType(FFTGrid::PARAMETER);
//  muRho->createRealGrid();
//  muRho->setAccessMode(FFTGrid::WRITE);
//
//  int i;
//  int rSize =  muRho->getrsize();
//  double ijkB, ijkR, compVal;
//  for(i=0; i  <  rSize; i++)
//  {
//    ijkB = Vs->getNextReal();
//    ijkR = Rho->getNextReal();
//    compVal = exp(2.0*(ijkB +ijkR)-13.81551); // -13.81551 in the exponent divides by 1e6=(1 000 000)
//    muRho->setNextReal(float( compVal));
//  }
//
//  Vp->endAccess();
//  Vs->endAccess();
//  Rho->endAccess();
//
//  muRho->endAccess();
//  writeToFile(simbox, modelGeneral, model_settings, muRho, file_name, "Mu rho");
//
//  delete muRho;
//}

void
ParameterOutput::ComputeMuRho(const Simbox        * simbox,
                              GridMapping         * time_depth_mapping,
                              const ModelSettings * model_settings,
                              StormContGrid       * vp,
                              StormContGrid       * vs,
                              StormContGrid       * rho,
                              const std::string   & file_name)
{
  StormContGrid * mu_rho;
  mu_rho = new StormContGrid(*vp);

  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;
  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_b    = vs->GetValue(i,j,k);
        ijk_r    = rho->GetValue(i,j,k);
        comp_val = static_cast<float>(exp(2.0*(ijk_b +ijk_r)-13.81551)); // -13.81551 in the exponent divides by 1e6=(1 000 000)

        mu_rho->SetValue(i,j,k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, mu_rho, file_name, "Mu rho");

  delete mu_rho;
}

//FFTGrid*
//ParameterOutput::createFFTGrid(FFTGrid * referenceGrid, bool file_grid)
//{
//  int nx  = referenceGrid->getNx();
//  int ny  = referenceGrid->getNy();
//  int nz  = referenceGrid->getNz();
//  int nxp = referenceGrid->getNxp();
//  int nyp = referenceGrid->getNyp();
//  int nzp = referenceGrid->getNzp();
//
//  FFTGrid * fftGrid;
//
//  if(file_grid)
//    fftGrid = new FFTFileGrid(nx,ny,nz,nxp,nyp,nzp);
//  else
//    fftGrid = new FFTGrid(nx,ny,nz,nxp,nyp,nzp);
//
//  return(fftGrid);
//}

//void
//ParameterOutput::WriteToFile(const Simbox        * simbox,
//                             GridMapping         * time_depth_mapping,
//                             const ModelSettings * model_settings,
//                             FFTGrid             * grid,
//                             const std::string   & file_name,
//                             const std::string   & sgri_label,
//                             bool                  padding)
//{
//  //GridMapping * timeDepthMapping = modelGeneral->GetTimeDepthMapping();
//  //GridMapping * timeCutMapping;//   = modelGeneral->getTimeCutMapping(); //Included in the new simbox format.
//  float seismic_start_time  = 0.0; //Hack for Sebastian, was: model->getModelSettings()->getSegyOffset();
//  TraceHeaderFormat *format = model_settings->getTraceHeaderFormatOutput();
//
//  grid->writeFile(file_name,
//                  IO::PathToInversionResults(),
//                  simbox,
//                  sgri_label,
//                  seismic_start_time,
//                  time_depth_mapping,
//                  *format,
//                  padding);
//}

void
ParameterOutput::WriteToFile(const Simbox        * simbox,
                             GridMapping         * time_depth_mapping,
                             const ModelSettings * model_settings,
                             StormContGrid       * grid,
                             const std::string   & file_name,
                             const std::string   & sgri_label,
                             bool                  padding)
{
  //GridMapping * timeDepthMapping = modelGeneral->GetTimeDepthMapping();
  //GridMapping * timeCutMapping;//   = modelGeneral->getTimeCutMapping(); //Included in the new simbox format.
  float         seismic_start_time = 0.0; //Hack for Sebastian, was: model->getModelSettings()->getSegyOffset();
  TraceHeaderFormat *format = model_settings->getTraceHeaderFormatOutput();

  //H-Writing
  WriteFile(model_settings,
            grid,
            file_name,
            IO::PathToInversionResults(),
            simbox,
            sgri_label,
            seismic_start_time,
            time_depth_mapping,
            //timeCutMapping,
            *format,
            padding);
}

void
ParameterOutput::ExpTransf(StormContGrid * grid)
{
  float value = 0.0f;
  for (size_t i = 0; i < grid->GetNI(); i++) {
    for (size_t j = 0; j < grid->GetNJ(); j++) {
      for (size_t k = 0; k < grid->GetNK(); k++) {

        value = grid->GetValue(i, j, k);

        if (value != RMISSING) {
          value = exp(value);
          grid->SetValue(i, j, k, value);
        }

      }
    }
  }

}


void
ParameterOutput::WriteFile(const ModelSettings     * model_settings,
                           StormContGrid           * storm_grid,
                           const std::string       & f_name,
                           const std::string       & sub_dir,
                           const Simbox            * simbox,
                           const std::string         label,
                           const float               z0,
                           const GridMapping       * depth_map,
                           const TraceHeaderFormat & thf,
                           bool padding)
{
  std::string file_name = IO::makeFullFileName(sub_dir, f_name);
  int format_flag       = model_settings->getOutputGridFormat();
  int domain_flag       = model_settings->getOutputGridDomain();

  if (format_flag > 0) //Output format specified.
  {
    if ((domain_flag & IO::TIMEDOMAIN) > 0) {
      //if(timeMap == NULL) { //No resampling of storm
        if ((format_flag & IO::STORM) > 0)
          storm_grid->WriteToFile(file_name, "", false);
          //FFTGrid::writeStormFile(file_name, simbox, false, padding); //H Must writeStormFile be changed to take in timeMap effects?
        if ((format_flag & IO::ASCII) > 0)
          storm_grid->WriteToFile(file_name, "", true);
          //FFTGrid::writeStormFile(file_name, simbox, true, padding);
      //}
      //else {
      //  FFTGrid::writeResampledStormCube(timeMap, fileName, simbox, formatFlag_);
      //}

      //SEGY, SGRI CRAVA are never resampled in time.
      if ((format_flag & IO::SEGY) > 0) {
        //FFTGrid::writeSegyFile(fileName, simbox, z0, thf);

        std::string file_name_segy = file_name + IO::SuffixSegy();
        LogKit::LogFormatted(LogKit::Low,"\nWriting SEGY file "+file_name_segy+"...");

        //SegY * segy = new SegY(storm_grid, z0, thf, simbox->getdz(), simbox->getIL0(), simbox->getXL0(),
        //                       simbox->getILStepX(), simbox->getILStepY(), simbox->getXLStepX(), simbox->getXLStepY());
        //segy->WriteAllTracesToFile();
        SegY * segy = new SegY(storm_grid, file_name_segy, true);

        LogKit::LogFormatted(LogKit::Low,"done\n");

        delete segy;

      }
      if ((format_flag & IO::SGRI) >0) {
        //FFTGrid::writeSgriFile(file_name, simbox, label);

        std::string file_name_sgri   = file_name + IO::SuffixSgri();
        std::string file_name_header = file_name + IO::SuffixSgriHeader();

        LogKit::LogFormatted(LogKit::Low,"\nWriting SGRI header file "+ file_name_header + "...");
        storm_grid->WriteToSgriFile(file_name_sgri, file_name_header, label, simbox->getdz());
        LogKit::LogFormatted(LogKit::Low,"done\n");

      }
      if ((format_flag & IO::CRAVA) > 0) {
        std::string crava_file_name = file_name + IO::SuffixCrava();
        storm_grid->WriteCravaFile(crava_file_name, simbox->getIL0(), simbox->getXL0(), simbox->getILStepX(), simbox->getILStepY(), simbox->getXLStepX(), simbox->getXLStepY());
      }
    }

    if (depth_map != NULL && (domain_flag & IO::DEPTHDOMAIN) > 0) { //Writing in depth. Currently, only stormfiles are written in depth.
      std::string depth_name = file_name+"_Depth";
      if (depth_map->getMapping() == NULL) {
        if (depth_map->getSimbox() == NULL) {
          LogKit::LogFormatted(LogKit::Warning,
            "WARNING: Depth interval lacking when trying to write %s. Write cancelled.\n",depth_name.c_str());
          return;
        }
        if ((format_flag & IO::STORM) > 0) {
          std::string header = depth_map->getSimbox()->getStormHeader(FFTGrid::PARAMETER, storm_grid->GetNI(), storm_grid->GetNJ(), storm_grid->GetNK(), false, false);
          storm_grid->WriteToFile(file_name, header, false);
          //FFTGrid::writeStormFile(depth_name, depth_map->getSimbox(), false);
        }
        if ((format_flag & IO::ASCII) > 0) {
          std::string header = depth_map->getSimbox()->getStormHeader(FFTGrid::PARAMETER, storm_grid->GetNI(), storm_grid->GetNJ(), storm_grid->GetNK(), false, true);
          storm_grid->WriteToFile(file_name, header, true);
          //FFTGrid::writeStormFile(depth_name, depth_map->getSimbox(), true);
        }
        if ((format_flag & IO::SEGY) >0) {
          //makeDepthCubeForSegy(depth_map->getSimbox(),depth_name);
          StormContGrid * storm_cube_depth = new StormContGrid(*depth_map->getSimbox(), storm_grid->GetNI(), storm_grid->GetNJ(), storm_grid->GetNK());

          for (size_t i = 0; i < storm_grid->GetNI(); i++) {
            for (size_t j = 0; j < storm_grid->GetNJ(); j++) {
              for (size_t k = 0; k < storm_grid->GetNK(); k++) {
                (*storm_cube_depth)(i, j, k) = storm_grid->GetValue(i, j, k);

              }
            }
          }

          std::string file_name_segy = file_name + IO::SuffixSegy();

          SegY * segy = new SegY(storm_cube_depth, file_name_segy, true);
          delete segy;
          delete storm_cube_depth;

        }
      }
      else {
        if (depth_map->getSimbox() == NULL) {
          LogKit::LogFormatted(LogKit::Warning,
            "WARNING: Depth mapping incomplete when trying to write %s. Write cancelled.\n",depth_name.c_str());
          return;
        }
        // Writes also segy in depth if required
        WriteResampledStormCube(storm_grid, depth_map, depth_name, simbox, format_flag);
      }
    }
  }
}

void
ParameterOutput::WriteResampledStormCube(const StormContGrid * storm_grid,
                                         const GridMapping   * gridmapping,
                                         const std::string   & file_name,
                                         const Simbox        * simbox,
                                         const int             format)
{
  // simbox is related to the cube we resample from. gridmapping contains simbox for the cube we resample to.

  float time, kindex;
  StormContGrid *mapping = gridmapping->getMapping();
  StormContGrid *outgrid = new StormContGrid(*mapping);

  double x,y;
  int nz = static_cast<int>(mapping->GetNK());
  for (int i = 0; i < storm_grid->GetNI(); i++) {
    for (int j = 0; j < storm_grid->GetNJ(); j++) {
      simbox->getXYCoord(i,j,x,y);
      for (int k = 0; k < nz; k++) {
        time = (*mapping)(i,j,k);
        kindex = float((time - static_cast<float>(simbox->getTop(x,y)))/simbox->getdz());

        double x_tmp, y_tmp, z_tmp;
        storm_grid->FindCenterOfCell(i, j, kindex, x_tmp, y_tmp, z_tmp);
        float value = storm_grid->GetValueZInterpolated(x_tmp, y_tmp, z_tmp);

        //float value = getRealValueInterpolated(i,j,kindex);
        (*outgrid)(i,j,k) = value;
      }
    }
  }

  std::string gf_name;
  std::string header;
  if ((format & IO::ASCII) > 0) { // ASCII
    gf_name = file_name + IO::SuffixGeneralData();
    header = gridmapping->getSimbox()->getStormHeader(FFTGrid::PARAMETER, storm_grid->GetNI(), storm_grid->GetNJ(), nz, 0, 1);
    outgrid->SetFormat(StormContGrid::STORM_ASCII);
    outgrid->WriteToFile(gf_name, header);
  }

  if ((format & IO::STORM) > 0) {
    gf_name =  file_name + IO::SuffixStormBinary();
    header = gridmapping->getSimbox()->getStormHeader(FFTGrid::PARAMETER,storm_grid->GetNI(), storm_grid->GetNJ(), nz, 0, 0);
    outgrid->SetFormat(StormContGrid::STORM_BINARY);
    outgrid->WriteToFile(gf_name,header);
  }
  if((format & IO::SEGY) > 0) {
    gf_name =  file_name + IO::SuffixSegy();

    SegY * segy = new SegY(outgrid, gf_name, true);
    delete segy;

    //writeSegyFromStorm(outgrid, gf_name);
  }
  delete outgrid;
}