#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/parameteroutput.h"
#include "src/modelsettings.h"
#include "src/simbox.h"
#include "src/model.h"
#include "src/io.h"

void
ParameterOutput::writeParameters(const Simbox  * simbox,
                                 Model         * model,
                                 FFTGrid       * alpha, 
                                 FFTGrid       * beta, 
                                 FFTGrid       * rho,
                                 int             outputFlag,
                                 bool            fileGrid,
                                 int             simNum,
                                 bool            kriged)
{
  std::string prefix;
  std::string suffix;
  std::string fileName;

  if(simNum >= 0)
  {
    prefix = IO::PrefixSimulations();
    suffix = "_"+NRLib::ToString(simNum+1);
  }
  else
  {
    prefix = IO::PrefixPredictions();
    suffix = "";
  }

  if(kriged)
    suffix = "_Kriged"+suffix;

  if((outputFlag & IO::MURHO) > 0)
  {
    fileName = prefix+"MuRho"+suffix;
    computeMuRho(simbox, model, alpha, beta, rho, fileGrid, fileName);
  }
  if((outputFlag & IO::LAMBDARHO) > 0)
  {
    fileName = prefix+"LambdaRho"+suffix;
    computeLambdaRho(simbox, model, alpha, beta, rho, fileGrid, fileName);
  }
  if((outputFlag & IO::LAMELAMBDA) > 0)
  {
    fileName = prefix+"LameLambda"+suffix;
    computeLameLambda(simbox, model, alpha, beta, rho, fileGrid, fileName);
  }
  if((outputFlag & IO::LAMEMU) > 0)
  {
    fileName = prefix+"LameMu"+suffix;
    computeLameMu(simbox, model,  beta, rho, fileGrid, fileName);
  }
  if((outputFlag & IO::POISSONRATIO) > 0)
  {
    fileName = prefix+"PoissonRatio"+suffix;
    computePoissonRatio(simbox, model, alpha, beta, fileGrid, fileName);
  }
  if((outputFlag & IO::AI) > 0)
  {
    fileName = prefix+"AI"+suffix;
    computeAcousticImpedance(simbox, model, alpha, rho, fileGrid, fileName);
  }
  if((outputFlag & IO::SI) > 0)
  {
    fileName = prefix+"SI"+suffix;
    computeShearImpedance(simbox, model, beta, rho, fileGrid, fileName);
  }
  if((outputFlag & IO::VPVSRATIO) > 0)
  {
    fileName = prefix+"VpVsRatio"+suffix;
    computeVpVsRatio(simbox, model, alpha, beta, fileGrid, fileName);
  }
  if((outputFlag & IO::VP) > 0)
  {
    fileName = prefix+"Vp"+suffix;
    alpha->setAccessMode(FFTGrid::RANDOMACCESS);
    alpha->expTransf();
    writeToFile(simbox, model, alpha, fileName, "Inverted Vp");
    if(simNum<0) //prediction, need grid unharmed.
      alpha->logTransf();
    alpha->endAccess();
  }
  if((outputFlag & IO::VS) > 0)
  {
    fileName = prefix+"Vs"+suffix;
    beta->setAccessMode(FFTGrid::RANDOMACCESS);
    beta->expTransf();
    writeToFile(simbox, model, beta, fileName, "Inverted Vs");
    if(simNum<0) //prediction, need grid unharmed.
      beta->logTransf();
    beta->endAccess();
  }
  if((outputFlag & IO::RHO) > 0)
  {
    fileName = prefix+"Rho"+suffix;
    rho->setAccessMode(FFTGrid::RANDOMACCESS);
    rho->expTransf();
    writeToFile(simbox, model, rho, fileName, "Inverted density");
    if(simNum<0) //prediction, need grid unharmed.
      rho->logTransf();
    rho->endAccess();
  }
}

void 
ParameterOutput::computeAcousticImpedance(const Simbox * simbox, Model * model, FFTGrid * Alpha, FFTGrid * Rho , 
                                          bool fileGrid, const std::string & fileName)
{   
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* prImpedance; 
  prImpedance = createFFTGrid(Alpha, fileGrid);
  prImpedance->setType(FFTGrid::PARAMETER);
  prImpedance->createRealGrid();
  prImpedance->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  prImpedance->getrsize(); 
  double ijkA, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkA + ijkR);
    prImpedance->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Rho->endAccess();

  prImpedance->endAccess();
  writeToFile(simbox, model, prImpedance, fileName, "Acoustic Impedance");
  delete prImpedance;
}

void
ParameterOutput::computeShearImpedance(const Simbox * simbox, Model * model, FFTGrid * Beta, FFTGrid * Rho,
                                       bool fileGrid, const std::string & fileName)
{

  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* shImpedance; 
  shImpedance  = createFFTGrid(Beta, fileGrid);
  shImpedance->setType(FFTGrid::PARAMETER);
  shImpedance->createRealGrid();
  shImpedance->setAccessMode(FFTGrid::WRITE);
  int i;
  int rSize =  shImpedance->getrsize(); 
  double ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkB + ijkR);
    shImpedance->setNextReal(float( compVal));
  }

  Beta->endAccess();
  Rho->endAccess();

  shImpedance->endAccess(); 
  writeToFile(simbox, model, shImpedance, fileName, "Shear impedance");
  delete shImpedance;
}


void
ParameterOutput::computeVpVsRatio(const Simbox * simbox, Model * model, FFTGrid * Alpha, FFTGrid * Beta,
                                  bool fileGrid, const std::string & fileName)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace(); 
  if(Beta->getIsTransformed())  Beta->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);

  FFTGrid* ratioVpVs; 
  ratioVpVs = createFFTGrid(Alpha, fileGrid);
  ratioVpVs->setType(FFTGrid::PARAMETER);
  ratioVpVs->createRealGrid();
  ratioVpVs->setAccessMode(FFTGrid::WRITE);
  int i;
  int rSize =  ratioVpVs->getrsize(); 
  double ijkA, ijkB, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkB = Beta->getNextReal();
    compVal = exp(ijkA - ijkB);
    ratioVpVs->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();

  ratioVpVs->endAccess();
  writeToFile(simbox, model, ratioVpVs, fileName, "Vp-Vs ratio");
  delete ratioVpVs;
}

void
ParameterOutput::computePoissonRatio(const Simbox * simbox, Model * model, FFTGrid * Alpha, FFTGrid * Beta,
                                     bool fileGrid, const std::string & fileName)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);

  FFTGrid* poiRat; 
  poiRat  = createFFTGrid(Alpha, fileGrid);
  poiRat->setType(FFTGrid::PARAMETER);
  poiRat->createRealGrid();
  poiRat->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  poiRat->getrsize(); 
  double ijkA, ijkB, compVal, vRatioSq;
  for(i=0; i  <  rSize; i++)
  {
    ijkA      = Alpha->getNextReal();
    ijkB      = Beta->getNextReal();
    vRatioSq  = exp(2*(ijkA-ijkB));
    compVal   = 0.5*(vRatioSq - 2)/(vRatioSq - 1);
    poiRat->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();

  poiRat->endAccess();
  writeToFile(simbox, model, poiRat, fileName, "Poisson ratio");
  delete poiRat;
}

void 
ParameterOutput::computeLameMu(const Simbox * simbox, Model * model, FFTGrid * Beta, FFTGrid * Rho,
                               bool fileGrid, const std::string & fileName )
{
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* mu; 
  mu  = createFFTGrid(Beta, fileGrid);
  mu->setType(FFTGrid::PARAMETER);
  mu->createRealGrid();
  mu->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  mu->getrsize(); 
  double ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkR+2*ijkB-13.81551); // -13.81551 in the exponent divides by 1 000 000
    mu->setNextReal(float( compVal));
  }

  Beta->endAccess();
  Rho->endAccess();
  mu->endAccess();
  writeToFile(simbox, model, mu, fileName, "Lame mu");

  delete mu;
}

void
ParameterOutput::computeLameLambda(const Simbox * simbox, Model * model, FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
                                   bool fileGrid, const std::string & fileName)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* lambda; 
  lambda  = createFFTGrid(Alpha, fileGrid);
  lambda->setType(FFTGrid::PARAMETER);
  lambda->createRealGrid();
  lambda->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  lambda->getrsize(); 
  double ijkA, ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(ijkR)*(exp(2*ijkA-13.81551)-exp(2*ijkB-13.81551)); // -13.81551 in the exponent divides by 1 000 000
    lambda->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  lambda->endAccess();
  writeToFile(simbox, model, lambda, fileName, "Lame lambda");

  delete lambda;
}

void
ParameterOutput::computeLambdaRho(const Simbox * simbox, Model * model, FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
                                  bool fileGrid, const std::string & fileName)
{
  if(Alpha->getIsTransformed()) Alpha->invFFTInPlace();
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Alpha->setAccessMode(FFTGrid::READ);
  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* lambdaRho; 
  lambdaRho  = createFFTGrid(Alpha, fileGrid);
  lambdaRho->setType(FFTGrid::PARAMETER);
  lambdaRho->createRealGrid();
  lambdaRho->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  lambdaRho->getrsize(); 
  double ijkA, ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkA = Alpha->getNextReal();
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(2.0*(ijkA +ijkR)-13.81551)-2.0*exp(2.0*(ijkB +ijkR)-13.81551); // -13.81551 in the exponent divides by 1e6=(1 000 000)
    lambdaRho->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  lambdaRho->endAccess();
 
  writeToFile(simbox, model, lambdaRho, fileName, "Lambda rho");
  delete lambdaRho;
}

void
ParameterOutput::computeMuRho(const Simbox * simbox, Model * model, FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho,
                              bool fileGrid, const std::string & fileName)
{
  if(Beta->getIsTransformed()) Beta->invFFTInPlace();
  if(Rho->getIsTransformed()) Rho->invFFTInPlace();

  Beta->setAccessMode(FFTGrid::READ);
  Rho->setAccessMode(FFTGrid::READ);

  FFTGrid* muRho; 
  muRho = createFFTGrid(Alpha, fileGrid);
  muRho->setType(FFTGrid::PARAMETER);
  muRho->createRealGrid();
  muRho->setAccessMode(FFTGrid::WRITE);

  int i;
  int rSize =  muRho->getrsize(); 
  double ijkB, ijkR, compVal;
  for(i=0; i  <  rSize; i++)
  {
    ijkB = Beta->getNextReal();
    ijkR = Rho->getNextReal();
    compVal = exp(2.0*(ijkB +ijkR)-13.81551); // -13.81551 in the exponent divides by 1e6=(1 000 000)
    muRho->setNextReal(float( compVal));
  }

  Alpha->endAccess();
  Beta->endAccess();
  Rho->endAccess();

  muRho->endAccess();
  writeToFile(simbox, model, muRho, fileName, "Mu rho");

  delete muRho;
}

FFTGrid*            
ParameterOutput::createFFTGrid(FFTGrid * referenceGrid, bool fileGrid)
{
  int nx  = referenceGrid->getNx();
  int ny  = referenceGrid->getNy();
  int nz  = referenceGrid->getNz(); 
  int nxp = referenceGrid->getNxp();
  int nyp = referenceGrid->getNyp();
  int nzp = referenceGrid->getNzp();

  FFTGrid * fftGrid;

  if(fileGrid)
    fftGrid = new FFTFileGrid(nx,ny,nz,nxp,nyp,nzp);
  else
    fftGrid = new FFTGrid(nx,ny,nz,nxp,nyp,nzp);

  return(fftGrid);
}

void 
ParameterOutput::writeToFile(const Simbox      * simbox, 
                             Model             * model, 
                             FFTGrid           * grid, 
                             const std::string & fileName, 
                             const std::string & sgriLabel) 
{
  GridMapping * timeDepthMapping = model->getTimeDepthMapping();
  GridMapping * timeCutMapping   = model->getTimeCutMapping();
  float         seismicStartTime = 0.0; //Hack for Sebastian, was: model->getModelSettings()->getSegyOffset();
  TraceHeaderFormat *format = model->getModelSettings()->getTraceHeaderFormatOutput();

  grid->writeFile(fileName, 
                  IO::PathToInversionResults(), 
                  simbox,
                  sgriLabel, 
                  seismicStartTime, 
                  timeDepthMapping, 
                  timeCutMapping, 
                  *format);
}
