#include "src/crava.h"

#include "fft/include/rfftw.h"

#include "src/wavelet.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/corr.h"
#include "src/modelgeneral.h"
#include "src/modelavostatic.h"
#include "src/modelavodynamic.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/vario.h"
#include "src/welldata.h"
#include "src/krigingdata3d.h"
#include "src/covgridseparated.h"
#include "src/krigingadmin.h"
#include "src/faciesprob.h"
#include "src/definitions.h"
#include "src/gridmapping.h"
#include "src/filterwelllogs.h"
#include "src/parameteroutput.h"
#include "src/timings.h"
#include "src/spatialwellfilter.h"
#include "src/qualitygrid.h"
#include "src/io.h"
#include "src/tasklist.h"
#include "lib/timekit.hpp"
#include "lib/random.h"
#include "lib/lib_matr.h"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/grid/grid2d.hpp"


#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>

Crava::Crava(ModelSettings * modelSettings,
             ModelGeneral * modelGeneral, ModelAVOStatic * modelAVOstatic, ModelAVODynamic * modelAVOdynamic,
             SpatialWellFilter *spatwellfilter)
{
  LogKit::WriteHeader("Building Stochastic Model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  modelSettings_     = modelSettings;
  modelGeneral_      = modelGeneral;
  modelAVOstatic_    = modelAVOstatic;
  modelAVOdynamic_   = modelAVOdynamic;

  nx_                = modelAVOdynamic_->getBackAlpha()->getNx();
  ny_                = modelAVOdynamic_->getBackAlpha()->getNy();
  nz_                = modelAVOdynamic_->getBackAlpha()->getNz();
  nxp_               = modelAVOdynamic_->getBackAlpha()->getNxp();
  nyp_               = modelAVOdynamic_->getBackAlpha()->getNyp();
  nzp_               = modelAVOdynamic_->getBackAlpha()->getNzp();
  lowCut_            = modelSettings_->getLowCut();
  highCut_           = modelSettings_->getHighCut();
  wnc_               = modelSettings_->getWNC();     // white noise component see crava.h
  energyTreshold_    = modelSettings_->getEnergyThreshold();
  ntheta_            = modelSettings_->getNumberOfAngles();
  fileGrid_          = modelSettings_->getFileGrid();
  outputGridsSeismic_= modelSettings_->getOutputGridsSeismic();
  outputGridsElastic_= modelSettings_->getOutputGridsElastic();
  writePrediction_   = modelSettings_->getWritePrediction();
  krigingParameter_  = modelSettings_->getKrigingParameter();
  nWells_            = modelSettings_->getNumberOfWells();
  nSim_              = modelSettings_->getNumberOfSimulations();
  wells_             = modelAVOstatic_->getWells();
  simbox_            = modelGeneral_->getTimeSimbox();
  meanAlpha_         = modelAVOdynamic_->getBackAlpha();
  meanBeta_          = modelAVOdynamic_->getBackBeta();
  meanRho_           = modelAVOdynamic_->getBackRho();
  correlations_      = modelAVOdynamic_->getCorrelations();
  random_            = modelGeneral_->getRandomGen();
  seisWavelet_       = modelAVOdynamic_->getWavelets();
  A_                 = modelAVOdynamic_->getAMatrix();
  postAlpha_         = meanAlpha_;         // Write over the input to save memory
  postBeta_          = meanBeta_;          // Write over the input to save memory
  postRho_           = meanRho_;           // Write over the input to save memory
  fprob_             = NULL;
  thetaDeg_          = new float[ntheta_];
  empSNRatio_        = new float[ntheta_];
  theoSNRatio_       = new float[ntheta_];
  modelVariance_     = new float[ntheta_];
  signalVariance_    = new float[ntheta_];
  errorVariance_     = new float[ntheta_];
  dataVariance_      = new float[ntheta_];
  scaleWarning_      = 0;
  scaleWarningText_  = "";
  errThetaCov_       = new double*[ntheta_];
  sigmamdnew_        = NULL;
  for(int i=0;i<ntheta_;i++) {
    errThetaCov_[i]  = new double[ntheta_];
    thetaDeg_[i]     = static_cast<float>(modelSettings_->getAngle(i)*180.0/NRLib::Pi);
  }

  fftw_real * corrT = NULL; // =  fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real));

  // Double-use grids to save memory
  FFTGrid * parSpatialCorr  = NULL;  // Parameter correlation
  FFTGrid * errCorrUnsmooth = NULL;  // Error correlation

  if(!modelSettings_->getForwardModeling())
  {
    seisData_       = modelAVOdynamic_->getSeisCubes();
    modelAVOdynamic_->releaseGrids();
    correlations_->createPostGrids(nx_,ny_,nz_,nxp_,nyp_,nzp_,fileGrid_);
    parPointCov_    = correlations_->getPriorVar0();
    parSpatialCorr  = correlations_->getPostCovAlpha();  // Double-use grids to save memory
    errCorrUnsmooth = correlations_->getPostCovBeta();   // Double-use grids to save memory
    // NBNB   nzp_*0.001*corr->getdt() = T    lowCut = lowIntCut*domega = lowIntCut/T
    int lowIntCut = int(floor(lowCut_*(nzp_*0.001*correlations_->getdt())));
    // computes the integer whis corresponds to the low cut frequency.
    float corrGradI, corrGradJ;
    modelGeneral_->getCorrGradIJ(corrGradI, corrGradJ);
    corrT = parSpatialCorr->fillInParamCorr(correlations_,lowIntCut,corrGradI, corrGradJ);
    if(spatwellfilter!=NULL)
    {
      parSpatialCorr->setAccessMode(FFTGrid::RANDOMACCESS);
      for(int i=0;i<nWells_;i++)
        spatwellfilter->setPriorSpatialCorr(parSpatialCorr, wells_[i], i);
      parSpatialCorr->endAccess();
    }
    correlations_->setPriorCorrTFiltered(corrT,nz_,nzp_); // Can has zeros in the middle
    errCorrUnsmooth->fillInErrCorr(correlations_,corrGradI,corrGradJ);
    if((modelSettings_->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
      correlations_->writeFilePriorCorrT(corrT,nzp_);     // No zeros in the middle
  }
  else
  {
    modelAVOdynamic_->releaseGrids();
  }

  // reality check: all dimensions involved match
  assert(meanBeta_->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));
  assert(meanRho_->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));

  for(int i=0 ; i< ntheta_ ; i++)
  {
    if(!modelSettings_->getForwardModeling())
      assert(seisData_[i]->consistentSize(nx_,ny_,nz_,nxp_,nyp_,nzp_));
  }

  if(!modelSettings_->getForwardModeling())
  {
    parSpatialCorr->fftInPlace();
    computeVariances(corrT,modelSettings_);
    scaleWarning_ = checkScale();  // fills in scaleWarningText_ if needed.
    fftw_free(corrT);
    if(simbox_->getIsConstantThick() == false)
      divideDataByScaleWavelet();
    errCorrUnsmooth->fftInPlace();
    for(int i = 0 ; i < ntheta_ ; i++)
    {
      seisData_[i]->setAccessMode(FFTGrid::RANDOMACCESS);
      seisData_[i]->fftInPlace();
      seisData_[i]->endAccess();
    }

    if ((modelSettings_->getEstimateFaciesProb() && modelSettings_->getFaciesProbRelative())
        || modelSettings_->getUseLocalNoise())
    {
      meanAlpha2_ = copyFFTGrid(meanAlpha_);
      meanBeta2_  = copyFFTGrid(meanBeta_);
      meanRho2_   = copyFFTGrid(meanRho_);
    }

    meanAlpha_->fftInPlace();
    meanBeta_ ->fftInPlace();
    meanRho_  ->fftInPlace();
  }

  Timings::setTimeStochasticModel(wall,cpu);
}

Crava::~Crava()
{
  delete [] thetaDeg_;
  delete [] empSNRatio_;
  delete [] theoSNRatio_;
  delete [] modelVariance_;
  delete [] signalVariance_;
  delete [] errorVariance_;
  delete [] dataVariance_;
  if(fprob_!=NULL) delete fprob_;

  for(int i = 0;i<ntheta_;i++)
    delete[] errThetaCov_[i];
  delete [] errThetaCov_;

  if(postAlpha_!=NULL) delete  postAlpha_ ;
  if(postBeta_!=NULL)  delete  postBeta_;
  if(postRho_!=NULL)   delete  postRho_ ;

  if(sigmamdnew_!=NULL)
  {
    for(int i=0;i<nx_;i++)
    {
      for(int j=0;j<ny_;j++)
      {
        if((*sigmamdnew_)(i,j)!=NULL)
        {
          for(int ii=0;ii<3;ii++)
            delete [] (*sigmamdnew_)(i,j)[ii];
          delete [] (*sigmamdnew_)(i,j);
        }
      }
    }
     delete sigmamdnew_;

  }
}

void
Crava::computeDataVariance(void)
{
  //
  // Compute variation in raw seismic
  //
  int rnxp = 2*(nxp_/2+1);
  for(int l=0 ; l < ntheta_ ; l++)
  {
    double  totvar = 0;
    long int ndata = 0;
    dataVariance_[l]=0.0;
    seisData_[l]->setAccessMode(FFTGrid::READ);
    for(int k=0 ; k<nzp_ ; k++)
    {
      double tmpvar1 = 0;
      for(int j=0;j<nyp_;j++)
      {
        double tmpvar2 = 0;
        for(int i=0; i <rnxp; i++)
        {
          float tmp=seisData_[l]->getNextReal();
          if(k < nz_ && j < ny_ &&  i < nx_ && tmp != 0.0)
          {
            tmpvar2 += double(tmp*tmp);
            ndata++;
          }
        }
        tmpvar1 += tmpvar2;
      }
      totvar += tmpvar1;
    }
    seisData_[l]->endAccess();
    dataVariance_[l] = static_cast<float>(totvar/static_cast<double>(ndata));
  }
}

void
Crava::setupErrorCorrelation(ModelSettings * modelSettings,
                             const std::vector<Grid2D *> & noiseScale)
{
  //
  //  Setup error correlation matrix
  //
  for(int l=0 ; l < ntheta_ ; l++)
  {
    empSNRatio_[l] = modelSettings->getSNRatio(l);
    if(modelSettings->getUseLocalNoise() == true) {
      double minScale = noiseScale[l]->FindMin(RMISSING);
      errorVariance_[l] = float(dataVariance_[l]*minScale/empSNRatio_[l]);
    }
    else
      errorVariance_[l] = dataVariance_[l]/empSNRatio_[l];

    if (empSNRatio_[l] < 1.1f)
    {
      LogKit::LogFormatted(LogKit::Low,"\nThe empirical signal-to-noise ratio for angle stack %d is %7.1e. Ratios smaller than\n",l+1,empSNRatio_[l]);
      LogKit::LogFormatted(LogKit::Low," 1 are illegal and CRAVA has to stop. CRAVA was for some reason not able to estimate\n");
      LogKit::LogFormatted(LogKit::Low," this ratio reliably, and you must give it as input to the model file\n\n");
      exit(1);
    }
  }

  Vario * angularCorr = modelSettings->getAngularCorr();

  for(int i = 0; i < ntheta_; i++)
    for(int j = 0; j < ntheta_; j++)
      {
        float dTheta = modelSettings->getAngle(i) - modelSettings->getAngle(j);
        errThetaCov_[i][j] = static_cast<float>(sqrt(errorVariance_[i])
                                                *sqrt(errorVariance_[j])
                                                *angularCorr->corr(dTheta,0));
      }

}

void
Crava::computeVariances(fftw_real     * corrT,
                        ModelSettings * modelSettings)
{
  computeDataVariance();

  setupErrorCorrelation(modelSettings, modelAVOdynamic_->getLocalNoiseScales());

  Wavelet1D ** errorSmooth = new Wavelet1D*[ntheta_];
  float    * paramVar    = new float[ntheta_] ;
  float    * WDCorrMVar  = new float[ntheta_];

  for(int i=0 ; i < ntheta_ ; i++)
  {
    Wavelet1D* wavelet1D=seisWavelet_[i]->getWavelet1DForErrorNorm();

    errorSmooth[i] = new Wavelet1D(wavelet1D,Wavelet::FIRSTORDERFORWARDDIFF);
    if (seisWavelet_[i]->getDim() != 1) // do not delete seisWavelet_[i]!!
      delete wavelet1D;

    std::string angle    = NRLib::ToString(thetaDeg_[i], 1);
    std::string fileName = IO::PrefixWavelet() + std::string("Diff_") + angle + IO::SuffixGeneralData();
    errorSmooth[i]->printToFile(fileName);
  }

  // Compute variation in parameters
  for(int i=0 ; i < ntheta_ ; i++)
  {
    paramVar[i]=0.0;
    for(int l=0; l<3 ; l++)
      for(int m=0 ; m<3 ; m++)
        paramVar[i] += parPointCov_[l][m]*A_[i][l]*A_[i][m];
  }

  // Compute variation in wavelet
  for(int l=0 ; l < ntheta_ ; l++)
  {
    WDCorrMVar[l] = computeWDCorrMVar(errorSmooth[l],corrT);
  }

  // Compute signal and model variance and theoretical signal-to-noise-ratio
  for(int l=0 ; l < ntheta_ ; l++)
  {
    modelVariance_[l]  = WDCorrMVar[l]*paramVar[l];
    signalVariance_[l] = errorVariance_[l] + modelVariance_[l];
  }

  for(int l=0 ; l < ntheta_ ; l++)
  {
    if (modelSettings->getMatchEnergies(l))
    {
      LogKit::LogFormatted(LogKit::Low,"Matching syntethic and empirical energies:\n");
      float gain = sqrt((errorVariance_[l]/modelVariance_[l])*(empSNRatio_[l] - 1.0f));
      seisWavelet_[l]->scale(gain);
      if((modelSettings->getWaveletOutputFlag() & IO::GLOBAL_WAVELETS) > 0 ||
          (modelSettings->getEstimationMode() && modelSettings->getEstimateWavelet(l)))
      {
        std::string angle    = NRLib::ToString(thetaDeg_[l], 1);
        std::string fileName = IO::PrefixWavelet() + std::string("EnergyMatched_") + angle;
        seisWavelet_[l]->writeWaveletToFile(fileName, 1.0); // dt_max = 1.0;
      }
      modelVariance_[l] *= gain*gain;
      signalVariance_[l] = errorVariance_[l] + modelVariance_[l];
    }
    theoSNRatio_[l] = signalVariance_[l]/errorVariance_[l];
  }

  delete [] paramVar ;
  delete [] WDCorrMVar;
  for(int i=0;i<ntheta_;i++)
    delete errorSmooth[i];
  delete [] errorSmooth;
}


int
Crava::checkScale(void)
{
  std::string scaleWarning1;
  std::string scaleWarning2;
  std::string scaleWarning3;
  std::string scaleWarning4;

  scaleWarning1 = "The observed variability in seismic data is much larger than in the model.\n   Variance of: \n";
  scaleWarning2 = "The observed variability in seismic data is much less than in the model.\n   Variance of: \n";
  scaleWarning3 = "Small signal to noise ratio detected";
  scaleWarning4 = "Large signal to noise ratio detected";

  bool thisThetaIsOk;
  int  isOk=0;

  for(int l=0 ; l < ntheta_ ; l++)
  {
    thisThetaIsOk=true;

    if(( dataVariance_[l] > 4 * signalVariance_[l]) && thisThetaIsOk) //1 var 4
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 1;
        scaleWarningText_ = "Model inconsistency in angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"for seismic data\n"+scaleWarning1+"\n";
      }
      else
      {
        isOk = 1;
        scaleWarningText_ += "Model inconsistency in angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"for seismic data\n"+scaleWarning1+"\n";
      }
    }
    if( (dataVariance_[l] < 0.1 * signalVariance_[l]) && thisThetaIsOk) //1 var 0.1
    {
      thisThetaIsOk=false;
      if(isOk==0)
      {
        isOk = 2;
        scaleWarningText_ = "Model inconsistency in angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"for seismic data\n"+scaleWarning2+"\n";
      }
      else
      {
        isOk = 2;
        scaleWarningText_ += "Model inconsistency in angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"for seismic data\n"+scaleWarning2+"\n";
      }
    }
    if( (modelVariance_[l] < 0.02 * errorVariance_[l]) && thisThetaIsOk)
    {
      thisThetaIsOk=false;
      if(isOk==0)
        {
        isOk = 3;
        scaleWarningText_ = scaleWarning3+" for angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"\n";
      }
      else
        {
        isOk = 3;
        scaleWarningText_ += scaleWarning3+" for angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"\n";
      }
    }
    if( (modelVariance_[l] > 50.0 * errorVariance_[l]) && thisThetaIsOk)
    {
      thisThetaIsOk=false;
      if(isOk==0)
        {
        isOk = 4;
        scaleWarningText_ = scaleWarning4+" for angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"\n";
      }
      else
      {
        isOk = 4;
        scaleWarningText_ += scaleWarning4+" for angle "+NRLib::ToString(int(thetaDeg_[l]+0.5))+"\n";
      }
    }
  }
  return isOk;
}

void
Crava:: divideDataByScaleWavelet()
{
  int i,j,k,l,flag;
  float modW, modScaleW;

  fftw_real*    rData;
  fftw_real     tmp;
  fftw_complex* cData ;

  fftw_complex scaleWVal;
  rfftwnd_plan plan1,plan2;

  rData  = static_cast<fftw_real*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)));
  cData  = reinterpret_cast<fftw_complex*>(rData);
  Wavelet* localWavelet ;

  flag   = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan1  = rfftwnd_create_plan(1,&nzp_,FFTW_REAL_TO_COMPLEX,flag);
  plan2  = rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);

  for(l=0 ; l< ntheta_ ; l++ )
  {
    std::string angle = NRLib::ToString(thetaDeg_[l], 1);
    if(ModelSettings::getDebugLevel() > 0) {
      std::string fileName = IO::PrefixOriginalSeismicData() + "With_Padding_" + angle;
      seisData_[l]->writeStormFile(fileName, simbox_, false, true, true);
    }

    seisWavelet_[l]->fft1DInPlace();
    modW = seisWavelet_[l]->getNorm();
    modW *= modW; // note the wavelet norm is in time domain. In frequency domain we have an additional factor float(nzp_);
                  // this is because we define the wavelet as an operator hence the fft is not norm preserving.


    //double maxfrequency = double((nzp_/2)*1000.0*nz_)/(simbox_->getlz()*nzp_);
    //modW *= maxfrequency/double(highCut_); // this is the mean squared sum over reklevant frequency band.
    //                                       // Makes the problem less sensitive to the padding size

    seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
    for(i=0; i < nxp_; i++)
      for(j=0; j< nyp_; j++)
      {
        int iInd=i;
        int jInd=j;

        if(iInd > 3*nx_-1  ){
          iInd = 0;
        }
        if(jInd > 3*ny_-1  ){
          jInd = 0;
        }

        if((iInd > (nxp_+nx_)/2))
          iInd = nxp_-iInd;
        if(iInd >= nx_ )
          iInd = 2*nx_-iInd-1;

        if(jInd > (nyp_+ny_)/2)
          jInd = nyp_-jInd;
        if(jInd >= ny_ )
          jInd = 2*ny_-jInd-1;

        localWavelet = seisWavelet_[l]->getLocalWavelet1D(iInd,jInd);  // NBNB causes difference ??


        for(k=0;k<nzp_;k++)
        {
          rData[k] = seisData_[l]->getRealValue(i,j,k, true)/static_cast<float>(sqrt(static_cast<float>(nzp_)));

          if(k > nz_)
          {
            float dist = seisData_[l]->getDistToBoundary( k, nz_, nzp_);
            rData[k] *= std::max<float>(1-dist*dist,0);
          }
        }

        rfftwnd_one_real_to_complex(plan1,rData ,cData);

        double sf     = simbox_->getRelThick(i,j)*seisWavelet_[l]->getLocalStretch(iInd,jInd);
        double relT   = simbox_->getRelThick(i,j);

        double deltaF = static_cast<double>(nz_)*1000.0/(relT*simbox_->getlz()*static_cast<double>(nzp_));

        for(k=0;k < (nzp_/2 +1);k++) // all complex values
        {
          scaleWVal =  localWavelet->getCAmp(k,static_cast<float>(sf));
          modScaleW =  scaleWVal.re * scaleWVal.re + scaleWVal.im * scaleWVal.im;
          // note scaleWVal is acctually the value of the complex conjugate
          // (see definition of getCAmp)
          // wVal      =  seisWavelet_[l]->getCAmp(k);
          // modW      =  wVal.re * wVal.re + wVal.im * wVal.im;
          //  Here we need only the modulus
          // ( see definition of getCAmp)
          if((modW > 0) && (deltaF*k < highCut_ ) && (deltaF*k > lowCut_ )) //NBNB frequency cleaning
          {
            float tolFac= 0.10f;
            if(modScaleW <  tolFac * modW)
              modScaleW =   float(sqrt(modScaleW*tolFac * modW));
              //modScaleW =   float(0.5*(modScaleW + tolFac * modW));
            if(modScaleW == 0)
              modScaleW = 1;
            tmp           = cData[k].re * (scaleWVal.re/modScaleW)
              -cData[k].im * (scaleWVal.im/modScaleW);
            cData[k].im   = cData[k].im * (scaleWVal.re/modScaleW)
              +cData[k].re * (scaleWVal.im/modScaleW);
            cData[k].re   = tmp;
          }
          else
          {
            cData[k].im = 0.0f;
            cData[k].re = 0.0f;
          }
        }
        delete localWavelet;
        rfftwnd_one_complex_to_real(plan2 ,cData ,rData);
        for(k=0;k<nzp_;k++)
        {
          seisData_[l]->setRealValue(i,j,k,rData[k]/static_cast<float>(sqrt(static_cast<float>(nzp_))),true);
        }
      }

      if(ModelSettings::getDebugLevel() > 0)
      {
        std::string fileName1 = IO::PrefixReflectionCoefficients() + angle;
        std::string fileName2 = IO::PrefixReflectionCoefficients() + "With_Padding_" + angle;
        std::string sgriLabel = "Reflection coefficients for incidence angle " + angle;
        seisData_[l]->writeFile(fileName1, IO::PathToDebug(), simbox_, sgriLabel);
        seisData_[l]->writeStormFile(fileName2, simbox_, false, true, true);
      }

      LogKit::LogFormatted(LogKit::Medium,"\nInterpolating reflections for angle stack "+angle+": ");
      seisData_[l]->interpolateSeismic(energyTreshold_);

      if(ModelSettings::getDebugLevel() > 0)
      {
        std::string sgriLabel = "Interpolated reflections for incidence angle "+angle;
        std::string fileName1 = IO::PrefixReflectionCoefficients() + "Interpolated_" + angle;
        std::string fileName2 = IO::PrefixReflectionCoefficients()  + "Interpolated_With_Padding_" + angle;
        seisData_[l]->writeFile(fileName1, IO::PathToDebug(), simbox_, sgriLabel);
        seisData_[l]->writeStormFile(fileName2, simbox_, false, true, true);
      }
      seisData_[l]->endAccess();
  }

  fftw_free(rData);
  fftwnd_destroy_plan(plan1);
  fftwnd_destroy_plan(plan2);
}


void
Crava::multiplyDataByScaleWaveletAndWriteToFile(const std::string & typeName)
{
  int i,j,k,l,flag;

  fftw_real*    rData;
  fftw_real     tmp;
  fftw_complex* cData;
  fftw_complex scaleWVal;
  rfftwnd_plan plan1,plan2;

  rData  = static_cast<fftw_real*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)));
  cData  = reinterpret_cast<fftw_complex*>(rData);

  flag   = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan1  = rfftwnd_create_plan(1, &nzp_ ,FFTW_REAL_TO_COMPLEX,flag);
  plan2  = rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);

  Wavelet1D* localWavelet;

  for(l=0 ; l< ntheta_ ; l++ )
  {
    seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
    seisData_[l]->invFFTInPlace();

    for(i=0; i < nx_; i++)
      for(j=0; j< ny_; j++)
      {
        float sf = static_cast<float>(simbox_->getRelThick(i,j))*seisWavelet_[l]->getLocalStretch(i,j);

        for(k=0;k<nzp_;k++)
        {
          rData[k] = seisData_[l]->getRealValue(i,j,k, true)/static_cast<float>(sqrt(static_cast<float>(nzp_)));
        }

        rfftwnd_one_real_to_complex(plan1,rData ,cData);
        localWavelet = seisWavelet_[l]->getLocalWavelet1D(i,j);

        for(k=0;k < (nzp_/2 +1);k++) // all complex values
        {
            scaleWVal    =  localWavelet->getCAmp(k,sf);    // NBNB change here
          //scaleWVal    =  localWavelet->getCAmp(k);
          // note scaleWVal is acctually the value of the complex conjugate
          // (see definition of getCAmp)
          tmp           = cData[k].re * scaleWVal.re + cData[k].im * scaleWVal.im;
          cData[k].im   = cData[k].im * scaleWVal.re - cData[k].re * scaleWVal.im;
          cData[k].re   = tmp;
        }
        delete localWavelet;
        rfftwnd_one_complex_to_real(plan2 ,cData ,rData);
        for(k=0;k<nzp_;k++)
        {
          seisData_[l]->setRealValue(i,j,k,rData[k]/static_cast<float>(sqrt(static_cast<double>(nzp_))),true);
        }

      }
      std::string angle     = NRLib::ToString(thetaDeg_[l],1);
      std::string sgriLabel = typeName + " for incidence angle "+angle;
      std::string fileName  = typeName + "_" + angle;
      seisData_[l]->writeFile(fileName, IO::PathToInversionResults(), simbox_, sgriLabel);
      seisData_[l]->endAccess();
  }

  fftw_free(rData);
  fftwnd_destroy_plan(plan1);
  fftwnd_destroy_plan(plan2);
}


int
Crava::computePostMeanResidAndFFTCov()
{
  LogKit::WriteHeader("Posterior model / Performing Inversion");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
  int i,j,k,l,m;

  fftw_complex * kW          = new fftw_complex[ntheta_];

  fftw_complex * errMult1    = new fftw_complex[ntheta_];
  fftw_complex * errMult2    = new fftw_complex[ntheta_];
  fftw_complex * errMult3    = new fftw_complex[ntheta_];

  fftw_complex * ijkData     = new fftw_complex[ntheta_];
  fftw_complex * ijkDataMean = new fftw_complex[ntheta_];
  fftw_complex * ijkRes      = new fftw_complex[ntheta_];
  fftw_complex * ijkMean     = new fftw_complex[3];
  fftw_complex * ijkAns      = new fftw_complex[3];
  fftw_complex   kD,kD3;
  fftw_complex   ijkParLam;
  fftw_complex   ijkErrLam;
  fftw_complex   ijkTmp;

  fftw_complex**  K  = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    K[i] = new fftw_complex[3];

  fftw_complex**  KS  = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    KS[i] = new fftw_complex[3];

  fftw_complex**  KScc  = new fftw_complex*[3]; // cc - complex conjugate (and transposed)
  for(i = 0; i < 3; i++)
    KScc[i] = new fftw_complex[ntheta_];

  fftw_complex**  parVar = new fftw_complex*[3];
  for(i = 0; i < 3; i++)
    parVar[i] = new fftw_complex[3];

  fftw_complex**  margVar = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    margVar[i] = new fftw_complex[ntheta_];

  fftw_complex**  errVar = new fftw_complex*[ntheta_];
  for(i = 0; i < ntheta_; i++)
    errVar[i] = new fftw_complex[ntheta_];

  fftw_complex** reduceVar = new fftw_complex*[3];
  for(i = 0; i < 3; i++)
    reduceVar[i]= new fftw_complex[3];

  Wavelet1D * diff1Operator = new Wavelet1D(Wavelet::FIRSTORDERFORWARDDIFF,nz_,nzp_);
  Wavelet1D * diff2Operator = new Wavelet1D(diff1Operator,Wavelet::FIRSTORDERBACKWARDDIFF);
  Wavelet1D * diff3Operator = new Wavelet1D(diff2Operator,Wavelet::FIRSTORDERCENTRALDIFF);

  diff1Operator->fft1DInPlace();
  delete diff2Operator;
  diff3Operator->fft1DInPlace();

  Wavelet1D ** errorSmooth  = new Wavelet1D*[ntheta_];
  Wavelet1D ** errorSmooth2 = new Wavelet1D*[ntheta_];
  Wavelet1D ** errorSmooth3 = new Wavelet1D*[ntheta_];

  int cnxp  = nxp_/2+1;

  for(l = 0; l < ntheta_ ; l++)
  {
    std::string angle = NRLib::ToString(thetaDeg_[l], 1);
    std::string fileName;
    seisData_[l]->setAccessMode(FFTGrid::READANDWRITE);

    Wavelet1D* wavelet1D = seisWavelet_[l]->getWavelet1DForErrorNorm(); //

    errorSmooth[l]  = new Wavelet1D(wavelet1D ,Wavelet::FIRSTORDERFORWARDDIFF);
    errorSmooth2[l] = new Wavelet1D(errorSmooth[l], Wavelet::FIRSTORDERBACKWARDDIFF);
    errorSmooth3[l] = new Wavelet1D(errorSmooth2[l],Wavelet::FIRSTORDERCENTRALDIFF);
    fileName = std::string("ErrorSmooth_") + angle + IO::SuffixGeneralData();
    errorSmooth3[l]->printToFile(fileName);
    errorSmooth3[l]->fft1DInPlace();

    fileName = IO::PrefixWavelet() + angle + IO::SuffixGeneralData();
    wavelet1D->printToFile(fileName);
    wavelet1D->fft1DInPlace();

    fileName = std::string("FourierWavelet_") + angle + IO::SuffixGeneralData();
    wavelet1D->printToFile(fileName);
    delete wavelet1D;
    delete errorSmooth[l];
    delete errorSmooth2[l];
  }

  delete[] errorSmooth;
  delete[] errorSmooth2;

  meanAlpha_->setAccessMode(FFTGrid::READANDWRITE);  //   Note
  meanBeta_ ->setAccessMode(FFTGrid::READANDWRITE);  //   the top five are over written
  meanRho_  ->setAccessMode(FFTGrid::READANDWRITE);  //   does not have the initial meaning.

  FFTGrid * parSpatialCorr     = correlations_->getPostCovAlpha(); // NB! Note double usage of postCovAlpha
  FFTGrid * errCorrUnsmooth    = correlations_->getPostCovBeta();  // NB! Note double usage of postCovBeta
  FFTGrid * postCovAlpha       = correlations_->getPostCovAlpha();
  FFTGrid * postCovBeta        = correlations_->getPostCovBeta();
  FFTGrid * postCovRho         = correlations_->getPostCovRho();
  FFTGrid * postCrCovAlphaBeta = correlations_->getPostCrCovAlphaBeta();
  FFTGrid * postCrCovAlphaRho  = correlations_->getPostCrCovAlphaRho();
  FFTGrid * postCrCovBetaRho   = correlations_->getPostCrCovBetaRho();
  parSpatialCorr    ->setAccessMode(FFTGrid::READANDWRITE);  //   after the prosessing
  errCorrUnsmooth   ->setAccessMode(FFTGrid::READANDWRITE);  //
  postCovRho        ->setAccessMode(FFTGrid::WRITE);
  postCrCovAlphaBeta->setAccessMode(FFTGrid::WRITE);
  postCrCovAlphaRho ->setAccessMode(FFTGrid::WRITE);
  postCrCovBetaRho  ->setAccessMode(FFTGrid::WRITE);

  // Computes the posterior mean first  below the covariance is computed
  // To avoid to many grids in mind at the same time
  double priorVarVp,justfactor;

  int cholFlag;
  //   long int timestart, timeend;
  //   time(&timestart);
  float realFrequency;

  Wavelet1D** seisWaveletForNorm = new Wavelet1D*[ntheta_];
  for(l = 0; l < ntheta_; l++)
  {
    seisWaveletForNorm[l]=seisWavelet_[l]->getWavelet1DForErrorNorm();
    seisWaveletForNorm[l]->fft1DInPlace();
  }

  LogKit::LogFormatted(LogKit::Low,"\nBuilding posterior distribution:");
  float monitorSize = std::max(1.0f, static_cast<float>(nzp_)*0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  for(k = 0; k < nzp_; k++)
  {
    realFrequency = static_cast<float>((nz_*1000.0f)/(simbox_->getlz()*nzp_)*std::min(k,nzp_-k)); // the physical frequency
    kD = diff1Operator->getCAmp(k);                   // defines content of kD
      if( simbox_->getIsConstantThick() == true)
      {
        // defines content of K=WDA
        fillkW(k,errMult1);                                    // errMult1 used as dummy
        lib_matrProdScalVecCpx(kD, errMult1, ntheta_);         // errMult1 used as dummy
        lib_matrProdDiagCpxR(errMult1, A_, ntheta_, 3, K);     // defines content of (WDA)     K

        // defines error-term multipliers

        fillkWNorm(k,errMult1,seisWaveletForNorm);               // defines input of  (kWNorm) errMult1
        fillkWNorm(k,errMult2,errorSmooth3);               // defines input of  (kWD3Norm) errMult2
        lib_matrFillOnesVecCpx(errMult3,ntheta_);          // defines content of errMult3

      }
      else //simbox_->getIsConstantThick() == false
      {
        kD3 = diff3Operator->getCAmp(k);          // defines  kD3

        // defines content of K = DA
        lib_matrFillValueVecCpx(kD, errMult1, ntheta_);      // errMult1 used as dummy
        lib_matrProdDiagCpxR(errMult1, A_, ntheta_, 3, K); // defines content of ( K = DA )

        // defines error-term multipliers
        lib_matrFillOnesVecCpx(errMult1,ntheta_);        // defines content of errMult1
        for(l=0; l < ntheta_; l++)
          errMult1[l].re /= seisWavelet_[l]->getNorm();  // defines content of errMult1

        lib_matrFillValueVecCpx(kD3,errMult2,ntheta_);   // defines content of errMult2
        for(l=0; l < ntheta_; l++)
        {
          errMult2[l].re  /= errorSmooth3[l]->getNorm(); // defines content of errMult2
          errMult2[l].im  /= errorSmooth3[l]->getNorm(); // defines content of errMult2
        }
        fillInverseAbskWRobust(k,errMult3,seisWaveletForNorm);// defines content of errMult3
      } //simbox_->getIsConstantThick()


    for( j = 0; j < nyp_; j++) {
      for( i = 0; i < cnxp; i++) {
        ijkMean[0] = meanAlpha_->getNextComplex();
        ijkMean[1] = meanBeta_ ->getNextComplex();
        ijkMean[2] = meanRho_  ->getNextComplex();

        for(l = 0; l < ntheta_; l++ )
        {
          ijkData[l] = seisData_[l]->getNextComplex();
          ijkRes[l]  = ijkData[l];
        }

        ijkTmp       = parSpatialCorr->getNextComplex();
        ijkParLam.re = float ( sqrt(ijkTmp.re * ijkTmp.re));
        ijkParLam.im = 0.0;

        for(l = 0; l < 3; l++ )
          for(m = 0; m < 3; m++ )
          {
            parVar[l][m].re = parPointCov_[l][m] * ijkParLam.re;
            parVar[l][m].im = 0.0;
            // if(l!=m)
            //   parVar[l][m].re *= 0.75;  //NBNB OK DEBUG TEST
          }

          priorVarVp = parVar[0][0].re;
          ijkTmp       = errCorrUnsmooth->getNextComplex();
          ijkErrLam.re = float( sqrt(ijkTmp.re * ijkTmp.re));
          ijkErrLam.im = 0.0;

          if(realFrequency > lowCut_*simbox_->getMinRelThick() &&  realFrequency < highCut_) // inverting only relevant frequencies
          {
            for(l = 0; l < ntheta_; l++ ) {
              for(m = 0; m < ntheta_; m++ )
              {        // Note we multiply kWNorm[l] and comp.conj(kWNorm[m]) hence the + and not a minus as in pure multiplication
                errVar[l][m].re  = float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * ( errMult1[l].re *  errMult1[m].re +  errMult1[l].im *  errMult1[m].im));
                errVar[l][m].re += float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * ( errMult2[l].re *  errMult2[m].re +  errMult2[l].im *  errMult2[m].im));
                if(l==m) {
                  errVar[l][m].re += float( wnc_*errThetaCov_[l][m] * errMult3[l].re  * errMult3[l].re);
                  errVar[l][m].im   = 0.0;
                }
                else {
                  errVar[l][m].im  = float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * (-errMult1[l].re * errMult1[m].im + errMult1[l].im * errMult1[m].re));
                  errVar[l][m].im += float( 0.5*(1.0-wnc_)*errThetaCov_[l][m] * ijkErrLam.re * (-errMult2[l].re * errMult2[m].im + errMult2[l].im * errMult2[m].re));
                }
              }
            }
            lib_matrProdCpx(K, parVar , ntheta_, 3 ,3, KS);              //  KS is defined here
            lib_matrProdAdjointCpx(KS, K, ntheta_, 3 ,ntheta_, margVar); // margVar = (K)S(K)' is defined here
            lib_matrAddMatCpx(errVar, ntheta_,ntheta_, margVar);         // errVar  is added to margVar = (WDA)S(WDA)'  + errVar

            cholFlag=lib_matrCholCpx(ntheta_,margVar);                   // Choleskey factor of margVar is Defined

            if(cholFlag==0)
            { // then it is ok else posterior is identical to prior

              lib_matrAdjoint(KS,ntheta_,3,KScc);                        //  WDAScc is adjoint of WDAS
              lib_matrAXeqBMatCpx(ntheta_, margVar, KS, 3);              // redefines WDAS
              lib_matrProdCpx(KScc,KS,3,ntheta_,3,reduceVar);            // defines reduceVar
              //double hj=1000000.0;
              //if(reduceVar[0][0].im!=0)
              // hj = MAXIM(reduceVar[0][0].re/reduceVar[0][0].im,-reduceVar[0][0].re/reduceVar[0][0].im); //NBNB DEBUG
              lib_matrSubtMatCpx(reduceVar,3,3,parVar);                  // redefines parVar as the posterior solution

              lib_matrProdMatVecCpx(K,ijkMean, ntheta_, 3, ijkDataMean); //  defines content of ijkDataMean
              lib_matrSubtVecCpx(ijkDataMean, ntheta_, ijkData);         //  redefines content of ijkData

              lib_matrProdAdjointMatVecCpx(KS,ijkData,3,ntheta_,ijkAns); // defines ijkAns

              lib_matrAddVecCpx(ijkAns, 3,ijkMean);                      // redefines ijkMean
              lib_matrProdMatVecCpx(K,ijkMean, ntheta_, 3, ijkData);     // redefines ijkData
              lib_matrSubtVecCpx(ijkData, ntheta_,ijkRes);               // redefines ijkRes
            }

            // quality control DEBUG
            if(priorVarVp*4 < ijkAns[0].re*ijkAns[0].re + ijkAns[0].re*ijkAns[0].re)
            {
              justfactor = sqrt(ijkAns[0].re*ijkAns[0].re + ijkAns[0].re*ijkAns[0].re)/sqrt(priorVarVp);
            }
          }
          postAlpha_->setNextComplex(ijkMean[0]);
          postBeta_ ->setNextComplex(ijkMean[1]);
          postRho_  ->setNextComplex(ijkMean[2]);
          postCovAlpha->setNextComplex(parVar[0][0]);
          postCovBeta ->setNextComplex(parVar[1][1]);
          postCovRho  ->setNextComplex(parVar[2][2]);
          postCrCovAlphaBeta->setNextComplex(parVar[0][1]);
          postCrCovAlphaRho ->setNextComplex(parVar[0][2]);
          postCrCovBetaRho  ->setNextComplex(parVar[1][2]);

          for(l=0;l<ntheta_;l++)
            seisData_[l]->setNextComplex(ijkRes[l]);
      }
    }
    // Log progress
    if (k+1 >= static_cast<int>(nextMonitor))
    {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  std::cout << "\n";

  //  time(&timeend);
  // LogKit::LogFormatted(LogKit::Low,"\n Core inversion finished after %ld seconds ***\n",timeend-timestart);
  // these does not have the initial meaning
  meanAlpha_      = NULL; // the content is taken care of by  postAlpha_
  meanBeta_       = NULL; // the content is taken care of by  postBeta_
  meanRho_        = NULL; // the content is taken care of by  postRho_
  parSpatialCorr  = NULL; // the content is taken care of by  postCovAlpha
  errCorrUnsmooth = NULL; // the content is taken care of by  postCovBeta

  postAlpha_->endAccess();
  postBeta_->endAccess();
  postRho_->endAccess();

  postCovAlpha->endAccess();
  postCovBeta->endAccess();
  postCovRho->endAccess();
  postCrCovAlphaBeta->endAccess();
  postCrCovAlphaRho->endAccess();
  postCrCovBetaRho->endAccess();

  postAlpha_->invFFTInPlace();
  postBeta_->invFFTInPlace();
  postRho_->invFFTInPlace();

  for(l=0;l<ntheta_;l++)
    seisData_[l]->endAccess();

  //Finish use of seisData_, since we need the memory.
  if((outputGridsSeismic_ & IO::RESIDUAL) > 0)
  {
    if(simbox_->getIsConstantThick() != true)
      multiplyDataByScaleWaveletAndWriteToFile("residuals");
    else
    {
      for(l=0;l<ntheta_;l++)
      {
        std::string angle     = NRLib::ToString(thetaDeg_[l],1);
        std::string sgriLabel = " Residuals for incidence angle "+angle;
        std::string fileName  = IO::PrefixResiduals() + angle;
        seisData_[l]->setAccessMode(FFTGrid::RANDOMACCESS);
        seisData_[l]->invFFTInPlace();
        seisData_[l]->writeFile(fileName, IO::PathToInversionResults(), simbox_, sgriLabel);
        seisData_[l]->endAccess();
      }
    }
  }
  for(l=0;l<ntheta_;l++)
    delete seisData_[l];
  LogKit::LogFormatted(LogKit::DebugLow,"\nDEALLOCATING: Seismic data\n");

  if(modelGeneral_->getVelocityFromInversion() == true) { //Conversion undefined until prediction ready. Complete it.
    postAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);
    postAlpha_->expTransf();
    GridMapping * tdMap = modelGeneral_->getTimeDepthMapping();
    const GridMapping * dcMap = modelGeneral_->getTimeCutMapping();
    const Simbox * timeSimbox = simbox_;
    if(dcMap != NULL)
      timeSimbox = dcMap->getSimbox();

    tdMap->setMappingFromVelocity(postAlpha_, timeSimbox);
    postAlpha_->logTransf();
    postAlpha_->endAccess();
  }

  //NBNB Anne Randi: Skaler traser ihht notat fra Hugo

  if(modelSettings_->getUseLocalNoise())
  {
    correlations_->invFFT();
    correlations_->createPostVariances();
    correlations_->FFT();
    correctAlphaBetaRho(modelSettings_);
  }

  if(writePrediction_ == true)
    ParameterOutput::writeParameters(simbox_, modelGeneral_, modelSettings_, postAlpha_, postBeta_, postRho_,
    outputGridsElastic_, fileGrid_, -1, false);

  writeBWPredicted();

  delete [] seisData_;
  delete [] kW;
  delete [] errMult1;
  delete [] errMult2;
  delete [] errMult3;
  delete [] ijkData;
  delete [] ijkDataMean;
  delete [] ijkRes;
  delete [] ijkMean ;
  delete [] ijkAns;
  delete    diff1Operator;
  delete    diff3Operator;


  for(i = 0; i < ntheta_; i++)
  {
    delete[]  K[i];
    delete[]  KS[i];
    delete[]  margVar[i] ;
    delete[] errVar[i] ;
    delete errorSmooth3[i];
    delete seisWaveletForNorm[i];
  }

  delete[] K;
  delete[] KS;
  delete[] margVar;
  delete[] errVar  ;
  delete[] errorSmooth3;
  delete[] seisWaveletForNorm;

  for(i = 0; i < 3; i++)
  {
    delete[] KScc[i];
    delete[] parVar[i] ;
    delete[] reduceVar[i];
  }
  delete[] KScc;
  delete[] parVar;
  delete[] reduceVar;

  Timings::setTimeInversion(wall,cpu);
  return(0);
}

void
Crava::doPredictionKriging()
{
  if(writePrediction_ == true) { //No need to do this if output not requested.
    double wall2=0.0, cpu2=0.0;
    TimeKit::getTime(wall2,cpu2);
    doPostKriging(*postAlpha_, *postBeta_, *postRho_);
    Timings::setTimeKrigingPred(wall2,cpu2);
    ParameterOutput::writeParameters(simbox_, modelGeneral_, modelSettings_, postAlpha_, postBeta_, postRho_,
                                     outputGridsElastic_, fileGrid_, -1, true);
  }
}

int
Crava::simulate(RandomGen * randomGen)
{
  LogKit::WriteHeader("Simulating from posterior model");

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  if(nSim_>0)
  {
    bool kriging = (krigingParameter_ > 0);
    FFTGrid * postCovAlpha       = correlations_->getPostCovAlpha();
    FFTGrid * postCovBeta        = correlations_->getPostCovBeta();
    FFTGrid * postCovRho         = correlations_->getPostCovRho();
    FFTGrid * postCrCovAlphaBeta = correlations_->getPostCrCovAlphaBeta();
    FFTGrid * postCrCovAlphaRho  = correlations_->getPostCrCovAlphaRho();
    FFTGrid * postCrCovBetaRho   = correlations_->getPostCrCovBetaRho();

    assert( postCovAlpha->getIsTransformed() );
    assert( postCovBeta->getIsTransformed() );
    assert( postCovRho->getIsTransformed() );
    assert( postCrCovAlphaBeta->getIsTransformed() );
    assert( postCrCovAlphaRho->getIsTransformed() );
    assert( postCrCovBetaRho->getIsTransformed() );

    int             simNr,i,j,k,l;
    fftw_complex ** ijkPostCov;
    fftw_complex *  ijkSeed;
    FFTGrid *       seed0;
    FFTGrid *       seed1;
    FFTGrid *       seed2;

    ijkPostCov = new fftw_complex*[3];
    for(l=0;l<3;l++)
      ijkPostCov[l]=new fftw_complex[3];

    ijkSeed = new fftw_complex[3];

    seed0 =  createFFTGrid();
    seed1 =  createFFTGrid();
    seed2 =  createFFTGrid();
    seed0->createComplexGrid();
    seed1->createComplexGrid();
    seed2->createComplexGrid();

    // long int timestart, timeend;

    for(simNr = 0; simNr < nSim_;  simNr++)
    {
      // time(&timestart);

      seed0->fillInComplexNoise(randomGen);
      seed1->fillInComplexNoise(randomGen);
      seed2->fillInComplexNoise(randomGen);

      postCovAlpha      ->setAccessMode(FFTGrid::READ);
      postCovBeta       ->setAccessMode(FFTGrid::READ);
      postCovRho        ->setAccessMode(FFTGrid::READ);
      postCrCovAlphaBeta->setAccessMode(FFTGrid::READ);
      postCrCovAlphaRho ->setAccessMode(FFTGrid::READ);
      postCrCovBetaRho  ->setAccessMode(FFTGrid::READ);
      seed0 ->setAccessMode(FFTGrid::READANDWRITE);
      seed1 ->setAccessMode(FFTGrid::READANDWRITE);
      seed2 ->setAccessMode(FFTGrid::READANDWRITE);

      int cnxp=nxp_/2+1;
      int cholFlag;
      for(k = 0; k < nzp_; k++)
        for(j = 0; j < nyp_; j++)
          for(i = 0; i < cnxp; i++)
          {
            ijkPostCov[0][0] = postCovAlpha      ->getNextComplex();
            ijkPostCov[1][1] = postCovBeta       ->getNextComplex();
            ijkPostCov[2][2] = postCovRho        ->getNextComplex();
            ijkPostCov[0][1] = postCrCovAlphaBeta->getNextComplex();
            ijkPostCov[0][2] = postCrCovAlphaRho ->getNextComplex();
            ijkPostCov[1][2] = postCrCovBetaRho  ->getNextComplex();

            ijkPostCov[1][0].re =  ijkPostCov[0][1].re;
            ijkPostCov[1][0].im = -ijkPostCov[0][1].im;
            ijkPostCov[2][0].re =  ijkPostCov[0][2].re;
            ijkPostCov[2][0].im = -ijkPostCov[0][2].im;
            ijkPostCov[2][1].re =  ijkPostCov[1][2].re;
            ijkPostCov[2][1].im = -ijkPostCov[1][2].im;

            ijkSeed[0]=seed0->getNextComplex();
            ijkSeed[1]=seed1->getNextComplex();
            ijkSeed[2]=seed2->getNextComplex();

            cholFlag = lib_matrCholCpx(3,ijkPostCov);  // Choleskey factor of posterior covariance write over ijkPostCov
            if(cholFlag == 0)
            {
              lib_matrProdCholVec(3,ijkPostCov,ijkSeed); // write over ijkSeed
            }
            else
            {
              for(l=0; l< 3;l++)
              {
                ijkSeed[l].re =0.0;
                ijkSeed[l].im = 0.0;
              }

            }
            seed0->setNextComplex(ijkSeed[0]);
            seed1->setNextComplex(ijkSeed[1]);
            seed2->setNextComplex(ijkSeed[2]);
          }

          postCovAlpha->endAccess();  //
          postCovBeta->endAccess();   //
          postCovRho->endAccess();
          postCrCovAlphaBeta->endAccess();
          postCrCovAlphaRho->endAccess();
          postCrCovBetaRho->endAccess();
          seed0->endAccess();
          seed1->endAccess();
          seed2->endAccess();

          // time(&timeend);
          // printf("Simulation in FFT domain in %ld seconds \n",timeend-timestart);
          // time(&timestart);

          seed0->setAccessMode(FFTGrid::RANDOMACCESS);
          seed0->invFFTInPlace();

          seed1->setAccessMode(FFTGrid::RANDOMACCESS);
          seed1->invFFTInPlace();


          seed2->setAccessMode(FFTGrid::RANDOMACCESS);
          seed2->invFFTInPlace();

          if(modelSettings_->getUseLocalNoise()==true)
          {
            float alpha,beta, rho;
            float alphanew, betanew, rhonew;

            for(j=0;j<ny_;j++)
              for(i=0;i<nx_;i++)
                for(k=0;k<nz_;k++)
                {
                  alpha = seed0->getRealValue(i,j,k);
                  beta = seed1->getRealValue(i,j,k);
                  rho = seed2->getRealValue(i,j,k);
                  alphanew = float((*sigmamdnew_)(i,j)[0][0]*alpha+ (*sigmamdnew_)(i,j)[0][1]*beta+(*sigmamdnew_)(i,j)[0][2]*rho);
                  betanew = float((*sigmamdnew_)(i,j)[1][0]*alpha+ (*sigmamdnew_)(i,j)[1][1]*beta+(*sigmamdnew_)(i,j)[1][2]*rho);
                  rhonew = float((*sigmamdnew_)(i,j)[2][0]*alpha+ (*sigmamdnew_)(i,j)[2][1]*beta+(*sigmamdnew_)(i,j)[2][2]*rho);
                  seed0->setRealValue(i,j,k,alphanew);
                  seed1->setRealValue(i,j,k,betanew);
                  seed2->setRealValue(i,j,k,rhonew);
                }
          }

          seed0->add(postAlpha_);
          seed0->endAccess();
          seed1->add(postBeta_);
          seed1->endAccess();
          seed2->add(postRho_);
          seed2->endAccess();

          if(kriging == true) {
            double wall2=0.0, cpu2=0.0;
            TimeKit::getTime(wall2,cpu2);
            doPostKriging(*seed0, *seed1, *seed2);
            Timings::addToTimeKrigingSim(wall2,cpu2);
          }
          ParameterOutput::writeParameters(simbox_, modelGeneral_, modelSettings_, seed0, seed1, seed2,
                                           outputGridsElastic_, fileGrid_, simNr, kriging);
          // time(&timeend);
          // printf("Back transform and write of simulation in %ld seconds \n",timeend-timestart);
    }

    delete seed0;
    delete seed1;
    delete seed2;

    for(l=0;l<3;l++)
      delete  [] ijkPostCov[l];
    delete [] ijkPostCov;
    delete [] ijkSeed;

  }
  Timings::setTimeSimulation(wall,cpu);
  return(0);
}

void
Crava::doPostKriging(FFTGrid & postAlpha,
                     FFTGrid & postBeta,
                     FFTGrid & postRho)
{

  LogKit::WriteHeader("Kriging to wells");

  CovGridSeparated covGridAlpha      (*correlations_->getPostCovAlpha()      );
  CovGridSeparated covGridBeta       (*correlations_->getPostCovBeta()       );
  CovGridSeparated covGridRho        (*correlations_->getPostCovRho()        );
  CovGridSeparated covGridCrAlphaBeta(*correlations_->getPostCrCovAlphaBeta());
  CovGridSeparated covGridCrAlphaRho (*correlations_->getPostCrCovAlphaRho() );
  CovGridSeparated covGridCrBetaRho  (*correlations_->getPostCrCovBetaRho()  );

  KrigingData3D kd(wells_, nWells_, 1); // 1 = full resolution logs

  std::string baseName = "Raw_" + IO::PrefixKrigingData() + IO::SuffixGeneralData();
  std::string fileName = IO::makeFullFileName(IO::PathToInversionResults(), baseName);
  kd.writeToFile(fileName);

  CKrigingAdmin pKriging(*simbox_,
                         kd.getData(), kd.getNumberOfData(),
                         covGridAlpha, covGridBeta, covGridRho,
                         covGridCrAlphaBeta, covGridCrAlphaRho, covGridCrBetaRho,
                         krigingParameter_);

  pKriging.KrigAll(postAlpha, postBeta, postRho, false, modelSettings_->getDebugFlag(), modelSettings_->getDoSmoothKriging());
}

FFTGrid *
Crava::computeSeismicImpedance(FFTGrid * alpha, FFTGrid * beta, FFTGrid * rho, int angle)
{
  FFTGrid * impedance = createFFTGrid();
  impedance->setType(FFTGrid::DATA);
  impedance->createRealGrid();
  impedance->setAccessMode(FFTGrid::WRITE);

  int rnxp  = alpha->getRNxp();
  alpha->setAccessMode(FFTGrid::READ);
  beta->setAccessMode(FFTGrid::READ);
  rho->setAccessMode(FFTGrid::READ);
  for(int k = 0; k < nzp_; k++) {
    for(int j = 0; j < nyp_; j++)
    {
      for(int i = 0; i < rnxp; i++)
      {
        float imp = 0;
        imp += alpha->getNextReal()*A_[angle][0];
        imp += beta->getNextReal()*A_[angle][1];
        imp += rho->getNextReal()*A_[angle][2];

        impedance->setNextReal(imp);
      }
    }
  }
  impedance->endAccess();
  alpha->endAccess();
  beta->endAccess();
  rho->endAccess();
  return(impedance);
}


void
Crava::computeSyntSeismic(FFTGrid * alpha, FFTGrid * beta, FFTGrid * rho)
{
  LogKit::WriteHeader("Compute Synthetic Seismic and Residuals");

  bool fftDomain = alpha->getIsTransformed();
  if(fftDomain == true) {
    alpha->invFFTInPlace();
    beta->invFFTInPlace();
    rho->invFFTInPlace();
  }

  for(int l=0;l<ntheta_;l++) {
    FFTGrid * imp = computeSeismicImpedance(alpha, beta, rho, l);
    imp->setAccessMode(FFTGrid::RANDOMACCESS);
    for(int i=0;i<nx_; i++) {
      for(int j=0;j<ny_;j++) {
        Wavelet1D impVec(0,nz_, nzp_);
        //impVec.setupAsVector();
        int k;
        for(k=0;k<nz_;k++){
          float value = imp->getRealValue(i, j, k, true);
          impVec.setRAmp(value, k);
        }
        //Tapering:
        float fac = 1.0f/static_cast<float>(nzp_-nz_-1);
        for(;k<nzp_;k++) {
          float value = fac*((k-nz_)*impVec.getRAmp(0)+(nzp_-k-1)*impVec.getRAmp(nz_-1));
          impVec.setRAmp(value, k);
        }
        Wavelet1D resultVec(&impVec, Wavelet::FIRSTORDERFORWARDDIFF);
        resultVec.fft1DInPlace();
        Wavelet1D * localWavelet = seisWavelet_[l]->getLocalWavelet1D(i,j);

        float sf = static_cast<float>(simbox_->getRelThick(i, j))*seisWavelet_[l]->getLocalStretch(i,j);

        for(int k=0;k<(nzp_/2 +1);k++) {
          fftw_complex r = resultVec.getCAmp(k);
          fftw_complex w = localWavelet->getCAmp(k,static_cast<float>(sf));
          fftw_complex s;
          s.re = r.re*w.re+r.im*w.im; //Use complex conjugate of w
          s.im = -r.re*w.im+r.im*w.re;
          resultVec.setCAmp(s,k);
        }
        resultVec.invFFT1DInPlace();
        for(int k=0;k<nzp_;k++){
          float value = resultVec.getRAmp(k);
          imp->setRealValue(i, j, k, value, true);
        }
      }
    }
    std::string angle     = NRLib::ToString(thetaDeg_[l],1);
    std::string sgriLabel = " Synthetic seismic for incidence angle "+angle;
    std::string fileName  = IO::PrefixSyntheticSeismicData() + angle;
    if(((modelSettings_->getOutputGridsSeismic() & IO::SYNTHETIC_SEISMIC_DATA) > 0) ||
      (modelSettings_->getForwardModeling() == true))
      imp->writeFile(fileName, IO::PathToSeismicData(), simbox_,sgriLabel);
    if((modelSettings_->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0) {
      FFTGrid seis(nx_, ny_, nz_, nxp_, nyp_, nzp_);
      std::string fileName = IO::FileTemporarySeismic()+NRLib::ToString(l)+IO::SuffixCrava();
      std::string errText;
      seis.readCravaFile(fileName, errText);
      if(errText == "") {
        seis.setAccessMode(FFTGrid::RANDOMACCESS);
        for(int k=0;k<nz_;k++) {
          for(int j=0;j<ny_;j++) {
            for(int i=0;i<nx_;i++) {
              float residual = seis.getRealValue(i, j, k) - imp->getRealValue(i,j,k);
              imp->setRealValue(i, j, k, residual);
            }
          }
        }
        sgriLabel = "Residual computed from synthetic seismic for incidence angle "+angle;
        fileName = IO::PrefixSyntheticResiduals() + angle;
        imp->writeFile(fileName, IO::PathToSeismicData(), simbox_,sgriLabel);
      }
      else {
        errText += "\nFailed to read temporary stored seismic data.\n";
        LogKit::LogMessage(LogKit::Error,errText);
      }
    }
    delete imp;
  }

  if(fftDomain == true) {
    alpha->fftInPlace();
    beta->fftInPlace();
    rho->fftInPlace();
  }
}


float
Crava::computeWDCorrMVar (Wavelet1D* WD ,fftw_real* corrT)
{
  float var = 0.0;
  int i,j,corrInd;

  for(i=0;i<nzp_;i++)
    for(j=0;j<nzp_;j++)
    {
      corrInd = std::max(i-j,j-i);
      var += WD->getRAmp(i)*corrT[corrInd]*WD->getRAmp(j);
    }
    return var;
}


void
Crava::fillkW(int k, fftw_complex* kW )
{
  int l;
  for(l = 0; l < ntheta_; l++)
  {
    kW[l].re  = float( seisWavelet_[l]->getCAmp(k).re );//
    kW[l].im  = float( seisWavelet_[l]->getCAmp(k).im );//
  }

}

void
Crava::fillkWNorm(int k, fftw_complex* kWNorm,Wavelet1D** wavelet )
{
  int l;
  for(l = 0; l < ntheta_; l++)
  {
    kWNorm[l].re   = float( wavelet[l]->getCAmp(k).re/wavelet[l]->getNorm());
    kWNorm[l].im   = float( wavelet[l]->getCAmp(k).im/wavelet[l]->getNorm());
  }
}

void
Crava::fillInverseAbskWRobust(int k, fftw_complex* invkW ,Wavelet1D** seisWaveletForNorm)
{
  int l;
  float modulus,modulusFine,maxMod;
  fftw_complex value;
  fftw_complex valueFine;
  for(l = 0; l < ntheta_; l++)
  {
    value  = seisWaveletForNorm[l]->getCAmp(k);
    valueFine = seisWaveletForNorm[l]->getCAmp(k,0.999f);// denser sampling of wavelet

    modulus      = value.re*value.re + value.im*value.im;
    modulusFine  = valueFine.re*valueFine.re + valueFine.im*valueFine.im;
    maxMod       = std::max(modulus,modulusFine);

    if(maxMod > 0.0)
    {
      invkW[l].re = float( 1.0/sqrt(maxMod) );
      invkW[l].im = 0.0f;
    }
    else
    {
      invkW[l].re  =  seisWaveletForNorm[l]->getNorm()*nzp_*nzp_*100.0f; // a big number
      invkW[l].im  =  0.0; // a big number
    }
  }
}


FFTGrid*
Crava::createFFTGrid()
{
  FFTGrid* fftGrid;

  if(fileGrid_)
    fftGrid =  new FFTFileGrid(nx_,ny_,nz_,nxp_,nyp_,nzp_);
  else
    fftGrid =  new FFTGrid(nx_,ny_,nz_,nxp_,nyp_,nzp_);

  return(fftGrid);
}

FFTGrid*
Crava::copyFFTGrid(FFTGrid * fftGridOld)
{
  FFTGrid* fftGrid;
  if(fileGrid_)
    fftGrid =  new FFTFileGrid(reinterpret_cast<FFTFileGrid*>(fftGridOld));
  else
    fftGrid =  new FFTGrid(fftGridOld);
  return(fftGrid);
}


FFTFileGrid*
Crava::copyFFTGrid(FFTFileGrid * fftGridOld)
{
  FFTFileGrid* fftGrid;
  fftGrid =  new FFTFileGrid(fftGridOld);
  return(fftGrid);
}

void
Crava::printEnergyToScreen()
{
  int i;
  LogKit::LogFormatted(LogKit::Low,"\n\n                       ");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"  Seismic %4.1f ",thetaDeg_[i]);
  LogKit::LogFormatted(LogKit::Low,"\n----------------------");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"---------------");
  LogKit::LogFormatted(LogKit::Low,"\nObserved data variance :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"    %1.3e  ",dataVariance_[i]);
  LogKit::LogFormatted(LogKit::Low,"\nModelled data variance :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"    %1.3e  ",signalVariance_[i]);
  LogKit::LogFormatted(LogKit::Low,"\nError variance         :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"    %1.3e  ",errorVariance_[i]);
  LogKit::LogFormatted(LogKit::Low,"\nWavelet scale          :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"    %2.3e  ",seisWavelet_[i]->getScale());
  LogKit::LogFormatted(LogKit::Low,"\nEmpirical S/N          :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"    %5.2f      ",empSNRatio_[i]);
  LogKit::LogFormatted(LogKit::Low,"\nModelled S/N           :");
  for(i=0;i < ntheta_; i++) LogKit::LogFormatted(LogKit::Low,"    %5.2f      ",theoSNRatio_[i]);
  LogKit::LogFormatted(LogKit::Low,"\n");
}


void
Crava::computeFaciesProb(SpatialWellFilter *filteredlogs, bool useFilter)
{
  ModelSettings * modelSettings = modelSettings_;

  if(modelSettings->getEstimateFaciesProb())
  {
    LogKit::WriteHeader("Facies probability volumes");

    double wall=0.0, cpu=0.0;
    TimeKit::getTime(wall,cpu);

    int nfac = modelSettings->getNumberOfFacies();

    LogKit::LogFormatted(LogKit::Low,"\nPrior facies probabilities:\n");
    LogKit::LogFormatted(LogKit::Low,"\n");
    LogKit::LogFormatted(LogKit::Low,"Facies         Probability\n");
    LogKit::LogFormatted(LogKit::Low,"--------------------------\n");
    const float * priorFacies = modelAVOstatic_->getPriorFacies();
    for(int i=0 ; i<nfac; i++) {
      LogKit::LogFormatted(LogKit::Low,"%-15s %10.4f\n",modelSettings->getFaciesName(i).c_str(),priorFacies[i]);
    }

    if (simbox_->getdz() > 4.01f) { // Require this density for estimation of facies probabilities
      LogKit::LogFormatted(LogKit::Low,"\nWARNING: The minimum sampling density is lower than 4.0. The FACIES PROBABILITIES\n");
      LogKit::LogFormatted(LogKit::Low,"         generated by CRAVA are not reliable. To get more reliable probabilities    \n");
      LogKit::LogFormatted(LogKit::Low,"         the number of layers must be increased.                                    \n");
      std::string text("");
      text += "Increase the number of layers to improve the quality of the facies probabilities.\n";
      text += "   The minimum sampling density is "+NRLib::ToString(simbox_->getdz())+", and it should be lower than 4.0.\n";
      text += "   To obtain the desired density, the number of layers should be at least "+NRLib::ToString(static_cast<int>(ceil(simbox_->GetLZ()/4.0)))+"\n";
      TaskList::addTask(text);
    }

    LogKit::LogFormatted(LogKit::Low,"\n");
    LogKit::LogFormatted(LogKit::Low,"Well                    Use    SyntheticVs    Deviated\n");
    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------------\n");
    for(int i=0 ; i<nWells_ ; i++) {
      LogKit::LogFormatted(LogKit::Low,"%-23s %3s        %3s          %3s\n",
                           wells_[i]->getWellname().c_str(),
                           ( wells_[i]->getUseForFaciesProbabilities() ? "yes" : " no" ),
                           ( wells_[i]->hasSyntheticVsLog()            ? "yes" : " no" ),
                           ( wells_[i]->isDeviated()                   ? "yes" : " no" ));
    }

    std::string baseName = IO::PrefixFaciesProbability();

    FFTGrid * likelihood = NULL;
    if((modelSettings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
      int nx = postAlpha_->getNx();
      int ny = postAlpha_->getNy();
      int nz = postAlpha_->getNz();
      if(postAlpha_->isFile()==1)
        likelihood = new FFTFileGrid(nx, ny, nz, nx, ny, nz);
      else
        likelihood = new FFTGrid(nx, ny, nz, nx, ny, nz);
      likelihood->createRealGrid(false);
    }

    if(modelSettings->getFaciesProbRelative())
    {
      meanAlpha2_->subtract(postAlpha_);
      meanAlpha2_->changeSign();
      meanBeta2_->subtract(postBeta_);
      meanBeta2_->changeSign();
      meanRho2_->subtract(postRho_);
      meanRho2_->changeSign();
      fprob_ = new FaciesProb(meanAlpha2_,
                              meanBeta2_,
                              meanRho2_,
                              nfac,
                              modelSettings->getPundef(),
                              modelAVOstatic_->getPriorFacies(),
                              modelAVOstatic_->getPriorFaciesCubes(),
                              filteredlogs->getSigmae(),
                              useFilter,
                              const_cast<const WellData **>(wells_),
                              nWells_,
                              modelAVOstatic_->getFaciesEstimInterval(),
                              simbox_->getdz(),
                              true,
                              modelSettings->getNoVsFaciesProb(),
                              this,
                              modelAVOdynamic_->getLocalNoiseScales(),
                              modelSettings_,
                              likelihood);
      delete meanAlpha2_;
      delete meanBeta2_;
      delete meanRho2_;
    }
    else
    {
      fprob_ = new FaciesProb(postAlpha_,
                              postBeta_,
                              postRho_,
                              nfac,
                              modelSettings->getPundef(),
                              modelAVOstatic_->getPriorFacies(),
                              modelAVOstatic_->getPriorFaciesCubes(),
                              filteredlogs->getSigmae(),
                              useFilter,
                              const_cast<const WellData **>(wells_),
                              nWells_,
                              modelAVOstatic_->getFaciesEstimInterval(),
                              simbox_->getdz(),
                              false,
                              modelSettings->getNoVsFaciesProb(),
                              this,
                              modelAVOdynamic_->getLocalNoiseScales(),
                              modelSettings_,
                              likelihood);
      baseName += "Absolute_";
    }
    fprob_->calculateConditionalFaciesProb(wells_,
                                           nWells_,
                                           modelAVOstatic_->getFaciesEstimInterval(),
                                           modelSettings->getFaciesNames(),
                                           simbox_->getdz());
    LogKit::LogFormatted(LogKit::Low,"\nProbability cubes done\n");


    if (modelSettings->getOutputGridsOther() & IO::FACIESPROB_WITH_UNDEF){
      for(int i=0;i<nfac;i++)
      {
        FFTGrid * grid = fprob_->getFaciesProb(i);
        std::string fileName = baseName +"With_Undef_"+ modelSettings->getFaciesName(i);
        ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, grid, fileName,"");
      }
      FFTGrid * grid = fprob_->getFaciesProbUndef();
      std::string fileName = baseName + "Undef";
      ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, grid, fileName,"");
    }

    fprob_->calculateFaciesProbGeomodel(modelAVOstatic_->getPriorFacies(),
                                        modelAVOstatic_->getPriorFaciesCubes());

    if (modelSettings->getOutputGridsOther() & IO::FACIESPROB){
      for(int i=0;i<nfac;i++)
      {
        FFTGrid * grid = fprob_->getFaciesProb(i);
        std::string fileName = baseName + modelSettings->getFaciesName(i);
        ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, grid, fileName,"");
      }
    }
    fprob_->writeBWFaciesProb(wells_, nWells_);
    std::vector<double> pValue = fprob_->calculateChiSquareTest(wells_, nWells_, modelAVOstatic_->getFaciesEstimInterval());

    if (modelSettings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID)
      QualityGrid qualityGrid(pValue, wells_, simbox_, modelSettings, modelGeneral_);

    if(likelihood != NULL) {
      for(int i=0;i<nfac;i++) {
        FFTGrid * grid = fprob_->createLHCube(likelihood, i,
                                              modelAVOstatic_->getPriorFacies(), modelAVOstatic_->getPriorFaciesCubes());
        std::string fileName = IO::PrefixLikelihood() + modelSettings->getFaciesName(i);
        ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, grid,fileName,"");
        delete grid;
      }
    }

    Timings::setTimeFaciesProb(wall,cpu);
  }
}


void
Crava::filterLogs(Simbox          * timeSimboxConstThick,
                  FilterWellLogs *& filterlogs)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);
  int relative;
  if(modelSettings_->getEstimateFaciesProb())
    relative = 1;
  else
    relative = 0;

  filterlogs = new FilterWellLogs(timeSimboxConstThick,
                                  simbox_,
                                  correlations_,
                                  nzp_, nz_,
                                  wells_, nWells_,
                                  lowCut_, highCut_,
                                  relative);
  Timings::setTimeFiltering(wall,cpu);
}

float**
Crava::getPriorVar0() const {return correlations_->getPriorVar0();}

float**
Crava::getPostVar0() const {return correlations_->getPostVar0();}



void Crava::computeG(double **G) const
{
  double **sigmam = new double*[3];
  int i,j;
  for(i=0;i<3;i++)
  {
    sigmam[i] = new double[3];
    for(j=0;j<3;j++)
      sigmam[i][j] = double(correlations_->getPriorVar0()[i][j]);
  }
  float **sigmamd = correlations_->getPostVar0();
  double **sigmadelta = new double*[3];
  for(i=0;i<3;i++)
    sigmadelta[i] = new double[3];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmadelta[i][j] = double(sigmam[i][j]-sigmamd[i][j]);

  double  * eigval      = new double[3];
  int     * error       = new int[1];
  double ** eigvec      = new double *[3];
  double ** eigvalmat   = new double *[3];
  double ** help        = new double *[3];
  double ** eigvectrans = new double *[3];
  for(i=0;i<3;i++)
  {
    eigvec[i]      = new double[3];
    eigvalmat[i]   = new double[3];
    help[i]        = new double[3];
    eigvectrans[i] = new double [3];
  }

  lib_matr_eigen(sigmadelta,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j] = eigval[i];
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmadelta);


  lib_matr_eigen(sigmam,3,eigvec,eigval,error);

  for(i=0 ; i<3 ; i++)
    delete [] sigmam[i];
  delete [] sigmam;

  double ** sigmaminv   = new double *[3];
  double ** A           = new double *[3];

  for(i=0;i<3;i++)
  {
    sigmaminv[i]   = new double[3];
    A[i]           = new double[3];
  }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0000001)
        eigvalmat[i][j] = 1.0/sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmaminv);
  lib_matr_prod(sigmaminv,sigmadelta,3,3,3,help);
  lib_matr_prod(help,sigmaminv,3,3,3,A);

  lib_matr_eigen(A,3,eigvec,eigval,error);
  lib_matr_sort3x3(eigval,eigvec);
  lib_matrTranspose(eigvec,3,3,eigvectrans);// V_A transponert

  double ** lambdag = new double*[ntheta_];
  for(i=0;i<ntheta_;i++)
    lambdag[i] = new double[3];

    for(i=0;i<ntheta_;i++)
      for(j=0;j<3;j++)
        if(i==j)
          lambdag[i][j] = sqrt(eigval[j]/(1.0-eigval[j]));
        else
          lambdag[i][j] = 0.0;


  double  * eigvale       = new double[ntheta_];
  double ** eigvece       = new double *[ntheta_];
  double ** eigvalmate    = new double*[ntheta_];
  double ** eigvecetrans  = new double *[ntheta_];
  double ** helpe         = new double *[ntheta_];


  for(i=0;i<ntheta_;i++)
  {
    eigvece[i]      = new double[ntheta_];
    eigvalmate[i]   = new double[ntheta_];
    eigvecetrans[i] = new double[ntheta_];
    helpe[i]        = new double[ntheta_];
  }
  lib_matr_eigen(errThetaCov_,ntheta_,eigvece,eigvale,error);
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      if(i==j && eigvale[i]>0.0)
        eigvalmate[i][j] = sqrt(eigvale[i]);
      else
        eigvalmate[i][j] = 0.0;

  lib_matrTranspose(eigvece,ntheta_,ntheta_,eigvecetrans);
  lib_matr_prod(eigvece,eigvalmate,ntheta_,ntheta_,ntheta_,helpe);
  lib_matr_prod(helpe,eigvecetrans,ntheta_,ntheta_,ntheta_,eigvece);
  double **help1= new double*[ntheta_];
  double **help2 = new double*[ntheta_];

  for(i=0;i<ntheta_;i++)
  {
    help1[i] = new double[3];
    help2[i] =new double[3];
  }
  lib_matr_prod(eigvece,lambdag,ntheta_,ntheta_,3,help1);
  lib_matr_prod(help1,eigvectrans,ntheta_,3,3,help2); //
  lib_matr_prod(help2,sigmaminv,ntheta_,3,3,G);

  for(i=0;i<3;i++)
  {
    delete [] eigvec[i];
    delete [] help[i];
    delete [] eigvectrans[i];
    delete [] eigvalmat[i];
    delete [] A[i];
    delete [] sigmadelta[i];
    delete [] sigmaminv[i];
  }
  delete [] eigval;
  delete [] eigvec;
  delete [] help;
  delete [] eigvectrans;
  delete [] eigvalmat;
  delete [] A;
  delete [] sigmadelta;
  delete [] sigmaminv;

  for(i=0;i<ntheta_;i++)
  {
    delete [] eigvecetrans[i];
    delete [] eigvalmate[i];
    delete [] eigvece[i];
    delete [] help1[i];
    delete [] help2[i];
    delete [] helpe[i];
    delete [] lambdag[i];
  }
  delete [] eigvecetrans;
  delete [] eigvalmate;
  delete [] eigvece;
  delete [] eigvale;
  delete [] help1;
  delete [] help2;
  delete [] helpe;
  delete [] error;
  delete [] lambdag;
}
void Crava::newPosteriorCovPointwise(double ** sigmanew, double **G, const std::vector<double> & scales,
                                     double ** sigmamdnew) const
{
  //  this function name is not suited... it returns not what we should think perhaps...
  //  sigmanew=  sqrt( (sigmaM - sigmaM|d_new )^-1 ) * sqrt( (sigmaM -s igmaM|d_old )^-1)
  //  sigmamdnew = Sqrt( Posterior covariance)

  double **sigmaenew = new double*[ntheta_];
  double **D         = new double*[ntheta_];
  double **help      = new double*[ntheta_];
  int i,j;
  for(i=0;i<ntheta_;i++)
  {
    sigmaenew[i] = new double[ntheta_];
    D[i]         = new double[ntheta_];
    help[i]      = new double[ntheta_];
  }

  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      if(i==j)
        D[i][j] = sqrt(scales[i]);
      else
        D[i][j] = 0.0;

  lib_matr_prod(D,errThetaCov_,ntheta_,ntheta_,ntheta_,help);
  lib_matr_prod(help,D,ntheta_,ntheta_,ntheta_,sigmaenew);

  double **GT    = new double *[3];
  double **help1 = new double *[ntheta_];
  double **help2 = new double *[3];
  double **help4 = new double *[3];
  for(i=0;i<3;i++)
  {
    GT[i] = new double[ntheta_];
    help2[i] = new double[ntheta_];
    help4[i] = new double[ntheta_];
  }
  for(i=0;i<ntheta_;i++)
    help1[i] = new double[3];
  lib_matrTranspose(G,ntheta_,3,GT);
  double **sigmam = new double*[3];
  for(i=0;i<3;i++)
  {
    sigmam[i] = new double[3];
    for(j=0;j<3;j++)
      sigmam[i][j] = double(correlations_->getPriorVar0()[i][j]);
  }
  lib_matr_prod(G,sigmam,ntheta_,3,3,help1);
  lib_matr_prod(help1,GT,ntheta_,3,ntheta_,help);
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      help[i][j]+=sigmaenew[i][j];
  // help = G*Sigmam*GT+sigmaE_New


  int     * error         = new int[1];
  double  * eigvale       = new double[ntheta_];
  double ** eigvece       = new double * [ntheta_];
  double ** eigvalmate    = new double * [ntheta_];
  double ** eigvecetrans  = new double * [ntheta_];

  for(i=0;i<ntheta_;i++)
  {
    eigvece[i]      = new double[ntheta_];
    eigvalmate[i]   = new double[ntheta_];
    eigvecetrans[i] = new double[ntheta_];
  }
  lib_matr_eigen(help,ntheta_,eigvece,eigvale,error);
  for(i=0;i<ntheta_;i++)
    for(j=0;j<ntheta_;j++)
      if(i==j && eigvale[i]>0.00000001)
        eigvalmate[i][j] = 1.0/eigvale[i];
      else
        eigvalmate[i][j] = 0.0;

  lib_matrTranspose(eigvece,ntheta_,ntheta_,eigvecetrans);
  lib_matr_prod(eigvece,eigvalmate,ntheta_,ntheta_,ntheta_,help);
  lib_matr_prod(help,eigvecetrans,ntheta_,ntheta_,ntheta_,eigvece); // eigvece = (G*Sigmam*GT+sigmaE_New)^-1

  double ** help3 = new double *[3];
  for(i=0;i<3;i++)
    help3[i] = new double[3];

  double **deltanew = new double *[3];
  for(i=0;i<3;i++)
    deltanew[i] = new double[3];

  lib_matr_prod(sigmam,GT,3,3,ntheta_,help4);
  lib_matr_prod(help4,eigvece,3,ntheta_,ntheta_,help2);// help2= SigmaM*GT*(G*Sigmam*GT+sigmaE_New)^-1
  lib_matr_prod(help2,G,3,ntheta_,3,help3);// help2= SigmaM*GT*(G*Sigmam*GT+sigmaE_New)^-1*G
  lib_matr_prod(help3,sigmam,3,3,3,deltanew);// delta new= SigmaM*GT*(G*Sigmam*GT+sigmaE_New)^-1*GSigmaM

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      sigmamdnew[i][j] = -deltanew[i][j];
      sigmamdnew[i][j]+=sigmam[i][j];
    }
   // sigmamdnew = SigmaM - SigmaM*GT*(G*Sigmam*GT+sigmaE_New)^-1*GSigmaM;  is the local posterior covariance.

  double  * eigval      = new double[3];
  double ** eigvalmat   = new double*[3];
  double ** eigvec      = new double *[3];
  double ** eigvectrans = new double *[3];
  for(i=0;i<3;i++)
  {
    eigvec[i] = new double[3];
    eigvalmat[i] = new double[3];
    eigvectrans[i] = new double[3];
  }

  lib_matr_eigen(sigmamdnew,3,eigvec,eigval,error); // take square root of sigmamdnew, because this is what is needed to save for later use in simulation.
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j]=sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help3);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help3,eigvectrans,3,3,3,sigmamdnew);
  // sigmamdnew = Sqrt( Posterior covariance)

  // delta new= SigmaM*GT*(G*Sigmam*GT+sigmaE_New)^-1*GSigmaM
  lib_matr_eigen(deltanew,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j]=sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;

  lib_matr_prod(eigvec,eigvalmat,3,3,3,help3);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help3,eigvectrans,3,3,3,deltanew);
  // deltanew = sqrt(  SigmaM*GT*(G*Sigmam*GT+sigmaE_New)^-1*GSigmaM)

  float **sigmamd = correlations_->getPostVar0();
  double **sigmadelta = new double*[3];

  for(i=0;i<3;i++)
    sigmadelta[i] = new double[3];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmadelta[i][j] = double(sigmam[i][j]-sigmamd[i][j]);
  // sigmadelta = sigmaM-sigmaM|d

  lib_matr_eigen(sigmadelta,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j&& eigval[i]>0.0000001)
        eigvalmat[i][j]=1.0/sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help3);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help3,eigvectrans,3,3,3,sigmadelta);
   // sigmadelta = sqrt( sigmaM-sigmaM|d ) (robustified )


  lib_matr_prod(deltanew,sigmadelta,3,3,3,sigmanew);
  // sigmanew=  sqrt( (sigmaM - sigmaM|d_new )^-1 ) * sqrt( (sigmaM -s igmaM|d_old )^-1)

  for(i=0;i<ntheta_;i++)
  {
    delete [] eigvecetrans[i];
    delete [] eigvalmate[i];
    delete [] eigvece[i];
    delete [] sigmaenew[i];
    delete [] help1[i];
    delete [] help[i];
    delete [] D[i];
  }
  delete [] eigvecetrans;
  delete [] eigvalmate;
  delete [] eigvece;
  delete [] sigmaenew;
  delete [] help1;
  delete [] help;
  delete [] D;

  for(i=0;i<3;i++)
  {
    delete [] eigvectrans[i];
    delete [] eigvalmat[i];
    delete [] eigvec[i];
    delete [] sigmadelta[i];
    delete [] sigmam[i];
    delete [] deltanew[i];
    delete [] GT[i];
    delete [] help2[i];
    delete [] help3[i];
    delete [] help4[i];
  }
  delete [] eigvectrans;
  delete [] eigvalmat;
  delete [] eigvec;
  delete [] eigval;
  delete [] sigmadelta;
  delete [] sigmam;
  delete [] deltanew;
  delete [] GT;
  delete [] help2;
  delete [] help3;
  delete [] help4;

  delete [] error;
  delete [] eigvale;
}

void
Crava::computeFilter( float ** priorCov, double ** posteriorCov,int n,double** filter) const
{
  double ** imat       = new double *[n];
  double ** priorCov2       = new double *[n];
  for(int i=0;i<n;i++)
  {
    imat[i] = new double [n];
    priorCov2[i] = new double [n];
    for(int j=0;j<n;j++)
    {
      priorCov2[i][j] =  priorCov[i][j];
      imat[i][j] = 0.0;
      if(i==j)
        imat[i][j] =1.0;
    }
  }

  lib_matrCholR(n, priorCov2);
  lib_matrAXeqBMatR(n, priorCov2, imat, n);
  lib_matr_prod(posteriorCov,imat,n,n,n,filter);
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      filter[i][j] *=-1.0;
      if(i==j)
        filter[i][j]+=1.0;
    }
    delete [] imat[i];
    delete [] priorCov2[i];
  }

  delete [] imat;
  delete [] priorCov2;
}


void Crava::correctAlphaBetaRho(ModelSettings *modelSettings)
{
  int i,j,k;
  double **G = new double*[ntheta_];
  for(i=0;i<ntheta_;i++)
    G[i] = new double[3];

  computeG(G);

  double **sigmanew = new double *[3];
  double **sigmamd = new double *[3];
  for(i=0;i<3;i++)
  {
    sigmanew[i] = new double[3];
    sigmamd[i] = new double[3];
  }
  double **sigmamdold = new double*[3];

  for(i=0;i<3;i++)
    sigmamdold[i] = new double[3];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sigmamdold[i][j] = correlations_->getPostVar0()[i][j];

  double  * eigval       = new double[3];
  double ** eigvalmat    = new double*[3];
  double ** eigvec       = new double*[3];
  double ** eigvectrans  = new double*[3];
  int     * error        = new int[1];
  double ** help         = new double*[3];
  for(i=0;i<3;i++)
  {
    eigvec[i]      = new double[3];
    eigvalmat[i]   = new double[3];
    eigvectrans[i] = new double[3];
    help[i]        = new double[3];
  }

  lib_matr_eigen(sigmamdold,3,eigvec,eigval,error);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(i==j && eigval[i]>0.0)
        eigvalmat[i][j]=1.0/sqrt(eigval[i]);
      else
        eigvalmat[i][j] = 0.0;
  lib_matr_prod(eigvec,eigvalmat,3,3,3,help);
  lib_matrTranspose(eigvec,3,3,eigvectrans);
  lib_matr_prod(help,eigvectrans,3,3,3,sigmamdold);

  float * alpha     = new float[nz_];
  float * beta      = new float[nz_];
  float * rho       = new float[nz_];
  float * meanalpha = new float[nz_];
  float * meanbeta  = new float[nz_];
  float * meanrho   = new float[nz_];

  if(modelSettings->getNumberOfSimulations()>0)
    sigmamdnew_ = new NRLib::Grid2D<double **>(nx_,ny_,NULL);
  else
    sigmamdnew_ = NULL;

  std::vector<double> minScale(modelSettings->getNumberOfAngles());

  for(int angle=0;angle<modelSettings->getNumberOfAngles();angle++)
    minScale[angle] = modelAVOdynamic_->getLocalNoiseScale(angle)->FindMin(RMISSING);

  postAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);
  postBeta_->setAccessMode(FFTGrid::RANDOMACCESS);
  postRho_->setAccessMode(FFTGrid::RANDOMACCESS);
  meanAlpha2_->setAccessMode(FFTGrid::RANDOMACCESS);
  meanBeta2_->setAccessMode(FFTGrid::RANDOMACCESS);
  meanRho2_->setAccessMode(FFTGrid::RANDOMACCESS);

  for(i=0;i<nx_;i++)
  {
    for(j=0;j<ny_;j++)
    {
      std::vector<double> scales(modelSettings->getNumberOfAngles());
      for(int angle=0;angle<modelSettings->getNumberOfAngles();angle++)
        scales[angle] = (*(modelAVOdynamic_->getLocalNoiseScale(angle)))(i, j)/minScale[angle];

      newPosteriorCovPointwise(sigmanew,G, scales, sigmamd);

      lib_matr_prod(sigmamd,sigmamdold,3,3,3,eigvec); // store product in eigvec

      if(sigmamdnew_!=NULL)
      {
        (*sigmamdnew_)(i,j) = new double*[3];
        for(int ii=0;ii<3;ii++)
        {
          (*sigmamdnew_)(i,j)[ii] = new double[3];
          for(int jj=0;jj<3;jj++)
            (*sigmamdnew_)(i,j)[ii][jj] = eigvec[ii][jj];
        }
      }

      postAlpha_->getRealTrace(alpha, i, j);
      postBeta_->getRealTrace(beta, i, j);
      postRho_->getRealTrace(rho, i, j);
      meanAlpha2_->getRealTrace(meanalpha, i, j);
      meanBeta2_->getRealTrace(meanbeta, i, j);
      meanRho2_->getRealTrace(meanrho, i, j);

      for(k=0;k<nz_;k++)
      {
        float alphadiff = alpha[k] - meanalpha[k];
        float betadiff  = beta[k]  - meanbeta[k];
        float rhodiff   = rho[k]   - meanrho[k];
        alpha[k]  = float(meanalpha[k]+sigmanew[0][0]*alphadiff + sigmanew[0][1]*betadiff + sigmanew[0][2]*rhodiff);
        beta[k]   = float(meanbeta[k] +sigmanew[1][0]*alphadiff + sigmanew[1][1]*betadiff + sigmanew[1][2]*rhodiff);
        rho[k]    = float(meanrho[k]  +sigmanew[2][0]*alphadiff + sigmanew[2][1]*betadiff + sigmanew[2][2]*rhodiff);
      }
      postAlpha_->setRealTrace(i,j, alpha);
      postBeta_->setRealTrace(i,j,beta);
      postRho_->setRealTrace(i,j,rho);
    }
  }

  postAlpha_->endAccess();
  postBeta_->endAccess();
  postRho_->endAccess();
  meanAlpha2_->endAccess();
  meanBeta2_->endAccess();
  meanRho2_->endAccess();

  if(!(modelSettings_->getEstimateFaciesProb() && modelSettings_->getFaciesProbRelative())) {
    //We do not need these, and they will not be deleted elsewhere in this case.
    delete meanAlpha2_;
    delete meanBeta2_;
    delete meanRho2_;
  }

  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] meanalpha;
  delete [] meanbeta;
  delete [] meanrho;

  for(i=0;i<3;i++)
  {
    delete [] sigmanew[i];
    delete [] sigmamd[i];
  }
  delete [] sigmanew;
  delete [] sigmamd;

  for(i=0;i<ntheta_;i++)
    delete [] G[i];
  delete [] G;

  for(i=0;i<3;i++)
  {
    delete [] eigvectrans[i];
    delete [] eigvalmat[i];
    delete [] eigvec[i];
    delete [] sigmamdold[i];
    delete [] help[i];
  }
  delete [] eigvectrans;
  delete [] eigvalmat;
  delete [] eigval;
  delete [] eigvec;
  delete [] sigmamdold;
  delete [] error;
  delete [] help;

}

void Crava::writeBWPredicted(void)
{
  int i;
  for (i=0; i<nWells_; i++)
  {
    BlockedLogs  * bw = wells_[i]->getBlockedLogsOrigThick();

    postAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);
    bw->setLogFromGrid(postAlpha_,0,1,"ALPHA_PREDICTED");
    postAlpha_->endAccess();

    postBeta_->setAccessMode(FFTGrid::RANDOMACCESS);
    bw->setLogFromGrid(postBeta_,0,1,"BETA_PREDICTED");
    postBeta_->endAccess();

    postRho_->setAccessMode(FFTGrid::RANDOMACCESS);
    bw->setLogFromGrid(postRho_,0,1,"RHO_PREDICTED");
    postRho_->endAccess();
   }
}
