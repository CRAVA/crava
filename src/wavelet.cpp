#include "wavelet.h"

#include <string.h>
#include <assert.h>
#include <math.h>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"
#include "lib/irapgrid.h"
#include "lib/log.h"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

Wavelet::Wavelet(Simbox         * simbox,
                 FFTGrid        * seisCube,
                 WellData      ** wells,
                 ModelSettings  * modelSettings,
                 float          * coeff) 
{
  LogKit::writeLog("\n  Estimating wavelet from seismic data and (nonfiltered) blocked wells\n");
  readtype_ = ESTIMATE;
  // initialization

  int   nWells        = modelSettings->getNumberOfWells();
  float waveletLength = modelSettings->getWaveletTaperingL();
  maxShift_           = modelSettings->getMaxWaveletShift();
  minRelativeAmp_     = modelSettings->getMinRelWaveletAmp();

  scale_=1.0f; 
  gridNI_=0;   
  gridNJ_=0;
  shiftGrid_=NULL;  
  gainGrid_=NULL; 
  errCode_=0;
  dz_        = static_cast<float>(simbox->getdz());
  nz_        = simbox->getnz();
  nzp_       = seisCube->getNzp();
  cnzp_      = nzp_/2+1;
  rnzp_	     = 2*cnzp_;
  theta_     = seisCube->getTheta();

  inFFTorder_= true;
  isReal_    = true; 
  cz_        = 0;

  float* dz = new float[nWells];
  float*  wellWeight=new float[nWells];
  int i, j, k;
  for(i=0;i<3;i++)
    coeff_[i] = coeff[i];

  int   nz   = simbox->getnz();
  float dz0  = static_cast<float>(simbox->getdz());
  int   nzp  = seisCube->getNzp();
  int   cnzp = (nzp/2+1);
  int   rnzp = 2*cnzp;

  float * alpha    = new float[nz];
  float * beta     = new float[nz];
  float * rho      = new float[nz];
  float * seisData = new float[nz];
  bool  * hasData  = new bool[nz];

  //Wavelet estimation
  fftw_real    ** cpp_r = new fftw_real*[nWells];
  fftw_complex ** cpp_c = (fftw_complex** ) cpp_r;
  
  fftw_real    ** seis_r = new fftw_real*[nWells];
  fftw_complex ** seis_c = (fftw_complex** ) seis_r;

  fftw_real    ** synt_seis_r = new fftw_real*[nWells];
  fftw_complex ** synt_seis_c = (fftw_complex** ) synt_seis_r;

  fftw_real    ** cor_cpp_r = new fftw_real*[nWells];
  fftw_complex ** cor_cpp_c = (fftw_complex** ) cor_cpp_r;
  
  fftw_real    ** ccor_seis_cpp_r = new fftw_real*[nWells];
  fftw_complex ** ccor_seis_cpp_c = (fftw_complex** ) ccor_seis_cpp_r;

  fftw_real    ** wavelet_r = new fftw_real*[nWells];
  fftw_complex ** wavelet_c = (fftw_complex** ) wavelet_r; 

  int maxBlocks = 0;
  for(i=0;i<nWells;i++)
  {
    cpp_r[i]           = new fftw_real[rnzp];
    seis_r[i]          = new fftw_real[rnzp];
    synt_seis_r[i]     = new fftw_real[rnzp];
    cor_cpp_r[i]       = new fftw_real[rnzp];
    ccor_seis_cpp_r[i] = new fftw_real[rnzp];
    wavelet_r[i]       = new fftw_real[rnzp];
    wellWeight[i]      = 0;
    const int * ipos = wells[i]->getBlockedLogsPropThick()->getIpos();
    const int * jpos = wells[i]->getBlockedLogsPropThick()->getJpos();
    dz[i] = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])*simbox->getdz());
    int nBlocks = wells[i]->getBlockedLogsPropThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++) 
  {
    if (!wells[w]->isDeviated())
    {
      BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
      //
      // Block seismic data for this well
      //
      bl->getBlockedGrid(seisCube,seisLog);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bl->getVerticalTrend(bl->getAlpha(), alpha);
      bl->getVerticalTrend(bl->getBeta(), beta);
      bl->getVerticalTrend(bl->getRho(), rho);
      bl->getVerticalTrend(seisLog, seisData);

      for (k = 0 ; k < nz ; k++) {
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;
      }

      int start,length;
      findContiniousPartOfData(hasData,nz,start,length);

      if(length*dz0 > waveletLength  ) // must have enough data
      {
        fillInCpp(alpha,beta,rho,start,length,cpp_r[w],nzp);
        fft(cpp_r[w],cpp_c[w],nzp);
        fillInSeismic(seisData,start,length,seis_r[w],nzp);
        fft(seis_r[w],seis_c[w],nzp);
        estimateCor(cpp_c[w],cpp_c[w],cor_cpp_c[w],cnzp);
        fftInv(cor_cpp_c[w],cor_cpp_r[w],nzp);
        estimateCor(seis_c[w],cpp_c[w],ccor_seis_cpp_c[w],cnzp);
        fftInv(ccor_seis_cpp_c[w],ccor_seis_cpp_r[w],nzp);
        fftInv(cpp_c[w],cpp_r[w],nzp);
        fftInv(seis_c[w],seis_r[w],nzp);
        wellWeight[w] = length*dz[w]*(cor_cpp_r[w][0]+cor_cpp_r[w][1]);// Gives most weight to long datasets with  
                                                                       // large reflection coefficients
      }
      else
        wellWeight[w] = 0; 
    }
  }
  delete [] seisLog;
  float* shiftWell = new float[nWells];
  float shiftAvg = shiftOptimal(ccor_seis_cpp_r,wellWeight,dz,nWells,nzp,shiftWell);

  multiplyPapolouis(ccor_seis_cpp_r,dz,nWells,nzp, waveletLength);
  multiplyPapolouis(cor_cpp_r,dz,nWells,nzp, waveletLength);
  getWavelet(ccor_seis_cpp_r,cor_cpp_r,wavelet_r,wellWeight, nWells, nzp);
  rAmp_      = averageWavelets(wavelet_r,nWells,nzp,wellWeight,dz,dz0); // wavelet centered
  cAmp_      = (fftw_complex*) rAmp_;
  char* fileName = new char[MAX_STRING];
  //sprintf(fileName,"wavelet");
  //printToFile(fileName,rAmp_, nzp);

  for(int w=0;w<nWells;w++) // gets syntetic seismic with estimated wavelet
  {
    fillInnWavelet(wavelet_r[w],nzp,dz[w]);
    shiftReal(shiftWell[w]/dz[w],wavelet_r[w],nzp);
    //sprintf(fileName,"waveletShift");
    //printToFile(fileName,wavelet_r[w], nzp);
    fft(wavelet_r[w],wavelet_c[w],nzp);
    //sprintf(fileName,"cpp");
    //printToFile(fileName,cpp_r[w], nzp);
    fft(cpp_r[w],cpp_c[w],nzp);
    convolve(wavelet_c[w],cpp_c[w],synt_seis_c[w],cnzp);
    fftInv(synt_seis_c[w],synt_seis_r[w],nzp); // 
    //sprintf(fileName,"syntSeis");
    //printToFile(fileName,synt_seis_r[w], nzp);
    //sprintf(fileName,"seis");
    //printToFile(fileName,seis_r[w], nzp);
    //printf("Test run\n");
  }

  float* scaleOptWell    = new float[nWells];
  float* errWellOptScale = new float[nWells];
  float* errWell = new float[nWells];
  float err;
  float scaleOpt = findOptimalWaveletScale(synt_seis_r,seis_r,nWells,nzp,wellWeight,err,errWell,scaleOptWell,errWellOptScale);

  delete [] scaleOptWell ;
  delete [] errWellOptScale;
  delete [] errWell ;
  
//  shiftReal(shiftAvg/dz0,rAmp_,nzp);
  
  shiftAndScale(shiftAvg,scaleOpt);//shifts wavelet average from wells
  invFFT1DInPlace();
  waveletLength_ = getWaveletLengthF();
  LogKit::writeLog("  Estimated wavelet length:  %.1fms\n",waveletLength_);

  if( LogKit::getDebugLevel() > 0 )
  {
    //flipUpDown();
    sprintf(fileName,"estimated_wavelet");
    float dzOut = 1.0; // sample at least as dense as this
    writeWaveletToFile(fileName, dzOut);
    //flipUpDown();
  }
  double norm2=0.0;
  for(i=0; i < nzp_; i++ )
      norm2 += rAmp_[i]*rAmp_[i];
  norm_= float( sqrt(norm2));
  //char* fileName=new char[MAX_STRING];


  if(LogKit::getDebugLevel() > 1) {
    sprintf(fileName,"estimated_wavelet_full_%d.txt",int(ceil((theta_*180/PI)-0.5)) );
    FILE* fid = fopen(fileName,"w");
    for(i=0;i<nzp_;i++)
      fprintf(fid,"%1.3e\n",rAmp_[i]);
    fclose(fid);
    
    for(j=0;j<nWells;j++)
    {
      sprintf(fileName,"seis_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",seis_r[0][i]);
      fclose(fid);
        
      sprintf(fileName,"cor_cpp_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",cor_cpp_r[0][i]);
      fclose(fid);
      
      sprintf(fileName,"ccor_seis_cpp_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",ccor_seis_cpp_r[0][i]);
      fclose(fid);

      sprintf(fileName,"cpp_%d_well_%d.txt",int(theta_*180/PI+0.5),j+1);
      fid = fopen(fileName,"w");
      for(i=0;i<nzp_;i++)
        fprintf(fid,"%1.3e\n",cpp_r[0][i]);
      fclose(fid);
    }
  }
  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;
  for(i=0;i<nWells;i++)
  {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
    delete [] cor_cpp_r[i] ;
    delete [] ccor_seis_cpp_r[i] ;
    delete [] wavelet_r[i];
  }
  //flipUpDown(); //NB ODD temporary fix
}

Wavelet::Wavelet(char * fileName, ModelSettings * modelSettings)
{
  maxShift_       = modelSettings->getMaxWaveletShift();
  minRelativeAmp_ = modelSettings->getMinRelWaveletAmp();

  isReal_     = true;
  inFFTorder_ = false;
  scale_=1; 
  gridNI_=0;   
  gridNJ_=0;
  shiftGrid_=NULL;  
  gainGrid_=NULL; 
  errCode_=0;
  FILE* file = fopen(fileName,"r");
 
  if(file == 0)
  {
    sprintf(errText_,"Error: Could not open file %s for reading.\n", fileName);
    errCode_=1;
    return;
  }
  else
    fclose(file);

  int fileFormat = getWaveletFileFormat(fileName);
  
  if(fileFormat < 0)
  {
    if (errCode_ == 0)  
    {
      sprintf(errText_,"Error: Unknown file format of file  %s.\n", fileName);
      errCode_=2;
      return;
    }
   }
   if(fileFormat==OLD)
      WaveletReadOld(fileName);
   if(fileFormat==JASON)
     WaveletReadJason(fileName);
   
   for(int i=0; i < rnzp_ ;i++)
    {  
      if(i < nzp_)
      {
        rAmp_[i]*=scale_;
      }
      else
      {
        rAmp_[i]=RMISSING;
      }// endif
    }//end for i 
  
}
void
Wavelet::WaveletReadJason(char * fileName)
{
  readtype_=JASON;
  FILE* file = fopen(fileName,"r");
  bool lineIsComment = true; 
  char* dummyStr= new char[MAX_STRING];
  while( lineIsComment ==true)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      readToEOL(file);
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
      return;
    } // endif
    else
    {
      readToEOL(file);
      if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
      {
        lineIsComment = false;
      }
    }
  }  // end while
  float shift=float(atof(dummyStr));
  if(fscanf(file,"%s",dummyStr) == EOF)
  {
    readToEOL(file);
    if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
    else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
    errCode_=2; 
    return;
  } // endif
  dz_=float(atof(dummyStr));
  if(fscanf(file,"%s",dummyStr) == EOF)
  {
    readToEOL(file);
    if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
    else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
    errCode_=2; 
    return;
  } // endif
  nz_=atoi(dummyStr);

  cz_ =  int(floor((fabs(shift/dz_))+0.5));
  nzp_         = nz_;
  cnzp_        = nzp_/2+1;
  rnzp_        = 2*cnzp_; 
  rAmp_        = (fftw_real* ) fftw_malloc(sizeof(float)*rnzp_);
  cAmp_        = (fftw_complex* ) rAmp_;
  norm_        = RMISSING;

  for(int i=0; i<nz_;i++)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      readToEOL(file);
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
      return;
    } // endif
    else
    {
      rAmp_[i] = (fftw_real) atof(dummyStr);
    }
  }
  fclose(file);
  waveletLength_ = getWaveletLengthF();
  LogKit::writeLog("\n  Estimated wavelet length:  %.1fms.",waveletLength_);
  delete [] dummyStr;
}
void
Wavelet::WaveletReadOld(char * fileName)
{
  readtype_=OLD;
  FILE* file = fopen(fileName,"r");
 
  int  maxWaveletL=10000;

  
  int   i,pos,shift,nSamples;
  char  headStr[MAX_STRING];
  char  tmpStr[MAX_STRING];
  char  targetString[MAX_STRING];
  char  number[MAX_STRING];
  float dz,ampMult;

  for(i = 0; i < 5; i++)
  {
    if(fscanf(file,"%s",headStr) == EOF)
    {
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
    } // endif
  }  // end for i


  strcpy(targetString,"CMX");
  pos = findEnd(headStr, 0, targetString);
  if(pos==-1) 
  {
    if (errCode_ == 0)  sprintf(errText_,"Error when reading wavelet amplitude from file  %s.\n", fileName);
    else sprintf(errText_,"%sError when reading wavelet amplitude from file  %s. \n", errText_,fileName);
    errCode_=3; 
    ampMult = RMISSING; // Dummy setting to avoid g++ warning 
  }
  else
  {
    readUntilStop(pos, headStr, number ,',');
    ampMult= float ( atof(number) );
  } //endif


  strcpy(targetString,"SI");
  pos = findEnd(headStr, 0, targetString);
  if(pos==-1) 
  {
    if (errCode_ == 0)  sprintf(errText_,"Error when reading sampling interval from file  %s.\n", fileName);
    else sprintf(errText_,"%sError when reading sampling interval from file  %s.\n",errText_,fileName);
    errCode_=4; 
    dz = RMISSING; // Dummy setting to avoid g++ warning
  }
  else
  {
    readUntilStop(pos, headStr, number ,',');
    dz= float ( atof(number) );
  } //endif

  float * tempWave= (float*) fftw_malloc(sizeof(float)* maxWaveletL );

  nSamples=0;

  while(fscanf(file,"%s",tmpStr) != EOF)
  {      
    if( maxWaveletL  > nSamples )
    {
      char * target = new char[1];
      strcpy(target,"F");
      pos = findEnd(tmpStr, 0, target);
      delete [] target;
      if(pos == -1)
      {
        tempWave[nSamples]=float( atof(tmpStr) );
      }
      else     
      {
        nSamples--;
      }
      nSamples++;
    }
    else
    {
      if (errCode_ == 0)  sprintf(errText_,"Error in memory use when reading wavelet from file  %s.\n", fileName);
      else  sprintf(errText_,"%sError in memory use when reading wavelet from file  %s.\n",errText_,fileName);
      errCode_=5;
    }//endif
  }//endwhile

  fclose(file);

  strcpy(targetString,"SHIFT");
  pos = findEnd(headStr, 0, targetString);
  if(pos==-1) 
  {
    shift=nSamples/2; // integer division
    if(shift*2 == nSamples)
    {
      if (errCode_ == 0)  
        sprintf(errText_,"Error when reading wavelet shift from file  %s.\n    --> No SHIFT and even number of data.\n", fileName);
      else  
        sprintf(errText_,"%sError when reading wavelet shift from file  %s.\n    --> No SHIFT and even number of data.\n",errText_,fileName); 
      errCode_=6; 
      //
      // NBNB-PAL: Temporary? hack since the errorCode is never returned from constructor
      //
      LogKit::writeLog("\nError when reading wavelet shift from file  %s.",fileName);
      LogKit::writeLog("\n    --> No SHIFT and even number of data.\n");
      exit(1);
    }
    //cz_ = shift;   // case no flip
    shift =-shift-1; // case flip
  }
  else
  {
    readUntilStop(pos, headStr, number ,',');
    shift= int ( atof(number) );

    //cz_=nSamples+shift;  // case no flip
  }//endif

  if(errCode_ == 0) 
  {
    cz_          = -shift-1; // case flip
    theta_       = RMISSING;
    dz_          = dz;
    nz_          = nSamples;
    nzp_         = nSamples;
    cnzp_        = nzp_/2+1;
    rnzp_        = 2*cnzp_; 
    inFFTorder_  = false;
    isReal_      = true;
    rAmp_        = (fftw_real* ) fftw_malloc(sizeof(float)*rnzp_);
    cAmp_        = (fftw_complex* ) rAmp_;
    norm_        = RMISSING;

    // Note the wavelet is fliped left to right  
    // since it is used for a convolution 

    for(i=0; i < rnzp_ ;i++)
    {  
      if(i < nzp_)
      {
        rAmp_[i]=ampMult*tempWave[nzp_-i-1];
      }
      else
      {
        rAmp_[i]=RMISSING;
      }// endif
    }//end for i 
  }

  fftw_free(tempWave);
  // LogKit::writeLog("\nReading wavelet file %s  ... done.\n",fileName);

  //
  // Estimate wavelet length
  //
  waveletLength_ = getWaveletLengthF();
  LogKit::writeLog("\n  Estimated wavelet length:  %.1fms.\n",waveletLength_);
}


Wavelet::Wavelet(Wavelet * wavelet)
{
  gridNI_    = 0;   
  gridNJ_    = 0;
  shiftGrid_ = NULL;  
  gainGrid_  = NULL; 

  theta_     = wavelet->theta_;
  dz_        = wavelet->dz_;
  nz_        = wavelet->nz_;
  nzp_       = wavelet->nzp_;
  cnzp_      = wavelet->cnzp_;
  rnzp_	     = wavelet->rnzp_;
  cz_        = wavelet->cz_;
  inFFTorder_= wavelet->inFFTorder_;
  isReal_    = wavelet->isReal_;  
  norm_      = wavelet->norm_;
  errCode_   = wavelet->errCode_;   
  if(errCode_ != 0)
    strcpy(errText_, wavelet->errText_) ;   
  rAmp_      = (fftw_real*) fftw_malloc(rnzp_*sizeof(fftw_real));  
  cAmp_      = (fftw_complex*) rAmp_;

  if(isReal_)
    for(int i = 0; i < rnzp_; i++)
    {
      rAmp_[i] = wavelet->rAmp_[i];  
    }
  else
    for(int i = 0; i < cnzp_; i++)
    {
      cAmp_[i].re = wavelet->cAmp_[i].re; 
      cAmp_[i].im = wavelet->cAmp_[i].im; 
    }
}

// Makes error correlations 
//
//
Wavelet::Wavelet(Wavelet * wavelet, int difftype)
{
  gridNI_    = 0;   
  gridNJ_    = 0;
  shiftGrid_ = NULL;  
  gainGrid_  = NULL; 
  if(! wavelet->getIsReal() ) wavelet->invFFT1DInPlace();
  theta_     = wavelet->theta_;
  dz_        = wavelet->dz_;
  nz_        = wavelet->nz_;
  nzp_       = wavelet->nzp_;
  cnzp_      = wavelet->cnzp_;
  rnzp_	     = wavelet->rnzp_;
  cz_        = wavelet->cz_;
  inFFTorder_= wavelet->inFFTorder_;
  isReal_    = wavelet->isReal_;  
  norm_      = wavelet->norm_;
  errCode_   = wavelet->errCode_;   
  if(errCode_ != 0)
    strcpy(errText_, wavelet->errText_) ;   
  rAmp_      = (fftw_real*) fftw_malloc(rnzp_*sizeof(fftw_real));  
  cAmp_      = (fftw_complex*) rAmp_;
  int i;

  double norm2 = 0.0;

  if(difftype != FOURIER)
  {
    for( i = 0; i < rnzp_; i++)
    {
      if(i < nzp_)
      {
        switch(difftype)
        {
        case FIRSTORDERFORWARDDIFF:
          if(i == nzp_-1 )
            rAmp_[i] = wavelet->rAmp_[0] - wavelet->rAmp_[i];
          else
            rAmp_[i] = wavelet->rAmp_[i+1]-wavelet->rAmp_[i];      
          break;
        case FIRSTORDERBACKWARDDIFF:
          if(i == 0 )
            rAmp_[i] = wavelet->rAmp_[0] - wavelet->rAmp_[nzp_-1];
          else
            rAmp_[i] = wavelet->rAmp_[i]-wavelet->rAmp_[i-1];         
          break;
        case FIRSTORDERCENTRALDIFF:
          if(i == 0 )
            rAmp_[i] = float( 0.5*(wavelet->rAmp_[i+1] - wavelet->rAmp_[nzp_-1]) );
          else
          {
            if(i == nzp_-1 )
              rAmp_[i] = float( 0.5*(wavelet->rAmp_[0] - wavelet->rAmp_[i-1]));
            else
              rAmp_[i] = float( 0.5*(wavelet->rAmp_[i+1] - wavelet->rAmp_[i-1]));
          }        
          break;
        }
        norm2 += rAmp_[i]*rAmp_[i];
      }
      else
        rAmp_[i] = RMISSING;
    }
    norm_=float( sqrt(norm2) );
  }
  else
  {
    fftw_complex  iValue;
    for(i=0; i < nzp_; i++ )
      rAmp_[i] = wavelet->rAmp_[i];
    fft1DInPlace();
    for(i=0; i < cnzp_; i++ )
    {
      iValue  =  cAmp_[i];
      cAmp_[i].re = float( - iValue.im * 2 * PI * float(i)/float(nzp_) );
      cAmp_[i].im = float(   iValue.re * 2 * PI * float(i)/float(nzp_) );
    }
    invFFT1DInPlace();
    for(i=0; i < nzp_; i++ )
      norm2 += rAmp_[i]*rAmp_[i];
    norm_= float( sqrt(norm2));
  }
}

Wavelet::Wavelet(int difftype, int nz, int nzp)
{
  gridNI_     = 0;   
  gridNJ_     = 0;
  shiftGrid_  = NULL;  
  gainGrid_   = NULL; 
  theta_      = RMISSING;
  dz_         = RMISSING;
  nz_         = nz;
  nzp_        = nzp;
  cnzp_       = nzp_/2+1;
  rnzp_	      = 2*cnzp_;
  cz_         = 0;
  inFFTorder_ = true;
  errCode_    = 0;
  rAmp_       = (fftw_real*) fftw_malloc(rnzp_*sizeof(fftw_real));  
  cAmp_       = (fftw_complex*) rAmp_;
  norm_       = RMISSING;
  int i;

  if(difftype != FOURIER)
  {
    isReal_    = true;

    for(i=0;i < rnzp_; i++)
    {
      if( i < nzp_)
        rAmp_[i] = 0.0;
      else
        rAmp_[i] = RMISSING;
    }

    switch(difftype)
    {
    case FIRSTORDERFORWARDDIFF:
      rAmp_[0] = -1.0; 
      rAmp_[1] = 1.0; 
      norm_    = float( sqrt(2.0) );
      break;
    case FIRSTORDERBACKWARDDIFF:
      rAmp_[0]      = 1.0;
      rAmp_[nzp_-1] =  -1.0;      
      norm_    = float( sqrt(2.0) ); 
      break;
    case FIRSTORDERCENTRALDIFF:
      rAmp_[1]      = 0.5;
      rAmp_[nzp_-1] = -0.5; 
      norm_    = float( sqrt(0.5) );
      break;
    }
  }
  else
  {
    double norm2 = 0.0;
    isReal_    = false;

    for(i=0;i < cnzp_; i++)
    {
      cAmp_[i].re = 0.0;
      cAmp_[i].im = float( 2.0 * PI * float( i ) / float( nzp_ ) );   
    }

    invFFT1DInPlace();

    for(i = 0; i < nzp_;i++) 
      norm2 +=  rAmp_[i]* rAmp_[i];

    norm_= float( sqrt( norm2 ) );
  }       
}


Wavelet::~Wavelet()
{
  fftw_free(rAmp_);

  if(shiftGrid_!=NULL)
    delete [] shiftGrid_;
  if(gainGrid_!=NULL)
    delete [] gainGrid_; 
}

int
Wavelet::getWaveletFileFormat(char * fileName)
{
  int fileformat=-1;
  char* dummyStr = new char[MAX_STRING];
  // test for old file format
  FILE* file = fopen(fileName,"r");
  for(int i = 0; i < 5; i++)
  {
    if(fscanf(file,"%s",dummyStr) == EOF)
    {
      if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
      else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
      errCode_=2; 
      return fileformat;
    } // endif
  }  // end for i
  fclose(file);
  char* targetString =new char[MAX_STRING];
  strcpy(targetString,"CMX");
  int  pos = findEnd(dummyStr, 0, targetString);
  if(pos>=0)
    fileformat= OLD;

  if(fileformat<0) // not old format
  {
    // test for jason file format
    file = fopen(fileName,"r");
    bool lineIsComment = true; 
    while( lineIsComment ==true)
    {
      if(fscanf(file,"%s",dummyStr) == EOF)
      {
        readToEOL(file);
        if (errCode_ == 0)  sprintf(errText_,"Error: End of file %s premature.\n", fileName);
        else sprintf(errText_,"%sError: End of file %s premature.\n", errText_,fileName);
        errCode_=2; 
        return fileformat;    
      } // endif
      else
      {
        readToEOL(file);
        if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
        {
          lineIsComment = false;
        }
      }
    }  // end while
    fclose(file);
    if(atof(dummyStr)!=0) // not convertable number
    {
      fileformat= JASON;
    }
  }

  delete[] dummyStr;
  delete [] targetString;
  return fileformat;
}
void
Wavelet::resample(float dz, int nz, float pz, float theta) 
{
  theta_=theta;
  if(errCode_==0){
    //LogKit::writeLog("  Resampling wavelet\n");
    assert(isReal_);
    assert(!inFFTorder_);
    int nzp,cnzp,rnzp,k;
    float z;
    fftw_real* wlet;

    nzp   =  findClosestFactorableNumber( (int) ceil( nz*(1.0f+pz) ) );
    cnzp  =  nzp/2 + 1;
    rnzp  =  2*cnzp;

    wlet  = (fftw_real *) fftw_malloc( sizeof(fftw_real)*rnzp );

    for(k=0; k < rnzp; k++)
    {
      if(k < nzp)
      {
        if(k < nzp/2+1)
        {
          z = float( dz*k );
        }
        else
        {
          z= float( dz*(k-nzp) );
        }
        wlet[k] = getWaveletValue(z, rAmp_ , cz_, nz_, dz_);
      }
      else
      {
        wlet[k] =RMISSING;
      }
    }
    fftw_free( rAmp_);

    float norm2 = 0.0; 
    for(k=0; k < nzp; k++) norm2 +=wlet[k]*wlet[k];

    rAmp_       = (fftw_real *)    wlet; // rAmp_ is not allocated 
    cAmp_       = (fftw_complex* ) rAmp_;
    nzp_        = nzp;
    rnzp_       = rnzp;
    cnzp_       = cnzp;
    cz_         = 0;
    nz_         = nz;
    dz_         = dz;
    norm_       = float( sqrt( norm2) );
    inFFTorder_ = true;
  }
  if(readtype_!=OLD)
     flipUpDown();
  if( LogKit::getDebugLevel() > 0 )
  {
    //flipUpDown();// ODD temporary debugfix
    char* fileName = new char[MAX_STRING];
    sprintf(fileName,"resampled_wavelet");
    float dzOut = 1.0; // sample at least as dense as this
    writeWaveletToFile(fileName, dzOut);
    //flipUpDown();// ODD temporary debugfix
    delete [] fileName;
  }
}

void           
Wavelet::flipUpDown()
{
  if(isReal_==true)
  {
    float tmp;
    for(int i=1;i<nzp_/2;i++)
    {
      tmp=rAmp_[i];
      rAmp_[i] = rAmp_[nzp_-i];
      rAmp_[nzp_-i] =tmp;
    }
  }
  else
  {
    for(int i=0;i<cnzp_;i++)
    {
      cAmp_[i].im *=-1.0;
    }
  }
}


Wavelet*  Wavelet::getLocalWavelet(int i, int j)
{
  Wavelet * localWavelet;
  localWavelet = new Wavelet(this);

  float  ts    = this->getLocalTimeshift(i,j);
  float gain   = this->getLocalGainFactor(i,j);

  localWavelet->shiftAndScale(ts,gain);

  return localWavelet;
}


float  Wavelet::getLocalTimeshift(int i, int j)
{
  float shift = 0.0f;

  int ind;
  if(shiftGrid_ != NULL && i < gridNI_ &&  i>=0 && j < gridNJ_ &&  j>=0)
  {
    ind = gridNI_*j + i;
    if(shiftGrid_[ind] != WELLMISSING)
      shift = shiftGrid_[ind];
  }


  return shift;
}

float  Wavelet::getLocalGainFactor(int i, int j)
{
  float gain = 1.0f;
  int ind;
  if(gainGrid_ != NULL && i < gridNI_ &&  i>=0 && j < gridNJ_ &&  j>=0)
  {
    ind = gridNI_*j + i;
    if(gainGrid_[ind] != WELLMISSING)
      gain = gainGrid_[ind];
  }

  return(gain);
}


void
Wavelet::shiftAndScale(float shift,float gain)
{
  int k;

  fftw_complex  ampMultiplier,tmp;

  if(isReal_) 
    fft1DInPlace();

  float iShift=shift/dz_;

  for(k=0;k < cnzp_; k++)
  {
    ampMultiplier.re = float(gain*cos(2.0*(PI*(iShift)*k)/float(nzp_)));
    ampMultiplier.im = float(gain*sin(-2.0*(PI*(iShift)*k)/float(nzp_)));

    tmp.re = ampMultiplier.re*cAmp_[k].re - ampMultiplier.im*cAmp_[k].im;
    tmp.im = ampMultiplier.im*cAmp_[k].re + ampMultiplier.re*cAmp_[k].im;

    cAmp_[k] =tmp;
  }
}

void
Wavelet::scale(float scale)
{
  if (scale != 1.0f)
    LogKit::writeLog("  Scaling wavelet with gain factor         : %.3e\n",scale);
  scale_ = scale;
  for(int i=0; i < rnzp_ ; i++)
    if(rAmp_[i] != RMISSING)
      rAmp_[i]=rAmp_[i]*scale;
  
  /*
    if(isReal_) 
    fft1DInPlace();
    
    for(int k=0 ; k < cnzp_ ; k++)
    {
    cAmp_[k].re = scale*cAmp_[k].re;
    cAmp_[k].im = scale*cAmp_[k].im;
    }
  */
}

float 
Wavelet::getWaveletValue(float z, float *Wavelet, int center, int nz, float dz)
{
  // returns the value of the vavelet in the location z. Wavelet have the length nz 
  // and the center value is Wavelet[center]
  // uses kriging with ricker 20Hz wavelet as correlation function.
  float value;
  // double a,b,c;
  //int    im1,i0,ip1,
  int    k,l;
  int*   ind=new int[6];// iL1,iL2,iL3,iR1,iR2,iR3;
  double* val=new double[6];//vL1,vL2,vL3,vR1,vR2,vR3;

  ind[2]= int( floor( (z/dz) ) );
  for(k=0;k<6;k++)
    ind[k]=  ind[2]+k-2;

  //i0    =  int( ceil( (z/dz) - 0.5) );
  //ip1   =  i0 + 1;
  //im1   =  i0 - 1;
  //hz    =  z - i0*dz;

  for(k=0;k<6;k++)
    val[k]=  getArrayValueOrZero(ind[k]+center , Wavelet,  nz); 

  double** Cov = new double*[6];
  double*  cov = new double[6];
  double   nu = 20;
  double   deltaT;

  for(k=0;k<6;k++)
  { 
    Cov[k] = new double[6];
    deltaT = (dz*ind[k]-z)*0.001;
    cov[k] = (1-2*nu*nu*PI*PI*(deltaT)*(deltaT))*exp(-nu*nu*PI*PI*(deltaT)*(deltaT));
    for(l=0;l<6;l++)
    {
      deltaT =(dz*ind[k]-dz*ind[l])*0.001;
      Cov[k][l] = (1-2*nu*nu*PI*PI*deltaT * deltaT )*exp(-nu*nu*PI*PI*deltaT*deltaT);
    }
  }
  //OK not very intellegent implementation since chol is done for each time step.
  lib_matrCholR(6,  Cov);
  lib_matrAxeqbR(6, Cov, cov); // cov contains now the kriging weigths; 

  value = 0.0;
  for(k=0;k<6;k++)
  {
    value+= float(val[k]*cov[k]);
    delete [] Cov[k];
  }
  delete [] Cov;
  delete [] cov;


  //v0    =  getArrayValueOrZero(i0+center  , Wavelet,  nz);    
  //vp1   =  getArrayValueOrZero(ip1+center , Wavelet,  nz);  
  //vm1   =  getArrayValueOrZero(im1+center , Wavelet,  nz);  

  //c     =  v0;
  //b     =  (vp1-vm1)/(2.0*dz);
  //a     =  (vp1+vm1-2.0*v0)/(2.0*dz*dz);
  //value2 =  float(a*hz*hz + b*hz + c);

  delete [] val;
  delete [] ind;
  return value;
}

float 
Wavelet::getArrayValueOrZero(int i  ,float * Wavelet, int nz)
{
  float value;

  if(i > -1 && i < nz)
    value = Wavelet[i];
  else
    value = 0.0; 
  return value;
}

void Wavelet::fft1DInPlace()
{
  // use the operator version of the fourier transform
  if(isReal_){
    int flag; 
    rfftwnd_plan plan;  
    flag    = FFTW_ESTIMATE | FFTW_IN_PLACE;
    plan    = rfftwnd_create_plan(1, &nzp_ ,FFTW_REAL_TO_COMPLEX,flag);
    rfftwnd_one_real_to_complex(plan,rAmp_,cAmp_);
    fftwnd_destroy_plan(plan);
    isReal_ = false;
  }
}

void Wavelet::invFFT1DInPlace()
{
  // use the operator version of the fourier transform
  if(!isReal_)
  {
    int flag;
    rfftwnd_plan plan; 

    flag = FFTW_ESTIMATE | FFTW_IN_PLACE;
    plan= rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);
    rfftwnd_one_complex_to_real(plan,cAmp_,rAmp_);
    fftwnd_destroy_plan(plan);
    isReal_=true;
    double scale= double(1.0/double(nzp_));
    for(int i=0; i < nzp_; i++)
      rAmp_[i] = (fftw_real) (rAmp_[i]*scale);  
  }
}


fftw_real  
Wavelet::getRAmp(int k)
{
  fftw_real value;

  if(isReal_)
  {
    if(k < nzp_)
      value = rAmp_[k];
    else
      value = 0.0;
  }
  else
  {
    invFFT1DInPlace();

    if(k < nzp_)
      value = rAmp_[k];
    else
      value = 0.0;

    fft1DInPlace();
  }
  return value;
}


fftw_complex   
Wavelet::getCAmp(int k)
{
  assert(!isReal_);
  fftw_complex  value;

  if(k < cnzp_)
  {
    value.re =  cAmp_[k].re;
    value.im = -cAmp_[k].im;
  }
  else
  {
    int refk =  nzp_-k;
    value.re =  cAmp_[refk].re;
    value.im =  cAmp_[refk].im;
  }
  return value;
}


fftw_complex   
Wavelet::getCAmp(int k, float scale)
{
  ///////////////////////////////////////////////////////////
  //
  // Get the  fourier transform of the streched wavelet.
  // scale is  in [0 1] wavelet
  // scale = 1 this is identical to getCAmp(int k)
  // Note we do not use the normal scale relation:  
  //  
  // FFT( w(s*t) ) = fw(w/s) * ( 1/s )   with  fw( w ) = FFT( w(t) )
  //
  // We return:  fw(w/s) 
  ///////////////////////////////////////////////////////////// 

  assert(!isReal_);

  fftw_complex  value;
  float omega,dOmega;  
  int   omU, omL;
  if( k < cnzp_)
  {
    omega  = float(k) / scale;
    if(omega >= cnzp_)
    {
      value.re =  0.0f;
      value.im =  0.0f;
    }
    else
    {
      omL        = int(floor( omega ));
      omU        = int( floor( omega )) + 1;
      if(omU >= cnzp_) 
        omU -= 1;    

      dOmega     = omega - float(omL);
      value.re =  cAmp_[omL].re * ( 1.0f - dOmega ) + cAmp_[omU].re * dOmega;
      value.im = -cAmp_[omL].im * ( 1.0f - dOmega ) - cAmp_[omU].im * dOmega;
      if(k < 0) //NBNB Ragnar
      {
        value.re = float(exp(log(cAmp_[omL].re) * ( 1.0f - dOmega ) + log(cAmp_[omU].re) * dOmega));
        value.re = float(exp(log(cAmp_[omL].im) * ( 1.0f - dOmega ) + log(cAmp_[omU].im) * dOmega));
      }
    }
  }
  else
  {
    int refk =  nzp_-k;
    omega = float(refk)/scale;

    if(omega >= cnzp_)
    {
      value.re =  0.0f;
      value.im =  0.0f;
    }
    else
    {
      omL        = int(floor( omega ));
      omU        = int( floor( omega )) + 1;
      if(omU >= cnzp_) 
        omU -= 1;
      dOmega     = omega - float(omL);
      value.re =  cAmp_[omL].re * ( 1.0f - dOmega ) + cAmp_[omU].re * dOmega;
      value.im =  cAmp_[omL].im * ( 1.0f - dOmega ) + cAmp_[omU].im * dOmega;
    }
  }
  return value;
}


void 
Wavelet::setShiftGrid(irapgrid * grid, Simbox * simbox)
{
  gridNI_ = simbox->getnx();
  gridNJ_ = simbox->getny();
  if(shiftGrid_ != NULL)
    delete [] shiftGrid_;
  shiftGrid_ = new float[gridNI_*gridNJ_];
  for(int j=0;j<gridNJ_;j++)
    for(int i=0;i<gridNI_;i++)
    {
      int outside;
      double x, y, z;
      simbox->getCoord(i, j, 0, x, y, z);
      shiftGrid_[i+gridNI_*j] = static_cast<float>(irapgridGetValue(x, y, grid, &outside));
    }
}

void 
Wavelet::setGainGrid(irapgrid * grid, Simbox * simbox)
{
  double sum = 0.0;
  int nData = 0;
  gridNI_ = simbox->getnx();
  gridNJ_ = simbox->getny();
  if(gainGrid_ != NULL)
    delete [] gainGrid_;
  gainGrid_ = new float[gridNI_*gridNJ_];
  for(int j=0;j<gridNJ_;j++)
    for(int i=0;i<gridNI_;i++)
    {
      int outside;
      double x, y, z;
      simbox->getCoord(i, j, 0, x, y, z);
      double value = irapgridGetValue(x, y, grid, &outside);
      gainGrid_[i+gridNI_*j] = static_cast<float>(value);
      if(value != WELLMISSING)
      {
        sum += log(value);
        nData++;
      }
    }
    float invGeoMean = float(exp(-sum/static_cast<double>(nData)));
    for(int j=0;j<gridNJ_;j++)
      for(int i=0;i<gridNI_;i++)
      {
        if(gainGrid_[i+gridNI_*j] != WELLMISSING)
          gainGrid_[i+gridNI_*j] *= invGeoMean;
      }
      float geoMean = 1/invGeoMean;
      norm_ *= geoMean;
      bool fftflag = isReal_;
      if(fftflag == true)
        invFFT1DInPlace();
      for(int i=0;i<nzp_;i++)
        rAmp_[i] *= geoMean;
      if(fftflag == true)
        fft1DInPlace();
}


void
Wavelet::printToFile(char* fileName, bool overrideDebug)
{
  if(overrideDebug == true || LogKit::getDebugLevel() > 0) {
      char * fName = LogKit::makeFullFileName(fileName, ".dat");
      FILE *file = fopen(fName,"wb");
      int i;
      for(i=0;i<nzp_;i++)
        fprintf(file,"%f\n",rAmp_[i]);

      fclose(file);
      delete [] fName;
    }
}

void
Wavelet::printToFile(char* fileName,fftw_real* vec, int nzp)
{
  if( LogKit::getDebugLevel() > 0) {
      char * fName = LogKit::makeFullFileName(fileName, ".dat");
      FILE *file = fopen(fName,"wb");
      int i;
      for(i=0;i<nzp;i++)
        fprintf(file,"%f\n",vec[i]);

      fclose(file);
      delete [] fName;
    }
}


void
Wavelet::writeWaveletToFile(char* fileName,float approxDzIn)
{
   flipUpDown(); //internal representation in CRAVA is UpDown (OKOK blame it on me OK)
   sprintf(fileName,"%s_%.1f_deg",fileName,theta_*180.0/PI);
  
   char * fName = LogKit::makeFullFileName(fileName, ".asc");
   LogKit::writeLog("  Writing Wavelet to file : \'%s\'\n", fName);
   int i;
   float approxDz;

   approxDz = MINIM(approxDzIn,floor(dz_*10)/10);
   approxDz = MINIM(approxDzIn,dz_);

   float T = nzp_*dz_;
   int nzpNew  = int(ceil(T/approxDz - 0.5));  
   float dznew   = T/float(nzpNew);
   int cnzpNew = (nzpNew/2)+1;

   fftw_real*     waveletNew_r = new fftw_real[2*cnzpNew];
   fftw_complex*  waveletNew_c =(fftw_complex*) waveletNew_r;
   fft1DInPlace();
   double multiplyer = double(nzpNew)/double(nzp_);

   for(i=0;i<cnzpNew;i++)
   {
     if(i < cnzp_)
     {
       waveletNew_c[i].re = (fftw_real) (cAmp_[i].re*multiplyer);
       waveletNew_c[i].im = (fftw_real) (cAmp_[i].im*multiplyer);
       if((i==(cnzp_-1)) & (2*((cnzp_-1)/2) != cnzp_-1)) //boundary effect in fft domain
         waveletNew_c[i].re*=0.5;
     }
     else
     { 
        waveletNew_c[i].re = 0;
        waveletNew_c[i].im = 0;
     }
   }
   invFFT1DInPlace();
   fftInv(waveletNew_c,waveletNew_r,nzpNew );// note might be n^2 algorithm for some nzpNew

   int wLength = int(floor(waveletLength_/dznew+0.5));
   int halfLength = wLength/2; // integer division
   wLength =  halfLength*2+1;// allways odd
   if( wLength>nzpNew)
   {  
     wLength=2*(nzpNew/2)-1;// allways odd
     halfLength=wLength/2;
   }

   float shift = -dznew*halfLength;
   FILE * file = fopen(fName,"wb");

   fprintf(file,"\"* Export format using Comma Separated Values\"\n");
   fprintf(file,"\"*Wavelet written from CRAVA\"\n");
   fprintf(file,"\"* Generated \"\n");
   fprintf(file,"\"*\"\n");
   fprintf(file,"\"* File format: \"\n");
   fprintf(file,"\"* - N lines starting with * are comment (such as this line)\"\n");
   fprintf(file,"\"* - 1 line with four fields (data type, data unit, depth type, depth unit) \"\n");
   fprintf(file,"\"* - 1 line with start time  \"\n");
   fprintf(file,"\"* - 1 line with sample interval \"\n");
   fprintf(file,"\"* - 1 line with number of data lines \"\n");
   fprintf(file,"\"* - N lines with trace data \"\n");
   fprintf(file,"\"* Data values are represented as floating point numbers,\"\n");
   fprintf(file,"\"* \"\n");
   fprintf(file,"\"wavelet\",\"none\",\"time\",\"ms\"\n");
   fprintf(file,"%1.0f\n",shift);
   fprintf(file,"%1.2f\n",dznew);
   fprintf(file,"%d\n", wLength);  

   for(i=halfLength;i > 0;i--)
     fprintf(file,"%f\n", waveletNew_r[nzpNew-i]);
   for(i=0;i<=halfLength;i++)
     fprintf(file,"%f\n", waveletNew_r[i]);
   
   fclose(file);
   flipUpDown();//internal representation in CRAVA is UpDown (OKOK blame it on me OK)
   delete [] fName;
   delete [] waveletNew_r;
}

void           
Wavelet::fftInv(fftw_complex* cAmp,fftw_real* rAmp,int nt)
{
  rfftwnd_plan p2 = rfftwnd_create_plan(1, &nt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_one_complex_to_real(p2, cAmp, rAmp);
  fftwnd_destroy_plan(p2);
  double sf = 1.0/double(nt);
  for(int i=0;i<nt;i++)
    rAmp[i]*=fftw_real(sf);
}

void
Wavelet::fft(fftw_real* rAmp,fftw_complex* cAmp,int nt)
{
  rfftwnd_plan p1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_one_real_to_complex(p1, rAmp, cAmp);
  fftwnd_destroy_plan(p1);
}

void 
Wavelet::shiftReal(float shift, fftw_real* rAmp,int nt)
{
  fftw_complex* cAmp = (fftw_complex*) rAmp;
  fft(rAmp,cAmp, nt);
  int cnzp= nt/2+1;
  float expo;
  fftw_complex tmp,mult;
  for(int i=0;i<cnzp;i++)
  {
    tmp     = cAmp[i];
    expo    = float(-2.0*shift*PI*float(i)/float(nt));
    mult.re = cos(expo);
    mult.im = sin(expo);
    cAmp[i].re = tmp.re*mult.re-tmp.im*mult.im;
    cAmp[i].im = tmp.re*mult.im+tmp.im*mult.re;
  }
  fftInv(cAmp,rAmp, nt);
}

void 
Wavelet::shiftReal(int shift, fftw_real* rAmp,int nt)
{
  float* tmp=new float[nt];
  int i,index;
  for(i=0;i<nt;i++)
    tmp[i]=rAmp[i];
  
  for(i=0;i<nt;i++)
  {
    index = i-shift;
    if(index<0)
      index+=nt;
    if(index>=nt)
      index-=nt;

    rAmp[i] = tmp[index];
  }
  delete [] tmp;
}

void
Wavelet::fillInCpp(float* alpha,float* beta,float* rho,int start,int length,fftw_real* cpp_r,int nzp)
{
  int i;

  for(i=0;i<nzp;i++)
    cpp_r[i]=0;

  for(i=start;i < start+length-1;i++)
  {
    float ei1 = computeElasticImpedance(alpha[i],beta[i],rho[i]);
    float ei2 = computeElasticImpedance(alpha[i+1],beta[i+1],rho[i+1]);
    cpp_r[i] =  ei2-ei1;
  } 
}


void
Wavelet::fillInSeismic(float* seisData,int start, int length,fftw_real* seis_r,int nzp)
{ 
  int i;
  for(i=0; i<nzp; i++)
    seis_r[i] = 0.0;

  for(i=start; i<start+length; i++)
  {
    seis_r[i] = seisData[i];
  }
/*
  int lTregion = 3;
  int* modify  = getIndexPrior(start,lTregion,nzp);
  int* conditionto = getIndexPost(start-1,lTregion,nzp);
  //NBNB Odd: interpolate endpoints?
*/

}

void
Wavelet::estimateCor(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* ccor_1_2_c,int cnzp)
{
  for(int i=0;i<cnzp;i++)
  {
    ccor_1_2_c[i].re = var1_c[i].re*var2_c[i].re+var1_c[i].im*var2_c[i].im;
    ccor_1_2_c[i].im = -var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }

}

void
Wavelet::convolve(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* out_c,int cnzp)
{
  for(int i=0;i<cnzp;i++)
  {
    out_c[i].re = var1_c[i].re*var2_c[i].re-var1_c[i].im*var2_c[i].im;
    out_c[i].im = var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }

}

float
Wavelet::computeElasticImpedance(float vp, float vs, float rho)
{
  // vp, vs, rho are logtransformed
  float angImp;

  angImp = float(coeff_[0]*vp+coeff_[1]*vs+coeff_[2]*rho );
  
  return(angImp); 
}


void
Wavelet::findContiniousPartOfData(bool* hasData,int nz,int &start, int &length)
{ 
  int i;
  int lPice=0;
  int lengthMaxPice=-1;
  int startLongestPice=0;
  bool previousHadData = false;

  for(i = 0; i < nz ;i++)
  {
    if(hasData[i])
    {
      if(! previousHadData)
        lPice=1;
      else
        lPice++;
      previousHadData = true;
    }
    else
    {
      if(previousHadData)
      {
        if(lengthMaxPice < lPice)
        {
          lengthMaxPice  = lPice;
          startLongestPice = i-lPice;
        }
      }
      previousHadData=false;
    }
  }

  if(previousHadData)
  {
    if(lengthMaxPice < lPice)
    {
      lengthMaxPice  = lPice;
      startLongestPice = i-lPice;
    }
  }

  start  = startLongestPice;
  length = lengthMaxPice; 
}

int* 
Wavelet::getIndexPrior(int start,int nInd,int nzp)
{
   //returns indexbefore start np cyclic;
  int* index= new int[nInd];
  int i,tmp;
  for(i=0;i<nInd;i++)
  {
    tmp=start-i-1;
    if(tmp<0)
      tmp = tmp+nzp;
    if(tmp>=nzp)
      tmp = tmp-nzp;
    index[i]= tmp;
  }
  return index;
}

int* 
Wavelet::getIndexPost(int start,int nInd,int nzp)
{
  //returns index after start np cyclic;
  int* index = new int[nInd];
  int i,tmp;
  for(i=0;i<nInd;i++)
  {
    tmp=start+i+1;
    if(tmp<0)
      tmp = tmp+nzp;
    if(tmp>=nzp)
      tmp = tmp-nzp;
    index[i]= tmp;
  }
  return index;
}

float
Wavelet::shiftOptimal(fftw_real** ccor_seis_cpp_r,float* wellWeight,float* dz,int nWells,int nzp,float* shiftWell)
{
  float shift=0.0f;
  float totalWeight=0;
  float sum=0;
  int w,i,polarity;
  // if the sum from -maxShift_ to maxShift_ ms is 
  // positive then polarity is positive
  for(w=0;w<nWells;w++)
  {
    shiftWell[w]=0;
    if(wellWeight[w] > 0) // only include those estimated
    {
      for(i=0;i<ceil(maxShift_/dz[w]);i++)//zero included
        sum+=ccor_seis_cpp_r[w][i];
      for(i=0;i<floor(maxShift_/dz[w]);i++)
        sum+=ccor_seis_cpp_r[w][nzp-i-1];
    }
  }

  polarity=-1;
  if(sum > 0)
    polarity=1;

  // gets optimal shift
  float maxValue;
  float shiftF;
  int shiftI;
  float f1,f2,f3;

  for(w=0;w<nWells;w++)
  {
    if(wellWeight[w]>0)
    {
      maxValue = 0.0f;
      shiftI=0;

      for(i=0;i<ceil(maxShift_/dz[w]);i++)
      {
        if(ccor_seis_cpp_r[w][i]*polarity > maxValue)
        {
          maxValue = ccor_seis_cpp_r[w][i]*polarity;
          shiftI = i;
        }
      }
      for(i=0;i<floor(maxShift_/dz[w]);i++)
      {
        if(ccor_seis_cpp_r[w][nzp-1-i]*polarity > maxValue)
        {
          maxValue = ccor_seis_cpp_r[w][nzp-1-i]*polarity;
          shiftI = -1-i;
        }
      }
      if(shiftI < 0)
      {
        if(ccor_seis_cpp_r[w][nzp+shiftI-1]*polarity < maxValue) //then local max
        {
          f1 = ccor_seis_cpp_r[w][nzp+shiftI-1];
          f2 = ccor_seis_cpp_r[w][nzp+shiftI];
          int ind3;
          if(shiftI==-1)
            ind3 = 0;
          else
            ind3=nzp+shiftI+1;
          f3 = ccor_seis_cpp_r[w][ind3];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI);
      }
      else //positive or zero shift
      {
        if(ccor_seis_cpp_r[w][shiftI+1]*polarity < maxValue) //then local max
        {
          f3 = ccor_seis_cpp_r[w][shiftI+1];
          f2 = ccor_seis_cpp_r[w][shiftI];
          int ind1;
          if(shiftI==0)
            ind1 = nzp-1;
          else
            ind1=shiftI-1;
          f1 = ccor_seis_cpp_r[w][ind1];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI);
      }
      shiftWell[w] = shiftF*dz[w];
      shiftReal(-shiftF, ccor_seis_cpp_r[w],nzp);// 
      shift += wellWeight[w]*shiftF*dz[w];//weigthing shift according to wellWeight
      totalWeight += wellWeight[w];
    }
  }
  shift/=totalWeight;
  return shift;
}

void
Wavelet::multiplyPapolouis(fftw_real** vec, float* dz,int nWells,int nzp, float waveletLength)
{
  int i,w;
  float wHL=float( waveletLength/2.0);
  float weight,dist;
  for(w=0;w<nWells;w++)
    for(i=1;i<nzp;i++)
    {
      dist =MINIM(i,nzp-i)*dz[w];
      if(dist < wHL)
      {
        weight  = float(1.0/PI*fabs(sin(PI*(dist)/wHL))); 
        weight += float((1-fabs(dist)/wHL)*cos(PI*dist/wHL));
      }
      else
        weight=0;

      vec[w][i]*=weight;
    }
}

void
Wavelet::getWavelet(fftw_real** ccor_seis_cpp_r,fftw_real** cor_cpp_r,fftw_real** wavelet_r,float* wellWeight,int nWells,int nt)
{
  fftw_complex* c_sc,*c_cc,*wav;

  int cnzp = nt/2+1;
  for(int w=0;w<nWells;w++)
  {
    if(wellWeight[w] > 0)
    {
      c_sc   = (fftw_complex*)ccor_seis_cpp_r[w];
      fft(ccor_seis_cpp_r[w],c_sc,nt);
      c_cc   = (fftw_complex*)cor_cpp_r[w];
      fft(cor_cpp_r[w],c_cc,nt);
      wav    = (fftw_complex*) wavelet_r[w];

      for(int i=0;i<cnzp;i++)
      {
        wav[i].re=c_sc[i].re/c_cc[i].re;//note c_cc[i].im =0
        wav[i].im=c_sc[i].im/c_cc[i].re;
      }
      fftInv(c_sc,ccor_seis_cpp_r[w],nt);
      fftInv(c_cc,cor_cpp_r[w],nt);
      fftInv(wav,wavelet_r[w],nt);
    }
  }
}

fftw_real* 
Wavelet::averageWavelets(fftw_real** wavelet_r,int nWells,int nzp,float* wellWeight,float* dz,float dzOut)
{
  // assumes dz[w] < dzOut for all w
  fftw_real* wave= (fftw_real*) fftw_malloc(rnzp_*sizeof(fftw_real));  

  int w,i;

  float* weight= new float[nWells];// weight is length of data interval in well
  float sum = 0;
  for(i=0;i<nWells;i++)
    sum+=wellWeight[i];

  for(i=0;i<nWells;i++)
    weight[i]=wellWeight[i]/sum;

  
  for(i=0;i<nzp;i++)
    wave[i] = 0.0; // initialize
  float ww;

  for( w=0;w<nWells;w++)
  {
    if(wellWeight[w]>0)
    {
      wave[0] += weight[w]* wavelet_r[w][0]; // central

      for(i=1;i <nzp/2-1;i++) // positive time
      {
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];  // fraction of distance to ind*dz[w]
        if(ind<nzp/2)
            ww= (1-t)* wavelet_r[w][ind]+t*(wavelet_r[w][ind+1]);
        else
          ww=0;
        wave[i]    += weight[w]*ww;

      }
      for(i=1;i < nzp/2-1;i++)// negative time Note dz[w] <= dzOut for all w
      {
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];
        if(ind<nzp/2)
          ww    = (1-t)* wavelet_r[w][nzp-ind]+t*(wavelet_r[w][nzp-(ind+1)]);
        else
          ww=0;

        wave[nzp-i] += weight[w]*ww;
      }
    }
  }
  
  char* fileName = new char[MAX_STRING];
  sprintf(fileName,"wavelet_%d_fftOrder_noshift",int(floor(theta_/PI*180+0.5)));
  printToFile(fileName,wave,nzp_);
  
  delete [] weight;
  return wave;
}
float         
Wavelet::getNoiseStandardDeviation(Simbox * simbox, FFTGrid * seisCube, WellData ** wells, int nWells)
{
  //LogKit::writeLog("\n  Estimating noise from seismic data and (nonfiltered) blocked wells");

  float errStd  = 0.0f;
  float dataVar = 0.0f;
  // initialization
  scale_=1; 
  gridNI_=0;   
  gridNJ_=0;
  shiftGrid_=NULL;  
  gainGrid_=NULL; 
  errCode_=0;
  
  float * dz = new float[nWells];
  int * nDataUsedInWell = new int[nWells];
  int i, k, w;

  int   nz   = simbox->getnz();
  float dz0  = static_cast<float>(simbox->getdz());
  int   nzp  = seisCube->getNzp();
  int   cnzp = (nzp/2+1);
  int   rnzp = 2*cnzp;

  float * alpha    = new float[nz];
  float * beta     = new float[nz];
  float * rho      = new float[nz];
  float * seisData = new float[nz];
  bool  * hasData  = new bool[nz];

  //Noise estimation
  fftw_real**       cpp_r = new fftw_real*[nWells];
  fftw_complex**    cpp_c = (fftw_complex** ) cpp_r;
  
  fftw_real**       seis_r = new fftw_real*[nWells];
  fftw_complex**    seis_c = (fftw_complex** ) seis_r; 

  fftw_real**       synt_r = new fftw_real*[nWells];
  fftw_complex**    synt_c = (fftw_complex** ) synt_r;

  fftw_real**       cor_seis_synt_r = new fftw_real*[nWells];
  fftw_complex**    cor_seis_synt_c = (fftw_complex** ) cor_seis_synt_r; 

  //fftw_real**       err_r = new fftw_real*[nWells];
//  fftw_complex**    err_c = (fftw_complex** ) err_r;  //Useful for debug

  fftw_real**      wavelet_r = new fftw_real*[nWells];
  fftw_complex**   wavelet_c = (fftw_complex**) wavelet_r;
  float* shiftWell   = new float[nWells];
  float* errVarWell  = new float[nWells];
  float* dataVarWell = new float[nWells];
  int*   nActiveData = new int[nWells];

  int maxBlocks = 0;
  //float wShift  = getWaveletShift();
  for(i=0;i<nWells;i++)
  {
    cpp_r[i]     = new fftw_real[rnzp];
    seis_r[i]    = new fftw_real[rnzp];
    cor_seis_synt_r[i]     = new fftw_real[rnzp];
    wavelet_r[i] = new fftw_real[rnzp];
    synt_r[i]    = new fftw_real[rnzp];
    nDataUsedInWell[i] = 0;
    const int * ipos = wells[i]->getBlockedLogsPropThick()->getIpos();
    const int * jpos = wells[i]->getBlockedLogsPropThick()->getJpos();
    dz[i] = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])*simbox->getdz());
    errVarWell[i]  = 0.0f;
    dataVarWell[i] =0.0;
    nActiveData[i] = 0;
    int nBlocks = wells[i]->getBlockedLogsPropThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (w = 0 ; w < nWells ; w++) 
  {
    if (!wells[w]->isDeviated())
    {
      BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
      //
      // Block seismic data for this well
      //
      bl->getBlockedGrid(seisCube,seisLog);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bl->getVerticalTrend(bl->getAlpha(), alpha);
      bl->getVerticalTrend(bl->getBeta(), beta);
      bl->getVerticalTrend(bl->getRho(), rho);
      bl->getVerticalTrend(seisLog, seisData);

      for (k = 0 ; k < nz ; k++) {
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;
      }

      int start,length;
      findContiniousPartOfData(hasData,nz,start,length);

      if(length*dz0 > waveletLength_) // must have enough data
      {
        nDataUsedInWell[w] = length; 
        fillInCpp(alpha,beta,rho,start,length,cpp_r[w],nzp);
        fft(cpp_r[w],cpp_c[w],nzp);
        fillInnWavelet(wavelet_r[w],nzp,dz[w]);
        fft(wavelet_r[w],wavelet_c[w],nzp);
        convolve(cpp_c[w],wavelet_c[w],synt_c[w],cnzp);
        fillInSeismic(seisData,start, length,seis_r[w],nzp);
        fft(seis_r[w],seis_c[w],nzp);
        estimateCor(synt_c[w],seis_c[w],cor_seis_synt_c[w],cnzp);
        fftInv(cor_seis_synt_c[w],cor_seis_synt_r[w],nzp);
        float shift=findBulkShift(cor_seis_synt_r[w],dz[w], nzp);
        shift = floor(shift*10.0f+0.5f)/10.0f;//rounds to nearest 0.1 ms (don't have more accuracy)
        shiftWell[w]=shift;
        fftInv(synt_c[w],synt_r[w],nzp);
        shiftReal(-shift/dz[w],synt_r[w],nzp);
        fillInSeismic(seisData,start, length,seis_r[w],nzp);
        if(LogKit::getDebugLevel() > 0)
        {
          char* fileName = new char[MAX_STRING];
          sprintf(fileName,"seismic_Well_%d_%d",w,int(floor(theta_/PI*180.0+0.5)));
          printToFile(fileName,seis_r[w], nzp);
          sprintf(fileName,"synthetic_seismic_Well_%d_%d",w,int(floor(theta_/PI*180.0+0.5)));
          printToFile(fileName,synt_r[w], nzp);
          sprintf(fileName,"wellTime_Well_%d",w);
          char * fName = LogKit::makeFullFileName(fileName, ".dat");
          FILE *file = fopen(fName,"wb");
          //fprint(file,"%f %f  %f \n",,dz[w],);
          fclose(file);
          delete [] fileName;
          delete fName;
        }
        for(i=start;i<start+length;i++)
        { 
          float err=(seis_r[w][i] - synt_r[w][i]);
          errVarWell[w]+=err*err;
          dataVarWell[w]+=seis_r[w][i] *seis_r[w][i] ;
        }
        nActiveData[w]=length;
      }
      else
      {
        nDataUsedInWell[w] = 0; 
        LogKit::writeLog("\n  Not using vertical well %s for error estimation (length=%.1fms  required length=%.1fms).",
                         wells[w]->getWellname(),length*dz0,waveletLength_);
      }
    }
  }
  float  errOptScale;
  float* errWell=new float[nWells];
  float* scaleOptWell=new float[nWells];
  float* errWellOptScale=new float[nWells];

  float optScale = findOptimalWaveletScale(synt_r,seis_r,nWells,nzp,dataVarWell,
                                           errOptScale,errWell,scaleOptWell,errWellOptScale);
  delete [] seisLog;

  int nData=0;
  for(i=0;i<nWells;i++)
  {
    nData   += nActiveData[i];
    errStd  += errVarWell[i];
    dataVar += dataVarWell[i];
    if(nActiveData[i]>0)
    {    
      errVarWell[i]  /= nActiveData[i];
      dataVarWell[i] /= nActiveData[i];
    }
  }

  if (nData == 0)
  {
    LogKit::writeLog("\nERROR: Cannot estimate signal-to-noiseratio. No legal well data available. CRAVA must stop\n");
    exit(1);
  }
  dataVar /= float(nData);
  errStd  /= float(nData);
  errStd   = sqrt(errStd);
  
  //LogKit::writeLog("\n  Reporting errors (as standard deviations) estimated in different ways:\n\n");
  LogKit::writeLog("                                     SeisData        ActuallyUsed       OptimalGlobal      OptimalLocal\n");
  LogKit::writeLog("  Well                  shift[ms]     StdDev          Gain   S/N         Gain   S/N         Gain   S/N \n");
  LogKit::writeLog("  --------------------------------------------------------------------------------------------------------\n");
  for(i=0;i<nWells;i++)
  {
    if(nActiveData[i]>0) {
      float SNActuallyUsed  = dataVarWell[i]/errVarWell[i];
      float SNOptimalGlobal = dataVarWell[i]/(errWell[i]*errWell[i]);
      float SNOptimalLocal  = dataVarWell[i]/(errWellOptScale[i]*errWellOptScale[i]);
      LogKit::writeLog("  %-20s   %6.2f     %9.2e        1.00 %7.2f     %6.2f %7.2f     %6.2f %7.2f\n", 
                       wells[i]->getWellname(),shiftWell[i],sqrt(dataVarWell[i]),
                       SNActuallyUsed,optScale,SNOptimalGlobal,scaleOptWell[i],SNOptimalLocal);
    }
    else
      LogKit::writeLog("  %-20s      -            -             -      -           -      -           -      -   \n",
                       wells[i]->getWellname()); 
  }
  float empSNRatio = dataVar/(errStd*errStd);
  LogKit::writeLog("\n  Signal to noise ratio used for this angle stack is: %6.2f\n", empSNRatio);

  if (empSNRatio < 1.1f) 
  {
    LogKit::writeLog("\nERROR: The empirical signal-to-noise ratio Var(seismic data)/Var(noise) is %.2f. Ratios\n",empSNRatio);
    LogKit::writeLog("       smaller than 1 are illegal and CRAVA has to stop. CRAVA was for some reason not able\n");
    LogKit::writeLog("       to estimate this ratio reliably, and you must give it as input to the model file.\n");
    LogKit::writeLog("\nNOTE: If the wavelet was estimated by CRAVA the solution may be to remove one or more wells\n");
    LogKit::writeLog("      from the wavelet estimation (compare shifts and SN-ratios for different wells).\n\n");
    exit(1);
  }
 
  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;
  for(i=0;i<nWells;i++)
  {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
    //delete [] err_r[i] ;
    delete [] synt_r[i] ;
    delete [] wavelet_r[i];
  }
  delete [] cpp_r;
  delete [] seis_r;
  //delete [] err_r;
  delete [] synt_r;
  delete [] wavelet_r;
      
  return errStd;
}

void
Wavelet::fillInnWavelet(fftw_real* wavelet_r,int nzp,float dz)
{
  //Note: Assumes that dz_ > dz
  // NBNB OddK 
  fftw_real previous1 = rAmp_[0];
  fftw_real previous2 = rAmp_[0];
  fftw_real current1 = 0.0f;
  fftw_real current2 = 0.0f;
  wavelet_r[0]        = previous1;
  //wavelet_r[nzp/2]    = 0;
  //wavelet_r[nzp/2+1]    = 0;
  int counterForWavelet   = 1;

  float w;
  int i;
  for(i=1;i<= nzp/2;i++)
  {
    if(counterForWavelet*dz_ < i*dz)
    {
      counterForWavelet++;
      previous1=current1;
      previous2=current2;
    }
    current1=rAmp_[counterForWavelet];
    current2=rAmp_[nzp-counterForWavelet];
    w = (counterForWavelet*dz_-i*dz) /dz_;
    wavelet_r[i]     = (1-w)*current1+ w*previous1;
    wavelet_r[nzp-i] = (1-w)*current2+ w*previous2;
  }
}


float          
Wavelet::findOptimalWaveletScale(fftw_real** synt_seis_r,fftw_real** seis_r,int nWells,int nzp,
                                 float* wellWeight,float& err,float* errWell,float* scaleOptWell,float* errWellOptScale)
{
  float optScale=1.0;
  int    nScales = 51; // should be odd to include 1.00
  float* scales  = new float[nScales];
  float* error  = new float[nScales];
  float  scaleLimit=3.0;

  for(int i=0;i<nScales;i++)
  {
    scales[i] = exp(-log(scaleLimit)+i*2*(log(scaleLimit))/(nScales-1));
    error[i]=0.0;
  }

  int* counter     = new int[nWells];
  float* seisNorm  = new float[nWells];
  float** resNorm  = new float*[nWells];
  for(int i=0;i<nWells;i++)
    resNorm[i] =   new float[nScales];

  for(int i=0;i<nWells;i++)
  {
    if(wellWeight[i]>0)
    {
      for(int j=0;j<nScales;j++)
      {
        resNorm[i][j]=0.0;
        counter[i]=0;
        for(int k=0;k<nzp;k++)
        {
          if( fabs(seis_r[i][k])> 1e-7)
          {
            seisNorm[i]  +=   seis_r[i][k]* seis_r[i][k];
            float     foo = scales[j]*synt_seis_r[i][k] - seis_r[i][k];
            resNorm[i][j] += foo*foo;
            counter[i]++;
          }
        }
      }
    }//if
  }
  int totCount=0;
  for(int i=0;i<nWells;i++)
  {
    for(int j=0;j<nScales;j++)
      error[j]+=resNorm[i][j];
    totCount+=counter[i];
  }

  int optInd=0; 
  float optValue=error[0];
  for(int i=1;i<nScales;i++)
    if(error[i]<optValue)
    {
      optValue=error[i];
      optInd=i;
    }

  err = sqrt(optValue/float(totCount));
  optScale = scales[optInd];
  for(int i=0;i<nWells;i++)
    errWell[i]=sqrt(resNorm[i][optInd]/counter[i]);

  for(int i=0;i<nWells;i++)
  {
    optValue = resNorm[i][0];
    optInd=0;
    for(int j=1;j<nScales;j++)
      if(resNorm[i][j]<optValue)
      {
        optValue=resNorm[i][j];
        optInd=j;
      }
    scaleOptWell[i]    = scales[optInd];
    errWellOptScale[i] = sqrt(optValue/float(counter[i]));
  }

 return optScale;
}


int Wavelet::getWaveletLengthI()
{
  bool trans=false;
  if(isReal_==false)
  {
    invFFT1DInPlace();
    trans=true;
  }
  float maxAmp =  fabs(rAmp_[0]); // gets max amp 
  for(int i=1;i <nzp_;i++)
    if(fabs(rAmp_[i]) > maxAmp)
      maxAmp = fabs(rAmp_[i]);

  float minAmp= maxAmp*minRelativeAmp_; // minimum relevant amplitude

  int wLength=nzp_;

  for(int i=nzp_/2;i>0;i--)
  {
    if(fabs(rAmp_[i]) >minAmp)
    {
      wLength= (i*2+1);// adds both sides 
      break;
    }
    if(fabs(rAmp_[nzp_-i]) > minAmp)
    {
      wLength= (2*i+1);// adds both sides 
      break;
    }
  }
  wLength =MINIM(wLength,2*(nzp_/2)-1); // allways even number
  if(trans==true)
    fft1DInPlace();

  return wLength;
}
float
Wavelet::getWaveletLengthF()
{
  return dz_*float( getWaveletLengthI() );
}

float
Wavelet::findBulkShift(fftw_real* vec_r,float dz,int nzp)
{
  
  
  float shift=0.0f;
  
  float sum=0;
  int i,polarity;
  // if the sum from -maxShift_ to maxShift_ ms is 
  // positive then polarity is positive   
  for(i=0;i<ceil(maxShift_/dz);i++)//zero included
    sum+=vec_r[i];
  for(i=0;i<floor(maxShift_/dz);i++)
    sum+=vec_r[nzp-i-1];
    
  polarity=-1;
  if(sum > 0)
    polarity=1;

  // gets optimal shift
  float maxValue;
  float shiftF;
  int shiftI;
  float f1,f2,f3;

 
  maxValue = 0.0f;
  shiftI=0;

  for(i=0;i<ceil(maxShift_/dz);i++)
  {
    if(vec_r[i]*polarity > maxValue)
    {
      maxValue = vec_r[i]*polarity;
      shiftI = i;
    }
  }
  for(i=0;i<floor(maxShift_/dz);i++)
  {
    if(vec_r[nzp-1-i]*polarity > maxValue)
    {
      maxValue = vec_r[nzp-1-i]*polarity;
      shiftI = -1-i;
    }
  }
  if(shiftI < 0)
  {
    if(vec_r[nzp+shiftI-1]*polarity < maxValue) //then local max
    {
      f1 = vec_r[nzp+shiftI-1];
      f2 = vec_r[nzp+shiftI];
      int ind3;
      if(shiftI==-1)
        ind3 = 0;
      else
        ind3=nzp+shiftI+1;
      f3 = vec_r[ind3];
      float x0=(f1-f3)/(2*(f1+f3-2*f2));
      shiftF=float(shiftI)+x0;
    }
    else  // do as good as we can
      shiftF=float(shiftI);
  }
  else //positive or zero shift
  {
    if(vec_r[shiftI+1]*polarity < maxValue) //then local max
    {
      f3 = vec_r[shiftI+1];
      f2 = vec_r[shiftI];
      int ind1;
      if(shiftI==0)
        ind1 = nzp-1;
      else
        ind1=shiftI-1;
      f1 = vec_r[ind1];
      float x0=(f1-f3)/(2*(f1+f3-2*f2));
      shiftF=shiftI+x0;
    }
    else  // do as good as we can
      shiftF=float(shiftI);
  }
  shift = shiftF*dz;
  shiftReal(-shiftF,vec_r,nzp);// for testing
 
  return shift;
}



