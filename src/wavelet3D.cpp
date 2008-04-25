#include "wavelet3D.h"

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
#include "lib/sgri.h"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

Wavelet3D::Wavelet3D(char * fileName, ModelSettings * modelSettings, int dim)
  :Wavelet(modelSettings, dim)
{
  WaveletReadSgri(fileName);
  printf("\nERROR: Only reading of Sgri wavelets implemented. No processing. Program terminates.\n");
}

Wavelet3D::Wavelet3D(Wavelet * wavelet, int dim)
  : Wavelet(wavelet, dim)
{
}

void
Wavelet3D::WaveletReadSgri(char * fileName)
{
  readtype_ = Wavelet::SGRI;
	Sgri sgri(fileName, errText_, errCode_);
  return;
}

void
Wavelet3D::resample(float dz, int nz, float pz, float theta) 
{
  theta_=theta;
  dz_ = dz;
  nz_ = nz;
  nzp_ = (int) pz;
}

bool           
Wavelet3D::consistentSize(int nzp, int nyp, int nxp) const 
{ 
  bool ok = true;
  if (nzp!=nzp_) {
    printf("nzp=%d  nzp_wavelet#D=%d\n",nzp,nzp_);
    ok = false;
  }
  if (nyp != nyp_) {
    printf("nyp=%d  nyp_wavelet3D=%d\n",nyp,nyp_);
    ok = false;
  }
  if (nxp != nxp_) {
    printf("nxp=%d  nxp_wavelet3D=%d\n",nxp,nxp_);
    ok = false;
  }
  return (ok);
}
