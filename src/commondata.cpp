/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/commondata.h"

CommonData::CommonData(void)
{
}

CommonData::~CommonData()
{
}

void
CommonData::createOuterTempSimbox()
{
  //Simbox: Flat topp og bunn, gitt av hhv. topp av øverste og bunn av nederste flate. Rund utover til nærmeste 4ms.
  //Hent segygeometri fra første AVO-vintage, første vinkel. (Få evt. med segyformatet).
  //Les inn ytterflatene for inversjonen.
  //Med dette: finn xy-utstrekninga (se makeTimeSimboxes)

  //Lag en simbox xy-utstrekninga, flat topp og flat bunn. -> esimation_simbox_
  //Lag et volum med samme xy-utstrekning. ->full_inversion_volume_
}

