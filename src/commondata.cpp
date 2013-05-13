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
  //Simbox: Flat topp og bunn, gitt av hhv. topp av �verste og bunn av nederste flate. Rund utover til n�rmeste 4ms.
  //Hent segygeometri fra f�rste AVO-vintage, f�rste vinkel. (F� evt. med segyformatet).
  //Les inn ytterflatene for inversjonen.
  //Med dette: finn xy-utstrekninga (se makeTimeSimboxes)

  //Lag en simbox xy-utstrekninga, flat topp og flat bunn. -> esimation_simbox_
  //Lag et volum med samme xy-utstrekning. ->full_inversion_volume_
}

