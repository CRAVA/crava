// moheres/moheres_lib/irapgrid.h, moheres_lib, storm1.1

/*Func: irapgrid 

---------------------------------------------------------------- 
Name:      irapgrid - Header file 
Syntax:     @irapgrid-syntaks
Description:  Structure declarations keeping irap grid and
prototype declarations.

Files:  Source file irapgrid.c
Author:     O.L 
Date:        1993
End:

----------------------------------------------------------------
*/

#ifndef IRAPGRID_H
#define IRAPGRID_H 1

/*<irapgrid-syntaks:*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib/global_def.h"

#define IRAPPRINT 0          /* 1 If printing MESSAGES, 0 - no MESSAGE */
#define IRAPMISSING WELLMISSING /* Missing code defined */
#define ABS_TOLERANCE 0.1    /* If two maps differ by less than ABS_TOLERANCE they are
said to be equal */
#define REL_TOLERANCE 0.01  /* If two maps differ by REL_TOLERANCE*(typical size)
where typical size is the average or maximum difference
between the two surfaces it is said to be equal */

class Point;

struct irapgrid
{int nx,ny;
double xinc,yinc;
double xmin,xmax,ymin,ymax;
double *grid;
double missingcode;
int bin;
char *filename;
double constValue;
};


struct firapgrid
{int nx,ny;
double xinc,yinc;
double xmin,xmax,ymin,ymax;
float *grid;
float missingcode;
int bin;
char *filename;
float constValue;
};


#ifdef __cplusplus
extern "C"
{
#endif

  /* Reading irap grid from file */
  struct irapgrid *irapgridRead(char *, int *); 

  /*Writing irap grid to file */
  int irapgridWrite(char *, struct irapgrid *);

  /* Reading irap grid from file in BINARY format */
  struct irapgrid *irapgridReadBin(char *, int *);

  /*Writing irap grid to file in BINARY format */
  int irapgridWriteBin(char *, struct irapgrid *);

  /*Writing  grid (without header info) to file in float BINARY format */
  int XitegridWriteBin(char *, struct irapgrid *);

  /* Reading irap grid from file specified by pointer*/
  struct irapgrid *irapgridReadpt(FILE *, int *); 

  /*Writing irap grid to file specified by pointer */
  int irapgridWritept(FILE *, struct irapgrid *);

  /* Reading irap grid from file specified by pointer in BINARY format */
  struct irapgrid *irapgridReadBinpt(FILE *, int *);

  /*Writing irap grid to file specified by pointer in BINARY format */
  int irapgridWriteBinpt(FILE *, struct irapgrid *);


  /* Get value in position (Ix,Iy) from grid */
  double irapgridGetValue(double ,double ,
  struct irapgrid *,int *);

  /* Get value in node nearest position (Ix,Iy) from grid */
  double irapgridGetNearestValue(double ,double ,
  struct irapgrid *,int *);

  /* Assign value to nodes around position (x,y) */
  int irapgridSetValue(double ,double ,
  struct irapgrid *,double );

  /* Free space allocated for struct irapgrid */
  void freeIrapgrid(struct irapgrid *);

  /* Add irap grids */
  struct irapgrid *irapgridAddGrid(struct irapgrid *,
  struct irapgrid *,int *);
  /* Subtract irap grids */
  struct irapgrid *irapgridSubGrid(struct irapgrid *,
  struct irapgrid *,int *);


  /* Multiply irap grids */
  struct irapgrid *irapgridMultGrid(struct irapgrid *,
  struct irapgrid *,
    int *);
  /* Divide irap grids */
  struct irapgrid *irapgridDivGrid(struct irapgrid *,
  struct irapgrid *,
    int *);

  /* Copy irap grid */
  struct irapgrid *irapgridCopyGrid(struct irapgrid *,int *);

  /* Erode irap grid */
  struct irapgrid *irapgridErodeGrid(struct irapgrid *,
  struct irapgrid *,int *);

  /* Erode irap grids */
  void irapgridErodeMultiGrid(struct irapgrid *,
  struct irapgrid **, int, int *);

  /* Aritmetic operations between irap grid and a constant */
  struct irapgrid *irapgridArithmeticConstant(struct irapgrid *,
    double,int, int *);

  /* Construct a constant grid */
  struct irapgrid *irapgridConstGrid(double xmin, double ymin, 
    double xinc, double yinc,
    int nx, int ny, 
    double constant,
    int *failure);

  /* Calculate max. of two grids node by node */
  struct irapgrid *irapgridMaxGrid(struct irapgrid *irap1,
  struct irapgrid *irap2,
    int *failure);

  /* Calculate min. of two grids node by node */
  struct irapgrid *irapgridMinGrid(struct irapgrid *irap1,
  struct irapgrid *irap2,
    int *failure);

  /* Calculate minimum and maximum value of a grid */
  void irapgridMinMaxValue(struct irapgrid *irap1,
    double *minvalue, double *maxvalue);

  /* Calculate thickness of formation after erosion */
  struct irapgrid *irapgridZoneThickness(struct irapgrid *irap1,
  struct irapgrid *irap2,
    double dz,
    int *failure);

  /* Calculate top erosion thickness */
  struct irapgrid *irapgridTopErosionThickness(struct irapgrid *iraptop,
  struct irapgrid *irapbot,
    double dz, 
    double erosionfactor,
    int *failure);

  /* Check whether irap1 cross irap2 or not 
  and use average value of these two surfaces cross eachother */
  void irapgridCheckCrossing(struct irapgrid *irap1,
  struct irapgrid *irap2, 
    int *failure);

  int irapgridCheckConsistency(struct irapgrid *top,
  struct irapgrid *bot,
  struct irapgrid *topEro,
  struct irapgrid *botEro,
    double dz);

  /* Check if two irapgrids are identical */
  int irapgridNotIdentical(struct irapgrid *irap1,struct irapgrid *irap2);

  /* check whether a rectangular area is within the grid */
  int rectangleWithinIrapgrid(Point *cornerPoint, struct irapgrid *grid);

  /* Reading irap grid from file specified by pointer. Floating point grid */
  struct firapgrid *firapgridReadpt(FILE *file, int *failure); 

  /*Writing irap grid to file specified by pointer. Floating point grid  */
  int firapgridWritept(FILE *file, struct firapgrid *irap);

  /* Reading irap grid from file specified by pointer in BINARY format.
  Floating point grid  */
  struct firapgrid *firapgridReadBinpt(FILE *, int *);

  /*Writing irap grid to file specified by pointer in BINARY format.
  Floating point grid  */
  int firapgridWriteBinpt(FILE *, struct firapgrid *);


  int get2DGridInfo(
    char *fileName,
    int *ncol,
    int *nrow,
    double * x1,
    double * x2,
    double * y1,
    double * y2,
    double * xinc,
    double * yinc);

  void fget2dgridinfo(
    char *fileName,
    int *ncol,
    int *nrow,
    double * x1,
    double * x2,
    double * y1,
    double * y2,
    double * xinc,
    double * yinc,
    int *ierr,
    int fileNameL);

#ifdef __cplusplus
}
#endif

/*>irapgrid-syntaks:*/
#endif

