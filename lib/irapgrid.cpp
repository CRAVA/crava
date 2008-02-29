#if (defined SGI || defined SUN)
#define fget2dgridinfo  fget2dgridinfo_
#endif

#include "lib/lib_misc.h"
#include "lib/irapgrid.h"
#include "lib/log.h"

// moheres/moheres_lib/irapgrid.C, moheres_lib, storm1.1

/*
________________________________________________________________

irapgridRead
________________________________________________________________

Description: 	Function for reading grid from file specified 
by a string, returning pointer to the grid struct 
and an error indicator. File format is ascii, but this
function can also be used to read binary format
defined by the internal write binary function
'irapgridWriteBin'. This function will identify
the format as ascii or binary depending on the keyword
written in the beginning of the file.

int *failure has value:\\
0: OK \\
-1: could not open file \\
-2: error when reading file \\
-3: error in allocating space for grid \\

Return value:  Pointer to irap grid. If it is not possible to read the file
the return value is NULL.

________________________________________________________________
*/

struct irapgrid *irapgridRead(char *file, int *failure)
{
  FILE *in;
  char text[MAX_STRING];
  int nread,nreadSum,i,npix;
  struct irapgrid *rIrap;
  float dummy[4];

  *failure = 0;
  in = fopen(file,"r");
  if (in == NULL || strcmp(file,"")==0) 
  {*failure = -1;
  return NULL;
  }
  nread = fscanf(in,"%s",text);

  if(strstr(text,"BINARY")){
    fclose(in);
    rIrap = irapgridReadBin(file,failure);
    return rIrap;
  }
  fclose(in);

  rIrap = (struct irapgrid *) calloc(1, sizeof(struct irapgrid));
  if (rIrap == NULL)
  {
    *failure = -3;
    nread = fclose(in);
    in = NULL;
    return NULL;
  }

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE reading irap file: %s\n",file);

  in = fopen(file,"r");

  nread = fscanf(in,"%d %d %e %e",&(rIrap->nx),&(rIrap->ny),
    &(dummy[0]),&(dummy[1]));
  if(rIrap->nx <= 0 || rIrap->ny <= 0) {
    /*     LogKit::writeLog("ERROR: number of nodes in irapgrid are < = 0 "); */
    *failure = -2;
    return NULL;
  }
  npix = rIrap->nx * rIrap->ny;
  rIrap->xinc = (double) dummy[0];
  rIrap->yinc = (double) dummy[1];
  if (nread != 4)
  {*failure = -2;
  nread = fclose(in);
  in = NULL;
  free (rIrap);
  rIrap = NULL;
  return NULL;
  };
  rIrap->grid = (double *) malloc(((unsigned) npix)*sizeof(double));
  if (rIrap->grid == NULL)
  {
    *failure = -3;
    nread = fclose(in);
    in = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };
  nread = fscanf(in,"%lf %lf %lf %lf",&(rIrap->xmin),&(rIrap->xmax),
    &(rIrap->ymin),&(rIrap->ymax));
  if (nread != 4)
  {
    *failure = -2;
    nread = fclose(in);
    in = NULL;
    free((void *)rIrap->grid);
    rIrap->grid = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };

  nreadSum = 0;
  for (i = 0; i < npix; i++)
  {
    nread = fscanf(in,"%lf ",&(rIrap->grid[i]));
    nreadSum +=nread;
  };

  if (nreadSum != npix)
  {
    *failure = -2;
    nread = fclose(in);
    in = NULL;
    free((void *)rIrap->grid);
    rIrap->grid = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };

  nread = fclose(in);
  in = NULL;
  rIrap->constValue = IRAPMISSING;
  rIrap->missingcode = IRAPMISSING;
  rIrap->filename = (char *) calloc(strlen(file)+1,sizeof(char));
  strcpy(rIrap->filename,file);
  rIrap->bin = 0;
  return rIrap;
}





/*F:irapgridWrite*

________________________________________________________________

irapgridWrite
________________________________________________________________

Name:		irapgridWrite
Syntax:		@irapgridWrite-syntax
Description: 	Function for writing grid to file specified 
by a string, returning an error indicator. 
File format is ascii.

Return value:  0 if writing succesfull \\
-1 if error when opening file \\
-2 if error under writing \\

Author:		Oddvar Lia, NR
Date:   Nov 1993

________________________________________________________________

*/

/*<irapgridWrite-syntax: */
int irapgridWrite(char *file, struct irapgrid *irap)
/*>irapgridWrite-syntax: */
{FILE *pf;
int i,nwritten,n,npix;

pf = fopen(file,"w");
if (pf == NULL || strcmp(file,"")==0) 
return -1;

if(irap == NULL){
  LogKit::writeLog("ERROR (irapgridWrite): Pointer irap = NULL\n"); 
  return -2;
}

if(irap->grid == NULL){
  LogKit::writeLog("ERROR (irapgridWrite): Pointer irap->grid = NULL\n"); 
  return -2;
}


if(IRAPPRINT)
LogKit::writeLog("MESSAGE writing irap file: %s\n",file);

npix = irap->nx*irap->ny;
nwritten = fprintf(pf,"%6d%6d%14.6f%14.6f\n",irap->nx,irap->ny,
                   irap->xinc,irap->yinc);
if (nwritten < 0)
return -2; 

nwritten = fprintf(pf,"%14.4f%14.4f%14.4f%14.4f\n",irap->xmin,irap->xmax,
                   irap->ymin,irap->ymax);
if (nwritten < 0) return -2; 

n = 0;
for ( i = 0; i < npix; i++) {
  n++;
  nwritten = fprintf(pf," %12.4f",irap->grid[i]);
  if(nwritten < 0) return -2;
  if(n == 6){
    n = 0;
    fprintf(pf,"\n");
  }
}
fprintf(pf,"\n");

nwritten = fclose(pf);
pf = NULL;

return 0; 
}








/*F:irapgridReadBin*

________________________________________________________________

irapgridReadBin
________________________________________________________________

Name:		irapgridReadBin
Syntax:		@irapgridReadBin-syntax
Description: 	Function for reading grid from file specified by a string, 
returning pointer to the grid struct 
and an error indicator. The file format is
binary.
The format is defined by the function 
'int irapgridWriteBin'.


int *failure  has the values:\\
0: OK \\
-1: could not open file\\
-2: error when reading file\\
-3: error in allocating space for grid\\

Return value: Pointer to an array with the irap structure. 
If read error occurs, the pointer is NULL.
Author:		Oddvar Lia, NR
Date:          Nov 1993
________________________________________________________________

*/

/*<irapgridReadBin-syntax: */
struct irapgrid *irapgridReadBin(char *file, int *failure) 
  /*>irapgridReadBin-syntax: */
{
  FILE *in;
  int nread,npix;
  struct irapgrid *rIrap;
  float dummy[4];
  char c;
  char text[MAX_STRING];

  *failure = 0;
  in = fopen(file,"rb");
  if (in == NULL || strcmp(file,"")==0) 
  {*failure = -1;
  return NULL;
  };
  rIrap = (struct irapgrid *) calloc(1, sizeof(struct irapgrid));
  if (rIrap == NULL)
  {*failure = -3;
  nread = fclose(in);
  in = NULL;
  return NULL;
  };

  nread = fscanf(in,"%s",text);

  if(!strstr(text,"BINARY")){
    fclose(in);
    *failure = -2;
    rIrap = NULL;
    return rIrap;
  }

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE reading irap file: %s\n",file);

  nread = fscanf(in,"%d %d %e %e",
    &(rIrap->nx),&(rIrap->ny),&(dummy[0]),&(dummy[1]));
  if(rIrap->nx <= 0 || rIrap->ny <= 0) {
    /*     LogKit::writeLog("ERROR: number of nodes in irapgrid are < = 0 "); */
    *failure = -2;
    return NULL;
  }

  rIrap->xinc = (double) dummy[0];
  rIrap->yinc = (double) dummy[1];
  if (nread != 4)
  {*failure = -2;
  nread = fclose(in);
  in = NULL;
  free (rIrap);
  rIrap = NULL;
  return NULL;
  };


  rIrap->grid = (double *) malloc(((unsigned)(rIrap->nx * rIrap->ny))* 
    sizeof(double));
  if (rIrap->grid == NULL)
  {*failure = -3;
  nread = fclose(in);
  in = NULL;
  free(rIrap);
  rIrap = NULL;
  return NULL;
  };


  nread = fscanf(in,"%lf %lf %lf %lf",&(rIrap->xmin),&(rIrap->xmax),
    &(rIrap->ymin),&(rIrap->ymax));
  if (nread != 4)
  {*failure = -2;
  nread = fclose(in);
  in = NULL;
  free((void *)rIrap->grid);
  rIrap->grid = NULL;
  free(rIrap);
  rIrap = NULL;
  return NULL;
  };

  /* Get new line */
  for (c = char(fgetc(in)); !((c == '\n') || (c == EOF) || (c == '\377'));
    c = char(fgetc(in)));

    npix = rIrap->nx*rIrap->ny;

  nread = fread(&(rIrap->grid[0]),sizeof(double),npix,in);
  if(nread != npix){
    *failure = -2;
    free(&(rIrap->grid[0]));
    free(rIrap);
    return NULL;
  }
#ifndef BIGENDIAN
  swapDoubles(&(rIrap->grid[0]),npix);
#endif

  fclose(in);
  in = NULL;
  rIrap->constValue = IRAPMISSING;
  rIrap->missingcode = IRAPMISSING;
  rIrap->bin = 1;
  rIrap->filename = (char *) calloc(strlen(file)+1,sizeof(char));
  strcpy(rIrap->filename,file);
  return rIrap;
}





/*F:irapgridWriteBin*

________________________________________________________________

irapgridWriteBin
________________________________________________________________

Name:		irapgridWriteBin
Syntax:		@irapgridWriteBin-syntax
Description: 	Function for writing grid to file specified by a string, 
returning an error indicator. The file format is binary.

Return value:  0 if writing succesfull\\
-1 if error when opening file\\
-2 if error under writing\\

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridWriteBin-syntax: */
int irapgridWriteBin(char *file, struct irapgrid *irap)
/*>irapgridWriteBin-syntax: */
{FILE *pf;
int npix,nwritten;

pf = fopen(file,"wb");
if (pf == NULL || strcmp(file,"")==0) 
return -1;

if(irap == NULL){
  LogKit::writeLog("ERROR (irapgridWriteBin): Pointer irap = NULL\n"); 
  return -2;
}

if(irap->grid == NULL){
  LogKit::writeLog("ERROR (irapgridWriteBin): Pointer irap->grid = NULL\n"); 
  return -2;
}


if(IRAPPRINT)
LogKit::writeLog("MESSAGE writing irap file: %s\n",file);

fprintf(pf,"STORMGRID_BINARY\n\n");
nwritten = fprintf(pf,"%d %d %f %f \n",irap->nx,irap->ny,
                   irap->xinc,irap->yinc);
if (nwritten < 0)
return -2; 

nwritten = fprintf(pf,"%f %f %f %f \n",irap->xmin,irap->xmax,
                   irap->ymin,irap->ymax);
if (nwritten < 0)
return -2; 

npix = irap->nx*irap->ny;
#ifndef BIGENDIAN
swapDoubles(&(irap->grid[0]),npix);
#endif
nwritten = fwrite(&(irap->grid[0]),sizeof(double),npix,pf);
#ifndef BIGENDIAN
swapDoubles(&(irap->grid[0]),npix);
#endif
if(nwritten != npix){
  return -2;
}

fclose(pf);
pf = NULL;

return 0; 
}






/*F:XitegridWriteBin*

________________________________________________________________

XitegridWriteBin
________________________________________________________________

Name:		XitegridWriteBin
Syntax:		@XitegridWriteBin-syntax
Description: Writes a irap grid without the header information in float
BINARY format. The output file is a rawdata file to be input
to the raw2biff program of the Xite system.
Note that the grid cells are rearranged in y direction.
(x-coord. runs fastest).

Side effects:
Return value:  0 if writing succesfull\\
-1 if error when opening file\\
-2 if error under writing\\

Author:		Rolf Clemetsen, NR
Date:         Nov 1993
________________________________________________________________

*/

/*<XitegridWriteBin-syntax: */
int XitegridWriteBin(char *file, struct irapgrid *irap)
/*>XitegridWriteBin-syntax: */
{FILE *pf;
int npix,nwritten,i,j,ind1,ind2;
float *floatgrid;

pf = fopen(file,"w");
if (pf == NULL || strcmp(file,"")==0) 
return -1;

if(IRAPPRINT)
LogKit::writeLog("MESSAGE writing Xite file: %s\n",file);

npix = irap->nx*irap->ny;

floatgrid = (float *) malloc(((unsigned)(npix))* sizeof(float));
if ( floatgrid == NULL ) return -3;

for ( j=0; j<irap->ny; j++ ) {
  ind1 = j * irap->nx;
  ind2 = (irap->ny -1 -j) * irap->nx; 
  for ( i=0; i<irap->nx; i++ ) 
    floatgrid[ind1 + i] = float(irap->grid[ind2 + i]);
}
nwritten = fwrite(&(floatgrid[0]),sizeof(float),npix,pf);
free ( floatgrid );
if(nwritten != npix)  return -2;

fclose(pf);
pf = NULL;

return 0; 
}






/*F:irapgridReadpt*

________________________________________________________________

irapgridReadpt
________________________________________________________________

Name:		irapgridReadpt
Syntax:		@irapgridReadpt-syntax
Description: 	Function for reading grid from file specified 
by a pointer, returning pointer to the grid 
struct and an error indicator. File format is ascii.

int *failure  has the value:\\
0: OK\\
-1: error file is not open\\
-2: error when reading file\\
-3: error in allocating space for grid\\

Side effects:
Return value: Pointer to irapgrid, NULL if the file is cannot be read.
Author:		Oddvar Lia, NR
Date:  Nov 1993
________________________________________________________________

*/

/*<irapgridReadpt-syntax: */
struct irapgrid *irapgridReadpt(FILE *file, int *failure) 
  /*>irapgridReadpt-syntax: */
{
  int nread,nreadSum,i,npix;
  struct irapgrid *rIrap;
  float dummy[4];

  *failure = 0;
  if(file == NULL){
    *failure = -1;
    return NULL;
  }

  rIrap = (struct irapgrid *) malloc(sizeof(struct irapgrid));
  if (rIrap == NULL)
  {*failure = -3;
  return NULL;
  };

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE reading irap grid\n");

  nread = fscanf(file,"%d %d %e %e",&(rIrap->nx),&(rIrap->ny),
    &(dummy[0]),&(dummy[1]));
  if(rIrap->nx <= 0 || rIrap->ny <= 0) {
    /*     LogKit::writeLog("ERROR: Number of nodes in irapgrid are < = 0 "); */
    *failure = -2;
    return NULL;
  }
  npix = rIrap->nx * rIrap->ny;
  rIrap->xinc = (double) dummy[0];
  rIrap->yinc = (double) dummy[1];
  if (nread != 4)
  {*failure = -2;
  free (rIrap);
  rIrap = NULL;
  return NULL;
  };
  rIrap->grid = (double *) malloc(((unsigned) npix)*sizeof(double));
  if (rIrap->grid == NULL)
  {
    *failure = -3;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };
  nread = fscanf(file,"%lf %lf %lf %lf",&(rIrap->xmin),&(rIrap->xmax),
    &(rIrap->ymin),&(rIrap->ymax));
  if (nread != 4)
  {
    *failure = -2;
    free((void *)rIrap->grid);
    rIrap->grid = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };

  nreadSum = 0;
  for (i = 0; i < npix; i++)
  {
    nread = fscanf(file,"%lf ",&(rIrap->grid[i]));
    nreadSum +=nread;
  };

  if (nreadSum != npix)
  {
    *failure = -2;
    free((void *)rIrap->grid);
    rIrap->grid = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };
  rIrap->constValue = IRAPMISSING;
  rIrap->missingcode = IRAPMISSING;
  rIrap->bin = 0;
  return rIrap;
}




/*F:irapgridWritept*

________________________________________________________________

irapgridWritept
________________________________________________________________

Name:		irapgridWritept
Syntax:		@irapgridWritept-syntax
Description: 	Function for writing grid to file specified by a 
pointer, returning an error indicator. 
The file format is ascii.


Return value:  0 if writing succesfull\\
-1 if file is not opened\\
-2 if error under writing\\

Author:		Oddvar Lia, NR
Date:   Nov 1993
________________________________________________________________

*/

/*<irapgridWritept-syntax: */
int irapgridWritept(FILE *file, struct irapgrid *irap)
/*>irapgridWritept-syntax: */
{
  int i,nwritten,n,npix;

  if (file == NULL) return -1;


  if(irap == NULL){
    LogKit::writeLog("ERROR (irapgridWritept): Pointer irap = NULL\n"); 
    return -2;
  }

  if(irap->grid == NULL){
    LogKit::writeLog("ERROR (irapgridWritept): Pointer irap->grid = NULL\n"); 
    return -2;
  }

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE writing irap grid\n");


  npix = irap->nx*irap->ny;
  nwritten = fprintf(file,"%d %d %f %f \n",irap->nx,irap->ny,
    irap->xinc,irap->yinc);
  if (nwritten < 0)
    return -2; 

  nwritten = fprintf(file,"%f %f %f %f \n",irap->xmin,irap->xmax,
    irap->ymin,irap->ymax);
  if (nwritten < 0) return -2; 

  n = 0;
  for ( i = 0; i < npix; i++) {
    n++;
    nwritten = fprintf(file,"%f ",irap->grid[i]);
    if(nwritten < 0) return -2;
    if(n == 6){
      n = 0;
      fprintf(file,"\n");
    }
  }
  fprintf(file,"\n");

  return 0; 
}





/*F:irapgridReadBinpt*

________________________________________________________________

irapgridReadBinpt
________________________________________________________________

Name:		irapgridReadBinpt
Syntax:		@irapgridReadBinpt-syntax
Description: 	Function for reading grid from file specified 
by a pointer, returning 
pointer to the grid struct and an error indicator. 
File format is binary.
The binary format is the format 'irapgridWriteBinpt' uses.

int *failure has the value:\\
0: OK\\
-1: file is not open\\
-2: error when reading file\\
-3: error in allocating space for grid\\

Return value: Pointer to an array with the irap structure. If read error, the
pointer is NULL.

Author:		Oddvar Lia, NR
Date:  Nov 1993
________________________________________________________________

*/

/*<irapgridReadBinpt-syntax: */
struct irapgrid *irapgridReadBinpt(FILE *file, int *failure) 
  /*>irapgridReadBinpt-syntax: */
{
  int nread,npix;
  struct irapgrid *rIrap;
  float dummy[4];
  char c;

  *failure = 0;
  if (file == NULL) 
  {*failure = -1;
  return NULL;
  };
  rIrap = (struct irapgrid *) malloc(sizeof(struct irapgrid));
  if (rIrap == NULL)
  {*failure = -3;
  return NULL;
  };

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE reading irap grid\n");

  nread = fscanf(file,"%d %d %e %e",
    &(rIrap->nx),&(rIrap->ny),&(dummy[0]),&(dummy[1]));
  if(rIrap->nx <= 0 || rIrap->ny <= 0) {
    LogKit::writeLog("%s \n","WARNING: number of nodes in irapgrid are < = 0 ");
    *failure = -2;
    return NULL;
  }

  rIrap->xinc = (double) dummy[0];
  rIrap->yinc = (double) dummy[1];
  if (nread != 4)
  {*failure = -2;
  free (rIrap);
  rIrap = NULL;
  return NULL;
  };


  rIrap->grid = (double *) malloc(((unsigned)(rIrap->nx * rIrap->ny))* 
    sizeof(double));
  if (rIrap->grid == NULL)
  {*failure = -3;
  free(rIrap);
  rIrap = NULL;
  return NULL;
  };


  nread = fscanf(file,"%lf %lf %lf %lf",&(rIrap->xmin),&(rIrap->xmax),
    &(rIrap->ymin),&(rIrap->ymax));
  if (nread != 4)
  {*failure = -2;
  free((void *)rIrap->grid);
  rIrap->grid = NULL;
  free(rIrap);
  rIrap = NULL;
  return NULL;
  };

  /* Get new line */
  for (c = char(fgetc(file)); !((c == '\n') || (c == EOF) || (c == '\377'));
    c = char(fgetc(file)));

    npix = rIrap->nx*rIrap->ny;

  nread = fread(&(rIrap->grid[0]),sizeof(double),npix,file);
  if(nread != npix){
    *failure = -2;
    free(&(rIrap->grid[0]));
    free(rIrap);
    return NULL;
  }

  rIrap->constValue = IRAPMISSING;
  rIrap->missingcode = IRAPMISSING;
  rIrap->bin = 1;
  return rIrap;
}




/*F:irapgridWriteBinpt*

________________________________________________________________

irapgridWriteBinpt
________________________________________________________________

Name:		irapgridWriteBinpt
Syntax:		@irapgridWriteBinpt-syntax
Description: 	Function for writing grid to file specified by a pointer, 
returning  an error indicator. The file format is binary.


Return value:  0 if writing succesfull\\
-1 if file is not open\\
-2 if error under writing\\

Author:		Oddvar Lia, NR
Date: Nov 1993
________________________________________________________________

*/

/*<irapgridWriteBinpt-syntax: */
int irapgridWriteBinpt(FILE *file, struct irapgrid *irap)
/*>irapgridWriteBinpt-syntax: */
{
  int npix,nwritten;

  if (file == NULL) 
    return -1;



  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE writing irap grid\n");

  if(irap == NULL){
    LogKit::writeLog("ERROR (irapgridWriteBinpt): Pointer irap = NULL\n"); 
    return -2;
  }

  if(irap->grid == NULL){
    LogKit::writeLog("ERROR (irapgridWriteBinpt): Pointer irap->grid = NULL\n"); 
    return -2;
  }


  nwritten = fprintf(file,"%d %d %f %f \n",irap->nx,irap->ny,
    irap->xinc,irap->yinc);
  if (nwritten < 0)
    return -2; 

  nwritten = fprintf(file,"%f %f %f %f \n",irap->xmin,irap->xmax,
    irap->ymin,irap->ymax);
  if (nwritten < 0)
    return -2; 

  npix = irap->nx*irap->ny;
  nwritten = fwrite(&(irap->grid[0]),sizeof(double),npix,file);
  if(nwritten != npix){
    return -2;
  }

  return 0; 
}





/*
________________________________________________________________

irapgridGetValue
________________________________________________________________

Description: 
Function that returns a value from irap grid at specified coordinate.
The value is calculated by bilinear interpolation of the 4 surrounding
node values. If any of the surrounding nodes have missing code a
pooled value of the existing non-missing values is returned.
If the coordinate is outside the grid area or if all 4 surrounding
nodes have missing value, the missingvalue is returned.
An error indicator returns 0 or 1 depending on whether the point is
within or outside the defined part of the grid.

Return value: The interpolated z-value.
________________________________________________________________
*/

double irapgridGetValue(double x, double y,
struct irapgrid *irap, int *outside)
{
  double	rz, sumfacs, dx, dy, dist;
  double	x0, y0, x1, y1, a, b, c, d;
  int	pix, i, j, ii, jj, nMissing;
  int	indx[5];
  double	z[5], fac[5];

  *outside = 0;

  // Return missing if coordinate is outside grid area
  if (x < irap->xmin || x > irap->xmax ||
    y < irap->ymin || y > irap->ymax)
  {
    *outside = 1;
    return irap->missingcode;
  }

  // Find the lower left surrounding grid node 
  i = (int) ((x - irap->xmin) / irap->xinc);
  j = (int) ((y - irap->ymin) / irap->yinc);

  // Find array indexes of four surrounding points 
  indx[1] = i + j*irap->nx;
  indx[2] = indx[1] + 1;
  indx[3] = indx[1] + irap->nx;
  indx[4] = indx[1] + 1 + irap->nx;

  // Check if point is exactly at the eastern edge 
  if (i == irap->nx - 1)
  {
    indx[2] = IMISSING;	
    indx[4] = IMISSING;
  }	

  // Check if point is exactly at the northern edge 
  if (j == irap->ny - 1)
  {
    indx[3] = IMISSING;	
    indx[4] = IMISSING;
  }	

  // Find grid values of the four points
  nMissing = 0;
  for (pix = 1; pix <= 4; pix++)
  {
    if (indx[pix] != IMISSING)
    {
      z[pix] = irap->grid[indx[pix]];
    }
    else
    {
      z[pix] = irap->missingcode;
    }

    if (z[pix] == irap->missingcode)
    {
      nMissing++;
    }             
  }

  if (nMissing == 4) // All surrounding values are missing
  {
    return irap->missingcode;
  }

  if (nMissing == 0) // All values present => Bilinear interpolation
  {
    x0 = irap->xmin + i * irap->xinc;
    y0 = irap->ymin + j * irap->yinc;

    x1 = (x - x0) / irap->xinc;
    y1 = (y - y0) / irap->yinc;

    a = z[1];
    b = z[2] - a;
    c = z[3] - a;
    d = z[4] - a - b - c;

    return a + b*x1 + c*y1 + d*x1*y1;
  }

  // Calculate pooled value of the 1, 2 or 3 non-missing nodes 

  sumfacs = 0.0;
  for (jj = 0; jj<= 1; jj++)
  {
    for (ii = 0; ii<= 1; ii++)
    {                         
      pix = 1 + ii + 2 * jj;
      dx = irap->xmin + (i+ii)*irap->xinc - x;
      dy = irap->ymin + (j+jj)*irap->yinc - y;
      dist = sqrt(dx*dx + dy*dy);

      if(dist < 0.000001)
      {        
        return z[pix];
      }

      if(z[pix] != irap->missingcode)
      {
        fac[pix] = 1.0/dist;
        sumfacs += fac[pix];
      }
      else
      {
        fac[pix] = 0.0;
        z[pix] = 0.0;
      }
    }
  }

  rz = 0.0;
  for (pix = 1; pix <= 4; pix++) 
  {
    rz += z[pix] * fac[pix];
  }

  return rz/sumfacs;
}


/*F:irapgridGetNearestValue*

________________________________________________________________

irapgridGetNearestValue
________________________________________________________________

Name:		irapgridGetNearestValue
Syntax:		@irapgridGetNearestValue-syntax
Description: 
Function that returns the value in the nearest node in the 
irap grid at specified coordinate. If the nearest node have 
missing code or if the point is outside the grid limits 
then the missing code is returned.

An error indicator returns 0 or 1 depending on whether the point is
within or outside the defined part of the grid.

Return value: The z-value in nearest node

Author:  Anne-Lise Hektoen, NR

Date:    Jan. 1994

________________________________________________________________

*/

/*<irapgridGetNearestValue-syntax: */
double irapgridGetNearestValue(double x,double y,
struct irapgrid *irap,int *outside)
  /*>irapgridGetNearestValue-syntax: */
{double rz;
int indx;
int i,j;

*outside = 0;

i = (int) ((x - irap->xmin) / irap->xinc + 0.5);
j = (int) ((y - irap->ymin) / irap->yinc + 0.5);

if (irap->nx == 1 && irap->ny == 1)
return irap->grid[0];

else if (i < 0 || j < 0 || i > irap->nx || j > irap->ny)
{*outside = 1;
return irap->missingcode;
}

indx = i + j*irap->nx ;

rz = irap->grid[indx];
if (rz == irap->missingcode) {
  *outside = 1;
  return irap->missingcode;
}

return rz;
}




/*F:irapgridSetValue*

________________________________________________________________

irapgridSetValue
________________________________________________________________

Name:		irapgridSetValue
Syntax:		@irapgridSetValue-syntax
Description: 
Function that locally modify the irap grid at a
specified coordinate. 4 nearest nodes (used in bilinear 
interpolation in irapgridGetValue()) is modified.

Return value:  0 if OK, 
-1 if position outside grid.

Author:	        Anne-Lise Hektoen, NR
Date:  May 1994

________________________________________________________________

*/

/*<irapgridSetValue-syntax: */
int irapgridSetValue(double x,double y,
struct irapgrid *irap,double value)
  /*>irapgridSetValue-syntax: */
{int indx[5],pix;
int i,j;

i = (int) ((x - irap->xmin) / irap->xinc) ;
j = (int) ((y - irap->ymin) / irap->yinc) ;

if (irap->nx == 1 && irap->ny == 1)
irap->grid[0] = value;

else if (i < 0 || j < 0 || i > irap->nx || j > irap->ny)
return -1;

indx[1] = i + j*irap->nx ;
indx[2] = indx[1] + 1;
indx[3] = indx[1] + irap->nx;
indx[4] = indx[1] + 1 + irap->nx;

for (pix = 1; pix <= 4; pix++) 
irap->grid[indx[pix]] = value;

return 0;
}




/*F:freeIrapgrid*

________________________________________________________________

freeIrapgrid
________________________________________________________________

Name:		freeIrapgrid
Syntax:		@freeIrapgrid-syntax
Description:
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Anne-Lise Hektoen, NR
Date:
Source file: 
________________________________________________________________

*/

/*<freeIrapgrid-syntax: */
void freeIrapgrid(struct irapgrid *map)
/*>freeIrapgrid-syntax: */
{
  if (map->grid) 
  {
    free(map->grid);
    map->grid = NULL;
  }

  if (map->filename) 
  {
    free(map->filename);
    map->filename = NULL;
  }

  return;
}	/* end of freeIrapgrid */





/*F:irapgridAddGrid*

________________________________________________________________

irapgridAddGrid
________________________________________________________________

Name:		irapgridAddGrid
Syntax:		@irapgridAddGrid-syntax
Description: Add two irap grids that are compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to 
be the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids can not add or undefined grids \\
-2 Allocation error\\

Return value: Pointer to a new irap grid that is the sum of the input
grids.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridAddGrid-syntax: */
struct irapgrid *irapgridAddGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridAddGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapAddGrid): Pointer irap = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapAddGrid): Pointer irap->grid = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapAddGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapAddGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapAddGrid): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
      rirap->constValue = irap1->constValue + irap2->constValue;
    else
      rirap->constValue = rirap->missingcode;


    /*  Add irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          rirap->grid[pix] = irap1->grid[pix] + irap2->grid[pix];
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    rirap->bin = irap1->bin;
    return rirap;
}





/*F:irapgridSubGrid*

________________________________________________________________

irapgridSubGrid
________________________________________________________________

Name:		irapgridSubGrid
Syntax:		@irapgridSubGrid-syntax
Description: 
Subtract two irap grids that are compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to be 
the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

The subtraction is:  irap1 - irap2

Error message through *failure:\\
0 Ok \\
-1 Incompatible grids can not add or undefined grids \\
-2 Allocation error \\

Return value: Pointer to a new irap grid that is the sum of the input
grids.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridSubGrid-syntax: */
struct irapgrid *irapgridSubGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridSubGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapSubGrid): Pointer irap = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapSubGrid): Pointer irap->grid = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapSubGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapSubGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapSubGrid): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
      rirap->constValue = irap1->constValue - irap2->constValue;
    else
      rirap->constValue = rirap->missingcode;


    /*  Subtract irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          rirap->grid[pix] = irap1->grid[pix] - irap2->grid[pix];
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    rirap->bin = irap1->bin;
    return rirap;
}


/*F:irapgridMultGrid*

________________________________________________________________

irapgridMultGrid
________________________________________________________________

Name:		irapgridMultGrid
Syntax:		@irapgridMultGrid-syntax
Description: Multiply two irap grids that are compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined 
to be the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids can not multiply or undefined grids \\
-2 Allocation error\\

Return value: Pointer to a new irap grid that is the product node by node
of the input grids.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridMultGrid-syntax: */
struct irapgrid *irapgridMultGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridMultGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapMultGrid): Pointer irap = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapMultGrid): Pointer irap->grid = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapMultGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapMultGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapMultGrid): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
      rirap->constValue = irap1->constValue * irap2->constValue;
    else
      rirap->constValue = rirap->missingcode;

    /*  Multiply irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          rirap->grid[pix] = irap1->grid[pix] * irap2->grid[pix];
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    rirap->bin = irap1->bin;
    return rirap;
}




/*F:irapgridDivGrid*

________________________________________________________________

irapgridDivGrid
________________________________________________________________

Name:		irapgridDivGrid
Syntax:		@irapgridDivGrid-syntax
Description: Divide two irap grids that are compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to 
be the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids can not divide node by node or
undefined grids \\
-2 Allocation error\\

Return value: Pointer to a new irap grid that is grid1 divided by grid2.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridDivGrid-syntax: */
struct irapgrid *irapgridDivGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridDivGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapDivGrid): Pointer irap = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapDivGrid): Pointer irap->grid = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapDivGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapDivGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapDivGrid): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
      rirap->constValue = irap1->constValue/irap2->constValue;
    else
      rirap->constValue = rirap->missingcode;

    /*  Divide irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode &&
        irap2->grid[pix] != 0.0){
          rirap->grid[pix] = (irap1->grid[pix])/(irap2->grid[pix]);
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    rirap->bin = irap1->bin;
    return rirap;
}




/*F:irapgridArithmeticConstant*

________________________________________________________________

irapgridArithmeticConstant
________________________________________________________________

Name:		irapgridArithmeticConstant
Syntax:		@irapgridArithmeticConstant-syntax
Description: Arithmetic operation between one irap grid and a constant.
The operation is specified by a number for the operation.
The missing code of the  grid is defined to be the missingcode
for the new grid. A missing code is returned  in all nodes
where the input grid has missing code.

Defined values for "int operation" is:
Addition : 1 \\
Subtraction : 2 \\
Multiplication : 3 \\
Division : 4 \\

Error message through *failure:\\
0 Ok\\
-1 Not defined operation \\
-2 Allocation error\\
-3 Internal error \\

Return value: Pointer to a new irap grid that is the arithmetic 
operation node by node between the grid and the 
specified input constant.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridArithmeticConstant-syntax: */
struct irapgrid *irapgridArithmeticConstant(struct irapgrid *irap1,
  double constant,int operation,
  int *failure)
  /*>irapgridArithmeticConstant-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL){
    LogKit::writeLog("ERROR (irapArithmeticConstant): Pointer irap = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->grid == NULL){
    LogKit::writeLog("ERROR (irapArithmeticConstant): Pointer irap->grid = NULL\n"); 
    *failure = -3;
    return rirap;
  }


  npix = irap1->nx*irap1->ny;

  rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
  if(rirap == NULL){
    /*    LogKit::writeLog("ERROR (irapArithmeticConstant): Allocation error\n"); */
    *failure = -2;
    return rirap;
  }


  rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
  if(rirap->grid == NULL){
    /*    LogKit::writeLog("ERROR (irapArithmeticConstant): Allocation error\n"); */
    *failure = -2;
    free(rirap);
    rirap = NULL;
    return rirap;
  }

  rirap->nx = irap1->nx;
  rirap->ny = irap1->ny;
  rirap->xmin = irap1->xmin;
  rirap->ymin = irap1->ymin;
  rirap->xinc = irap1->xinc;
  rirap->yinc = irap1->yinc;
  rirap->xmax = irap1->xmax;
  rirap->ymax = irap1->ymax;
  rirap->missingcode = irap1->missingcode;

  if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5))
    rirap->constValue = irap1->constValue + constant;
  else
    rirap->constValue = rirap->missingcode;

  /*  Arithmetic operation node by node between grid and constant */
  if(operation == 1){
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode){
        rirap->grid[pix] = irap1->grid[pix] + constant;
      }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }

  }
  else if(operation == 2){
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode){
        rirap->grid[pix] = irap1->grid[pix] - constant;
      }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }

  }
  else if(operation == 3){
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode){
        rirap->grid[pix] = irap1->grid[pix] * constant;
      }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }

  }
  else if(operation ==4){
    if(constant != 0.0){
      for(pix=0;pix < npix; pix++){
        if(irap1->grid[pix] != irap1->missingcode){
          rirap->grid[pix] = irap1->grid[pix] / constant;
        }
        else{
          rirap->grid[pix] = rirap->missingcode;
        }
      }
    }
    else {
      /*     LogKit::writeLog("ERROR (irapArithmeticConstant): Constant = 0 in division\n"); */
      *failure = -1;
      free(rirap);
      rirap = NULL;
      return rirap;
    }
  }
  else {
    /*  LogKit::writeLog("ERROR (irapArithmeticConstant): Undefined operation for irapgrid\n"); */
    *failure = -3;
    free(rirap);
    rirap = NULL;
    return rirap;
  }
  rirap->bin = irap1->bin;  
  return rirap;
}








/*F:irapgridCopyGrid*

________________________________________________________________

irapgridCopyGrid
________________________________________________________________

Name:		irapgridCopyGrid
Syntax:		@irapgridCopyGrid-syntax
Description:
Copy one irapgrid to another

Error message through *failure:\\
0 Ok \\
-1 Undefined grid \\
-2 Allocation error \\
Return value: Pointer to a new irap grid that is the copy of the input grids.

Author:		 NR
Date: 
________________________________________________________________

*/

/*<irapgridCopyGrid-syntax: */
struct irapgrid *irapgridCopyGrid(struct irapgrid *irap,
  int *failure)
  /*>irapgridCopyGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap == NULL){
    LogKit::writeLog("ERROR (irapCopyGrid): Pointer irap = NULL\n"); 
    *failure = -1;
    return rirap;
  }

  if(irap->grid == NULL){
    LogKit::writeLog("ERROR (irapCopyGrid): Pointer irap->grid = NULL\n"); 
    *failure = -1;
    return rirap;
  }

  npix = irap->nx*irap->ny;

  rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
  if(rirap == NULL){
    /*    LogKit::writeLog("ERROR (irapCopyGrid): Allocation error\n"); */
    *failure = -2;
    return rirap;
  }


  rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
  if(rirap->grid == NULL){
    /*    LogKit::writeLog("ERROR (irapCopyGrid): Allocation error\n"); */
    *failure = -2;
    free(rirap);
    rirap = NULL;
    return rirap;
  }

  rirap->nx = irap->nx;
  rirap->ny = irap->ny;
  rirap->xmin = irap->xmin;
  rirap->ymin = irap->ymin;
  rirap->xinc = irap->xinc;
  rirap->yinc = irap->yinc;
  rirap->xmax = irap->xmax;
  rirap->ymax = irap->ymax;
  rirap->missingcode = irap->missingcode;
  rirap->constValue = irap->constValue;

  /*  Copy irap grid */
  for(pix=0;pix < npix; pix++){
    rirap->grid[pix] = irap->grid[pix];
  }
  rirap->bin = irap->bin;  
  return rirap;
}





/*F:irapgridErodeGrid*

________________________________________________________________

irapgridErodeGrid
________________________________________________________________

Name:		irapgridErodeGrid
Syntax:		@irapgridErodeGrid-syntax
Description:	Irap grid, grid1, erode down in grid2 if the 
value of grid1 is larger 
than the value of grid2. Return a new pointer to a new grid that is 
the eroded surface.
If one or both grids have missing code the resulting value
is defined to be missing code. A grid can only erode another
if they are of same size and same grid increments and are located
at the same coordinates.


Error message through *failure:\\
0 Ok \\
-1 Incompatible grids can not erode or undefined grids \\
-2 Allocation error \\

Return value: Pointer to a new irap grid that is the resulting surface 
after erosion.

Author:		Anne-Lise Hektoen, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridErodeGrid-syntax: */
struct irapgrid *irapgridErodeGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridErodeGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapErodeGrid): Pointer irap = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapErodeGrid): Pointer irap->grid = NULL\n"); 
    *failure = -2;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapErodeGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapErodeGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapErodeGrid): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;
    rirap->constValue = rirap->missingcode;

    /*  Erode irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          if (irap1->grid[pix] >irap2->grid[pix])
            rirap->grid[pix] = irap1->grid[pix];
          else
            rirap->grid[pix] = irap2->grid[pix];
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    rirap->bin = irap1->bin;  
    return rirap;
}





/*F:irapgridErodeMultiGrid*

________________________________________________________________

irapgridErodeMultiGrid
________________________________________________________________

Name:		irapgridErodeMultiGrid
Syntax:		@irapgridErodeMultiGrid-syntax
Description: 	Irap grid, grid0, erode down in grid1...gridN 
if the value of grid0 is larger than the 
value of grid1...gridN. 
If one or both grids (grid0 and grid1...gridN) 
have missing code the resulting value
is defined to be missing code. 
A grid can only erode another
if they are of same size and same grid increments 
and are located at the same coordinates. 
The resulting grids are returned in grid1...gridN.

Error message through *failure:\\
0 Ok \\
-1 Incompatible grids can not erode or undefined grids \\
-2 Allocation error \\

Author:		Anne-Lise Hektoen, NR

Date: Nov 1993

________________________________________________________________

*/

/*<irapgridErodeMultiGrid-syntax: */
void irapgridErodeMultiGrid(struct irapgrid *irapErode,
struct irapgrid **irap,
  int nGrids,int *failure)
  /*>irapgridErodeMultiGrid-syntax: */
{
  int pix,npix;
  int n;
  *failure = 0;


  if(irapErode == NULL){
    LogKit::writeLog("ERROR (irapErodeMultiGrid): Pointer irap = NULL\n"); 
    *failure = -2;
    return;
  }

  if(irapErode->grid == NULL) {
    LogKit::writeLog("ERROR (irapErodeMultiGrid): Pointer irap->grid NULL\n"); 
    *failure = -2;
    return;
  }



  for (n=0; n<nGrids; n++) {
    if(irap[n] != NULL){
      if(irap[n]->grid == NULL) {
        printf("ERROR (irapErodeMultiGrid): Undefined irap grid \n"); 
        *failure = -2;
        return;
      }

      if(irapErode->nx != irap[n]->nx || irapErode->ny != irap[n]->ny ||
        irapErode->xmin != irap[n]->xmin || 
        irapErode->ymin != irap[n]->ymin ||
        irapErode->xinc != irap[n]->xinc || irapErode->yinc != irap[n]->yinc){
          printf("ERROR (irapErodeMultiGrid): Incompatible irap grids\n"); 
          *failure = -1;
          return;
        }
    }
  }

  npix = irapErode->nx*irapErode->ny;


  /*  Erode irap grids */
  for (n=0; n<nGrids; n++) {
    if (irap[n] != NULL) {
      for(pix=0;pix < npix; pix++){
        if(irapErode->grid[pix] != irapErode->missingcode && 
          irap[n]->grid[pix] != irap[n]->missingcode){
            if (irapErode->grid[pix] >irap[n]->grid[pix])
              irap[n]->grid[pix] = irapErode->grid[pix];
            else
              irap[n]->grid[pix] = irap[n]->grid[pix];
          }
        else{
          irap[n]->grid[pix] = irap[n]->missingcode;
        }
      }
    }
  }

  return;
}








/*F:irapgridConstGrid*

________________________________________________________________

irapgridConstGrid
________________________________________________________________

Name:		irapgridConstGrid
Syntax:		@irapgridConstGrid-syntax
Description:	Return a constant irapgrid.


Error message through *failure:\\
0 Ok \\
-2 Allocation error \\

Return value: Pointer to a new irap grid with a constant value. 

Author:		O.Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridConstGrid-syntax: */
struct irapgrid *irapgridConstGrid(double xmin, double ymin, 
  double xinc, double yinc,
  int nx, int ny, 
  double constant,
  int *failure)
  /*>irapgridConstGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  npix = nx*ny;

  rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
  if(rirap == NULL){
    /*    LogKit::writeLog("ERROR (irapConstGrid): Allocation error\n");  */
    *failure = -2;
    return rirap;
  }


  rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
  if(rirap->grid == NULL){
    /*    LogKit::writeLog("ERROR (irapConstGrid): Allocation error\n"); */
    *failure = -2;
    free(rirap);
    rirap = NULL;
    return rirap;
  }

  rirap->nx = nx;
  rirap->ny = ny;
  rirap->xmin = xmin;
  rirap->ymin = ymin;
  rirap->xinc = xinc;
  rirap->yinc = yinc;
  rirap->xmax = xmin + (nx-1)*xinc;
  rirap->ymax = ymin + (ny-1)*yinc;
  rirap->missingcode = IRAPMISSING;
  rirap->constValue = constant;

  /*  Assign value to the grid */
  for(pix=0;pix < npix; pix++){
    rirap->grid[pix] = constant;
  }
  rirap->bin = 0;
  return rirap;
}






/*F:irapgridMaxGrid*

________________________________________________________________

irapgridMaxGrid
________________________________________________________________

Name:		irapgridMaxGrid
Syntax:		@irapgridMaxGrid-syntax
Description: Calculate a new grid that is the maximum of two grids.
The two grids must be compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to be 
the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids can not add or undefined grids \\
-2 Allocation error\\
-3 Internal error \\

Return value: Pointer to a new irap grid that is the maximum of the input
grids.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridMaxGrid-syntax: */
struct irapgrid *irapgridMaxGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridMaxGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapMaxGrid): Pointer irap= NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapMaxGrid): Pointer irap->grid = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapMaxGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapMaxGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapMaxGrid): Allocation error\n");  */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
    {
      if (irap1->constValue >= irap2->constValue)
        rirap->constValue = irap1->constValue;
      else
        rirap->constValue = irap2->constValue;
    }
    else
      rirap->constValue = rirap->missingcode;


    /*  Maximum of irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          if(irap1->grid[pix] >= irap2->grid[pix]){
            rirap->grid[pix] = irap1->grid[pix];	
          }
          else {
            rirap->grid[pix] = irap2->grid[pix];	
          }
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    rirap->bin = irap1->bin;
    return rirap;
}


/*F:irapgridMinGrid*

________________________________________________________________

irapgridMinGrid
________________________________________________________________

Name:		irapgridMinGrid
Syntax:		@irapgridMinGrid-syntax
Description: Calculate a new grid that is the minimum of two grids.
The two grids must be compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to be 
the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids or undefined grids \\
-2 Allocation error\\
-3 Internal error \\

Return value: Pointer to a new irap grid that is the minimum of the input
grids.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridMinGrid-syntax: */
struct irapgrid *irapgridMinGrid(struct irapgrid *irap1,
struct irapgrid *irap2,
  int *failure)
  /*>irapgridMinGrid-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapMinGrid): Pointer irap = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapMinGrid): Pointer irap->grid = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapMinGrid): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapMinGrid): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapMinGrid): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
    {
      if (irap1->constValue <= irap2->constValue)
        rirap->constValue = irap1->constValue;
      else
        rirap->constValue = irap2->constValue;
    }
    else
      rirap->constValue = rirap->missingcode;

    /*  Minimum of irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          if(irap1->grid[pix] <= irap2->grid[pix]){
            rirap->grid[pix] = irap1->grid[pix];	
          }
          else {
            rirap->grid[pix] = irap2->grid[pix];	
          }
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }

    rirap->bin = irap1->bin;
    return rirap;
}


/*F:irapgridZoneThickness*

________________________________________________________________

irapgridZoneThickness
________________________________________________________________

Name:		irapgridZoneThickness
Syntax:		@irapgridZoneThickness-syntax
Description: Calculate a new grid that is the 
thickness of a zone after removing top and bottom erosion.
If the thickness is less than zero an error is returned in
'int *failure'.

$ The thickness is defined by:
$ \[
$    \mbox{Thickness} = \Delta Z - e_t(x,y) - e_b(x,y)
$ \]
$ where $e_t(x,y)$ and $e_b(x,y)$ is the top and bottom
$ erosion and $\Delta Z$ is the original uncompressed and
$ uneroded formation thickness.

The two grids must be compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to be 
the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids or undefined grids \\
-2 Allocation error\\
-3 Internal error \\

Warning message:\\
1 Erosion is larger than original thickness.\\

Return value: Pointer to a new irap grid that is the thickness after erosion.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridZoneThickness-syntax: */
struct irapgrid *irapgridZoneThickness(struct irapgrid *irap1,
struct irapgrid *irap2,
  double dz,
  int *failure)
  /*>irapgridZoneThickness-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;
  int nneg=0;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapZoneThickness): Pointer irap = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapZoneThickness): Pointer irap->grid = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapZoneThickness): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    npix = irap1->nx*irap1->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapZoneThickness): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapZoneThickness): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = irap1->nx;
    rirap->ny = irap1->ny;
    rirap->xmin = irap1->xmin;
    rirap->ymin = irap1->ymin;
    rirap->xinc = irap1->xinc;
    rirap->yinc = irap1->yinc;
    rirap->xmax = irap1->xmax;
    rirap->ymax = irap1->ymax;
    rirap->missingcode = irap1->missingcode;

    if (!equalRealDoubles(irap1->constValue,irap1->missingcode,5) &&
      !equalRealDoubles(irap2->constValue,irap2->missingcode,5))
      rirap->constValue = dz - irap1->constValue - irap2->constValue;
    else
      rirap->constValue = rirap->missingcode;

    /*  Calculate thickness */
    nneg=0;
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          rirap->grid[pix] = dz - irap1->grid[pix] - irap2->grid[pix];	
          if(rirap->grid[pix] < 0.0){
            nneg++;
          }
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    if(nneg>0){
      /*    LogKit::writeLog("ERROR: Formation thickness after erosion is negative\n"); */
      /*    LogKit::writeLog("       in %d nodes\n",nneg); */
      *failure = 1;
    }
    rirap->bin = irap1->bin;
    return rirap;
}



/*F:irapgridCheckConsistency*

________________________________________________________________

irapgridCheckConsistency
________________________________________________________________

Name:		irapgridCheckConsistency
Syntax:		@irapgridCheckConsistency-syntax
Description: Check that the four surfaces and the original thickness 'dz' 
defining are transformation from global coordinate system 
to simulation box coordinate system are consistent.


Side effects: The points of the surfaces that do not fullfill the consistency 
checks but with error less than the tolerance, these points are
corrected such that the surfaces returned fullfill the consistency.
If however the error is larger than the tolerance, error code is returned.

Return value:  0 Ok\\
-1 Error: Incompatible grids or undefined grids or internal error \\
-2 Error: Difference between bottom and top surface is negative \\
-3 Error:Thickness of zone after transforming to reservoir box coordinates \\
is negative which means that dz is too small relative to the \\
erosion surfaces. \\
-4 Error: Difference between bottom and top surface is positive but \\
thickness after transformation is 0.0 in some point whick means \\
that the erosion is too large or dz is too small. \\
-5 Error: Top erosion is negative at some points.	 \\
-6 Error: Bottom erosion is negative at some points. \\
-7 Warning: dz < max(bot -top) which means expansion not compaction \\
-8 Warning: There are nodes where thickness after transformation \\
to reservoir box coordinates is less than thickness \\
of the zone in world coordinate system. This means that \\
there are positions where the thickness will expand from \\
reservoir box coordinate system to the thickness as it \\
is observed in world coordinate system.\\

Author:		Oddvar Lia, NR
Date: Nov 1993, Updated: Feb 1994

________________________________________________________________
*/

/*<irapgridCheckConsistency-syntax: */
int irapgridCheckConsistency(struct irapgrid *top,
struct irapgrid *bot,
struct irapgrid *topEro,
struct irapgrid *botEro,
  double dz)
  /*>irapgridCheckConsistency-syntax: */
{
  int pix,npix;
  int expansion=0;
  int nofNotMissings;
  int negative_untransformed = 0;
  int negative_transformed = 0;
  int untransformed_transformed = 0;
  int negative_eTop = 0;
  int negative_eBot = 0;
  int missingNonCorrespond = 0;

  double thickness_transformed,thickness_untransformed;
  double maxthickness;
  double tolerance;
  double eTop,eBot,zTop,zBot;
  double average;

  if(top == NULL || bot == NULL){
    LogKit::writeLog("ERROR (irapgridCheckConsistency): Pointers are NULL\n"); 
    return -1;
  }

  if(top->grid == NULL || bot->grid == NULL){
    LogKit::writeLog("ERROR (irapgridCheckConsistency): Pointers are NULL\n"); 
    return -1;
  }

  if(top->nx != bot->nx || top->ny != bot->ny ||
    top->xmin != bot->xmin || top->ymin != bot->ymin ||
    top->xinc != bot->xinc || top->yinc != bot->yinc){
      return  -1;
    }

    if(top->nx != topEro->nx || top->ny != topEro->ny ||
      top->xmin != topEro->xmin || top->ymin != topEro->ymin ||
      top->xinc != topEro->xinc || top->yinc != topEro->yinc){
        return  -1;
      }

      if(top->nx != botEro->nx || top->ny != botEro->ny ||
        top->xmin != botEro->xmin || top->ymin != botEro->ymin ||
        top->xinc != botEro->xinc || top->yinc != botEro->yinc){
          return  -1;
        }

        npix = top->nx*top->ny;

        tolerance = REL_TOLERANCE*dz;
        tolerance = MINIM(tolerance,ABS_TOLERANCE);

        /*  Calculate thickness */
        maxthickness = -10000000;
        for(pix=0;pix < npix; pix++){
          zTop = top->grid[pix];
          zBot = bot->grid[pix];
          eTop = topEro->grid[pix];
          eBot = botEro->grid[pix];
          if(zTop != top->missingcode && 
            zBot != bot->missingcode &&
            eTop != topEro->missingcode &&
            eBot != botEro->missingcode)
          {
            /* Check that erosion surfaces are non negative */

            if(eTop < -tolerance){
              negative_eTop = 1;
            }
            else if(eTop < 0.0){
              topEro->grid[pix] = 0.0;
              eTop = 0.0;
            }

            if(eBot < -tolerance){
              negative_eBot = 1;
            }
            else if(eBot < 0.0){
              botEro->grid[pix] = 0.0;
              eBot = 0.0;
            }

            /* Check that dz >= (top erosion + bottom erosion)  */

            thickness_transformed = dz - eTop - eBot;	
            if(thickness_transformed < -tolerance){
              negative_transformed = 1;
            }
            else if(thickness_transformed < 0.0){
              average = (dz - eBot + eTop)/2.0;
              topEro->grid[pix] = average;
              botEro->grid[pix] = dz - average;
              eTop = average;
              eBot = dz - average;
              thickness_transformed = 0.0;
            }

            /* Check that thickness of zone is non-negative */

            thickness_untransformed = zBot - zTop;
            if(thickness_untransformed < -tolerance){
              negative_untransformed = 1;
            }
            else if(thickness_untransformed < 0){
              thickness_untransformed = 0.0;
              average = (zTop + zBot)/2.0;
              bot->grid[pix] = average;
              top->grid[pix] = average;
            }

            if(thickness_untransformed > 0.0 && thickness_transformed == 0){
              untransformed_transformed = 1;
            }
            if(thickness_transformed < thickness_untransformed &&
              !(equalRealDoubles(thickness_transformed,thickness_untransformed,5))){
                expansion = 1;
              }

              if(thickness_untransformed > maxthickness){
                maxthickness = thickness_untransformed;
              }

          }/* end when no missing codes */
          else /*at least one node contains missing values*/
          {
            nofNotMissings = 0;
            if(top->filename == NULL)
            {
              top->grid[pix] = top->missingcode;
            }
            else if (zTop != top->missingcode)
            {
              nofNotMissings++;
              top->grid[pix] = top->missingcode;
            }
            if(bot->filename == NULL)
            {
              bot->grid[pix] = bot->missingcode;
            }
            else if (zBot != bot->missingcode)
            {
              nofNotMissings++;
              bot->grid[pix] = bot->missingcode;
            }
            if(topEro->filename == NULL)
            {
              topEro->grid[pix] = topEro->missingcode;
            }
            else if (eTop != topEro->missingcode)
            {
              nofNotMissings++;
              topEro->grid[pix] = topEro->missingcode;
            }
            if(botEro->filename == NULL)
            {
              botEro->grid[pix] = botEro->missingcode;
            }
            else if (eBot != botEro->missingcode)
            {
              nofNotMissings++;
              botEro->grid[pix] = botEro->missingcode;
            }

            if(nofNotMissings > 0)
            {
              missingNonCorrespond = 1; 
            }
          }
        }

        if(missingNonCorrespond == 1)
        {
          LogKit::writeLog("%s\n%s\n%s%s\n",
            "WARNING (irapgridCheckConsistency):",
            "- Missing codes in maps do not correspond.",
            "There are at least one node with missingcode in some ",
            "maps and value in others.");
        }

        if(negative_eTop){
          return -5;
        }
        if(negative_eBot){
          return -6;
        }
        if(negative_untransformed){
          return -2;
        }
        if(negative_transformed){
          return -3;
        }
        if(untransformed_transformed){
          return -4;
        }

        if(maxthickness > dz && !equalRealDoubles(maxthickness,dz,4)){
          return -7;
        }
        if(expansion){
          return -8;
        }

        return 0;
}


/*F:irapgridNotIdentical*

________________________________________________________________

irapgridNotIdentical
________________________________________________________________

Name:		irapgridNotIdentical
Syntax:		@irapgridNotIdentical-syntax
Description:    Check if two irapgrids are identical. 
Return value:   0 if identical grids, 1 otherwise
Author:		Anne-Lise Hektoen, NR
Date:           October 1994
________________________________________________________________

*/

/*<irapgridNotIdentical-syntax: */
int irapgridNotIdentical(struct irapgrid *irap1,
struct irapgrid *irap2)
  /*>irapgridNotIdentical-syntax: */
{ int i, pix,npix; 
double value1, value2;
double constant1 = True;
double constant2 = True;

/* Grids may be constant, but covering different area --> still identical */

i = 0;
value1 = irap1->constValue;
if (value1 == IRAPMISSING)
{ /* May still be constant, but may have missing values */
  while (value1 == IRAPMISSING)
    value1 = irap1->grid[i++];

  npix = irap1->nx*irap1->ny;
  for (pix=i; pix<npix; pix++)
  {
    if (!equalRealDoubles(irap1->grid[pix],IRAPMISSING,5))
    {
      if (!(equalRealDoubles(value1,irap1->grid[pix],5)))
        constant1 = False;
    }
  }
}

i=0;
value2 = irap2->constValue;
if (value2 == IRAPMISSING)
{
  while (value2 == IRAPMISSING)
    value2 = irap2->grid[i++];
  npix = irap2->nx*irap2->ny;
  for (pix=i; pix<npix; pix++)
  {
    if (!equalRealDoubles(irap2->grid[pix],IRAPMISSING,5))
    {
      if (!(equalRealDoubles(value2,irap2->grid[pix],5)))
        constant2 = False;
    }
  }
}

if (constant1 && constant2)
{ /* Both grid constants - check value */
  if (!equalRealDoubles(value1,value2,5))
    return 1;
}
else
{ /* Not constant grids, require that grid cover same area 
  and have same grid density */
  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc)
    return 1;

  npix = irap1->nx*irap1->ny;
  for(pix=0;pix < npix; pix++)
  {
    if (!equalRealDoubles(irap1->grid[pix],irap2->grid[pix],5))
      return 1;
  }
}

return 0;
}	/* end of irapgridNotIdentical */



/*F:irapgridTopErosionThickness*

________________________________________________________________

irapgridTopErosionThickness
________________________________________________________________

Name:		irapgridTopErosionThickness
Syntax:		@irapgridTopErosionThickness-syntax
Description: Calculate a new grid that is the 
thickness of a top erosion map.
The top erosion thickness map is defined by
the difference between an original constant zone
thickness 'dz' minus the difference between bottom and top
surface maps for the zone. The erosion thickness calculated
is the rescaled by a factor between 0 and 1. If this factor
'erosionfactor' is 1.0, the returned erosion thickness map
is constructed such that all thickness variation of the
zone is taken out as erosion, but if the factor is 0.0
no thickness variation is due to erosion all all thickness
variation is due to compaction and the returned erosion map
is 0.0. The original zone thickness specified must at least
be as large as the maximum difference between the bottom
and top surface for the zone.


$ The erosion thickness is defined by:
$ \[
$    \mbox{Thickness} = \alpha(\Delta Z - (z_b(x,y) - z_t(x,y)))
$ \]
$ where $z_t(x,y)$ and $z_b(x,y)$ is the top and bottom
$ surfaces and $\Delta Z$ is the original uncompressed and
$ uneroded formation thickness. The factor $\alpha$ is the
$ factor {\em erosionfactor} defining the degree of erosion and
$ compaction.

The two grids must be compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to be 
the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids or undefined grids \\
-2 Allocation error\\
-3 Internal error \\
-4 dz <= 0.0 or erosionfactor is not between 0.0 and 1.0 \\
-5 top grid is deeper than bottom grid

Warning message:\\
1  Zone thickness originally is less than the difference between
top and bottom surfaces.

Return value: Pointer to a new irap grid that is the erosion thickness 
of top erosion.

Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridTopErosionThickness-syntax: */
struct irapgrid *irapgridTopErosionThickness(struct irapgrid *iraptop,
struct irapgrid *irapbot,
  double dz, 
  double erosionfactor,
  int *failure)
  /*>irapgridTopErosionThickness-syntax: */
{
  struct irapgrid *rirap=NULL;
  int pix,npix;
  int nneg=0;
  int nnegthick=0;
  double diff;

  *failure = 0;

  if(iraptop == NULL || irapbot == NULL){
    LogKit::writeLog("ERROR (irapTopErosionThickness): Pointer irap = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(iraptop->grid == NULL || irapbot->grid == NULL){
    LogKit::writeLog("ERROR (irapTopErosionThickness): Pointer irap->grid = NULL\n"); 
    *failure = -3;
    return rirap;
  }

  if(iraptop->nx != irapbot->nx || iraptop->ny != irapbot->ny ||
    iraptop->xmin != irapbot->xmin || iraptop->ymin != irapbot->ymin ||
    iraptop->xinc != irapbot->xinc || iraptop->yinc != irapbot->yinc){
      /*    LogKit::writeLog("ERROR (irapTopErosionThickness): Incompatible irap grids\n"); */
      *failure = -1;
      return rirap;
    }

    if(erosionfactor < 0.0 || erosionfactor > 1.0 || dz <= 0.0){
      *failure = -4;
      return rirap;
    }

    npix = iraptop->nx*iraptop->ny;

    rirap = (struct irapgrid *) calloc(1,sizeof(struct irapgrid));
    if(rirap == NULL){
      /*    LogKit::writeLog("ERROR (irapTopErosionThickness): Allocation error\n"); */
      *failure = -2;
      return rirap;
    }


    rirap->grid = (double *) malloc(((unsigned)npix)*sizeof(double));
    if(rirap->grid == NULL){
      /*    LogKit::writeLog("ERROR (irapTopErosionThickness): Allocation error\n"); */
      *failure = -2;
      free(rirap);
      rirap = NULL;
      return rirap;
    }

    rirap->nx = iraptop->nx;
    rirap->ny = iraptop->ny;
    rirap->xmin = iraptop->xmin;
    rirap->ymin = iraptop->ymin;
    rirap->xinc = iraptop->xinc;
    rirap->yinc = iraptop->yinc;
    rirap->xmax = iraptop->xmax;
    rirap->ymax = iraptop->ymax;
    rirap->missingcode = iraptop->missingcode;
    rirap->constValue = rirap->missingcode;


    /*  Calculate thickness */
    nneg=0;
    nnegthick=0;
    for(pix=0;pix < npix; pix++){
      if(iraptop->grid[pix] != iraptop->missingcode && 
        irapbot->grid[pix] != irapbot->missingcode){
          diff = irapbot->grid[pix] - iraptop->grid[pix];
          if(diff < 0.0){
            nnegthick++;
          }
          else {
            rirap->grid[pix] = 
              erosionfactor*(dz - (irapbot->grid[pix] - iraptop->grid[pix]));
            if(rirap->grid[pix] < 0.0){
              nneg++;
              rirap->grid[pix] = 0.0;
            }
          }
        }
      else{
        rirap->grid[pix] = rirap->missingcode;
      }
    }
    if(nnegthick>0){
      /*    LogKit::writeLog("ERROR: Thickness of zone is negative\n"); */
      /*    LogKit::writeLog("       in %d nodes\n",nnegthick); */
      *failure = -5;
    }
    if(nneg>0){
      /*    LogKit::writeLog("ERROR: Erosion thickness  is truncatet from negative to 0\n"); */
      /*    LogKit::writeLog("       in %d nodes\n",nneg); */
      *failure = 1;
    }
    rirap->bin = iraptop->bin;
    return rirap;
}



/*F:irapgridCheckCrossing*

________________________________________________________________

irapgridCheckCrossing
________________________________________________________________

Name:		irapgridCheckCrossing
Syntax:		@irapgridCheckCrossing-syntax
Description: Check that two grids 'grid1' and 'grid2' do
not cross each other.
If the grids cross each other, 
the function returns an 'updated' version of each grid.

The updated version of 'grid1' is equal to the original value of
'grid1' where 
this grid has value less than 'grid2' 
and is equal to the average
of 'grid1' and 'grid2' where 'grid1' is larger than 'grid2'.

The updated version of 'grid2' is equal to the original value of
'grid2' where this grid has value larger than than 'grid1' 
and is equal to the average
of 'grid1' and 'grid2' where 'grid1' is larger than 'grid2'.

The two grids must be compatible which means they are
defined to have the same number of nodes,xmin,ymin,xinc and yinc.
The missing code of the first grid is defined to be 
the missingcode
for the new grid. A missing code is returned  in all nodes
where at least one of the two grids have  missing code.

Error message through *failure:\\
0 Ok\\
-1 Incompatible grids or undefined grids \\
-2 Allocation error\\
-3 Internal error \\

Warning message throug *failure:\\
1 if the grids cross each other.


Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridCheckCrossing-syntax: */
void irapgridCheckCrossing(struct irapgrid *irap1,
struct irapgrid *irap2, 
  int *failure)
  /*>irapgridCheckCrossing-syntax: */
{
  int pix,npix;

  *failure = 0;

  if(irap1 == NULL || irap2 == NULL){
    LogKit::writeLog("ERROR (irapCheckCrossing): Pointer irap = NULL\n"); 
    *failure = -3;
    return;
  }

  if(irap1->grid == NULL || irap2->grid == NULL){
    LogKit::writeLog("ERROR (irapCheckCrossing): Pointer irap->grid = NULL\n"); 
    *failure = -3;
    return;
  }

  if(irap1->nx != irap2->nx || irap1->ny != irap2->ny ||
    irap1->xmin != irap2->xmin || irap1->ymin != irap2->ymin ||
    irap1->xinc != irap2->xinc || irap1->yinc != irap2->yinc){
      /*    LogKit::writeLog("ERROR (irapCheckCrossing): Incompatible irap grids\n"); */
      *failure = -1;
      return;
    }

    npix = irap1->nx*irap1->ny;


    /*  Check crossing of irap grids */
    for(pix=0;pix < npix; pix++){
      if(irap1->grid[pix] != irap1->missingcode && 
        irap2->grid[pix] != irap2->missingcode){
          if(irap1->grid[pix] > irap2->grid[pix]){
            irap1->grid[pix] = 
              (irap1->grid[pix] + irap2->grid[pix])/2.0;	
            irap2->grid[pix] = irap1->grid[pix];
            *failure = 1;
          }
        }
    }

    return;
}


/*F:irapgridMinMaxValue*

________________________________________________________________

irapgridMinMaxValue
________________________________________________________________

Name:		irapgridMinMaxValue
Syntax:		@irapgridMinMaxValue-syntax
Description: Calculate minimum and maximum value of all defined node values.
Returns the minimum and maximum value of the grid.


Author:		Oddvar Lia, NR
Date: Nov 1993

________________________________________________________________

*/

/*<irapgridMinMaxValue-syntax: */
void irapgridMinMaxValue(struct irapgrid *irap1,
                         double *minvalue, double *maxvalue)
                         /*>irapgridMinMaxValue-syntax: */
{
  int pix,npix;
  double minval,maxval;

  if(irap1 == NULL){
    LogKit::writeLog("ERROR (irapMinMaxValue): Pointer irap = NULL\n"); 
    return;
  }

  if(irap1->grid == NULL){
    LogKit::writeLog("ERROR (irapMinMaxValue): Pointer irap->grid = NULL\n"); 
    return;
  }


  npix = irap1->nx*irap1->ny;

  /*  Minimum value of irap grids */
  minval = 1.0e+30;
  maxval = - minval;
  for(pix=0;pix < npix; pix++){
    if(irap1->grid[pix] != irap1->missingcode ){
      if(irap1->grid[pix] < minval){
        minval = irap1->grid[pix];	
      }
      if(irap1->grid[pix] > maxval){
        maxval = irap1->grid[pix];	
      }
    }
  }
  *maxvalue = maxval;
  *minvalue = minval;
  return;
}




/*F:firapgridReadpt*

________________________________________________________________

firapgridReadpt
________________________________________________________________

Name:		firapgridReadpt
Syntax:		@firapgridReadpt-syntax
Description: 	Function for reading grid from file specified 
by a pointer, returning pointer to the grid 
struct and an error indicator. File format is ascii.
The grid is floating point.

int *failure  has the value:\\
0: OK\\
-1: error file is not open\\
-2: error when reading file\\
-3: error in allocating space for grid\\

Side effects:
Return value: Pointer to firapgrid, NULL if the file is cannot be read.
Author:		Oddvar Lia, NR
Date:  March 1994
________________________________________________________________

*/

/*<firapgridReadpt-syntax: */
struct firapgrid *firapgridReadpt(FILE *file, int *failure) 
  /*>firapgridReadpt-syntax: */
{
  int nread,nreadSum,i,npix;
  struct firapgrid *rIrap;
  float dummy[4];

  *failure = 0;
  if(file == NULL){
    *failure = -1;
    return NULL;
  }

  rIrap = (struct firapgrid *) malloc(sizeof(struct firapgrid));
  if (rIrap == NULL)
  {*failure = -3;
  return NULL;
  };

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE reading irap grid\n");

  nread = fscanf(file,"%d %d %e %e",&(rIrap->nx),&(rIrap->ny),
    &(dummy[0]),&(dummy[1]));
  if(rIrap->nx <= 0 || rIrap->ny <= 0) {
    /*     LogKit::writeLog("ERROR: Number of nodes in irapgrid are < = 0 "); */
    *failure = -2;
    return NULL;
  }
  npix = rIrap->nx * rIrap->ny;
  rIrap->xinc = (double) dummy[0];
  rIrap->yinc = (double) dummy[1];
  if (nread != 4)
  {*failure = -2;
  free (rIrap);
  rIrap = NULL;
  return NULL;
  };
  rIrap->grid = (float *) malloc(((unsigned) npix)*sizeof(float));
  if (rIrap->grid == NULL)
  {
    *failure = -3;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };
  nread = fscanf(file,"%lf %lf %lf %lf",&(rIrap->xmin),&(rIrap->xmax),
    &(rIrap->ymin),&(rIrap->ymax));
  if (nread != 4)
  {
    *failure = -2;
    free((void *)rIrap->grid);
    rIrap->grid = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };

  nreadSum = 0;
  for (i = 0; i < npix; i++)
  {
    nread = fscanf(file,"%f ",&(rIrap->grid[i]));
    nreadSum +=nread;
  };

  if (nreadSum != npix)
  {
    *failure = -2;
    free((void *)rIrap->grid);
    rIrap->grid = NULL;
    free(rIrap);
    rIrap = NULL;
    return NULL;
  };

  rIrap->constValue = IRAPMISSING;
  rIrap->missingcode = IRAPMISSING;
  rIrap->bin = 0;
  return rIrap;
}




/*F:firapgridWritept*

________________________________________________________________

firapgridWritept
________________________________________________________________

Name:		firapgridWritept
Syntax:		@firapgridWritept-syntax
Description: 	Function for writing grid to file specified by a 
pointer, returning an error indicator. 
The file format is ascii.
The grid is floating point.

Return value:  0 if writing succesfull\\
-1 if file is not opened\\
-2 if error under writing\\

Author:		Oddvar Lia, NR
Date:   March 1994
________________________________________________________________

*/

/*<firapgridWritept-syntax: */
int firapgridWritept(FILE *file, struct firapgrid *irap)
/*>firapgridWritept-syntax: */
{
  int i,nwritten,n,npix;

  if (file == NULL) return -1;


  if(irap == NULL){
    LogKit::writeLog("ERROR (firapgridWritept): Pointer irap = NULL\n"); 
    return -2;
  }

  if(irap->grid == NULL){
    LogKit::writeLog("ERROR (firapgridWritept): Pointer irap->grid = NULL\n"); 
    return -2;
  }

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE writing irap grid\n");


  npix = irap->nx*irap->ny;
  nwritten = fprintf(file,"%d %d %f %f \n",irap->nx,irap->ny,
    irap->xinc,irap->yinc);
  if (nwritten < 0)
    return -2; 

  nwritten = fprintf(file,"%f %f %f %f \n",irap->xmin,irap->xmax,
    irap->ymin,irap->ymax);
  if (nwritten < 0) return -2; 

  n = 0;
  for ( i = 0; i < npix; i++) {
    n++;
    nwritten = fprintf(file,"%f ",irap->grid[i]);
    if(nwritten < 0) return -2;
    if(n == 6){
      n = 0;
      fprintf(file,"\n");
    }
  }
  fprintf(file,"\n");

  return 0; 
}





/*F:firapgridReadBinpt*

________________________________________________________________

firapgridReadBinpt
________________________________________________________________

Name:		firapgridReadBinpt
Syntax:		@firapgridReadBinpt-syntax
Description: 	Function for reading grid from file specified 
by a pointer, returning 
pointer to the grid struct and an error indicator. 
File format is binary.
The binary format is the format 'firapgridWriteBinpt' uses.
The grid is floating points.

int *failure has the value:\\
0: OK\\
-1: file is not open\\
-2: error when reading file\\
-3: error in allocating space for grid\\

Return value: Pointer to an array with the irap structure. If read error, the
pointer is NULL.

Author:		Oddvar Lia, NR
Date:  Nov 1993
________________________________________________________________

*/

/*<firapgridReadBinpt-syntax: */
struct firapgrid *firapgridReadBinpt(FILE *file, int *failure) 
  /*>firapgridReadBinpt-syntax: */
{
  int nread,npix;
  struct firapgrid *rIrap;
  float dummy[4];
  char c;

  *failure = 0;
  if (file == NULL) 
  {*failure = -1;
  return NULL;
  };
  rIrap = (struct firapgrid *) malloc(sizeof(struct firapgrid));
  if (rIrap == NULL)
  {*failure = -3;
  return NULL;
  };

  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE reading irap grid\n");

  nread = fscanf(file,"%d %d %e %e",
    &(rIrap->nx),&(rIrap->ny),&(dummy[0]),&(dummy[1]));
  if(rIrap->nx <= 0 || rIrap->ny <= 0) {
    LogKit::writeLog("%s \n","WARNING: number of nodes in irapgrid are < = 0 ");
    *failure = -2;
    return NULL;
  }

  rIrap->xinc = (double) dummy[0];
  rIrap->yinc = (double) dummy[1];
  if (nread != 4)
  {*failure = -2;
  free (rIrap);
  rIrap = NULL;
  return NULL;
  };


  rIrap->grid = (float *) malloc(((unsigned)(rIrap->nx * rIrap->ny))* 
    sizeof(float));
  if (rIrap->grid == NULL)
  {*failure = -3;
  free(rIrap);
  rIrap = NULL;
  return NULL;
  };


  nread = fscanf(file,"%lf %lf %lf %lf",&(rIrap->xmin),&(rIrap->xmax),
    &(rIrap->ymin),&(rIrap->ymax));
  if (nread != 4)
  {*failure = -2;
  free((void *)rIrap->grid);
  rIrap->grid = NULL;
  free(rIrap);
  rIrap = NULL;
  return NULL;
  };

  /* Get new line */
  for (c = char(fgetc(file)); !((c == '\n') || (c == EOF) || (c == '\377'));
    c = char(fgetc(file)));

    npix = rIrap->nx*rIrap->ny;

  nread = fread(&(rIrap->grid[0]),sizeof(float),npix,file);
  if(nread != npix){
    *failure = -2;
    free(&(rIrap->grid[0]));
    free(rIrap);
    return NULL;
  }

  rIrap->constValue = IRAPMISSING;
  rIrap->missingcode = IRAPMISSING;
  rIrap->bin = 1;
  return rIrap;
}




/*F:firapgridWriteBinpt*

________________________________________________________________

firapgridWriteBinpt
________________________________________________________________

Name:		firapgridWriteBinpt
Syntax:		@firapgridWriteBinpt-syntax
Description: 	Function for writing grid to file specified by a pointer, 
returning  an error indicator. The file format is binary.


Return value:  0 if writing succesfull\\
-1 if file is not open\\
-2 if error under writing\\

Author:		Oddvar Lia, NR
Date: Nov 1993
________________________________________________________________

*/

/*<firapgridWriteBinpt-syntax: */
int firapgridWriteBinpt(FILE *file, struct firapgrid *irap)
/*>firapgridWriteBinpt-syntax: */
{
  int npix,nwritten;

  if (file == NULL) 
    return -1;



  if(IRAPPRINT)
    LogKit::writeLog("MESSAGE writing irap grid\n");

  if(irap == NULL){
    LogKit::writeLog("ERROR (firapgridWriteBinpt): Pointer irap = NULL\n"); 
    return -2;
  }

  if(irap->grid == NULL){
    LogKit::writeLog("ERROR (firapgridWriteBinpt): Pointer irap->grid = NULL\n"); 
    return -2;
  }


  nwritten = fprintf(file,"%d %d %f %f \n",irap->nx,irap->ny,
    irap->xinc,irap->yinc);
  if (nwritten < 0)
    return -2; 

  nwritten = fprintf(file,"%f %f %f %f \n",irap->xmin,irap->xmax,
    irap->ymin,irap->ymax);
  if (nwritten < 0)
    return -2; 

  npix = irap->nx*irap->ny;
  nwritten = fwrite(&(irap->grid[0]),sizeof(float),npix,file);
  if(nwritten != npix){
    return -2;
  }

  return 0; 
}

int
get2DGridInfo(
              char *fileName,
              int *ncol,
              int *nrow,
              double * x1,
              double * x2,
              double * y1,
              double * y2,
              double * xinc,
              double * yinc)
{
  FILE *	file;
  int	status;
  int	err = 1;
  char	text[MAX_STRING];

  file = fopen(fileName, "r");
  if (!file)
  {
    return(1);
  }

  status = fscanf(file, "%s", text);
  if (status != 1) goto error_exit;

  status = fscanf(file, "%d %d %lf %lf", ncol, nrow, xinc, yinc);
  if (status != 4) goto error_exit;

  status = fscanf(file, "%lf %lf %lf %lf", x1, x2, y1, y2);
  if (status != 4) goto error_exit;

  err = 0;

error_exit:

  fclose(file);
  return(err);
}

void
fget2dgridinfo(
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
               int fileNameL)
{
  char *	file;
  int	i;

  file = (char *) calloc(fileNameL+1, sizeof(char));
  strncpy(file, fileName, fileNameL);

  for (i = 0; i < fileNameL; i++)
  {
    if (file[i] == ' ')
    {
      file[i] = 0;
      break;
    }
  }

  file[fileNameL] = 0;

  *ierr = get2DGridInfo(file, ncol, nrow, x1, x2, y1, y2, xinc, yinc);

  free(file);
  return;
}
