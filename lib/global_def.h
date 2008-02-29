/*
(C) Copyright ODIN Reservoir Software & Services AS
P.O.Box 585, Madla, N-4040 Hafrsfjord, Norway

PURPOSE 
Define macros 
*/

#if !defined(GLOBAL_DEF_H)
#define GLOBAL_DEF_H 1 
#ifndef PI 
#define PI  3.14159265358979323846
#endif

#ifndef RMISSING
#define RMISSING -99999.000
#endif

#ifndef WELLMISSING
#define WELLMISSING -999.00000
#endif


#define IMISSING -99999
#define CHNMISSING -888 /* Must differ from IMISSING!! */

#define DUMMYZONE 999999 /* Used in lib_obs/GenWell to identify zone border
where the zonenumber above/below is unknown    */

#define MAX_PATH  2000  /* Max length of filenames */

#define MAX_STRING 2400  /* Max length of string parameters */

#define MAXIM(A,B) ((A) > (B) ? (A) : (B))
#define MINIM(A,B) ((A) < (B) ? (A) : (B))
#define MINMAX(a,b,c) ((a) > (b) ? (a) : ((b) > (c) ? (c) : (b)))

#define MALLOC(type) (type *) calloc(1,sizeof(type))
#define CALLOC(size,type) (type *) calloc(size + 1,sizeof(type))
#define REALLOC(p,size,type) p = (type *) realloc(p,(size + 1) * sizeof(type))
#define FREE(p) {free(p);p=NULL;}

#ifndef True
#define True 1
#endif

#ifndef False
#define False 0
#endif

#endif
