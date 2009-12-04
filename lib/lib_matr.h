#ifndef LIB_MATR_H
#define LIB_MATR_H

#include "fft/include/fftw.h"
#ifdef __cplusplus
extern "C"
{
#endif
  extern int lib_matrCholR(int i_dim, double **x_mat);
  extern void lib_matrAxeqbR(int i_dim, double **i_mat, double *x_vec);
  extern void lib_matrAXeqBMatR(int i_dim, double ** i_mat, double **x_mat, int n);
  extern int lib_matrCholCpx(int i_dim, fftw_complex **x_mat);
  extern void lib_matrAxeqbCpx(int i_dim, fftw_complex **i_mat,fftw_complex *x_vec);
  extern void lib_matrAXeqBMatCpx(int i_dim, fftw_complex **i_mat, fftw_complex **x_mat, int n);
  extern void lib_matrLXeqBMatCpx(int i_dim, fftw_complex **i_mat, fftw_complex **x_mat, int n);
  extern void lib_matrLXeqMatR(int i_dim, double **i_mat, double **x_mat, int n);

  extern void lib_matrProdCholVec(int n, fftw_complex ** mat, fftw_complex * vec);

  extern void lib_matrProdScalMatCpx(fftw_complex s, fftw_complex **mat,int n1, int n2);
  extern void lib_matrProdScalVecCpx(fftw_complex s, fftw_complex *vec,int n);
  extern void lib_matrProdScalVecRCpx(float s, fftw_complex *vec,int n);

  extern void lib_matrProdCpx(fftw_complex **mat1, fftw_complex **mat2, int n1, int n2, int n3, fftw_complex **outmat);
  extern void lib_matrProdAdjointCpx(fftw_complex **mat1, fftw_complex **mat2, int n1, int n2, int n3, fftw_complex **outmat);
  extern void lib_matrProdCpxR(fftw_complex **mat1, float **mat2, int n1, int n2, int n3, fftw_complex **outmat);
  extern void lib_matrProdRCpx(float **mat1, fftw_complex **mat2, int n1, int n2, int n3, fftw_complex **outmat);
  extern void lib_matrProdDiagCpxR(fftw_complex *mat1, float **mat2, int n1,int n2, fftw_complex **outmat);
  extern void lib_matrProd2Cpx(fftw_complex **mat1, fftw_complex **mat2, int n1, int n2, fftw_complex **outmat);
  extern void lib_matrProdMatVecCpx(fftw_complex **mat, fftw_complex *vec, int n1, int n2, fftw_complex *outvec);
  extern void lib_matrProdMatRVecCpx(double **mat, fftw_complex *vec, int n1, int n2, fftw_complex *outvec);
  extern void lib_matrProdAdjointMatVecCpx(fftw_complex **mat, fftw_complex *vec, int n1, int n2, fftw_complex *outvec);

  extern void lib_matrAddMatCpx(fftw_complex **x, int n1, int n2, fftw_complex **y);
  extern void lib_matrSubtMatCpx(fftw_complex **x, int n1, int n2, fftw_complex **y);
  extern void lib_matrAddVecCpx(fftw_complex *x, int n,fftw_complex *y);
  extern void lib_matrSubtVecCpx(fftw_complex *x, int n, fftw_complex *y);

  extern void lib_matrAdjoint(fftw_complex **mat, int n1, int n2,fftw_complex **outmat);
  extern void lib_matrConj(fftw_complex **mat, int n1, int n2, fftw_complex **outmat);
  extern void lib_matrCopyCpx(fftw_complex **mat, int n1, int n2, fftw_complex **outmat);

  extern void lib_matrAddMat(double **x, int n1, int n2, double **y);
  extern void lib_matrSubtMat(double **x, int n1, int n2, double **y);
  extern void lib_matrCopy(double **mat, int n1, int n2, double **outmat);

  extern void lib_matrFillOnesVecCpx(fftw_complex* v1, int n1);
  extern void lib_matrFillValueVecCpx(fftw_complex value,fftw_complex* v1, int n1); 

  extern void lib_matrPrintR(float ** mat, int n1, int n2);
  extern void lib_matrPrintCpx(fftw_complex ** mat, int n1, int n2);
  extern void  lib_matr_eigen(double **i_mat,int i_n,double **o_eigvec,double *o_eigval,int *o_error);
  extern void lib_matr_prod(double **i_mat1, double **i_mat2,int i_n1,int i_n2,int i_n3,double **o_mat);
  extern void lib_matrTranspose(double **mat, int n1, int n2, double **outmat);
  extern void lib_matrLtXeqBR(int i_dim, double **i_mat, double **x_mat, int n);
  extern void lib_matr_sort3x3(double *eigenval, double **eigenvec);

  extern void lib_matrDump(const char * fName, double ** mat, int n1, int n2);
#ifdef __cplusplus
}
#endif

#endif
