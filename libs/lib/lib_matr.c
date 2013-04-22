/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include "lib/lib_matr.h"


/*#include <moheres/moheres_lib/utl.h>*/
static void lib_matr_tred2(int i_n, double **x_a, double o_d[], double o_e[]);
static void lib_matr_tqli(double i_e[], int i_n, double **x_z, double x_d[],
                 int *o_error);

/* Module spesific definition */
# define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
# define TOL 1e-20
int lib_matrCholR(
                  int i_dim, /*The dimension of the input matrix. */
                  double **x_mat) /*The matrix beeing decomposed on input, on output the lower
                                  triangel contains the L-matrix. */
                                  /*FUNC********************************************************************

                                  DESCRIPTION:

                                  Factorizes the positive definite, symmetric matrix 'x_mat' into
                                  L * L(transposed), when L is a lower triangular matrix. L is stored in the
                                  lower part of 'x_mat' on output.

                                  HOW TO USE THE FUNCTION:
                                  lib_matr_CholR(i_dim, x_mat);

                                  SIDE-EFFECTS: The input matrix 'x_mat' is destroyed on output.

                                  RETURN VALUE: 0 if everything O.K.
                                  1 if illegal matrix

                                  ************************************************************************/
{
  int l_i, l_j, l_k, r_value = 0;
  double l_r;

  for (l_i=0; l_i < i_dim; l_i++)
    if (x_mat[l_i][l_i] <= TOL) {
      /*
        fprintf(stderr, "Matrix to be Cholesky factorized not valid! %s\n",
        "Diagonal-element near 0.0");
        lib_matrPrintR( (float**)x_mat, i_dim, i_dim);
      */
      r_value = 1;
      return(r_value);
    } else {
      for (l_j=0; l_j < l_i; l_j++) {
        l_r = 0.0;
        for (l_k = 0; l_k < l_j; l_k++)
          l_r += x_mat[l_i][l_k] * x_mat[l_j][l_k];
        x_mat[l_i][l_j] = (x_mat[l_i][l_j] - l_r) / x_mat[l_j][l_j];
      }
      l_r = 0.0;
      for (l_k=0; l_k < l_i; l_k++)
        l_r += x_mat[l_i][l_k] * x_mat[l_i][l_k];
      l_r = x_mat[l_i][l_i] - l_r;
      if (l_r <= TOL) {
        /*
          fprintf(stderr, "Matrix to be Cholesky factorized not valid !\n");
        */
        r_value = 1;
        return(r_value);
      }
      x_mat[l_i][l_i] = sqrt(l_r);
    }
    return(r_value);
}

int lib_matrCholCpx(
                    int i_dim, /*The dimension of the input matrix. */
                    fftw_complex **x_mat) /*The matrix beeing decomposed on input, on output the lower
                                          triangel contains the L-matrix. */
                                          /*FUNC********************************************************************

                                          DESCRIPTION:

                                          Factorizes the positive definite, hermittian complex matrix 'x_mat' into
                                          L * L(transposed), when L is a lower triangular matrix. L is stored in the
                                          lower part of 'x_mat' on output.

                                          HOW TO USE THE FUNCTION:
                                          lib_matrCholCpx(i_dim, x_mat);

                                          SIDE-EFFECTS: The input matrix 'x_mat' is destroyed on output.

                                          RETURN VALUE: 0 if everything O.K.
                                          1 if illegal matrix

                                          ************************************************************************/
{
  int l_i, l_j, l_k, r_value = 0;
  fftw_complex l_r;
  float help;
  float factor = x_mat[0][0].re;
  if(factor <= 0)
    return(1);
  for(l_i=0;l_i<i_dim;l_i++)
    for(l_j=0;l_j<i_dim;l_j++)
    {
      x_mat[l_i][l_j].re = x_mat[l_i][l_j].re/factor;
      x_mat[l_i][l_j].im = x_mat[l_i][l_j].im/factor;
    }

    for (l_i=0; l_i < i_dim; l_i++)
      if (x_mat[l_i][l_i].re <= TOL) {
        /*
          fprintf(stderr, "Matrix to be Cholesky factorized not valid! %s\n",
          "Diagonal-element near 0.0");
          lib_matrPrintCpx( x_mat, i_dim, i_dim);
        */
        r_value = 1;
        return(r_value);
      } else {
        for (l_j=0; l_j < l_i; l_j++) {
          l_r.re = 0.0;
          l_r.im = 0.0;
          for (l_k = 0; l_k < l_j; l_k++){
            l_r.re += (x_mat[l_i][l_k].re * x_mat[l_j][l_k].re)+(x_mat[l_i][l_k].im * x_mat[l_j][l_k].im);
            l_r.im += -(x_mat[l_i][l_k].re * x_mat[l_j][l_k].im)+(x_mat[l_i][l_k].im * x_mat[l_j][l_k].re);
          }
          help = x_mat[l_j][l_j].re*x_mat[l_j][l_j].re + x_mat[l_j][l_j].im*x_mat[l_j][l_j].im;
          x_mat[l_i][l_j].re = ((x_mat[l_i][l_j].re - l_r.re)*x_mat[l_j][l_j].re
            +(x_mat[l_i][l_j].im - l_r.im)*x_mat[l_j][l_j].im)/help;

          x_mat[l_i][l_j].im = ((x_mat[l_i][l_j].im - l_r.im)*x_mat[l_j][l_j].re
            -(x_mat[l_i][l_j].re - l_r.re)*x_mat[l_j][l_j].im)/help;
        }
        l_r.re = 0.0;
        l_r.im = 0.0;
        for (l_k=0; l_k < l_i; l_k++){
          l_r.re += (x_mat[l_i][l_k].re * x_mat[l_i][l_k].re) + (x_mat[l_i][l_k].im * x_mat[l_i][l_k].im);
          /* l_r.im -= 2*x_mat[l_i][l_k].re * x_mat[l_i][l_k].im;*/
        }
        l_r.re = x_mat[l_i][l_i].re - l_r.re;
        /* l_r.im = x_mat[l_i][l_i].im - l_r.im;*/
        if (l_r.re <= TOL) {
          /*
          fprintf(stderr, "Matrix to be Cholesky factorized not valid !\n");
          lib_matrPrintCpx( x_mat, i_dim, i_dim);
          */
          r_value = 1;
          return(r_value);
        }

        assert(l_r.im<0.00001);
        x_mat[l_i][l_i].re = (float) (sqrt(l_r.re));
        x_mat[l_i][l_i].im = 0.0;
      }

      for(l_i=0;l_i<i_dim;l_i++)
        for(l_j=0;l_j<=l_i;l_j++)
        {
          x_mat[l_i][l_j].re *= (float) (sqrt(factor));
          x_mat[l_i][l_j].im *= (float) (sqrt(factor));
        }

        return(r_value);
}



void lib_matrAxeqbR(
                    int i_dim, /*The dimension of the equation system to solve. */
                    double **i_mat, /*The LU-decomp. of the A-matrix in the equation-system.*/
                    double *x_vec) /*On input the vector b, on output the solution x.*/
                    /*FUNC*****************************************************************
                    DESCRIPTION:

                    Solves the set of i_dim linear equations A * x = b. Here the matrix A is input
                    as its LU-decomposition. The function 'lib_matrCholR' should be used on
                    A before this function is called.

                    HOW TO USE THE FUNCTION:
                    lib_matrAxeqbR(i_dim, i_mat, x_vec);

                    SIDE-EFFECTS: The input-vector b is destroyed on output.

                    RETURN VALUE: void.

                    ************************************************************************/
{
  int l_i, l_j;
  double l_x;

  for (l_i = 0; l_i < i_dim; l_i++) {
    l_x = x_vec[l_i];
    for (l_j = 0; l_j < l_i; l_j++)
      l_x -= x_vec[l_j] * i_mat[l_i][l_j];
    x_vec[l_i] = l_x / i_mat[l_i][l_i];
  }

  for (l_i = i_dim - 1; l_i >= 0; l_i--) {
    l_x = x_vec[l_i];
    for (l_j = i_dim - 1; l_j > l_i; l_j--)
      l_x = l_x - x_vec[l_j] * i_mat[l_j][l_i];
    x_vec[l_i] = l_x / i_mat[l_i][l_i];
  }
}

void lib_matrAxeqbCpx(
                      int i_dim, /*The dimension of the equation system to solve. */
                      fftw_complex **i_mat, /*The LU-decomp. of the A-matrix in the equation-system.*/
                      fftw_complex *x_vec) /*On input the vector b, on output the solution x.*/
                      /*FUNC*****************************************************************
                      DESCRIPTION:

                      Solves the set of i_dim linear equations A * x = b for complex
                      numbers. Here the matrix A is input
                      as its LU-decomposition. The function 'lib_matrCholCpx' should be used on
                      A before this function is called.

                      HOW TO USE THE FUNCTION:
                      lib_matrAxeqbCpx(i_dim, i_mat, x_vec);

                      SIDE-EFFECTS: The input-vector b is destroyed on output.

                      RETURN VALUE: void.

                      ************************************************************************/
{
  int l_i, l_j;
  fftw_complex l_x;
  float help;

  for (l_i = 0; l_i < i_dim; l_i++) {
    l_x.re = x_vec[l_i].re;
    l_x.im = x_vec[l_i].im;
    for (l_j = 0; l_j < l_i; l_j++) {
      l_x.re -= (x_vec[l_j].re * i_mat[l_i][l_j].re - x_vec[l_j].im * i_mat[l_i][l_j].im);
      l_x.im -= (x_vec[l_j].im * i_mat[l_i][l_j].re + x_vec[l_j].re * i_mat[l_i][l_j].im);
    }
    help = (i_mat[l_i][l_i].re*i_mat[l_i][l_i].re+i_mat[l_i][l_i].im*i_mat[l_i][l_i].im);
    x_vec[l_i].re = (l_x.re*i_mat[l_i][l_i].re+l_x.im*i_mat[l_i][l_i].im)/help;
    x_vec[l_i].im = (l_x.im*i_mat[l_i][l_i].re-l_x.re*i_mat[l_i][l_i].im)/help;
  }

  for (l_i = i_dim - 1; l_i >= 0; l_i--) {
    l_x.re = x_vec[l_i].re;
    l_x.im = x_vec[l_i].im;
    for (l_j = i_dim - 1; l_j > l_i; l_j--) {
      l_x.re = l_x.re - (x_vec[l_j].re * i_mat[l_j][l_i].re + x_vec[l_j].im * i_mat[l_j][l_i].im) ;
      l_x.im = l_x.im - (-x_vec[l_j].re * i_mat[l_j][l_i].im + x_vec[l_j].im * i_mat[l_j][l_i].re);
    }
    help = (i_mat[l_i][l_i].re*i_mat[l_i][l_i].re+i_mat[l_i][l_i].im*i_mat[l_i][l_i].im);
    x_vec[l_i].re = (l_x.re*i_mat[l_i][l_i].re-l_x.im*i_mat[l_i][l_i].im)/help;
    x_vec[l_i].im = (l_x.im*i_mat[l_i][l_i].re+l_x.re*i_mat[l_i][l_i].im)/help;
  }

}

/*FUNC*****************************************************************
DESCRIPTION:

Solves the n set of i_dim linear equations A * X = B for real
numbers. Here the matrix A is input
as its LU-decomposition. The function 'lib_matrCholR' should be used on
A before this function is called.

B is i_dim*n, where n is number of columns.


HOW TO USE THE FUNCTION:
lib_matrAXeqBMatCpx(i_dim, i_mat, x_mat, n);

SIDE-EFFECTS: The input-matrix B is destroyed on output.

RETURN VALUE: void.

************************************************************************/

void lib_matrAXeqBMatR(int i_dim, double **i_mat, double **x_mat, int n)
{

  int l_i, l_j;
  double l_x;
  double help;
  int i;
  for(i=0;i<n;i++)
  {
    for (l_i = 0; l_i < i_dim; l_i++) {
      l_x = x_mat[l_i][i];
      for (l_j = 0; l_j < l_i; l_j++) {
        l_x -= (x_mat[l_j][i] * i_mat[l_i][l_j]);
      }
      help = (i_mat[l_i][l_i]*i_mat[l_i][l_i]);
      x_mat[l_i][i] = (l_x*i_mat[l_i][l_i])/help;
    }

    for (l_i = i_dim - 1; l_i >= 0; l_i--) {
      l_x = x_mat[l_i][i];

      for (l_j = i_dim - 1; l_j > l_i; l_j--) {
        l_x = l_x - (x_mat[l_j][i] * i_mat[l_j][l_i] ) ;
      }
      help = (i_mat[l_i][l_i]*i_mat[l_i][l_i]);
      x_mat[l_i][i] = (l_x*i_mat[l_i][l_i])/help;
    }
  }
}

/*FUNC*****************************************************************
DESCRIPTION:

Solves the n set of i_dim linear equations A * X = B for complex
numbers. Here the matrix A is input
as its LU-decomposition. The function 'lib_matrCholCpx' should be used on
A before this function is called.

B is i_dim*n, where n is number of columns.


HOW TO USE THE FUNCTION:
lib_matrAXeqBMatCpx(i_dim, i_mat, x_mat, n);

SIDE-EFFECTS: The input-matrix B is destroyed on output.

RETURN VALUE: void.

************************************************************************/

void lib_matrAXeqBMatCpx(int i_dim, fftw_complex **i_mat, fftw_complex **x_mat, int n)
{

  int l_i, l_j;
  fftw_complex l_x;
  float help;
  int i;
  for(i=0;i<n;i++)
  {
    for (l_i = 0; l_i < i_dim; l_i++) {
      l_x.re = x_mat[l_i][i].re;
      l_x.im = x_mat[l_i][i].im;
      for (l_j = 0; l_j < l_i; l_j++) {
        l_x.re -= (x_mat[l_j][i].re * i_mat[l_i][l_j].re - x_mat[l_j][i].im * i_mat[l_i][l_j].im);
        l_x.im -= (x_mat[l_j][i].im * i_mat[l_i][l_j].re + x_mat[l_j][i].re * i_mat[l_i][l_j].im);
      }
      help = (i_mat[l_i][l_i].re*i_mat[l_i][l_i].re+i_mat[l_i][l_i].im*i_mat[l_i][l_i].im);
      x_mat[l_i][i].re = (l_x.re*i_mat[l_i][l_i].re+l_x.im*i_mat[l_i][l_i].im)/help;
      x_mat[l_i][i].im = (l_x.im*i_mat[l_i][l_i].re-l_x.re*i_mat[l_i][l_i].im)/help;
    }

    for (l_i = i_dim - 1; l_i >= 0; l_i--) {
      l_x.re = x_mat[l_i][i].re;
      l_x.im = x_mat[l_i][i].im;
      for (l_j = i_dim - 1; l_j > l_i; l_j--) {
        l_x.re = l_x.re - (x_mat[l_j][i].re * i_mat[l_j][l_i].re + x_mat[l_j][i].im * i_mat[l_j][l_i].im) ;
        l_x.im = l_x.im - (-x_mat[l_j][i].re * i_mat[l_j][l_i].im + x_mat[l_j][i].im * i_mat[l_j][l_i].re);
      }
      help = (i_mat[l_i][l_i].re*i_mat[l_i][l_i].re+i_mat[l_i][l_i].im*i_mat[l_i][l_i].im);
      x_mat[l_i][i].re = (l_x.re*i_mat[l_i][l_i].re-l_x.im*i_mat[l_i][l_i].im)/help;
      x_mat[l_i][i].im = (l_x.im*i_mat[l_i][l_i].re+l_x.re*i_mat[l_i][l_i].im)/help;
    }
  }


}


/*FUNC*****************************************************************
DESCRIPTION:

Solves the n set of i_dim linear equations L * X = B for complex
numbers. Here the matrix L is  the L from Cholesky decomposition of a matrix A.
The function 'lib_matrCholCpx' should be used on
A before this function is called, and goes as input here.

This function does forward substitution as inlib_matrAXeqBMatCpx, but no backward substitution

B is i_dim*n, where n is number of columns.


HOW TO USE THE FUNCTION:
lib_matrAXeqBMatCpx(i_dim, i_mat, x_mat, n);

SIDE-EFFECTS: The input-matrix B is destroyed on output.

RETURN VALUE: void.

************************************************************************/

void lib_matrLXeqBMatCpx(int i_dim, fftw_complex **i_mat, fftw_complex **x_mat, int n)
{

  int l_i, l_j;
  fftw_complex l_x;
  float help;
  int i;
  for(i=0;i<n;i++)
  {
    for (l_i = 0; l_i < i_dim; l_i++) {
      l_x.re = x_mat[l_i][i].re;
      l_x.im = x_mat[l_i][i].im;
      for (l_j = 0; l_j < l_i; l_j++) {
        l_x.re -= (x_mat[l_j][i].re * i_mat[l_i][l_j].re - x_mat[l_j][i].im * i_mat[l_i][l_j].im);
        l_x.im -= (x_mat[l_j][i].im * i_mat[l_i][l_j].re + x_mat[l_j][i].re * i_mat[l_i][l_j].im);
      }
      help = (i_mat[l_i][l_i].re*i_mat[l_i][l_i].re+i_mat[l_i][l_i].im*i_mat[l_i][l_i].im);
      x_mat[l_i][i].re = (l_x.re*i_mat[l_i][l_i].re+l_x.im*i_mat[l_i][l_i].im)/help;
      x_mat[l_i][i].im = (l_x.im*i_mat[l_i][l_i].re-l_x.re*i_mat[l_i][l_i].im)/help;
    }
  }
}

void lib_matrLXeqMatR(
                    int i_dim, /*The dimension of the equation system to solve. */
                    double **i_mat, /*The LU-decomp. of the A-matrix in the equation-system.*/
                    double **x_mat, int n) /*On input the matrix B, on output the solution X.*/
                    /*FUNC*****************************************************************
                    DESCRIPTION:

                    Solves the set of i_dim linear equations L * X = B. L is from LU decomposition of a matrix A.
                    Here the matrix A is input
                    as its LU-decomposition. The function 'lib_matrCholR' should be used on
                    A before this function is called.
                    B is i_dim*n, where n is number of columns.

                    HOW TO USE THE FUNCTION:
                    lib_matrLXeqMatR(i_dim, i_mat, x_mat, n);

                    SIDE-EFFECTS: The input-matrix B is destroyed on output.

                    RETURN VALUE: void.

                    ************************************************************************/
{
  int l_i, l_j;
  double l_x;
  int i;
  for(i=0;i<n;i++)
  {
    for (l_i = 0; l_i < i_dim; l_i++) {
      l_x = x_mat[l_i][i];
      for (l_j = 0; l_j < l_i; l_j++)
        l_x -= x_mat[l_j][i] * i_mat[l_i][l_j];
      x_mat[l_i][i] = l_x / i_mat[l_i][l_i];
    }
  }
}
/*
Calculate  L x from a lover triangular matrix, mat, and the vector vec
mat is typically the result of an  cholesky factoring (see lib_matrCholCpx above)
it write over vec
*/
void lib_matrProdCholVec(int n, fftw_complex ** mat, fftw_complex * vec)
{
  int i,j;
  fftw_complex * inseed= (fftw_complex *) fftw_malloc( sizeof( fftw_complex) * n );

  for(i=0; i < n; i++) { inseed[i].re = vec[i].re; inseed[i].im = vec[i].im;}
  for(i=0; i < n; i++) { vec[i].re=0.0; vec[i].im=0.0;}

  for(i=0; i < n; i++)
    for(j=0; j < i+1 ; j++)
    {
      vec[i].re += mat[i][j].re * inseed[j].re - mat[i][j].im * inseed[j].im;
      vec[i].im += mat[i][j].re * inseed[j].im + mat[i][j].im * inseed[j].re;
    }

    fftw_free(inseed);

}



/*
Calculate  product of a scalar and n complex vector
*/
void lib_matrProdScalVecCpx(fftw_complex s, fftw_complex *vec, int n )
{
  fftw_complex  inseed;

  int i;
  for(i=0;i<n;i++)
  {
    inseed=vec[i];
    vec[i].re = s.re*inseed.re - s.im*inseed.im;
    vec[i].im = s.re*inseed.im + s.im*inseed.re;
  }
}

/*
Calculate  product of a scalar and n complex vector
*/
void lib_matrProdScalVecRCpx(float s, fftw_complex *vec,int n)
{
  int i;
  for(i=0;i<n;i++)
  {
    vec[i].re *= s;
    vec[i].im *= s;
  }
}

/*
Calculate  product of a scalar and n1 x n2 complex matrix
*/
void lib_matrProdScalMatCpx(fftw_complex s, fftw_complex **mat, int n1,
                            int n2)
{
  fftw_complex  inseed;
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)

    {
      inseed= mat[i][j];
      mat[i][j].re = s.re*inseed.re - s.im*inseed.im;
      mat[i][j].im = s.re*inseed.im + s.im*inseed.re;
    }
}


/*
Calculate matrix product of a n1 x n2 and the adjoint of an n3 x n2 complex matrix
*/
void lib_matrProdAdjointCpx(fftw_complex **mat1, fftw_complex **mat2, int n1,
                            int n2, int n3, fftw_complex **outmat)
{
  int i, j, k;
  fftw_complex x;
  for(i=0;i<n1;i++)
    for(j=0;j<n3;j++)
    {
      x.re = 0.0;
      x.im = 0.0;
      for(k=0;k<n2;k++)
      {
        x.re += mat1[i][k].re*mat2[j][k].re + mat1[i][k].im*mat2[j][k].im;
        x.im += mat1[i][k].im*mat2[j][k].re - mat1[i][k].re*mat2[j][k].im;
      }
      outmat[i][j].re = x.re;
      outmat[i][j].im = x.im;
    }
}

/*
Calculate matrix product of a n1 x n2 and n2 x n3 complex matrix
*/
void lib_matrProdCpx(fftw_complex **mat1, fftw_complex **mat2, int n1,
                     int n2, int n3, fftw_complex **outmat)
{
  int i, j, k;
  fftw_complex x;
  for(i=0;i<n1;i++)
    for(j=0;j<n3;j++)
    {
      x.re = 0.0;
      x.im = 0.0;
      for(k=0;k<n2;k++)
      {
        x.re += mat1[i][k].re*mat2[k][j].re - mat1[i][k].im*mat2[k][j].im;
        x.im += mat1[i][k].im*mat2[k][j].re + mat1[i][k].re*mat2[k][j].im;
      }
      outmat[i][j].re = x.re;
      outmat[i][j].im = x.im;
    }
}




/*
Calculate matrix product of a n1 x n2 complex and n2 x n3 real matrix
*/

void lib_matrProdCpxR(fftw_complex **mat1, float **mat2, int n1,
                      int n2, int n3, fftw_complex **outmat)
{
  int i, j, k;
  fftw_complex x;
  for(i=0;i<n1;i++)
    for(j=0;j<n3;j++)
    {
      x.re = 0.0;
      x.im = 0.0;
      for(k=0;k<n2;k++)
      {
        x.re += mat1[i][k].re*mat2[k][j];
        x.im += mat1[i][k].im*mat2[k][j];
      }
      outmat[i][j].re = x.re;
      outmat[i][j].im = x.im;
    }


}
/*
Calculate matrix product of a n1 x n2 real and n2 x n3 complex matrix
*/
void lib_matrProdRCpx(float **mat1, fftw_complex **mat2, int n1,
                      int n2, int n3, fftw_complex **outmat)
{
  int i, j, k;
  fftw_complex x;
  for(i=0;i<n1;i++)
    for(j=0;j<n3;j++)
    {
      x.re = 0.0;
      x.im = 0.0;
      for(k=0;k<n2;k++)
      {
        x.re += mat1[i][k]*mat2[k][j].re ;
        x.im += mat1[i][k]*mat2[k][j].im ;
      }
      outmat[i][j].re = x.re;
      outmat[i][j].im = x.im;
    }
}



/*
Calculate matrix product of a n1 x n1 diagonal complex matrix and n1 x n2 real matrix
*/
void lib_matrProdDiagCpxR(fftw_complex *mat1, float **mat2, int n1,
                          int n2, fftw_complex **outmat)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
    {
      outmat[i][j].re = mat1[i].re * mat2[i][j];
      outmat[i][j].im = mat1[i].im * mat2[i][j];
    }
}



/*
calculate matrix product FAF^* for complex matrices F: n1 x n2, A: n2 x n2
*/

void lib_matrProd2Cpx(fftw_complex **mat1, fftw_complex **mat2, int n1,
                      int n2, fftw_complex **outmat)
{
  int i, j, k, l;
  fftw_complex x;
  fftw_complex help;
  for(i=0;i<n1;i++)
    for(j=0;j<n1;j++)
    {
      x.re = 0.0;
      x.im = 0.0;
      for(k=0;k<n2;k++)
      {
        for(l=0;l<n2;l++)
        {
          help.re = (mat1[i][k].re*mat2[k][l].re - mat1[i][k].im*mat2[k][l].im);
          help.im = mat1[i][k].im*mat2[k][l].re + mat1[i][k].re*mat2[k][l].im;
          x.re +=help.re*mat1[j][l].re + help.im*mat1[j][l].im;
          x.im +=help.im*mat1[j][l].re - help.re*mat1[j][l].im;
        }
      }
      outmat[i][j].re = x.re;
      outmat[i][j].im = x.im;
    }
}



/*
Calculate the product of a n1 x n2 complex matrix and a n2 x 1 complex
vector, return in the n1 x 1 vector outvec.
*/
void lib_matrProdMatVecCpx(fftw_complex **mat, fftw_complex *vec, int n1, int n2,
                           fftw_complex *outvec)
{
  int i, j;
  fftw_complex x;
  for(i=0;i<n1;i++)
  {
    x.re = 0.0;
    x.im = 0.0;
    for(j=0;j<n2;j++)
    {
      x.re += mat[i][j].re*vec[j].re - mat[i][j].im*vec[j].im;
      x.im += mat[i][j].im*vec[j].re + mat[i][j].re*vec[j].im;
    }
    outvec[i].re = x.re;
    outvec[i].im = x.im;
  }
}

/*
Calculate the product of a n1 x n2 real matrix and a n2 x 1 complex
vector, return in the n1 x 1 vector outvec.
*/
void lib_matrProdMatRVecCpx(double **mat, fftw_complex *vec, int n1, int n2,
                           fftw_complex *outvec)
{
  int i, j;
  fftw_complex x;
  for(i=0;i<n1;i++)
  {
    x.re = 0.0;
    x.im = 0.0;
    for(j=0;j<n2;j++)
    {
      x.re += (float) (mat[i][j])*vec[j].re;
      x.im += (float) (mat[i][j])*vec[j].im;
    }
    outvec[i].re = x.re;
    outvec[i].im = x.im;
  }
}

/*
Calculate the product of the adjoint of a n2 x n1 complex matrix and a n2 x 1 complex
vector, return in the n1 x 1 vector outvec.
*/
void lib_matrProdAdjointMatVecCpx(fftw_complex **mat, fftw_complex *vec, int n1, int n2,
                                  fftw_complex *outvec)
{
  int i, j;
  fftw_complex x;
  for(i=0;i<n1;i++)
  {
    x.re = 0.0;
    x.im = 0.0;
    for(j=0;j<n2;j++)
    {
      x.re += mat[j][i].re*vec[j].re + mat[j][i].im*vec[j].im;
      x.im += -mat[j][i].im*vec[j].re + mat[j][i].re*vec[j].im;
    }
    outvec[i].re = x.re;
    outvec[i].im = x.im;
  }
}



/*
Addition of two complex matrices x and y. Result returned in y.
*/

void lib_matrAddMatCpx(fftw_complex **x, int n1, int n2, fftw_complex **y)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
    {
      y[i][j].re +=x[i][j].re;
      y[i][j].im +=x[i][j].im;
    }
}

void lib_matrSubtMatCpx(fftw_complex **x, int n1, int n2, fftw_complex **y)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
    {
      y[i][j].re -=x[i][j].re;
      y[i][j].im -=x[i][j].im;
    }
}

void lib_matrAddVecCpx(fftw_complex *x, int n,fftw_complex *y)
{
  int i;
  for(i=0;i<n;i++)
  {
    y[i].re +=x[i].re;
    y[i].im +=x[i].im;
  }
}

void lib_matrSubtVecCpx(fftw_complex *x, int n, fftw_complex *y)
{
  int i;
  for(i=0;i<n;i++)
  {
    y[i].re -=x[i].re;
    y[i].im -=x[i].im;
  }
}

/*
Find the adjoint(=conjungate transpose) to mat and return in output
*/

void lib_matrAdjoint(fftw_complex **mat, int n1, int n2, fftw_complex **outmat)
{
  int i, j;
  for(i = 0; i  < n1; i++)
    for(j = 0;j < n2; j++)
    {
      outmat[j][i].re = mat[i][j].re;
      outmat[j][i].im = -mat[i][j].im;
    }
}


/*
  Conjungates the matrix mat and returns it in
*/
void lib_matrConj(fftw_complex **mat, int n1, int n2, fftw_complex **outmat)
{
  int i, j;
  for(i = 0; i  < n1; i++)
    for(j = 0;j < n2; j++)
    {
      outmat[i][j].re =  mat[i][j].re;
      outmat[i][j].im = -mat[i][j].im;
    }

}

void lib_matrCopyCpx(fftw_complex **mat, int n1, int n2, fftw_complex **outmat)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
    {
      outmat[i][j].re = mat[i][j].re;
      outmat[i][j].im = mat[i][j].im;
    }

}

void lib_matrFillOnesVecCpx(fftw_complex* v1, int n1)
{
  int i;
  for(i=0;i<n1;i++)
  {
    v1[i].re = 1.0f;
    v1[i].im = 0.0f;
  }
}

void lib_matrFillValueVecCpx(fftw_complex value,fftw_complex* v1, int n1)
{
  int i;
  for(i=0;i<n1;i++)
  {
    v1[i].re = value.re;
    v1[i].im = value.im;
  }
}


/*
Addition of two matrices x and y. Result returned in y.
*/

void lib_matrAddMat(double **x, int n1, int n2, double **y)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      y[i][j] +=x[i][j];
}

void lib_matrSubtMat(double **x, int n1, int n2, double **y)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      y[i][j] -=x[i][j];
}

void lib_matrCopy(double **mat, int n1, int n2, double **outmat)
{
  int i, j;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      outmat[i][j] = mat[i][j];

}


void lib_matrPrintR(float ** mat, int n1, int n2)
{
  int i, j;
  for(i=0;i<n1;i++)
  {
    printf("\n ");
    for(j=0;j<n2;j++)
    {
      printf("%3e ",mat[i][j]);
    }
  }
  printf("\n ");
}


void lib_matrPrintCpx(fftw_complex ** mat, int n1, int n2)
{
  int i, j;
  printf("\nREAL:");
  for(i=0;i<n1;i++)
  {
    printf("\n ");
    for(j=0;j<n2;j++)
    {
      printf("%3e ",mat[i][j].re);
    }
  }

  printf("\nIMAG:");
  for(i=0;i<n1;i++)
  {
    printf("\n ");
    for(j=0;j<n2;j++)
    {
      printf("%3e ",mat[i][j].im);
    }
  }
  printf("\n ");
}

void  lib_matr_eigen(
   double      **i_mat,
   int         i_n,
   double      **o_eigvec,
   double       *o_eigval,
   int         *o_error) /* 0 if ok, 1 if not all eigenvalues determined.
                            input to the function tqli. */
/*FUNC******************************************************************

DESCRIPTION:

Calculate the eigenvectors and eigenvalues of a real, symmetric
matrix using the Householders algorithm described in 'Numerical recipes in C',
p 353-381. To do this we use the two local functions 'tred2' and 'tqli'.

HOW TO USE THE FUNCTION:

lib_matr_eigen(i_mat, i_n, o_eigvec, o_eigval, o_error);

SIDE-EFFECTS:

RETURN VALUE:

*********************************************************************/
{
 /* local variables */

  double *l_d, *l_e, **l_mat;
  int i,j;

/* Since the input to the functions tred2 and tqli will be destroyed, we
   make a copy of i_mat.  */

/*
  l_mat = (double **) Mmatrix_2d(0, i_n-1, 0, i_n-1, sizeof(double), 1);
  l_mat = new double *[i_n];
*/
  l_mat = malloc(sizeof(double)*i_n);
  for(i=0;i<i_n;i++)
    l_mat[i]= malloc(sizeof(double)*i_n);
  /*
    l_mat[i] = new double[i_n];
  */
  l_d = malloc(sizeof(double)*i_n);
  l_e = malloc(sizeof(double)*i_n);
  /*
    l_d = new double[i_n];
    l_e = new double[i_n];
    l_d = (double *) Mmatrix_1d(0,i_n-1, sizeof(double), 1);
    l_e = (double *) Mmatrix_1d(0,i_n-1, sizeof(double), 1);
  */
  for (i=0; i<=i_n-1; i++)
    for (j=0; j<=i_n-1; j++)
      l_mat[i][j] = i_mat[i][j];

  lib_matr_tred2(i_n, l_mat, l_d, l_e);
  lib_matr_tqli(l_e, i_n, l_mat, l_d, o_error);

  if (*o_error > 0) {
    for(i=0;i<i_n;i++)
      free(l_mat[i]);
    free(l_mat);
    free(l_d);
    free(l_e);
    return;
  }

  for (i=0; i<=i_n-1; i++) {
    o_eigval[i] = l_d[i];
    for (j=0; j<=i_n-1; j++)
      o_eigvec[i][j] = l_mat[i][j];
  }

  for(i=0;i<i_n;i++)
    free(l_mat[i]);
  /* delete [] l_mat; */
  free(l_mat);
  free(l_d);
  free(l_e);
  /*
     delete [] l_d;
     delete [] l_e;
     l_mat = (double **)Fmatrix_2d((char **)&l_mat[0][0], (char *)&l_mat[0]);
     l_d = (double *)Fmatrix_1d((char *)&l_d[0]);
     l_e = (double *)Fmatrix_1d((char *)&l_e[0]);
  */
}

static void lib_matr_tred2(int i_n, double **x_a, double o_d[], double o_e[])
/*FUNC**********************************************************************

DESCRIPTION:

Householder reduction of a real, symmetric, matrix a[0,...,n-1][0,...,n-1].
On output, x_a is replaced by the orthogonal matrix Q effecting the
transformation. d[0,...,n-1] returns the diagonal elements of the
tridiagonal matrix, and e[0,...,n-1] the off diagonal elements, e[0]=0.
Several statements, as noted in comments, can be omitted if only eigenvalues
are to be found, in which case x_a contains no useful information on output.
Otherwise they are to be included.
This implementation is the same as stated in 'Numerical recipes in C' p. 373,
except that here the arrays and matrices starts in 0.

HOW TO USE THE FUNCTION:

tred2(i_n, x_a, o_d, o_e);

SIDE-EFFECTS:

RETURN VALUE:

****************************************************************************/
{
  int l, k, j, i;
  double scale, hh, h, g, f;

  for (i = i_n - 1; i>=1; i--) {
    l = i-1;
    h = scale = 0.0;
    if (l > 0) {
      for (k=0; k<=l; k++)
        scale += fabs(x_a[i][k]);
      if (scale == 0)
        o_e[i] = x_a[i][l];
      else {
        for (k=0; k<= l; k++) {
          x_a[i][k] /= scale;
          h += x_a[i][k]*x_a[i][k];
        }
        f = x_a[i][l];
        g = f>0 ? -sqrt(h) : sqrt(h);
        o_e[i]=scale*g;
        h -= f*g;
        x_a[i][l] = f-g;
        f = 0.0;
        for (j=0; j<=l; j++) {
        /* Next statement can be omitted if eigenvectors not wanted */
          x_a[j][i] = x_a[i][j]/h;
          g = 0.0;
          for (k=0; k<=j; k++)
            g += x_a[j][k]*x_a[i][k];
          for (k=j+1; k<=l; k++)
            g += x_a[k][j]*x_a[i][k];
          o_e[j] = g/h;
          f += o_e[j]*x_a[i][j];
        }
        hh = f/(h+h);
        for (j=0; j<=l; j++) {
          f = x_a[i][j];
          o_e[j] = g = o_e[j]-hh*f;
          for (k=0; k<=j; k++)
            x_a[j][k] -= (f*o_e[k] + g*x_a[i][k]);
        }
      }
    }
    else
      o_e[i] = x_a[i][l];
    o_d[i] = h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
  o_d[0] = 0.0;
  o_e[0] = 0.0;
  /* Contents of this loop can be omitted if eigenvectors not wanted except
     for statement o_d[i] = x_a[i][i];  */
  for (i=0; i <= i_n-1; i++) {
    l = i-1;
    if (o_d[i]) {
      for (j=0; j <= l; j++) {
        g = 0.0;
        for (k=0; k<=l; k++)
          g += x_a[i][k]*x_a[k][j];
        for (k=0; k<=l; k++)
          x_a[k][j] -= g*x_a[k][i];
      }
    }
    o_d[i] = x_a[i][i];
    x_a[i][i] = 1.0;
    for (j=0; j<=l; j++)
      x_a[j][i]=x_a[i][j]=0.0;
  }
}



static void lib_matr_tqli(double i_e[], int i_n, double **x_z, double x_d[],
                 int *o_error)

/*FUNC************************************************************************

DESCRIPTION:

QL-algorithm with implicit shifts, to determine the eigenvalues and
eigenvectors of a real, symmetric, tridiagonal matrix, or of a real, symmetric
matrix previously reduced by 'tred2'. On input x_d[0,...i_n-1] contains the
diagonal elements of the tridiagonal matrix. On output, it returns the
eigenvalues. The vector i_e[0,...,i_n-1] inputs the subdiagonal elements of the
tridiagonal matrix, with i_e[0] arbitrary. On output i_e is destroyed. When
finding only the eigenvalues, several lines may be omitted, as noted in the
comments. If the eigenvectors of a tridiagonal matrix are desired, the matrix
x_z[0,...,i_n-1][0,...,i_n-1] is input as the identity matrix.
If the eigenvectors of a matrix that has been reduced by 'tred2' are
required, then x_z is input as the matrix output by 'tred2'.
In either case, the kth column of x_z
returns the normalized eigenvector corresponding to d[k].

HOW TO USE THE FUNCTION:

tqli(i_n, i_e, x_d, x_z);

SIDE-EFFECTS:

RETURN VALUE:

***************************************************************************/
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  *o_error = 0;
  for (i=1; i<=i_n-1; i++) i_e[i-1] = i_e[i];
  i_e[i_n-1] = 0.0;
  for (l=0; l<=i_n-1; l++) {
    iter=0;
    do {
      for (m=l; m<=i_n-2; m++) {
        dd = fabs(x_d[m]) + fabs(x_d[m+1]);
        if (fabs(i_e[m])+dd == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {*o_error = 1; return;}
        g = (x_d[l+1] - x_d[l]) / (2.0*i_e[l]);
        r = sqrt((g*g) + 1.0);
        g = x_d[m] - x_d[l] + i_e[l] /(g + SIGN(r,g));
        s = c = 1.0;
        p = 0.0;
        for (i=m-1; i>=l; i--) {
          f = s*i_e[i];
          b = c*i_e[i];
          if (fabs(f) >= fabs(g)) {
            c = g/f;
            r = sqrt((c*c) + 1.0);
            i_e[i+1] = f*r;
            c *= (s=1.0/r);
          } else {
              s = f/g;
              r = sqrt((s*s) + 1.0);
              i_e[i+1] = g*r;
              s *= (c=1.0/r);
            }
          g = x_d[i+1] - p;
          r = (x_d[i]-g) * s + 2.0 * c * b;
          p = s*r;
          x_d[i+1] = g + p;
          g = c*r-b;
          /* Next loop can be omitted if eigenvectors not wanted */
          for (k=0; k<=i_n-1; k++) {
            f = x_z[k][i+1];
            x_z[k][i+1] = s*x_z[k][i] + c*f;
            x_z[k][i] = c*x_z[k][i]-s*f;
          }
        }
        x_d[l] = x_d[l]-p;
        i_e[l] = g;
        i_e[m] = 0.0;
      }
    } while (m != l);
  }
}

void lib_matr_prod(
   double      **i_mat1,
   double      **i_mat2,
   int         i_n1,
   int         i_n2,
   int         i_n3,
   double      **o_mat)
/*FUNC******************************************************************

DESCRIPTION:

Calculate the matrix product of a n1 x n2 matrix and a n2 x n3
matrix

HOW TO USE THE FUNCTION:

SIDE-EFFECTS:

RETURN VALUE:

***********************************************************************/

{
 /* local variables */

  int         l_i,l_j,l_k;
  double      l_x;

  for(l_i = 0 ; l_i <= i_n1-1 ; l_i++)
    for(l_j = 0 ; l_j <= i_n3-1 ; l_j++)
    {
      l_x = 0.0000000000000;
      for(l_k = 0 ; l_k <= i_n2-1 ; l_k++)
        l_x = l_x + i_mat1[l_i][l_k] * i_mat2[l_k][l_j];
      o_mat[l_i][l_j] = l_x;
    }
}


void lib_matr_prod_sym(
   double      **i_mat1,
   double      **i_mat2,
   int         i_n1,
   int         i_n2,
   int         i_n3,
   double      **o_mat)
/*FUNC******************************************************************

DESCRIPTION:

Calculate the matrix product of a n1 x n2 matrix and a n2 x n3
matrix

HOW TO USE THE FUNCTION:

SIDE-EFFECTS:

RETURN VALUE:

***********************************************************************/

{
 /* local variables */

  int         l_i,l_j,l_k;
  double      l_x;

  printf("Using matrix multiplication for symmetric matrices\n");

  for(l_i = 0 ; l_i <= i_n1-1 ; l_i++)
    for(l_j = 0 ; l_j <= i_n3-1 ; l_j++)
    {
      l_x = 0.0000000000000;
      for(l_k = 0 ; l_k <= i_n2-1 ; l_k++)
        l_x = l_x + i_mat1[l_i][l_k] * i_mat2[l_j][l_k];  /* Not swap of i_mat2 indices*/
      o_mat[l_i][l_j] = l_x;
    }
}

/*
Find the transpose to mat and return in output
*/

void lib_matrTranspose(double **mat, int n1, int n2, double **outmat)
{
  int i, j;
  for(i = 0; i  < n1; i++)
    for(j = 0;j < n2; j++)
    {
      outmat[j][i] = mat[i][j];
    }
}

void lib_matrLtXeqBR(
                    int i_dim, /*The dimension of the equation system to solve. */
                    double **i_mat, /*The LU-decomp. of the A-matrix in the equation-system.*/
                    double **x_mat, int n) /*On input the vector b, on output the solution x.*/
                    /*FUNC*****************************************************************
                    DESCRIPTION:

                    Solves the set of i_dim linear equations L^T * X = B. Here the matrix A is input
                    as its LU-decomposition. The function 'lib_matrCholR' should be used on
                    A before this function is called. Only the backsubstitution is performed.

                    HOW TO USE THE FUNCTION:
                    lib_matrLtXeqBR(i_dim, i_mat, x_mat, n);

                    SIDE-EFFECTS: The input-matrix B is destroyed on output.

                    RETURN VALUE: void.

                    ************************************************************************/
{
  int l_i, l_j;
  double l_x;

  int i;
  for(i=0;i<n;i++)
  {
 /* for (l_i = 0; l_i < i_dim; l_i++) {
    l_x = x_mat[l_i][i];
    for (l_j = 0; l_j < l_i; l_j++)
      l_x -= x_mat[l_j][i] * i_mat[l_i][l_j];
    x_mat[l_i][i] = l_x / i_mat[l_i][l_i];
  }
*/
  for (l_i = i_dim - 1; l_i >= 0; l_i--) {
    l_x = x_mat[l_i][i];
    for (l_j = i_dim - 1; l_j > l_i; l_j--)
      l_x = l_x - x_mat[l_j][i] * i_mat[l_j][l_i];
    x_mat[l_i][i] = l_x / i_mat[l_i][l_i];
  }
  }
}

void lib_matr_sort3x3(double *eigenval, double **eigenvec)
{
  double *help1;
  int i,j,k;
  int *index;
  double **help2;
  double max;
  help1 = malloc(sizeof(double)*3);
  help2 = malloc(sizeof(double)*3);
  for(i=0;i<3;i++)
    help2[i]= malloc(sizeof(double)*3);
  index = malloc(sizeof(int)*3);

  for(j=0;j<3;j++)
  {
    max = -100000.0;
    for(i=0;i<3;i++)
    {
      if(j==0 ||(j==1 && i!=index[0]) ||(j==2 && i!=index[0] && i!=index[1]))
      {
        if(eigenval[i]>max)
        {
          max = eigenval[i];
          index[j] = i;
        }
      }
    }
  help1[j] = max;
  for(k=0;k<3;k++)
    help2[k][j] = eigenvec[k][index[j]];
  }

  for(i=0;i<3;i++)
  {
    eigenval[i] = help1[i];
    for(j=0;j<3;j++)
      eigenvec[i][j] = help2[i][j];
  }

  for(i=0;i<3;i++)
    free (help2[i]);
  free (help2);
  free (help1);
  free (index);
}

void lib_matrDump(const char * fName, double ** mat, int n1, int n2)
{
  int i,j;
  FILE * dump = fopen(fName,"w");
  if(dump != NULL) {
    for(i=0;i<n1;i++) {
      for(j=0;j<n2;j++)
        fprintf(dump, "%f ", mat[i][j]);
      fprintf(dump,"\n");
    }
    fclose(dump);
  }
}

 extern void lib_matrDumpCpx(const char * fName, fftw_complex **  mat, int n1, int n2)
 {
   int i,j;
   FILE * dump = fopen(fName,"w");
   if(dump != NULL) {
    for(i=0;i<n1;i++) {
      for(j=0;j<n2;j++)
        fprintf(dump, "%f ", mat[i][j].re);
      fprintf(dump,"\n");
    }
    for(i=0;i<n1;i++) {
      for(j=0;j<n2;j++)
        fprintf(dump, "%f ", mat[i][j].im);
      fprintf(dump,"\n");
    }

    fclose(dump);
  }

 }

