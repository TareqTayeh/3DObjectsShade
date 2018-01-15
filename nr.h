/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #3
 * Date: 9 November 2017
 * Program purpose: The purpose of this program is to design and implementat a program that allows a user to shade
 *                  objects (spheres, cones, and torii) in various colours. The program should scale the parametric
 *                  objects within a scene (using 3D homogeneous transformation matrices) and render them using
 *                  a simple Lambertian shading model that assumes diffuse light reflection from objects.
 * File: nr.h
 */


#ifndef NR_H
#define NR_H

#include <iostream>

void nrerror(const char error_text[])

{
  printf("Run-time error...\n") ;
  printf("%s\n",error_text) ;
  exit(1) ;
}


float *vector(int nl, int nh)

{ float *v ;

  v = new float[nh - nl + 1];
  if (!v) {
    nrerror("VECTOR: allocation failure") ;
  }
  return (v - nl) ;
}


int *ivector(int nl, int nh)

{ int *v ;

  v = new int[nh - nl + 1];
  if (!v) {
    nrerror("IVECTOR: allocation failure") ;
  }
  return (v - nl) ;
}


double *dvector(int nl,int nh)

{ double *v ;

  v = new double[nh -nl + 1];
  if (!v) {
    nrerror("DVECTOR: allocation failure") ;
  }
  return (v - nl) ;
}


float **matrix(int nrl, int nrh, int ncl, int nch)

{ int i ;
  float **m ;

  m = new float*[nrh - nrl + 1];
  if (!m) {
    nrerror("MATRIX: allocation failure") ;
  }
  m -= nrl ;

  for (i = nrl ; i <= nrh ; i++) {
    m[i] = new float[nch - ncl + 1];
    if (!m[i]) {
      nrerror("MATRIX: allocation failure") ;
    }
    m[i] -= ncl ;
  }
  return (m) ;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)

{ int i ;
  double **m ;

  m = new double*[nrh - nrl + 1];
  if (!m) {
    nrerror("DMATRIX: allocation failure") ;
  }
  m -= nrl ;

  for (i = nrl ; i <= nrh ; i++) {
    m[i] = new double[nch - ncl + 1];
    if (!m[i]) {
      nrerror("DMATRIX: allocation failure") ;
    }
    m[i] -= ncl ;
  }
  return (m) ;
}


int **imatrix(int nrl, int nrh, int ncl, int nch)

{ int i ;
  int **m ;

  m = new int*[nrh - nrl + 1];
  if (!m) {
    nrerror("IMATRIX: allocation failure") ;
  }
  m -= nrl ;

  for (i = nrl ; i <= nrh ; i++) {
    m[i] = new int(nch - ncl + 1);
    if (!m[i]) {
      nrerror("IMATRIX: allocation failure") ;
    }
    m[i] -= ncl ;
  }
  return (m) ;
}


float **submatrix(float **a, int oldrl,int oldrh, int oldcl, int oldch, int newrl, int newcl)

{ int i, j ;

  float **m ;

  m = new float*[oldrh - oldrl + 1];
  if (!m) {
    nrerror("SUBMATRIX: allocation failure") ;
  }
  m -= newrl ;

  for (i = oldrl, j = newrl ; i <= oldrh ; i++, j++) {
    m[j] = a[i] + oldcl - newcl ;
  }
  return (m) ;
}


void free_vector(float *v, int nl, int nh)

{ delete (v+nl); }


void free_ivector(int *v, int nl, int nh)

{ delete (v+nl); }


void free_dvector(double *v, int nl, int nh)

{ delete (v+nl); }


void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

  for (i = nrh ; i >= nrl ; i--) {
    delete (m[i] + nrl);
  }
  delete(m + nrl) ;
}


void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

  for (i = nrh ; i >= nrl ; i--) {
    delete(m[i] + ncl);
  }
  delete(m + nrl);
}


void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

  for (i = nrh ; i >= nrl ; i--) {
    delete(m[i] + ncl) ;
  }
  delete(m + nrl) ;
}


void free_submatrix(float **b, int nrl, int nrh,int ncl, int nch)

{ delete(b + nrl); }


float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch)

{ int i, j, nrow, ncol ;
  float **m ;

  nrow = nrh - nrl + 1 ;
  ncol = nch - ncl + 1 ;

  m = new float*[(nrow)*sizeof(float *)];
  if (!m) {
    nrerror("CONVERT_MATRIX: allocation failure") ;
  }
  m -= nrl ;
  for (i = 0, j = nrl ; i <= nrow - 1 ; i++, j++) {
    m[j] = a + ncol*i - ncl ;
  }
  return(m) ;
}

void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch)

{ delete(b + nrl) ; }

#endif // NR_H
