#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "haar.h"
#define SQRT1_2 sqrt(1.0/2)

// -----------------------------------------------------------------------------
// diego domenzain
// spring 2021
// 
// code from the book,
// 
// @book{rostamian2014programming,
//   title={Programming Projects in C for Students of Engineering, Science, and Mathematics},
//   author={Rostamian, Rouben},
//   year={2014},
//   publisher={SIAM}
// }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//    1d
// -----------------------------------------------------------------------------

// forward haar transform (vector)
static void haar1d_fwd(double *v, int n)
{
  double h = sqrt(n);
  for (int i = 0; i < n; i++)
    v[i] /= h;
  for (int d = 1; d < n; d *= 2)
    for (int i=0; i < n; i += 2*d){
      double x = SQRT1_2 * (v[i] + v[i+d]); 
      double y = SQRT1_2 * (v[i] - v[i+d]); 
      v[i]  = x;
      v[i+d]= y;
    }
}

// inverse haar transform (vector)
static void haar1d_inv(double *v, int n)
{
  double h = sqrt(n);
  for (int i = 0; i < n; i++)
    v[i] *= h;
  for (int d = n/2; d >= 1; d /= 2){  
    for (int i=0; i < n; i += 2*d){
      double x = SQRT1_2 * (v[i] + v[i+d]); 
      double y = SQRT1_2 * (v[i] - v[i+d]);
      v[i]  = x;
      v[i+d]= y;
    }
  }
}

// -----------------------------------------------------------------------------
//    2d
// -----------------------------------------------------------------------------

// forward haar transform (matrix)
static void haar2d_fwd(double **a, int m, int n)
{
  // transform rows
  for (int i = 0; i < m; i++)
    haar1d(a[i], n, HA_FWD);
  
  // transform columns
  double h = sqrt(m);
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++)
      a[i][j] /= h;
    for (int d = 1; d < m; d*=2)
      for (int i = 0; i < m; i += 2*d) {
        double x = SQRT1_2 * (a[i][j] + a[i+d][j]);
        double y = SQRT1_2 * (a[i][j] - a[i+d][j]);
        a[i][j] = x;
        a[i+d][j] = y;
      }
  }
}

// inverse haar transform (matrix)
static void haar2d_inv(double **a, int m, int n)
{
  // transform rows
  for (int i = 0; i < m; i++)
    haar1d(a[i], n, HA_INV);
  
  // transform columns
  double h = sqrt(m);
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++)
      a[i][j] *= h;
    for (int d = m/2; d >= 1; d /= 2)
      for (int i = 0; i < m; i += 2*d) {
        double x = SQRT1_2 * (a[i][j] + a[i+d][j]);
        double y = SQRT1_2 * (a[i][j] - a[i+d][j]);
        a[i][j] = x;
        a[i+d][j] = y;
      }
  }
}

// -----------------------------------------------------------------------------
//    1d and 2d public
// -----------------------------------------------------------------------------

// forward or inverse Haar transform (vector)
void haar1d(double *v, int n, int dir)
{
  if (dir == HA_FWD) 
    haar1d_fwd(v, n);
  else if (dir == HA_INV)
    haar1d_inv(v, n);
  else { // shouldn’t be here!
    fprintf(stderr, "*** error in haar1d(): "
      "the third argument should be one of "
      "HA_FWD or HA_INV\n");
  exit(EXIT_FAILURE);
  }
}

// forward or inverse Haar transform (matrix)
void haar2d(double **a, int m, int n, int dir)
{
  if (dir == HA_FWD) 
    haar2d_fwd(a, m, n);
  else if (dir == HA_INV)
    haar2d_inv(a, m, n);
  else { // shouldn’t be here!
    fprintf(stderr, "*** error in haar1d(): "
      "the third argument should be one of "
      "HA_FWD or HA_INV\n");
  exit(EXIT_FAILURE);
  }
}
