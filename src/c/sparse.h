#ifndef H_SPARSE_H
#define H_SPARSE_H
#include "array.h"
// -----------------------------------------------------------------------------
// diego domenzain
// sometime after fall 2017
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
// sparse "a" to full (m x n)
void sparse_unpack(double **a, int m, int n, int *Ap, int *Ai, double *Ax);

// full "a" to sparse (m x n)
void sparse_pack(double **a, int m, int n, int *Ap, int *Ai, double *Ax);

// get non-zero elements of a (m x n)
int nonz(double **a, int m, int n);

// Aj to Ap
int *sparse_Aj2Ap(int *Aj, int n, int nnz);

#endif // H_SPARSE_H
