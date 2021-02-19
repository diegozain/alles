#ifndef H_WAVELET_H
#define H_WAVELET_H
#define HA_FWD 1
#define HA_INV -1
void haar1d(double *v, int n, int dir);
void haar2d(double **a, int m, int n, int dir);
#endif /* H_WAVELET_H */

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
