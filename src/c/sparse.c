#include "sparse.h"
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
void sparse_unpack(double **a, int m, int n, int *Ap, int *Ai, double *Ax)
{
    // populate zeros
    int i,j,k;
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            a[i][j] = 0.0;
    
    // populate non-zeros
    for(j=0; j<n; j++)
        for(k=Ap[j]; k<Ap[j+1]; k++){
            i=Ai[k];
            a[i][j] = Ax[k];
        }
}

int nonz(double **a, int m, int n)
{
    // find # of non-zero entries in a
    int i,j;
    int nnz = 0;
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            if (a[i][j] != 0.0){
                nnz++;
            }
        }
    }
    return nnz;
}

// full "a" to sparse (m x n)
void sparse_pack(double **a, int m, int n, int *Ap, int *Ai, double *Ax)
{
    // populate Ai, Ap, Ax
    int nz=0;
    int ap = 0;
    for(int j=0; j<n; j++){
        for(int i=0; i<m; i++){
            if (a[i][j] != 0.0){
                Ai[nz] = i;
                Ax[nz] = a[i][j];
                nz++;
                ap++;
            }
            Ap[j+1] = ap;
        }
    }
}

// Aj to Ap
int *sparse_Aj2Ap(int *Aj, int n, int nnz)
{
//    int nnz = sizeof Aj / Aj[0];
   int *Ap;
   make_vector( Ap, n+1 );
   Ap[0] = 0;
   Ap[n] = nnz;
   int i = 1;
   
   for(int j=0; j<nnz-1; j++){
       if (Aj[j] != Aj[j+1]){
           Ap[i] = j+1;
           i++;
       }
   }
   return Ap;
   
   free_vector(Ap);
}














