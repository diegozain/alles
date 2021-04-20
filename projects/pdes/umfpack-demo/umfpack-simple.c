#include <stdlib.h>
#include <stdio.h>
#include "../../../src/c/sparse.h"
#include "../../../src/c/xmalloc.h"
#include "../../../../SuiteSparse/UMFPACK/Include/umfpack.h"
// -----------------------------------------------------------------------------
// solve Ax = b using LU decomposition with the the umfpack lib
// -----------------------------------------------------------------------------
int main(void)
{
    void *Symbolic, *Numeric;
    double **A;
    double *b, *x, *Ax;
    int *Ap, *Ai;
    int n = 5;
    int i, j, nz;
    
    // matrix A
    //
    make_matrix(A, n, n);
    
    // populate zeros
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            A[i][j] = 0.0;
        
    // populate non-zeros of A
    // A:
    // 2   3   0   0   0
    // 3   0   4   0   6
    // 0  -1  -3   2   0
    // 0   0   1   0   0
    // 0   4   2   0   1
    A[0][0] = 2; A[0][1] = 3;
    A[1][0] = 3; A[1][2] = 4; A[1][4] = 6;
    A[2][1] = -1; A[2][2] = -3; A[2][3] = 2;
    A[3][2] = 1;
    A[4][1] = 4; A[4][2] = 2; A[4][4] = 1;
    
    // number of non-zero entries of A
    nz = nonz(A,n,n);
    
    // sparsify A
    make_vector( Ap, n+1 );
    Ap[0] = 0;
    make_vector( Ai, nz );
    make_vector( Ax, nz );
    
    sparse_pack(A, n, n, Ap, Ai, Ax);
    
    // vector b
    make_vector(b, n);
    
    // populate b
    // b:
    // 8  45  -3   3  19
    b[0] = 8; b[1] = 45; b[2] = -3; b[3] = 3; b[4] = 19;
    
    // ----------------   UMFPACK magic ----------------------------------------
    make_vector(x,n);
    // get permutation matrices
    umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    // LU decomposition
    umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    // get x with simple matrix-multiply from LU
    umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
    // ----------------   end of UMFPACK magic ---------------------------------
    
    // print
    printf("\nA:\n");
    print_matrix("%7.2f", A, n, n);
    printf("b:\n");
    print_vector("%7.2f", b, n);
    printf("\n ----- solve Ax = b -----\n\n");
    printf("x \n");
    print_vector("%7.2f \n", x, n);
    printf("\n ----- sparse A -----\n\n");
    printf("Ap \n");
    print_vector("%7i", Ap, n+1);
    printf("Ai \n");
    print_vector("%7i", Ai, nz);
    printf("Ax \n");
    print_vector("%7.2f", Ax, nz);
    
    // clean up
    free_matrix(A);
    free_vector(b);
    free_vector(Ap);
    free_vector(Ai);
    free_vector(Ax);
    free_vector(x);
    
    // clean up
    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);
    
    return 0;
}
