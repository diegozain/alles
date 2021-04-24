#include <stdlib.h>
#include <stdio.h>
#include "../../../src/c/xmalloc.h"
#include "../../../src/c/array.h"
#include "../../../../SuiteSparse/UMFPACK/Include/umfpack.h"
// -----------------------------------------------------------------------------
// 
// solve Ax = b using LU decomposition with the the umfpack lib
// 
// -----------------------------------------------------------------------------
// some error nonsense
static void error_and_exit(int status, const char *file, int line)
{
    fprintf(stderr, "*** file %s, line %d: ", file, line); 
    switch (status) {
        case UMFPACK_ERROR_out_of_memory: 
            fprintf(stderr, "out of memory!\n"); 
            break;
        case UMFPACK_WARNING_singular_matrix: 
            fprintf(stderr, "matrix is singular!\n"); 
            break;
        default: 
            fprintf(stderr, "UMFPACK error code %d\n", status);
    }
    exit(EXIT_FAILURE);
}
// -----------------------------------------------------------------------------
int main(void)
{
    void *Symbolic, *Numeric;
    double *b, *x;
    int n = 5;
    
    // nz_=20 is an overestimate, nz_=12 would do fine.
    // however, vector I and J *can* be larger than nz.
    // if so, the idea is to sum values in repeated entries.
    // clever way for defining PDE matrices!
    int nz_ = 20; 
    int *Ti, *Tj;
    double *Tx;
    
    int *Ap, *Ai;
    double *Ax;
    
    int nz;
    
    int status;
    // ----------------   matrix A ---------------------------------------------
    make_vector(Ti, nz_);
    make_vector(Tj, nz_);
    make_vector(Tx, nz_);
    Ti[0] = 0; Tj[0] = 0; Tx[0] = 2;
    Ti[1] = 1; Tj[1] = 0; Tx[1] = 3;
    Ti[2] = 0; Tj[2] = 1; Tx[2] = 3;
    Ti[3] = 2; Tj[3] = 1; Tx[3] = -1;
    Ti[4] = 4; Tj[4] = 1; Tx[4] = 4;
    Ti[5] = 1; Tj[5] = 2; Tx[5] = 4;
    Ti[6] = 2; Tj[6] = 2; Tx[6] = -3;
    Ti[7] = 3; Tj[7] = 2; Tx[7] = 1;
    Ti[8] = 4; Tj[8] = 2; Tx[8] = 2;
    Ti[9] = 2; Tj[9] = 3; Tx[9] = 2;
    Ti[10]= 1; Tj[10]= 4; Tx[10]= 6;
    Ti[11]= 4; Tj[11]= 4; Tx[11]= 1;
    // ----------------   get CSS  ---------------------------------------------
    make_vector(Ai, nz_);
    make_vector(Ap, n+1);
    make_vector(Ax, nz_);
    umfpack_di_triplet_to_col(n, n, nz_, Ti, Tj, Tx, Ap, Ai, Ax, NULL);
    
    // get true number of nonzero elements
    nz = Ap[n];
    
    // clean - we dont need these anymore
    free_vector(Ti);
    free_vector(Tj);
    free_vector(Tx);
    // ----------------   vector b ---------------------------------------------
    make_vector(b, n);
    
    // populate b
    // b:
    // 8  45  -3   3  19
    b[0] = 8; b[1] = 45; b[2] = -3; b[3] = 3; b[4] = 19;
    // ----------------   UMFPACK magic ----------------------------------------
    make_vector(x,n);
    // get permutation matrices
    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    if (status != UMFPACK_OK)
          error_and_exit(status, __FILE__, __LINE__);
    // LU decomposition
    umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    // get x with simple matrix-multiply from LU
    umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
    // ----------------   end of UMFPACK magic ---------------------------------
    // print
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
    // ----------------   clean up  --------------------------------------------
    // clean up
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
