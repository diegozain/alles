#include <cuda_runtime.h>
#include <cusparse.h>
#include <cstdint>

void Lxysolver(
  int32_t nl,
  int32_t nnz,
  int32_t n,
  const double* Vs_,   // [nnz * n], device
  const int32_t* J_,   // [nnz], device
  const int32_t* I_,   // [nl + 1], device
  const double* ys_,   // [nl * n], device
  double* xs_          // [nl * n], device
)
{
  cusparseHandle_t handle;
  cusparseCreate(&handle);

  // ------------------------------------
  // Sparse matrix descriptor (structure)
  // ------------------------------------
  cusparseSpMatDescr_t matA;
  cusparseCreateCsr(
    &matA,
    nl, nl, nnz,
    (void*)I_, (void*)J_, nullptr,   // values set later
    CUSPARSE_INDEX_32I,
    CUSPARSE_INDEX_32I,
    CUSPARSE_INDEX_BASE_ZERO,
    CUDA_R_64F
  );

  // ------------------------------------
  // Dense vector descriptors
  // ------------------------------------
  cusparseDnVecDescr_t vecY, vecZ, vecX;

  double* z_;
  cudaMalloc(&z_, nl * sizeof(double));

  cusparseCreateDnVec(&vecY, nl, nullptr, CUDA_R_64F);
  cusparseCreateDnVec(&vecZ, nl, z_,    CUDA_R_64F);
  cusparseCreateDnVec(&vecX, nl, nullptr, CUDA_R_64F);

  // ------------------------------------
  // SpILU setup
  // ------------------------------------
  cusparseSpILUDescr_t iluDescr;
  cusparseSpILU_createDescr(&iluDescr);

  size_t bufferSizeILU = 0;
  cusparseSpILU_bufferSize(
    handle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    matA,
    iluDescr,
    CUDA_R_64F,
    &bufferSizeILU
  );

  void* bufferILU;
  cudaMalloc(&bufferILU, bufferSizeILU);

  // Symbolic ILU analysis (once)
  cusparseSpILU_analysis(
    handle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    matA,
    iluDescr,
    CUDA_R_64F,
    bufferILU
  );

  // ------------------------------------
  // SpSV (triangular solve) setup
  // ------------------------------------
  cusparseSpSVDescr_t spsvL, spsvU;
  cusparseSpSV_createDescr(&spsvL);
  cusparseSpSV_createDescr(&spsvU);

  size_t bufferSizeL = 0, bufferSizeU = 0;
  const double alpha = 1.0;

  // Lower triangular (L)
  cusparseSpSV_bufferSize(
    handle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha,
    matA,
    vecY,
    vecZ,
    CUDA_R_64F,
    CUSPARSE_SPSV_ALG_DEFAULT,
    spsvL,
    &bufferSizeL
  );

  // Upper triangular (U)
  cusparseSpSV_bufferSize(
    handle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha,
    matA,
    vecZ,
    vecX,
    CUDA_R_64F,
    CUSPARSE_SPSV_ALG_DEFAULT,
    spsvU,
    &bufferSizeU
  );

  size_t bufferSizeSolve = std::max(bufferSizeL, bufferSizeU);
  void* bufferSolve;
  cudaMalloc(&bufferSolve, bufferSizeSolve);

  // Analysis (once)
  cusparseSpSV_analysis(
    handle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha,
    matA,
    vecY,
    vecZ,
    CUDA_R_64F,
    CUSPARSE_SPSV_ALG_DEFAULT,
    spsvL,
    bufferSolve
  );

  cusparseSpSV_analysis(
    handle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    &alpha,
    matA,
    vecZ,
    vecX,
    CUDA_R_64F,
    CUSPARSE_SPSV_ALG_DEFAULT,
    spsvU,
    bufferSolve
  );

  // ------------------------------------
  // Solve loop (required)
  // ------------------------------------
  for (int32_t i = 0; i < n; ++i)
  {
    const double* Vsrc = Vs_ + i * nnz;
    const double* y  = ys_ + i * nl;
    double* x = xs_ + i * nl;

    // Update matrix values
    cusparseCsrSetValues(matA, (void*)Vsrc);

    cusparseDnVecSetValues(vecY, (void*)y);
    cusparseDnVecSetValues(vecX, (void*)x);

    // Numeric ILU
    cusparseSpILU_factorize(
      handle,
      CUSPARSE_OPERATION_NON_TRANSPOSE,
      matA,
      iluDescr,
      CUDA_R_64F,
      bufferILU
    );

    // Check zero pivot
    int position;
    cusparseSpILU_zeroPivot(handle, iluDescr, &position);

    // Forward solve: L z = y
    cusparseSpSV_solve(
      handle,
      CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha,
      matA,
      vecY,
      vecZ,
      CUDA_R_64F,
      CUSPARSE_SPSV_ALG_DEFAULT,
      spsvL
    );

    // Backward solve: U x = z
    cusparseSpSV_solve(
      handle,
      CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha,
      matA,
      vecZ,
      vecX,
      CUDA_R_64F,
      CUSPARSE_SPSV_ALG_DEFAULT,
      spsvU
    );
  }

  // ------------------------------------
  // Cleanup
  // ------------------------------------
  cudaFree(z_);
  cudaFree(bufferILU);
  cudaFree(bufferSolve);

  cusparseSpSV_destroyDescr(spsvL);
  cusparseSpSV_destroyDescr(spsvU);
  cusparseSpILU_destroyDescr(iluDescr);

  cusparseDestroySpMat(matA);
  cusparseDestroyDnVec(vecY);
  cusparseDestroyDnVec(vecZ);
  cusparseDestroyDnVec(vecX);

  cusparseDestroy(handle);
}


int32_t main() 
{
  // 🧠
  int32_t nl = 5;
  int32_t nnz = 13;
  int32_t n = 2; // Number of systems to solve

  double* Vs = (double*)malloc(nnz * n * sizeof(double));
  int32_t* J = (int32_t*)malloc(nnz * sizeof(int32_t));
  int32_t* I = (int32_t*)malloc((nl + 1) * sizeof(int32_t));
  double* ys = (double*)malloc(nl * n * sizeof(double));
  double* xs = (double*)malloc(nl * n * sizeof(double));

  // 🖊️
  // system 1:
  // [ 10  0  0  0  0 ]
  // [  3  9  0  0  0 ]
  // [  0  7  8  7  0 ]
  // [  3  0  8  7  5 ]
  // [  0  8  0  9  9 ]
  // 
  // system 2:
  // [  5  0  0  0  0 ]
  // [  2  6  0  0  0 ]
  // [  0  4  7  3  0 ]
  // [  2  0  6  4  3 ]
  // [  0  7  0  8  6 ]
  // 
  double Vs_example[] = {
    10, 3, 9, 7, 8, 7, 3, 8, 7, 5, 8, 9, 9, // system 1
    5, 2, 6, 4, 7, 3, 2, 6, 4, 3, 7, 8, 6  //  system 2
  };
  int32_t J_example[] = {0, 0, 1, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4};
  int32_t I_example[] = {0, 1, 3, 6, 10, 13};
  double ys_example[] = {
    1, 2, 3, 4, 5, // system 1
    5, 4, 3, 2, 1  // system 2
  };

  memcpy(Vs, Vs_example, nnz * n * sizeof(double));
  memcpy(J, J_example, nnz * sizeof(int32_t));
  memcpy(I, I_example, (nl + 1) * sizeof(int32_t));
  memcpy(ys, ys_example, nl * n * sizeof(double));

  // 🚀🐦
  Lxysolver(nl, nnz, n, Vs, J, I, ys, xs);

  // 🖨️
  std::cout << "solutions xs:" << std::endl;
  for (int32_t i = 0; i < n; i++) 
  {
    std::cout << "system " << i + 1 << ": ";
    for (int32_t j = 0; j < nl; j++) 
    {
      std::cout << xs[i * nl + j] << " ";
    }
    std::cout << std::endl;
  }

  // 🚿
  free(Vs);
  free(J);
  free(I);
  free(ys);
  free(xs);

  return 0;
}