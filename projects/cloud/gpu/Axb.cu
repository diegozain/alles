#include <iostream>
#include <cuda_runtime.h>
#include <cublas_v2.h>
// --------------------------------------------------------------------
#define CHECK_CUDA(call)                                            \
  {                                                                 \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
      std::cerr << "CUDA error in " << __FILE__ << " at line "      \
            << __LINE__ << ": " << cudaGetErrorString(err) << "\n"; \
      exit(EXIT_FAILURE);                                           \
    }                                                               \
  }

#define CHECK_CUBLAS(call)                                          \
  {                                                                 \
    cublasStatus_t err = call;                                      \
    if (err != CUBLAS_STATUS_SUCCESS) {                             \
      std::cerr << "cuBLAS error in " << __FILE__ << " at line "    \
            << __LINE__ << ": " << err << "\n";                     \
      exit(EXIT_FAILURE);                                           \
    }                                                               \
  }
// --------------------------------------------------------------------
void solveBatchedLinearSystems(int32_t n, int32_t m, double* As, double* bs, double* xs) {
  double *As_, *bs_;
  double **As_array_, **bs_array_;
  int32_t *pivots_, *info_;

  cublasHandle_t handle_;
  CHECK_CUBLAS(cublasCreate(&handle_));

  // Allocate GPU memory
  CHECK_CUDA(cudaMalloc((void**)&As_, m * n * n * sizeof(double)));
  CHECK_CUDA(cudaMalloc((void**)&bs_, m * n * sizeof(double)));
  CHECK_CUDA(cudaMalloc((void**)&pivots_, m * n * sizeof(int32_t)));
  CHECK_CUDA(cudaMalloc((void**)&info_, m * sizeof(int32_t)));

  // Copy data from pinned host memory to device memory
  CHECK_CUDA(cudaMemcpy(As_, As, m * n * n * sizeof(double), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(bs_, bs, m * n * sizeof(double), cudaMemcpyHostToDevice));

  // Allocate and set up batched matrix pointers
  CHECK_CUDA(cudaMalloc((void**)&As_array_, m * sizeof(double*)));
  CHECK_CUDA(cudaMalloc((void**)&bs_array_, m * sizeof(double*)));

  double **As_array, **bs_array;
  CHECK_CUDA(cudaMallocHost((void**)&As_array, m * sizeof(double*))); // Pinned memory for faster transfer
  CHECK_CUDA(cudaMallocHost((void**)&bs_array, m * sizeof(double*))); // Pinned memory for faster transfer

  for (int32_t i = 0; i < m; i++) {
    As_array[i] = As_ + i * n * n;
    bs_array[i] = bs_ + i * n;
  }

  CHECK_CUDA(cudaMemcpy(As_array_, As_array, m * sizeof(double*), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(bs_array_, bs_array, m * sizeof(double*), cudaMemcpyHostToDevice));

  // Perform LU factorization
  CHECK_CUBLAS(cublasDgetrfBatched(handle_, n, As_array_, n, pivots_, info_, m));

  // Solve Ax = b using LU factors
  CHECK_CUBLAS(cublasDgetrsBatched(handle_, CUBLAS_OP_N, n, 1, As_array_, n, pivots_, bs_array_, n, info_, m));

  // Copy solution from device to pinned host memory
  CHECK_CUDA(cudaMemcpy(xs, bs_, m * n * sizeof(double), cudaMemcpyDeviceToHost));

  // Free GPU memory
  cudaFree(As_);
  cudaFree(bs_);
  cudaFree(pivots_);
  cudaFree(info_);
  cudaFree(As_array_);
  cudaFree(bs_array_);

  // Free pinned host memory
  cudaFreeHost(As_array);
  cudaFreeHost(bs_array);

  cublasDestroy(handle_);
}
// --------------------------------------------------------------------
int32_t main() {
  int32_t n = 3; // Matrix size
  int32_t m = 2; // Number of systems

  double *As, *bs, *xs;

  // Allocate pinned host memory
  CHECK_CUDA(cudaMallocHost((void**)&As, m * n * n * sizeof(double)));
  CHECK_CUDA(cudaMallocHost((void**)&bs, m * n * sizeof(double)));
  CHECK_CUDA(cudaMallocHost((void**)&xs, m * n * sizeof(double)));

  // Initialize A and B
  double A_data[] = {
    4, 1, 2,   // System 1: A
    1, 3, 2,
    2, 2, 3,

    5, 2, 1,   // System 2: A
    2, 6, 3,
    1, 3, 7
  };

  double B_data[] = {
    1, 2, 3,  // System 1: b
    4, 5, 6   // System 2: b
  };

  memcpy(As, A_data, m * n * n * sizeof(double));
  memcpy(bs, B_data, m * n * sizeof(double));

  // Solve the linear systems
  solveBatchedLinearSystems(n, m, As, bs, xs);

  // Print solutions
  std::cout << "Solutions (X):\n";
  for (int32_t i = 0; i < m; i++) {
    std::cout << "System " << i + 1 << ": ";
    for (int32_t j = 0; j < n; j++) {
      std::cout << xs[i * n + j] << " ";
    }
    std::cout << "\n";
  }

  // Free pinned host memory
  cudaFreeHost(As);
  cudaFreeHost(bs);
  cudaFreeHost(xs);

  return 0;
}
// --------------------------------------------------------------------