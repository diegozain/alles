#include <iostream>
#include <string>
#include <cuda_runtime.h>
#include <cublas_v2.h>
// --------------------------------------------------------------------
// 
// 🦠💻🧠
// ** pinned memory **
// double *As;
// cudaMallocHost(&As, m * n * n * sizeof(double));
// cudaFreeHost(As);
// 
// 💻🧠
// ** heap memory **
// double* As = (double*)malloc(m * n * n *  sizeof(double));
// free(As);
// 
// --------------------------------------------------------------------
// 🔵🔵 windows 🔵🔵
// nvcc -o Axb.exe Axb.cu -I. -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\include" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\lib\x64" -lcudart -lcublas --gpu-architecture=sm_75 --expt-relaxed-constexpr -O3 -dlto -use_fast_math
// --------------------------------------------------------------------
void Axbsolver(int32_t n, int32_t m, double* As, double* bs, double* xs) {
  // ⌚
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  // --------------------------------------------------------------------
  // 🐦‍⬛🧠
  cublasHandle_t handle_;
  cublasCreate(&handle_);

  double *As_, *bs_;
  double **As_array_, **bs_array_;
  int32_t *info_; // for cublasDgetrfBatched
  int32_t *pivots_;

  cudaMalloc(&As_, m * n * n * sizeof(double));
  cudaMalloc(&bs_, m * n * sizeof(double));
  cudaMalloc(&As_array_, m * sizeof(double*));
  cudaMalloc(&bs_array_, m * sizeof(double*));
  cudaMalloc(&info_, m * sizeof(int32_t));
  cudaMalloc(&pivots_, m * n * sizeof(int32_t));
  // --------------------------------------------------------------------
  // 🦠💻🧠 pinned memory
  double **As_array, **bs_array;
  cudaMallocHost(&As_array, m * sizeof(double*));
  cudaMallocHost(&bs_array, m * sizeof(double*));
  int32_t *info; // for cublasDgetrsBatched
  cudaMallocHost(&info, m * sizeof(int32_t*));
  // --------------------------------------------------------------------
  // 🦠 ⟵ 💻
  cudaMemcpy(As_, As, m * n * n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(bs_, bs, m * n * sizeof(double), cudaMemcpyHostToDevice);

  // 🐦‍⬛ ⟵ 💻
  for (int32_t i = 0; i < m; i++) {
    As_array[i] = As_ + i * n * n;
    bs_array[i] = bs_ + i * n;
  }
  cudaMemcpy(As_array_, As_array, m * sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(bs_array_, bs_array, m * sizeof(double*), cudaMemcpyHostToDevice);
  // --------------------------------------------------------------------
  // 🚀🐦‍⬛
  
  // ∘∘ LU ∘∘
  cublasDgetrfBatched(handle_, n, As_array_, n, pivots_, info_, m);
  cublasDgetrsBatched(handle_, CUBLAS_OP_N, n, 1, As_array_, n, pivots_, bs_array_, n, info, m);
  // --------------------------------------------------------------------
  // 💻 ⟵ 🐦‍⬛
  cudaMemcpy(xs, bs_, m * n * sizeof(double), cudaMemcpyDeviceToHost);
  // --------------------------------------------------------------------
  // 🚿🦠🐦‍⬛
  cublasDestroy(handle_);
  
  cudaFree(As_);
  As_ = nullptr;
  cudaFree(bs_);
  bs_ = nullptr;
  cudaFree(info_);
  info_ = nullptr;
  cudaFree(As_array_);
  As_array_ = nullptr;
  cudaFree(bs_array_);
  bs_array_ = nullptr;
  cudaFree(pivots_);
  pivots_ = nullptr;

  // 🚿🦠🐦‍⬛💻
  cudaFreeHost(As_array);
  cudaFreeHost(bs_array);
  cudaFreeHost(info);
  // --------------------------------------------------------------------
  // ⌚
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTime;
  cudaEventElapsedTime(&elapsedTime, start, stop);
  // ⌚🚿
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  std::cout << "  Axbsolver elapsed time: " << elapsedTime/1e3 << " seconds" << "\n";
}
// --------------------------------------------------------------------
int32_t main() {
  // 🧠
  int32_t n = 5;
  int32_t m = 3;

  // 🦠💻🧠 pinned memory
  double *As, *bs, *xs;
  cudaMallocHost(&As, m * n * n * sizeof(double));
  cudaMallocHost(&bs, m * n * sizeof(double));
  cudaMallocHost(&xs, m * n * sizeof(double));

  // 🖊️
  double* As_data = (double*)malloc(m * n * n *  sizeof(double));
  double* bs_data = (double*)malloc(m * n *  sizeof(double));
  for (int32_t i = 0; i < m * n * n; i++) {
    As_data[i] = static_cast<double>(rand() % 10);
  }
  for (int32_t i = 0; i < m * n; i++) {
    bs_data[i] = static_cast<double>(rand() % 10);
  }

  memcpy(As, As_data, m * n * n * sizeof(double));
  memcpy(bs, bs_data, m * n * sizeof(double));

  // 🖨️
  std::cout << "\n\n ::: solving now :::" << "\n\n";
  for (int32_t j = 0; j < m; j++) {
    std::cout << "A" << j + 1 << ": \n";
    for (int32_t ii = 0; ii < n; ii++) {
      for (int32_t i = 0; i < n; i++) {
        std::cout << As[j*n*n + ii * n + i] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  for (int32_t i = 0; i < m; i++) {
    std::cout << "bs" << i + 1 << ": ";
    for (int32_t j = 0; j < n; j++) {
      std::cout << bs[i * n + j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n\n";

  // 🚀🐦‍⬛
  Axbsolver(n, m, As, bs, xs);

  // 🖨️
  std::cout << "\n\n ::: result :::" << "\n\n";
  for (int32_t i = 0; i < m; i++) {
    std::cout << "xs" << i + 1 << ": ";
    for (int32_t j = 0; j < n; j++) {
      std::cout << xs[i * n + j] << " ";
    }
    std::cout << "\n";
  }

  // 🚿🦠💻
  cudaFreeHost(As);
  cudaFreeHost(bs);
  cudaFreeHost(xs);

  // 🚿
  free(As_data);
  As_data = nullptr;
  free(bs_data);
  bs_data = nullptr;

  return 0;
}
// --------------------------------------------------------------------