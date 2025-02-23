// --------------------------------------------------------------------
// m systems of equations of the form:
//      Ax=b, A of size n×n
// 
// * store in memory
//   all As and bs in continuous arrays
// 
// * factorize
//   cholesky (positive definite?)
//     cublasDpotrfBatched()
//   LU
//     cublasDgetrfBatched()
// 
// * solve 
//   cholesky
//     cublasDpotrsBatched()
//   LU
//     cublasDgetrsBatched()
// 
// * retreive solution
// 
// * clean gpu
// --------------------------------------------------------------------
// 🐦‍⬛🧠
double* As_; // Batched matrix A
double* bs_; // Batched RHS vector B
int32_t* info_; // Status info for factorization
double** As_array_;
double** bs_array_;

cudaMalloc((void**)&As_, m * n * n * sizeof(double));
cudaMalloc((void**)&bs_, m * n * sizeof(double));
cudaMalloc((void**)&info_, m * sizeof(int32_t));
cudaMalloc((void**)&As_array_, m * sizeof(double*));
cudaMalloc((void**)&bs_array_, m * sizeof(double*));

double **As_array, **bs_array;
cudaMallocHost((void**)&As_array, m * sizeof(double*));
cudaMallocHost((void**)&bs_array, m * sizeof(double*)); 

cublasHandle_t handle_;
cublasCreate(&handle_);

// 🐦‍⬛ ⟵ 💻
cudaMemcpy(As_, As, m * n * n * sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(bs_, bs, m * n * sizeof(double), cudaMemcpyHostToDevice);

for (int32_t i = 0; i < m; i++) {
  As_array[i] = As_ + i * n * n;
  bs_array[i] = bs_ + i * n;
}

cudaMemcpy(As_array_, As_array, m * sizeof(double*), cudaMemcpyHostToDevice);
cudaMemcpy(bs_array_, bs_array, m * sizeof(double*), cudaMemcpyHostToDevice);

// 🚀🐦‍⬛

// cholesky factorization in batch
cublasDpotrfBatched(handle_, CUBLAS_FILL_MODE_LOWER, n, As_array_, n, info_, m);
// solve Ax = b using Cholesky factors
cublasDpotrsBatched(handle_, CUBLAS_FILL_MODE_LOWER, n, 1, As_array_, n, bs_array, n, info_, m);

// LU
int32_t* d_pivots; // Pivot array for LU
cudaMalloc((void**)&d_pivots, m * n * sizeof(int32_t));
cublasDgetrfBatched(handle_, n, As_array_, n, d_pivots, info_, m);
cublasDgetrsBatched(handle_, CUBLAS_OP_N, n, 1, As_array_, n, d_pivots, bs_array, n, info_, m);

// 💻 ⟵ 🐦‍⬛
cudaMemcpy(xs, bs_, m * n * sizeof(double), cudaMemcpyDeviceToHost);

// 🚿🐦‍⬛
cudaFree(As_);
cudaFree(bs_);
cudaFree(info_);
cudaFree(As_array_);
cudaFree(bs_array);
cudaFree(d_pivots);
cublasDestroy(handle_);




