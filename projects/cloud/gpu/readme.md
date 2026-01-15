# gpu trips

diego domenzain

12.2025

---

## opencl

* ``kernelios.cl`` 
  * hosts all the device functions,
  * it is compiled with ``ocloc`` intel, or ``cloc`` amd.
* ``file.c``
  * hosts all the host functions,
  * it is compiled with ``gcc -lOpenCL``,
  * the host code reads by name the binary of ``kernelios.cl``.

get all the stuff (*Intel*),
```shell
sudo apt update
sudo apt install ocl-icd-opencl-dev opencl-headers clinfo
sudo apt install intel-gpu-tools

clinfo | grep 'Device Name'
```

compile,
```shell
ocloc compile -file kernelios.cl -device acm-g10
gcc openclex_.c -o openclex_ -lOpenCL && ./openclex_
```

check if it's working (*Intel*),
```shell
sudo intel_gpu_top
```

### ocloc

needs "device name". the way to figure that out is as follows,
1. run ``clinfo | grep 'Device Name'`` and see something like ``Intel(R) Graphics [0x46a8]``.
1. run ``ocloc compile --help`` and get the names of available devices.
1. ask an llm about it using info from previous steps.

``opencl`` needs special enabling for *floating point* support. you can see which ``fp`` are available in your device by running,
```shell
clinfo | grep cl_khr_fp
```

## cuda & opencl

⬟  host cuda
```c++
int nthreinbloc= 1024;
int nbloingrid = 2;
kernel<<<nbloingrid, nthreinbloc>>>();
```

⌲ device cuda
```c++
__global__ void kernel()
{
  int iblock = blockIdx.x;
  int nthr = blockDim.x;
  int ithread = threadIdx.x;
  int ithread_= blockIdx.x * blockDim.x + threadIdx.x;
}
```

⬟  host opencl
```c
size_t nthreinbloc= 1024;
size_t nbloingrid = 2;
size_t nglobalwork= nthreinbloc * nbloingrid;

clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &nglobalwork, &nthreinbloc, 0, NULL, NULL);
```

⌲ device opencl
```c
__kernel void kernel()
{
  int iblock = get_group_id(0);
  int nthr = get_local_size(0);
  int ithread = get_local_id(0);;
  int ithread_= get_global_id(0);
}
```

### in more detail

```text
CUDA	                           OpenCL
blocks ........................... work-gropus
threads .......................... work-items

grid of blocks ................... global work size (NDRange)
threads per block ................ local work size
thread ID (threadIdx) ............ local ID (get_local_id)
block ID (blockIdx) .............. group ID (get_group_id)
block size (blockDim) ............ local size (get_local_size)
grid size (gridDim) .............. number of work groups (get_num_groups)

blockDim=nthreinbloc

blockIdx * blockDim + threadIdx .. get_global_id(dim)
threadIdx ........................ get_local_id(dim)
blockIdx ......................... get_group_id(dim)
blockDim ......................... get_local_size(dim)
gridDim .......................... get_num_groups(dim)
                                 . dim=0,1,2.
```

cuda:
```c++
dim3 nbloingrid(32, 32);
dim3 nthreinbloc(16, 16);
kernel<<<nbloingrid, nthreinbloc>>>();
// kernel<<<nfwd, nthreinbloc>>>();
```
opencl:
```c
// threads per work-group (blockDim)
size_t nthreinbloc[2] = {16, 16};
// total threads (gridDim * blockDim)
size_t global_work_size[2] = {32 * 16, 32 * 16};
// queue kernel
clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global_work_size, nthreinbloc, 0, NULL, NULL);
```

opencl:
```c
// dim=1
size_t global_work_size = 1024;
clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
// dim=2
size_t global_work_size[2] = {1024, 512};
size_t nthreinbloc[2] = {16, 16};
clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global_work_size, nthreinbloc, 0, NULL, NULL);
```

## cuda

host,
```c++
// ⌲ 𝔹
double *var_;
cudaMalloc(&var_, nvar * sizeof(double));
// ⬟ ⟶ ⌲
cudaMemcpy(var_, var, nvar * sizeof(double), cudaMemcpyHostToDevice);
// ⌲ ℧
// 1 << a ⟶ 2ᵃ
// a=9 ⟶ 512, a=10 ⟶ 1024
nthreinbloc= 1 << 10;
if (nthingstodo < nthreinbloc)
{
  nthreinbloc = nthingstodo;
}
// nbloingrid = ⌈ nthingstodo ÷ nthreinbloc ⌉
nbloingrid = static_cast<int32_t>((nthingstodo + nthreinbloc - 1) / nthreinbloc);
// ∗ ⌲
kernel<<<nbloingrid, nthreinbloc>>>(var_, nvar);
// ⬟ ⟵ ⌲
cudaMemcpy(var, var_, nvar * sizeof(double), cudaMemcpyDeviceToHost);
// ⌲ ⫤
cudaFree(var_);
var_ = nullptr;
```

device,
```c++
__global__ void kernel(double *var_, int32_t nvar)
{
  // blockIdx.x  • block index in the x-dimension
  // blockDim.x  • number of threads per block in the x-dimension
  // threadIdx.x • thread index within the current block in the x-dimension
  // gridDim.x   • number of blocks in the x-dimension

  // ℧ get thread index
  int32_t ithread_ = blockIdx.x * blockDim.x + threadIdx.x;
  // break down ithread_ to the indicies you want.

  // ⌲ are loco 𝒞
  if (indexnew is bad)
  {
    return;
  }

  // ℧
  double *varr_ = var_ + indexnew*noffset;

  // ∗
  // do stuff with varr_
}
``
