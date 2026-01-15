#define CL_TARGET_OPENCL_VERSION 300

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
// ---------------------------------------------------------------------
// sudo apt update
// sudo apt install ocl-icd-opencl-dev opencl-headers clinfo
// sudo apt install intel-gpu-tools
// 
// clinfo | grep 'Device Name'
// ---------------------------------------------------------------------
// gcc openclex.c -o openclex -lOpenCL && ./openclex
// 
// check 🦠:
//   sudo intel_gpu_top
// ---------------------------------------------------------------------
// benchmark 🦠:
//   glmark2
// ---------------------------------------------------------------------
#define ARRAY_SIZE 5

// 🦠 device kernel 
const char *kernelio = 
  "__kernel void doubleArray(__global int *arr) {"
  "   int id = get_global_id(0);"
  "   arr[id] *= 2;"
  "}";

int main()
{
  // 🖊️
  int data[ARRAY_SIZE] = {1, 2, 3, 4, 5};
  size_t dataSize = sizeof(data);
  // -------------------------------------------------------------------
  // get 🦠
  cl_platform_id platform;
  cl_device_id device;
  if (clGetPlatformIDs(1, &platform, NULL) != CL_SUCCESS)
  {
    printf("failed to get platform id 😥\n");
    return 1;
  }
  if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL) != CL_SUCCESS)
  {
    printf("failed to get device id 😥\n");
    return 1;
  }
  char devicename[256];
  clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(devicename), devicename, NULL);
  char deviceversion[256];
  clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(deviceversion), deviceversion, NULL);

  // 🖨️🦠
  printf("device name .................... %s \n", devicename);
  printf("device version ................. %s \n", deviceversion);

  // 🧠🦠
  cl_ulong global_mem = 0;
  cl_ulong max_alloc = 0;
  cl_ulong local_mem = 0;
  clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem), &global_mem, NULL);
  clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(max_alloc), &max_alloc, NULL);
  clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem), &local_mem, NULL);

  // 🖨️🦠
  printf("global memory .................. %llu bytes (%.2f gb)\n", (unsigned long long)global_mem, (double)global_mem / (1024.0 * 1024.0 * 1000));
  printf("max allocable single ........... %llu bytes (%.2f gb)\n", (unsigned long long)max_alloc, (double)max_alloc / (1024.0 * 1024.0 * 1000));
  printf("local memory (per work group) .. %llu bytes\n", (unsigned long long)local_mem);
  // -------------------------------------------------------------------
  // 🦠🦠 opencl context & command queue 🦠🦠
  cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, NULL);
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, 0, NULL);

  // 🦠🧠 & 💻 ⟶ 🦠
  cl_mem buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, dataSize, data, NULL);

  // 🦠⬛
  cl_program program = clCreateProgramWithSource(context, 1, &kernelio, NULL, NULL);
  clBuildProgram(program, 1, &device, NULL, NULL, NULL);
  cl_kernel kernel = clCreateKernel(program, "doubleArray", NULL);
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &buffer);

  // 🚀🦠
  size_t globalSize = ARRAY_SIZE;
  clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalSize, NULL, 0, NULL, NULL);
  
  // 💻 ⟵ 🦠
  clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, dataSize, data, 0, NULL, NULL);

  // 🚿🦠
  clReleaseMemObject(buffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
  // -------------------------------------------------------------------
  // 🖨️
  printf("output: ");
  for (int i = 0; i < ARRAY_SIZE; i++)
    printf("%d ", data[i]);
  printf("\n");

  return 0;
}
