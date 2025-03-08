#define CL_TARGET_OPENCL_VERSION 300

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
// ---------------------------------------------------------------------
// sudo apt update
// sudo apt install ocl-icd-opencl-dev opencl-headers clinfo
// 
// clinfo | grep 'Device Name'
// ---------------------------------------------------------------------
// gcc opencex.c -o opencex -lOpenCL
// ---------------------------------------------------------------------
#define ARRAY_SIZE 5

// 🦠 device kernel 
const char *kernelio = 
  "__kernel void doubleArray(__global int *arr) {"
  "   int id = get_global_id(0);"
  "   arr[id] *= 2;"
  "}";

int main() {
  // 🖊️
  int data[ARRAY_SIZE] = {1, 2, 3, 4, 5};
  size_t dataSize = sizeof(data);

  // get 🦠
  cl_platform_id platform;
  cl_device_id device;
  clGetPlatformIDs(1, &platform, NULL);
  clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);

  // opencl context & command queue
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

  // 🖨️
  printf("Output: ");
  for (int i = 0; i < ARRAY_SIZE; i++)
    printf("%d ", data[i]);
  printf("\n");

  return 0;
}
