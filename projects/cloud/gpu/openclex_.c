#define CL_TARGET_OPENCL_VERSION 300

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <CL/cl.h>
// ---------------------------------------------------------------------
// sudo apt update
// sudo apt install ocl-icd-opencl-dev opencl-headers clinfo
// sudo apt install intel-gpu-tools
// 
// clinfo | grep 'Device Name'
// ---------------------------------------------------------------------
// ocloc .. intel
// cloc ... amd
// ---------------------------------------------------------------------
// ocloc compile -file kernelios.cl -device acm-g10
// gcc openclex_.c -o openclex_ -lOpenCL && ./openclex_
// 
// check 🦠:
//   sudo intel_gpu_top
// ---------------------------------------------------------------------
#define ARRAY_SIZE 5

// ℕ ⟶ ℝ
float ui16fp16(uint16_t half)
{
  uint32_t sign = (half & 0x8000) << 16;
  uint32_t exponent = (half & 0x7C00) >> 10;
  uint32_t mantissa = (half & 0x03FF);

  if (exponent == 0x1F) // Inf or NaN
    exponent = 0xFF;
  else if (exponent == 0) // Subnormal or zero
    exponent = 0;
  else
    exponent = exponent + (127 - 15);

  uint32_t float_bits = sign | (exponent << 23) | (mantissa << 13);
  float result;
  memcpy(&result, &float_bits, sizeof(result));
  return result;
}
// ℝ ⟶ ℕ
uint16_t fp32ui16(float value)
{
  uint32_t float_bits;
  memcpy(&float_bits, &value, sizeof(float_bits));

  uint32_t sign = (float_bits >> 16) & 0x8000; // Extract sign bit
  uint32_t exponent = ((float_bits >> 23) & 0xFF) - 127 + 15; // Adjust exponent bias
  uint32_t mantissa = (float_bits & 0x007FFFFF) >> 13; // Truncate mantissa to 10 bits

  if (exponent <= 0) // Subnormal or zero
  {
    if (exponent < -10)
    {
      // Too small, becomes zero
      return sign;
    }
    // Subnormal number
    mantissa = (mantissa | 0x0800) >> (1 - exponent);
    return sign | mantissa;
  }
  else if (exponent >= 31) // Overflow, becomes Inf or NaN
  {
    return sign | 0x7C00 | (mantissa ? 0x01 : 0x00); // Inf or NaN
  }
  // Normalized number
  return sign | (exponent << 10) | mantissa;
}
// ---------------------------------------------------------------------
int main()
{
  // 🖊️
  size_t narray = ARRAY_SIZE;

  int dataint[ARRAY_SIZE] = {1, 2, 3, 4, 5};
  size_t ndataint = sizeof(dataint);

  uint16_t datafp16[ARRAY_SIZE] = {0x3C00, 0x4000, 0x4200, 0x4400, 0x4500};
  size_t ndatafp16 = sizeof(datafp16);

  size_t nthreinbloo, nbloingrii, nglobalwork;

  // ℝ ⟶ ℕ
  float datafp32[ARRAY_SIZE] = {1.0, 2.0, 3.0, 4.0, 5.0};
  for (int i=0; i<narray; i++)
  {
    datafp16[i] = fp32ui16(datafp32[i]);
  }
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
  // -------------------------------------------------------------------
  // 🦠🦠 opencl context & command queue 🦠🦠
  cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, NULL);
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, 0, NULL);
  // -------------------------------------------------------------------
  // 👁️🍍
  FILE *fp = fopen("kernelios_dg2.bin", "rb");
  fseek(fp, 0, SEEK_END);
  size_t binary_size = ftell(fp);
  rewind(fp);
  unsigned char *binary = (unsigned char *)malloc(binary_size);
  fread(binary, 1, binary_size, fp);
  fclose(fp);

  // 🦠🍍
  cl_int status;
  cl_program program = clCreateProgramWithBinary(context, 1, &device, &  binary_size, (const unsigned char**)&binary, NULL, &status);
  if (status != CL_SUCCESS)
  {
    printf("failed to create program with binary 😥\n");
    return 1;
  }
  free(binary);
  clBuildProgram(program, 1, &device, NULL, NULL, NULL);
  // -------------------------------------------------------------------
  // 🦠🔪🍍🧠 & 💻 ⟶ 🦠
  cl_mem dataint_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, ndataint, dataint, NULL);

  // 🦠🔪🍍
  cl_kernel doublearrayintker = clCreateKernel(program, "doublearrayint", NULL);
  clSetKernelArg(doublearrayintker, 0, sizeof(cl_mem), &dataint_);

  // 🚀🔪🍍⬛
  nthreinbloo= narray;
  nbloingrii = 1;
  nglobalwork= nthreinbloo * nbloingrii;
  clEnqueueNDRangeKernel(queue, doublearrayintker, 1, NULL, &nglobalwork, &nthreinbloo, 0, NULL, NULL);

  // 💻 ⟵ 🦠🔪🍍
  clEnqueueReadBuffer(queue, dataint_, CL_TRUE, 0, ndataint, dataint, 0, NULL, NULL);

  // 🚿🦠🔪🍍
  clReleaseMemObject(dataint_);
  clReleaseKernel(doublearrayintker);
  // -------------------------------------------------------------------
  // 🦠🔪🍍🧠 & 💻 ⟶ 🦠
  cl_mem datafp16_ = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, ndatafp16, datafp16, NULL);

  // 🦠🔪🍍
  cl_kernel halfarrayfp16ker = clCreateKernel(program, "halfarrayfp16", NULL);
  clSetKernelArg(halfarrayfp16ker, 0, sizeof(cl_mem), &datafp16_);

  // 🚀🔪🍍⬛
  nthreinbloo= narray;
  nbloingrii = 1;
  nglobalwork= nthreinbloo * nbloingrii;
  clEnqueueNDRangeKernel(queue, halfarrayfp16ker, 1, NULL, &nglobalwork, &nthreinbloo, 0, NULL, NULL);

  // 💻 ⟵ 🦠🔪🍍
  clEnqueueReadBuffer(queue, datafp16_, CL_TRUE, 0, ndatafp16, datafp16, 0, NULL, NULL);

  // 🚿🦠🔪🍍
  clReleaseMemObject(datafp16_);
  clReleaseKernel(halfarrayfp16ker);
  // -------------------------------------------------------------------
  // 🚿🦠
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
  // -------------------------------------------------------------------
  // 🖨️
  printf("output int: ");
  for (int i = 0; i < ARRAY_SIZE; i++)
  {
    printf("%d ", dataint[i]);
  }
  printf("\n");
  // 🖨️
  printf("output fp16: ");
  for (int i = 0; i < ARRAY_SIZE; i++)
  {
    printf("%f ", ui16fp16(datafp16[i]));
  }
  printf("\n");
  // -------------------------------------------------------------------
  // 🏠
  return 0;
}