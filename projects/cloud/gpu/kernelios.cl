// ---------------------------------------------------------------------
// ocloc .. intel
// cloc ... amd
// ---------------------------------------------------------------------
// finding the device name for ocloc:
// 
// clinfo | grep 'Device Name'
//   Intel(R) Graphics [0x46a8]
// 
// get the names of available devices in ocloc:
//   ocloc compile --help
// 
// ask an llm about it.
// 
// turns out that gpu is Intel Arc (alchemist) 
//   ---> device name = acm-g10
// ---------------------------------------------------------------------
// ocloc compile -file kernelios.cl -device acm-g10
// ---------------------------------------------------------------------
#pragma OPENCL EXTENSION cl_khr_fp16 : enable
// ---------------------------------------------------------------------
inline half addonefp16(half value) 
{
  return value + (half)1;
}
// ---------------------------------------------------------------------
__kernel void doublearrayint(__global int *arr) 
{
  int id = get_global_id(0);
  arr[id] *= 2;
}
__kernel void halfarrayfp16(__global half *arr) 
{
  int id = get_global_id(0);
  arr[id] = addonefp16(arr[id]);
  arr[id] /= (half)2;
}
// ---------------------------------------------------------------------
