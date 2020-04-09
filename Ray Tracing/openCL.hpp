////If you want to build the file directly at the command prompt then use the following commands.
////AMD commands
////cl /c saxpy.cpp /I"%AMDAPPSDKROOT%\include"
////link  /OUT:"saxpy.exe" "%AMDAPPSDKROOT%\lib\x86_64\OpenCL.lib" saxpy.obj
//
////nVIDIA commands
////cl /c saxpy.cpp /I"%NVSDKCOMPUTE_ROOT%\OpenCL\common\inc"
////link  /OUT:"saxpy.exe" "%NVSDKCOMPUTE_ROOT%\OpenCL\common\lib\x64\OpenCL.lib" saxpy.obj
//
//#include <iostream>
//#include <array>
//#include <string>
//#include <fstream>
//#include <OpenCL/cl.hpp>
//
//#define VECTOR_SIZE 4096
//
////OpenCL kernel which is run for every work item created.
////The below const char string is compiled by the runtime complier
////when a program object is created with clCreateProgramWithSource
////and built with clBuildProgram.
//const std::string saxpy_kernel =
//"__kernel                                   \n"
//"void saxpy_kernel(                         \n"
//"                  __global float *A,       \n"
//"                  __global float *C)       \n"
//"{                                          \n"
//"    //Get the index of the work-item       \n"
//"    int index = get_global_id(0);          \n"
//"    C[index] = 2 * A[index]; \n"
//"}                                          \n";
//
////const std::string mykernel =
////"__kernel void NumericalReduction(__global int* data, __local int* localData, __global int *outData)\n"
////"{\n"
////"    size_t globalId = get_global_id(0);\n"
////"    size_t localSize = get_local_size(0);\n"
////"    size_t localId = get_local_id(0);\n"
////"\n"
////"    localData[localId] = data[globalId];\n"
////"    barrier(CLK_LOCAL_MEM_FENCE);\n"
////"\n"
////"    for (int i = localSize >> 1; i > 0 ; i >>= 1)\n"
////"    {\n"
////"        if (localId < i)\n"
////"            localData[localId] += localData[localId + i];\n"
////"        barrier(CLK_LOCAL_MEM_FENCE);\n"
////"    }\n"
////"    outData[get_group_id(0)] = 20;\n"
////"}\n";
//
//void init_openCL(){
//    std::vector<cl::Platform> platforms;
//    cl::Platform::get(&platforms);
//
//    auto platform = platforms.front();
//    std::vector<cl::Device> devices;
//    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
//
//    auto device = devices.front();
//
//    std::ifstream helloworldFile("kernel.cl");
//    std::string src(std::istreambuf_iterator<char>(helloworldFile), (std::istreambuf_iterator<char>()));
//
//    cl::Program::Sources sources(1, std::make_pair(saxpy_kernel.c_str(), saxpy_kernel.length() + 1));
//
//    cl::Context context(device);
//    cl::Program program(context, sources);
//
//    auto err = program.build("-cl-std=CL1.2");
//
//    std::vector<int> vec(1024);
//
//    for (int i = 0; i < vec.size(); ++i)
//    {
//        vec[i] = i;
//    }
//
//    cl::Kernel kernel(program, "saxpy_kernel");
//
//    auto maxWorkGroupSize = kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device, &err);
//
//    auto workGroupSize = maxWorkGroupSize;
//
//    workGroupSize = 1024;
//
//    auto numWorkGroups = vec.size() / workGroupSize;
//
//    cl::Buffer buf(context, CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR, sizeof(int) * vec.size(), vec.data());
//    cl::Buffer outBuf(context, CL_MEM_WRITE_ONLY, sizeof(int) * numWorkGroups);
//
//    err = kernel.setArg(0, buf);
//    err = kernel.setArg(1, outBuf);
//
//    std::vector<int> outVec(numWorkGroups);
//
//    cl::CommandQueue queue(context, device);
//    err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(vec.size()));
//    err = queue.enqueueReadBuffer(outBuf, CL_FALSE, 0, sizeof(int) * outVec.size(), outVec.data());
//
//
//    cl::finish();
//
//    printf("%d\n", outVec[0]);
//
//}


#include <iostream>
#include <vector>
#include <OpenCL/cl.hpp> // main OpenCL include file

#define NUM_ELEMENTS

using namespace std;

typedef struct __attribute__((packed))_foo
{
    float x;
    float y;
}foo;

void init_openCL()
{
    // Find all available OpenCL platforms (e.g. AMD, Nvidia, Intel)
    vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    
    // Choose and create an OpenCL platform
    unsigned int input = 0;
    // Handle incorrect user input
    
    cl::Platform platform = platforms[input];
    
    // Print the name of chosen OpenCL platform
    cout << "Using OpenCL platform: \t" << platform.getInfo<CL_PLATFORM_NAME>() << endl;
    
    // Find all available OpenCL devices (e.g. CPU, GPU or integrated GPU)
    vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    
    input = 0;
    
    cl::Device device = devices[input];
    
    // Print the name of the chosen OpenCL device
    cout << endl << "Using OpenCL device: \t" << device.getInfo<CL_DEVICE_NAME>() << endl << endl;
    
    // Create an OpenCL context on that device.
    // the context manages all the OpenCL resources
    cl::Context context = cl::Context(device);
    
    ///////////////////
    // OPENCL KERNEL //
    ///////////////////
    
    const char* source_string =
    "typedef struct __attribute__((packed))_foo"
    "{"
    "    float x;"
    "    float y;"
    "}foo;"
    ""
    "__kernel void parallel_add(__global foo* arr, __global float* z){ "
    " const int i = get_global_id(0); " // get a unique number identifying the work item in the global pool
    " z[i] = arr[i].x + arr[i].y;    " // add two arrays
    "}";
    
    // Create an OpenCL program by performing runtime source compilation
    cl::Program program = cl::Program(context, source_string);
    
    // Build the program and check for compilation errors
    cl_int result = program.build({ device }, "");
    if (result) cout << "Error during compilation! (" << result << ")" << endl;
    
    // Create a kernel (entry point in the OpenCL source program)
    // kernels are the basic units of executable code that run on the OpenCL device
    // the kernel forms the starting point into the OpenCL program, analogous to main() in CPU code
    // kernels can be called from the host (CPU)
    cl::Kernel kernel = cl::Kernel(program, "parallel_add");
    
    // Create input data arrays on the host (= CPU)
    foo cpuArrayA[NUM_ELEMENTS];
    for(int i = 0; i < NUM_ELEMENTS; ++i)
    {
        cpuArrayA[i].x = i * 1.0f;
        cpuArrayA[i].y = i * 0.1f;
    }
    float cpuOutput[NUM_ELEMENTS] = {}; // empty array for storing the results of the OpenCL program
    
    // Create buffers (memory objects) on the OpenCL device, allocate memory and copy input data to device.
    // Flags indicate how the buffer should be used e.g. read-only, write-only, read-write
    cl::Buffer clBufferA = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, NUM_ELEMENTS * sizeof(foo), cpuArrayA);
//    cl::Buffer clBufferB = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, numElements * sizeof(cl_int), cpuArrayB);
    cl::Buffer clOutput = cl::Buffer(context, CL_MEM_WRITE_ONLY, NUM_ELEMENTS * sizeof(cl_int), NULL);
    
    // Specify the arguments for the OpenCL kernel
    // (the arguments are __global float* x, __global float* y and __global float* z)
    kernel.setArg(0, clBufferA); // first argument
    //kernel.setArg(1, clBufferB); // second argument
    kernel.setArg(1, clOutput);  // third argument
    
    // Create a command queue for the OpenCL device
    // the command queue allows kernel execution commands to be sent to the device
    cl::CommandQueue queue = cl::CommandQueue(context, device);
    
    // Determine the global and local number of "work items"
    // The global work size is the total number of work items (threads) that execute in parallel
    // Work items executing together on the same compute unit are grouped into "work groups"
    // The local work size defines the number of work items in each work group
    // Important: global_work_size must be an integer multiple of local_work_size
    std::size_t global_work_size = NUM_ELEMENTS;
    std::size_t local_work_size = 10; // could also be 1, 2 or 5 in this example
    // when local_work_size equals 10, all ten number pairs from both arrays will be added together in one go
    
    // Launch the kernel and specify the global and local number of work items (threads)
    queue.enqueueNDRangeKernel(kernel, NULL, global_work_size, local_work_size);
    
    // Read and copy OpenCL output to CPU
    // the "CL_TRUE" flag blocks the read operation until all work items have finished their computation
    queue.enqueueReadBuffer(clOutput, CL_TRUE, 0, NUM_ELEMENTS * sizeof(cl_float), cpuOutput);
    
    // Print results to console
    for (int i = 0; i < NUM_ELEMENTS; i++)
        cout << cpuArrayA[i].x << " + " << cpuArrayA[i].y << " = " << cpuOutput[i] << endl;
    
}
