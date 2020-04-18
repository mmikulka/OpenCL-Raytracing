#include <iostream>
#include <vector>
#include <OpenCL/cl.hpp> // main OpenCL include file

#include "CLStructs.hpp"
#include "gfxraytrace.hpp"
#include "CLStructs.hpp"

#define NUM_GLOBAL_WITEMS 1024

const int n = 8*32*512;             // size of vectors
const int k = 10000;                // number of loop iterations

cl_float2* cl_uv(cl_float3 positions[], vp &viewport, int numPixels)
{
    //setup and get platforms of computer
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);

    if (all_platforms.size()==0) {
        std::cout<<" No platforms found. Check OpenCL installation!\n";
        exit(1);
    }
    //grab only first platform
    cl::Platform default_platform=all_platforms[0];
    // std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";

    // get default device (CPUs, GPUs) of the default platform
    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if(all_devices.size()==0){
        std::cout<<" No devices found. Check OpenCL installation!\n";
        exit(1);
    }

    // use device[1] because that's a GPU; device[0] is the CPU
    cl::Device default_device=all_devices[1];
    // std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";

    cl::Context context({default_device});
    cl::Program::Sources sources;

    // calculates for each element; C = A + B
    char * kernel_code = nullptr;
       int filesize = 0;
           std::ifstream myfile("kernel.cl");
           myfile.seekg(0, myfile.end);
           filesize = myfile.tellg();
           myfile.seekg(0, std::ios::beg);
           
           kernel_code = new char[filesize];
           myfile.read(kernel_code, filesize);
           myfile.close()
    ;
    sources.push_back({kernel_code, filesize});
    

    cl::Program program(context, sources);
    if (program.build({default_device}) != CL_SUCCESS) {
        std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    // set up kernels and vectors for GPU code
    cl::CommandQueue queue(context, default_device);
    cl::Kernel uv_kernel = cl::Kernel(program, "uv");

    // construct vectors
    cl_float2 C[n];

    //put aside memory for buffer
    cl::Buffer buffer_A2(context, CL_MEM_READ_WRITE, sizeof(cl_float3)*numPixels);
    cl::Buffer buffer_B2(context, CL_MEM_READ_WRITE, sizeof(vp));
    cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(cl_float3)*numPixels);
    
    //write arrays to buffer
    queue.enqueueWriteBuffer(buffer_A2, CL_TRUE, 0, sizeof(cl_float3)*numPixels, positions);
    queue.enqueueWriteBuffer(buffer_B2, CL_TRUE, 0, sizeof(vp), &viewport);

    uv_kernel.setArg(0, buffer_A2);
    uv_kernel.setArg(1, buffer_B2);
    uv_kernel.setArg(2, buffer_C2);
    
    queue.enqueueNDRangeKernel(uv_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(cl_float2)*n, C);
    
    queue.finish();
    
    delete[]kernel_code;
    
    for (int i = 0; i < 10; ++i)
    {
        
        std::cout << viewport.left << ", " << C[i].s[1] << std::endl;
    }
    
    return C;
    
}

