#include <iostream>
#include <vector>
#include <OpenCL/cl.hpp> // main OpenCL include file

#include "CLStructs.hpp"
#include "gfxraytrace.hpp"

cl::Program init_openCL()
{
    // Find all available OpenCL platforms (e.g. AMD, Nvidia, Intel)
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    
    // Choose and create an OpenCL platform
    unsigned int input = 0;
    
    cl::Platform platform = platforms[input];
    
    // Print the name of chosen OpenCL platform
    std::cout << "Using OpenCL platform: \t" << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
    
    // Find all available OpenCL devices (e.g. CPU, GPU or integrated GPU)
    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    
    if (devices.size() > 1) // set device to descrite graphics card if there is one
        input = 1;
    
    cl::Device device = devices[input];
    
    // Print the name of the chosen OpenCL device
    std::cout << std::endl << "Using OpenCL device: \t" << device.getInfo<CL_DEVICE_NAME>() << std::endl << std::endl;
    
    // Create an OpenCL context on that device.
    cl::Context context = cl::Context(device);
    
    //load in kernel programs
    char* source_string = nullptr;
    std::ifstream myfile ("kernel.cl");
    if (myfile.is_open())
    {
        myfile.seekg(0, myfile.end);
        int filesize = myfile.tellg();
        myfile.seekg(0, std::ios::beg);
       
        source_string = new char[filesize];
        myfile.read(source_string, filesize);
        myfile.close();
                
        // Create an OpenCL program by performing runtime source compilation
        cl::Program program = cl::Program(context, source_string);
        
        // Build the program and check for compilation errors
        cl_int result = program.build({ device }, "");
        if (result) std::cout << "Error during compilation! (" << result << ")" << std::endl;
        
        delete[] source_string;
        
        return program;
    }
    else{
        std::cout << "Could not open file: \"kernel.cl\"" << std::endl;
        return cl::Program();
    }
}

cl_float2* cl_uv(cl_float3 * positions, vp* viewport, int numPixels)
{
    cl::Program program = init_openCL();
    
    auto context = program.getInfo<CL_PROGRAM_CONTEXT>();
    auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
    auto& device = devices.front();
    
    // Create a kernel
    cl::Kernel kernel = cl::Kernel(program, "uv");

    cl_float2* output = new cl_float2[numPixels]; // empty array for storing the results of the OpenCL program
    for (int i = 0; i < 10; ++i)
    {
        output[i].s[0] = 1.0f;
        output[i].s[1] = 2.0f;
    }
    
    /*for (int i = 0; i < 10; ++i)
    {
        std::cout << "pre CL operation: " << output[i].s[0] << ", " << output[i].s[1] << std::endl;
    }*/
    
    // Create buffers (memory objects) on the OpenCL device, allocate memory and copy input data to device.
    cl::Buffer clBufferA = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, numPixels * sizeof(cl_float3), positions);
    cl::Buffer clBufferB = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(vp), viewport);
    cl::Buffer clOutput = cl::Buffer(context, CL_MEM_WRITE_ONLY, numPixels * sizeof(cl_float2), output);
    
    // Specify the arguments for the OpenCL kernel
    kernel.setArg(0, clBufferA); // first argument
    kernel.setArg(1, clBufferB); // second argument
    kernel.setArg(2, clOutput);  // third argument
    
    // Create a command queue for the OpenCL device
    cl::CommandQueue queue = cl::CommandQueue(context, device);
    
    // Determine the global and local number of "work items"
    // The global work size is the total number of work items (threads) that execute in parallel
    // Work items executing together on the same compute unit are grouped into "work groups"
    // The local work size defines the number of work items in each work group
    // Important: global_work_size must be an integer multiple of local_work_size
    std::size_t global_work_size = numPixels;
    std::size_t local_work_size = 200;
    
    // Launch the kernel and specify the global and local number of work items (threads)
    queue.enqueueNDRangeKernel(kernel, NULL, global_work_size, local_work_size);
    
    // Read and copy OpenCL output to CPU
    // the "CL_TRUE" flag blocks the read operation until all work items have finished their computation
    queue.enqueueReadBuffer(clOutput, CL_TRUE, 0, numPixels * sizeof(cl_float2), output);
    
    return output;
    
}

viewRay* cl_ortho_viewrays(cam * camera, cl_float2* uv, int numPixels)
{
    cl::Program program = init_openCL();
    
    auto context = program.getInfo<CL_PROGRAM_CONTEXT>();
    auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
    auto& device = devices.front();
    
    // Create a kernel
    cl::Kernel kernel = cl::Kernel(program, "ortho_viewrays");

    viewRay* viewRays = new viewRay[numPixels]; // empty array for storing the results of the OpenCL program
    
    // Create buffers (memory objects) on the OpenCL device, allocate memory and copy input data to device.
    cl::Buffer clBufferCam = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cam), camera);
    cl::Buffer clBufferUV = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, numPixels * sizeof(UV), uv);
    cl::Buffer clOutput = cl::Buffer(context, CL_MEM_WRITE_ONLY, numPixels * sizeof(viewRay), viewRays);
    
    // Specify the arguments for the OpenCL kernel
    kernel.setArg(0, clBufferCam); // first argument
    kernel.setArg(1, clBufferUV); // second argument
    kernel.setArg(2, clOutput);  // third argument
    
    // Create a command queue for the OpenCL device
    cl::CommandQueue queue = cl::CommandQueue(context, device);
    
    // Determine the global and local number of "work items"
    // The global work size is the total number of work items (threads) that execute in parallel
    // Work items executing together on the same compute unit are grouped into "work groups"
    // The local work size defines the number of work items in each work group
    // Important: global_work_size must be an integer multiple of local_work_size
    std::size_t global_work_size = numPixels;
    std::size_t local_work_size = 10;
    
    // Launch the kernel and specify the global and local number of work items (threads)
    queue.enqueueNDRangeKernel(kernel, NULL, global_work_size, local_work_size);
    
    // Read and copy OpenCL output to CPU
    // the "CL_TRUE" flag blocks the read operation until all work items have finished their computation
    queue.enqueueReadBuffer(clOutput, CL_TRUE, 0, numPixels * sizeof(cl_float3), viewRays);
    
    return viewRays;
}

viewRay* cl_persp_viewrays(cam * camera, cl_float2* uv, int numPixels, float focal_length)
{
    cl::Program program = init_openCL();
    
    auto context = program.getInfo<CL_PROGRAM_CONTEXT>();
    auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
    auto& device = devices.front();
    
    // Create a kernel
    cl::Kernel kernel = cl::Kernel(program, "ortho_viewrays");

    viewRay* viewRays = new viewRay[numPixels]; // empty array for storing the results of the OpenCL program
    
    // Create buffers (memory objects) on the OpenCL device, allocate memory and copy input data to device.
    cl::Buffer clBufferCam = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cam), camera);
    cl::Buffer clBufferUV = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, numPixels * sizeof(UV), uv);
    cl::Buffer clOutput = cl::Buffer(context, CL_MEM_WRITE_ONLY, numPixels * sizeof(viewRay), viewRays);
    
    // Specify the arguments for the OpenCL kernel
    kernel.setArg(0, clBufferCam); // first argument
    kernel.setArg(1, clBufferUV); // second argument
    kernel.setArg(2, clOutput);  // third argument
    
    // Create a command queue for the OpenCL device
    cl::CommandQueue queue = cl::CommandQueue(context, device);
    
    // Determine the global and local number of "work items"
    // The global work size is the total number of work items (threads) that execute in parallel
    // Work items executing together on the same compute unit are grouped into "work groups"
    // The local work size defines the number of work items in each work group
    // Important: global_work_size must be an integer multiple of local_work_size
    std::size_t global_work_size = numPixels;
    std::size_t local_work_size = 10;
    
    // Launch the kernel and specify the global and local number of work items (threads)
    queue.enqueueNDRangeKernel(kernel, NULL, global_work_size, local_work_size);
    
    // Read and copy OpenCL output to CPU
    // the "CL_TRUE" flag blocks the read operation until all work items have finished their computation
    queue.enqueueReadBuffer(clOutput, CL_TRUE, 0, numPixels * sizeof(cl_float3), viewRays);
    
    return viewRays;
}

