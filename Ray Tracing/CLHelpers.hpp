#include <iostream>
#include <vector>
#include <OpenCL/cl.hpp> // main OpenCL include file

#include "CLStructs.hpp"
#include "gfxraytrace.hpp"
#include "CLStructs.hpp"

#define ARRAY_SPLIT 160000
#define NUM_GLOBAL_WITEMS ARRAY_SPLIT

const int n = 160000;             // size of arrays

void printfloat3(cl_float3 debug)
{
    std::cout << debug.s[0] << std::endl << debug.s[1] << std::endl <<  debug.s[2] << std::endl;
}

cl_float2* cl_uv(cl_float3* positions, vp &viewport, int numPixels)
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
    default_platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
    if(all_devices.size()==0){
        std::cout<<" No devices found. Check OpenCL installation!\n";
        exit(1);
    }
    
    // use device[1] because that's a GPU; device[0] is the CPU
    cl::Device default_device=all_devices[1];
    std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
    
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
    cl_float2* C = new cl_float2[n];
    
    //put aside memory for buffer
    cl::Buffer buffer_A2(context, CL_MEM_READ_WRITE, sizeof(cl_float3)*numPixels);
    cl::Buffer buffer_B2(context, CL_MEM_READ_WRITE, sizeof(vp));
    cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(cl_float2)*numPixels);
    
    //write arrays to buffer
    queue.enqueueWriteBuffer(buffer_A2, CL_TRUE, 0, sizeof(cl_float3)*numPixels, positions);
    queue.enqueueWriteBuffer(buffer_B2, CL_TRUE, 0, sizeof(vp), &viewport);
    
    uv_kernel.setArg(0, buffer_A2);
    uv_kernel.setArg(1, buffer_B2);
    uv_kernel.setArg(2, buffer_C2);
    
    queue.enqueueNDRangeKernel(uv_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(cl_float2)*numPixels, C);
    
    queue.finish();
    
    delete[]kernel_code;
    
    return C;
    
}

viewRay* cl_ortho_viewrays(cam& camera, cl_float2 * uV, int numPixels)
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
    
    // read kernel from kernel.cl
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
    cl::Kernel ortho_kernel = cl::Kernel(program, "ortho_viewrays");
    
    // construct arrays
    viewRay* C = new viewRay[n];
    
    //put aside memory for buffer
    //    cl::Buffer buffer_A2(context, CL_MEM_READ_WRITE, sizeof(cam));
    cl::Buffer buffer_B2(context, CL_MEM_READ_WRITE, sizeof(cl_float2) * numPixels);
    cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(viewRay)*numPixels);
    
    
    //write arrays to buffer
    //    queue.enqueueWriteBuffer(buffer_A2, CL_TRUE, 0, sizeof(cam), &camera);
    queue.enqueueWriteBuffer(buffer_B2, CL_TRUE, 0, sizeof(cl_float2) * numPixels, uV);
    
    ortho_kernel.setArg(0, camera);
    ortho_kernel.setArg(1, buffer_B2);
    ortho_kernel.setArg(2, buffer_C2);
    
    queue.enqueueNDRangeKernel(ortho_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(viewRay)*n, C);
    
    queue.finish();
    
    delete[]kernel_code;
    
    return C;
}

viewRay* cl_persp_viewrays(cam& camera, cl_float2 * uV, int numPixels, float focal_length)
{
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
    
    // read kernel from kernel.cl
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
    cl::Kernel persp_kernel = cl::Kernel(program, "persp_viewrays");
    
    // construct arrays
    viewRay* C = new viewRay[n];
    
    //put aside memory for buffer
    cl::Buffer buffer_B2(context, CL_MEM_READ_WRITE, sizeof(cl_float2) * numPixels);
    cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(viewRay)*numPixels);
    
    
    //write arrays to buffer
    queue.enqueueWriteBuffer(buffer_B2, CL_TRUE, 0, sizeof(cl_float2) * numPixels, uV);
    
    persp_kernel.setArg(0, camera);
    persp_kernel.setArg(1, buffer_B2);
    persp_kernel.setArg(2, buffer_C2);
    persp_kernel.setArg(3, focal_length);
    
    
    queue.enqueueNDRangeKernel(persp_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(viewRay)*n, C);
    
    queue.finish();
    
    delete[]kernel_code;
    
    //testing only
    //        for (int i = 0; i < numPixels; ++i)
    //        {
    //
    //            std::cout << "i: " << i << "\tRay direction: " << C[i].direction.s[0] << ", " << C[i].direction.s[1] << ", " << C[i].direction.s[2] << std::endl;
    //        }
    
    return C;
}

intersect* cl_intersect (object * objects, int numObjects, const viewRay* rays, int numRays, float t_upper_bound, float t_lower_bound)
{
    typedef float bug;
    
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
    cl::Kernel intersect_kernel = cl::Kernel(program, "intersections");
    
    // construct vectors
    intersect* C = new intersect[n];
    bug* debug  = new bug[n];
    
    //put aside memory for buffer
    cl::Buffer buffer_A2(context, CL_MEM_READ_WRITE, sizeof(object) * numObjects * 2);
    cl::Buffer buffer_B2(context, CL_MEM_READ_WRITE, sizeof(viewRay) * numRays);
    cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(intersect) * numRays);
    cl::Buffer buffer_debug(context, CL_MEM_READ_WRITE, sizeof(bug) * numRays);
    
    //write arrays to buffer
    queue.enqueueWriteBuffer(buffer_A2, CL_TRUE, 0, sizeof(object) * numObjects, objects);
    queue.enqueueWriteBuffer(buffer_B2, CL_TRUE, 0, sizeof(viewRay) * numRays, rays);
    
    intersect_kernel.setArg(0, buffer_A2);
    intersect_kernel.setArg(2, buffer_B2);
    intersect_kernel.setArg(3, buffer_C2);
    intersect_kernel.setArg(1, numObjects);
    intersect_kernel.setArg(4, t_upper_bound);
    intersect_kernel.setArg(5, t_lower_bound);
    intersect_kernel.setArg(6, buffer_debug);
    
    
    queue.enqueueNDRangeKernel(intersect_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(intersect)*n, C);
    queue.enqueueReadBuffer(buffer_debug, CL_TRUE, 0, sizeof(bug)*n, debug);
    
    queue.finish();
    
    delete[]kernel_code;
    for (int i = 0; i < numRays; ++i)
    {

        std::cout << std::endl << debug[0] << std::endl;
    }
    
    //
//        std::cout << "intersection T val" << std::endl;
//        for (int i = 0; i < numRays; ++i)
//        {
//            if (C[i].intersects > 0)
//            {
////            {
////            std::cout << "i = " << i; //<< "\t" <<debug[i] << std::endl;
////            printfloat3(C[i].normal);
////            }
//
//            std::cout <<  C[i].t_ << std::endl;
//            }
//        }
    //
    return C;
}

rgbColor* cl_flat_shader(const intersect* xsect, int numIntersections, rgbColor background)
{
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
    cl::Kernel flat_kernel = cl::Kernel(program, "flat_shader");
    
    // construct vectors
    rgbColor* C = new rgbColor[n];
    
    //put aside memory for buffer
    cl::Buffer buffer_A2(context, CL_MEM_READ_WRITE, sizeof(intersect) * numIntersections);
    cl::Buffer buffer_background(context, CL_MEM_READ_WRITE, sizeof(rgbColor));
    cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(rgbColor) * numIntersections);
    //write arrays to buffer
    queue.enqueueWriteBuffer(buffer_A2, CL_TRUE, 0, sizeof(intersect) * numIntersections, xsect);
    queue.enqueueWriteBuffer(buffer_background, CL_TRUE, 0, sizeof(rgbColor), &background);
    
    flat_kernel.setArg(0, buffer_A2);
    flat_kernel.setArg(1, buffer_background);
    flat_kernel.setArg(2, buffer_C2);
    
    queue.enqueueNDRangeKernel(flat_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(rgbColor)*n, C);
    
    queue.finish();
    
    delete[]kernel_code;
    
    return C;
}

rgbColor* cl_phong_shader(const intersect* xsect, phong& phongInfo, int numIntersections, rgbColor background, const cam camera, light* lights, int numLights)
{
    std::cout << "normals: " ;
//    for (int i = 0; i < numIntersections; ++i)
//    {
//        if (xsect[i].intersects)
//            printfloat3(xsect[i].normal);
//    }
    
    
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
    myfile.close();
    
    sources.push_back({kernel_code, filesize});
    
    cl::Program program(context, sources);
    if (program.build({default_device}) != CL_SUCCESS) {
        std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }
    
    // set up kernels and vectors for GPU code
    cl::CommandQueue queue(context, default_device);
    cl::Kernel persp_kernel = cl::Kernel(program, "phong_shader");
    
    // construct vectors
    rgbColor* C = new rgbColor[n];
    
    //put aside memory for buffer
    cl::Buffer buffer_intersects(context, CL_MEM_READ_WRITE, sizeof(intersect) * numIntersections);
    cl::Buffer buffer_result(context, CL_MEM_READ_WRITE, sizeof(rgbColor) * numIntersections);
    cl::Buffer buffer_lights(context, CL_MEM_READ_WRITE, sizeof(light) * numLights);
    
    //write arrays to buffer
    queue.enqueueWriteBuffer(buffer_intersects, CL_TRUE, 0, sizeof(intersect) * numIntersections, xsect);
    queue.enqueueWriteBuffer(buffer_lights, CL_TRUE, 0, sizeof(light) * numLights, lights);
    
    persp_kernel.setArg(0, phongInfo);
    persp_kernel.setArg(1, buffer_intersects);
    persp_kernel.setArg(2, background);
    persp_kernel.setArg(3, buffer_result);
    persp_kernel.setArg(4, camera);
    persp_kernel.setArg(5, buffer_lights);
    persp_kernel.setArg(6, numLights);
    
    queue.enqueueNDRangeKernel(persp_kernel, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
    queue.enqueueReadBuffer(buffer_result, CL_TRUE, 0, sizeof(rgbColor)*n, C);
    
    queue.finish();
    
    delete[]kernel_code;
    
    return C;
}
