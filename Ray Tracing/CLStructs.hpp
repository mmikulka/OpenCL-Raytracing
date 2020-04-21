/*
 Structures that are needed since OpenCL does not deal with classes.
 
 With these you pass less variables into the CL programs.
 */
#pragma once

typedef struct __attribute__((packed))_vec3
{
    float x;
    float y;
    float z;
}vec3;

typedef struct __attribute__((packed))_UV
{
    float u;
    float v;
}UV;

typedef struct __attribute__((packed))_rgbColor
{
    float r;
    float g;
    float b;
    float a;
}rgbColor;

typedef struct __attribute__((packed))_light
{
    cl_float3 location;
    rgbColor color;
    float intensity;
    //extra dummy vars to make sure size of struct is multiple of float4
    float dummy1;
    float dummy2;
    float dummy3;
}light;

typedef struct __attribute__((packed))__cam
{
    cl_float3 eye;
    cl_float3 u;
    cl_float3 v;
    cl_float3 w;
}cam;

typedef struct __attribute__((packed))_vp
{
    float x_resolution;
    float y_resolution;
    float left;
    float right;
    float top;
    float bottom;
    //extra dummy vars to make sure size of struct is multiple of float4
    float dummy1;
    float dummy2;
}vp;

typedef struct __attribute__((packed))_viewRay
{
    cl_float3 origin;
    cl_float3 direction;
}viewRay;

typedef struct __attribute__((packed))_phong
{
    float ambient_coef;
    float diffuse_coef;
    float specular_coef;
    rgbColor ambient_color;
    //extra dummy vars to make sure size of struct is multiple of float4
    float dummy;
}phong;

typedef struct __attribute__((packed))_triangle
{
    rgbColor color;
    float shininess;
    cl_float3 a;
    cl_float3 b;
    cl_float3 c;
    //extra dummy vars to make sure size of struct is multiple of float4
    float dummy1;
    float dummy2;
    float dummy3;
}triangle;

typedef struct __attribute__((packed))_circle
{
    rgbColor color;
    float shininess;
    cl_float3 center;
    float radius;
    //extra dummy vars to make sure size of struct is multiple of float4
    float dummy1;
    float dummy2;
}circle;

typedef struct __attribute__((packed))_object
{
    triangle triObj;
    circle circleObj;
    bool is_triangle;
    bool is_circle;
    //extra dummy vars to make sure size of struct is multiple of float4
    bool dummy1;
    bool dummy2;
    float dummy3;
    float dummy4;
    float dummy5;
}object;

typedef struct __attribute__((packed))_intersect
{
    object obj;
    cl_float3 location;
    cl_float3 normal;
    bool intersects;
    float t_;
    //extra dummy vars to make sure size of struct is multiple of float4
    bool dummy1;
    bool dummy2;
    bool dummy3;
    float dummy4;
    float dummy5;
}intersect;
