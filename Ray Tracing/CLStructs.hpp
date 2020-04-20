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
}rgbColor;

typedef struct __attribute__((packed))_light
{
    cl_float3 location;
    rgbColor color;
    float intensity;
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
}phong;

typedef struct __attribute__((packed))_triangle
{
    rgbColor color;
    float shininess;
    cl_float3 a;
    cl_float3 b;
    cl_float3 c;
}triangle;

typedef struct __attribute__((packed))_circle
{
    rgbColor color;
    float shininess;
    cl_float3 center;
    float radius;
}circle;

typedef struct __attribute__((packed))_object
{
    triangle triObj;
    circle circleObj;
    bool is_triangle;
    bool is_circle;
}object;

typedef struct __attribute__((packed))_intersect
{
    object obj;
    cl_float3 location;
    cl_float3 normal;
    float t_;
}intersect;
