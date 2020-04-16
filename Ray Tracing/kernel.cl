typedef struct __attribute__((packed))_rgbColor
{
    float r;
    float g;
    float b;
}rgbColor;

typedef struct __attribute__((packed))_light
{
    float3 location;
    rgbColor color;
    float intensity;
}lght;

typedef struct __attribute__((packed))__cam
{
    float3 eye;
    float3 u;
    float3 v;
    float3 w;
}cam;

typedef struct __attribute__((packed))_vp
{
    size_t x_resolution;
    size_t y_resolution;
    float left;
    float right;
    float top;
    float bottom;
}vp;

typedef struct __attribute__((packed))_viewRay
{
    float3 orgin;
    float3 direction;
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
    float3 a;
    float3 b;
    float3 c;
}triangle;

typedef struct __attribute__((packed))_circle
{
    rgbColor color;
    float shininess;
    float3 center;
    float radius;
}circle;

typedef struct __attribute__((packed))_intersect
{
    triangle *triangle_obj;
    circle *circle_obj;
    float3 location;
    float3 normal;
    float t_;
}intersect;

__kernel void uv(__global float3* pos, __global vp * viewPort, __global float3* out)
{
 const int i = get_global_id(0);
    out[i].x = viewPort[i].left + ((viewPort[i].right - viewPort[i].left) * (pos[i].x + 0.5)) / viewPort[i].x_resolution;
    out[i].y = viewPort[i].bottom + ((viewPort[i].top - viewPort[i].bottom) * (pos[i].y + 0.5)) / viewPort[i].y_resolution;
    out[i].z = 0.0;
}

__kernel void parallel_add(__global float3* arr, __global float3* z){
 const int i = get_global_id(0);
    z[i].x = arr[i].x + arr[i].y;
    z[i].y = 0.0;
    z[i].z = 0.0;
}
