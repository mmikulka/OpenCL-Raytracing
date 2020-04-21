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
    float x_resolution;
    float y_resolution;
    float left;
    float right;
    float top;
    float bottom;
}vp;

typedef struct __attribute__((packed))_viewRay
{
    float3 origin;
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
    float3 location;
    float3 normal;
    bool intersects;
    float t_;
}intersect;

__kernel void uv(__global float3* pos, __global vp * viewPort, __global float2* out)
{
 const int i = get_global_id(0);
    out[i].x = viewPort->left + ((viewPort->right - viewPort->left) * (pos[i].x + 0.5)) / viewPort->x_resolution;
    out[i].y = viewPort->bottom + ((viewPort->top - viewPort->bottom) * (pos[i].y + 0.5)) / viewPort->y_resolution;
}

__kernel void ortho_viewrays(cam c, __global float2* uv, __global viewRay * rays)
{
    const int i = get_global_id(0);
    rays[i].origin = c.eye;
    rays[i].origin.x = (c.u.x * uv[i].x) + (c.v.x * uv[i].y);
    rays[i].origin.y = (c.u.y * uv[i].x) + (c.v.y * uv[i].y);
    rays[i].origin.z = (c.u.z * uv[i].x) + (c.v.z * uv[i].y);
    rays[i].direction = -c.w;
}

static float determinant3x3(float3 system1, float3 system2, float3 system3)
{
    float det = system1.x * (system2.y * system3.z - system2.z * system3.y);
    det -= system1.y * (system2.x * system3.z - system2.z * system3.x);
    det += system1.z * (system2.x * system3.y - system2.y * system3.x);
    return det;
}

static float3 cl_cross ( float3 a, float3 b)
{
    float3 crossProduct;
    crossProduct.x = a.y * b.z - a.z * a.y;
    crossProduct.y = a.z * b.x - a.x * b.z;
    crossProduct.z = a.x * b.y - a.y * b.x;
    return crossProduct;
}

static float cl_dot (float3 a, float3 b)
{
    float product = a.x * b.x;
    product += a.y * b.y;
    product += a.z * b.z;
    return product;
}

static float cl_magnitude(float3 a)
{
    float squaredMag = dot(a, a);
    return sqrt(squaredMag);
}

static float3 cl_normalize (float3 a)
{
    return a / cl_magnitude(a);
}

static float3 solve(float3 system1, float3 system2, float3 system3, float3 solution)
{
    float3 linear_solution;
    float determinant = determinant3x3(system1, system2, system3);
    linear_solution.x = determinant3x3(solution, system2, system3) / determinant;
    linear_solution.y = determinant3x3(system1, solution, system3) / determinant;
    linear_solution.z = determinant3x3(system1, system2, solution) / determinant;
    return linear_solution;
}

static intersect tri_intersect(object triangle, viewRay ray, float t_upper_bound, float t_lower_bound)
{
    intersect intersection;
    intersection.t_ = -1;
    intersection.intersects = false;
    //setup variables to solve linear solution
    float3 system1 = triangle.a - triangle.b;
    float3 system2 = triangle.a - triangle.c;
    float3 system3 = ray.direction;
    float3 solution = triangle.a - ray.origin;
    
    float3 linearSolution = solve(system1, system2, system3, solution);
    
    
    //analize solution to make sure ray hits the triangle.
    if (linearSolution.z < t_lower_bound || linearSolution.z > t_upper_bound)
    {  return intersection;}
    if (linearSolution.y < 0 || linearSolution.y > 1)
    {  return intersection;}
    if (linearSolution.x < 0 || linearSolution.x > 1 - linearSolution.y)
    {  return intersection;}
    //ray hits the triangle, so calculate intersection, return the intersection object.
    intersection.obj.triObj = triangle;
    intersection.obj.is_triangle = true;
    intersection.obj.is_circle = false;
    intersection.location = ray.origin + ray.direction * linearSolution.z;
    float3 a_b = triangle.a - triangle.b;
    float3 a_c = triangle.a - triangle.c;
    intersection.normal = cl_cross(a_b, a_c);
    intersection.intersects = true;
    intersection.t_ = linearSolution.z;
    return intersection;
}

static intersect circle_intersect(circle circ, viewRay ray, float t_upper_bound, float t_lower_bound)
{
    intersect xsect;
    xsect.intersects = false;
    xsect.t_ = -1.0f;
    float3 origin_min_center = ray.origin - circ.center;
    float a = cl_dot(ray.direction, ray.direction);
    float b = cl_dot(ray.direction, origin_min_center);
    float radsqr = circ.radius * circ.radius;
    float c = cl_dot(origin_min_center, origin_min_center);
    c = c - radsqr;
    float descriminate = b * b - a * c;
    if (descriminate < 0) // no intersection
    {
//        return 1.0f;
       return xsect;
    }
    else
    {
        float t =  -b;
        if (descriminate > 0)//has 2 intersection
        {
          t -= sqrt(descriminate);
          t /= a;
          float3 location = ray.origin + ray.direction * t;
            float3 normal = cl_normalize((location - circ.center) / circ.radius);
          //test to make sure t is within the uper and lower bounds
          if(t < t_upper_bound && t > t_lower_bound)
          {
              xsect.obj.circleObj = circ;
              xsect.obj.is_circle = true;
              xsect.obj.is_triangle = false;
              xsect.location = location;
              xsect.normal = normal;
              xsect.intersects = true;
              xsect.t_ = t;
//              return 2.0f;
              return xsect;
          }
          //t is not withing upper and lower bounds, so the ray does not cross object.
//          else return 3.0f;
          else return xsect;
        }
        else // has only one intersecton
        {
          t /= a;
          float3 location = ray.origin + ray.direction * t;
            float3 normal = cl_normalize((location - circ.center) / circ.radius);
          //make sure t is within upper and lower bounds.
          if (t < t_upper_bound && t > t_lower_bound)
          {
              xsect.obj.circleObj = circ;
              xsect.obj.is_circle = true;
              xsect.obj.is_triangle = false;
              xsect.location = location;
              xsect.normal = normal;
              xsect.intersects = true;
              xsect.t_ = t;
//              return 4.0f;
              return xsect;;
          }
          //t is not within upper and lower bounds, so we act as though the ray does not cross with object.
//          else return 5.0f;
           else return xsect;
        }
    }
}

__kernel void intersections(__global object* objects, int numObjects, __global viewRay* rays, __global intersect* intersections, float t_upper_bound, float t_lower_bound, __global float3* debug)
{
    float temp_max = t_upper_bound;
    intersect closest;
    closest.t_ = -1;
    closest.intersects = false;
    closest.obj.is_triangle = false;
    closest.obj.is_circle = false;
    const int i = get_global_id(0);
    
    for (int j = 0; j < 1; ++j)
    {
        if (objects[j].is_triangle)
        {
            intersect temp = tri_intersect(objects[j].triObj, rays[i], temp_max, t_lower_bound);
            
            if (temp.intersects && temp.t_ > t_lower_bound && temp.t_ < temp_max)
            {
                    closest = temp;
                    temp_max = closest.t_;
            }
        }
        else
        {
            circle temp_circle = objects[i].circleObj;
            intersect temp = circle_intersect(temp_circle, rays[i], temp_max, t_lower_bound);
           
            if (temp.intersects && temp.t_ > t_lower_bound && temp.t_ < temp_max)
            {
                closest = temp;
                temp_max = closest.t_;
            }
        }
    }
    if (!closest.intersects)
    {
        closest.t_ = -1.0f;
    }
    intersections[i] = closest;
}

__kernel void flat_shader(__global intersect* xsects, __global rgbColor* background, __global rgbColor* result)
{
    int i = get_global_id(0);
    if (xsects[i].t_ < 0)
    {
        result[i].r = background->r;
        result[i].g = background->g;
        result[i].b = background->b;
    }
    else if(xsects[i].obj.is_circle)
    {
        result[i] = xsects[i].obj.circleObj.color;
    }
    else
    {
        result[i] = xsects[i].obj.triObj.color;
    }
}



//test function not neede by program
__kernel void parallel_add(__global float3* arr, __global float3* z){
 const int i = get_global_id(0);
    z[i].x = arr[i].x + arr[i].y;
    z[i].y = 0.0;
    z[i].z = 0.0;
}
