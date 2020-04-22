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
        float3 location;
        rgbColor color;
        float intensity;
        //extra dummy vars to make sure size of struct is multiple of float4
        float dummy1;
        float dummy2;
        float dummy3;
    }light;

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
        //extra dummy vars to make sure size of struct is multiple of float4
        float dummy1;
        float dummy2;
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
        //extra dummy vars to make sure size of struct is multiple of float4
        float dummy;
    }phong;

typedef struct __attribute__((packed))_triangle
    {
        rgbColor color;
        float shininess;
        float3 a;
        float3 b;
        float3 c;
        //extra dummy vars to make sure size of struct is multiple of float4
        float dummy1;
        float dummy2;
        float dummy3;
    }triangle;

typedef struct __attribute__((packed))_circle
    {
        rgbColor color;
        float shininess;
        float3 center;
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
        float3 location;
        float3 normal;
        bool intersects;
        float t_;
        //extra dummy vars to make sure size of struct is multiple of float4
        bool dummy1;
        bool dummy2;
        bool dummy3;
        float dummy4;
        float dummy5;
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

__kernel void persp_viewrays(cam c, __global float2* uv, __global viewRay * rays, float focal_length)
{
    const int i = get_global_id(0);
    rays[i].origin = c.eye;
    rays[i].direction = -c.w * focal_length;
    rays[i].direction += c.u * uv[i].x;
    rays[i].direction += c.v * uv[i].y;
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
    crossProduct.x = a.y * b.z - a.z * b.y;
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
    float squaredMag = a.x * a.x;
    squaredMag += a.y * a.y;
    squaredMag += a.z * a.z;
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

static intersect tri_intersect(triangle triangle, viewRay ray, float t_upper_bound, float t_lower_bound)
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
            float3 normal = cl_normalize((location - circ.center));
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

__kernel void intersections(__global object* objects, int numObjects, __global viewRay* rays, __global intersect* intersections, float t_upper_bound, float t_lower_bound, __global float* debug)
{
    float temp_max = t_upper_bound;
    intersect closest;
    closest.t_ = -1;
    closest.intersects = false;
    closest.obj.is_triangle = false;
    closest.obj.is_circle = false;
    const int i = get_global_id(0);
    
    for (int j = 0; j < numObjects; ++j)
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
            intersect temp = circle_intersect(objects[j].circleObj, rays[i], temp_max, t_lower_bound);
            
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

__kernel void phong_shader(phong phongInfo, __global intersect* xsects, rgbColor background, __global rgbColor* color, cam camera, __global light * lights, int numLights)
{
    int i = get_global_id(0);
    if (xsects[i].intersects)
    {
    color[i].r = (phongInfo.ambient_coef * phongInfo.ambient_color.r);
    color[i].g = (phongInfo.ambient_coef* phongInfo.ambient_color.g);
    color[i].b = (phongInfo.ambient_coef * phongInfo.ambient_color.b);


    for (int j = 0; j < numLights; ++j)
    {
        //compute light vector from object.
        float3 lightDirection = cl_normalize(lights[j].location - xsects[i].location);
        float3 viewDirection = cl_normalize(camera.eye - xsects[i].location);
        //computre bisector
        float3 bisector = cl_normalize(viewDirection + lightDirection);
        //compute normal dot light location and choose either that or 0
        float ndl = cl_dot(xsects[i].normal, lightDirection);
        if (ndl < 0) ndl = 0;
        //set variable to the object's color to make code cleaner
        rgbColor objColor;
        if (xsects[i].obj.is_triangle)
        {
            objColor = xsects[i].obj.triObj.color;
        }
        else
        {
            objColor = xsects[i].obj.circleObj.color;
        }
        // add diffuse light output to intensity of pixel.
        color[i].r += objColor.r * phongInfo.diffuse_coef * ndl;
        color[i].g += objColor.g * phongInfo.diffuse_coef * ndl;
        color[i].b += objColor.b * phongInfo.diffuse_coef * ndl;
        // compute normal and bisector and choose 0 or that.
        float ndh =  dot(bisector, xsects[i].normal);
        if (ndh < 0) ndh = 0;
        // raise ndh to phong exponent
        if (xsects[i].obj.is_triangle)
        {
            ndh = powr(ndh, xsects[i].obj.triObj.shininess);
        }
        else
        {
            ndh = powr(ndh, xsects[i].obj.circleObj.shininess);
        }        //add specular coefficient to intensity of pixel
        color[i].r += phongInfo.specular_coef * ndh * lights[j].color.r;
        color[i].g += phongInfo.specular_coef * ndh * lights[j].color.g;
        color[i].b += phongInfo.specular_coef * ndh * lights[j].color.b;
    }
    //loop though each color and make sure it is not less than 0 or greatre than 1

    color[i].r = (color[i].r < 0)? 0 : color[i].r;
    color[i].r = (color[i].r > 1)? 1 : color[i].r;
    color[i].g = (color[i].g < 0)? 0 : color[i].g;
    color[i].g = (color[i].g > 1)? 1 : color[i].g;
    color[i].b = (color[i].b < 0)? 0 : color[i].b;
    color[i].b = (color[i].b > 1)? 1 : color[i].b;
    }
    else
    {
        color[i] = background;
    }
}
