
//
//  gfxraytrace.hpp
//  raytracer
//
//  Created by Theatre Floater on 4/5/20.
//  Copyright © 2020 Theatre Floater. All rights reserved.
//

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "gfxalgebra.hpp"
#include "gfximage.hpp"
#include "rayson.hpp"
#include "CLHelpers.hpp"

namespace gfx {

// Forward declarations of new class types, in alphabetical order.

class abstract_projection;
class abstract_scene_object;
class abstract_shader;
class blinn_phong_shader;
class camera;
class flat_shader;
class intersection;
class orthographic_projection;
class perspective_projection;
class point_light;
class raytracer;
class scene;
class scene_read_exception;
class scene_mesh;
class scene_sphere;
class scene_triangle;
class view_ray;
class viewport;


// Class declarations, in the order that classes are introduced in
// the comment above.

// A scene represents all the geometric information necessary to
// render an image. That includes:
//
// - a background color, used to fill a pixel whose view ray does
//   not intersect any scene object;
//
// - a vector of point lights; and
//
// - a vector of scene objects.
//
// Ordinarily you need at least one light, and many scene objects,
// to make an interesting image. However this is not enforced with
// assertions.
class scene {
public:
    using light_storage_type = std::vector<std::unique_ptr<point_light>>;
    using object_storage_type = std::vector<std::unique_ptr<abstract_scene_object>>;
    
private:
    std::unique_ptr<camera> camera_;
    std::unique_ptr<viewport> viewport_;
    std::unique_ptr<abstract_projection> projection_;
    std::unique_ptr<abstract_shader> shader_;
    hdr_rgb background_;
    light_storage_type lights_;
    object_storage_type objects_;
    
public:
    
    scene(scene&&) noexcept = default;
    
    scene() noexcept
    : background_(BLACK) {
        assert(!complete());
    }
    
    // Constructor.
    scene(std::unique_ptr<camera> camera,
          std::unique_ptr<viewport> viewport,
          std::unique_ptr<abstract_projection> projection,
          std::unique_ptr<abstract_shader> shader,
          const hdr_rgb& background) noexcept
    : camera_(std::move(camera)),
    viewport_(std::move(viewport)),
    projection_(std::move(projection)),
    shader_(std::move(shader)),
    background_(background) {
        
        assert(camera_);
        assert(viewport_);
        assert(projection_);
        assert(shader_);
        
        assert(complete());
    }
    
    constexpr bool complete() const noexcept {
        return (camera_ && viewport_ && projection_ && shader_);
    }
    
    // Accessors.
    const camera& camera() const noexcept {
        assert(camera_);
        return *camera_;
    }
    const viewport& viewport() const noexcept {
        assert(viewport_);
        return *viewport_;
    }
    const abstract_projection& projection() const noexcept {
        assert(projection_);
        return *projection_;
    }
    const abstract_shader& shader() const noexcept {
        assert(shader_);
        return *shader_;
    }
    constexpr const hdr_rgb& background() const noexcept {
        return background_;
    }
    constexpr const light_storage_type& lights() const noexcept {
        return lights_;
    }
    constexpr const object_storage_type& objects() const noexcept {
        return objects_;
    }
    
    // Mutators.
    void camera(std::unique_ptr<gfx::camera> camera) noexcept {
        camera_ = std::move(camera);
    }
    void viewport(std::unique_ptr<gfx::viewport> viewport) noexcept {
        viewport_ = std::move(viewport);
    }
    void projection(std::unique_ptr<abstract_projection> projection) noexcept {
        projection_ = std::move(projection);
    }
    void shader(std::unique_ptr<abstract_shader> shader) noexcept {
        shader_ = std::move(shader);
    }
    void background(const hdr_rgb& background) noexcept {
        background_ = background;
    }
    
    // Adding lights and objects.
    void add_light(std::unique_ptr<point_light> light) noexcept {
        lights_.emplace_back(std::move(light));
    }
    void add_object(std::unique_ptr<abstract_scene_object> object) noexcept {
        objects_.emplace_back(std::move(object));
    }
    
    // Trace a ray and find the closest intersecting scene object.
    //
    // If no such intersection exists, return an empty optional.
    //
    // If there is an intersection within that t range, return an optional that
    // contains that intersection object.
    _intersect * intersect(const viewRay* rays, int numRays) const noexcept;
    
    // Render the given scene and return the resulting image.
    hdr_image render() const noexcept;
    
    // Load a scene from a JSON file; throws scene_read_exception on error.
    static scene read_json(const std::string& path) noexcept(false);
    
};

// An error encountered while trying to read and parse a scene file.
class scene_read_exception {
private:
    std::string path_;
    std::string message_;
    
public:
    scene_read_exception(const std::string& path,
                         const std::string& message) noexcept
    : path_(path), message_(message) { }
    
    const std::string& path() const noexcept { return path_; }
    const std::string& message() const noexcept { return message_; }
};

// The location and orientation of the camera. This is defined by a
// 3D point for the eye location; and a basis defined by three 3D
// normalized vectors.
class camera {
private:
    vector3<double> eye_, u_, v_, w_;
    
public:
    
    // Constructor that provides the eye location and three basis
    // vectors directly. u, v, and w must each be normalized (magnitude
    // 1.0).
    constexpr camera(const vector3<double>& eye,
                     const vector3<double>& u,
                     const vector3<double>& v,
                     const vector3<double>& w) noexcept
    : eye_(eye), u_(u), v_(v), w_(w) {
        
        assert(approx_equal(u.magnitude(), 1.0, .01));
        assert(approx_equal(v.magnitude(), 1.0, .01));
        assert(approx_equal(w.magnitude(), 1.0, .01));
    }
    
    // Constructor that computes the basis in terms of a given
    // view-direction and up vector.
    constexpr camera(const vector3<double>& eye,
                     const vector3<double>& view_direction,
                     const vector3<double>& up) noexcept;
    
    // Accessors and mutators.
    constexpr const vector3<double>& eye() const noexcept {
        return eye_;
    }
    constexpr const vector3<double>& u() const noexcept {
        return u_;
    }
    constexpr const vector3<double>& v() const noexcept {
        return v_;
    }
    constexpr const vector3<double>& w() const noexcept {
        return w_;
    }
    
    cam to_struct() const noexcept;
};

// A viewport defines the boundary of the viewing window. It stores
// the width and height of the image, in screen coordinates; and the
// left, right, top, and bottom of the view window in world
// coordinates.
class viewport {
private:
    size_t x_resolution_, y_resolution_;
    double left_, top_, right_, bottom_;
    
public:
    
    // Constructor. The following inequalities must hold:
    //
    // x_resolution, y_resolution > 0
    // left < 0 < right
    // bottom < 0 < top
    //
    constexpr viewport(size_t x_resolution,
                       size_t y_resolution,
                       double left,
                       double top,
                       double right,
                       double bottom) noexcept
    : x_resolution_(x_resolution),
    y_resolution_(y_resolution),
    left_(left),
    top_(top),
    right_(right),
    bottom_(bottom) {
        assert(x_resolution > 0);
        assert(y_resolution > 0);
        assert(left < right);
        assert(bottom < top);
        assert(left < 0.0);
        assert(right > 0.0);
        assert(top > 0.0);
        assert(bottom < 0.0);
    }
    
    // Accessors.
    constexpr size_t x_resolution() const noexcept {
        return x_resolution_;
    }
    constexpr size_t y_resolution() const noexcept {
        return y_resolution_;
    }
    constexpr double left() const noexcept {
        return left_;
    }
    constexpr double right() const noexcept {
        return right_;
    }
    constexpr double top() const noexcept {
        return top_;
    }
    constexpr double bottom() const noexcept {
        return bottom_;
    }
    
    // Map an (x, y) screen coordinate to a (u, v) coordinate in
    // [0, 1]^2. Return a 2D vector x, with u in x[0] and v in x[1];
    vector2<double> uv(size_t x, size_t y) const noexcept;
    
    vp to_struct()const noexcept;
};

// Abstract class defining a projection algorithm.
class abstract_projection {
public:
    
    // Given a camera and (u, v) coordinate within the viewport,
    // create and return the corresponding viewing ray. (u, v) are
    // expected to come from the camera::uv function.
    virtual viewRay* compute_view_ray(cam& c,
                                      cl_float2* uv, int numPixels) const noexcept = 0;
    
    virtual ~abstract_projection() noexcept = default;
};

// Orthographic implementation of abstract_projection.
class orthographic_projection : public abstract_projection {
public:
    
    constexpr orthographic_projection() noexcept = default;
    
    virtual viewRay* compute_view_ray(cam& c,
                                      cl_float2* uv, int numPixels) const noexcept;
};

// Perspective implementation of abstract_projection.
class perspective_projection : public abstract_projection {
private:
    double focal_length_;
    
public:
    
    // The perspective projection algorithm needs to know the
    // focal_length of the camera, which is the distance between the
    // eye and the view plane, and must be positive.
    constexpr perspective_projection(double focal_length) noexcept
    : focal_length_(focal_length) {
        assert(focal_length > 0.0);
    }
    
    // Accessor.
    constexpr double focal_length() const noexcept {
        return focal_length_;
    }
    
    virtual viewRay* compute_view_ray(cam& c,
                                      cl_float2* uv, int numPixels) const noexcept;
};

// Abstract class defining a shading algorithm.
class abstract_shader {
public:
    
    // Given a scene, camera, and particular ray-object intersection,
    // compute the color of the pixel corresponding to the view
    // ray. The pixel's color is returned.
    virtual hdr_rgb* shade(const scene& scene,
                          const cam& camera,
                          const _intersect* xsect,
                          int numIntersections) const noexcept = 0;
    
    virtual ~abstract_shader() noexcept = default;
};

// Flat-shader implementation of abstract_shader.
class flat_shader : public abstract_shader {
public:
    
    virtual hdr_rgb* shade(const scene& scene,
                          const cam& camera,
                          const _intersect* xsect,
                          int numIntersections) const noexcept;
};

// Blin-Phong implementation of abstract_shader.
class blinn_phong_shader : public abstract_shader {
private:
    double ambient_coefficient_;
    hdr_rgb ambient_color_;
    double diffuse_coefficient_, specular_coefficient_;
    
public:
    
    // The Blinn-Phong model depends on the following parameters:
    //
    // ambient_coefficient is the multiplier for ambient light, which
    // must be non-negative. When zero there will be no ambient light.
    //
    // ambient_color is the color of ambient light; usually white in
    // daylight.
    //
    // diffuse_coefficient is the multiplier for diffuse light (object
    // color), which must be non-negative. When zero there is no
    // diffuse light, so only ambient and specular light would be
    // visible.
    //
    // specular_coefficient is the multiplier for specular light
    // (speckles/gloss/glare), which must be non-negative. When zero
    // there are no specular highlights so all objects appear matte.
    //
    blinn_phong_shader(double ambient_coefficient,
                       const hdr_rgb& ambient_color,
                       double diffuse_coefficient,
                       double specular_coefficient)
    : ambient_coefficient_(ambient_coefficient),
    ambient_color_(ambient_color),
    diffuse_coefficient_(diffuse_coefficient),
    specular_coefficient_(specular_coefficient) {
        
        assert(ambient_coefficient >= 0.0);
        assert(diffuse_coefficient >= 0.0);
        assert(specular_coefficient >= 0.0);
    }
    
    // Accessors.
    constexpr double ambient_coefficient () const noexcept {
        return ambient_coefficient_ ;
    }
    constexpr const hdr_rgb& ambient_color() const noexcept {
        return ambient_color_;
    }
    constexpr double diffuse_coefficient () const noexcept {
        return diffuse_coefficient_;
    }
    constexpr double specular_coefficient() const noexcept {
        return specular_coefficient_;
    }
    
    virtual hdr_rgb* shade(const scene& scene,
                          const cam& camera,
                          const _intersect* xsect,
                          int numIntersections) const noexcept;
    
    phong to_struct() const noexcept;
};

// A view ray represents a ray traveling from the viewer out into
// the scene. It is defined by an origin, and direction, each of
// which is a 3D vector.
class view_ray {
private:
    vector3<double> origin_, direction_;
    
public:
    
    // Constructor with an explicit origin and direction. Direction
    // must be normalized (magnitude 1).
    constexpr view_ray(const vector3<double>& origin,
                       const vector3<double>& direction) noexcept
    : origin_(origin),
    direction_(direction) { }
    
    // Accessors.
    constexpr const vector3<double>& origin() const noexcept {
        return origin_;
    }
    constexpr const vector3<double>& direction() const noexcept {
        return direction_;
    }
};

// Abstract class for some kind of scene object.
class abstract_scene_object {
private:
    hdr_rgb color_;
    double shininess_;
    
public:
    
    // Construct an object with the given diffuse color and shininesss
    // value (Phong exponent). shininess must be positive.
    constexpr abstract_scene_object(const hdr_rgb& color,
                                    double shininess) noexcept
    : color_(color),
    shininess_(shininess) {
        
        assert(shininess > 0.0);
    }
    
    virtual ~abstract_scene_object() noexcept = default;
    
    // Accessors.
    constexpr const hdr_rgb& color() const noexcept {
        return color_;
    }
    constexpr double shininess() const noexcept {
        return shininess_;
    }
    
    virtual object to_struct() const noexcept = 0;
};

// A scene object that is a 3D sphere.
class scene_sphere : public abstract_scene_object {
private:
    vector3<double> center_;
    double radius_;
    
public:
    
    // Create a sphere with the given color, shininess, center
    // location, and radius. radius must be positive.
    constexpr scene_sphere(const hdr_rgb& color,
                           double shininess,
                           const vector3<double>& center,
                           double radius) noexcept
    : abstract_scene_object(color, shininess),
    center_(center),
    radius_(radius) {
        
        assert(radius > 0.0);
    }
    
    // Accessors.
    constexpr const vector3<double>& center() const noexcept {
        return center_;
    }
    constexpr double radius() const noexcept {
        return radius_;
    }
    virtual object to_struct() const noexcept;
};

class scene_triangle : public abstract_scene_object {
private:
    vector3<double> a_, b_, c_;
    
public:
    
    // The three vertices of the triangle are called a, b, c. Each is
    // a 3D location.
    constexpr scene_triangle(const hdr_rgb& color,
                             double shininess,
                             const vector3<double>& a,
                             const vector3<double>& b,
                             const vector3<double>& c) noexcept
    : abstract_scene_object(color, shininess),
    a_(a),
    b_(b),
    c_(c) { }
    
    // Accessors.
    constexpr const vector3<double>& a() const noexcept {
        return a_;
    }
    constexpr const vector3<double>& b() const noexcept {
        return b_;
    }
    constexpr const vector3<double>& c() const noexcept {
        return c_;
    }
    
    virtual object to_struct() const noexcept;
};

// A point_light represents a light source that gives off the same
// amount of light in all directions. The sun, or an idealized light
// bulb, can be modeled as a point light.
class point_light {
private:
    vector3<double> location_;
    hdr_rgb color_;
    double intensity_;
    
public:
    
    // Construct a point light at the given location, that emits light
    // of the given color, with the given scalar intensity. Intensity must
    // be positive.
    constexpr point_light(const vector3<double>& location,
                          const hdr_rgb& color,
                          double intensity) noexcept
    : location_(location),
    color_(color),
    intensity_(intensity) {
        
        assert(intensity > 0.0);
    }
    
    // Accessors.
    constexpr const vector3<double>& location () const noexcept {
        return location_;
    }
    constexpr const hdr_rgb& color() const noexcept {
        return color_;
    }
    constexpr double intensity() const noexcept {
        return intensity_;
    }
    
    light to_struct() const noexcept;
};

// An intersection represents a place where a view ray hits a
// scene object. It is defined by:
//
// - a non-owning pointer to the object that was hit;
//
// - the 3D point where the hit occurs;
//
// - a normal vector, that is perpendicular to the object at the hit
//   location; and
//
// - the t value where the hit happened relative to the view ray
//   direction, i.e.
//       location == ray.origin + (t * ray.direction)
//class intersection {
//private:
//    const abstract_scene_object *object_; // non-owning pointer
//    vector3<double> location_, normal_;
//    double t_;
//
//public:
//
//    // Construct an intersection.
//    // The object pointer must not be nullptr.
//    // The normal must be normalized (magnitude 1).
//    constexpr intersection(const abstract_scene_object* object,
//                           const vector3<double>& location,
//                           const vector3<double>& normal,
//                           double t) noexcept
//    : object_(object),
//    location_(location),
//    normal_(normal),
//    t_(t) {
//
//        assert(object != nullptr);
//        assert(approx_equal(normal.magnitude(), 1.0, .01));
//    }
//
//    constexpr const abstract_scene_object& object() const noexcept {
//        return *object_;
//    }
//    constexpr const vector3<double>& location() const noexcept {
//        return location_;
//    }
//    constexpr const vector3<double>& normal() const noexcept {
//        return normal_;
//    }
//    constexpr double t() const noexcept {
//        return t_;
//    }
//};

scene scene::read_json(const std::string& path) noexcept(false) {
    
    auto import_vector = [](const rayson::vector3& v) noexcept {
        return vector3<double>{v.x(), v.y(), v.z()};
    };
    
    auto import_color = [](const rayson::color& c) noexcept {
        return hdr_rgb{float(c.r()), float(c.g()), float(c.b())};
    };
    
    try {
        
        rayson::scene loaded = rayson::read_file(path);
        
        std::unique_ptr<abstract_projection> the_projection;
        if (loaded.projection().type == rayson::Projection::ORTHO) {
            the_projection = std::make_unique<orthographic_projection>();
        } else if (loaded.projection().type == rayson::Projection::PERSP) {
            auto persp = &loaded.projection().perspProjection;
            the_projection = std::make_unique<perspective_projection>(persp->focal_length());
        } else {
            assert(false); // unknown projection type, should be unreachable
        }
        
        std::unique_ptr<abstract_shader> the_shader;
        if (loaded.shader().type == rayson::Shader::FLAT) {
            the_shader = std::make_unique<flat_shader>();
        } else if (loaded.shader().type == rayson::Shader::PHONG) {
            auto phong = &loaded.shader().phongShader;
            the_shader = std::make_unique<blinn_phong_shader>(phong->ambient_coeff(),
                                                              import_color(phong->ambient_color()),
                                                              phong->diffuse_coeff(),
                                                              phong->specular_coeff());
        } else {
            assert(false); // unknown shader type, should be unreachable
        }
        
        scene result(std::make_unique<::gfx::camera>(import_vector(loaded.camera().eye()),
                                                     import_vector(loaded.camera().view()),
                                                     import_vector(loaded.camera().up())),
                     std::make_unique<::gfx::viewport>(loaded.viewport().x_resolution(),
                                                       loaded.viewport().y_resolution(),
                                                       loaded.viewport().left(),
                                                       loaded.viewport().top(),
                                                       loaded.viewport().right(),
                                                       loaded.viewport().bottom()),
                     std::move(the_projection),
                     std::move(the_shader),
                     import_color(loaded.background()));
        
        for (auto& light : loaded.point_lights()) {
            result.add_light(std::make_unique<point_light>(import_vector(light.location()),
                                                           import_color(light.color()),
                                                           light.intensity()));
        }
        
        for (auto& sphere : loaded.spheres()) {
            result.add_object(std::make_unique<scene_sphere>(import_color(sphere.material().color()),
                                                             sphere.material().shininess(),
                                                             import_vector(sphere.center()),
                                                             sphere.radius()));
        }
        
        for (auto& tri : loaded.triangles()) {
            result.add_object(std::make_unique<scene_triangle>(import_color(tri.material().color()),
                                                               tri.material().shininess(),
                                                               import_vector(tri.a()),
                                                               import_vector(tri.b()),
                                                               import_vector(tri.c())));
        }
        
        return result;
        
    } catch (rayson::read_exception e) {
        throw scene_read_exception(path, e.message());
    }
}

//implimentation

//conversion functions to clean up rest of code

cam camera::to_struct() const noexcept
{
    cam cl_camera = {{{static_cast<float>(eye_[0]), static_cast<float>(eye_[1]), static_cast<float>(eye_[2])}},
        {{static_cast<float>(u_[0]), static_cast<float>(u_[1]), static_cast<float>(u_[2])}},
        {{static_cast<float>(v_[0]), static_cast<float>(v_[1]), static_cast<float>(v_[2])}},
        {{static_cast<float>(w_[0]), static_cast<float>(w_[1]), static_cast<float>(w_[2])}}};
    return cl_camera;
}

vp viewport::to_struct() const noexcept
{
    vp viewPort = {static_cast<float>(x_resolution()), static_cast<float>(y_resolution()), static_cast<float>(left()), static_cast<float>(right()), static_cast<float>(top()), static_cast<float>(bottom())};
    return viewPort;
}

object scene_triangle::to_struct() const noexcept
{
    object tri;
    tri.triObj.color.s[0] = color().r();
    tri.triObj.color.s[1] = color().g();
    tri.triObj.color.s[2] = color().b();
    tri.triObj.shininess = static_cast<float>(shininess());
    tri.triObj.a.s[0] = static_cast<float>(a_[0]);
    tri.triObj.a.s[1] = static_cast<float>(a_[1]);
    tri.triObj.a.s[2] = static_cast<float>(a_[2]);
    tri.triObj.b.s[0] = static_cast<float>(b_[0]);
    tri.triObj.b.s[1] = static_cast<float>(b_[1]);
    tri.triObj.b.s[2] = static_cast<float>(b_[2]);
    tri.triObj.c.s[0] = static_cast<float>(c_[0]);
    tri.triObj.c.s[1] = static_cast<float>(c_[1]);
    tri.triObj.c.s[2] = static_cast<float>(c_[2]);
    
    tri.is_triangle = true;
    tri.is_circle = false;
    return tri;
}

object scene_sphere::to_struct() const noexcept
{
    object circle;
    circle.circleObj.color.s[0] = color().r();
    circle.circleObj.color.s[1] = color().g();
    circle.circleObj.color.s[2] = color().b();
    circle.circleObj.shininess = static_cast<float>(shininess());
    circle.circleObj.center.s[0] = static_cast<float>(center_[0]);
    circle.circleObj.center.s[1] = static_cast<float>(center_[1]);
    circle.circleObj.center.s[2] = static_cast<float>(center_[2]);
    circle.circleObj.radius = static_cast<float>(radius_);
    
    circle.is_triangle = false;
    circle.is_circle = true;
    return circle;
}

light point_light::to_struct() const noexcept
{
    _light light;
    
    light.location.s[0] = location_[0];
    light.location.s[1] = location_[1];
    light.location.s[2] = location_[2];
    
    light.color.s[0] = color_.r();
    light.color.s[1] = color_.g();
    light.color.s[2] = color_.b();
    
    light.intensity = intensity_;
    
    return light;
}

phong blinn_phong_shader::to_struct() const noexcept
{
    phong shaderInfo;
    shaderInfo.ambient_coef = ambient_coefficient_;
    shaderInfo.diffuse_coef = diffuse_coefficient_;
    shaderInfo.specular_coef = specular_coefficient_;
    shaderInfo.ambient_color.s[0] = ambient_color_.r();
    shaderInfo.ambient_color.s[1] = ambient_color_.g();
    shaderInfo.ambient_color.s[2] = ambient_color_.b();
    return shaderInfo;
}
    
_intersect* scene::intersect (const viewRay* rays, int numRays) const noexcept
    {
    float tmin = 0;
    float tmax = std::numeric_limits<float>::max();
    
    _intersect * xsects = new _intersect[numRays];
    
    int numObjects = objects_.size();
    
    object objects[numObjects];
    
    for (int i = 0; i < numObjects; ++i)
    {
        objects[i] = objects_[i]->to_struct();
    }
    xsects = cl_intersect(objects, numObjects, rays, numRays, tmax, tmin);
    
    return xsects;
}

hdr_image scene::render() const noexcept {
    
    assert(camera_);
    assert(viewport_);
    assert(projection_);
    assert(shader_);
    
    assert(viewport_->x_resolution() > 0);
    assert(viewport_->y_resolution() > 0);
    
    size_t w = viewport_->x_resolution(),
    h = viewport_->y_resolution();
    
    
    hdr_image result(w, h, background_);
    assert(!result.is_empty());
    
    int numPixels = h * w;
    cl_float3 pixels[h * w];
    
    for (size_t i = 0; i < numPixels; ++i) {
        pixels[i] = {{static_cast<float>(i % w) + 1.0f, static_cast<float>(i / w) + 1.0f, 0.0f}};
    }
    
    
    vp cl_viewport = viewport_->to_struct();
    //compute uv
    cl_float2 * uV = cl_uv(pixels, cl_viewport, numPixels);
    
    //convert gfx::camera to a struct.
    cam cl_camera = camera_->to_struct();
    
   viewRay* rays = projection_->compute_view_ray(cl_camera, uV, numPixels);
    
    // if ray hits object then evaluate shading model
    _intersect* xsects = intersect(rays, numPixels);
    
    //shade each pixel
    hdr_rgb* colors;
    colors = shader_->shade(*this, cl_camera, xsects, numPixels);
    
    //transfer colors to hdr_image variable.
    for (int i = 0; i < numPixels; ++i)
    {
        result.pixel(i % w, i / w, colors[i]);
    }
    
    return result;
}


constexpr camera::camera(const vector3<double>& eye,
                         const vector3<double>& view_direction,
                         const vector3<double>& up) noexcept {
    
    eye_ = eye;
    w_ = -(view_direction / view_direction.magnitude());
    u_ = (up.cross(w_)) / (up.cross(w_).magnitude());
    v_ = w_.cross(u_);
}

viewRay* orthographic_projection::compute_view_ray(cam& c,
                                                   cl_float2 * uv, int numPixels) const noexcept {
    
    viewRay* rays = cl_ortho_viewrays(c, uv, numPixels);

    return rays;
}

viewRay*  perspective_projection::compute_view_ray(cam& c,
                                                  cl_float2 * uv, int numPixels) const noexcept {
    viewRay* rays = cl_persp_viewrays(c, uv, numPixels, static_cast<float>(focal_length_));
    return rays;
}

hdr_rgb* flat_shader::shade(const scene& scene,
                           const cam& camera,
                           const _intersect* xsect, int numIntersections) const noexcept {
    
    rgbColor background;
    background.r = scene.background().r();
    background.g = scene.background().g();
    background.b = scene.background().b();
    
    rgbColor* color = cl_flat_shader(xsect, numIntersections, background);

    hdr_rgb* hdrColors = new hdr_rgb[numIntersections];
    for (int i = 0; i < numIntersections; ++i)
    {
        hdrColors[i] = hdr_rgb(color[i].r, color[i].g, color[i].b);
    }

    
    return hdrColors;
}

hdr_rgb* blinn_phong_shader::shade(const scene& scene,
                                  const cam& camera,
                                  const _intersect* xsect, int numIntersections) const noexcept {
    
    //convert light objects to struct form
    const int numLights = scene.lights().size();
    light lights[numLights];
    for(int i = 0; i < numLights; ++i)
    {
        lights[i] = scene.lights()[i]->to_struct();
    }
    phong phongInfo = to_struct();
    
    rgbColor background;
    background.r = scene.background().r();
    background.g = scene.background().g();
    background.b = scene.background().b();
    
    rgbColor* color = cl_phong_shader(xsect, phongInfo, numIntersections, background, camera, lights, numLights);
    
    hdr_rgb* hdrColors = new hdr_rgb[numIntersections];
    for (int i = 0; i < numIntersections; ++i)
    {
        hdrColors[i] = hdr_rgb(color[i].r, color[i].g, color[i].b);
    }
    return hdrColors;
}


}
