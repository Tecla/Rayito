////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RAYITO_H__
#define __RAYITO_H__

#include "RMath.h"
#include "RRay.h"
#include "RMaterial.h"
#include "RScene.h"
#include "RLight.h"


namespace Rayito
{


//
// Image (collection of colored pixels with a width x height)
//

class Image
{
public:
    Image(size_t width, size_t height)
        : m_width(width), m_height(height), m_pixels(new Color[width * height]) { }
    
    virtual ~Image() { delete[] m_pixels; }
    
    size_t width()  const { return m_width; }
    size_t height() const { return m_height; }
    
    Color& pixel(size_t x, size_t y)
    {
        return m_pixels[y * m_width + x];
    }
    
protected:
    size_t m_width, m_height;
    Color *m_pixels;
};


//
// Cameras
//

class Camera
{
public:
    Camera() { }
    
    virtual ~Camera() { }
    
    // Generate a ray origin+direction for the camera, possibly with depth-of-field
    virtual Ray makeRay(float xScreen, float yScreen, float lensU, float lensV) const = 0;
};


class PerspectiveCamera : public Camera
{
public:
    // Create a perspective camera, with a given field-of-view in degrees,
    // look-at parameters, and depth-of-field parameters (lensRadius=0 to disable DOF)
    PerspectiveCamera(float fieldOfViewInDegrees,
                      const Point& origin,
                      const Vector& target,
                      const Vector& targetUpDirection,
                      float focalDistance,
                      float lensRadius);
    
    virtual ~PerspectiveCamera() { }
    
    virtual Ray makeRay(float xScreen, float yScreen, float lensU, float lensV) const;
    
protected:
    Point m_origin;
    Vector m_forward;
    Vector m_right;
    Vector m_up;
    float m_tanFov;
    float m_focalDistance;
    float m_lensRadius;
};


//
// Ray tracing
//

// Path trace through the scene, starting with an initial ray.
// Pass along scene information and various samplers so that we can reduce noise
// along the way.
Color pathTrace(const Ray& ray,
                ShapeSet& scene,
                std::list<Shape*>& lights,
                Rng& rng,
                size_t lightSamplesHint,
                size_t maxRayDepth,
                size_t pixelSampleIndex,
                Sampler** bounceSamplers);

// Generate a ray-traced image of the scene, with the given camera, resolution,
// and sample settings
Image* raytrace(ShapeSet& scene,
                const Camera& cam,
                size_t width,
                size_t height,
                size_t pixelSamplesHint,
                size_t lightSamplesHint,
                size_t maxRayDepth);


} // namespace Rayito


#endif // __RAYITO_H__
