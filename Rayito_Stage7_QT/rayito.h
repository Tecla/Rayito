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
    Camera(float shutterOpen = 0.0f, float shutterClose = 0.0f)
        : m_shutterOpen(shutterOpen), m_shutterClose(shutterClose) { }
    
    virtual ~Camera() { }
    
    // Generate a ray origin+direction for the camera, possibly with depth-of-field
    virtual Ray makeRay(float xScreen, float yScreen, float lensU, float lensV, float timeU) const = 0;
    
protected:
    float m_shutterOpen;
    float m_shutterClose;
    
    float time(float timeU) const { return m_shutterOpen + (m_shutterClose - m_shutterOpen) * timeU; }
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
                      float lensRadius,
                      float shutterOpen,
                      float shutterClose);
    
    virtual ~PerspectiveCamera() { }
    
    virtual Ray makeRay(float xScreen, float yScreen, float lensU, float lensV, float timeU) const;
    
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
// Sampler container (for a given pixel, holds the samplers for all random features and bounces)
//

struct SamplerContainer
{
    // These are sampled once per pixel sample to start a path
    Sampler* m_lensSampler;
    Sampler* m_subpixelSampler;
    Sampler* m_timeSampler;
    // This is sampled once per bounce to determine the next leg of the path
    std::vector<Sampler*> m_bounceSamplers;
    // These are sampled N times per bounce (once for each light sample)
    std::vector<Sampler*> m_lightSelectionSamplers;
    std::vector<Sampler*> m_lightElementSamplers;
    std::vector<Sampler*> m_lightSamplers;
    std::vector<Sampler*> m_brdfSamplers;
    
    unsigned int m_numLightSamples;
    unsigned int m_maxRayDepth;
};


//
// Ray tracing
//

// Path trace through the scene, starting with an initial ray.
// Pass along scene information and various samplers so that we can reduce noise
// along the way.
Color pathTrace(const Ray& ray,
                ShapeSet& scene,
                std::vector<Shape*>& lights,
                Rng& rng,
                SamplerContainer& samplers,
                unsigned int pixelSampleIndex);

// Generate a ray-traced image of the scene, with the given camera, resolution,
// and sample settings
Image* raytrace(ShapeSet& scene,
                const Camera& cam,
                size_t width,
                size_t height,
                unsigned int pixelSamplesHint,
                unsigned int lightSamplesHint,
                unsigned int maxRayDepth);


} // namespace Rayito


#endif // __RAYITO_H__
