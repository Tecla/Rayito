#include <string>
#include <iostream>
#include <sstream>

#include "rayito.h"

#include <QThread>


using namespace Rayito;


// A few debug print helpers for Color and Vector, in case we need them.

std::ostream& operator <<(std::ostream& stream, const Color& c)
{
    stream << '(' << c.m_r << ", " << c.m_g << ", " << c.m_b << ')';
    return stream;
}

std::ostream& operator <<(std::ostream& stream, const Vector& v)
{
    stream << '[' << v.m_x << ", " << v.m_y << ", " << v.m_z << ']';
    return stream;
}


// Marsaglia multiply-with-carry psuedo random number generator.  It's very fast
// and has good distribution properties.  Has a period of 2^60. See
// http://groups.google.com/group/sci.crypt/browse_thread/thread/ca8682a4658a124d/
struct Rng
{
    unsigned int m_z, m_w;
    
    Rng(unsigned int z = 362436069, unsigned int w = 521288629) : m_z(z), m_w(w) { }
    
    
    // Returns a 'canonical' float from [0,1)
    float nextFloat()
    {
        unsigned int i = nextUInt32();
        return i * 2.328306e-10f;
    }
 
    // Returns an int with random bits set
    unsigned int nextUInt32()
    {
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
        return (m_z << 16) + m_w;  /* 32-bit result */
    }
};


// Set up a camera ray given the look-at spec, FOV, and screen position to aim at.
Ray makeCameraRay(float fieldOfViewInDegrees,
                  const Point& origin,
                  const Vector& target,
                  const Vector& targetUpDirection,
                  float xScreenPos0To1,
                  float yScreenPos0To1)
{
    Vector forward = (target - origin).normalized();
    Vector right = cross(forward, targetUpDirection).normalized();
    Vector up = cross(right, forward).normalized();
    
    // Convert to radians, as that is what the math calls expect
    float tanFov = std::tan(fieldOfViewInDegrees * M_PI / 180.0f);
    
    Ray ray;

    // Set up ray info
    ray.m_origin = origin;
    ray.m_direction = forward +
                      right * ((xScreenPos0To1 - 0.5f) * tanFov) +
                      up * ((yScreenPos0To1 - 0.5f) * tanFov);
    ray.m_direction.normalize();
    
    return ray;
}


Color trace(const Ray& ray, ShapeSet& scene, std::list<Shape*>& lights, Rng& rng, size_t lightSamplesHint)
{
    Color result = Color(0.0f, 0.0f, 0.0f);
    
    // Trace the initial ray to see if we hit anything
    Intersection intersection(ray);
    if (!scene.intersect(intersection))
    {
        // No hit, return black (background)
        return result;
    }
    
    const size_t numLightSamplesU = lightSamplesHint;
    const size_t numLightSamplesV = lightSamplesHint;
    
    // Add in emission at intersection
    result += intersection.m_pMaterial->emittance();
    
    // Find out what lights the intersected point can see
    Point position = intersection.position();
    for (std::list<Shape*>::iterator iter = lights.begin();
         iter != lights.end();
         ++iter)
    {
        // Sample the light (with stratified random sampling to reduce noise)
        Color lightResult = Color(0.0f, 0.0f, 0.0f);
        for (size_t lsv = 0; lsv < numLightSamplesV; ++lsv)
        {
            for (size_t lsu = 0; lsu < numLightSamplesU; ++lsu)
            {
                // Ask the light for a random position/normal we can use
                // for lighting
                Point lightPoint;
                Vector lightNormal;
                Light *pLightShape = dynamic_cast<Light*>(*iter);
                pLightShape->sampleSurface((lsu + rng.nextFloat()) / float(numLightSamplesU),
                                           (lsv + rng.nextFloat()) / float(numLightSamplesV),
                                           position,
                                           lightPoint,
                                           lightNormal);
                
                // Fire a shadow ray to make sure we can actually see
                // that light position
                Vector toLight = lightPoint - position;
                float lightDistance = toLight.normalize();
                Ray shadowRay(position, toLight, lightDistance - kRayTMin);
                Intersection shadowIntersection(shadowRay);
                bool intersected = scene.intersect(shadowIntersection);
                
                if (!intersected || shadowIntersection.m_pShape == pLightShape)
                {
                    // The light point is visible, so let's add that
                    // lighting contribution
                    lightResult += pLightShape->emitted() *
                        intersection.m_colorModifier *
                        intersection.m_pMaterial->shade(position,
                                                        intersection.m_normal,
                                                        ray.m_direction,
                                                        toLight);
                }
            }
        }
        lightResult /= numLightSamplesU * numLightSamplesV;
        
        result += lightResult;
    }
    
    return result;
}


//
// RenderThread works on a small chunk of the image
//
class RenderThread : public QThread
{
public:
    RenderThread(size_t xstart, size_t xend, size_t ystart, size_t yend,
                 Image *pImage,
                 ShapeSet& masterSet,
                 std::list<Shape*>& lights,
                 size_t pixelSamplesHint, size_t lightSamplesHint)
        : m_xstart(xstart), m_xend(xend), m_ystart(ystart), m_yend(yend),
          m_pImage(pImage), m_masterSet(masterSet), m_lights(lights),
          m_pixelSamplesHint(pixelSamplesHint), m_lightSamplesHint(lightSamplesHint) { }
    
protected:
    virtual void run()
    {
        // Random number generator (for random pixel positions, light positions, etc)
        // We seed the generator for this render thread based on something that
        // doesn't change, but gives us a good variable seed for each thread.
        Rng rng(static_cast<unsigned int>(((m_xstart << 16) | m_xend) ^ m_xstart),
                static_cast<unsigned int>(((m_ystart << 16) | m_yend) ^ m_ystart));
        
        const size_t numPixelSamplesU = m_pixelSamplesHint;
        const size_t numPixelSamplesV = m_pixelSamplesHint;
        
        float aspectRatioXToY = float(m_pImage->width()) / float(m_pImage->height());
        
        // For each pixel row...
        for (size_t y = m_ystart; y < m_yend; ++y)
        {
            // For each pixel across the row...
            for (size_t x = m_xstart; x < m_xend; ++x)
            {
                // For each sample in the pixel...
                Color pixelColor(0.0f, 0.0f, 0.0f);
                for (size_t vsi = 0; vsi < numPixelSamplesV; ++vsi)
                {
                    for (size_t usi = 0; usi < numPixelSamplesU; ++usi)
                    {
                        // Calculate a stratified random position within the pixel
                        // to hide aliasing.  Also, PPMs are top-down, and we're
                        // bottom up.  Flip pixel row to be in screen space.
                        float yu = 1.0f - ((y + (vsi + rng.nextFloat()) / float(numPixelSamplesV)) / float(m_pImage->height()));
                        float xu = (x + (usi + rng.nextFloat()) / float(numPixelSamplesU)) / float(m_pImage->width());
                        
                        // Find where this pixel sample hits in the scene
                        Ray ray = makeCameraRay(45.0f,
                                                Point(0.0f, 5.0f, 15.0f),
                                                Point(0.0f, 0.0f, 0.0f),
                                                Point(0.0f, 1.0f, 0.0f),
                                                (xu - 0.5f) * aspectRatioXToY + 0.5f,
                                                yu);
                        
                        pixelColor += trace(ray, m_masterSet, m_lights, rng, m_lightSamplesHint);
                    }
                }
                // Divide by the number of pixel samples (a box filter, essentially)
                pixelColor /= numPixelSamplesU * numPixelSamplesV;
                
                // Store off the computed pixel in a big buffer
                m_pImage->pixel(x, y) = pixelColor;
            }
        }
    }
    
    size_t m_xstart, m_xend, m_ystart, m_yend;
    Image *m_pImage;
    ShapeSet& m_masterSet;
    std::list<Shape*>& m_lights;
    size_t m_pixelSamplesHint, m_lightSamplesHint;
};


Image* raytrace(size_t width, size_t height, size_t pixelSamplesHint, size_t lightSamplesHint)
{
    // Available materials
    Lambert blueishLambert(Color(0.8f, 0.8f, 1.0f));
    Lambert purplishLambert(Color(0.9f, 0.5f, 0.7f));
    Phong greenishPhong(Color(0.3f, 0.9f, 0.3f), 16.0f);
    
    // The 'scene'
    ShapeSet masterSet;
    
    // Put a ground plane in (with bullseye texture!)
    Plane plane(Point(0.0f, -2.0f, 0.0f),
                Vector(0.0f, 1.0f, 0.0f),
                &blueishLambert,
                true);
    masterSet.addShape(&plane);
    
    Sphere sphere1(Point(3.0f, -1.0f, 0.0f),
                   1.0f,
                   &purplishLambert);
    masterSet.addShape(&sphere1);
    
    Sphere sphere2(Point(-3.0f, 0.0f, -2.0f),
                   2.0f,
                   &greenishPhong);
    masterSet.addShape(&sphere2);
    
    // Add an area light
    RectangleLight areaLight(Point(-2.5f, 4.0f, -2.5f),
                             Vector(5.0f, 0.0f, 0.0f),
                             Vector(0.0f, 0.0f, 5.0f),
                             Color(1.0f, 1.0f, 1.0f),
                             1.0f);
    masterSet.addShape(&areaLight);

    Sphere sphereForLight(Point(0.0f, 0.0f, 2.0f),
                          1.0f,
                          &blueishLambert);
    ShapeLight sphereLight(&sphereForLight, Color(1.0f, 1.0f, 0.3f), 4.0f);
    masterSet.addShape(&sphereLight);
    
    // Get light list from the scene
    std::list<Shape*> lights;
    masterSet.findLights(lights);
    
    // Set up the output image
    Image *pImage = new Image(width, height);
    
    // Set up render threads; we make as much as 16 chunks of the image that
    // can render in parallel.
    
    // Chunk size is the number of pixels per image chunk (we have to take care
    // to deal with tiny images)
    size_t xChunkSize = width >= 4 ? width / 4 : 1;
    size_t yChunkSize = height >= 4 ? height / 4 : 1;
    // Chunks are the number of chunks in each dimension we can chop the image
    // into (again, taking care to deal with tiny images, and also images that
    // don't divide clealy into 4 chunks)
    size_t xChunks = width > 4 ? width / xChunkSize : 1;
    size_t yChunks = height > 4 ? height / yChunkSize : 1;
    if (xChunks * xChunkSize < width) xChunks++;
    if (yChunks * yChunkSize < height) yChunks++;
    
    size_t numRenderThreads = xChunks * yChunks;
    RenderThread **renderThreads = new RenderThread*[numRenderThreads];
    
    // Launch render threads
    for (size_t yc = 0; yc < yChunks; ++yc)
    {
        // Get the row start/end (making sure the last chunk doesn't go off the end)
        size_t yStart = yc * yChunkSize;
        size_t yEnd = std::min((yc + 1) * yChunkSize, height);
        for (size_t xc = 0; xc < xChunks; ++xc)
        {
            // Get the column start/end (making sure the last chunk doesn't go off the end)
            size_t xStart = xc * xChunkSize;
            size_t xEnd = std::min((xc + 1) * xChunkSize, width);
            // Render the chunk!
            renderThreads[yc * xChunks + xc] = new RenderThread(xStart,
                                                                xEnd,
                                                                yStart,
                                                                yEnd,
                                                                pImage,
                                                                masterSet,
                                                                lights,
                                                                pixelSamplesHint,
                                                                lightSamplesHint);
            renderThreads[yc * xChunks + xc]->start();
        }
    }
    
    // Wait until the render finishes
    bool stillRunning;
    do
    {
        // See if any render thread is still going...
        stillRunning = false;
        for (size_t i = 0; i < numRenderThreads; ++i)
        {
            if (renderThreads[i]->isRunning())
            {
                stillRunning = true;
                break;
            }
        }
        if (stillRunning)
        {
            // Give up the CPU so the render threads can do their thing
            QThread::yieldCurrentThread();
        }
    } while (stillRunning);
    
    // Clean up render thread objects
    for (size_t i = 0; i < numRenderThreads; ++i)
    {
        delete renderThreads[i];
    }
    delete[] renderThreads;
    
    // We made a picture!
    return pImage;
}

