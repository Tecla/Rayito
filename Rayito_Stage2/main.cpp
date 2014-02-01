#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "rayito.h"


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


// Turn on PFM writing if you can find a viewer that will read it.  This will
// output an HDR image.
//#define WRITE_PFM 1


// TODO: these should probably be read in as commandline parameters.
const size_t kWidth = 512;
const size_t kHeight = 512;
const size_t kNumPixelSamples = 64;


int main(int argc, char **argv)
{
    // The 'scene'
    ShapeSet masterSet;
    
    // Put a ground plane in (with bullseye texture!)
    Plane plane(Point(0.0f, -2.0f, 0.0f),
                Vector(0.0f, 1.0f, 0.0f),
                Color(1.0f, 1.0f, 1.0f),
                true);
    masterSet.addShape(&plane);
    
    // Add an area light
    RectangleLight areaLight(Point(-2.5f, 2.0f, -2.5f),
                             Vector(5.0f, 0.0f, 0.0f),
                             Vector(0.0f, 0.0f, 5.0f),
                             Color(1.0f, 0.5f, 1.0f),
                             3.0f);
    masterSet.addShape(&areaLight);
    
    // Add another area light below it, darker, that will make a shadow too.
    RectangleLight smallAreaLight(Point(-2.0f, -1.0f, -2.0f),
                                  Vector(4.0f, 0.0f, 0.0f),
                                  Vector(0.0f, 0.0f, 4.0f),
                                  Color(1.0f, 1.0f, 0.5f),
                                  0.75f);
    masterSet.addShape(&smallAreaLight);
    
    // Get light list from the scene
    std::list<Shape*> lights;
    masterSet.findLights(lights);
    
    // Random number generator (for random pixel positions, light positions, etc)
    Rng rng;
    
    
    // Set up the output file (TODO: the filename should probably be a commandline parameter)
    std::ostringstream headerStream;
#if WRITE_PFM
    headerStream << "PF\n";
    headerStream << kWidth << ' ' << kHeight << '\n';
    headerStream << "-1.0\n";
    std::ofstream fileStream("out.pfm", std::ios::out | std::ios::binary);
#else
    headerStream << "P6\n";
    headerStream << kWidth << ' ' << kHeight << '\n';
    headerStream << "255\n";
    std::ofstream fileStream("out.ppm", std::ios::out | std::ios::binary);
#endif
    fileStream << headerStream.str();
    
    // For each row...
    for (size_t y = 0; y < kHeight; ++y)
    {
        // For each pixel across the row...
        for (size_t x = 0; x < kWidth; ++x)
        {
            // For each sample in the pixel...
            Color pixelColor(0.0f, 0.0f, 0.0f);
            for (size_t si = 0; si < kNumPixelSamples; ++si)
            {
                // Calculate a random position within the pixel to hide aliasing.
                // PPMs are top-down, and we're bottom up.  Flip pixel row to be in screen space.
                float yu = 1.0f - ((y + rng.nextFloat()) / float(kHeight - 1));
                float xu = (x + rng.nextFloat()) / float(kWidth - 1);
                
                // Find where this pixel sample hits in the scene
                Ray ray = makeCameraRay(45.0f,
                                        Point(0.0f, 5.0f, 15.0f),
                                        Point(0.0f, 0.0f, 0.0f),
                                        Point(0.0f, 1.0f, 0.0f),
                                        xu,
                                        yu);
                Intersection intersection(ray);
                if (masterSet.intersect(intersection))
                {
                    // Add in emission at intersection
                    pixelColor += intersection.m_emitted;
                    
                    // Find out what lights the intersected point can see
                    Point position = intersection.position();
                    for (std::list<Shape*>::iterator iter = lights.begin();
                         iter != lights.end();
                         ++iter)
                    {
                        // Ask the light for a random position/normal we can use
                        // for lighting
                        Point lightPoint;
                        Vector lightNormal;
                        Light *pLightShape = dynamic_cast<Light*>(*iter);
                        pLightShape->sampleSurface(rng.nextFloat(),
                                                   rng.nextFloat(),
                                                   position,
                                                   lightPoint,
                                                   lightNormal);
                        
                        // Fire a shadow ray to make sure we can actually see
                        // that light position
                        Vector toLight = lightPoint - position;
                        float lightDistance = toLight.normalize();
                        Ray shadowRay(position, toLight, lightDistance);
                        Intersection shadowIntersection(shadowRay);
                        bool intersected = masterSet.intersect(shadowIntersection);
                        
                        if (!intersected || shadowIntersection.m_pShape == pLightShape)
                        {
                            // The light point is visible, so let's add that
                            // lighting contribution
                            float lightAttenuation = std::max(0.0f, dot(intersection.m_normal, toLight));
                            pixelColor += intersection.m_color * pLightShape->emitted() * lightAttenuation;
                        }
                    }
                }
            }
            // Divide by the number of pixel samples (a box filter, essentially)
            pixelColor /= kNumPixelSamples;
            
#if WRITE_PFM
            fileStream << pixelColor.m_r << pixelColor.m_g << pixelColor.m_b;
#else
            // We're writing LDR pixel values, so clamp to 0..1 range first
            pixelColor.clamp();
            // Get 24-bit pixel sample and write it out
            unsigned char r, g, b;
            r = static_cast<unsigned char>(pixelColor.m_r * 255.0f);
            g = static_cast<unsigned char>(pixelColor.m_g * 255.0f);
            b = static_cast<unsigned char>(pixelColor.m_b * 255.0f);
            fileStream << r << g << b;
#endif
        }
    }
    
    // Tidy up (probably unnecessary)
    fileStream.flush();
    fileStream.close();
    
    return 0;
}

