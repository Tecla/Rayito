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
const size_t kNumPixelSamplesU = 4;
const size_t kNumPixelSamplesV = 4;
const size_t kNumLightSamplesU = 4;
const size_t kNumLightSamplesV = 4;


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


int main(int argc, char **argv)
{
    // Change these to suit yourself, or as an exercise grab them from the commandline
    size_t pixelSamplesHint = 4;
    size_t lightSamplesHint = 4;
    
    // Available materials
    Lambert blueishLambert(Color(0.9f, 0.9f, 1.0f));
    Lambert purplishLambert(Color(0.9f, 0.7f, 0.8f));
    Phong greenishPhong(Color(0.7f, 0.9f, 0.7f), 16.0f);
    
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
    ShapeLight sphereLight(&sphereForLight, Color(1.0f, 1.0f, 0.1f), 4.0f);
    masterSet.addShape(&sphereLight);
    
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
            for (size_t vsi = 0; vsi < pixelSamplesHint; ++vsi)
            {
                for (size_t usi = 0; usi < pixelSamplesHint; ++usi)
                {
                    // Calculate a stratified random position within the pixel
                    // to hide aliasing.  Also, PPMs are top-down, and we're
                    // bottom up.  Flip pixel row to be in screen space.
                    float yu = 1.0f - ((y + (vsi + rng.nextFloat()) / float(kNumPixelSamplesV)) / float(kHeight));
                    float xu = (x + (usi + rng.nextFloat()) / float(kNumPixelSamplesU)) / float(kWidth);
                    
                    // Find where this pixel sample hits in the scene
                    Ray ray = makeCameraRay(45.0f,
                                            Point(0.0f, 5.0f, 15.0f),
                                            Point(0.0f, 0.0f, 0.0f),
                                            Point(0.0f, 1.0f, 0.0f),
                                            xu,
                                            yu);
                    
                    pixelColor += trace(ray, masterSet, lights, rng, lightSamplesHint);
                }
            }
            // Divide by the number of pixel samples (a box filter, essentially)
            pixelColor /= pixelSamplesHint * pixelSamplesHint;
            
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

