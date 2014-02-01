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


int main(int argc, char **argv)
{
    // The 'scene'
    ShapeSet masterSet;
    
    // Put a ground plane in (with bullseye texture!)
    Plane plane(Point(0.0f, -2.0f, 0.0f),
                Vector(0.0f, 1.0f, 0.0f),
                Color(1.0f, 0.5f, 0.8f));
    masterSet.addShape(&plane);
    
    
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
        // PPMs are top-down, and we're bottom up.  Flip pixel row to be in screen space.
        float yu = 1.0f - (float(y) / float(kHeight - 1));
        
        // For each pixel across the row...
        for (size_t x = 0; x < kWidth; ++x)
        {
            float xu = float(x) / float(kWidth - 1);
            
            // Find where this pixel sample hits in the scene
            Ray ray = makeCameraRay(30.0f,
                                    Point(0.0f, 0.0f, 0.0f),
                                    Point(0.0f, 0.0f, 1.0f),
                                    Point(0.0f, 1.0f, 0.0f),
                                    xu,
                                    yu);
            Intersection intersection(ray);
            bool intersected = masterSet.intersect(intersection);
            
            Color pixelColor(0.0f, 0.0f, 0.0f);
            if (intersected)
            {
                pixelColor = intersection.m_color;
            }
            
            // Li'l debug output if you need
//            std::cout << "Ray: " << ray.m_origin << " -> " << ray.m_direction << " :: " << pixelColor << '\n';
            
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

