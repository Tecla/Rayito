////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RRAY_H__
#define __RRAY_H__

#include "RMath.h"


namespace Rayito
{


//
// Ray (directed line segment)
//

// Don't ever start a ray exactly where you previously hit; you must offset it
// a little bit so you don't accidentally 'self-intersect'.
const float kRayTMin = 0.0001f;
// Unless otherwise specified, rays are defined to be able to hit anything as
// far as the computer can see.  You can limit a ray's max if you need to though
// as is done often when calculating shadows, so you only check the range from
// the point on the surface to the point on the light.
const float kRayTMax = 1.0e30f;


struct Ray
{
    Point m_origin;
    Vector m_direction;
    float m_tMax;
    
    // Some sane defaults
    Ray()
        : m_origin(),
          m_direction(0.0f, 0.0f, 1.0f),
          m_tMax(kRayTMax)
    {
        
    }
    
    Ray(const Ray& r)
        : m_origin(r.m_origin),
          m_direction(r.m_direction),
          m_tMax(r.m_tMax)
    {
        
    }
    
    Ray(const Point& origin, const Vector& direction, float tMax = kRayTMax)
        : m_origin(origin),
          m_direction(direction),
          m_tMax(tMax)
    {
        
    }
    
    Ray& operator =(const Ray& r)
    {
        m_origin = r.m_origin;
        m_direction = r.m_direction;
        m_tMax = r.m_tMax;
        return *this;
    }
    
    Point calculate(float t) const { return m_origin + t * m_direction; }
};


//
// Intersection (results from casting a ray)
//

class Shape;
class Material;

struct Intersection
{
    Ray m_ray;
    float m_t;
    Shape *m_pShape;
    Material *m_pMaterial;
    Color m_colorModifier;
    Vector m_normal;
    
    
    Intersection()
        : m_ray(),
          m_t(kRayTMax),
          m_pShape(NULL),
          m_pMaterial(NULL),
          m_colorModifier(1.0f, 1.0f, 1.0f),
          m_normal()
    {
        
    }
    
    Intersection(const Intersection& i)
        : m_ray(i.m_ray),
          m_t(i.m_t),
          m_pShape(i.m_pShape),
          m_pMaterial(i.m_pMaterial),
          m_colorModifier(i.m_colorModifier),
          m_normal(i.m_normal)
    {
        
    }
    
    Intersection(const Ray& ray)
         : m_ray(ray),
           m_t(ray.m_tMax),
           m_pShape(NULL),
           m_pMaterial(NULL),
           m_colorModifier(1.0f, 1.0f, 1.0f),
           m_normal()
    {
        
    }
    
    Intersection& operator =(const Intersection& i)
    {
        m_ray = i.m_ray;
        m_t = i.m_t;
        m_pShape = i.m_pShape;
        m_pMaterial = i.m_pMaterial;
        m_colorModifier = i.m_colorModifier;
        m_normal = i.m_normal;
        return *this;
    }
    
    bool intersected() const { return (m_pShape == NULL) ? false : true; }
    
    Point position() const { return m_ray.calculate(m_t); }
};


} // namespace Rayito


#endif // __RRAY_H__
