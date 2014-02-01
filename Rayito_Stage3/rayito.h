////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RAYITO_H__
#define __RAYITO_H__


#include <cmath>
#include <list>
#include <algorithm>


#ifndef M_PI
    // For some reason, MSVC doesn't define this when <cmath> is included
    #define M_PI 3.14159265358979
#endif


namespace Rayito
{


//
// RGB color class
//
// Operators supported:
//     color = color
//     color + color, color += color
//     color - color, color -= color
//     color * color, color *= color, float * color, color * float, color *= float
//     color / color, color /= color, color / float
//

struct Color
{
    float m_r, m_g, m_b;
    
    Color()                          : m_r(0.0f), m_g(0.0f), m_b(0.0f)    { }
    Color(const Color& c)            : m_r(c.m_r), m_g(c.m_g), m_b(c.m_b) { }
    Color(float r, float g, float b) : m_r(r), m_g(g), m_b(b)             { }
    explicit Color(float f)          : m_r(f), m_g(f), m_b(f)             { }
    
    
    void clamp(float min = 0.0f, float max = 1.0f)
    {
        m_r = std::max(min, std::min(max, m_r));
        m_g = std::max(min, std::min(max, m_g));
        m_b = std::max(min, std::min(max, m_b));
    }
    
    
    Color& operator =(const Color& c)
    {
        m_r = c.m_r;
        m_g = c.m_g;
        m_b = c.m_b;
        return *this;
    }
    
    Color& operator +=(const Color& c)
    {
        m_r += c.m_r;
        m_g += c.m_g;
        m_b += c.m_b;
        return *this;
    }
    
    Color& operator -=(const Color& c)
    {
        m_r -= c.m_r;
        m_g -= c.m_g;
        m_b -= c.m_b;
        return *this;
    }
    
    Color& operator *=(const Color& c)
    {
        m_r *= c.m_r;
        m_g *= c.m_g;
        m_b *= c.m_b;
        return *this;
    }
    
    Color& operator /=(const Color& c)
    {
        m_r /= c.m_r;
        m_g /= c.m_g;
        m_b /= c.m_b;
        return *this;
    }
    
    Color& operator *=(float f)
    {
        m_r *= f;
        m_g *= f;
        m_b *= f;
        return *this;
    }
    
    Color& operator /=(float f)
    {
        m_r /= f;
        m_g /= f;
        m_b /= f;
        return *this;
    }
};


inline Color operator +(const Color& c1, const Color& c2)
{
    return Color(c1.m_r + c2.m_r,
                 c1.m_g + c2.m_g,
                 c1.m_b + c2.m_b);
}


inline Color operator -(const Color& c1, const Color& c2)
{
    return Color(c1.m_r - c2.m_r,
                 c1.m_g - c2.m_g,
                 c1.m_b - c2.m_b);
}


inline Color operator *(const Color& c1, const Color& c2)
{
    return Color(c1.m_r * c2.m_r,
                 c1.m_g * c2.m_g,
                 c1.m_b * c2.m_b);
}


inline Color operator /(const Color& c1, const Color& c2)
{
    return Color(c1.m_r / c2.m_r,
                 c1.m_g / c2.m_g,
                 c1.m_b / c2.m_b);
}


inline Color operator *(const Color& c, float f)
{
    return Color(f * c.m_r,
                 f * c.m_g,
                 f * c.m_b);
}


inline Color operator *(float f, const Color& c)
{
    return Color(f * c.m_r,
                 f * c.m_g,
                 f * c.m_b);
}


inline Color operator /(const Color& c, float f)
{
    return Color(c.m_r / f,
                 c.m_g / f,
                 c.m_b / f);
}


//
// 3D vector class (and associated operations)
//
// Operators supported:
//     vector = vector
//     vector + vector, vector += vector
//     vector - vector, vector -= vector
//     vector * float, vector *= float, float * vector
//     vector / float, vector /= float
//

struct Vector
{
    float m_x, m_y, m_z;
    
    Vector()                          : m_x(0.0f), m_y(0.0f), m_z(0.0f)    { }
    Vector(const Vector& v)           : m_x(v.m_x), m_y(v.m_y), m_z(v.m_z) { }
    Vector(float x, float y, float z) : m_x(x), m_y(y), m_z(z)             { }
    explicit Vector(float f)          : m_x(f), m_y(f), m_z(f)             { }
    
    
    float length2() const { return m_x * m_x + m_y * m_y + m_z * m_z; }
    float length()  const { return std::sqrt(length2()); }
    
    // Returns old length from before normalization (ignore the return value if you don't need it)
    float  normalize()        { float len = length(); *this /= len; return len; }
    // Return a vector in this same direction, but normalized
    Vector normalized() const { Vector r(*this); r.normalize(); return r; }
    
    
    Vector& operator =(const Vector& v)
    {
        m_x = v.m_x;
        m_y = v.m_y;
        m_z = v.m_z;
        return *this;
    }
    
    Vector& operator +=(const Vector& v)
    {
        m_x += v.m_x;
        m_y += v.m_y;
        m_z += v.m_z;
        return *this;
    }
    
    Vector& operator -=(const Vector& v)
    {
        m_x -= v.m_x;
        m_y -= v.m_y;
        m_z -= v.m_z;
        return *this;
    }
    
    Vector& operator *=(float f)
    {
        m_x *= f;
        m_y *= f;
        m_z *= f;
        return *this;
    }
    
    Vector& operator /=(float f)
    {
        m_x /= f;
        m_y /= f;
        m_z /= f;
        return *this;
    }
};


inline Vector operator +(const Vector& v1, const Vector& v2)
{
    return Vector(v1.m_x + v2.m_x,
                  v1.m_y + v2.m_y,
                  v1.m_z + v2.m_z);
}


inline Vector operator -(const Vector& v1, const Vector& v2)
{
    return Vector(v1.m_x - v2.m_x,
                  v1.m_y - v2.m_y,
                  v1.m_z - v2.m_z);
}


inline Vector operator *(const Vector& v, float f)
{
    return Vector(f * v.m_x,
                  f * v.m_y,
                  f * v.m_z);
}


inline Vector operator *(float f, const Vector& v)
{
    return Vector(f * v.m_x,
                  f * v.m_y,
                  f * v.m_z);
}


// dot(v1, v2) = length(v1) * length(v2) * cos(angle between v1, v2)
inline float dot(const Vector& v1, const Vector& v2)
{
    // In cartesian coordinates, it simplifies to this simple calculation:
    return v1.m_x * v2.m_x + v1.m_y * v2.m_y + v1.m_z * v2.m_z;
}


// cross(v1, v2) = length(v1) * length(v2) * sin(angle between v1, v2);
// result is perpendicular to both v1, v2.
inline Vector cross(const Vector& v1, const Vector& v2)
{
    // In cartesian coordinates, it simplifies down to this calculation:
    return Vector(v1.m_y * v2.m_z - v1.m_z * v2.m_y,
                  v1.m_z * v2.m_x - v1.m_x * v2.m_z,
                  v1.m_x * v2.m_y - v1.m_y * v2.m_x);
}


// Oh, by the way, a point can be thought of as just a vector, but where you
// refrain from doing dot/cross/normalize operations on it.
typedef Vector Point;


//
// Ray (directed line segment)
//

// Don't ever start a ray exactly where you previously hit; you must offset it
// a little bit so you don't accidentally 'self-intersect'.
const float kRayTMin = 0.00001f;
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


//
// Material
//

class Material
{
public:
    virtual ~Material() { }
    
    // If the material is emitting, override this
    virtual Color emittance() { return Color(); }
    
    virtual Color shade(const Point& position,
                        const Vector& normal,
                        const Vector& incomingRayDirection,
                        const Vector& lightDirection) = 0;
};


// Lambertian diffuse material
class Lambert : public Material
{
public:
    Lambert(const Color& color) : m_color(color) { }
    
    virtual ~Lambert() { }
    
    virtual Color shade(const Point& position,
                        const Vector& normal,
                        const Vector& incomingRayDirection,
                        const Vector& lightDirectionNorm)
    {
        return std::max(0.0f, dot(lightDirectionNorm, normal)) * m_color;
    }
    
protected:
    Color m_color;
};


// Phong glossy material
class Phong : public Material
{
public:
    Phong(const Color& color, float exponent) : m_color(color), m_exponent(exponent) { }
    
    virtual ~Phong() { }
    
    virtual Color shade(const Point& position,
                        const Vector& normal,
                        const Vector& incomingRayDirection,
                        const Vector& lightDirectionNorm)
    {
        Vector halfVec = (lightDirectionNorm - incomingRayDirection).normalized();
        return std::pow(std::max(0.0f, dot(halfVec, normal)), m_exponent) * m_color;
    }
    
protected:
    Color m_color;
    float m_exponent;
};


// Emitter (light) material
class Emitter : public Material
{
public:
    Emitter(const Color& color, float power) : m_color(color), m_power(power) { }
    
    virtual ~Emitter() { }
    
    virtual Color emittance() { return m_color * m_power; }
    
    virtual Color shade(const Point& position,
                        const Vector& normal,
                        const Vector& incomingRayDirection,
                        const Vector& lightDirectionNorm)
    {
        // Let the emittance take care of business
        return Color();
    }
    
protected:
    Color m_color;
    float m_power;
};


//
// Shapes (scene hierarchy)
//

class Shape
{
public:
    virtual ~Shape() { }
    
    // Subclasses must implement this; this is the meat of ray tracing
    virtual bool intersect(Intersection& intersection) = 0;
    
    // Usually for lights: given two random numbers between 0.0 and 1.0, find a
    // location + surface normal on the surface.  Return false if not a surface.
    virtual bool sampleSurface(float u1,
                               float u2,
                               const Point& referencePosition,
                               Point& outPosition,
                               Vector& outNormal)
    {
        return false;
    }
    
    // Find all lights in the scene starting with this shape
    virtual void findLights(std::list<Shape*>& outLightList) { }
};


// List of shapes, so you can aggregate a pile of them
class ShapeSet : public Shape
{
public:
    virtual ~ShapeSet() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        bool intersectedAny = false;
        for (std::list<Shape*>::iterator iter = m_shapes.begin();
             iter != m_shapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            bool intersected = pShape->intersect(intersection);
            if (intersected)
            {
                intersectedAny = true;
            }
        }
        return intersectedAny;
    }
    
    virtual void findLights(std::list<Shape*>& outLightList)
    {
        for (std::list<Shape*>::iterator iter = m_shapes.begin();
             iter != m_shapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            pShape->findLights(outLightList);
        }
    }
    
    void addShape(Shape *pShape) { m_shapes.push_back(pShape); }
    
    void clearShapes() { m_shapes.clear(); }
    
protected:
    std::list<Shape*> m_shapes;
};


// Light base class, making it easy to find all the lights in the scene.
class Light : public Shape
{
public:
    Light(const Color& c, float power) : m_color(c), m_power(power), m_material(c, power) { }
    
    virtual ~Light() { }
    
    // This *is* a light, so we put ourself on the list
    virtual void findLights(std::list<Shape*>& outLightList) { outLightList.push_back(this); }
    
    virtual Color emitted() const { return m_color * m_power; }
    
protected:
    Color m_color;
    float m_power;
    Emitter m_material;
};


// Area light with a corner and two sides to define a rectangular/parallelipiped shape
class RectangleLight : public Light
{
public:
    RectangleLight(const Point& pos,
                   const Vector& side1,
                   const Vector& side2,
                   const Color& color,
                   float power)
        : Light(color, power), m_position(pos), m_side1(side1), m_side2(side2)
    {
        
    }
    
    virtual ~RectangleLight() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        // This is much like a plane intersection, except we also range check it
        // to make sure it's within the rectangle.  Please see the plane shape
        // intersection method for a little more info.
        
        Vector normal = cross(m_side1, m_side2).normalized();
        float nDotD = dot(normal, intersection.m_ray.m_direction);
        if (nDotD == 0.0f)
        {
            return false;
        }
        
        float t = (dot(m_position, normal) - dot(intersection.m_ray.m_origin, normal)) /
                  dot(intersection.m_ray.m_direction, normal);
        
        // Make sure t is not behind the ray, and is closer than the current
        // closest intersection.
        if (t >= intersection.m_t || t < kRayTMin)
        {
            return false;
        }
        
        // Take the intersection point on the plane and transform it to a local
        // space where we can really easily check if it's in range or not.
        Vector side1Norm = m_side1;
        Vector side2Norm = m_side2;
        float side1Length = side1Norm.normalize();
        float side2Length = side2Norm.normalize();
        
        Point worldPoint = intersection.m_ray.calculate(t);
        Point worldRelativePoint = worldPoint - m_position;
        Point localPoint = Point(dot(worldRelativePoint, side1Norm),
                                 dot(worldRelativePoint, side2Norm),
                                 0.0f);
        
        // Do the actual range check
        if (localPoint.m_x < 0.0f || localPoint.m_x > side1Length ||
            localPoint.m_y < 0.0f || localPoint.m_y > side2Length)
        {
            return false;
        }
        
        // This intersection is the closest so far, so record it.
        intersection.m_t = t;
        intersection.m_pShape = this;
        intersection.m_pMaterial = &m_material;
        intersection.m_colorModifier = Color(1.0f, 1.0f, 1.0f);
        intersection.m_normal = normal;
        // Hit the back side of the light?  We'll count it, so flip the normal
        // to effectively make a double-sided light.
        if (dot(intersection.m_normal, intersection.m_ray.m_direction) > 0.0f)
        {
            intersection.m_normal *= -1.0f;
        }
        
        return true;
    }
    
    // Given two random numbers between 0.0 and 1.0, find a location + surface
    // normal on the surface of the *light*.
    virtual bool sampleSurface(float u1, float u2, const Point& referencePosition, Point& outPosition, Vector& outNormal)
    {
        outNormal = cross(m_side1, m_side2).normalized();
        outPosition = m_position + m_side1 * u1 + m_side2 * u2;
        // Reference point out in back of the light?  That's okay, we'll flip
        // the normal to have a double-sided light.
        if (dot(outNormal, outPosition - referencePosition) > 0.0f)
        {
            outNormal *= -1.0f;
        }
        return true;
    }
    
protected:
    Point m_position;
    Vector m_side1, m_side2;
};


// Area light based on an arbitrary shape from the scene
class ShapeLight : public Light
{
public:
    ShapeLight(Shape *pShape,
               const Color& color,
               float power)
        : Light(color, power), m_pShape(pShape)
    {
        
    }
    
    virtual ~ShapeLight() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        // Forward intersection test on to the shape, but patch in the light material
        if (m_pShape->intersect(intersection))
        {
            intersection.m_pMaterial = &m_material;
            return true;
        }
        return false;
    }
    
    // Given two random numbers between 0.0 and 1.0, find a location + surface
    // normal on the surface of the *light*.
    virtual bool sampleSurface(float u1, float u2, const Point& referencePosition, Point& outPosition, Vector& outNormal)
    {
        // Forward surface sampling on to the shape
        return m_pShape->sampleSurface(u1, u2, referencePosition, outPosition, outNormal);
    }
    
protected:
    Shape *m_pShape;
};


// Infinite-extent plane, with option bullseye texturing to make it interesting.
class Plane : public Shape
{
public:
    Plane(const Point& position, const Vector& normal, Material *pMaterial, bool bullseye = false)
        : m_position(position),
          m_normal(normal.normalized()),
          m_pMaterial(pMaterial),
          m_bullseye(bullseye)
    {
        
    }
    
    virtual ~Plane() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        // Plane eqn: ax+by+cz+d=0; another way of writing it is: dot(n, p-p0)=0
        // where n=normal=(a,b,c), and p=(x,y,z), and p0 is position.  Now, p is
        // the ray equation (the intersection point is along the ray): p=origin+t*direction
        // So the plane-ray intersection eqn is dot(n, origin+t*direction-p0)=0.
        // Distribute, and you get:
        //     dot(n, origin) + t*dot(n, direction) - dot(n, p0) = 0
        // Solve for t, and you get:
        //    t = (dot(n, p0) - dot(n, origin)) / dot(n, direction)
        
        // Check if it's even possible to intersect
        float nDotD = dot(m_normal, intersection.m_ray.m_direction);
        if (nDotD >= 0.0f)
        {
            return false;
        }
        
        float t = (dot(m_position, m_normal) - dot(intersection.m_ray.m_origin, m_normal)) /
                  dot(intersection.m_ray.m_direction, m_normal);
        
        // Make sure t is not behind the ray, and is closer than the current
        // closest intersection.
        if (t >= intersection.m_t || t < kRayTMin)
        {
            return false;
        }
        
        // This intersection is closer, so record it.
        intersection.m_t = t;
        intersection.m_pShape = this;
        intersection.m_pMaterial = m_pMaterial;
        intersection.m_normal = m_normal;
        intersection.m_colorModifier = Color(1.0f, 1.0f, 1.0f);
        
        // Hack bullseye pattern to get some variation
        if (m_bullseye && std::fmod((intersection.position() - m_position).length() * 0.25f, 1.0f) > 0.5f)
        {
            intersection.m_colorModifier *= 0.2f;
        }
        
        return true;
    }

protected:
    Point m_position;
    Vector m_normal;
    Material *m_pMaterial;
    bool m_bullseye;
};


// Sphere (heh, what else?)
class Sphere : public Shape
{
public:
    Sphere(const Point& position, float radius, Material* pMaterial)
        : m_position(position),
          m_radius(radius),
          m_pMaterial(pMaterial)
    {
        
    }
    
    virtual ~Sphere() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        Ray localRay = intersection.m_ray;
        localRay.m_origin -= m_position;
        
        // Ray-sphere intersection can result in either zero, one or two points
        // of intersection.  It turns into a quadratic equation, so we just find
        // the solution using the quadratic formula.  Note that there is a
        // slightly more stable form of it when computing it on a computer, and
        // we use that method to keep everything accurate.
        
        // Calculate quadratic coeffs
        float a = localRay.m_direction.length2();
        float b = 2.0f * dot(localRay.m_direction, localRay.m_origin);
        float c = localRay.m_origin.length2() - m_radius * m_radius;
        
        float t0, t1, discriminant;
        discriminant = b * b - 4.0f * a * c;
        if (discriminant < 0.0f)
        {
            // Discriminant less than zero?  No solution => no intersection.
            return false;
        }
        discriminant = std::sqrt(discriminant);
        
        // Compute a more stable form of our param t (t0 = q/a, t1 = c/q)
        // q = -0.5 * (b - sqrt(b * b - 4.0 * a * c)) if b < 0, or
        // q = -0.5 * (b + sqrt(b * b - 4.0 * a * c)) if b >= 0
        float q;
        if (b < 0.0f)
        {
            q = -0.5f * (b - discriminant);
        }
        else
        {
            q = -0.5f * (b + discriminant);
        }
        
        // Get our final parametric values
        t0 = q / a;
        if (q != 0.0f)
        {
            t1 = c / q;
        }
        else
        {
            t1 = intersection.m_t;
        }
        
        // Swap them so they are ordered right
        if (t0 > t1)
        {
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }
        
        // Check our intersection for validity against this ray's extents
        if (t0 >= intersection.m_t || t1 < kRayTMin)
        {
            return false;
        }
        
        if (t0 >= kRayTMin)
        {
            intersection.m_t = t0;
        }
        else if (t1 < intersection.m_t)
        {
            intersection.m_t = t1;
        }
        else
        {
            return false;
        }
        
        // Create our intersection data
        Point localPos = localRay.calculate(intersection.m_t);
        Vector worldNorm = localPos.normalized();
        
        intersection.m_pShape = this;
        intersection.m_pMaterial = m_pMaterial;
        intersection.m_normal = worldNorm;
        intersection.m_colorModifier = Color(1.0f, 1.0f, 1.0f);
        
        return true;
    }
    
    // Given two random numbers between 0.0 and 1.0, find a location + surface
    // normal on the surface of the *light*.
    virtual bool sampleSurface(float u1, float u2, const Point& referencePosition, Point& outPosition, Vector& outNormal)
    {
        // Pick a random point on the whole sphere surface
        outNormal = uniformToSphere(u1, u2);
        outPosition = outNormal * m_radius + m_position;
        if (dot(outNormal, referencePosition - outPosition) < 0.0f)
        {
            // Point was on the opposite side?  Flip around to make this
            // more reasonable.
            outNormal *= -1.0f;
            outPosition = outNormal * m_radius + m_position;
        }
        return true;
    }

protected:
    Point m_position;
    float m_radius;
    Material *m_pMaterial;
    
    // Helper method for finding a random point on the sphere
    Vector uniformToSphere(float u1, float u2)
    {
        // Find a height uniformly distributed on the sphere
        float z = 1.0f - 2.0f * u1;
        // Find the radius based on that height that sits on the sphere surface
        float radius = std::sqrt(std::max(0.0f, 1.0f - z * z));
        // Find a random angle around the sphere's equator
        float phi = M_PI * 2.0f * u2;
        // And put it all together...
        return Vector(radius * std::cos(phi),
                      radius * std::sin(phi),
                      z);
    }
};


} // namespace Rayito


#endif // __RAYITO_H__
