////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RSCENE_H__
#define __RSCENE_H__

#include <list>
#include <vector>
#include <algorithm>

#include "RMath.h"
#include "RMaterial.h"
#include "RRay.h"
#include "RSampling.h"
#include "RAccel.h"


namespace Rayito
{


//
// Shapes (scene hierarchy)
//

class Shape
{
public:
    Shape() : m_transform() { }
    
    virtual ~Shape() { }
    
    const Transform& transform() const { return m_transform; }
    Transform&       transform()       { return m_transform; }
    
    // Subclasses must implement this; this is the meat of ray tracing.
    // The first version finds the nearest intersection, the second just tells
    // us if the ray hits anything at all (generally used for shadow rays).
    virtual bool intersect(Intersection& intersection) = 0;
    virtual bool doesIntersect(const Ray& ray) = 0;
    
    // Get bbox of this shape (and its children)
    virtual BBox bbox() = 0;
    // Is the bbox of this shape infinitely big in at least one dimension?
    virtual bool infiniteExtent() const { return false; }
    
    // Prep bboxes, transforms, areas, whatever the shape needs. The renderer
    // calls this on the scene root shape which will prep all of the shapes.
    virtual void prepare() { m_transform.prepare(); }
    
    // Usually for lights: given two random numbers between 0.0 and 1.0, find a
    // location + surface normal on the surface, and return the PDF for how
    // likely the sample was (with respect to solid angle).  Return false if not
    // a surface.
    virtual bool sampleSurface(const Point& refPosition,
                               const Vector& refNormal,
                               float refTime,
                               float u1,
                               float u2,
                               float u3,
                               Point& outPosition,
                               Vector& outNormal,
                               float& outPdf)
    {
        outPdf = 0.0f;
        return false;
    }
    
    // Given a reference position and normal, and a position and normal on the
    // surface of this shape, return the PDF (with respect to solid angle) of
    // having chosen that surface position.
    virtual float pdfSA(const Point& refPosition,
                        const Vector& refNormal,
                        float refTime,
                        const Point& surfPosition,
                        const Vector& surfNormal) const
    {
        Vector incoming = surfPosition - refPosition;
        float dist = incoming.normalize();
        return dist * dist * surfaceAreaPdf() / std::fabs(dot(surfNormal, incoming));
    }
    
    // Return the likelihood of sampling a point on a surface
    virtual float surfaceAreaPdf() const
    {
        return 0.0f;
    }
    
    // Find all lights in the scene starting with this shape
    virtual void findLights(std::vector<Shape*>& outLightList) { }
    
    // Is this shape itself actually a light?
    virtual bool isLight() const { return false; }
    
    // Methods for BVH build
    virtual unsigned int numElements()             const { return 0; }
    virtual BBox         elementBBox(unsigned int) const { return BBox(); }
    virtual float        elementArea(unsigned int) const { return 0; }
    
    // Methods for BVH intersection
    virtual bool intersect(Intersection&, unsigned int)      { return false; }
    virtual bool doesIntersect(const Ray& ray, unsigned int) { return false; }
    
protected:
    Transform m_transform;
};


// List of shapes, so you can aggregate a pile of them
class ShapeSet : public Shape
{
public:
    ShapeSet() : Shape(), m_shapes(), m_infiniteShapes(), m_bvh(*this) { }
    
    virtual ~ShapeSet() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        // Transform ray to local space for intersection
        Ray nonLocalRay = intersection.m_ray;
        intersection.m_ray = intersection.m_ray.transformToLocal(m_transform);
        bool intersectedAny = false;
        for (std::vector<Shape*>::iterator iter = m_infiniteShapes.begin();
             iter != m_infiniteShapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            if (pShape->intersect(intersection))
                intersectedAny = true;
        }
        
        if (m_shapes.size() > 2)
        {
            if (m_bvh.intersect(intersection))
                intersectedAny = true;
        }
        else
        {
            for (std::vector<Shape*>::iterator iter = m_shapes.begin();
                 iter != m_shapes.end();
                 ++iter)
            {
                Shape *pShape = *iter;
                if (pShape->intersect(intersection))
                    intersectedAny = true;
            }
        }
        // Put ray back in non-local space, and patch up normal if we intersected
        if (intersectedAny)
            intersection.m_normal = m_transform.fromLocalNormal(nonLocalRay.m_time, intersection.m_normal);
        intersection.m_ray = nonLocalRay;
        return intersectedAny;
    }
    
    virtual bool doesIntersect(const Ray& ray)
    {
        // Put ray in local space for intersection test
        Ray localRay = ray.transformToLocal(m_transform);
        for (std::vector<Shape*>::iterator iter = m_infiniteShapes.begin();
             iter != m_infiniteShapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            if (pShape->doesIntersect(localRay))
                return true;
        }
        
        if (m_shapes.size() > 2)
        {
            return m_bvh.doesIntersect(localRay);
        }
        for (std::vector<Shape*>::iterator iter = m_shapes.begin();
             iter != m_shapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            if (pShape->doesIntersect(localRay))
                return true;
        }
        return false;
    }
    
    virtual void prepare()
    {
        Shape::prepare();
        for (std::vector<Shape*>::iterator iter = m_infiniteShapes.begin();
             iter != m_infiniteShapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            pShape->prepare();
        }
        for (std::vector<Shape*>::iterator iter = m_shapes.begin();
             iter != m_shapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            pShape->prepare();
        }
        if (m_shapes.size() > 2)
            m_bvh.build();
    }
    
    virtual BBox bbox()
    {
        // Compute combined bbox (in non-local space)
        BBox totalBBox;
        for (size_t ti = 0; ti < m_transform.numKeys(); ++ti)
        {
            float time = m_transform.keyTime(ti);
            for (std::vector<Shape*>::iterator iter = m_shapes.begin();
                 iter != m_shapes.end();
                 ++iter)
            {
                totalBBox = totalBBox.combined((*iter)->bbox().transformFromLocal(time, m_transform));
            }
        }
        return totalBBox;
    }
    
    virtual float surfaceAreaPdf() const
    {
        // TODO: this does not account for scaling!
        float areaTotal = 0.0f;
        for (unsigned int i = 0; i < numElements(); ++i)
        {
            areaTotal += elementArea(i);
        }
        return areaTotal > 0.0f ? 1.0f / areaTotal : 0.0f;
    }
    
    virtual void findLights(std::vector<Shape*>& outLightList)
    {
        for (std::vector<Shape*>::iterator iter = m_shapes.begin();
             iter != m_shapes.end();
             ++iter)
        {
            Shape *pShape = *iter;
            pShape->findLights(outLightList);
        }
    }
    
    void addShape(Shape *pShape)
    {
        if (pShape->infiniteExtent())
            m_infiniteShapes.push_back(pShape);
        else
            m_shapes.push_back(pShape);
    }
    
    void clearShapes() { m_shapes.clear(); m_infiniteShapes.clear(); }
    
    // Methods for BVH build
    virtual unsigned int numElements()                   const { return m_shapes.size(); }
    virtual BBox         elementBBox(unsigned int index) const { return m_shapes[index]->bbox(); }
    virtual float        elementArea(unsigned int index) const { return 1.0f / m_shapes[index]->surfaceAreaPdf(); }
    
    // Methods for BVH intersection
    virtual bool intersect(Intersection& intersection, unsigned int index) { return m_shapes[index]->intersect(intersection); }
    virtual bool doesIntersect(const Ray& ray, unsigned int index)         { return m_shapes[index]->doesIntersect(ray); }
    
protected:
    std::vector<Shape*> m_shapes;
    std::vector<Shape*> m_infiniteShapes;
    Bvh<ShapeSet> m_bvh;
};


// Infinite-extent plane, with option bullseye texturing to make it interesting.
class Plane : public Shape
{
public:
    Plane(const Point& position, const Vector& normal, Material *pMaterial, bool bullseye = false)
        : Shape(),
          m_position(position),
          m_normal(normal.normalized()),
          m_pMaterial(pMaterial),
          m_bullseye(bullseye)
    {
        
    }
    
    virtual ~Plane() { }
    
    virtual bool intersect(Intersection& intersection)
    {
        Ray localRay = intersection.m_ray.transformToLocal(m_transform);
        
        // Plane eqn: ax+by+cz+d=0; another way of writing it is: dot(n, p-p0)=0
        // where n=normal=(a,b,c), and p=(x,y,z), and p0 is position.  Now, p is
        // the ray equation (the intersection point is along the ray): p=origin+t*direction
        // So the plane-ray intersection eqn is dot(n, origin+t*direction-p0)=0.
        // Distribute, and you get:
        //     dot(n, origin) + t*dot(n, direction) - dot(n, p0) = 0
        // Solve for t, and you get:
        //    t = (dot(n, p0) - dot(n, origin)) / dot(n, direction)
        
        // Check if it's even possible to intersect
        float nDotD = dot(m_normal, localRay.m_direction);
        if (nDotD >= 0.0f)
        {
            return false;
        }
        
        float t = (dot(m_position, m_normal) - dot(localRay.m_origin, m_normal)) / nDotD;
        
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
        intersection.m_normal = m_transform.fromLocalNormal(localRay.m_time, m_normal);
        intersection.m_colorModifier = Color(1.0f, 1.0f, 1.0f);
        
        // Hack bullseye pattern to get some variation
        if (m_bullseye && std::fmod((localRay.calculate(t) - m_position).length() * 0.25f, 1.0f) > 0.5f)
        {
            intersection.m_colorModifier *= 0.2f;
        }
        
        return true;
    }

    virtual bool doesIntersect(const Ray& ray)
    {
        Ray localRay = ray.transformToLocal(m_transform);
        
        // Plane eqn: ax+by+cz+d=0; another way of writing it is: dot(n, p-p0)=0
        // where n=normal=(a,b,c), and p=(x,y,z), and p0 is position.  Now, p is
        // the ray equation (the intersection point is along the ray): p=origin+t*direction
        // So the plane-ray intersection eqn is dot(n, origin+t*direction-p0)=0.
        // Distribute, and you get:
        //     dot(n, origin) + t*dot(n, direction) - dot(n, p0) = 0
        // Solve for t, and you get:
        //    t = (dot(n, p0) - dot(n, origin)) / dot(n, direction)
        
        // Check if it's even possible to intersect
        float nDotD = dot(m_normal, localRay.m_direction);
        if (nDotD >= 0.0f)
        {
            return false;
        }
        
        float t = (dot(m_position, m_normal) - dot(localRay.m_origin, m_normal)) / nDotD;
        
        // Make sure t is not behind the ray, and is closer than the current
        // closest intersection.
        if (t >= ray.m_tMax || t < kRayTMin)
        {
            return false;
        }
        
        return true;
    }
    
    virtual BBox bbox()
    {
        return BBox();
    }
    
    virtual bool infiniteExtent() const { return true; }

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
    Sphere(const Point& position = Point(), float radius = 1.0f, Material* pMaterial = NULL)
        : Shape(),
          m_position(position),
          m_radius(radius),
          m_pMaterial(pMaterial)
    {
        
    }
    
    virtual ~Sphere() { }
    
    void setMaterial(Material* pMaterial) { m_pMaterial = pMaterial; }
    
    virtual bool intersect(Intersection& intersection)
    {
        // Transform ray to local space.  Beyond the tranform, we have to move the
        // sphere center to the origin (and the ray along with it).   This makes
        // the intersection logic easier to follow to intersect a sphere at the origin.
        Ray localRay = intersection.m_ray.transformToLocal(m_transform);
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
        
        float discriminant = b * b - 4.0f * a * c;
        if (discriminant < 0.0f)
        {
            // Discriminant less than zero?  No solution => no intersection.
            return false;
        }
        discriminant = std::sqrt(discriminant);
        
        // Compute a more stable form of our param t (t0 = q/a, t1 = c/q)
        // q = -0.5 * (b - sqrt(b * b - 4.0 * a * c)) if b < 0, or
        // q = -0.5 * (b + sqrt(b * b - 4.0 * a * c)) if b >= 0
        float q = (b < 0.0f) ? (-0.5f * (b - discriminant)) : (-0.5f * (b + discriminant));
        
        // Get our final parametric values
        float t0 = q / a;
        float t1 = (q != 0.0f) ? (c / q) : intersection.m_t;
        
        // Swap them so they are ordered closest-to-farthest
        if (t0 > t1)
        {
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }
        
        // Check our intersection for validity against this ray's extents
        if (t0 >= kRayTMin && t0 < intersection.m_t)
        {
            intersection.m_t = t0;
        }
        else if (t1 >= kRayTMin && t1 < intersection.m_t)
        {
            intersection.m_t = t1;
        }
        else
        {
            // Both intersections are outside of the ray's extent
            return false;
        }
        
        // Create our intersection data
        Vector localNorm = localRay.calculate(intersection.m_t);
        Vector worldNorm = m_transform.fromLocalNormal(localRay.m_time, localNorm).normalized();
        
        intersection.m_pShape = this;
        intersection.m_pMaterial = m_pMaterial;
        intersection.m_normal = worldNorm;
        intersection.m_colorModifier = Color(1.0f, 1.0f, 1.0f);
        
        return true;
    }
    
    virtual bool doesIntersect(const Ray& ray)
    {
        // Transform ray to local space.  In this case it's just moving the
        // sphere center to the origin (and the ray along with it).   This makes
        // the intersection logic easier to follow.
        Ray localRay = ray.transformToLocal(m_transform);
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
        
        float discriminant = b * b - 4.0f * a * c;
        if (discriminant < 0.0f)
        {
            // Discriminant less than zero?  No solution => no intersection.
            return false;
        }
        discriminant = std::sqrt(discriminant);
        
        // Compute a more stable form of our param t (t0 = q/a, t1 = c/q)
        // q = -0.5 * (b - sqrt(b * b - 4.0 * a * c)) if b < 0, or
        // q = -0.5 * (b + sqrt(b * b - 4.0 * a * c)) if b >= 0
        float q = (b < 0.0f) ? (-0.5f * (b - discriminant)) : (-0.5f * (b + discriminant));
        
        // Get our final parametric values, see if any intersect
        float t0 = q / a;
        if (t0 >= kRayTMin && t0 < ray.m_tMax)
        {
            return true;
        }
        float t1 = c / q;
        if (q != 0.0f && t1 < ray.m_tMax && t1 >= kRayTMin)
        {
            return true;
        }
        return false;
    }
    
    virtual BBox bbox()
    {
        BBox result;
        for (size_t ti = 0; ti < m_transform.numKeys(); ++ti)
        {
            float time = m_transform.keyTime(ti);
            result = result.combined(BBox(m_position - Point(m_radius),
                                          m_position + Point(m_radius)).transformFromLocal(time, m_transform));
        }
        return result;
    }
    
    // Given two random numbers between 0.0 and 1.0, find a location + surface
    // normal on the surface of the *light*.
    virtual bool sampleSurface(const Point& refPosition,
                               const Vector& refNormal,
                               float refTime,
                               float u1,
                               float u2,
                               float u3,
                               Point& outPosition,
                               Vector& outNormal,
                               float& outPdf)
    {
        // Take care which calculations must be done in local space, and which should be non-local
        Vector localRefPosition = m_transform.toLocalPoint(refTime, refPosition);
        Vector toCenter = m_position - localRefPosition;
        float dist2 = toCenter.length2();
        if (dist2 < m_radius * m_radius * 1.00001f)
        {
            // Inside or on the surface of the sphere?  Just pick a spot to sample.
            outNormal = uniformToSphere(u1, u2);
            outPosition = m_position + outNormal * m_radius;
            outNormal = m_transform.fromLocalNormal(refTime, outNormal);
            outPosition = m_transform.fromLocalPoint(refTime, outPosition);
            Vector toSurf = refPosition - outPosition;
            outPdf = toSurf.length2() * surfaceAreaPdf() / std::fabs(dot(toSurf.normalized(), outNormal));
            return true;
        }
        // Outside the surface of the sphere; fit a cone around the sphere to
        // sample more efficiently within the cone
        float sinThetaMax2 = m_radius * m_radius / dist2;
        float cosThetaMax = std::sqrt(std::max(0.0f, 1.0f - sinThetaMax2));
        Vector x, y, z;
        makeCoordinateSpace(toCenter, x, y, z);
        Vector localCone = uniformToCone(u1, u2, cosThetaMax);
        Vector cone = transformFromLocalCoordinateSpace(localCone, x, y, z).normalized();
        // Make sure the generated direction actually hits the sphere; if not, adjust
        Ray localRay(localRefPosition, cone);
        Ray ray = localRay.transformFromLocal(m_transform);
        Intersection isect(ray);
        if (!intersect(isect))
        {
            isect.m_t = dot(toCenter, cone);
        }
        outPosition = localRay.calculate(isect.m_t);
        outNormal = (outPosition - m_position).normalized();
        outNormal = m_transform.fromLocalNormal(refTime, outNormal);
        outPosition = m_transform.fromLocalPoint(refTime, outPosition);
        outPdf = uniformConePdf(cosThetaMax);
        return true;
    }
    
    virtual float pdfSA(const Point &refPosition,
                        const Vector &refNormal,
                        float refTime,
                        const Point &surfPosition,
                        const Vector &surfNormal) const
    {
        // Take care which calculations must be done in local space, and which should be non-local
        Vector localRefPosition = m_transform.toLocalPoint(refTime, refPosition);
        Vector toCenter = m_position - localRefPosition;
        float dist2 = toCenter.length2();
        if (dist2 < m_radius * m_radius * 1.00001f)
        {
            // Inside or on the surface of the sphere?  Just pick a spot to sample.
            Vector toSurf = refPosition - surfPosition;
            return toSurf.length2() * surfaceAreaPdf() / std::fabs(dot(toSurf.normalized(), surfNormal));
        }
        float sinThetaMax2 = m_radius * m_radius / dist2;
        float cosThetaMax = std::sqrt(std::max(0.0f, 1.0f - sinThetaMax2));
        return uniformConePdf(cosThetaMax);
    }
    
    virtual float surfaceAreaPdf() const
    {
        return 3.0f / (4.0f * M_PI * m_radius * m_radius);
    }

protected:
    Point m_position;
    float m_radius;
    Material *m_pMaterial;
};


} // namespace Rayito


#endif // __RSCENE_H__
