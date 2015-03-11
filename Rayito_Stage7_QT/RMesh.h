////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RMESH_H__
#define __RMESH_H__

#include <list>
#include <vector>
#include <algorithm>

#include "RMath.h"
#include "RMaterial.h"
#include "RRay.h"
#include "RSampling.h"
#include "RAccel.h"
#include "RScene.h"


namespace Rayito
{


// Polygon face
struct Face
{
    // Both of these must be the same size (or else m_normalIndices must be empty)
    // These indices point into the vertices and normals in the owner mesh
    std::vector<unsigned int> m_vertexIndices;
    std::vector<unsigned int> m_normalIndices;
};


// Polygon mesh.  Faces may have 3 or more sides, but each face must be convex
// (no holes or edges going back inside the hull at all).  Faces are triangulated
// by making a triangle fan out from the first vertex.
class Mesh : public Shape
{
public:
    Mesh(const std::vector<Point>& verts,
         const std::vector<Vector>& normals,
         const std::vector<Face>& faces,
         Material* pMaterial)
        : m_vertices(verts),
          m_normals(normals),
          m_faces(faces),
          m_pMaterial(pMaterial),
          m_bbox(),
          m_bvh(*this),
          m_faceAreaCDF(),
          m_totalArea(0.0f)
    {
        
    }
    
    virtual ~Mesh() { }
    
    void setMaterial(Material* pMaterial) { m_pMaterial = pMaterial; }
    
    virtual bool intersect(Intersection& intersection)
    {
        // Transform ray to the local space of our transformation
        Ray nonLocalRay = intersection.m_ray;
        intersection.m_ray = intersection.m_ray.transformToLocal(m_transform);
        // Let the BVH do the work of finding the intersection quickly
        bool intersected = m_bvh.intersect(intersection);
        // Patch ray back to non-local space, and fix up the normal if we intersected
        if (intersected)
            intersection.m_normal = m_transform.fromLocalNormal(nonLocalRay.m_time, intersection.m_normal);
        intersection.m_ray = nonLocalRay;
        return intersected;
    }
    
    virtual bool doesIntersect(const Ray& ray)
    {
        // Let the BVH do the work of finding the intersection quickly (in local space)
        Ray localRay = ray.transformToLocal(m_transform);
        return m_bvh.doesIntersect(localRay);
    }
    
    virtual BBox bbox()
    {
        // This is only valid after prepare() is called; note bbox is already in non-local space
        return m_bbox;
    }
    
    virtual void prepare()
    {
        Shape::prepare();
        
        // Calculate the bounding box (in non-local space!)
        m_bbox = BBox();
        for (size_t ti = 0; ti < m_transform.numKeys(); ++ti)
        {
            float time = m_transform.keyTime(ti);
            for (size_t i = 0; i < m_vertices.size(); ++i)
            {
                m_bbox.expand(m_transform.fromLocalPoint(time, m_vertices[i]));
            }
        }
        
        // Calculate total surface area, and the running total of per-face area.
        // This is used to create the "cumulative distribution function" of the
        // face areas we can use to quickly choose a face proportional to its
        // area based on a random number (this means you can use meshes as area
        // lights).
        m_faceAreaCDF.clear();
        m_faceAreaCDF.reserve(m_faces.size() + 1);
        m_totalArea = 0.0f;
        for (size_t faceIndex = 0; faceIndex < m_faces.size(); ++faceIndex)
        {
            float faceArea = 0.0f;
            for (size_t tri = 0; tri < m_faces[faceIndex].m_vertexIndices.size() - 2; ++tri)
            {
                Point p0 = m_vertices[m_faces[faceIndex].m_vertexIndices[0]];
                Point p1 = m_vertices[m_faces[faceIndex].m_vertexIndices[tri + 1]];
                Point p2 = m_vertices[m_faces[faceIndex].m_vertexIndices[tri + 2]];
                faceArea += cross(p1 - p0, p2 - p0).length() * 0.5f;
            }
            m_faceAreaCDF.push_back(m_totalArea);
            m_totalArea += faceArea;
        }
        m_faceAreaCDF.push_back(m_totalArea);
        
        // Build the BVH so ray intersections are nice and fast
        m_bvh.build();
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
        // Select a face based on a random number (u3), proportional to face
        // surface area; a face with double the surface area of another is twice
        // as likely to be selected.
        std::vector<float>::iterator iter = std::upper_bound(m_faceAreaCDF.begin(),
                                                             m_faceAreaCDF.end(),
                                                             u3 * m_totalArea);
        // Get the face index, taking care to make sure we get a face index in range
        size_t faceIndex;
        if (iter == m_faceAreaCDF.end())
            faceIndex = m_faceAreaCDF.size() - 1;
        else if (iter == m_faceAreaCDF.begin())
            faceIndex = 0;
        else
            faceIndex = std::distance(m_faceAreaCDF.begin(), iter) - 1;
        // Find the actual triangle on the face we are choosing
        float faceArea = m_faceAreaCDF[faceIndex + 1] - m_faceAreaCDF[faceIndex];
        float triangleSelector = (u3 * m_totalArea - m_faceAreaCDF[faceIndex]) / faceArea;
        float triangleAreaSoFar = 0.0f;
        for (size_t tri = 0; tri < m_faces[faceIndex].m_vertexIndices.size() - 2; ++tri)
        {
            Point p0 = m_vertices[m_faces[faceIndex].m_vertexIndices[0]];
            Point p1 = m_vertices[m_faces[faceIndex].m_vertexIndices[tri + 1]];
            Point p2 = m_vertices[m_faces[faceIndex].m_vertexIndices[tri + 2]];
            triangleAreaSoFar += cross(p1 - p0, p2 - p0).length() * 0.5f;
            if (triangleSelector * faceArea < triangleAreaSoFar)
            {
                // We found our triangle!  Now, find out which point on the
                // triangle we selected, and put it in non-local space.
                float alpha = 0.0f, beta = 0.0f;
                uniformToBarycentricTriangle(u1, u2, alpha, beta);
                float gamma = 1.0f - alpha - beta;
                outPosition = p0 * alpha + p1 * beta + p2 * gamma;
                outPosition = m_transform.fromLocalPoint(refTime, outPosition);
                // Calculate normal and the PDF of having selected this position (w.r.t. solid angle)
                outNormal = cross(p1 - p0, p2 - p0);
                outNormal = m_transform.fromLocalNormal(refTime, outNormal).normalized();
                Vector toSurf = refPosition - outPosition;
                outPdf = toSurf.length2() * surfaceAreaPdf() / std::fabs(dot(toSurf.normalized(), outNormal));
                return true;
            }
        }
        return false;
    }
    
    virtual float pdfSA(const Point &refPosition,
                        const Vector &refNormal,
                        float refTime,
                        const Point &surfPosition,
                        const Vector &surfNormal) const
    {
        // Likelihood of having selected this position (w.r.t. solid angle)
        Vector toSurf = refPosition - surfPosition;
        return toSurf.length2() * surfaceAreaPdf() / std::fabs(dot(toSurf.normalized(), surfNormal));
    }
    
    virtual float surfaceAreaPdf() const
    {
        // TODO: this does not account for scaling
        return 1.0f / m_totalArea;
    }
    
    // Methods for BVH build
    
    virtual unsigned int numElements() const { return m_faces.size(); }
    
    virtual BBox elementBBox(unsigned int index) const
    {
        // Build a bbox around the face
        BBox bbox;
        for (size_t i = 0; i < m_faces[index].m_vertexIndices.size(); ++i)
        {
            bbox.expand(m_vertices[m_faces[index].m_vertexIndices[i]]);
        }
        return bbox;
    }
    
    virtual float elementArea(unsigned int index) const
    {
        return m_faceAreaCDF[index + 1] - m_faceAreaCDF[index];
    }
    
    // Methods for BVH intersection
    
    virtual bool intersect(Intersection& intersection, unsigned int index)
    {
        // Intersect the triangles of the face, bailing if we find one; we can
        // to this because we assume and hope the triangles are coplanar with
        // each other.  The first one to claim the intersection thus wins.
        bool intersectAny = false;
        for (size_t i = 0; i < m_faces[index].m_vertexIndices.size() - 2; ++i)
        {
            if (intersectTri(index, i, intersection))
                intersectAny = true;
        }
        return intersectAny;
    }
    
    virtual bool doesIntersect(const Ray& ray, unsigned int index)
    {
        // Intersect the triangles of the face
        for (size_t i = 0; i < m_faces[index].m_vertexIndices.size() - 2; ++i)
        {
            if (doesIntersectTri(index, i, ray))
                return true;
        }
        return false;
    }

protected:
    std::vector<Point> m_vertices;
    std::vector<Vector> m_normals;
    std::vector<Face> m_faces;
    Material *m_pMaterial;
    BBox m_bbox;
    Bvh<Mesh> m_bvh;
    std::vector<float> m_faceAreaCDF;
    float m_totalArea;
    
    bool intersectTri(unsigned int faceIndex, unsigned int tri, Intersection& intersection)
    {
        unsigned int v0 = m_faces[faceIndex].m_vertexIndices[0];
        unsigned int v1 = m_faces[faceIndex].m_vertexIndices[tri + 1];
        unsigned int v2 = m_faces[faceIndex].m_vertexIndices[tri + 2];
        
        // Moller-Trumbore ray-triangle intersection test.  The point here is to
        // find the barycentric coordinates of the triangle where the ray hits
        // the plane the triangle lives in.  If the barycentric coordinates
        // alpha, beta, gamma all add up to 1 (and each is in the 0.0 to 1.0
        // range) then we have a valid intersection in the triangle.  Then, each
        // of alpha, beta, gamma are the amounts of influence each vertex has
        // on the values at the intersection.  So if we store things at the
        // vertices (like normals, UVs, colors, etc) we can just weight them
        // with the barycentric coordinates to get the interpolated result.
        
        Vector v0To1 = m_vertices[v1] - m_vertices[v0];
        Vector v0To2 = m_vertices[v2] - m_vertices[v0];
        Vector gnormal = cross(v0To1, v0To2);
        float det = -dot(intersection.m_ray.m_direction, gnormal);
        if (det == 0.0f)
            return false;
        
        Vector rOriginToV0 = m_vertices[v0] - intersection.m_ray.m_origin;
        Vector rayVertCross = cross(intersection.m_ray.m_direction, rOriginToV0);
        Vector rOriginToV1 = m_vertices[v1] - intersection.m_ray.m_origin;
        float invDet = 1.0f / det;
        
        // Calculate barycentric gamma coord
        float gamma = -dot(rOriginToV1, rayVertCross) * invDet;
        if (gamma < 0.0f || gamma > 1.0f)
            return false;
        
        Vector rOriginToV2 = m_vertices[v2] - intersection.m_ray.m_origin;
        
        // Calculate barycentric beta coord
        float beta = dot(rOriginToV2, rayVertCross) * invDet;
        if (beta < 0.0f || beta + gamma > 1.0f)
            return false;
        
        float t = -dot(rOriginToV0, gnormal) * invDet;
        if (t < kRayTMin || t >= intersection.m_t)
            return false;
        
        float alpha = 1.0f - beta - gamma;
    
        // Calculate shading normal...
        Vector shadingNormal;
        if (!m_faces[faceIndex].m_normalIndices.empty())
        {
            // We have normals stored at the vertices, so use them.
            unsigned int n0 = m_faces[faceIndex].m_normalIndices[0];
            unsigned int n1 = m_faces[faceIndex].m_normalIndices[tri + 1];
            unsigned int n2 = m_faces[faceIndex].m_normalIndices[tri + 2];
            
            // Weight normals at each vertex by barycentric coords to create
            // the interpolated normal at the intersection point.
            shadingNormal = (m_normals[n0] * alpha) +
                            (m_normals[n1] * beta) +
                            (m_normals[n2] * gamma);
            shadingNormal.normalize();
        }
        else
        {
            // Use the geometric (flat-shaded) normal
            shadingNormal = gnormal.normalized();
        }
        
        intersection.m_t = t;
        intersection.m_pShape = this;
        intersection.m_pMaterial = m_pMaterial;
        intersection.m_normal = shadingNormal;
        intersection.m_colorModifier = Color(1.0f);
        
        return true;
    }
    
    bool doesIntersectTri(unsigned int faceIndex, unsigned int tri, const Ray& ray)
    {
        unsigned int v0 = m_faces[faceIndex].m_vertexIndices[0];
        unsigned int v1 = m_faces[faceIndex].m_vertexIndices[tri + 1];
        unsigned int v2 = m_faces[faceIndex].m_vertexIndices[tri + 2];
        
        // Moller-Trumbore ray-triangle intersection test.  The point here is to
        // find the barycentric coordinates of the triangle where the ray hits
        // the plane the triangle lives in.  If the barycentric coordinates
        // alpha, beta, gamma all add up to 1 (and each is in the 0.0 to 1.0
        // range) then we have a valid intersection in the triangle.
        
        Vector v0To1 = m_vertices[v1] - m_vertices[v0];
        Vector v0To2 = m_vertices[v2] - m_vertices[v0];
        Vector gnormal = cross(v0To1, v0To2);
        float det = -dot(ray.m_direction, gnormal);
        if (det == 0.0f)
            return false;
        
        Vector rOriginToV0 = m_vertices[v0] - ray.m_origin;
        Vector rayVertCross = cross(ray.m_direction, rOriginToV0);
        Vector rOriginToV1 = m_vertices[v1] - ray.m_origin;
        float invDet = 1.0f / det;
        
        // Calculate barycentric gamma coord
        float gamma = -dot(rOriginToV1, rayVertCross) * invDet;
        if (gamma < 0.0f || gamma > 1.0f)
            return false;
        
        Vector rOriginToV2 = m_vertices[v2] - ray.m_origin;
        
        // Calculate barycentric beta coord
        float beta = dot(rOriginToV2, rayVertCross) * invDet;
        if (beta < 0.0f || beta + gamma > 1.0f)
            return false;
        
        float t = -dot(rOriginToV0, gnormal) * invDet;
        if (t < kRayTMin || t >= ray.m_tMax)
            return false;
        
        return true;
    }
};


Mesh* createFromOBJFile(const char* filename);


} // namespace Rayito


#endif // __RMESH_H__
