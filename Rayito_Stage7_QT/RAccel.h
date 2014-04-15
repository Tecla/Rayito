////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RACCEL_H__
#define __RACCEL_H__

#include <limits>
#include <algorithm>

#include "RMath.h"
#include "RRay.h"


namespace Rayito
{


// Axis-aligned bounding box, with plenty of handy utilities inside it
struct BBox
{
    Point m_min, m_max;
    
    BBox() : m_min(std::numeric_limits<float>::max()), m_max(-std::numeric_limits<float>::max()) { }
    BBox(const BBox& bbox) : m_min(bbox.m_min), m_max(bbox.m_max) { }
    BBox(const Point& minCorner, const Point& maxCorner) : m_min(minCorner), m_max(maxCorner) { }
    
    BBox& operator =(const BBox& bbox)
    {
        m_min = bbox.m_min;
        m_max = bbox.m_max;
        return *this;
    }
    
    bool valid() const { return m_min.m_x < m_max.m_x && m_min.m_y < m_max.m_y && m_min.m_z < m_max.m_z; }
    bool empty() const { return !valid(); }
    
    bool intersects(const Ray& ray, float& inout_t0, float& inout_t1) const
    {
        // Ray-box intersection, recording the distances along the ray it enters/exits
        Vector invDir = 1.0f / ray.m_direction;
        return intersects(ray.m_origin, invDir, inout_t0, inout_t1);
    }
    
    bool intersects(const Point& origin, const Vector& invDir, float& inout_t0, float& inout_t1) const
    {
        // Ray-box intersection, recording the distances along the ray it enters/exits
        Vector vt0 = (m_min - origin) * invDir;
        Vector vt1 = (m_max - origin) * invDir;
        Vector vtNear = min(vt0, vt1);
        Vector vtFar = max(vt0, vt1);
        float btMin = vtNear.maxComponent();
        float btMax = vtFar.minComponent();
        inout_t0 = std::max(btMin, inout_t0);
        inout_t1 = std::min(btMax, inout_t1);
        return inout_t0 <= inout_t1;
    }
    
    BBox combined(const BBox& bbox) const
    {
        // Union of the two bboxes
        return BBox(min(m_min, bbox.m_min), max(m_max, bbox.m_max));
    }

    void expand(const Point& p)
    {
        // Expand the bbox to include the point
        m_min = min(m_min, p);
        m_max = max(m_max, p);
    }
    
    bool overlaps(const BBox& bbox) const
    {
        return intersection(bbox).valid();
    }
    
    bool contains(const Point& p) const
    {
        // Is the point inside the bbox?
        return m_min.m_x <= p.m_x && m_max.m_x >= p.m_x &&
               m_min.m_y <= p.m_y && m_max.m_y >= p.m_y &&
               m_min.m_z <= p.m_z && m_max.m_z >= p.m_z;
    }
    
    BBox intersection(const BBox& bbox) const
    {
        // Get the bbox that represents the overlap of these two (if any)
        return BBox(max(m_min, bbox.m_min), min(m_max, bbox.m_max));
    }
    
    BBox transformFromLocal(float time, const Transform& txform)
    {
        // Transform each corner of the local box, and then put a box around
        // those points.  This is the best we can do with an axis-aligned bbox.
        Point corners[8] =
        {
            m_min,
            Point(m_min.m_x, m_min.m_y, m_max.m_z),
            Point(m_min.m_x, m_max.m_y, m_min.m_z),
            Point(m_min.m_x, m_max.m_y, m_max.m_z),
            Point(m_max.m_x, m_min.m_y, m_min.m_z),
            Point(m_max.m_x, m_min.m_y, m_max.m_z),
            Point(m_max.m_x, m_max.m_y, m_min.m_z),
            m_max
        };
        BBox result;
        for (int i = 0; i < 8; i++)
        {
            result.expand(txform.fromLocalPoint(time, corners[i]));
        }
        return result;
    }
};


// BVH node flags: split axis takes up the first two bits, and the leaf vs interior takes the 3rd bit
typedef unsigned int BvhNodeFlags;
const BvhNodeFlags kSplitX = 0;
const BvhNodeFlags kSplitY = 1;
const BvhNodeFlags kSplitZ = 2;
const BvhNodeFlags kSplitFlags = 0x3;
const BvhNodeFlags kLeafNode = 0x4;
// 29 bits left over for # of prims if we ever get around to that


// BVH node: it has a bounding box around the contents of the node, flags that
// indicate if it's a leaf node (has no child nodes, holds a primitive) or is an
// interior node (has two child nodes, and has a splitting axis).  Note that the
// children nodes will always be stored consecutively, so we only have to store
// the index to the first child node.  Also note that leaf nodes will store the
// primitive index, but don't need the child node index (and vice-versa), so we
// stick them in a union because the node uses either the child index or the
// primitive index, but not both at the same time (ever).
struct BvhNode
{
    BBox m_bbox;
    union
    {
        unsigned int m_firstChild;
        unsigned int m_prim;
    };
    BvhNodeFlags m_flags;
    
    BvhNode() { }
    BvhNode(const BvhNode& n) : m_bbox(n.m_bbox), m_prim(n.m_prim), m_flags(n.m_flags) { }
    
    BvhNode& operator =(const BvhNode& n)
    {
        m_bbox = n.m_bbox;
        m_prim = n.m_prim;
        m_flags = n.m_flags;
        return *this;
    }
    
    bool leafNode()     const { return (m_flags & kLeafNode) != 0; }
    bool interiorNode() const { return (m_flags & kLeafNode) == 0; }
    
    // Splitting axis
    BvhNodeFlags split() const { return m_flags & kSplitFlags; }
    
    // Children nodes stored consecutively
    unsigned int leftChildIndex()  const { return m_firstChild; }
    unsigned int rightChildIndex() const { return m_firstChild + 1; }
    
    unsigned int prim() const { return m_prim; }
};


/*
 * BVH (bounding volume hierarchy).  This is a binary tree data spatial data
 * structure used to find ray intersections much more quickly (algorithmically
 * it does so in O(log N) time, instead of O(N) time if we didn't have a BVH).
 * Each node in the tree either stores a primitive (a leaf node) or a pointer
 * to two child BVH nodes.  Each node has a bounding box, which *may* overlap
 * with its sibling node.
 * 
 * This particular BVH uses spatial splits, so the trees it generates are not
 * amazingly efficient, but they're way, WAY better than nothing.  A more
 * advanced implementation would use SAH (surface-area hueristic) to pick better
 * splitting axis locations.  Those trees generally take longer to build, but
 * the time to actually trace rays through them is faster.
 * 
 * The template param type for the BVH must have the following methods:
 *     unsigned int numElements() const;
 *     BBox elementBBox(unsigned int index) const;
 *     float elementArea(unsigned int index) const;
 *     bool intersect(Intersection& intersection, unsigned int elementIndex);
 *     bool doesIntersect(const Ray& ray, unsigned int elementIndex);
 * The first two methods are used during building, the second two during tracing.
 */
template<typename T>
class Bvh
{
public:
    Bvh(T& object);
    
    ~Bvh();
    
    // Call this before tracing any rays through the BVH!
    bool build();
    
    // Trace rays, forwarding final ray intersection logic to the object
    bool intersect(Intersection& intersection);
    bool doesIntersect(const Ray& ray);
    
private:
    T& m_object;
    BvhNode *m_nodes;
    unsigned int m_numNodes;
    
    // A couple of helper structs for building the BVH
    
    // The bbox and actual primitive index for each primitive are needed during the build
    struct BuildElement
    {
        unsigned int m_prim;
        BBox m_bbox;
    };
    
    // At each step of the build, we have to divide the primitives so that those
    // to each side of the splitting axis get put in the right part of the list.
    // This helper struct is used to decide which part of the list each primitive
    // goes in.
    struct BuildElementPredicate
    {
        float m_splitAxis;
        BvhNodeFlags m_split;
        
        BuildElementPredicate(float splitAxis, BvhNodeFlags split)
            : m_splitAxis(splitAxis), m_split(split) { }
        
        bool operator ()(const BuildElement& elem)
        {
            return (m_split == kSplitX && m_splitAxis < (elem.m_bbox.m_max.m_x + elem.m_bbox.m_min.m_x) * 0.5f) ||
                   (m_split == kSplitY && m_splitAxis < (elem.m_bbox.m_max.m_y + elem.m_bbox.m_min.m_y) * 0.5f) ||
                   (m_split == kSplitZ && m_splitAxis < (elem.m_bbox.m_max.m_z + elem.m_bbox.m_min.m_z) * 0.5f);
        }
    };
    
    // At each step of the build, this is called recursively to fill out a BVH node
    bool buildRange(BuildElement *permutedElements,
                    unsigned int begin, unsigned int end,
                    unsigned int nodeIndex, const BBox& nodeBBox);
};


template<typename T>
Bvh<T>::Bvh(T& object)
    : m_object(object), m_nodes(NULL), m_numNodes(0)
{
    
}

template<typename T>
Bvh<T>::~Bvh()
{
    if (m_nodes != NULL) delete[] m_nodes;
}

template<typename T>
bool Bvh<T>::build()
{
    // Prep for the build: get primitive bboxes, indices, and set up the actual
    // BVH node storage so we can start filling it out.
    unsigned int numElems = m_object.numElements();
    if (numElems == 0)
        return true;
    
    BuildElement *elems = new BuildElement[numElems];
    BBox totalBBox;
    for (unsigned int i = 0; i < numElems; ++i)
    {
        elems[i].m_prim = i;
        elems[i].m_bbox = m_object.elementBBox(i);
        totalBBox = totalBBox.combined(elems[i].m_bbox);
    }
    // There can be exactly this many BVH nodes total.  It just works.
    m_nodes = new BvhNode[numElems * 2 - 1];
    // We start with one node already set aside (the root node)
    m_numNodes = 1;
    // Start building (with the root node)
    bool built = buildRange(elems, 0, numElems, 0, totalBBox);
    // Clean up temp help for building and get outta here
    delete[] elems;
    return built;
}

template<typename T>
bool Bvh<T>::buildRange(BuildElement *permutedElements,
                        unsigned int begin, unsigned int end,
                        unsigned int nodeIndex, const BBox& nodeBBox)
{
    // Is there only one primitive?  If so, make this a leaf node.
    if (end - begin <= 1)
    {
        m_nodes[nodeIndex].m_flags = kLeafNode;
        m_nodes[nodeIndex].m_bbox = nodeBBox;
        m_nodes[nodeIndex].m_prim = permutedElements[begin].m_prim;
        return true;
    }
    
    // Interior node...
    
    // Pick split axis
    Vector extents = nodeBBox.m_max - nodeBBox.m_min;
    BvhNodeFlags split;
    if (extents.m_x > extents.m_y)
    {
        if (extents.m_x > extents.m_z)
            split = kSplitX;
        else
            split = kSplitZ;
    }
    else if (extents.m_y > extents.m_z)
        split = kSplitY;
    else
        split = kSplitZ;
    
    // Pick split axis location (this is a vanilla spatial split, a SAH tree
    // build would do something more sophisticated here).
    float splitAxis;
    if (split == kSplitX)
        splitAxis = (nodeBBox.m_max.m_x + nodeBBox.m_min.m_x) * 0.5f;
    else if (split == kSplitY)
        splitAxis = (nodeBBox.m_max.m_y + nodeBBox.m_min.m_y) * 0.5f;
    else
        splitAxis = (nodeBBox.m_max.m_z + nodeBBox.m_min.m_z) * 0.5f;
    
    m_nodes[nodeIndex].m_bbox = nodeBBox;
    m_nodes[nodeIndex].m_flags = split;
    
    // Separate primitives such that those on the left of the split are in the
    // earlier part of the list (for the range we're dealing with) and those on
    // the right part of the split are later part of the list.
    BuildElementPredicate pred(splitAxis, split);
    BuildElement* partitionIter = std::partition(&permutedElements[begin], (&permutedElements[0]) + end, pred);
    unsigned int splitIndex = (unsigned int)(partitionIter - (&permutedElements[0]));
    
    // Peel off half of the elements if one side of the partition was empty
    // Note: doing this makes *crappy* BVH nodes at this part of the tree, but
    // it keeps us from generating pathologically-stupid trees instead in some
    // difficult cases.  Better to be merely crappy than pathologically stupid.
    if (splitIndex <= begin || splitIndex >= end)
    {
        splitIndex = begin + (end - begin) / 2;
        if (splitIndex < begin + 1)
            splitIndex = begin + 1;
        else if (splitIndex > end - 1)
            splitIndex = end - 1;
    }
    
    // Calculate bboxes for each new child we're about to create
    BBox leftBBox, rightBBox;
    for (unsigned int i = begin; i < splitIndex; ++i)
    {
        leftBBox = leftBBox.combined(permutedElements[i].m_bbox);
    }
    for (unsigned int i = splitIndex; i < end; ++i)
    {
        rightBBox = rightBBox.combined(permutedElements[i].m_bbox);
    }
    
    // Create children nodes, recurse to keep building
    m_nodes[nodeIndex].m_firstChild = m_numNodes;
    m_numNodes += 2;
    if (!buildRange(permutedElements, begin, splitIndex, m_nodes[nodeIndex].m_firstChild, leftBBox))
        return false;
    if (!buildRange(permutedElements, splitIndex, end, m_nodes[nodeIndex].m_firstChild + 1, rightBBox))
        return false;
    
    return true;
}

// Arbitrary limit on tree depth; there can be 2^32 nodes, or 2^31 prims implying
// a max depth of 32, but the trees are not perfectly balanced, so we add some
// slack that hopefully will suffice.
const unsigned int kMaxTraversalSteps = 50;

// Temporary data used during traversal to remember a node we need to potentially
// still visit and examine for intersection.
struct TraversalStep
{
    unsigned int m_nodeIndex;
    float m_t0, m_t1;
};

template<typename T>
bool Bvh<T>::doesIntersect(const Ray& ray)
{
    // Ray-bbox intersection uses the inverse direction (for performance reasons)
    Vector invDir(1.0f / ray.m_direction);
    
    // In order to find which child is "closer" along the ray we have to know
    // which direction the ray is going relative to each BVH node's spliting axis
    bool dirSigns[3] =
    {
        invDir.m_x < 0.0f,
        invDir.m_y < 0.0f,
        invDir.m_z < 0.0f
    };
    
    // Maintain a list of nodes we need to examine, and the enter/exit distances
    // along the ray they live in.
    TraversalStep steps[kMaxTraversalSteps];
    // Start with the root node (if we have one)
    unsigned int numSteps = (m_nodes != NULL && m_numNodes > 0) ? 1 : 0;
    steps[0].m_nodeIndex = 0;
    steps[0].m_t0 = kRayTMin;
    steps[0].m_t1 = ray.m_tMax;

    // Process pending nodes until we run out
    while (numSteps > 0 && numSteps <= kMaxTraversalSteps)
    {
        unsigned int step = numSteps - 1;
        const BvhNode& node = m_nodes[steps[step].m_nodeIndex];
        
        // Test prim if this is a prim node
        if (node.leafNode())
        {
            if (m_object.doesIntersect(ray, node.m_prim))
            {
                return true;
            }
            // Done with this prim node
            numSteps--;
            continue;
        }
        
        // Test ray against node bbox, adjusting ranges back if possible based
        // on previous near intersections
        float t0 = steps[step].m_t0;
        float t1 = steps[step].m_t1;
        if (!node.m_bbox.intersects(ray.m_origin, invDir, t0, t1))
        {
            // Ray misses the bbox, skip the node
            numSteps--;
            continue;
        }
        
        // Find which child node is closest
        // NOTE: it's not unreasonable to skip this check and just pick
        // an order, since we only care if *something* intersected at all, but
        // we leave the ordering in the hopes that closer objects will be more
        // frequently hit along the ray.
        unsigned int closestNode, furthestNode;
        if (dirSigns[node.split()] == false)
        {
            furthestNode = node.leftChildIndex();
            closestNode = node.rightChildIndex();
        }
        else
        {
            closestNode = node.leftChildIndex();
            furthestNode = node.rightChildIndex();
        }
        
        // Replace current step with furthest child
        steps[step].m_nodeIndex = furthestNode;
        steps[step].m_t0 = t0;
        steps[step].m_t1 = t1;
        // Push closest child as the next step to evaluate
        numSteps++;
        step++;
        steps[step].m_nodeIndex = closestNode;
        steps[step].m_t0 = t0;
        steps[step].m_t1 = t1;
    }
    return false;
}

template<typename T>
bool Bvh<T>::intersect(Intersection& intersection)
{
    // Ray-bbox intersection uses the inverse direction (for performance reasons)
    Vector invDir(1.0f / intersection.m_ray.m_direction);
    
    // In order to find which child is "closer" along the ray we have to know
    // which direction the ray is going relative to each BVH node's spliting axis
    bool dirSigns[3] =
    {
        invDir.m_x < 0.0f,
        invDir.m_y < 0.0f,
        invDir.m_z < 0.0f
    };
    
    // Maintain a list of nodes we need to examine, and the enter/exit distances
    // along the ray they live in.  We use the enter/exit information as we go
    // to find out if a node to be examined goes out of range, since a nearer
    // intersection may have already been found.  It allows us to skip nodes
    // quickly as they get out of range.
    TraversalStep steps[kMaxTraversalSteps];
    // Start with the root node (if we have one)
    unsigned int numSteps = (m_nodes != NULL && m_numNodes > 0) ? 1 : 0;
    steps[0].m_nodeIndex = 0;
    steps[0].m_t0 = kRayTMin;
    steps[0].m_t1 = intersection.m_t;

    // Process pending nodes until we run out
    bool intersected = false;
    while (numSteps > 0 && numSteps <= kMaxTraversalSteps)
    {
        unsigned int step = numSteps - 1;
        const BvhNode& node = m_nodes[steps[step].m_nodeIndex];
        
        // Test prim if this is a prim node
        if (node.leafNode())
        {
            if (m_object.intersect(intersection, node.m_prim))
            {
                intersected = true;
            }
            // Done with this prim node
            numSteps--;
            continue;
        }
        
        // Test ray against node bbox, adjusting ranges back if possible based
        // on previous near intersections
        float t0 = steps[step].m_t0;
        float t1 = steps[step].m_t1;
        if (t0 >= intersection.m_t)
        {
            // Previous near intersection was closer than this entire node, skip it
            numSteps--;
            continue;
        }
        if (t1 > intersection.m_t)
            t1 = intersection.m_t;
        if (!node.m_bbox.intersects(intersection.m_ray.m_origin, invDir, t0, t1))
        {
            // Ray misses the bbox, skip the node
            numSteps--;
            continue;
        }
        
        // Find which child node is closest
        unsigned int closestNode, furthestNode;
        if (dirSigns[node.split()] == false)
        {
            furthestNode = node.leftChildIndex();
            closestNode = node.rightChildIndex();
        }
        else
        {
            closestNode = node.leftChildIndex();
            furthestNode = node.rightChildIndex();
        }
        
        // Replace current step with furthest child
        steps[step].m_nodeIndex = furthestNode;
        steps[step].m_t0 = t0;
        steps[step].m_t1 = t1;
        // Push closest child as the next step to evaluate
        numSteps++;
        step++;
        steps[step].m_nodeIndex = closestNode;
        steps[step].m_t0 = t0;
        steps[step].m_t1 = t1;
    }
    return intersected;
}


} // namespace Rayito


#endif // __RACCEL_H__
