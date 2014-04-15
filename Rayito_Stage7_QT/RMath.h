////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RMATH_H__
#define __RMATH_H__

#include <cmath>
#include <algorithm>
#include <vector>


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
//     -vector
//     vector + vector, vector += vector
//     vector - vector, vector -= vector
//     vector * float, vector *= float, float * vector
//     vector / float, vector /= float, float / vector
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
    float  normalize()        { float len = length(); if (len > 0) *this /= len; return len; }
    // Return a vector in this same direction, but normalized
    Vector normalized() const { Vector r(*this); r.normalize(); return r; }
    
    float maxComponent() const { return std::max(std::max(m_x, m_y), m_z); }
    float minComponent() const { return std::min(std::min(m_x, m_y), m_z); }
    
    
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
    
    Vector& operator *=(const Vector& v)
    {
        m_x *= v.m_x;
        m_y *= v.m_y;
        m_z *= v.m_z;
        return *this;
    }
    
    Vector& operator *=(float f)
    {
        m_x *= f;
        m_y *= f;
        m_z *= f;
        return *this;
    }
    
    Vector& operator /=(const Vector& v)
    {
        m_x /= v.m_x;
        m_y /= v.m_y;
        m_z /= v.m_z;
        return *this;
    }
    
    Vector& operator /=(float f)
    {
        m_x /= f;
        m_y /= f;
        m_z /= f;
        return *this;
    }
    
    Vector operator -() const
    {
        return Vector(-m_x, -m_y, -m_z);
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


inline Vector operator *(const Vector& v1, const Vector& v2)
{
    return Vector(v1.m_x * v2.m_x,
                  v1.m_y * v2.m_y,
                  v1.m_z * v2.m_z);
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


inline Vector operator /(const Vector& v1, const Vector& v2)
{
    return Vector(v1.m_x / v2.m_x,
                  v1.m_y / v2.m_y,
                  v1.m_z / v2.m_z);
}


inline Vector operator /(float f, const Vector& v)
{
    return Vector(f / v.m_x,
                  f / v.m_y,
                  f / v.m_z);
}


inline Vector operator /(const Vector& v, float f)
{
    return Vector(v.m_x / f,
                  v.m_y / f,
                  v.m_z / f);
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


inline Vector max(const Vector& v1, const Vector& v2)
{
    return Vector(std::max(v1.m_x, v2.m_x),
                  std::max(v1.m_y, v2.m_y),
                  std::max(v1.m_z, v2.m_z));
}

inline Vector min(const Vector& v1, const Vector& v2)
{
    return Vector(std::min(v1.m_x, v2.m_x),
                  std::min(v1.m_y, v2.m_y),
                  std::min(v1.m_z, v2.m_z));
}


// Oh, by the way, a point can be thought of as just a vector, but where you
// refrain from doing dot/cross/normalize operations on it.
typedef Vector Point;


//
// 3D quaternion (axis-angle rotation) class (and associated operations)
//
// Operators supported:
//     quat = quat
//     -quat
//     ~quat (the conjugate)
//     quat + quat, quat += quat
//     quat - quat, quat -= quat
//     quat * quat (composition, not commutative)
//     quat * float, quat *= float, float * quat
//     quat / quat (composition w/ inverse, not commutative)
//     quat / float, quat /= float, float / quat
//     quat * vector (rotate the vector or point)
//

struct Quaternion
{
    float m_w;
    Vector m_v;
    
    Quaternion()                                   : m_w(1.0f), m_v(0.0f)   { }
    Quaternion(const Quaternion& q)                : m_w(q.m_w), m_v(q.m_v) { }
    Quaternion(float w, float x, float y, float z) : m_w(w), m_v(x, y, z)   { }
    Quaternion(float w, const Vector& v)           : m_w(w), m_v(v)         { }
    
    // From axis and angle
    Quaternion(const Vector& axis, float angle)
        : m_w(std::cos(angle * 0.5f)), m_v(axis * std::sin(angle * 0.5f))   { }
    
    // From three Euler angles, assuming order is ZYX
    Quaternion(float ex, float ey, float ez)
    {
        float cosX2 = std::cos(ex * 0.5f);
        float cosY2 = std::cos(ey * 0.5f);
        float cosZ2 = std::cos(ez * 0.5f);
        
        float sinX2 = std::sin(ex * 0.5f);
        float sinY2 = std::sin(ey * 0.5f);
        float sinZ2 = std::sin(ez * 0.5f);
        
        m_w     = cosZ2 * cosY2 * cosX2 + sinZ2 * sinY2 * sinX2;
        m_v.m_x = cosZ2 * cosY2 * sinX2 - sinZ2 * sinY2 * cosX2;
        m_v.m_y = cosZ2 * sinY2 * cosX2 + sinZ2 * cosY2 * sinX2;
        m_v.m_z = sinZ2 * cosY2 * cosX2 - cosZ2 * sinY2 * sinX2;
    }
    
   void toAxisAngle(float& outAngle, Vector& outAxis) const
    {
        outAngle = std::acos(m_w);
        float invSinAngle = 1.0f / std::sin(outAngle);
        outAxis = m_v * invSinAngle;
        outAngle *= 2.0f;
    }
    
    
    float length2() const { return m_w * m_w + m_v.length2(); }
    float length()  const { return std::sqrt(length2()); }
    
    // Returns old length from before normalization (ignore the return value if you don't need it)
    float      normalize()        { float len = length(); if (len > 0) *this /= len; return len; }
    // Return a vector in this same direction, but normalized
    Quaternion normalized() const { Quaternion q(*this); q.normalize(); return q; }
    
    // q * q.inverse() = q_identity
    Quaternion inverse() const
    {
        float len2 = length2();
        return Quaternion(m_w / len2, m_v / len2);
    }
    
    
    Quaternion& operator =(const Quaternion& q)
    {
        m_w = q.m_w;
        m_v = q.m_v;
        return *this;
    }
    
    Quaternion& operator +=(const Quaternion& q)
    {
        m_w += q.m_w;
        m_v += q.m_v;
        return *this;
    }
    
    Quaternion& operator -=(const Quaternion& q)
    {
        m_w -= q.m_w;
        m_v -= q.m_v;
        return *this;
    }
    
    Quaternion& operator *=(const Quaternion& q)
    {
        // Expansion of terms in: q1*q2 = (w1*w2 - dot(v1,v2), w1*v2 + w2*v1 + cross(v1, v2))
        m_w     = m_w * q.m_w     - m_v.m_x * q.m_v.m_x - m_v.m_y * q.m_v.m_y - m_v.m_z * q.m_v.m_z;
        m_v.m_x = m_w * q.m_v.m_x + m_v.m_x * q.m_w     + m_v.m_y * q.m_v.m_z - m_v.m_z * q.m_v.m_y;
        m_v.m_y = m_w * q.m_v.m_y - m_v.m_x * q.m_v.m_z + m_v.m_y * q.m_w     + m_v.m_z * q.m_v.m_x;
        m_v.m_z = m_w * q.m_v.m_z + m_v.m_x * q.m_v.m_y - m_v.m_y * q.m_v.m_x + m_v.m_z * q.m_w;
        return *this;
    }
    
    Quaternion& operator *=(float f)
    {
        m_w *= f;
        m_v *= f;
        return *this;
    }
    
    Quaternion& operator /=(const Quaternion& q)
    {
        *this *= q.inverse();
        return *this;
    }
    
    Quaternion& operator /=(float f)
    {
        m_w /= f;
        m_v /= f;
        return *this;
    }
    
    Quaternion operator -() const
    {
        return Quaternion(-m_w, -m_v);
    }
    
    Quaternion operator ~() const
    {
        return Quaternion(m_w, -m_v);
    }
};


inline Quaternion operator +(const Quaternion& q1, const Quaternion& q2)
{
    return Quaternion(q1.m_w + q2.m_w, q1.m_v + q2.m_v);
}


inline Quaternion operator -(const Quaternion& q1, const Quaternion& q2)
{
    return Quaternion(q1.m_w - q2.m_w, q1.m_v - q2.m_v);
}


inline Quaternion operator *(const Quaternion& q1, const Quaternion& q2)
{
    // Expansion of terms in: q1*q2 = (w1*w2 - dot(v1,v2), w1*v2 + w2*v1 + cross(v1, v2))
    return Quaternion(q1.m_w * q2.m_w     - q1.m_v.m_x * q2.m_v.m_x - q1.m_v.m_y * q2.m_v.m_y - q1.m_v.m_z * q2.m_v.m_z,
                      q1.m_w * q2.m_v.m_x + q1.m_v.m_x * q2.m_w     + q1.m_v.m_y * q2.m_v.m_z - q1.m_v.m_z * q2.m_v.m_y,
                      q1.m_w * q2.m_v.m_y - q1.m_v.m_x * q2.m_v.m_z + q1.m_v.m_y * q2.m_w     + q1.m_v.m_z * q2.m_v.m_x,
                      q1.m_w * q2.m_v.m_z + q1.m_v.m_x * q2.m_v.m_y - q1.m_v.m_y * q2.m_v.m_x + q1.m_v.m_z * q2.m_w);
}


inline Quaternion operator *(const Quaternion& q, float f)
{
    return Quaternion(f * q.m_w, f * q.m_v);
}


inline Quaternion operator *(float f, const Quaternion& q)
{
    return Quaternion(f * q.m_w, f * q.m_v);
}

inline Vector operator *(const Quaternion& q, const Vector& v)
{
    // Canonical quaternion-vector multiply:
    //     q * v = q * quat(0, v.x, v.y, v.z) * ~q
    // In matrix form:
    //     [ 1-2(qy*qy+qz*qz)  2(qx*qy-qw*qz)  2(qx*qz+qw*qy)  ][ vx ]
    //     [  2(qx*qy+qw*qz)  1-2(qx*qx+qz*qz) 2(qy*qz-qw*qx)  ][ vy ]
    //     [  2(qx*qz-qw*qy)   2(qy*qz+qw*qx) 1-2(qx*qx+qy*qy) ][ vz ]
    // But can be optimized to:
    //     t = 2*cross(qv, v)
    //     result = v + qw*t + cross(qv, t)
    Vector t = 2.0f * cross(q.m_v, v);
    return v + t * q.m_w + cross(q.m_v, t);
}

inline Quaternion operator /(const Quaternion& q1, const Quaternion& q2)
{
    return q1 * q2.inverse();
}


inline Quaternion operator /(float f, const Quaternion& q)
{
    return f * q.inverse();
}


inline Quaternion operator /(const Quaternion& q, float f)
{
    return q * (1.0f / f);
}


inline float dot(const Quaternion& q1, const Quaternion& q2)
{
    // In cartesian coordinates, it simplifies to this simple calculation:
    return q1.m_w * q2.m_w + dot(q1.m_v, q2.m_v);
}

// Linearly interpolate quaternions; assumes they are normalized
inline Quaternion lerp(const Quaternion& q1, const Quaternion& q2, float t)
{
    
    return (q1 * (1.0f - t) + q2 * t).normalized();
}

// Spherical linear interpolation of quaternions (slower, but higher quality),
// assumes they are normalized
inline Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t)
{
    Quaternion q2Adjusted;
    float dotq = dot(q1, q2);
    // Adjust to find the shortest arc between q1 and q2
    if (dotq < 0.0f)
    {
        q2Adjusted = -q2;
        dotq = -dotq;
    }
    else
    {
        q2Adjusted = q2;
    }
    if (dotq < 0.95f)
    {
        // High quality rotation
        float angle = std::acos(dotq);
        return (q1 * std::sin(angle * (1.0f - t)) + q2Adjusted * std::sin(angle * t)) / std::sin(angle);
    }
    else
    {
        // Rotation angle is small, linear interp will be just as good
        return lerp(q1, q2, t);
    }
}


//
// Transformations
//


// Transformation class (*not* a matrix, but instead for simplicity of motion blur
// it encodes a scale, rotate, and then translate (applied in that order))
class Transform
{
public:
    Transform() : m_time(), m_scale(), m_rotate(), m_translate() { }
    Transform(const Transform& t)
                : m_time(t.m_time), m_scale(t.m_scale), m_rotate(t.m_rotate), m_translate(t.m_translate) { }
    
    Transform& operator =(const Transform& t)
    {
        m_time = t.m_time;
        m_scale = t.m_scale;
        m_rotate = t.m_rotate;
        m_translate = t.m_translate;
        return *this;
    }
    
    
    // There is always at least ONE key, the default (identity transform) case
    size_t numKeys()     const { return m_time.empty() ? 1 : m_time.size(); }
    size_t numSegments() const { return m_time.size() < 1 ? 0 : m_time.size() - 1; }
    
    float keyTime(size_t keyIndex) const { return keyIndex < m_time.size() ? m_time[keyIndex] : 0.0f; }
    
    
    void clear()
    {
        m_time.clear();
        m_scale.clear();
        m_rotate.clear();
        m_translate.clear();
    }
    

    // Current accumulated aspects of the transformation
    
    Vector translationKey(size_t keyIndex) const
    {
        if (m_translate.empty())
            return Vector(0.0f);
        if (keyIndex >= m_translate.size())
            keyIndex = m_translate.size() - 1;
        return m_translate[keyIndex];
    }
    
    Vector scalingKey(size_t keyIndex) const
    {
        if (m_scale.empty())
            return Vector(1.0f);
        if (keyIndex >= m_scale.size())
            keyIndex = m_scale.size() - 1;
        return m_scale[keyIndex];
    }
    
    Quaternion rotationKey(size_t keyIndex) const
    {
        if (m_rotate.empty())
            return Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
        if (keyIndex >= m_rotate.size())
            keyIndex = m_rotate.size() - 1;
        return m_rotate[keyIndex];
    }
    
    Vector translation(float time) const
    {
        if (m_time.empty())
            return Vector(0.0f);
        float t;
        size_t index = timeIndex(time, t);
        if (t == 0.0f)
            return m_translate[index];
        else
            return m_translate[index] * (1.0f - t) + m_translate[index + 1] * t;
    }
    
    Vector scaling(float time) const
    {
        if (m_time.empty())
            return Vector(1.0f);
        float t;
        size_t index = timeIndex(time, t);
        if (t == 0.0f)
            return m_scale[index];
        else
            return m_scale[index] * (1.0f - t) + m_scale[index + 1] * t;
    }
    
    Quaternion rotation(float time) const
    {
        if (m_time.empty())
            return Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
        float t;
        size_t index = timeIndex(time, t);
        if (t == 0.0f)
            return m_rotate[index];
        else
            return lerp(m_rotate[index], m_rotate[index + 1], t);
    }
    
    // Overwrite aspects of the transformation
    
    void setTranslationKey(size_t keyIndex, const Vector& trans)
    {
        if (keyIndex >= m_translate.size())
            return;
        m_translate[keyIndex] = trans;
    }
    
    void setScalingKey(size_t keyIndex, const Vector& scaling)
    {
        if (keyIndex >= m_scale.size())
            return;
        m_scale[keyIndex] = scaling;
    }
    
    void setRotationKey(size_t keyIndex, const Quaternion& rot)
    {
        if (keyIndex >= m_rotate.size())
            return;
        m_rotate[keyIndex] = rot;
    }
    
    void setTranslation(float time, const Vector& trans)
    {
        size_t index = findOrInsertKey(time);
        m_translate[index] = trans;
    }
    
    void setScaling(float time, const Vector& scaling)
    {
        size_t index = findOrInsertKey(time);
        m_scale[index] = scaling;
    }
    
    void setRotation(float time, const Quaternion& rot)
    {
        size_t index = findOrInsertKey(time);
        m_rotate[index] = rot;
    }
    
    // Concatenate additional transformations
    
    void translateKey(size_t keyIndex, const Vector& trans)
    {
        if (keyIndex >= m_translate.size())
            return;
        m_translate[keyIndex] += trans;
    }
    
    void scaleKey(size_t keyIndex, const Vector& scaling)
    {
        if (keyIndex >= m_scale.size())
            return;
        m_scale[keyIndex] *= scaling;
    }
    
    void rotateKey(size_t keyIndex, const Quaternion& rot)
    {
        if (keyIndex >= m_rotate.size())
            return;
        m_rotate[keyIndex] *= rot;
    }
    
    void translate(float time, const Vector& trans)
    {
        size_t index = findOrInsertKey(time);
        m_translate[index] += trans;
    }
    
    void scale(float time, const Vector& scaling)
    {
        size_t index = findOrInsertKey(time);
        m_scale[index] *= scaling;
    }
    
    void rotate(float time, const Quaternion& rot)
    {
        size_t index = findOrInsertKey(time);
        m_rotate[index] *= rot;
    }
    
    
    void prepare()
    {
        // For the rotations to behave properly, they need to be normalized
        for (size_t i = 0; i < m_rotate.size(); ++i)
        {
            m_rotate[i].normalize();
        }
    }
    
    
    // Take various quantities into and out of the space defined by this transformation.
    // From local: apply scale, then rotation, then translation.
    // To local: back out translation, unapply rotation, and unapply scale.
    
    Point toLocalPoint(float time, const Point& p) const
    {
        return ((~rotation(time)) * (p - translation(time))) / scaling(time);
    }
    
    Point fromLocalPoint(float time, const Point& p) const
    {
        return rotation(time) * (p * scaling(time)) + translation(time);
    }
    
    Vector toLocalVector(float time, const Vector& v) const
    {
        return ((~rotation(time)) * v) / scaling(time);
    }
    
    Vector fromLocalVector(float time, const Vector& v) const
    {
        return rotation(time) * (v * scaling(time));
    }
    
    Vector toLocalNormal(float time, const Vector& n) const
    {
        return (~rotation(time)) * n;
    }
    
    Vector fromLocalNormal(float time, const Vector& n) const
    {
        return rotation(time) * n;
    }
    
private:
    std::vector<float>      m_time;
    std::vector<Vector>     m_scale;
    std::vector<Quaternion> m_rotate;
    std::vector<Vector>     m_translate;
    
    size_t timeIndex(float time, float& outT) const
    {
        // Find the index in m_time that is just before the specified time.
        // If the time is in between two times, we output the proportion of the
        // time at the index vs the time at index + 1 so values can be mixed
        // between the two indexes.
        
        // Times are strictly increasing, so we can do a quick little binary
        // search here to find the keys bracketing this time
        size_t lower = 0;
        size_t upper = m_time.size() - 1;
        if (m_time[upper] <= time)
            lower = upper;
        else if (m_time[lower] >= time)
            upper = lower;
        while (upper - lower > 0)
        {
            size_t mid = (lower + upper) / 2;
            if (time < m_time[mid])
                upper = mid;
            else if (mid > lower)
                lower = mid;
            else
                break;
        }
        size_t index = lower;
        // Compute the 0..1 mix of the keys surrounding the time
        if (index == m_time.size() - 1)
            outT = 0.0f;
        else if (m_time[index] >= time)
            outT = 0.0f; // Peg it to the beginning of the range
        else
            outT = (time - m_time[index]) / (m_time[index + 1] - m_time[index]);
        return index;
    }
    
    size_t findOrInsertKey(float time)
    {
        // Here's the deal.  If a key at the time slot exists, we just return it.
        // If it's empty or the time is before or after existing time keys, we
        // insert the time for all parts of the transform.  If the time is
        // between two existing keys, we insert an interpolated key.
        size_t index;
        if (m_time.empty())
        {
            // No times yet, make a key
            m_translate.push_back(Vector(0.0f));
            m_scale.push_back(Vector(1.0f));
            m_rotate.push_back(Quaternion(1.0f, 0.0f, 0.0f, 0.0f));
            m_time.push_back(time);
            index = 0;
        }
        else if (time > m_time.back())
        {
            // Time is past the end, add a new key at the end
            m_translate.push_back(m_translate.back());
            m_scale.push_back(m_scale.back());
            m_rotate.push_back(m_rotate.back());
            m_time.push_back(time);
            index = m_time.size() - 1;
        }
        else if (time < m_time[0])
        {
            // Time is before the first, insert a new key at the beginning
            m_translate.insert(m_translate.begin(), m_translate.front());
            m_scale.insert(m_scale.begin(), m_scale.front());
            m_rotate.insert(m_rotate.begin(), m_rotate.front());
            m_time.insert(m_time.begin(), time);
            index = 0;
        }
        else
        {
            // Find the index of time just before the time specified
            float t;
            index = timeIndex(time, t);
            if (t != 0.0f && t != 1.0f && index < m_time.size() - 1)
            {
                // Time is in between two keys, insert one in between them,
                // mixing the two values.
                index++;
                m_translate.insert(m_translate.begin() + index,
                                   m_translate[index - 1] * (1.0f - t) + m_translate[index] * t);
                m_scale.insert(m_scale.begin() + index,
                               m_scale[index - 1] * (1.0f - t) + m_scale[index] * t);
                m_rotate.insert(m_rotate.begin() + index,
                                lerp(m_rotate[index - 1], m_rotate[index], t));
                m_time.insert(m_time.begin() + index, time);
            }
        }
        return index;
    }
};


// Given a single direction, generate a coordinate space where the direction
// becomes the Z axis, and the X and Y axes are consistently made to match
inline void makeCoordinateSpace(const Vector& normalRef,
                                Vector& outXAxis, Vector& outYAxis, Vector& outZAxis)
{
    outZAxis = normalRef.normalized();
    Vector v2 = (outZAxis.m_x != 0.0f || outZAxis.m_z != 0.0f) ?
                    Vector(0.0f, 1.0f, 0.0f) :
                    Vector(1.0f, 0.0f, 0.0f);
    outXAxis = cross(v2, outZAxis).normalized();
    outYAxis = cross(outZAxis, outXAxis);
}

// Given a two directions, generate a coordinate space where the normal becomes
// the Z axis, and the X and Y axes are consistently made to match and align the
// X axis as much as possible with the tangent direction
inline void makeCoordinateSpace(const Vector& normalRef, const Vector& tangent,
                                Vector& outXAxis, Vector& outYAxis, Vector& outZAxis)
{
    outZAxis = normalRef.normalized();
    outYAxis = cross(tangent, outZAxis).normalized();
    outXAxis = cross(outZAxis, outYAxis);
}

// Transform a vector into the local space defined by X,Y,Z orthonormal axes
inline Vector transformToLocalCoordinateSpace(const Vector& v,
                                              const Vector& xAxis,
                                              const Vector& yAxis,
                                              const Vector& zAxis)
{
    return Vector(dot(v, xAxis), dot(v, yAxis), dot(v, zAxis));
}

// Transform a vector out of the local space defined by X,Y,Z orthonormal axes
inline Vector transformFromLocalCoordinateSpace(const Vector& v,
                                                const Vector& xAxis,
                                                const Vector& yAxis,
                                                const Vector& zAxis)
{
    return Vector(v.m_x * xAxis.m_x + v.m_y * yAxis.m_x + v.m_z * zAxis.m_x,
                  v.m_x * xAxis.m_y + v.m_y * yAxis.m_y + v.m_z * zAxis.m_y,
                  v.m_x * xAxis.m_z + v.m_y * yAxis.m_z + v.m_z * zAxis.m_z);
}


} // namespace Rayito


#endif // __RMATH_H__
