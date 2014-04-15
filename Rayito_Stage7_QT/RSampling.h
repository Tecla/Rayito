////////////////////////////////////////////////////////////////////////////////
//
// Very simple ray tracing example
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __RSAMPLING_H__
#define __RSAMPLING_H__

#include <cmath>
#include <algorithm>

#include "RMath.h"


namespace Rayito
{


//
// Random number and 1D/2D sample pattern generators
//

// Marsaglia multiply-with-carry psuedo random number generator.  It's very fast
// and has good distribution properties.  Has a period of 2^60. See
// http://groups.google.com/group/sci.crypt/browse_thread/thread/ca8682a4658a124d/
struct Rng
{
    unsigned int m_z, m_w;
    
    Rng(unsigned int z = 362436069, unsigned int w = 521288629)
        : m_z(z), m_w(w) { }
    
    Rng(const Rng& rng) : m_z(rng.m_z), m_w(rng.m_w) { }
    
    Rng& operator =(const Rng& rng)
    {
        m_z = rng.m_z;
        m_w = rng.m_w;
        return *this;
    }
    
    
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
        return (m_z << 16) + m_w;  // 32-bit result
    }
};


const unsigned int kUnlimitedSamples = 0;


class Sampler
{
public:
    Sampler(Rng& rng) : m_rng(rng), m_currentSampleIndex(0) { }
    
    virtual ~Sampler() { }
    
    
    float nextSample1D()
    {
        if (m_currentSampleIndex < total1DSamplesAvailable())
        {
            float result = sample1D(m_currentSampleIndex);
            ++m_currentSampleIndex;
            return result;
        }
        return m_rng.nextFloat();
    }

    virtual float        sample1D(unsigned int index)    = 0;
    virtual unsigned int total1DSamplesAvailable() const = 0;
    
    
    void nextSample2D(float& outD1, float& outD2)
    {
        if (m_currentSampleIndex < total1DSamplesAvailable())
        {
            sample2D(m_currentSampleIndex, outD1, outD2);
            ++m_currentSampleIndex;
        }
        else
        {
            outD1 = m_rng.nextFloat();
            outD2 = m_rng.nextFloat();
        }
    }

    virtual void         sample2D(unsigned int index, float& outD1, float& outD2) = 0;
    virtual unsigned int total2DSamplesAvailable() const                          = 0;
    
    virtual void refill(unsigned int permutation) = 0;
    
protected:
    Rng& m_rng;
    unsigned int m_currentSampleIndex;
};


class RandomSampler : public Sampler
{
public:
    RandomSampler(unsigned int xSamples, unsigned int ySamples, Rng& rng)
        : Sampler(rng), m_xSamples(xSamples), m_ySamples(ySamples) { }
    
    RandomSampler(unsigned int samples, Rng& rng)
        : Sampler(rng), m_xSamples(samples), m_ySamples(0) { }
    
    virtual ~RandomSampler() { }
    
    
    virtual float sample1D(unsigned int index)
    {
        return m_rng.nextFloat();
    }

    virtual unsigned int total1DSamplesAvailable() const
    {
        return m_xSamples;
    }

    virtual void sample2D(unsigned int index, float& outD1, float& outD2)
    {
        outD1 = m_rng.nextFloat();
        outD2 = m_rng.nextFloat();
    }

    virtual unsigned int total2DSamplesAvailable() const
    {
        return m_xSamples * m_ySamples;
    }
    
    virtual void refill(unsigned int) { m_currentSampleIndex = 0; }
    
protected:
    unsigned int m_xSamples, m_ySamples;
};


class StratifiedRandomSampler : public Sampler
{
public:
    StratifiedRandomSampler(unsigned int xSamples,
                            unsigned int ySamples,
                            Rng& rng,
                            unsigned int permutation)
        : Sampler(rng), m_permutation(permutation), m_xSamples(xSamples),
          m_ySamples(ySamples), m_samples(new float[xSamples * ySamples * 2]),
          m_is2D(true) { generateSamples(); }
    
    StratifiedRandomSampler(unsigned int samples, Rng& rng, unsigned int permutation)
        : Sampler(rng), m_permutation(permutation), m_xSamples(samples),
          m_ySamples(0), m_samples(new float[samples]), m_is2D(false) { generateSamples(); }
    
    virtual ~StratifiedRandomSampler() { delete[] m_samples; }
    
    
    virtual float sample1D(unsigned int index)
    {
        if (m_is2D || index >= m_xSamples)
            return m_rng.nextFloat(); // Bad, bad programmer!  Using a 1D sample from a 2D pattern?
        // Not the best, but good enough: rotate the samples per permutation
        unsigned int permutedIndex = (index + m_permutation) % m_xSamples;
        return m_samples[permutedIndex];
    }

    virtual unsigned int total1DSamplesAvailable() const
    {
        if (m_is2D)
            return 0;
        return m_xSamples;
    }

    virtual void sample2D(unsigned int index, float& outD1, float& outD2)
    {
        if (!m_is2D || index >= m_xSamples * m_ySamples)
        {
            // Bad, bad programmer!  Using a 2D sample from a 1D pattern?
            outD1 = m_rng.nextFloat(); 
            outD2 = m_rng.nextFloat();
        }
        else
        {
            unsigned int permutedIndex = (index + m_permutation) % (m_xSamples * m_ySamples);
            outD1 = m_samples[permutedIndex * 2];
            outD2 = m_samples[permutedIndex * 2 + 1];
        }
    }

    virtual unsigned int total2DSamplesAvailable() const
    {
        if (!m_is2D)
            return 0;
        return m_xSamples * m_ySamples;
    }
    
    virtual void refill(unsigned int permutation)
    {
        m_permutation = permutation;
        generateSamples();
        m_currentSampleIndex = 0;
    }
    
protected:
    unsigned int m_permutation, m_xSamples, m_ySamples;
    float *m_samples;
    bool m_is2D;
    
    void generateSamples()
    {
        if (m_is2D)
        {
            for (unsigned int i = 0; i < m_ySamples; ++i)
            {
                for (unsigned int j = 0; j < m_xSamples; ++j)
                {
                    float xPos = (j + m_rng.nextFloat()) / float(m_xSamples);
                    float yPos = (i + m_rng.nextFloat()) / float(m_ySamples);
                    m_samples[(i * m_xSamples + j) * 2] = xPos;
                    m_samples[(i * m_xSamples + j) * 2 + 1] = yPos;
                }
            }
        }
        else
        {
            for (unsigned int i = 0; i < m_xSamples; ++i)
            {
                float pos = (i + m_rng.nextFloat()) / float(m_xSamples);
                m_samples[i] = pos;
            }
        }
    }
};


// Correlated multi-jitter sampler, which has some very nice properties.  It has
// less clumping than just stratified sampling, and doesn't require any storage
// of samples due to how it computes its permutations.
//
// See "Correlated Multi-Jittered Sampling", Andrew Kensler, Pixar Tech Memo 13-01.
class CorrelatedMultiJitterSampler : public Sampler
{
public:
    CorrelatedMultiJitterSampler(unsigned int xSamples,
                                 unsigned int ySamples,
                                 Rng& rng,
                                 unsigned int permutation)
        : Sampler(rng), m_permutation(permutation), m_xSamples(xSamples),
          m_ySamples(ySamples), m_is2D(true) { }
    
    CorrelatedMultiJitterSampler(unsigned int samples,
                                 Rng& rng,
                                 unsigned int permutation)
        : Sampler(rng), m_permutation(permutation), m_xSamples(samples),
          m_ySamples(0), m_is2D(false) { }
    
    virtual ~CorrelatedMultiJitterSampler() { }
    
    
    virtual float sample1D(unsigned int index)
    {
        if (m_is2D || index >= m_xSamples)
            return m_rng.nextFloat(); // Bad, bad programmer!  Using a 1D sample from a 2D pattern?
        unsigned int permutedIndex = permute(index, m_xSamples, m_permutation * 0x8ff3cd11);
        float sx = randFloat01(permutedIndex, m_permutation * 0xa399d265);
        return (permutedIndex + sx) / float(m_xSamples);
    }

    virtual unsigned int total1DSamplesAvailable() const
    {
        if (m_is2D)
            return 0;
        return m_xSamples;
    }

    virtual void sample2D(unsigned int index, float& outD1, float& outD2)
    {
        if (!m_is2D || index >= m_xSamples * m_ySamples)
        {
            // Bad, bad programmer!  Using a 2D sample from a 1D pattern?
            outD1 = m_rng.nextFloat(); 
            outD2 = m_rng.nextFloat();
        }
        else
        {
            unsigned int permutedIndex = permute(index, m_xSamples * m_ySamples, m_permutation * 0xc2d3c8fb);
            int ix = permute(permutedIndex % m_xSamples, m_xSamples, m_permutation * 0xa511e9b3);
            int iy = permute(permutedIndex / m_xSamples, m_ySamples, m_permutation * 0x63d83595);
            float sx = randFloat01(permutedIndex, m_permutation * 0xa399d265);
            float sy = randFloat01(permutedIndex, m_permutation * 0x711ad6a5);
            outD1 = (ix + (iy + sx) / float(m_ySamples)) / float(m_xSamples);
            outD2 = (permutedIndex + sy) / float(m_xSamples * m_ySamples);
        }
    }

    virtual unsigned int total2DSamplesAvailable() const
    {
        if (!m_is2D)
            return 0;
        return m_xSamples * m_ySamples;
    }
    
    virtual void refill(unsigned int permutation)
    {
        m_permutation = permutation;
        m_currentSampleIndex = 0;
    }
    
protected:
    unsigned int m_permutation, m_xSamples, m_ySamples;
    bool m_is2D;
    
    // Permute a value in pseudo-random fashion via hashing.  Achieves avalanche
    // and works on values up to 2^27.  Uses "cycle-walking" when len is not a
    // power-of-two -- this will have a cost of up to 2 loop iterations.
    unsigned int permute(unsigned int i, unsigned int num, unsigned int permutation)
    {
        unsigned int w = num - 1;
        w |= w >> 1;
        w |= w >> 2;
        w |= w >> 4;
        w |= w >> 8;
        w |= w >> 16;
        do
        {
            i ^= permutation;
            i *= 0xe170893d;
            i ^= permutation >> 16;
            i ^= (i & w) >> 4;
            i ^= permutation >> 8;
            i *= 0x0929eb3f;
            i ^= permutation >> 23;
            i ^= (i & w) >> 1;
            i *= 1 | permutation >> 27;
            i *= 0x6935fa69;
            i ^= (i & w) >> 11;
            i *= 0x74dcb303;
            i ^= (i & w) >> 2;
            i *= 0x9e501cc3;
            i ^= (i & w) >> 2;
            i *= 0xc860a3df;
            i &= w;
            i ^= i >> 5;
        } while (i >= num);
        return (i + permutation) % num;
    }
    
    // Generated a canonical pseudo-random float via hashing.  Achieves avalanche.
    static float randFloat01(unsigned int i, unsigned int permutation)
    {
        i ^= permutation;
        i ^= i >> 17;
        i ^= i >> 10;
        i *= 0xb36534e5;
        i ^= i >> 12;
        i ^= i >> 21;
        i *= 0x93fc4795;
        i ^= 0xdf6e307f;
        i ^= i >> 17;
        i *= 1 | permutation >> 18;
        return i * 2.328306e-10f;
    }
};


//
// Multiple importance sampling weightings
//

inline float balanceHeuristic(unsigned int numSamples1, float pdf1, unsigned int numSamples2, float pdf2)
{
    return numSamples1 * pdf1 / (numSamples1 * pdf1 + numSamples2 * pdf2);
}

inline float powerHeuristic(unsigned int numSamples1, float pdf1, unsigned int numSamples2, float pdf2)
{
    float weighted1 = numSamples1 * pdf1;
    float weighted2 = numSamples2 * pdf2;
    return weighted1 * weighted1 / (weighted1 * weighted1 + weighted2 * weighted2);
}


//
// Sample space transformations
//


inline void concentricSampleDisk(float u1, float u2, float& outDx, float& outDy)
{
    float r, theta;
    // Map uniform random numbers to [-1, 1]
    float sx = 2.0f * u1 - 1.0f;
    float sy = 2.0f * u2 - 1.0f;

    // Map square to $(r,\theta)$

    // Handle degeneracy at the origin
    if (sx == 0.0f && sy == 0.0f)
    {
        outDx = 0.0f;
        outDy = 0.0f;
        return;
    }
    if (sx >= -sy)
    {
        if (sx > sy)
        {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0f)
                theta = sy / r;
            else
                theta = 8.0f + sy / r;
        }
        else
        {
            // Handle second region of disk
            r = sy;
            theta = 2.0f - sx / r;
        }
    }
    else
    {
        if (sx <= sy)
        {
            // Handle third region of disk
            r = -sx;
            theta = 4.0f - sy / r;
        }
        else
        {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0f + sx / r;
        }
    }
    theta *= M_PI / 4.0f;
    outDx = r * std::cos(theta);
    outDy = r * std::sin(theta);
}

    
// Method for finding a random point on the sphere
inline Vector uniformToSphere(float u1, float u2)
{
    // Find a height uniformly distributed on the sphere
    float z = 1.0f - 2.0f * u1;
    // Find the radius based on that height that sits on the sphere surface
    float radius = std::sqrt(std::max(0.0f, 1.0f - z * z));
    // Find a random angle around the sphere's equator
    float phi = M_PI * 2.0f * u2;
    // And put it all together...
    return Vector(radius * std::cos(phi), radius * std::sin(phi), z);
}


// Method for finding a random point on a disk
inline void uniformToUniformDisk(float u1, float u2, float& uDx, float& uDy)
{
    // We are mapping from a unit square (0.0 to 1.0 in each dimension)
    // to a disk, accounting for the fact that samples will bunch up
    // in the center if we don't spread them out.  We do this by taking
    // the square root of the first sample, which drives the values
    // towards a radius of 1.0 and yields uniform samples (area-wise)
    // on the disk.  The other uniform parameter allows us to evenly
    // distribute around the disk angle-wise, so we just spread it out
    // over 2*Pi radians.
    float radius = std::sqrt(u1);
    float theta = M_PI * 2.0f * u2;
    // Convert back to x,y coordinates
    uDx = radius * std::cos(theta);
    uDy = radius * std::sin(theta);
}


// Evenly distribute random points over a hemisphere
inline Vector uniformToHemisphere(float u1, float u2)
{
    float radius = std::sqrt(std::max(0.0f, 1.0f - u1 * u1));
    float phi = M_PI * 2.0f * u2;
    return Vector(radius * std::cos(phi),
                  radius * std::sin(phi),
                  u1);
}


// Distribute points over a hemisphere, with more congregating at the pole
inline Vector uniformToCosineHemisphere(float u1, float u2)
{
    float diskX, diskY;
    // Could use uniformToUniformDisk() too
    //uniformToUniformDisk(u1, u2, diskX, diskY);
    concentricSampleDisk(u1, u2, diskX, diskY);
    float z = std::sqrt(std::max(0.0f, 1.0f - diskX * diskX - diskY * diskY));
    return Vector(diskX, diskY, z);
}


// Distribute points over a cone
inline Vector uniformToCone(float u1, float u2, float cosThetaMax)
{
    float cosTheta = u1 * (cosThetaMax - 1.0f) + 1.0f;
    float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = u2 * M_PI * 2.0f;
    return Vector(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
}

inline float uniformConePdf(float cosThetaMax)
{
    return cosThetaMax >= 1.0f ? 0.0 : 1.0f / (2.0f * M_PI * (1.0f - cosThetaMax));
}


// Random point on a triangle, converted to barycentric alpha and beta (gamma is 1 - alpha - beta)
inline void uniformToBarycentricTriangle(float u1, float u2, float& btu, float& btv)
{
    float u1_sqrt = std::sqrt(u1);
    btu = 1.0f - u1_sqrt;
    btv = u2 * u1_sqrt;
}


} // namespace Rayito


#endif // __RSAMPLING_H__
