#include <string>
#include <iostream>
#include <sstream>

#include "rayito.h"

#include <QThread>


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


namespace
{


//
// RenderThread works on a small chunk of the image
//
class RenderThread : public QThread
{
public:
    RenderThread(size_t xstart, size_t xend, size_t ystart, size_t yend,
                 Image *pImage,
                 ShapeSet& masterSet,
                 const Camera& cam,
                 std::list<Shape*>& lights,
                 size_t pixelSamplesHint, size_t lightSamplesHint,
                 size_t maxRayDepth)
        : m_xstart(xstart), m_xend(xend), m_ystart(ystart), m_yend(yend),
          m_pImage(pImage), m_masterSet(masterSet), m_camera(cam), m_lights(lights),
          m_pixelSamplesHint(pixelSamplesHint), m_lightSamplesHint(lightSamplesHint),
          m_maxRayDepth(maxRayDepth) { }
    
protected:
    virtual void run()
    {
        // Random number generator (for random pixel positions, light positions, etc)
        // We seed the generator for this render thread based on something that
        // doesn't change, but gives us a good variable seed for each thread.
        Rng rng(static_cast<unsigned int>(((m_xstart << 16) | m_xend) ^ m_xstart),
                static_cast<unsigned int>(((m_ystart << 16) | m_yend) ^ m_ystart));
        
        // The aspect ratio is used to make the image only get more zoomed in when
        // the height changes (and not the width)
        float aspectRatioXToY = float(m_pImage->width()) / float(m_pImage->height());
        
        // Set up samplers for each of the ray bounces.  Each bounce will use
        // the same sampler for all pixel samples in the pixel to reduce noise.
        Sampler **bounceSamplers = new Sampler*[m_maxRayDepth];
        for (size_t i = 0; i < m_maxRayDepth; ++i)
        {
            bounceSamplers[i] = new StratifiedRandomSampler(m_pixelSamplesHint,
                                                            m_pixelSamplesHint,
                                                            rng);
        }
        // Set up samplers for each pixel sample
        StratifiedRandomSampler lensSampler(m_pixelSamplesHint, m_pixelSamplesHint, rng);
        StratifiedRandomSampler sampler(m_pixelSamplesHint, m_pixelSamplesHint, rng);
        size_t totalPixelSamples = sampler.total2DSamplesAvailable();

        // For each pixel row...
        for (size_t y = m_ystart; y < m_yend; ++y)
        {
            // For each pixel across the row...
            for (size_t x = m_xstart; x < m_xend; ++x)
            {
                // Accumulate pixel color
                Color pixelColor(0.0f, 0.0f, 0.0f);
                // For each sample in the pixel...
                for (size_t psi = 0; psi < totalPixelSamples; ++psi)
                {
                    // Calculate a stratified random position within the pixel
                    // to hide aliasing
                    float pu, pv;
                    sampler.sample2D(psi, pu, pv);
                    float xu = (x + pu) / float(m_pImage->width());
                    // Flip pixel row to be in screen space (images are top-down)
                    float yu = 1.0f - (y + pv) / float(m_pImage->height());
                    
                    // Calculate a stratified random variation for depth-of-field
                    float lensU, lensV;
                    lensSampler.sample2D(psi, lensU, lensV);
                    
                    // Find where this pixel sample hits in the scene
                    Ray ray = m_camera.makeRay((xu - 0.5f) * aspectRatioXToY + 0.5f,
                                               yu,
                                               lensU,
                                               lensV);
                    
                    // Trace a path out, gathering estimated radiance along the path
                    pixelColor += pathTrace(ray,
                                            m_masterSet,
                                            m_lights,
                                            rng,
                                            m_lightSamplesHint,
                                            m_maxRayDepth,
                                            psi,
                                            bounceSamplers);
                }
                // Divide by the number of pixel samples (a box pixel filter, essentially)
                pixelColor /= totalPixelSamples;
                
                // Store off the computed pixel in a big buffer
                m_pImage->pixel(x, y) = pixelColor;
                
                // Reset samplers for the next pixel sample
                for (size_t i = 0; i < m_maxRayDepth; ++i)
                {
                    bounceSamplers[i]->refill();
                }
                lensSampler.refill();
                sampler.refill();
            }
        }
        
        // Deallocate bounce samplers
        for (size_t i = 0; i < m_maxRayDepth; ++i)
        {
            delete bounceSamplers[i];
        }
        delete[] bounceSamplers;
    }
    
    size_t m_xstart, m_xend, m_ystart, m_yend;
    Image *m_pImage;
    ShapeSet& m_masterSet;
    const Camera& m_camera;
    std::list<Shape*>& m_lights;
    size_t m_pixelSamplesHint, m_lightSamplesHint;
    size_t m_maxRayDepth;
};


} // namespace


namespace Rayito
{


// Construct a perspective camera, precomputing a few things to ray trace faster
PerspectiveCamera::PerspectiveCamera(float fieldOfViewInDegrees,
                                     const Point& origin,
                                     const Vector& target,
                                     const Vector& targetUpDirection,
                                     float focalDistance,
                                     float lensRadius)
    : m_origin(origin),
      m_forward((target - origin).normalized()),
      m_tanFov(std::tan(fieldOfViewInDegrees * M_PI / 180.0f)),
      m_focalDistance(focalDistance),
      m_lensRadius(lensRadius)
{
    m_right = cross(m_forward, targetUpDirection);
    m_up = cross(m_right, m_forward); // no need to normalize, it already is
}

Ray PerspectiveCamera::makeRay(float xScreen, float yScreen, float lensU, float lensV) const
{
    // Set up a camera ray given the look-at spec, FOV, and screen position to aim at
    
    // Set up ray info
    Ray ray;
    ray.m_origin = m_origin;
    ray.m_direction = m_forward +
                      m_right * ((xScreen - 0.5f) * m_tanFov) +
                      m_up * ((yScreen - 0.5f) * m_tanFov);
    ray.m_direction.normalize();
    
    if (m_lensRadius > 0)
    {
        // Modify it for DOF if necessary
        float horizontalShift = 0.0f, verticalShift = 0.0f;
        float nearDistance = 0.0f; // TODO: be able to set a custom near distance?

        // Compute shifts on the near plane due to DOF
        uniformToUniformDisk(lensU, lensV, horizontalShift, verticalShift);
        horizontalShift *= m_lensRadius;
        verticalShift *= m_lensRadius;

        // Compute local direction rays for computing focal parameters
        Vector localRayDirection((xScreen - 0.5f) * m_tanFov,
                                 (yScreen - 0.5f) * m_tanFov,
                                 1.0f);
        localRayDirection.normalize();

        // Compute primary ray focal plane intersection
        float focusT = (m_focalDistance - nearDistance) / localRayDirection.m_z;
        Point focusPoint = m_origin + ray.m_direction * focusT;

        // Modify primary ray origin
        ray.m_origin += m_right * horizontalShift + m_up * verticalShift;

        // Compute primary ray direction
        ray.m_direction = focusPoint - ray.m_origin;
        ray.m_direction.normalize();
    }
    
    return ray;
}


Color pathTrace(const Ray& ray,
                ShapeSet& scene,
                std::list<Shape*>& lights,
                Rng& rng,
                size_t lightSamplesHint,
                size_t maxRayDepth,
                size_t pixelSampleIndex,
                Sampler **bounceSamplers)
{
    // Accumulate total incoming radiance in 'result'
    Color result = Color(0.0f, 0.0f, 0.0f);
    // As we get through more and more bounces, we track how much the light is
    // diminished through each bounce
    Color throughput = Color(1.0f, 1.0f, 1.0f);
    
    // Start with the initial ray from the camera
    Ray currentRay = ray;
    
    // While we have bounces left we can still take...
    size_t numBounces = 0;
    while (numBounces < maxRayDepth)
    {
        // Trace the ray to see if we hit anything
        Intersection intersection(currentRay);
        if (!scene.intersect(intersection))
        {
            // No hit, return black (background)
            break;
        }
        
        // Add in emission when directly visible or via a specular bounce
        if (numBounces == 0)
        {
            result += throughput * intersection.m_pMaterial->emittance();
        }
        
        // Evaluate the material and intersection information at this bounce
        Point position = intersection.position();
        Vector normal = intersection.m_normal;
        Vector outgoing = -currentRay.m_direction;
        Brdf* pBrdf = NULL;
        float brdfWeight = 1.0f;
        Color matColor = intersection.m_pMaterial->evaluate(position,
                                                            normal,
                                                            outgoing,
                                                            pBrdf,
                                                            brdfWeight);
        // No BRDF?  We can't evaluate lighting, so bail.
        if (pBrdf == NULL)
        {
            return result;
        }
        
        // Evaluate direct lighting at this bounce.
        // For each light...
        for (std::list<Shape*>::iterator iter = lights.begin();
             iter != lights.end();
             ++iter)
        {
            // Set up samplers for evaluating this light
            StratifiedRandomSampler lightSampler(lightSamplesHint, lightSamplesHint, rng);
            StratifiedRandomSampler brdfSampler(lightSamplesHint, lightSamplesHint, rng);
            
            // Sample the light (with stratified random sampling to reduce noise)
            Color lightResult = Color(0.0f, 0.0f, 0.0f);
            size_t numLightSamples = lightSampler.total2DSamplesAvailable();
            for (size_t lightSampleIndex = 0; lightSampleIndex < numLightSamples; ++lightSampleIndex)
            {
                // Sample the light using MIS between the light and the BRDF.
                // This means we ask the light for a direction, and the likelihood
                // of having sampled that direction (the PDF).  Then we ask the
                // BRDF what it thinks of that direction (its PDF), and weight
                // the light sample with MIS.
                //
                // Then, we ask the BRDF for a direction, and the likelihood of
                // having sampled that direction (the PDF).  Then we ask the
                // light what it thinks of that direction (its PDF, and whether
                // that direction even runs into the light at all), and weight
                // the BRDF sample with MIS.
                //
                // By doing both samples and asking both the BRDF and light for
                // their PDF for each one, we can combine the strengths of both
                // sampling methods and get the best of both worlds.  It does
                // cost an extra shadow ray and evaluation, though, but it is
                // generally such an improvement in quality that it is very much
                // worth the overhead.
                
                // Ask the light for a random position/normal we can use for lighting
                Light *pLightShape = dynamic_cast<Light*>(*iter);
                float lsu, lsv;
                lightSampler.sample2D(lightSampleIndex, lsu, lsv);
                Point lightPoint;
                Vector lightNormal;
                float lightPdf = 0.0f;
                pLightShape->sampleSurface(position,
                                           normal,
                                           lsu, lsv, rng.nextFloat(),
                                           lightPoint,
                                           lightNormal,
                                           lightPdf);
                
                if (lightPdf > 0.0f)
                {   
                    // Ask the BRDF what it thinks of this light position (for MIS)
                    Vector lightIncoming = position - lightPoint;
                    float lightDistance = lightIncoming.normalize();
                    float brdfPdf = 0.0f;
                    float brdfResult = pBrdf->evaluateSA(lightIncoming,
                                                         outgoing,
                                                         normal,
                                                         brdfPdf);
                    if (brdfResult > 0.0f && brdfPdf > 0.0f)
                    {
                        // Fire a shadow ray to make sure we can actually see the light position
                        Ray shadowRay(position, -lightIncoming, lightDistance - kRayTMin);
                        if (!scene.doesIntersect(shadowRay))
                        {
                            // The light point is visible, so let's add that
                            // contribution (mixed by MIS)
                            float misWeightLight = powerHeuristic(1, lightPdf, 1, brdfPdf);
                            lightResult += pLightShape->emitted() *
                                           intersection.m_colorModifier * matColor *
                                           brdfResult *
                                           std::fabs(dot(-lightIncoming, normal)) *
                                           misWeightLight / (lightPdf * brdfWeight);
                        }
                    }
                }
                
                // Ask the BRDF for a sample direction
                float bsu, bsv;
                brdfSampler.sample2D(lightSampleIndex, bsu, bsv);
                Vector brdfIncoming;
                float brdfPdf = 0.0f;
                float brdfResult = pBrdf->sampleSA(brdfIncoming,
                                                   outgoing,
                                                   normal,
                                                   bsu,
                                                   bsv,
                                                   brdfPdf);
                if (brdfPdf > 0.0f && brdfResult > 0.0f)
                {
                    Intersection shadowIntersection(Ray(position, -brdfIncoming));
                    bool intersected = scene.intersect(shadowIntersection);
                    if (intersected && shadowIntersection.m_pShape == pLightShape)
                    {
                        // Ask the light what it thinks of this direction (for MIS)
                        lightPdf = pLightShape->intersectPdf(shadowIntersection);
                        if (lightPdf > 0.0f)
                        {
                            // BRDF chose the light, so let's add that
                            // contribution (mixed by MIS)
                            float misWeightBrdf = powerHeuristic(1, brdfPdf, 1, lightPdf);
                            lightResult += pLightShape->emitted() * 
                                           intersection.m_colorModifier * matColor * brdfResult *
                                           std::fabs(dot(-brdfIncoming, normal)) * misWeightBrdf /
                                           (brdfPdf * brdfWeight);
                        }
                    }
                }
            }
            // Average light samples
            if (numLightSamples)
            {
                lightResult /= numLightSamples;
            }
            
            // Add direct lighting at this bounce (modified by how much the
            // previous bounces have dimmed it)
            result += throughput * lightResult;
        }
        
        // Sample the BRDF to find the direction the next leg of the path goes in
        float brdfSampleU, brdfSampleV;
        bounceSamplers[numBounces]->sample2D(pixelSampleIndex, brdfSampleU, brdfSampleV);
        Vector incoming;
        float incomingBrdfPdf = 0.0f;
        float incomingBrdfResult = pBrdf->sampleSA(incoming,
                                                   outgoing,
                                                   normal,
                                                   brdfSampleU,
                                                   brdfSampleV,
                                                   incomingBrdfPdf);

        if (incomingBrdfPdf > 0.0f)
        {
            currentRay.m_origin = position;
            currentRay.m_direction = -incoming;
            currentRay.m_tMax = kRayTMax;
            // Reduce lighting effect for the next bounce based on this bounce's BRDF
            throughput *= intersection.m_colorModifier * matColor * incomingBrdfResult *
                          (std::fabs(dot(-incoming, normal)) /
                          (incomingBrdfPdf * brdfWeight));
        }
        else
        {
            break; // BRDF is zero, stop bouncing
        }
        
        numBounces++;
    }
    
    // This represents an estimate of the total light coming in along the path
    return result;
}


Image* raytrace(ShapeSet& scene,
                const Camera& cam,
                size_t width,
                size_t height,
                size_t pixelSamplesHint,
                size_t lightSamplesHint,
                size_t maxRayDepth)
{
    // Get light list from the scene
    std::list<Shape*> lights;
    scene.findLights(lights);
    
    scene.prepare();
    
    // Set up the output image
    Image *pImage = new Image(width, height);
    
    // Set up render threads; we make as much as 16 chunks of the image that
    // can render in parallel.
    const size_t kChunkDim = 4;
    
    // Chunk size is the number of pixels per image chunk (we have to take care
    // to deal with tiny images)
    size_t xChunkSize = width >= kChunkDim ? width / kChunkDim : 1;
    size_t yChunkSize = height >= kChunkDim ? height / kChunkDim : 1;
    // Chunks are the number of chunks in each dimension we can chop the image
    // into (again, taking care to deal with tiny images, and also images that
    // don't divide clealy into 4 chunks)
    size_t xChunks = width > kChunkDim ? width / xChunkSize : 1;
    size_t yChunks = height > kChunkDim ? height / yChunkSize : 1;
    if (xChunks * xChunkSize < width) xChunks++;
    if (yChunks * yChunkSize < height) yChunks++;
    
    // Set up render threads
    size_t numRenderThreads = xChunks * yChunks;
    RenderThread **renderThreads = new RenderThread*[numRenderThreads];
    
    // Launch render threads
    for (size_t yc = 0; yc < yChunks; ++yc)
    {
        // Get the row start/end (making sure the last chunk doesn't go off the end)
        size_t yStart = yc * yChunkSize;
        size_t yEnd = std::min((yc + 1) * yChunkSize, height);
        for (size_t xc = 0; xc < xChunks; ++xc)
        {
            // Get the column start/end (making sure the last chunk doesn't go off the end)
            size_t xStart = xc * xChunkSize;
            size_t xEnd = std::min((xc + 1) * xChunkSize, width);
            // Render the chunk!
            renderThreads[yc * xChunks + xc] = new RenderThread(xStart,
                                                                xEnd,
                                                                yStart,
                                                                yEnd,
                                                                pImage,
                                                                scene,
                                                                cam,
                                                                lights,
                                                                pixelSamplesHint,
                                                                lightSamplesHint,
                                                                maxRayDepth);
            renderThreads[yc * xChunks + xc]->start();
        }
    }
    
    // Wait until the render finishes
    bool stillRunning;
    do
    {
        // See if any render thread is still going...
        stillRunning = false;
        for (size_t i = 0; i < numRenderThreads; ++i)
        {
            if (renderThreads[i]->isRunning())
            {
                stillRunning = true;
                break;
            }
        }
        if (stillRunning)
        {
            // Give up the CPU so the render threads can do their thing
            QThread::yieldCurrentThread();
        }
    } while (stillRunning);
    
    // Clean up render thread objects
    for (size_t i = 0; i < numRenderThreads; ++i)
    {
        delete renderThreads[i];
    }
    delete[] renderThreads;
    
    // We made a picture!
    return pImage;
}


} // namespace Rayito
