Stage 7 (Qt):

Everything in stage 6, plus the following changes...

* Transformations on shapes (featuring quaternion rotations)
* Transform motion blur, camera shutter settings
* Improved samplers (correlated multi-jitter), fixed number of samples for all lights
* Perfect specular reflection BRDF
* A second scene to render showing off more motion blur

Please see the code comments, they offer explanations of each feature.


Stage 6 (Qt):

Everything in stage 5, plus the following changes...

* BVH acceleration structure for both the scene and for polygon meshes
* Polygon mesh object
* OBJ reader

Please see the code comments, they offer explanation of BVHs and OBJ files
and ray-mesh intersections.


Stage 5 (Qt):

Everything in stage 4, plus the following changes...

* Path tracing!  Rayito stage 5 can do a full GI solution.
* Sampler classes (only stratified random for now)
* Camera classes (only perspective for now, but with depth-of-field!)
* Broke up rayito.h into many other headers for better organization
* Lots of little additions and organization fixes, more comments

Please see the code comments, they are meant to help the reader understand how
to do path tracing for real.  It is a bit simplified, but it does it right.


Stage 4 (Qt):

Everything in stage 3, plus the following changes...

* Qt app wrapper so you can easily visualize the results, change samples/pixel, samples/light, etc
* QThread implementation of multithreaded rendering (so we use more than one CPU core now)


Stage 3:

Everything in stage 2, plus the following changes...

* Shape light, which supports using any finite shape as a light (a sphere is demoed)
* Material class, with Lambert, Phong, and Emitter materials available
* Sphere shape class added
* Stratified random sampling of pixel samples
* Stratified random sampling of light samples


Stage 2:

Everything in stage 1, plus the following changes...

* Light shape subclass, with (double-sided) rectangle area light subclass
* Method to find the lights in a scene
* Shape sampling method to find a random point on a shape (used for lights currently)
* Emission added to intersection results
* Infinite plane has a procedural bullseye texture option
* Multiple samples per pixel supported (antialiasing)
* Fast random number stream added
* Shadow rays cast
* Simple Lambertian diffuse shading


Stage 1:

* Math support: RGB colors, 3D vectors/points, ray definition
* Scene graph: abstract shape base class with ray intersectio method, infinite plane subclass
* Camera ray construction
* Single-sample-per-pixel ray-traced image

