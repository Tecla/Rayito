Stage 4:

Everything in stage 3, plus the following changes...

* Qt app wrapper so you can easily visualize the results, change samples/pixel, samples/light, etc
* QThread implementation of multithreaded rendering (so we use more than one CPU core now)

MAKE YOUR OWN CHANGES!
...


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

