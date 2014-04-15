#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "rayito.h"
#include "RMesh.h"

#include <QGraphicsScene>


using namespace Rayito;


MainWindow::MainWindow(QWidget *pParent)
    : QMainWindow(pParent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    
    // We have to allocate this, because Qt doesn't do it for us by default.
    // Oftentimes you would create your own QGraphicsScene subclass, but we don't.
    QGraphicsScene *pScene = new QGraphicsScene(this);
    ui->renderGraphicsView->setScene(pScene);
}


MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionQuit_triggered()
{
    qApp->quit();
}


void MainWindow::displayImage(Image *pImage)
{
    // Convert from floating-point RGB to 32-bit ARGB format,
    // applying exposure and gamma along the way
    
    // A gamma'd value is value^(1/gamma)
    float gammaExponent = 1.0f / (float)ui->gammaSpinBox->value();
    // An exposure'd value is value*2^exposure (applied before gamma)
    float exposure = std::pow(2.0f, (float)ui->exposureSpinBox->value());
    uchar *argbPixels = new uchar[pImage->width() * pImage->height() * 4];
    for (size_t y = 0; y < pImage->height(); ++y)
    {
        for (size_t x = 0; x < pImage->width(); ++x)
        {
            size_t pixelOffset = (y * pImage->width() + x) * 4;
            Color color = pImage->pixel(x, y);
            // Check for negative values (we don't like those).  Make them green.
            if (color.m_r < 0.0f || color.m_g < 0.0f || color.m_b < 0.0f)
            {
                color = Color(0.0f, 1.0f, 0.0f);
            }
            else
            {
                // Combined gamma and exposure: result = (value*(2^exposure))^(1/gamma)
                color.m_r = std::pow(color.m_r * exposure, gammaExponent);
                color.m_g = std::pow(color.m_g * exposure, gammaExponent);
                color.m_b = std::pow(color.m_b * exposure, gammaExponent);
                // Check for NaNs (not-a-number), we HATE those.  Make them blue.
                if (color.m_r != color.m_r || color.m_g != color.m_g || color.m_b != color.m_b)
                {
                    color = Color(0.0f, 0.0f, 1.0f);
                }
            }
            // We're displaying in LDR in the end (you've got an exposure control, after all)
            color.clamp();
            // Stuff actual 32-bit ARGB pixel data in
            argbPixels[pixelOffset + 3] = 0xFF;
            argbPixels[pixelOffset + 2] = static_cast<uchar>(color.m_r * 255.0f);
            argbPixels[pixelOffset + 1] = static_cast<uchar>(color.m_g * 255.0f);
            argbPixels[pixelOffset + 0] = static_cast<uchar>(color.m_b * 255.0f);
        }
    }
    
    // Make an image, then make a pixmap for the graphics scene
    QImage image(argbPixels,
                 static_cast<int>(pImage->width()),
                 static_cast<int>(pImage->height()),
                 QImage::Format_ARGB32_Premultiplied);
    
    ui->renderGraphicsView->scene()->clear();
    ui->renderGraphicsView->scene()->addPixmap(QPixmap::fromImage(image));
    
    // Clean up the image and pixel conversion buffers
    delete[] argbPixels;
}

inline Mesh* makeCube()
{
    std::vector<Face> faces;
    std::vector<Vector> normals;
    std::vector<Point> vertices;
    vertices.push_back(Point(0.0f, 0.0f, 0.0f));
    vertices.push_back(Point(1.0f, 0.0f, 0.0f));
    vertices.push_back(Point(1.0f, 1.0f, 0.0f));
    vertices.push_back(Point(0.0f, 1.0f, 0.0f));
    vertices.push_back(Point(0.0f, 0.0f, 1.0f));
    vertices.push_back(Point(1.0f, 0.0f, 1.0f));
    vertices.push_back(Point(1.0f, 1.0f, 1.0f));
    vertices.push_back(Point(0.0f, 1.0f, 1.0f));
    faces.push_back(Face());
    faces.back().m_vertexIndices.push_back(0);
    faces.back().m_vertexIndices.push_back(1);
    faces.back().m_vertexIndices.push_back(2);
    faces.back().m_vertexIndices.push_back(3);
    faces.push_back(Face());
    faces.back().m_vertexIndices.push_back(1);
    faces.back().m_vertexIndices.push_back(5);
    faces.back().m_vertexIndices.push_back(6);
    faces.back().m_vertexIndices.push_back(2);
    faces.push_back(Face());
    faces.back().m_vertexIndices.push_back(5);
    faces.back().m_vertexIndices.push_back(4);
    faces.back().m_vertexIndices.push_back(7);
    faces.back().m_vertexIndices.push_back(6);
    faces.push_back(Face());
    faces.back().m_vertexIndices.push_back(4);
    faces.back().m_vertexIndices.push_back(0);
    faces.back().m_vertexIndices.push_back(3);
    faces.back().m_vertexIndices.push_back(7);
    faces.push_back(Face());
    faces.back().m_vertexIndices.push_back(3);
    faces.back().m_vertexIndices.push_back(2);
    faces.back().m_vertexIndices.push_back(6);
    faces.back().m_vertexIndices.push_back(7);
    faces.push_back(Face());
    faces.back().m_vertexIndices.push_back(3);
    faces.back().m_vertexIndices.push_back(2);
    faces.back().m_vertexIndices.push_back(6);
    faces.back().m_vertexIndices.push_back(7);
    return new Mesh(vertices, normals, faces, NULL);
}

void MainWindow::on_renderButton_clicked()
{
    // Make a picture...
    
    // Available materials
    DiffuseMaterial blueishLambert(Color(0.6f, 0.6f, 0.9f));
    DiffuseMaterial purplishLambert(Color(0.8f, 0.3f, 0.7f));
    DiffuseMaterial reddishLambert(Color(0.8f, 0.3f, 0.1f));
    GlossyMaterial bluishGlossy(Color(0.5f, 0.3f, 0.8f), 0.3);
    GlossyMaterial greenishGlossy(Color(0.3f, 0.9f, 0.3f), 0.1f);
    GlossyMaterial reddishGlossy(Color(0.8f, 0.1f, 0.1f), 0.3f);
    ReflectionMaterial reflective(Color(0.7f, 0.7f, 0.2f));
    
    // The 'scene'
    ShapeSet masterSet;
    
    // Put a ground plane in (with bullseye texture!)
    // Last parameter is whether to do the bullseye texture or not
    Plane plane(Point(), Vector(0.0f, 1.0f, 0.0f), &blueishLambert, true);
    plane.transform().translate(0.0f, Vector(0.0f, -2.0f, 0.0f));
    masterSet.addShape(&plane);
    
    // Add a pile-o-spheres with a few interesting materials
    
    Sphere sphere1(Point(), 1.0f, &purplishLambert);
    sphere1.transform().setTranslation(0.0f, Vector(2.0f, -1.0f, 0.0f));
    sphere1.transform().setTranslation(1.0f, Vector(3.0f, -1.0f, 0.0f));
    masterSet.addShape(&sphere1);
    
    Sphere sphere2(Point(), 2.0f, &greenishGlossy);
    sphere2.transform().translate(0.0f, Vector(-3.0f, 0.0f, -2.0f));
    masterSet.addShape(&sphere2);
    
    Sphere sphere3(Point(), 0.5f, &bluishGlossy);
    sphere3.transform().translate(0.0f, Vector(1.5f, -1.5f, 2.5f));
    masterSet.addShape(&sphere3);
    
    Sphere sphere4(Point(), 0.5f, &reflective);
    sphere4.transform().translate(0.0f, Vector(-2.0, -1.5f, 1.0f));
    masterSet.addShape(&sphere4);
    
    // Add a manually created mesh (a box), and read an OBJ file into a mesh
    
    Mesh *cubeMesh = makeCube();
    cubeMesh->setMaterial(&reddishLambert);
    cubeMesh->transform().translate(0.0f, Vector(0.0f, -2.0f, -2.0f));
    cubeMesh->transform().rotate(1.0f, Quaternion(Vector(0.0f, 1.0f, 0.0f), M_PI / 4.0f));
    masterSet.addShape(cubeMesh);

    Mesh* pOBJMesh = createFromOBJFile("../models/bumpy.obj");
    pOBJMesh->setMaterial(&reddishGlossy);
    pOBJMesh->transform().setTranslation(0.0f, Vector(0.2f, 0.0f, 0.0f));
    pOBJMesh->transform().rotate(0.5f, Quaternion(Vector(0.0f, 1.0f, 0.0f), M_PI / 4.0f));
    pOBJMesh->transform().rotate(1.0f, Quaternion(Vector(0.0f, 1.0f, 0.0f), M_PI / 2.0f));
#if MAKE_OBJ_A_MESH_LIGHT
    // For some fun, you can turn the OBJ mesh into a light (it's a bit noisy, though)
    ShapeLight meshLight(pOBJMesh, Color(1.0f, 1.0f, 1.0f), 10.0f);
    masterSet.addShape(&meshLight);
#else
    masterSet.addShape(pOBJMesh);
#endif

    // Add an area light
    RectangleLight areaLight(Point(),
                             Vector(3.0f, 0.0f, 0.0f),
                             Vector(0.0f, 0.0f, 3.0f),
                             Color(1.0f, 1.0f, 1.0f),
                             5.0f);
    areaLight.transform().setTranslation(0.0f, Vector(-1.5f, 4.0f, -1.5f));
    // Uncomment this to have the rect light hinge-swing downward
//    areaLight.transform().setRotation(1.0f, Quaternion(Vector(0.0f, 0.0f, 1.0f), -M_PI / 4.0f));
    masterSet.addShape(&areaLight);

    // Add an area light based on a shape (a sphere)
    Sphere sphereForLight(Point(), 0.1f, &blueishLambert);
    sphereForLight.transform().setTranslation(0.0f, Vector(0.0f, 0.5f, 4.0f));
    sphereForLight.transform().setTranslation(0.33f, Vector(0.0f, 1.5f, 4.0f));
    sphereForLight.transform().setTranslation(0.67f, Vector(1.0f, 1.5f, 4.0f));
    sphereForLight.transform().setTranslation(1.0f, Vector(1.0f, 0.5f, 4.0f));
    ShapeLight sphereLight(&sphereForLight, Color(1.0f, 1.0f, 0.3f), 100.0f);
    masterSet.addShape(&sphereLight);
    
    // Create the camera based on user settings
    PerspectiveCamera cam((float)ui->camFovSpinBox->value(),
                          Point(-4.0f, 5.0f, 15.0f),
                          Point(0.0f, 0.0f, 0.0f),
                          Point(0.0f, 1.0f, 0.0f),
                          (float)ui->focalDistanceSpinBox->value(),
                          (float)ui->lensRadiusSpinBox->value(),
                          (float)ui->shutterOpenSpinBox->value(),
                          (float)ui->shutterCloseSpinBox->value());
    
    // Ray trace!
    Image *pImage = raytrace(masterSet,
                             cam,
                             (size_t)ui->widthSpinBox->value(),
                             (size_t)ui->heightSpinBox->value(),
                             (unsigned int)ui->pixelSamplesSpinBox->value(),
                             (unsigned int)ui->lightSamplesSpinBox->value(),
                             (unsigned int)ui->rayDepthSpinBox->value());
    
    displayImage(pImage);
    
    // Clean up the scene and render
    delete pImage;
    delete pOBJMesh;
}

// Let's have some fun, and make a scene where objects fall and bounce off the
// ground using real physics kinematics.
inline Point kinematicPosition(const Point& start,
                               const Vector& velocity,
                               float time,
                               const Vector& gravity = Vector(0.0f, -9.8f, 0.0f),
                               float groundHeight = 0.0f)
{
    // Kinematic equations:
    //     p = p0 + v*t
    //     p = t*(vi + vf)/2
    //     vf = v0 + a*t
    //     p = p0 + (v0 + v0 + a*t)*t/2 = p0 + v0*t + a*t*t/2
    // Solving for time to find when we smack the plane:
    //     t = (-v0 + sqrt(v0*v0 - 4*a*p0/2))/(2*a/2) = (-v0 + sqrt(v0*v0 - 2*a*p0))/a
    
    // Gravity direction
    Vector up = -gravity.normalized();
    float vUp = dot(velocity, up);
    float pUp = dot(start, up);
    float aUp = -gravity.length();
    
    float discriminant = vUp * vUp - 2.0f * aUp * pUp;
    if (discriminant > 0.0f)
    {
        float intersectionTime = (-vUp - std::sqrt(discriminant)) / aUp;
        if (intersectionTime < time)
        {
            // Find intersection position
            Point isect = start + velocity * intersectionTime +
                          gravity * intersectionTime * intersectionTime * 0.5f;
            Vector isectVelocity = (velocity + gravity * intersectionTime);
            Vector reboundVelocity = isectVelocity - 2.0f * up * dot(isectVelocity, up);
            float reboundTime = time - intersectionTime;
            return isect + reboundVelocity * reboundTime +
                   gravity * reboundTime * reboundTime * 0.5f;
        }
    }
    // Do the usual kinematic result
    return start + velocity * time + gravity * time * time * 0.5f;
}

void MainWindow::on_renderButton2_clicked()
{
    // Make a picture...
    
    // Available materials
    DiffuseMaterial blueishLambert(Color(0.6f, 0.6f, 0.9f));
    GlossyMaterial yellowishGlossy(Color(0.9f, 0.9f, 0.3f), 0.3f);
    DiffuseMaterial redLambert(Color(1.0f, 0.2f, 0.2f));
    
    // The 'scene'
    ShapeSet masterSet;
    
    // Put a ground plane in (with bullseye texture!)
    // Last parameter is whether to do the bullseye texture or not
    Plane plane(Point(), Vector(0.0f, 1.0f, 0.0f), &redLambert, true);
    masterSet.addShape(&plane);
    
    // Add a pile-o-spheres at various stages of falling/bouncing
    
    Sphere spheres[10];
    Point start(-10.0f, 10.0f, 0.0f);
    Vector velocity(4.5f, 0.0f, 0.0f);
    float timeOffset = 0.0f;
    const float timeDelta = 0.2f;
    for (unsigned int i = 0; i < 10; ++i)
    {
        Point position0 = kinematicPosition(start, velocity, timeOffset);
        Point position1 = kinematicPosition(start, velocity, timeOffset + timeDelta);
        
        spheres[i].transform().setTranslation(0.0f, position0);
        spheres[i].transform().setTranslation(1.0f, position1);
        spheres[i].setMaterial(&blueishLambert);
        
        masterSet.addShape(&spheres[i]);
        
        timeOffset += timeDelta * 2.0f;
    }
    
    // Add a pile-o-cubes at various stages of falling/bouncing, and rotating
    
    Mesh* cubes[10];
    start = Point(10.0f, 10.0f, 2.0f);
    velocity = Vector(-4.5f, 0.0f, 0.0f);
    timeOffset = 0.0f;
    for (unsigned int i = 0; i < 10; ++i)
    {
        Point position0 = kinematicPosition(start, velocity, timeOffset);
        Point position1 = kinematicPosition(start, velocity, timeOffset + timeDelta);
        float rotation0 = timeOffset * M_PI * 0.5;
        if (rotation0 > M_PI * 2.0f)
            rotation0 -= M_PI * 2.0f;
        float rotation1 = rotation0 + timeDelta * M_PI * 0.5;
        
        cubes[i] = makeCube();
        cubes[i]->transform().setTranslation(0.0f, position0);
        cubes[i]->transform().setRotation(0.0f, Quaternion(Vector(1.0f, 0.0f, 1.0f).normalized(), rotation0));
        cubes[i]->transform().setTranslation(1.0f, position1);
        cubes[i]->transform().setRotation(1.0f, Quaternion(Vector(1.0f, 0.0f, 1.0f).normalized(), rotation1));
        cubes[i]->setMaterial(&yellowishGlossy);
        
        masterSet.addShape(cubes[i]);
        
        timeOffset += timeDelta * 2.0f;
    }
    
    // Add an area light
    RectangleLight areaLight(Point(),
                             Vector(2.0f, 0.0f, 0.0f),
                             Vector(0.0f, 0.0f, 2.0f),
                             Color(1.0f, 1.0f, 1.0f),
                             50.0f);
    areaLight.transform().setTranslation(0.0f, Vector(-1.0f, 15.0f, 1.0f));
    masterSet.addShape(&areaLight);

    // Create the camera based on user settings
    PerspectiveCamera cam((float)ui->camFovSpinBox->value(),
                          Point(-4.0f, 10.0f, 30.0f),
                          Point(0.0f, 5.0f, 0.0f),
                          Point(0.0f, 1.0f, 0.0f),
                          (float)ui->focalDistanceSpinBox->value(),
                          (float)ui->lensRadiusSpinBox->value(),
                          (float)ui->shutterOpenSpinBox->value(),
                          (float)ui->shutterCloseSpinBox->value());
    
    // Ray trace!
    Image *pImage = raytrace(masterSet,
                             cam,
                             (size_t)ui->widthSpinBox->value(),
                             (size_t)ui->heightSpinBox->value(),
                             (unsigned int)ui->pixelSamplesSpinBox->value(),
                             (unsigned int)ui->lightSamplesSpinBox->value(),
                             (unsigned int)ui->rayDepthSpinBox->value());
    
    displayImage(pImage);
    
    // Clean up the scene and render
    delete pImage;
    for (unsigned int i = 0; i < 10; ++i)
    {
        delete cubes[i];
    }
}

void MainWindow::on_actionRender_Scene_triggered()
{
   
}
