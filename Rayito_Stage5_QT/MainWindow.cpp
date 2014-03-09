#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "rayito.h"

#include <QGraphicsScene>


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


void MainWindow::on_renderButton_clicked()
{
    // Make a picture...
    
    // Available materials
    Rayito::DiffuseMaterial blueishLambert(Rayito::Color(0.7f, 0.7f, 0.9f));
    Rayito::DiffuseMaterial purplishLambert(Rayito::Color(0.8f, 0.3f, 0.7f));
    Rayito::DiffuseMaterial yellowishLambert(Rayito::Color(0.7f, 0.7f, 0.2f));
    Rayito::GlossyMaterial bluishGlossy(Rayito::Color(0.5f, 0.3f, 0.8f), 0.3);
    Rayito::GlossyMaterial greenishGlossy(Rayito::Color(0.3f, 0.9f, 0.3f), 0.1f);
    
    // The 'scene'
    Rayito::ShapeSet masterSet;
    
    // Put a ground plane in (with bullseye texture!)
    Rayito::Plane plane(Rayito::Point(0.0f, -2.0f, 0.0f),
                        Rayito::Vector(0.0f, 1.0f, 0.0f),
                        &blueishLambert,
                        true); // Last parameter is whether to do the bullseye texture or not
    masterSet.addShape(&plane);
    
    // Add a pile-o-spheres with a few interesting materials
    
    Rayito::Sphere sphere1(Rayito::Point(3.0f, -1.0f, 0.0f),
                           1.0f,
                           &purplishLambert);
    masterSet.addShape(&sphere1);
    
    Rayito::Sphere sphere2(Rayito::Point(-3.0f, 0.0f, -2.0f),
                           2.0f,
                           &greenishGlossy);
    masterSet.addShape(&sphere2);
    
    Rayito::Sphere sphere3(Rayito::Point(1.5f, -1.5f, 2.5f),
                           0.5f,
                           &bluishGlossy);
    masterSet.addShape(&sphere3);
    
    Rayito::Sphere sphere4(Rayito::Point(-2.0f, -1.5f, 1.0f),
                           0.5f,
                           &yellowishLambert);
    masterSet.addShape(&sphere4);
    
    // Add an area light
    Rayito::RectangleLight areaLight(Rayito::Point(-1.5f, 4.0f, -1.5f),
                                     Rayito::Vector(3.0f, 0.0f, 0.0f),
                                     Rayito::Vector(0.0f, 0.0f, 3.0f),
                                     Rayito::Color(1.0f, 1.0f, 1.0f),
                                     5.0f);
    masterSet.addShape(&areaLight);

    // Add an area light based on a shape (a sphere)
    Rayito::Sphere sphereForLight(Rayito::Point(0.0f, 0.5f, 2.0f),
                                  0.5f,
                                  &blueishLambert);
    Rayito::ShapeLight sphereLight(&sphereForLight, Rayito::Color(1.0f, 1.0f, 0.3f), 10.0f);
    masterSet.addShape(&sphereLight);
    
    // Create the camera based on user settings
    Rayito::PerspectiveCamera cam((float)ui->camFovSpinBox->value(),
                                  Rayito::Point(0.0f, 5.0f, 15.0f),
                                  Rayito::Point(0.0f, 0.0f, 0.0f),
                                  Rayito::Point(0.0f, 1.0f, 0.0f),
                                  (float)ui->focalDistanceSpinBox->value(),
                                  (float)ui->lensRadiusSpinBox->value());
    
    // Ray trace!
    Rayito::Image *pImage = raytrace(masterSet,
                                     cam,
                                     ui->widthSpinBox->value(),
                                     ui->heightSpinBox->value(),
                                     ui->pixelSamplesSpinBox->value(),
                                     ui->lightSamplesSpinBox->value(),
                                     ui->rayDepthSpinBox->value());
    
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
            Rayito::Color color = pImage->pixel(x, y);
            // Check for negative values (we don't like those).  Make them green.
            if (color.m_r < 0.0f || color.m_g < 0.0f || color.m_b < 0.0f)
            {
                color = Rayito::Color(0.0f, 1.0f, 0.0f);
            }
            else
            {
                // Combined gamma+exposure: result = (value*(2^exposure))^(1/gamma)
                color.m_r = std::pow(color.m_r * exposure, gammaExponent);
                color.m_g = std::pow(color.m_g * exposure, gammaExponent);
                color.m_b = std::pow(color.m_b * exposure, gammaExponent);
                // Check for NaNs (not-a-number), we HATE those.  Make them blue.
                if (color.m_r != color.m_r || color.m_g != color.m_g || color.m_b != color.m_b)
                {
                    color = Rayito::Color(0.0f, 0.0f, 1.0f);
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
    delete pImage;
}

void MainWindow::on_actionRender_Scene_triggered()
{
   
}
