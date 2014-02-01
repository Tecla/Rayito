#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "rayito.h"

#include <QGraphicsScene>


// Declaration for the function that does the actual work
Rayito::Image* raytrace(size_t width, size_t height, size_t pixelSamplesHint, size_t lightSamplesHint);


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
    Rayito::Image *pImage = raytrace(ui->widthSpinBox->value(),
                                     ui->heightSpinBox->value(),
                                     ui->pixelSamplesSpinBox->value(),
                                     ui->lightSamplesSpinBox->value());
    
    // Convert from floating-point RGB to 32-bit ARGB format,
    // applying exposure and gamma along the way
    
    // A gamma'd value is value^(1/gamma)
    float gammaExponent = 1.0f / (float)ui->gammaSpinBox->value();
    // An exposure'd value is value*2^exposure
    float exposure = std::pow(2.0f, (float)ui->exposureSpinBox->value());
    uchar *argbPixels = new uchar[pImage->width() * pImage->height() * 4];
    for (size_t y = 0; y < pImage->height(); ++y)
    {
        for (size_t x = 0; x < pImage->width(); ++x)
        {
            size_t pixelOffset = (y * pImage->width() + x) * 4;
            Rayito::Color color = pImage->pixel(x, y);
            // Combined gamma+exposure: result = (value*(2^exposure))^(1/gamma)
            color.m_r = std::pow(color.m_r * exposure, gammaExponent);
            color.m_g = std::pow(color.m_g * exposure, gammaExponent);
            color.m_b = std::pow(color.m_b * exposure, gammaExponent);
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
