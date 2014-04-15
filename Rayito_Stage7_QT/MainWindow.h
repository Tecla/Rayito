#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


namespace Ui {
class MainWindow;
}


namespace Rayito
{
    class Image;
}


class MainWindow : public QMainWindow
{
   Q_OBJECT
   
public:
   explicit MainWindow(QWidget *pParent = NULL);
   ~MainWindow();
   
    void displayImage(Rayito::Image *pImage);
    
private slots:
   void on_actionQuit_triggered();
   
   void on_renderButton_clicked();
   void on_renderButton2_clicked();
   
   void on_actionRender_Scene_triggered();
   
private:
   Ui::MainWindow *ui;
};


#endif // MAINWINDOW_H
