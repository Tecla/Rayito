#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


namespace Ui {
class MainWindow;
}


class MainWindow : public QMainWindow
{
   Q_OBJECT
   
public:
   explicit MainWindow(QWidget *pParent = NULL);
   ~MainWindow();
   
private slots:
   void on_actionQuit_triggered();
   
   void on_renderButton_clicked();
   
private:
   Ui::MainWindow *ui;
};


#endif // MAINWINDOW_H
