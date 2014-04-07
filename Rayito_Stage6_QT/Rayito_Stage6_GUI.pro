#-------------------------------------------------
#
# Project created by QtCreator 2012-08-23T16:57:47
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Rayito_Stage6_GUI
TEMPLATE = app


SOURCES += main.cpp\
        MainWindow.cpp \
    RaytraceMain.cpp \
    OBJMesh.cpp

HEADERS  += MainWindow.h \
    rayito.h \
    RMath.h \
    RRay.h \
    RMaterial.h \
    RLight.h \
    RScene.h \
    RSampling.h \
    RAccel.h \
    RMesh.h

FORMS    += MainWindow.ui

QMAKE_CXXFLAGS += -Wno-unused-parameter
