#-------------------------------------------------
#
# Project created by QtCreator 2012-08-23T16:57:47
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Rayito_stage4_GUI
TEMPLATE = app


SOURCES += main.cpp\
        MainWindow.cpp \
    RaytraceMain.cpp

HEADERS  += MainWindow.h \
    rayito.h

FORMS    += MainWindow.ui

QMAKE_CXXFLAGS += -Wno-unused-parameter

