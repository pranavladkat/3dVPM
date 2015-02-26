TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \
    src/surface.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    src/surface.hpp


QMAKE_CXXFLAGS += -std=c++11
