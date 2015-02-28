TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \
    src/surface.cpp \
    src/plot3d.cpp \
    src/vtk_writer.cpp \
    src/solver.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    src/surface.hpp \
    src/plot3d.hpp \
    src/vector3d.h \
    src/parameters.h \
    src/vtk_writer.hpp \
    src/solver.hpp


QMAKE_CXXFLAGS += -std=c++11
