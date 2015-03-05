TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += -I /home/pranav/Libraries/petsc/include -I /home/pranav/Libraries/petsc/ubuntu/include -I /home/pranav/Libraries/openmpi-1.8.4/openmpi-1.8.4/include \
               /usr/local/MATLAB/R2013a/extern/include \

LIBS += -llapack -lblas -Wl,-rpath,/home/pranav/Libraries/petsc/ubuntu/lib -Wl,-rpath,/home/pranav/Libraries/petsc/ubuntu/lib -L/home/pranav/Libraries/petsc/ubuntu/lib -lpetsc -Wl,-rpath,/home/pranav/Libraries/petsc/ubuntu/lib -lf2clapack -lf2cblas -lm -lpthread -lssl -lcrypto -lm -Wl,-rpath,/home/pranav/Libraries/openmpi-1.8.4/openmpi-1.8.4/lib -L/home/pranav/Libraries/openmpi-1.8.4/openmpi-1.8.4/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -Wl,-rpath,/home/pranav/Libraries/openmpi-1.8.4/openmpi-1.8.4/lib -L/home/pranav/Libraries/openmpi-1.8.4/openmpi-1.8.4/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -Wl,-rpath,/home/pranav/Libraries/openmpi-1.8.4/openmpi-1.8.4/lib -lmpi -lgcc_s -lpthread -ldl \
        -L /usr/local/MATLAB/R2013a/bin/glnxa64 -leng -lmx

unix:QMAKE_RPATHDIR += /usr/local/MATLAB/R2013a/bin/glnxa64

QMAKE_CXX = g++
QMAKE_CXXFLAGS += -std=c++11

SOURCES += \
    src/main.cpp \
    src/surface.cpp \
    src/plot3d.cpp \
    src/vtk_writer.cpp \
    src/solver.cpp \
    src/wake.cpp \
    src/parameters.cpp \
    src/domain.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    src/surface.hpp \
    src/plot3d.hpp \
    src/vector3d.h \
    src/vtk_writer.hpp \
    src/solver.hpp \
    src/wake.hpp \
    src/parameters.hpp \
    src/domain.hpp
