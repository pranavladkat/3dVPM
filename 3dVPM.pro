TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \

include(deployment.pri)
qtcAddDeployment()

HEADERS += \


QMAKE_CXXFLAGS += -std=c++11
