TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp
LIBS += -LC:/Armadillo/newblas -llibblas
LIBS += -LC:/Armadillo/newblas -lliblapack

INCLUDEPATH += C:/Armadillo/include
DEPENDPATH += C:/Armadillo/include
