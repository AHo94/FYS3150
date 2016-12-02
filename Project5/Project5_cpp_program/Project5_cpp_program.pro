TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    vec3.cpp \
    metropolis_quantum.cpp \
    wavefunctions.cpp

HEADERS += \
    vec3.h \
    metropolis_quantum.h \
    wavefunctions.h
