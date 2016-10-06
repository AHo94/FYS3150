TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp\
     vec3.cpp \
    celestials.cpp \
    solarsystem.cpp \
    odesolvers.cpp

HEADERS += \
    vec3.h \
    celestials.h \
    solarsystem.h \
    odesolvers.h
