TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp
#INCLUDEPATH += C:\Users\Alex\Documents\Armadillo\armadillo-7.400.2\armadillo-7.400.2\include
#INCLUDEPATH += C:\Users\Alex\Documents\Armadillo\armadillo-7.400.2\armadillo-7.400.2\examples\lib_win64
#INCLUDEPATH += -L..\examples\lib_win64
#INCLUDEPATH += -L..\include

#INCLUDEPATH += C:\Armadillo\include
#INCLUDEPATH += C:\Armadillo\examples\lib_win64
#LIBS += -lblas_win64_MT -llapack_win64_MT -larmadillo
INCLUDEPATH += C:\Armadillo\include
#INCLUDEPATH += C:\Armadillo\examples\lib_win64
LIBS += \
    -LC:\Armadillo\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT

#LIBS += -LC:\Armadillo\include\
#    -larmadillo
#LIBS += -LC:\Armadillo\examples\lib_win64\
#    -llapack_win64_MT\
#    -lblas_win64_MT
