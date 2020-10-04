TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        eno_advection.cpp \
        grid2d.cpp \
        main.cpp \
        setup.cpp

## Stuff for running in parallel
macx:
{
QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
}
macx:
{
QMAKE_LFLAGS += -lomp
}
macx:
{
LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib
}

HEADERS += \
    cf_2.h \
    eno_advection.h \
    grid2d.h \
    setup.h
