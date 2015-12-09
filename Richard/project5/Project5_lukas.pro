TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    planet.cpp \
    solarsystem.cpp \
    gaussian_random.cpp

HEADERS += \
    planet.h \
    solarsystem.h \
    gaussian_random.h

INCLUDEPATH += C:\Qt2\include
# LIBS += -L C:\Qt2
