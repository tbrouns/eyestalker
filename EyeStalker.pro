#  Copyright (C) 2016  Terence Brouns

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>


QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = EyeTrackingTB
TEMPLATE = app

QT       += core gui

CONFIG += c++11

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

QMAKE_LFLAGS += -Wl,-rpath,"'\$$ORIGIN'"

# LibUSB
 INCLUDEPATH += /usr/local/include/libusb-1.0/
 LIBS += -L/usr/local/lib -lusb-1.0

# OpenCV
INCLUDEPATH += "/usr/local/include/opencv2"
LIBS += `pkg-config --cflags --libs opencv`

# Eigen
INCLUDEPATH += "/usr/local/include/eigen3"

# Boost
LIBS += -L$$PWD/../../../../usr/local/lib/ -lboost_filesystem
INCLUDEPATH += $$PWD/../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../usr/local/include

# Qwt
unix:!macx: LIBS += -L$$PWD/../../../../usr/local/qwt-6.1.2/lib/ -lqwt
INCLUDEPATH += $$PWD/../../../../usr/local/qwt-6.1.2/include
DEPENDPATH += $$PWD/../../../../usr/local/qwt-6.1.2/include

# UEye
LIBS += -L$$PWD/../../../../../usr/lib/ -lueye_api
INCLUDEPATH += $$PWD/../../../../../usr/include
DEPENDPATH += $$PWD/../../../../../usr/include


SOURCES += main.cpp\
        mainwindow.cpp \
    source/drawfunctions.cpp \
    source/eyetracking.cpp \
    source/mainwindowexperiment.cpp \
    source/mainwindowparameters.cpp \
    source/mainwindowreview.cpp \
    source/parameters.cpp \
    source/qimageopencv.cpp \
    source/sliderdouble.cpp \
    source/startupwindow.cpp \
    source/ueyeopencv.cpp \
    source/confirmationwindow.cpp

HEADERS  += \
    source/constants.h \
    source/drawfunctions.h \
    source/eyetracking.h \
    source/parameters.h \
    source/qimageopencv.h \
    source/sliderdouble.h \
    source/startupwindow.h \
    source/structures.h \
    source/ueyeopencv.h \
    source/confirmationwindow.h

RESOURCES += \
    resources/qdarkstyle/style.qrc
