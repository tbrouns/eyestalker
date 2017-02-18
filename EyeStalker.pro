#  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

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

QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = EyeStalker
TEMPLATE = app

CONFIG += c++11

unix {

# optimization

#QMAKE_CXXFLAGS += -O3
#QMAKE_CXXFLAGS += -DEIGEN_NO_DEBUG
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS   += -fopenmp

# default

QMAKE_LFLAGS += -Wl,-rpath,"'\$$ORIGIN'"

# libraries

# LibUSB
INCLUDEPATH += /usr/local/include/libusb-1.0/
LIBS        += -L/usr/local/lib -lusb-1.0

# OpenCV
INCLUDEPATH += "/usr/local/include/opencv2"
LIBS        += `pkg-config --cflags --libs opencv`

# Eigen
INCLUDEPATH += "/usr/local/include/eigen3"

# Boost
LIBS        += -L$$PWD/../../../../usr/local/lib/ -lboost_filesystem -lboost_system
INCLUDEPATH += $$PWD/../../../../usr/local/include
DEPENDPATH  += $$PWD/../../../../usr/local/include

# Qwt
LIBS        += -L$$PWD/../../../../usr/local/qwt-6.1.3/lib/ -lqwt
INCLUDEPATH += $$PWD/../../../../usr/local/qwt-6.1.3/include
DEPENDPATH  += $$PWD/../../../../usr/local/qwt-6.1.3/include

# UEye
LIBS        += -L$$PWD/../../../../../usr/lib/ -lueye_api
INCLUDEPATH += $$PWD/../../../../../usr/include
DEPENDPATH  += $$PWD/../../../../../usr/include
}

win32 {
# OpenCV
INCLUDEPATH += C:/libs/OpenCV-2.4.9/include
LIBS += -LC:\\libs\\OpenCV-2.4.9\\build-qt\\lib \
    -lopencv_calib3d249d \
    -lopencv_contrib249d \
    -lopencv_core249d \
    -lopencv_features2d249d \
    -lopencv_flann249d \
    -lopencv_gpu249d \
    -lopencv_highgui249d \
    -lopencv_imgproc249d \
    -lopencv_legacy249d \
    -lopencv_ml249d \
    -lopencv_nonfree249d \
    -lopencv_objdetect249d \
    -lopencv_ocl249d \
    -lopencv_photo249d \
    -lopencv_stitching249d \
    -lopencv_superres249d \
    -lopencv_ts249d \
    -lopencv_video249d \
    -lopencv_videostab249d

# Eigen
INCLUDEPATH += C://libs//eigen

# Boost
INCLUDEPATH += C://libs//boost

LIBS += -L$$PWD/../../../../libs/boost/stage/lib/ -lboost_filesystem-mgw49-mt-1_62
LIBS += -L$$PWD/../../../../libs/boost/stage/lib/ -lboost_system-mgw49-mt-1_62
LIBS += -L$$PWD/../../../../libs/boost/stage/lib/ -lboost_filesystem-mgw49-mt-d-1_62
LIBS += -L$$PWD/../../../../libs/boost/stage/lib/ -lboost_system-mgw49-mt-d-1_62

# UEye
INCLUDEPATH += C://libs//UEye//include
LIBS        += C://libs//UEye//Lib//uEye_api.lib
LIBS        += C://libs//UEye//Lib//uEye_api_64.lib
LIBS        += C://libs//UEye//Lib//uEye_tools.lib
LIBS        += C://libs//UEye//Lib//ueye_tools_64.lib
}

SOURCES += source/main.cpp\
        source/mainwindow.cpp \
    source/drawfunctions.cpp \
    source/mainwindowexperiment.cpp \
    source/mainwindowparameters.cpp \
    source/mainwindowoffline.cpp \
    source/parameters.cpp \
    source/qimageopencv.cpp \
    source/sliderdouble.cpp \
    source/ueyeopencv.cpp \
    source/confirmationwindow.cpp \
    source/parameterwidget.cpp \
    source/variablewidget.cpp \
    source/eyestalker.cpp \
    source/qwtplotwidget.cpp

HEADERS  += \
    source/confirmationwindow.h \
    source/constants.h \
    source/drawfunctions.h \
    source/mainwindow.h \
    source/parameters.h \
    source/qimageopencv.h \
    source/sliderdouble.h \
    source/structures.h \
    source/ueyeopencv.h \
    source/parameterwidget.h \
    source/variablewidget.h \
    source/eyestalker.h \
    source/qwtplotwidget.h

RESOURCES += \
    resources/qdarkstyle/style.qrc
