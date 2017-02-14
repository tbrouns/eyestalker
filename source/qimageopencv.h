//  Copyright (C) 2016  Terence Brouns

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#ifndef QIMAGEOPENCV_H
#define QIMAGEOPENCV_H

#include "parameters.h"

#include <iostream>     // std::ofstream
#include <mutex>

// Libraries

// OpenCV

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

// QT

#include <QColor>
#include <QFont>
#include <QGridLayout>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>
#include <QPen>
#include <QSize>
#include <QWidget>

class QImageOpenCV : public QLabel
{
    Q_OBJECT

public:

    explicit QImageOpenCV(int type = 0, QWidget *parent = 0);
    ~QImageOpenCV();

    QSize sizeHint() const;

    void clearImage();
    void loadImage(const cv::Mat& cvimage);
    void resizeImage();
    void setImage();

    void setFindingCamera();

    void setAOIError();
    void setAOIBead (AOIProperties beadAOINew);
    void setAOIEye  (AOIProperties eyeAOINew);
    void setAOIFlash(AOIProperties flashAOINew);
    void drawAOI(QPixmap &img, AOIProperties mAOI, QColor col);

    void showAOIBead(bool);

    void setSize(int, int);
    void setSpinner();

private:

    bool SHOW_BEAD_AOI;

    double aspectRatio;
    double imageScaleFactorX;
    double imageScaleFactorY;

    AOIProperties beadAOI;
    AOIProperties eyeAOI;
    AOIProperties flashAOI;

    int imageHght;
    int imageHghtScaled;
    int imageWdthScaled;
    int imageType;
    int imageWdth;
    int widgetHght;
    int widgetWdth;
    int spinnerDegrees;
    QColor backgroundColour;
    QColor textColour;
    QPixmap image;
    QPixmap imageScaled;


protected:

    void mousePressEvent(QMouseEvent *event);

signals:

    void imageMouseClick(double, double);
    void updateImage(int);

};

#endif // QIMAGEOPENCV_H
