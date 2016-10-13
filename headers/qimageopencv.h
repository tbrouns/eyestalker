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

#include <headers/parameters.h>

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
    void drawAOI(QPixmap &img, int x0, int y0, int W, int H, QColor col);
    void loadImage(const cv::Mat& cvimage);
    void resizeImage();
    void setAOIError();
    void setEyeAOI(int x, int y, int w, int h);
    void setFindingCamera();
    void setFlashAOI(int x, int y, int w, int h);
    void setImage();
    void setSize(int, int);
    void setSpinner();

private:

    double aspectRatio;
    double imageScaleFactorX;
    double imageScaleFactorY;
    int eyeHghtAOI;
    int eyeWdthAOI;
    int eyeXPosAOI;
    int eyeYPosAOI;
    int flashHghtAOI;
    int flashWdthAOI;
    int flashXPosAOI;
    int flashYPosAOI;
    int imageHght;
    int imageScaledHght;
    int imageScaledWdth;
    int imageType;
    int imageWdth;
    int sizeH;
    int sizeW;
    int spinnerDegrees;
    QColor backgroundColour;
    QColor textColour;
    QPixmap image;
    QPixmap imageScaled;


protected:

    void mousePressEvent(QMouseEvent *event);

signals:

    void imageMouseClick(double, double);
    void updateImage();

};

#endif // QIMAGEOPENCV_H
