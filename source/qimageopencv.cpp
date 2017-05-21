//  EyeStalker: robust video-based eye tracking
//  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

//  EyeStalker is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  EyeStalker is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#include "qimageopencv.h"

QImageOpenCV::QImageOpenCV(QWidget *parent) : QLabel(parent)
{
    spinnerDegrees = 0;
    imageWdth  = 0;
    imageHght  = 0;
    widgetWdth = 0;
    widgetHght = 0;

    backgroundColour = QColor( 48,  47,  47);
    textColour       = QColor(177, 177, 177);

    SHOW_BEAD_AOI = false;

    this->setLineWidth(2);
    this->setAlignment(Qt::AlignCenter);
    this->setFrameStyle(QFrame::Panel | QFrame::Raised);
}

QImageOpenCV::~QImageOpenCV()
{

}

void QImageOpenCV::setSize(int W, int H)
{
    widgetWdth  = W;
    widgetHght  = H;
    aspectRatio = W / (double) H;
}

QSize QImageOpenCV::sizeHint() const
{
    return QSize(widgetWdth, widgetHght);
}

void QImageOpenCV::loadImage(const cv::Mat& cvimage)
{
    this->clear();

    // Copy input Mat
    const uchar *qImageBuffer = (const uchar*)cvimage.data;

    // Create QImage with same dimensions as input Mat
    QImage img(qImageBuffer, cvimage.cols, cvimage.rows, cvimage.step, QImage::Format_RGB888);

    image = QPixmap::fromImage(img.rgbSwapped());

    imageWdth = image.width();
    imageHght = image.height();

    resizeImage();
}

void QImageOpenCV::resizeImage()
{
    if (imageWdth > widgetWdth || imageHght > widgetHght)
    {
        imageScaled         = image.scaled(QSize(widgetWdth, widgetHght), Qt::KeepAspectRatio);
        imageWdthScaled     = imageScaled.width();
        imageHghtScaled     = imageScaled.height();
        imageScaleFactorX   = imageWdthScaled / (double) imageWdth;
        imageScaleFactorY   = imageHghtScaled / (double) imageHght;
    }
    else
    {
        imageScaled         = image;
        imageWdthScaled     = imageWdth;
        imageHghtScaled     = imageHght;
        imageScaleFactorX   = 1.0;
        imageScaleFactorY   = 1.0;
    }
}

void QImageOpenCV::drawAOI(QPixmap& img, AOIProperties mAOI, QColor col)
{
    if (imageWdthScaled > 0 && mAOI.wdth > 0 && mAOI.hght > 0)
    {
        mAOI.xPos = round(imageScaleFactorX * mAOI.xPos);
        mAOI.wdth = round(imageScaleFactorX * mAOI.wdth);
        mAOI.yPos = round(imageScaleFactorY * mAOI.yPos);
        mAOI.hght = round(imageScaleFactorY * mAOI.hght);

        QPainter painter(&img);

        int lineWidth = round(0.003 * (widgetWdth + widgetHght));

        painter.setPen(QPen(col, round(lineWidth)));
        painter.drawLine(QPoint(mAOI.xPos, mAOI.yPos),             QPoint(mAOI.xPos + mAOI.wdth, mAOI.yPos));
        painter.drawLine(QPoint(mAOI.xPos, mAOI.yPos + mAOI.hght), QPoint(mAOI.xPos + mAOI.wdth, mAOI.yPos + mAOI.hght));

        painter.setPen(QPen(col, round(lineWidth)));
        painter.drawLine(QPoint(mAOI.xPos,             mAOI.yPos), QPoint(mAOI.xPos,             mAOI.yPos + mAOI.hght));
        painter.drawLine(QPoint(mAOI.xPos + mAOI.wdth, mAOI.yPos), QPoint(mAOI.xPos + mAOI.wdth, mAOI.yPos + mAOI.hght));
    }
}

void QImageOpenCV::setImage()
{
    if (imageWdthScaled > 0)
    {
        QPixmap imageEdited = imageScaled;
        drawAOI(imageEdited,   eyeAOI, QColor(255,   0,   0));
        drawAOI(imageEdited, flashAOI, QColor(  0,   0, 255));
        if (SHOW_BEAD_AOI) { drawAOI(imageEdited,  beadAOI, QColor(  0, 255,   0)); }
        this->setPixmap(imageEdited);
    }
}

void QImageOpenCV::clearImage() { this->clear(); }

void QImageOpenCV::setAOIEye (AOIProperties eyeAOINew)  { eyeAOI  = eyeAOINew; }
void QImageOpenCV::setAOIBead(AOIProperties beadAOINew) { beadAOI = beadAOINew; }

void QImageOpenCV::setAOIFlash(AOIProperties flashAOINew)
{
    flashAOI = flashAOINew;
    flashAOI.xPos = flashAOI.xPos - Parameters::camAOI.xPos;
    flashAOI.yPos = flashAOI.yPos - Parameters::camAOI.yPos;

    if (flashAOI.xPos < 0)
    {
        if (flashAOI.xPos + flashAOI.wdth < 0) { flashAOI.wdth = 0; }
        else
        {
            flashAOI.wdth = flashAOI.xPos + flashAOI.wdth;
            flashAOI.xPos = 0;
        }
    }

    if (flashAOI.yPos < 0)
    {
        if (flashAOI.yPos + flashAOI.hght < 0) { flashAOI.hght = 0; }
        else
        {
            flashAOI.hght = flashAOI.yPos + flashAOI.hght;
            flashAOI.yPos = 0;
        }
    }
}

void QImageOpenCV::showAOIBead(bool state)
{
    SHOW_BEAD_AOI = state;
    this->setImage();
}

void QImageOpenCV::setAOIError()
{
    QPixmap pic(widgetWdth, widgetHght);
    pic.fill(backgroundColour);

    QRect rec(0, 0, round(0.9 * widgetWdth), round(0.5 * widgetHght));
    rec.moveCenter(rect().center());

    QPainter painter(&pic);

    QFont font = painter.font();
    font.setPointSize(round(0.032 * widgetWdth));
    font.setWeight(QFont::DemiBold);
    painter.setFont(font);

    painter.setPen(textColour);
    painter.drawText(rec, Qt::AlignCenter, QString("AREA OF INTEREST TOO SMALL"));
    this->setPixmap(pic);
}

void QImageOpenCV::setSpinner()
{
    QPixmap pic(widgetWdth, widgetHght);
    pic.fill(backgroundColour);

    QRect rec(0, 0, round(0.4 * widgetHght), round(0.4 * widgetHght));
    rec.moveCenter(rect().center());

    QPainter painter(&pic);

    QConicalGradient gradient;
    gradient.setCenter(rec.center());
    gradient.setAngle(round(spinnerDegrees / 16));
    gradient.setColorAt(0, backgroundColour);
    gradient.setColorAt(1, Qt::white);

    QPen pen(QBrush(gradient), 2);
    pen.setCapStyle(Qt::RoundCap);
    painter.setPen(pen);

    painter.drawArc(rec, spinnerDegrees, 16 * 360); // each step is 1/16th of a degree

    this->setPixmap(pic);

    spinnerDegrees = (spinnerDegrees % (360 * 16)) + 180;
}

void QImageOpenCV::mousePressEvent(QMouseEvent *event)
{
    int imageScaledXOffset = 0.5 * (widgetWdth - imageWdthScaled);
    int imageScaledYOffset = 0.5 * (widgetHght - imageHghtScaled);

    if (event->button() == Qt::LeftButton)
    {
        std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
        std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);

        double mouseXPos = (event->x()) - imageScaledXOffset;

        Parameters::eyeAOI.xPos = round(mouseXPos * (imageWdth / (double) imageWdthScaled) - 0.5 * Parameters::eyeAOI.wdth);
        if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth >= imageWdth)
        {   Parameters::eyeAOI.xPos = imageWdth - Parameters::eyeAOI.wdth; }
        else if (Parameters::eyeAOI.xPos < 0)
        {        Parameters::eyeAOI.xPos = 0; }
        Parameters::eyeAOIRatio.xPos = Parameters::eyeAOI.xPos / (double) imageWdth;

        double mouseYPos = (event->y()) - imageScaledYOffset;

        Parameters::eyeAOI.yPos = round(mouseYPos * (imageHght / (double) imageHghtScaled) - 0.5 * Parameters::eyeAOI.hght);
        if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght >= imageHght)
        {   Parameters::eyeAOI.yPos = imageHght - Parameters::eyeAOI.hght; }
        else if (Parameters::eyeAOI.yPos < 0)
        {        Parameters::eyeAOI.yPos = 0; }
        Parameters::eyeAOIRatio.yPos = Parameters::eyeAOI.yPos / (double) imageHght;

        setAOIEye(Parameters::eyeAOI);
        setImage();

        emit updateImage(-1);
    }
    else if (event->button() == Qt::RightButton)
    {
        std::lock_guard<std::mutex> AOICamLock (Parameters::AOICamMutex);
        std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);

        double mouseXPos = (event->x()) - imageScaledXOffset;

        Parameters::beadAOI.xPos = round(mouseXPos * (imageWdth / (double) imageWdthScaled) - 0.5 * Parameters::beadAOI.wdth);
        if (Parameters::beadAOI.xPos + Parameters::beadAOI.wdth >= imageWdth)
        {   Parameters::beadAOI.xPos = imageWdth - Parameters::beadAOI.wdth; }
        else if (Parameters::beadAOI.xPos < 0)
        {        Parameters::beadAOI.xPos = 0; }
        Parameters::beadAOIRatio.xPos = Parameters::beadAOI.xPos / (double) imageWdth;

        double mouseYPos = (event->y()) - imageScaledYOffset;

        Parameters::beadAOI.yPos = round(mouseYPos * (imageHght / (double) imageHghtScaled) - 0.5 * Parameters::beadAOI.hght);
        if (Parameters::beadAOI.yPos + Parameters::beadAOI.hght >= imageHght)
        {   Parameters::beadAOI.yPos = imageHght - Parameters::beadAOI.hght; }
        else if (Parameters::beadAOI.yPos < 0)
        {        Parameters::beadAOI.yPos = 0; }
        Parameters::beadAOIRatio.yPos = Parameters::beadAOI.yPos / (double) imageHght;

        setAOIBead(Parameters::beadAOI);
        setImage();

        emit updateImage(-1);
    }
}
