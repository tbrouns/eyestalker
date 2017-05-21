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

#ifndef PARAMETERWIDGET_H
#define PARAMETERWIDGET_H

#include <QGridLayout>
#include <QLabel>
#include <QSlider>
#include <QObject>
#include <QSettings>
#include <QScrollArea>
#include <QVBoxLayout>
#include <QWidget>

#include "constants.h"
#include "structures.h"
#include "sliderdouble.h"

class ParameterWidget : public QWidget
{
    Q_OBJECT

public:

    void reset();
    detectionParameters getStructure();
    bool getState();
    void setState(bool);
    void setStructure(detectionParameters);

    explicit ParameterWidget(QWidget *parent = 0);
    ~ParameterWidget();

private:

    QLabel      *GainAveragesLabel;
    SliderDouble*GainAveragesSlider;
    QLabel      *GainCertaintyLabel;
    SliderDouble*GainCertaintySlider;
    QLabel      *GainAppearanceLabel;
    SliderDouble*GainAppearanceSlider;
    QLabel      *GainPositionLabel;
    SliderDouble*GainPositionSlider;
    QLabel      *AspectRatioMinLabel;
    SliderDouble*AspectRatioMinSlider;
    QLabel      *CannyBlurLevelLabel;
    QSlider     *CannyBlurLevelSlider;
    QLabel      *CannyKernelSizeLabel;
    QSlider     *CannyKernelSizeSlider;
    QLabel      *CannyThresholdHighLabel;
    SliderDouble*CannyThresholdHighSlider;
    QLabel      *CannyThresholdLowLabel;
    SliderDouble*CannyThresholdLowSlider;
    QLabel      *CircumferenceMaxLabel;
    SliderDouble*CircumferenceMaxSlider;
    QLabel      *CircumferenceMinLabel;
    SliderDouble*CircumferenceMinSlider;
    QLabel      *CurvatureOffsetLabel;
    SliderDouble*CurvatureOffsetSlider;
    QLabel      *FitEdgeFractionLabel;
    SliderDouble*FitEdgeFractionSlider;
    QLabel      *FitEdgeMaximumLabel;
    QSlider     *FitEdgeMaximumSlider;
    QLabel      *FitMaximumLabel;
    QSlider     *FitMaximumSlider;
    QLabel      *GlintSizeLabel;
    QSlider     *GlintSizeSlider;
    QLabel      *ThresholdAspectRatioUpperLabel;
    SliderDouble*ThresholdAspectRatioUpperSlider;
    QLabel      *ThresholdCircumferenceUpperLabel;
    SliderDouble*ThresholdCircumferenceUpperSlider;
    QLabel      *ThresholdDisplacementUpperLabel;
    SliderDouble*ThresholdDisplacementUpperSlider;
    QLabel      *ThresholdAspectRatioLowerLabel;
    SliderDouble*ThresholdAspectRatioLowerSlider;
    QLabel      *ThresholdCircumferenceLowerLabel;
    SliderDouble*ThresholdCircumferenceLowerSlider;
    QLabel      *ThresholdDisplacementLowerLabel;
    SliderDouble*ThresholdDisplacementLowerSlider;
    QLabel      *ThresholdScoreEdgeLabel;
    SliderDouble*ThresholdScoreEdgeSlider;
    QLabel      *ThresholdScoreFitLabel;
    SliderDouble*ThresholdScoreFitSlider;
    QLabel      *ThresholdScoreDiffEdgeLabel;
    SliderDouble*ThresholdScoreDiffEdgeSlider;
    QLabel      *ThresholdScoreDiffFitLabel;
    SliderDouble*ThresholdScoreDiffFitSlider;
    QLabel      *ThresholdFitErrorLabel;
    SliderDouble*ThresholdFitErrorSlider;
    QLabel      *WindowLengthEdgeLabel;
    QSlider     *WindowLengthEdgeSlider;

    detectionParameters mDetectionParameters;

signals:

private slots:

    void setGainAppearance                    (double);
    void setGainAverages                    (double);
    void setGainPosition                    (double);
    void setGainCertainty                   (double);
    void setCannyBlurLevel                  (int);
    void setCannyKernelSize                 (int);
    void setCannyThresholdLow               (double);
    void setCannyThresholdHigh              (double);
    void setCurvatureOffset                 (double);
    void setFitEdgeFraction                 (double);
    void setFitEdgeMaximum                  (int);
    void setFitMaximum                      (int);
    void setThresholdFitError               (double);
    void setGlintSize                       (int);
    void setCircumferenceMin                (double);
    void setCircumferenceMax                (double);
    void setAspectRatioMin                  (double);
    void setThresholdCircumferenceLower     (double);
    void setThresholdCircumferenceUpper     (double);
    void setThresholdAspectRatioLower       (double);
    void setThresholdAspectRatioUpper       (double);
    void setThresholdDisplacementLower      (double);
    void setThresholdDisplacementUpper      (double);
    void setThresholdScoreEdge              (double);
    void setThresholdScoreFit               (double);
    void setThresholdScoreDiffEdge          (double);
    void setThresholdScoreDiffFit           (double);
    void setWindowLengthEdge                (int);

public slots:
};

#endif // PARAMETERWIDGET_H
