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

    QLabel      *AlphaAveragesLabel;
    SliderDouble*AlphaAveragesSlider;
    QLabel      *AlphaCertaintyLabel;
    SliderDouble*AlphaCertaintySlider;
    QLabel      *AlphaFeaturesLabel;
    SliderDouble*AlphaFeaturesSlider;
    QLabel      *AlphaPositionLabel;
    SliderDouble*AlphaPositionSlider;
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
    QLabel      *ThresholdAspectRatioLabel;
    SliderDouble*ThresholdAspectRatioSlider;
    QLabel      *ThresholdCircumferenceLabel;
    SliderDouble*ThresholdCircumferenceSlider;
    QLabel      *ThresholdDisplacementLabel;
    SliderDouble*ThresholdDisplacementSlider;
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

    void setAlphaFeatures           (double);
    void setAlphaAverages           (double);
    void setAlphaPosition           (double);
    void setAlphaCertainty          (double);
    void setCannyBlurLevel          (int);
    void setCannyKernelSize         (int);
    void setCannyThresholdLow       (double);
    void setCannyThresholdHigh      (double);
    void setCurvatureOffset         (double);
    void setFitEdgeFraction         (double);
    void setFitEdgeMaximum          (int);
    void setFitMaximum              (int);
    void setThresholdFitError       (double);
    void setGlintSize               (int);
    void setCircumferenceMin        (double);
    void setCircumferenceMax        (double);
    void setAspectRatioMin          (double);
    void setThresholdCircumference  (double);
    void setThresholdAspectRatio    (double);
    void setThresholdDisplacement   (double);
    void setThresholdScoreEdge      (double);
    void setThresholdScoreFit       (double);
    void setThresholdScoreDiffEdge  (double);
    void setThresholdScoreDiffFit   (double);
    void setWindowLengthEdge        (int);

public slots:
};

#endif // PARAMETERWIDGET_H
