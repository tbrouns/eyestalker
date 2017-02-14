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
    QLabel      *CircumferenceOffsetLabel;
    SliderDouble*CircumferenceOffsetSlider;
    QLabel      *CurvatureOffsetLabel;
    SliderDouble*CurvatureOffsetSlider;
    QLabel      *EdgeLengthFractionLabel;
    SliderDouble*EdgeLengthFractionSlider;
    QLabel      *EllipseFitErrorMaximumLabel;
    SliderDouble*EllipseFitErrorMaximumSlider;
    QLabel      *EllipseFitNumberMaximumLabel;
    QSlider     *EllipseFitNumberMaximumSlider;
    QLabel      *GlintSizeLabel;
    QSlider     *GlintSizeSlider;
    QLabel      *ThresholdAspectRatioLabel;
    SliderDouble*ThresholdAspectRatioSlider;
    QLabel      *ThresholdCircumferenceLabel;
    SliderDouble*ThresholdCircumferenceSlider;
    QLabel      *ThresholdDisplacementLabel;
    SliderDouble*ThresholdDisplacementSlider;
    QLabel      *ThresholdScoreLabel;
    SliderDouble*ThresholdScoreSlider;

    detectionParameters mDetectionParameters;

signals:

private slots:

    void setAlphaFeatures           (double);
    void setAlphaAverages            (double);
    void setAlphaPosition           (double);
    void setAlphaCertainty          (double);
    void setCannyBlurLevel          (int);
    void setCannyKernelSize         (int);
    void setCannyThresholdLow       (double);
    void setCannyThresholdHigh      (double);
    void setCurvatureOffset         (double);
    void setEdgeLengthFraction      (double);
    void setEllipseFitNumberMaximum (int);
    void setEllipseFitErrorMaximum  (double);
    void setGlintSize               (int);
    void setCircumferenceMin        (double);
    void setCircumferenceMax        (double);
    void setCircumferenceOffset     (double);
    void setAspectRatioMin          (double);
    void setThresholdCircumference  (double);
    void setThresholdAspectRatio    (double);
    void setThresholdDisplacement   (double);
    void setThresholdScore          (double);

public slots:
};

#endif // PARAMETERWIDGET_H
