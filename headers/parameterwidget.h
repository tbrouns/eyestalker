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

#include "headers/constants.h"
#include "headers/structures.h"
#include "headers/sliderdouble.h"

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

    QLabel *AlphaFeaturesLabel;
    QLabel *AlphaMomentumLabel;
    QLabel *AlphaAverageLabel;
    QLabel *AlphaPositionLabel;
    QLabel *CannyBlurLevelLabel;
    QLabel *CannyKernelSizeLabel;
    QLabel *CannyThresholdLowLabel;
    QLabel *CannyThresholdHighLabel;
    QLabel *CurvatureOffsetLabel;
    QLabel *EdgeLengthFractionLabel;
    QLabel *EllipseFitNumberMaximumLabel;
    QLabel *EllipseFitErrorMaximumLabel;
    QLabel *GlintSizeLabel;
    QLabel *CircumferenceMinLabel;
    QLabel *CircumferenceMaxLabel;
    QLabel *AspectRatioMinLabel;
    QLabel *HaarOffsetLabel;
    QLabel *ThresholdCircumferenceLabel;
    QLabel *ThresholdAspectRatioLabel;
    QLabel *ThresholdDisplacementLabel;
    QLabel *ThresholdScoreLabel;

    QSlider *CannyBlurLevelSlider;
    QSlider *CannyKernelSizeSlider;
    QSlider *EllipseFitNumberMaximumSlider;
    QSlider *GlintSizeSlider;
    QSlider *HaarOffsetSlider;

    SliderDouble *AlphaFeaturesSlider;
    SliderDouble *AlphaMomentumSlider;
    SliderDouble *AlphaAverageSlider;
    SliderDouble *AlphaPositionSlider;
    SliderDouble *CannyThresholdLowSlider;
    SliderDouble *CannyThresholdHighSlider;
    SliderDouble *CurvatureOffsetSlider;
    SliderDouble *EdgeLengthFractionSlider;
    SliderDouble *EllipseFitErrorMaximumSlider;
    SliderDouble *CircumferenceMinSlider;
    SliderDouble *CircumferenceMaxSlider;
    SliderDouble *AspectRatioMinSlider;
    SliderDouble *ThresholdCircumferenceSlider;
    SliderDouble *ThresholdAspectRatioSlider;
    SliderDouble *ThresholdDisplacementSlider;
    SliderDouble *ThresholdScoreSlider;

    detectionParameters mDetectionParameters;

signals:

private slots:

    void setAlphaFeatures      (double);
    void setAlphaMomentum           (double);
    void setAlphaAverage            (double);
    void setAlphaPosition         (double);
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
    void setAspectRatioMin          (double);
    void setHaarOffset              (int);
    void setThresholdCircumference  (double);
    void setThresholdAspectRatio    (double);
    void setThresholdDisplacement   (double);
    void setThresholdScore          (double);

public slots:
};

#endif // PARAMETERWIDGET_H
