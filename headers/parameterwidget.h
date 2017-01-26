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

    QLabel *AlphaMiscellaneousLabel;
    QLabel *AlphaMomentumLabel;
    QLabel *AlphaAverageLabel;
    QLabel *AlphaPredictionLabel;
    QLabel *CannyBlurLevelLabel;
    QLabel *CannyKernelSizeLabel;
    QLabel *CannyThresholdLowLabel;
    QLabel *CannyThresholdHighLabel;
    QLabel *CurvatureFactorLabel;
    QLabel *CurvatureOffsetLabel;
    QLabel *EdgeIntensityOffsetLabel;
    QLabel *EdgeLengthMinimumLabel;
    QLabel *EllipseFitNumberMaximumLabel;
    QLabel *EllipseFitErrorMaximumLabel;
    QLabel *GlintSizeLabel;
    QLabel *CircumferenceMinLabel;
    QLabel *CircumferenceMaxLabel;
    QLabel *AspectRatioMinLabel;
    QLabel *HaarOffsetLabel;
    QLabel *ThresholdCircumferenceLabel;
    QLabel *ThresholdAspectRatioLabel;

    QSlider *CannyBlurLevelSlider;
    QSlider *CannyKernelSizeSlider;
    QSlider *CannyThresholdLowSlider;
    QSlider *CannyThresholdHighSlider;
    QSlider *EllipseFitNumberMaximumSlider;
    QSlider *GlintSizeSlider;
    QSlider *HaarOffsetSlider;

    SliderDouble *AlphaMiscellaneousSlider;
    SliderDouble *AlphaMomentumSlider;
    SliderDouble *AlphaAverageSlider;
    SliderDouble *AlphaPredictionSlider;
    SliderDouble *CurvatureFactorSlider;
    SliderDouble *CurvatureOffsetSlider;
    SliderDouble *EdgeIntensityOffsetSlider;
    SliderDouble *EdgeLengthMinimumSlider;
    SliderDouble *EllipseFitErrorMaximumSlider;
    SliderDouble *CircumferenceMinSlider;
    SliderDouble *CircumferenceMaxSlider;
    SliderDouble *AspectRatioMinSlider;
    SliderDouble *ThresholdCircumferenceSlider;
    SliderDouble *ThresholdAspectRatioSlider;

    detectionParameters mDetectionParameters;

signals:

private slots:

    void setAlphaMiscellaneous      (double);
    void setAlphaMomentum           (double);
    void setAlphaAverage            (double);
    void setAlphaPrediction         (double);
    void setCannyBlurLevel          (int);
    void setCannyKernelSize         (int);
    void setCannyThresholdLow       (int);
    void setCannyThresholdHigh      (int);
    void setCurvatureFactor         (double);
    void setCurvatureOffset         (double);
    void setEdgeIntensityOffset     (double);
    void setEdgeLengthMinimum       (double);
    void setEllipseFitNumberMaximum (int);
    void setEllipseFitErrorMaximum  (double);
    void setGlintSize               (int);
    void setCircumferenceMin   (double);
    void setCircumferenceMax   (double);
    void setAspectRatioMin     (double);
    void setHaarOffset         (int);
    void setThresholdCircumference  (double);
    void setThresholdAspectRatio    (double);

public slots:
};

#endif // PARAMETERWIDGET_H
