#include "parameterwidget.h"

ParameterWidget::ParameterWidget(QWidget *parent) : QWidget(parent)
{
    // Averages and Limits/Thresholds

    QLabel *ParametersTextBox = new QLabel;
    ParametersTextBox->setText("<b>Eye-tracking parameters</b>");
    ParametersTextBox->setAlignment(Qt::AlignCenter);

    // Pupil circumference

    QLabel *CircumferenceMinTextBox = new QLabel;
    CircumferenceMinTextBox->setText("<b>Circumference min:</b>");

    CircumferenceMinLabel  = new QLabel;
    CircumferenceMinSlider = new SliderDouble;
    CircumferenceMinSlider->setPrecision(1);
    CircumferenceMinSlider->setDoubleRange(1, 500);
    CircumferenceMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CircumferenceMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCircumferenceMin(double)));

    QLabel *CircumferenceMaxTextBox = new QLabel;
    CircumferenceMaxTextBox->setText("<b>Circumference max:</b>");

    CircumferenceMaxLabel  = new QLabel;
    CircumferenceMaxSlider = new SliderDouble;
    CircumferenceMaxSlider->setPrecision(1);
    CircumferenceMaxSlider->setDoubleRange(1, 500);
    CircumferenceMaxSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CircumferenceMaxSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCircumferenceMax(double)));

    // Pupil aspect ratio

    QLabel *AspectRatioMinTextBox = new QLabel;
    AspectRatioMinTextBox->setText("<b>Aspect ratio min:</b>");

    AspectRatioMinLabel  = new QLabel;
    AspectRatioMinSlider = new SliderDouble;
    AspectRatioMinSlider->setPrecision(2);
    AspectRatioMinSlider->setDoubleRange(0.0, 1.0);
    AspectRatioMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(AspectRatioMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAspectRatioMin(double)));

    // Sliders for canny edge parameters

    QLabel *CannyEdgeTextBox = new QLabel;
    CannyEdgeTextBox->setText("<b>Canny edge detection</b>");
    CannyEdgeTextBox->setAlignment(Qt::AlignCenter);

    QLabel *CannyThresholdLowTextBox = new QLabel;
    CannyThresholdLowTextBox->setText("<b>Low threshold:</b>");

    CannyThresholdLowLabel  = new QLabel;
    CannyThresholdLowSlider = new SliderDouble;
    CannyThresholdLowSlider->setPrecision(1);
    CannyThresholdLowSlider->setDoubleRange(0, 1000);
    CannyThresholdLowSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyThresholdLowSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCannyThresholdLow(double)));

    QLabel *CannyThresholdHighTextBox = new QLabel;
    CannyThresholdHighTextBox->setText("<b>High threshold:</b>");

    CannyThresholdHighLabel  = new QLabel;
    CannyThresholdHighSlider = new SliderDouble;
    CannyThresholdHighSlider->setPrecision(1);
    CannyThresholdHighSlider->setDoubleRange(0, 1000);
    CannyThresholdHighSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyThresholdHighSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCannyThresholdHigh(double)));

    QLabel *CannyBlurLevelTextBox = new QLabel;
    CannyBlurLevelTextBox->setText("<b>Blur level:</b>");

    CannyBlurLevelLabel  = new QLabel;
    CannyBlurLevelSlider = new QSlider;
    CannyBlurLevelSlider->setRange(0, 10);
    CannyBlurLevelSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyBlurLevelSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyBlurLevel(int)));

    QLabel *CannyKernelSizeTextBox = new QLabel;
    CannyKernelSizeTextBox->setText("<b>Kernel size:</b>");

    CannyKernelSizeLabel  = new QLabel;
    CannyKernelSizeSlider = new QSlider;
    CannyKernelSizeSlider->setRange(2, 4);
    CannyKernelSizeSlider->setOrientation(Qt::Horizontal);
    CannyKernelSizeSlider->setSingleStep(1);
    QObject::connect(CannyKernelSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyKernelSize(int)));

    // Learning rate parameters

    QLabel *LearningRatesTextBox = new QLabel;
    LearningRatesTextBox->setText("<b>Learning rates</b>");
    LearningRatesTextBox->setAlignment(Qt::AlignCenter);

    QLabel *AlphaAveragesTextBox  = new QLabel;
    AlphaAveragesTextBox->setText("<b>Averages:</b>");

    AlphaAveragesLabel  = new QLabel;
    AlphaAveragesSlider = new SliderDouble;
    AlphaAveragesSlider->setPrecision(3);
    AlphaAveragesSlider->setDoubleRange(0, 0.1);
    AlphaAveragesSlider->setOrientation(Qt::Horizontal);
    QObject::connect(AlphaAveragesSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaAverages(double)));

    QLabel *AlphaPositionTextBox = new QLabel;
    AlphaPositionTextBox->setText("<b>Position:</b>");

    AlphaPositionLabel  = new QLabel;
    AlphaPositionSlider = new SliderDouble;
    AlphaPositionSlider->setPrecision(2);
    AlphaPositionSlider->setDoubleRange(0, 1.0);
    AlphaPositionSlider->setOrientation(Qt::Horizontal);
    QObject::connect(AlphaPositionSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaPosition(double)));

    QLabel *AlphaFeaturesTextBox = new QLabel;
    AlphaFeaturesTextBox->setText("<b>Features:</b>");

    AlphaFeaturesLabel  = new QLabel;
    AlphaFeaturesSlider = new SliderDouble;
    AlphaFeaturesSlider->setPrecision(2);
    AlphaFeaturesSlider->setDoubleRange(0, 1.0);
    AlphaFeaturesSlider->setOrientation(Qt::Horizontal);
    QObject::connect(AlphaFeaturesSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaFeatures(double)));

    QLabel *AlphaCertaintyTextBox = new QLabel;
    AlphaCertaintyTextBox->setText("<b>Certainty:</b>");

    AlphaCertaintyLabel  = new QLabel;
    AlphaCertaintySlider = new SliderDouble;
    AlphaCertaintySlider->setPrecision(2);
    AlphaCertaintySlider->setDoubleRange(0, 1.0);
    AlphaCertaintySlider->setOrientation(Qt::Horizontal);
    QObject::connect(AlphaCertaintySlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaCertainty(double)));

    // Threshold parameters

    QLabel *ThresholdCircumferenceTextBox = new QLabel;
    ThresholdCircumferenceTextBox->setText("<b>Circumference change</b>");

    QLabel *ThresholdAspectRatioTextBox = new QLabel;
    ThresholdAspectRatioTextBox->setText("<b>Aspect ratio change</b>");

    QLabel *ThresholdDisplacementTextBox = new QLabel;
    ThresholdDisplacementTextBox->setText("<b>Displacement change</b>");

    ThresholdAspectRatioLabel  = new QLabel;
    ThresholdAspectRatioSlider = new SliderDouble;
    ThresholdAspectRatioSlider->setPrecision(2);
    ThresholdAspectRatioSlider->setDoubleRange(0.0, 0.5);
    ThresholdAspectRatioSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdAspectRatioSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdAspectRatio(double)));

    ThresholdCircumferenceLabel  = new QLabel;
    ThresholdCircumferenceSlider = new SliderDouble;
    ThresholdCircumferenceSlider->setPrecision(2);
    ThresholdCircumferenceSlider->setDoubleRange(0.0, 0.5);
    ThresholdCircumferenceSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdCircumferenceSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdCircumference(double)));

    ThresholdDisplacementLabel  = new QLabel;
    ThresholdDisplacementSlider = new SliderDouble;
    ThresholdDisplacementSlider->setPrecision(1);
    ThresholdDisplacementSlider->setDoubleRange(0, 25);
    ThresholdDisplacementSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdDisplacementSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdDisplacement(double)));

    QLabel *ThresholdScoreEdgeTextBox = new QLabel;
    ThresholdScoreEdgeTextBox->setText("<b>Edge score</b>");

    ThresholdScoreEdgeLabel  = new QLabel;
    ThresholdScoreEdgeSlider = new SliderDouble;
    ThresholdScoreEdgeSlider->setPrecision(2);
    ThresholdScoreEdgeSlider->setDoubleRange(0, 0.5);
    ThresholdScoreEdgeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreEdgeSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreEdge(double)));

    QLabel *ThresholdScoreFitTextBox = new QLabel;
    ThresholdScoreFitTextBox->setText("<b>Fit score</b>");

    ThresholdScoreFitLabel  = new QLabel;
    ThresholdScoreFitSlider = new SliderDouble;
    ThresholdScoreFitSlider->setPrecision(2);
    ThresholdScoreFitSlider->setDoubleRange(0, 0.5);
    ThresholdScoreFitSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreFitSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreFit(double)));

    QLabel *ThresholdScoreDiffEdgeTextBox = new QLabel;
    ThresholdScoreDiffEdgeTextBox->setText("<b>Edge score difference</b>");

    ThresholdScoreDiffEdgeLabel  = new QLabel;
    ThresholdScoreDiffEdgeSlider = new SliderDouble;
    ThresholdScoreDiffEdgeSlider->setPrecision(2);
    ThresholdScoreDiffEdgeSlider->setDoubleRange(0, 1.0);
    ThresholdScoreDiffEdgeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreDiffEdgeSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreDiffEdge(double)));

    QLabel *ThresholdScoreDiffFitTextBox = new QLabel;
    ThresholdScoreDiffFitTextBox->setText("<b>Fit score difference</b>");

    ThresholdScoreDiffFitLabel  = new QLabel;
    ThresholdScoreDiffFitSlider = new SliderDouble;
    ThresholdScoreDiffFitSlider->setPrecision(2);
    ThresholdScoreDiffFitSlider->setDoubleRange(0, 0.5);
    ThresholdScoreDiffFitSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreDiffFitSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreDiffFit(double)));

    QLabel *ThresholdFitErrorTextBox = new QLabel;
    ThresholdFitErrorTextBox->setText("<b>Fit error:</b>");

    ThresholdFitErrorLabel  = new QLabel;
    ThresholdFitErrorSlider = new SliderDouble;
    ThresholdFitErrorSlider->setPrecision(2);
    ThresholdFitErrorSlider->setDoubleRange(0, 1.0);
    ThresholdFitErrorSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdFitErrorSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdFitError(double)));

    // Miscellaneous parameters

    QLabel *MiscParametersTextBox = new QLabel;
    MiscParametersTextBox->setText("<b>Miscellaneous</b>");
    MiscParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *GlintSizeTextBox = new QLabel;
    GlintSizeTextBox->setText("<b>Glint size</b>");

    GlintSizeLabel  = new QLabel;
    GlintSizeSlider = new QSlider;
    GlintSizeSlider->setRange(0, 10);
    GlintSizeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GlintSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setGlintSize(int)));

    QLabel *CurvatureOffsetTextBox = new QLabel;
    CurvatureOffsetTextBox->setText("<b>Curvature offset:</b>");

    CurvatureOffsetLabel  = new QLabel;
    CurvatureOffsetSlider = new SliderDouble;
    CurvatureOffsetSlider->setPrecision(1);
    CurvatureOffsetSlider->setDoubleRange(0, 40);
    CurvatureOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureOffset(double)));

    QLabel *WindowLengthEdgeTextBox = new QLabel;
    WindowLengthEdgeTextBox->setText("<b>Edge window length:</b>");

    WindowLengthEdgeLabel  = new QLabel;
    WindowLengthEdgeSlider = new QSlider;
    WindowLengthEdgeSlider->setRange(5, 11);
    WindowLengthEdgeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(WindowLengthEdgeSlider, SIGNAL(valueChanged(int)), this, SLOT(setWindowLengthEdge(int)));

    QLabel *FitEdgeFractionTextBox = new QLabel;
    FitEdgeFractionTextBox->setText("<b>Fit error edge fraction:");

    FitEdgeFractionLabel  = new QLabel;
    FitEdgeFractionSlider = new SliderDouble;
    FitEdgeFractionSlider->setPrecision(2);
    FitEdgeFractionSlider->setDoubleRange(0.0, 0.2);
    FitEdgeFractionSlider->setOrientation(Qt::Horizontal);
    QObject::connect(FitEdgeFractionSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setFitEdgeFraction(double)));

    QLabel *FitEdgeMaximumTextBox = new QLabel;
    FitEdgeMaximumTextBox->setText("<b>Maximum number of edges:</b>");

    FitEdgeMaximumLabel  = new QLabel;
    FitEdgeMaximumSlider = new QSlider;
    FitEdgeMaximumSlider->setRange(1, 7);
    FitEdgeMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(FitEdgeMaximumSlider, SIGNAL(valueChanged(int)), this, SLOT(setFitEdgeMaximum(int)));

    QLabel *FitMaximumTextBox = new QLabel;
    FitMaximumTextBox->setText("<b>Maximum number of fits:</b>");

    FitMaximumLabel  = new QLabel;
    FitMaximumSlider = new QSlider;
    FitMaximumSlider->setRange(1, 10);
    FitMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(FitMaximumSlider, SIGNAL(valueChanged(int)), this, SLOT(setFitMaximum(int)));

    QLabel *TitleLimitTextBox  = new QLabel;
    QLabel *TitleCannyTextBox  = new QLabel;
    QLabel *TitleLearnTextBox  = new QLabel;
    QLabel *TitleChangeTextBox = new QLabel;
    QLabel *TitleMiscTextBox   = new QLabel;

    TitleLimitTextBox ->setText("<b>Variable limits</b>");
    TitleCannyTextBox ->setText("<b>Canny edge detection</b>");
    TitleLearnTextBox ->setText("<b>Learning rates</b>");
    TitleChangeTextBox->setText("<b>Thresholds</b>");
    TitleMiscTextBox  ->setText("<b>Miscellaneous</b>");

    QWidget *MainWidget = new QWidget;
    QGridLayout *MainLayout = new QGridLayout(MainWidget);

    // Text labels

    MainLayout->addWidget(CircumferenceMaxTextBox,          1, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CircumferenceMinTextBox,          2, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AspectRatioMinTextBox,            3, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(CannyThresholdHighTextBox,        6, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyThresholdLowTextBox,         7, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyKernelSizeTextBox,           8, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyBlurLevelTextBox,            9, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(AlphaAveragesTextBox,             11, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaPositionTextBox,             12, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaFeaturesTextBox,             13, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaCertaintyTextBox,            15, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(ThresholdAspectRatioTextBox,      18, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdCircumferenceTextBox,    19, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdDisplacementTextBox,     20, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreEdgeTextBox,        21, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreFitTextBox,         22, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreDiffEdgeTextBox,    23, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreDiffFitTextBox,     24, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdFitErrorTextBox,         25, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(GlintSizeTextBox,                 27, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CurvatureOffsetTextBox,           28, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(WindowLengthEdgeTextBox,          29, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(FitEdgeFractionTextBox,           30, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(FitEdgeMaximumTextBox,            31, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(FitMaximumTextBox,                32, 0, 1, 1, Qt::AlignRight);

    // Sliders and titles

    MainLayout->addWidget(TitleLimitTextBox,                0, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(CircumferenceMaxSlider,           1, 1);
    MainLayout->addWidget(CircumferenceMinSlider,           2, 1);
    MainLayout->addWidget(AspectRatioMinSlider,             3, 1);

    MainLayout->addWidget(TitleCannyTextBox,                5, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(CannyThresholdHighSlider,         6, 1);
    MainLayout->addWidget(CannyThresholdLowSlider,          7, 1);
    MainLayout->addWidget(CannyKernelSizeSlider,            8, 1);
    MainLayout->addWidget(CannyBlurLevelSlider,             9, 1);

    MainLayout->addWidget(TitleLearnTextBox,               10, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(AlphaAveragesSlider,             11, 1);
    MainLayout->addWidget(AlphaPositionSlider,             12, 1);
    MainLayout->addWidget(AlphaFeaturesSlider,             13, 1);
    MainLayout->addWidget(AlphaCertaintySlider,            15, 1);

    MainLayout->addWidget(TitleChangeTextBox,              17, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(ThresholdAspectRatioSlider,      18, 1);
    MainLayout->addWidget(ThresholdCircumferenceSlider,    19, 1);
    MainLayout->addWidget(ThresholdDisplacementSlider,     20, 1);
    MainLayout->addWidget(ThresholdScoreEdgeSlider,        21, 1);
    MainLayout->addWidget(ThresholdScoreFitSlider,         22, 1);

    MainLayout->addWidget(ThresholdScoreDiffEdgeSlider,    23, 1);
    MainLayout->addWidget(ThresholdScoreDiffFitSlider,     24, 1);
    MainLayout->addWidget(ThresholdFitErrorSlider,         25, 1);

    MainLayout->addWidget(TitleMiscTextBox,                26, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(GlintSizeSlider,                 27, 1);
    MainLayout->addWidget(CurvatureOffsetSlider,           28, 1);
    MainLayout->addWidget(WindowLengthEdgeSlider,          29, 1);
    MainLayout->addWidget(FitEdgeFractionSlider,           30, 1);
    MainLayout->addWidget(FitEdgeMaximumSlider,            31, 1);
    MainLayout->addWidget(FitMaximumSlider,                32, 1);

    // Value labels

    MainLayout->addWidget(CircumferenceMaxLabel,         1, 2);
    MainLayout->addWidget(CircumferenceMinLabel,         2, 2);
    MainLayout->addWidget(AspectRatioMinLabel,           3, 2);

    MainLayout->addWidget(CannyThresholdHighLabel,       6, 2);
    MainLayout->addWidget(CannyThresholdLowLabel,        7, 2);
    MainLayout->addWidget(CannyKernelSizeLabel,          8, 2);
    MainLayout->addWidget(CannyBlurLevelLabel,           9, 2);

    MainLayout->addWidget(AlphaAveragesLabel,            11, 2);
    MainLayout->addWidget(AlphaPositionLabel,            12, 2);
    MainLayout->addWidget(AlphaFeaturesLabel,            13, 2);
    MainLayout->addWidget(AlphaCertaintyLabel,           15, 2);

    MainLayout->addWidget(ThresholdAspectRatioLabel,     18, 2);
    MainLayout->addWidget(ThresholdCircumferenceLabel,   19, 2);
    MainLayout->addWidget(ThresholdDisplacementLabel,    20, 2);
    MainLayout->addWidget(ThresholdScoreEdgeLabel,       21, 2);
    MainLayout->addWidget(ThresholdScoreFitLabel,        22, 2);
    MainLayout->addWidget(ThresholdScoreDiffEdgeLabel,   23, 2);
    MainLayout->addWidget(ThresholdScoreDiffFitLabel,    24, 2);
    MainLayout->addWidget(ThresholdFitErrorLabel,        25, 2);

    MainLayout->addWidget(GlintSizeLabel,                27, 2);
    MainLayout->addWidget(CurvatureOffsetLabel,          28, 2);
    MainLayout->addWidget(WindowLengthEdgeLabel,         29, 2);
    MainLayout->addWidget(FitEdgeFractionLabel,          30, 2);
    MainLayout->addWidget(FitEdgeMaximumLabel,           31, 2);
    MainLayout->addWidget(FitMaximumLabel,               32, 2);

    MainLayout->setColumnStretch(0,1);
    MainLayout->setColumnStretch(1,3);
    MainLayout->setColumnStretch(2,1);

    setLayout(MainLayout);
}

ParameterWidget::~ParameterWidget()
{

}

detectionParameters ParameterWidget::getStructure()
{
    return mDetectionParameters;
}

bool ParameterWidget::getState()
{
  return mDetectionParameters.DETECTION_ON;
}

void ParameterWidget::setState(bool state)
{
    mDetectionParameters.DETECTION_ON = state;
}

void ParameterWidget::setStructure(detectionParameters mEyePropertiesParametersNew)
{
    mDetectionParameters = mEyePropertiesParametersNew;
    this->reset();
}

void ParameterWidget::reset()
{
    // Size and shape limits

    CircumferenceMaxSlider->setDoubleValue(mDetectionParameters.thresholdCircumferenceMax);
    CircumferenceMaxLabel ->setText(QString::number(mDetectionParameters.thresholdCircumferenceMax, 'f', 1));

    CircumferenceMinSlider->setDoubleValue(mDetectionParameters.thresholdCircumferenceMin);
    CircumferenceMinLabel ->setText(QString::number(mDetectionParameters.thresholdCircumferenceMin, 'f', 1));

    AspectRatioMinSlider->setDoubleValue(mDetectionParameters.thresholdAspectRatioMin);
    AspectRatioMinLabel ->setText(QString::number(mDetectionParameters.thresholdAspectRatioMin, 'f', 2));

    // Canny edge detection

    CannyThresholdHighSlider->setDoubleValue(mDetectionParameters.cannyThresholdHigh);
    CannyThresholdHighLabel ->setText(QString::number(mDetectionParameters.cannyThresholdHigh, 'f', 1));

    CannyThresholdLowSlider->setDoubleValue(mDetectionParameters.cannyThresholdLow);
    CannyThresholdLowLabel ->setText(QString::number(mDetectionParameters.cannyThresholdLow, 'f', 1));

    CannyBlurLevelSlider->setValue(mDetectionParameters.cannyBlurLevel);
    CannyBlurLevelLabel ->setText(QString::number(mDetectionParameters.cannyBlurLevel));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mDetectionParameters.cannyKernelSize));
    CannyKernelSizeLabel ->setText(QString::number(mDetectionParameters.cannyKernelSize));

    // Learning rates

    AlphaAveragesSlider->setDoubleValue(mDetectionParameters.alphaAverages);
    AlphaAveragesLabel ->setText(QString::number(mDetectionParameters.alphaAverages, 'f', 3));

    AlphaPositionSlider->setDoubleValue(mDetectionParameters.alphaPosition);
    AlphaPositionLabel ->setText(QString::number(mDetectionParameters.alphaPosition, 'f', 2));

    AlphaFeaturesSlider->setDoubleValue(mDetectionParameters.alphaFeatures);
    AlphaFeaturesLabel ->setText(QString::number(mDetectionParameters.alphaFeatures, 'f', 2));

    AlphaCertaintySlider->setDoubleValue(mDetectionParameters.alphaCertainty);
    AlphaCertaintyLabel->setText(QString::number(mDetectionParameters.alphaCertainty, 'f', 2));

    // Thresholds

    ThresholdCircumferenceSlider->setDoubleValue(mDetectionParameters.thresholdChangeCircumference);
    ThresholdCircumferenceLabel ->setText(QString::number(mDetectionParameters.thresholdChangeCircumference, 'f', 2));

    ThresholdAspectRatioSlider->setDoubleValue(mDetectionParameters.thresholdChangeAspectRatio);
    ThresholdAspectRatioLabel ->setText(QString::number(mDetectionParameters.thresholdChangeAspectRatio, 'f', 2));

    ThresholdDisplacementSlider->setDoubleValue(mDetectionParameters.thresholdChangePosition);
    ThresholdDisplacementLabel ->setText(QString::number(mDetectionParameters.thresholdChangePosition, 'f', 1));

    ThresholdScoreEdgeSlider->setDoubleValue(mDetectionParameters.thresholdScoreEdge);
    ThresholdScoreEdgeLabel ->setText(QString::number(mDetectionParameters.thresholdScoreEdge, 'f', 2));

    ThresholdScoreFitSlider->setDoubleValue(mDetectionParameters.thresholdScoreEdge);
    ThresholdScoreFitLabel ->setText(QString::number(mDetectionParameters.thresholdScoreEdge, 'f', 2));

    ThresholdScoreDiffEdgeSlider->setDoubleValue(mDetectionParameters.thresholdScoreDiffEdge);
    ThresholdScoreDiffEdgeLabel ->setText(QString::number(mDetectionParameters.thresholdScoreDiffEdge, 'f', 2));

    ThresholdScoreDiffFitSlider->setDoubleValue(mDetectionParameters.thresholdScoreDiffFit);
    ThresholdScoreDiffFitLabel ->setText(QString::number(mDetectionParameters.thresholdScoreDiffFit, 'f', 2));

    ThresholdFitErrorSlider->setDoubleValue(mDetectionParameters.thresholdFitError);
    ThresholdFitErrorLabel ->setText(QString::number(mDetectionParameters.thresholdFitError, 'f', 2));

    // Misc

    GlintSizeSlider->setValue(round(0.5 * mDetectionParameters.glintWdth));
    GlintSizeLabel->setText(QString::number(mDetectionParameters.glintWdth));

    CurvatureOffsetSlider->setDoubleValue(mDetectionParameters.curvatureOffset);
    CurvatureOffsetLabel ->setText(QString::number(mDetectionParameters.curvatureOffset, 'f', 1));

    WindowLengthEdgeSlider->setValue(mDetectionParameters.windowLengthEdge);
    WindowLengthEdgeLabel ->setText(QString::number(mDetectionParameters.windowLengthEdge));

    FitEdgeFractionSlider->setDoubleValue(mDetectionParameters.fitEdgeFraction);
    FitEdgeFractionLabel ->setText(QString::number(mDetectionParameters.fitEdgeFraction, 'f', 2));

    FitEdgeMaximumSlider->setValue(mDetectionParameters.fitEdgeMaximum);
    FitEdgeMaximumLabel ->setText(QString::number(mDetectionParameters.fitEdgeMaximum));

    FitMaximumSlider->setValue(mDetectionParameters.fitMaximum);
    FitMaximumLabel ->setText(QString::number(mDetectionParameters.fitMaximum));
}

void ParameterWidget::setCircumferenceMin(double value)
{
    if (mDetectionParameters.thresholdCircumferenceMax < value)
    {
        CircumferenceMaxSlider->setDoubleValue(value);
    }

    mDetectionParameters.thresholdCircumferenceMin = value;
    CircumferenceMinLabel->setText(QString::number(mDetectionParameters.thresholdCircumferenceMin, 'f', 1));
}

void ParameterWidget::setCircumferenceMax(double value)
{
    if (mDetectionParameters.thresholdCircumferenceMin > value)
    {
        CircumferenceMinSlider->setDoubleValue(value);
    }

    mDetectionParameters.thresholdCircumferenceMax = value;
    CircumferenceMaxLabel->setText(QString::number(mDetectionParameters.thresholdCircumferenceMax, 'f', 1));
}

void ParameterWidget::setAspectRatioMin(double value)
{
    mDetectionParameters.thresholdAspectRatioMin = value;
    AspectRatioMinLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setFitEdgeFraction(double value)
{
    mDetectionParameters.fitEdgeFraction = value;
    FitEdgeFractionLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setCannyThresholdLow(double value)
{
    if (value > mDetectionParameters.cannyThresholdHigh)
    {
        CannyThresholdHighSlider->setDoubleValue(value);
    }

    mDetectionParameters.cannyThresholdLow = value;
    CannyThresholdLowLabel->setText(QString::number(value, 'f', 1));

}

void ParameterWidget::setCannyThresholdHigh(double value)
{
    if (value < mDetectionParameters.cannyThresholdLow)
    {
        CannyThresholdLowSlider->setDoubleValue(value);
    }

    mDetectionParameters.cannyThresholdHigh = value;
    CannyThresholdHighLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setCannyKernelSize(int value)
{
    int newValue = 2 * value - 1;
    mDetectionParameters.cannyKernelSize = newValue;
    CannyKernelSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCannyBlurLevel(int value)
{
    mDetectionParameters.cannyBlurLevel = value;
    CannyBlurLevelLabel->setText(QString::number(value));
}

void ParameterWidget::setAlphaAverages(double value)
{
    mDetectionParameters.alphaAverages = value;
    AlphaAveragesLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setAlphaPosition(double value)
{
    mDetectionParameters.alphaPosition = value;
    AlphaPositionLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setAlphaFeatures(double value)
{
    mDetectionParameters.alphaFeatures = value;
    AlphaFeaturesLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setAlphaCertainty(double value)
{
    mDetectionParameters.alphaCertainty = value;
    AlphaCertaintyLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdCircumference(double value)
{
    mDetectionParameters.thresholdChangeCircumference = value;
    ThresholdCircumferenceLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdAspectRatio(double value)
{
    mDetectionParameters.thresholdChangeAspectRatio = value;
    ThresholdAspectRatioLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdDisplacement(double value)
{
    mDetectionParameters.thresholdChangePosition = value;
    ThresholdDisplacementLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setThresholdScoreEdge(double value)
{
    mDetectionParameters.thresholdScoreEdge = value;
    ThresholdScoreEdgeLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreFit(double value)
{
    mDetectionParameters.thresholdScoreFit = value;
    ThresholdScoreFitLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreDiffEdge(double value)
{
    mDetectionParameters.thresholdScoreDiffEdge = value;
    ThresholdScoreDiffEdgeLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreDiffFit(double value)
{
    mDetectionParameters.thresholdScoreDiffFit = value;
    ThresholdScoreDiffFitLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setGlintSize(int value)
{
    int newValue = 2 * value;
    mDetectionParameters.glintWdth = newValue;
    GlintSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCurvatureOffset(double value)
{
    mDetectionParameters.curvatureOffset = value;
    CurvatureOffsetLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setWindowLengthEdge(int value)
{
    mDetectionParameters.windowLengthEdge = value;
    WindowLengthEdgeLabel->setText(QString::number(value));
}

void ParameterWidget::setFitEdgeMaximum(int value)
{
    mDetectionParameters.fitEdgeMaximum = value;
    FitEdgeMaximumLabel->setText(QString::number(value));
}

void ParameterWidget::setFitMaximum(int value)
{
    mDetectionParameters.fitMaximum = value;
    FitMaximumLabel->setText(QString::number(value));
}

void ParameterWidget::setThresholdFitError(double value)
{
    mDetectionParameters.thresholdFitError = value;
    ThresholdFitErrorLabel->setText(QString::number(value, 'f', 2));
}


