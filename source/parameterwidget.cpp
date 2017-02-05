#include "headers/parameterwidget.h"

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
    CircumferenceMinSlider->setDoubleRange(0, 500);
    CircumferenceMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CircumferenceMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCircumferenceMin(double)));

    QLabel *CircumferenceMaxTextBox = new QLabel;
    CircumferenceMaxTextBox->setText("<b>Circumference max:</b>");

    CircumferenceMaxLabel  = new QLabel;
    CircumferenceMaxSlider = new SliderDouble;
    CircumferenceMaxSlider->setPrecision(1);
    CircumferenceMaxSlider->setDoubleRange(0, 500);
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

    QLabel *AlphaAverageTextBox  = new QLabel;
    QLabel *AlphaPositionTextBox = new QLabel;
    QLabel *AlphaFeaturesTextBox = new QLabel;
    QLabel *AlphaMomentumTextBox = new QLabel;

    AlphaAverageTextBox ->setText("<b>Averages:</b>");
    AlphaPositionTextBox->setText("<b>Position:</b>");
    AlphaFeaturesTextBox->setText("<b>Features:</b>");
    AlphaMomentumTextBox->setText("<b>Momentum:</b>");

    AlphaAverageLabel  = new QLabel;
    AlphaPositionLabel = new QLabel;
    AlphaFeaturesLabel = new QLabel;
    AlphaMomentumLabel = new QLabel;

    AlphaAverageSlider  = new SliderDouble;
    AlphaPositionSlider = new SliderDouble;
    AlphaFeaturesSlider = new SliderDouble;
    AlphaMomentumSlider = new SliderDouble;

    AlphaAverageSlider ->setPrecision(3);
    AlphaPositionSlider->setPrecision(2);
    AlphaFeaturesSlider->setPrecision(2);
    AlphaMomentumSlider->setPrecision(2);

    AlphaAverageSlider ->setDoubleRange(0, 0.1);
    AlphaPositionSlider->setDoubleRange(0, 1.0);
    AlphaFeaturesSlider->setDoubleRange(0, 1.0);
    AlphaMomentumSlider->setDoubleRange(0, 1.0);

    AlphaAverageSlider ->setOrientation(Qt::Horizontal);
    AlphaPositionSlider->setOrientation(Qt::Horizontal);
    AlphaFeaturesSlider->setOrientation(Qt::Horizontal);
    AlphaMomentumSlider->setOrientation(Qt::Horizontal);

    QObject::connect(AlphaAverageSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaAverage(double)));
    QObject::connect(AlphaPositionSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaPosition(double)));
    QObject::connect(AlphaFeaturesSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaFeatures(double)));
    QObject::connect(AlphaMomentumSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaMomentum(double)));

    // Threshold parameters

    QLabel *ThresholdParametersTextBox = new QLabel;
    ThresholdParametersTextBox->setText("<b>Change thresholds</b>");
    ThresholdParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *ThresholdCircumferenceTextBox = new QLabel;
    ThresholdCircumferenceTextBox->setText("<b>Circumference:</b>");

    QLabel *ThresholdAspectRatioTextBox = new QLabel;
    ThresholdAspectRatioTextBox->setText("<b>Aspect ratio:</b>");

    QLabel *ThresholdDisplacementTextBox = new QLabel;
    ThresholdDisplacementTextBox->setText("<b>Displacement:</b>");

    ThresholdCircumferenceLabel  = new QLabel;
    ThresholdCircumferenceSlider = new SliderDouble;
    ThresholdCircumferenceSlider->setPrecision(1);
    ThresholdCircumferenceSlider->setDoubleRange(0, 25);
    ThresholdCircumferenceSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdCircumferenceSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdCircumference(double)));

    ThresholdDisplacementLabel  = new QLabel;
    ThresholdDisplacementSlider = new SliderDouble;
    ThresholdDisplacementSlider->setPrecision(1);
    ThresholdDisplacementSlider->setDoubleRange(0, 25);
    ThresholdDisplacementSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdDisplacementSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdDisplacement(double)));

    ThresholdAspectRatioLabel  = new QLabel;
    ThresholdAspectRatioSlider = new SliderDouble;
    ThresholdAspectRatioSlider->setPrecision(2);
    ThresholdAspectRatioSlider->setDoubleRange(0, 1.0);
    ThresholdAspectRatioSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdAspectRatioSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdAspectRatio(double)));

    // Miscellaneous parameters

    QLabel *MiscParametersTextBox = new QLabel;
    MiscParametersTextBox->setText("<b>Miscellaneous</b>");
    MiscParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *GlintSizeTextBox = new QLabel;
    GlintSizeTextBox->setText("<b>Glint size:</b>");

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
    CurvatureOffsetSlider->setDoubleRange(0, 50);
    CurvatureOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureOffset(double)));

    QLabel *ThresholdScoreTextBox = new QLabel;
    ThresholdScoreTextBox->setText("<b>Score threshold:</b>");

    ThresholdScoreLabel  = new QLabel;
    ThresholdScoreSlider = new SliderDouble;
    ThresholdScoreSlider->setPrecision(2);
    ThresholdScoreSlider->setDoubleRange(0, 0.5);
    ThresholdScoreSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScore(double)));


    QLabel *EdgeLengthFractionTextBox = new QLabel;
    EdgeLengthFractionTextBox->setText("<b>Edge minimum length:");

    EdgeLengthFractionLabel  = new QLabel;
    EdgeLengthFractionSlider = new SliderDouble;
    EdgeLengthFractionSlider->setPrecision(2);
    EdgeLengthFractionSlider->setDoubleRange(0, 1.0);
    EdgeLengthFractionSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EdgeLengthFractionSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEdgeLengthFraction(double)));

    QLabel *EllipseFitNumberMaximumTextBox = new QLabel;
    EllipseFitNumberMaximumTextBox->setText("<b>Edge maximum fit number:</b>");

    EllipseFitNumberMaximumLabel  = new QLabel;
    EllipseFitNumberMaximumSlider = new QSlider;
    EllipseFitNumberMaximumSlider->setRange(0, 7);
    EllipseFitNumberMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EllipseFitNumberMaximumSlider, SIGNAL(valueChanged(int)), this, SLOT(setEllipseFitNumberMaximum(int)));

    QLabel *EllipseFitErrorMaximumTextBox = new QLabel;
    EllipseFitErrorMaximumTextBox->setText("<b>Ellipse maximum fit error:</b>");

    EllipseFitErrorMaximumLabel  = new QLabel;
    EllipseFitErrorMaximumSlider = new SliderDouble;
    EllipseFitErrorMaximumSlider->setPrecision(1);
    EllipseFitErrorMaximumSlider->setDoubleRange(0, 160);
    EllipseFitErrorMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EllipseFitErrorMaximumSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEllipseFitErrorMaximum(double)));

    QLabel *TitleLimitTextBox  = new QLabel;
    QLabel *TitleCannyTextBox  = new QLabel;
    QLabel *TitleLearnTextBox  = new QLabel;
    QLabel *TitleChangeTextBox = new QLabel;
    QLabel *TitleMiscTextBox   = new QLabel;

    TitleLimitTextBox ->setText("<b>Variable limits</b>");
    TitleCannyTextBox ->setText("<b>Canny edge detection</b>");
    TitleLearnTextBox ->setText("<b>Learning rates</b>");
    TitleChangeTextBox->setText("<b>Change thresholds</b>");
    TitleMiscTextBox  ->setText("<b>Miscellaneous</b>");

    QWidget *MainWidget = new QWidget;
    QGridLayout *MainLayout = new QGridLayout(MainWidget);

    MainLayout->addWidget(CircumferenceMaxTextBox,          1, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CircumferenceMinTextBox,          2, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AspectRatioMinTextBox,            3, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(CannyThresholdHighTextBox,        6, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyThresholdLowTextBox,         7, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyKernelSizeTextBox,           8, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyBlurLevelTextBox,            9, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(AlphaPositionTextBox,             11, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaAverageTextBox,              12, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaMomentumTextBox,             13, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaFeaturesTextBox,             14, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(ThresholdDisplacementTextBox,     16, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdCircumferenceTextBox,    17, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdAspectRatioTextBox,      18, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(GlintSizeTextBox,                 21, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CurvatureOffsetTextBox,           22, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreTextBox,            23, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EdgeLengthFractionTextBox,        24, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EllipseFitNumberMaximumTextBox,   25, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EllipseFitErrorMaximumTextBox,    26, 0, 1, 1, Qt::AlignRight);

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
    MainLayout->addWidget(AlphaPositionSlider,             11, 1);
    MainLayout->addWidget(AlphaAverageSlider,              12, 1);
    MainLayout->addWidget(AlphaMomentumSlider,             13, 1);
    MainLayout->addWidget(AlphaFeaturesSlider,             14, 1);
    MainLayout->addWidget(TitleChangeTextBox,              15, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(ThresholdDisplacementSlider,     16, 1);
    MainLayout->addWidget(ThresholdCircumferenceSlider,    17, 1);
    MainLayout->addWidget(ThresholdAspectRatioSlider,      18, 1);
    MainLayout->addWidget(TitleMiscTextBox,                19, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(GlintSizeSlider,                 21, 1);
    MainLayout->addWidget(CurvatureOffsetSlider,           22, 1);
    MainLayout->addWidget(ThresholdScoreSlider,            23, 1);
    MainLayout->addWidget(EdgeLengthFractionSlider,        24, 1);
    MainLayout->addWidget(EllipseFitNumberMaximumSlider,   25, 1);
    MainLayout->addWidget(EllipseFitErrorMaximumSlider,    26, 1);

    MainLayout->addWidget(CircumferenceMaxLabel,         1, 2);
    MainLayout->addWidget(CircumferenceMinLabel,         2, 2);
    MainLayout->addWidget(AspectRatioMinLabel,           3, 2);

    MainLayout->addWidget(CannyThresholdHighLabel,       6, 2);
    MainLayout->addWidget(CannyThresholdLowLabel,        7, 2);
    MainLayout->addWidget(CannyKernelSizeLabel,          8, 2);
    MainLayout->addWidget(CannyBlurLevelLabel,           9, 2);

    MainLayout->addWidget(AlphaPositionLabel,            11, 2);
    MainLayout->addWidget(AlphaAverageLabel,             12, 2);
    MainLayout->addWidget(AlphaMomentumLabel,            13, 2);
    MainLayout->addWidget(AlphaFeaturesLabel,            14, 2);

    MainLayout->addWidget(ThresholdDisplacementLabel,    16, 2);
    MainLayout->addWidget(ThresholdCircumferenceLabel,   17, 2);
    MainLayout->addWidget(ThresholdAspectRatioLabel,     18, 2);

    MainLayout->addWidget(GlintSizeLabel,                21, 2);
    MainLayout->addWidget(CurvatureOffsetLabel,          22, 2);
    MainLayout->addWidget(ThresholdScoreLabel,           23, 2);
    MainLayout->addWidget(EdgeLengthFractionLabel,       24, 2);
    MainLayout->addWidget(EllipseFitNumberMaximumLabel,  25, 2);
    MainLayout->addWidget(EllipseFitErrorMaximumLabel,   26, 2);

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
    CircumferenceMaxSlider->setDoubleValue(mDetectionParameters.circumferenceMax);
    CircumferenceMaxLabel ->setText(QString::number(mDetectionParameters.circumferenceMax, 'f', 1));

    CircumferenceMinSlider->setDoubleValue(mDetectionParameters.circumferenceMin);
    CircumferenceMinLabel ->setText(QString::number(mDetectionParameters.circumferenceMin, 'f', 1));

    AspectRatioMinSlider->setDoubleValue(mDetectionParameters.aspectRatioMin);
    AspectRatioMinLabel ->setText(QString::number(mDetectionParameters.aspectRatioMin, 'f', 2));

    CannyThresholdHighSlider->setDoubleValue(mDetectionParameters.cannyThresholdHigh);
    CannyThresholdHighLabel ->setText(QString::number(mDetectionParameters.cannyThresholdHigh, 'f', 1));

    CannyThresholdLowSlider->setDoubleValue(mDetectionParameters.cannyThresholdLow);
    CannyThresholdLowLabel ->setText(QString::number(mDetectionParameters.cannyThresholdLow, 'f', 1));

    CannyBlurLevelSlider->setValue(mDetectionParameters.cannyBlurLevel);
    CannyBlurLevelLabel ->setText(QString::number(mDetectionParameters.cannyBlurLevel));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mDetectionParameters.cannyKernelSize));
    CannyKernelSizeLabel ->setText(QString::number(mDetectionParameters.cannyKernelSize));

    AlphaAverageSlider ->setDoubleValue(mDetectionParameters.alphaAverage);
    AlphaPositionSlider->setDoubleValue(mDetectionParameters.alphaPosition);
    AlphaFeaturesSlider->setDoubleValue(mDetectionParameters.alphaFeatures);
    AlphaMomentumSlider->setDoubleValue(mDetectionParameters.alphaMomentum);

    AlphaAverageLabel ->setText(QString::number(mDetectionParameters.alphaAverage, 'f', 3));
    AlphaPositionLabel->setText(QString::number(mDetectionParameters.alphaPosition, 'f', 2));
    AlphaFeaturesLabel->setText(QString::number(mDetectionParameters.alphaFeatures, 'f', 2));
    AlphaMomentumLabel->setText(QString::number(mDetectionParameters.alphaMomentum, 'f', 2));

    ThresholdCircumferenceSlider->setDoubleValue(mDetectionParameters.changeThresholdCircumference);
    ThresholdCircumferenceLabel ->setText(QString::number(mDetectionParameters.changeThresholdCircumference, 'f', 1));

    ThresholdAspectRatioSlider->setDoubleValue(mDetectionParameters.changeThresholdAspectRatio);
    ThresholdAspectRatioLabel ->setText(QString::number(mDetectionParameters.changeThresholdAspectRatio, 'f', 2));

    ThresholdDisplacementSlider->setDoubleValue(mDetectionParameters.changeThresholdPosition);
    ThresholdDisplacementLabel ->setText(QString::number(mDetectionParameters.changeThresholdPosition, 'f', 1));

    ThresholdScoreSlider->setDoubleValue(mDetectionParameters.scoreThreshold);
    ThresholdScoreLabel ->setText(QString::number(mDetectionParameters.scoreThreshold, 'f', 2));

    GlintSizeSlider->setValue(round(0.5 * mDetectionParameters.glintWdth));
    GlintSizeLabel->setText(QString::number(mDetectionParameters.glintWdth));

    CurvatureOffsetSlider->setDoubleValue(mDetectionParameters.curvatureOffset);
    CurvatureOffsetLabel ->setText(QString::number(mDetectionParameters.curvatureOffset, 'f', 1));

    EdgeLengthFractionSlider->setDoubleValue(mDetectionParameters.edgeLengthFraction);
    EdgeLengthFractionLabel ->setText(QString::number(mDetectionParameters.edgeLengthFraction, 'f', 2));

    EllipseFitNumberMaximumSlider->setValue(mDetectionParameters.ellipseFitNumberMaximum);
    EllipseFitNumberMaximumLabel ->setText(QString::number(mDetectionParameters.ellipseFitNumberMaximum));

    EllipseFitErrorMaximumSlider->setDoubleValue(mDetectionParameters.ellipseFitErrorMaximum);
    EllipseFitErrorMaximumLabel ->setText(QString::number(mDetectionParameters.ellipseFitErrorMaximum, 'f', 1));
}

void ParameterWidget::setCircumferenceMin(double value)
{
    if (mDetectionParameters.circumferenceMax < value)
    {
        CircumferenceMaxSlider->setDoubleValue(value);
    }

    mDetectionParameters.circumferenceMin = value;
    CircumferenceMinLabel->setText(QString::number(mDetectionParameters.circumferenceMin, 'f', 1));
}

void ParameterWidget::setCircumferenceMax(double value)
{
    if (mDetectionParameters.circumferenceMin > value)
    {
        CircumferenceMinSlider->setDoubleValue(value);
    }

    mDetectionParameters.circumferenceMax = value;
    CircumferenceMaxLabel->setText(QString::number(mDetectionParameters.circumferenceMax, 'f', 1));
}

void ParameterWidget::setAspectRatioMin(double value)
{
    mDetectionParameters.aspectRatioMin = value;
    AspectRatioMinLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setEdgeLengthFraction(double value)
{
    mDetectionParameters.edgeLengthFraction = value;
    EdgeLengthFractionLabel->setText(QString::number(value, 'f', 2));
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

void ParameterWidget::setAlphaAverage(double value)
{
    mDetectionParameters.alphaAverage = value;
    AlphaAverageLabel->setText(QString::number(value, 'f', 3));
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

void ParameterWidget::setAlphaMomentum(double value)
{
    mDetectionParameters.alphaMomentum = value;
    AlphaMomentumLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdCircumference(double value)
{
    mDetectionParameters.changeThresholdCircumference = value;
    ThresholdCircumferenceLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setThresholdAspectRatio(double value)
{
    mDetectionParameters.changeThresholdAspectRatio = value;
    ThresholdAspectRatioLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdDisplacement(double value)
{
    mDetectionParameters.changeThresholdPosition = value;
    ThresholdDisplacementLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setThresholdScore(double value)
{
    mDetectionParameters.scoreThreshold = value;
    ThresholdScoreLabel->setText(QString::number(value, 'f', 2));
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

void ParameterWidget::setEllipseFitNumberMaximum(int value)
{
    mDetectionParameters.ellipseFitNumberMaximum = value;
    EllipseFitNumberMaximumLabel->setText(QString::number(value));
}

void ParameterWidget::setEllipseFitErrorMaximum(double value)
{
    mDetectionParameters.ellipseFitErrorMaximum = value;
    EllipseFitErrorMaximumLabel->setText(QString::number(value, 'f', 1));
}


