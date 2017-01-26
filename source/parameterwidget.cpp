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

    // Pupil fraction

    QLabel *AspectRatioMinTextBox = new QLabel;
    AspectRatioMinTextBox->setText("<b>Fraction minimum:</b>");

    AspectRatioMinLabel  = new QLabel;
    AspectRatioMinSlider = new SliderDouble;
    AspectRatioMinSlider->setPrecision(2);
    AspectRatioMinSlider->setDoubleRange(0.0, 1.0);
    AspectRatioMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(AspectRatioMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAspectRatioMin(double)));

    // Edge intensity

    QLabel *EdgeIntensityOffsetTextBox = new QLabel;
    EdgeIntensityOffsetTextBox->setText("<b>Edge intensity offset:</b>");

    EdgeIntensityOffsetLabel  = new QLabel;
    EdgeIntensityOffsetSlider = new SliderDouble();
    EdgeIntensityOffsetSlider->setPrecision(1);
    EdgeIntensityOffsetSlider->setDoubleRange(0.0, 255.0);
    EdgeIntensityOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EdgeIntensityOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEdgeIntensityOffset(double)));

    // Sliders for canny edge parameters

    QLabel *CannyEdgeTextBox = new QLabel;
    CannyEdgeTextBox->setText("<b>Canny edge detection</b>");
    CannyEdgeTextBox->setAlignment(Qt::AlignCenter);

    QLabel *CannyThresholdLowTextBox = new QLabel;
    CannyThresholdLowTextBox->setText("<b>Low threshold:</b>");

    CannyThresholdLowLabel  = new QLabel;
    CannyThresholdLowSlider = new QSlider;
    CannyThresholdLowSlider->setRange(0, 300);
    CannyThresholdLowSlider->setOrientation(Qt::Horizontal);
    CannyThresholdLowSlider->setSingleStep(1);
    QObject::connect(CannyThresholdLowSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyThresholdLow(int)));

    QLabel *CannyThresholdHighTextBox = new QLabel;
    CannyThresholdHighTextBox->setText("<b>High threshold:</b>");

    CannyThresholdHighLabel  = new QLabel;
    CannyThresholdHighSlider = new QSlider;
    CannyThresholdHighSlider->setRange(0, 300);
    CannyThresholdHighSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyThresholdHighSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyThresholdHigh(int)));

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

    QLabel *AlphaAverageTextBox       = new QLabel;
    QLabel *AlphaPredictionTextBox    = new QLabel;
    QLabel *AlphaMiscellaneousTextBox = new QLabel;
    QLabel *AlphaMomentumTextBox      = new QLabel;

    AlphaAverageTextBox      ->setText("<b>Average:</b>");
    AlphaPredictionTextBox   ->setText("<b>Prediction:</b>");
    AlphaMiscellaneousTextBox->setText("<b>Miscellaneous:</b>");
    AlphaMomentumTextBox     ->setText("<b>Momentum:</b>");

    AlphaAverageLabel       = new QLabel;
    AlphaPredictionLabel    = new QLabel;
    AlphaMiscellaneousLabel = new QLabel;
    AlphaMomentumLabel      = new QLabel;

    AlphaAverageSlider       = new SliderDouble;
    AlphaPredictionSlider    = new SliderDouble;
    AlphaMiscellaneousSlider = new SliderDouble;
    AlphaMomentumSlider      = new SliderDouble;

    AlphaAverageSlider      ->setPrecision(3);
    AlphaPredictionSlider   ->setPrecision(2);
    AlphaMiscellaneousSlider->setPrecision(2);
    AlphaMomentumSlider     ->setPrecision(2);

    AlphaAverageSlider      ->setDoubleRange(0, 0.1);
    AlphaPredictionSlider   ->setDoubleRange(0, 1.0);
    AlphaMiscellaneousSlider->setDoubleRange(0, 1.0);
    AlphaMomentumSlider     ->setDoubleRange(0, 1.0);

    AlphaAverageSlider      ->setOrientation(Qt::Horizontal);
    AlphaPredictionSlider   ->setOrientation(Qt::Horizontal);
    AlphaMiscellaneousSlider->setOrientation(Qt::Horizontal);
    AlphaMomentumSlider     ->setOrientation(Qt::Horizontal);

    QObject::connect(AlphaAverageSlider,       SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaAverage(double)));
    QObject::connect(AlphaPredictionSlider,    SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaPrediction(double)));
    QObject::connect(AlphaMiscellaneousSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaMiscellaneous(double)));
    QObject::connect(AlphaMomentumSlider,      SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaMomentum(double)));

    // Threshold parameters

    QLabel *ThresholdParametersTextBox = new QLabel;
    ThresholdParametersTextBox->setText("<b>Change thresholds</b>");
    ThresholdParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *ThresholdCircumferenceTextBox = new QLabel;
    ThresholdCircumferenceTextBox->setText("<b>Circumference threshold:</b>");

    QLabel *ThresholdAspectRatioTextBox = new QLabel;
    ThresholdAspectRatioTextBox->setText("<b>Aspect ratio threshold:</b>");

    ThresholdCircumferenceLabel  = new QLabel;
    ThresholdCircumferenceSlider = new SliderDouble;
    ThresholdCircumferenceSlider->setPrecision(1);
    ThresholdCircumferenceSlider->setDoubleRange(0, 50);
    ThresholdCircumferenceSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdCircumferenceSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdCircumference(double)));

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

    QLabel *HaarOffsetTextBox = new QLabel;
    HaarOffsetTextBox->setText("<b>Pupil Haar-offset:</b>");

    HaarOffsetLabel  = new QLabel;
    HaarOffsetSlider = new QSlider;
    HaarOffsetSlider->setRange(0, 50);
    HaarOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(HaarOffsetSlider, SIGNAL(valueChanged(int)), this, SLOT(setHaarOffset(int)));

    QLabel *GlintSizeTextBox = new QLabel;
    GlintSizeTextBox->setText("<b>Glint size:</b>");

    GlintSizeLabel  = new QLabel;
    GlintSizeSlider = new QSlider;
    GlintSizeSlider->setRange(0, 10);
    GlintSizeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GlintSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setGlintSize(int)));

    QLabel *CurvatureFactorTextBox = new QLabel;
    CurvatureFactorTextBox->setText("<b>Curvature factor:</b>");

    CurvatureFactorLabel  = new QLabel;
    CurvatureFactorSlider = new SliderDouble;
    CurvatureFactorSlider->setPrecision(2);
    CurvatureFactorSlider->setDoubleRange(1.0, 1.2);
    CurvatureFactorSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureFactorSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureFactor(double)));

    QLabel *CurvatureOffsetTextBox = new QLabel;
    CurvatureOffsetTextBox->setText("<b>Curvature offset:</b>");

    CurvatureOffsetLabel  = new QLabel;
    CurvatureOffsetSlider = new SliderDouble;
    CurvatureOffsetSlider->setPrecision(1);
    CurvatureOffsetSlider->setDoubleRange(0, 180);
    CurvatureOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureOffset(double)));

    QLabel *EdgeLengthMinimumTextBox = new QLabel;
    EdgeLengthMinimumTextBox->setText("<b>Edge minimum length:");

    EdgeLengthMinimumLabel  = new QLabel;
    EdgeLengthMinimumSlider = new SliderDouble;
    EdgeLengthMinimumSlider->setPrecision(2);
    EdgeLengthMinimumSlider->setDoubleRange(0, 1.0);
    EdgeLengthMinimumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EdgeLengthMinimumSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEdgeLengthMinimum(double)));

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

    MainLayout->addWidget(CircumferenceMaxTextBox,     1, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CircumferenceMinTextBox,     2, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AspectRatioMinTextBox,       3, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EdgeIntensityOffsetTextBox,       4, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(CannyThresholdHighTextBox,        6, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyThresholdLowTextBox,         7, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyKernelSizeTextBox,           8, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyBlurLevelTextBox,            9, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(AlphaPredictionTextBox,          11, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaAverageTextBox,             12, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaMomentumTextBox,            13, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AlphaMiscellaneousTextBox,       14, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(ThresholdCircumferenceTextBox,   16, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdAspectRatioTextBox,     17, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(HaarOffsetTextBox,          19, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(GlintSizeTextBox,                20, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CurvatureFactorTextBox,          21, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CurvatureOffsetTextBox,          22, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EdgeLengthMinimumTextBox,        23, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EllipseFitNumberMaximumTextBox,  24, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(EllipseFitErrorMaximumTextBox,   25, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(TitleLimitTextBox,                0, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(CircumferenceMaxSlider,      1, 1);
    MainLayout->addWidget(CircumferenceMinSlider,      2, 1);
    MainLayout->addWidget(AspectRatioMinSlider,        3, 1);
    MainLayout->addWidget(EdgeIntensityOffsetSlider,        4, 1);
    MainLayout->addWidget(TitleCannyTextBox,                5, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(CannyThresholdHighSlider,         6, 1);
    MainLayout->addWidget(CannyThresholdLowSlider,          7, 1);
    MainLayout->addWidget(CannyKernelSizeSlider,            8, 1);
    MainLayout->addWidget(CannyBlurLevelSlider,             9, 1);
    MainLayout->addWidget(TitleLearnTextBox,               10, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(AlphaPredictionSlider,           11, 1);
    MainLayout->addWidget(AlphaAverageSlider,              12, 1);
    MainLayout->addWidget(AlphaMomentumSlider,             13, 1);
    MainLayout->addWidget(AlphaMiscellaneousSlider,        14, 1);
    MainLayout->addWidget(TitleChangeTextBox,              15, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(ThresholdCircumferenceSlider,    16, 1);
    MainLayout->addWidget(ThresholdAspectRatioSlider,      17, 1);
    MainLayout->addWidget(TitleMiscTextBox,                18, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(HaarOffsetSlider,           19, 1);
    MainLayout->addWidget(GlintSizeSlider,                 20, 1);
    MainLayout->addWidget(CurvatureFactorSlider,           21, 1);
    MainLayout->addWidget(CurvatureOffsetSlider,           22, 1);
    MainLayout->addWidget(EdgeLengthMinimumSlider,         23, 1);
    MainLayout->addWidget(EllipseFitNumberMaximumSlider,   24, 1);
    MainLayout->addWidget(EllipseFitErrorMaximumSlider,    25, 1);

    MainLayout->addWidget(CircumferenceMaxLabel,    1, 2);
    MainLayout->addWidget(CircumferenceMinLabel,    2, 2);
    MainLayout->addWidget(AspectRatioMinLabel,      3, 2);
    MainLayout->addWidget(EdgeIntensityOffsetLabel,      4, 2);

    MainLayout->addWidget(CannyThresholdHighLabel,       6, 2);
    MainLayout->addWidget(CannyThresholdLowLabel,        7, 2);
    MainLayout->addWidget(CannyKernelSizeLabel,          8, 2);
    MainLayout->addWidget(CannyBlurLevelLabel,           9, 2);

    MainLayout->addWidget(AlphaPredictionLabel,          11, 2);
    MainLayout->addWidget(AlphaAverageLabel,             12, 2);
    MainLayout->addWidget(AlphaMomentumLabel,            13, 2);
    MainLayout->addWidget(AlphaMiscellaneousLabel,       14, 2);

    MainLayout->addWidget(ThresholdCircumferenceLabel,   16, 2);
    MainLayout->addWidget(ThresholdAspectRatioLabel,     17, 2);

    MainLayout->addWidget(HaarOffsetLabel,          19, 2);
    MainLayout->addWidget(GlintSizeLabel,                20, 2);
    MainLayout->addWidget(CurvatureFactorLabel,          21, 2);
    MainLayout->addWidget(CurvatureOffsetLabel,          22, 2);
    MainLayout->addWidget(EdgeLengthMinimumLabel,        23, 2);
    MainLayout->addWidget(EllipseFitNumberMaximumLabel,  24, 2);
    MainLayout->addWidget(EllipseFitErrorMaximumLabel,   25, 2);

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
    CircumferenceMinSlider->setDoubleValue(mDetectionParameters.circumferenceMin);
    CircumferenceMinLabel ->setText(QString::number(mDetectionParameters.circumferenceMin, 'f', 1));

    CircumferenceMaxSlider->setDoubleValue(mDetectionParameters.circumferenceMax);
    CircumferenceMaxLabel ->setText(QString::number(mDetectionParameters.circumferenceMax, 'f', 1));

    AspectRatioMinSlider->setDoubleValue(mDetectionParameters.aspectRatioMin);
    AspectRatioMinLabel ->setText(QString::number(mDetectionParameters.aspectRatioMin, 'f', 2));

    EdgeIntensityOffsetSlider->setDoubleValue(mDetectionParameters.edgeIntensityOffset);
    EdgeIntensityOffsetLabel ->setText(QString::number(mDetectionParameters.edgeIntensityOffset, 'f', 1));

    CannyThresholdLowSlider->setValue(mDetectionParameters.cannyThresholdLow);
    CannyThresholdLowLabel ->setText(QString::number(mDetectionParameters.cannyThresholdLow));

    CannyThresholdHighSlider->setValue(mDetectionParameters.cannyThresholdHigh);
    CannyThresholdHighLabel ->setText(QString::number(mDetectionParameters.cannyThresholdHigh));

    CannyBlurLevelSlider->setValue(mDetectionParameters.cannyBlurLevel);
    CannyBlurLevelLabel ->setText(QString::number(mDetectionParameters.cannyBlurLevel));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mDetectionParameters.cannyKernelSize));
    CannyKernelSizeLabel ->setText(QString::number(mDetectionParameters.cannyKernelSize));

    AlphaAverageSlider      ->setDoubleValue(mDetectionParameters.alphaAverage);
    AlphaPredictionSlider   ->setDoubleValue(mDetectionParameters.alphaPrediction);
    AlphaMiscellaneousSlider->setDoubleValue(mDetectionParameters.alphaMiscellaneous);
    AlphaMomentumSlider     ->setDoubleValue(mDetectionParameters.alphaMomentum);

    AlphaAverageLabel      ->setText(QString::number(mDetectionParameters.alphaAverage, 'f', 3));
    AlphaPredictionLabel   ->setText(QString::number(mDetectionParameters.alphaPrediction, 'f', 2));
    AlphaMiscellaneousLabel->setText(QString::number(mDetectionParameters.alphaMiscellaneous, 'f', 2));
    AlphaMomentumLabel     ->setText(QString::number(mDetectionParameters.alphaMomentum, 'f', 2));

    ThresholdCircumferenceSlider->setDoubleValue(mDetectionParameters.circumferenceChangeThreshold);
    ThresholdCircumferenceLabel ->setText(QString::number(mDetectionParameters.circumferenceChangeThreshold, 'f', 1));

    ThresholdAspectRatioSlider->setDoubleValue(mDetectionParameters.aspectRatioChangeThreshold);
    ThresholdAspectRatioLabel ->setText(QString::number(mDetectionParameters.aspectRatioChangeThreshold, 'f', 2));

    HaarOffsetSlider->setValue(mDetectionParameters.pupilOffset);
    HaarOffsetLabel ->setText(QString::number(mDetectionParameters.pupilOffset));

    GlintSizeSlider->setValue(round(0.5 * mDetectionParameters.glintSize));
    GlintSizeLabel->setText(QString::number(mDetectionParameters.glintSize));

    CurvatureFactorSlider->setDoubleValue(mDetectionParameters.curvatureFactor);
    CurvatureFactorLabel ->setText(QString::number(mDetectionParameters.curvatureFactor, 'f', 2));

    CurvatureOffsetSlider->setDoubleValue(mDetectionParameters.curvatureOffsetMin);
    CurvatureOffsetLabel ->setText(QString::number(mDetectionParameters.curvatureOffsetMin, 'f', 1));

    EdgeLengthMinimumSlider->setDoubleValue(mDetectionParameters.edgeLengthMinimum);
    EdgeLengthMinimumLabel ->setText(QString::number(mDetectionParameters.edgeLengthMinimum, 'f', 2));

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

void ParameterWidget::setEdgeIntensityOffset(double value)
{
    mDetectionParameters.edgeIntensityOffset = value;
    EdgeIntensityOffsetLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setEdgeLengthMinimum(double value)
{
    mDetectionParameters.edgeLengthMinimum = value;
    EdgeLengthMinimumLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setCannyThresholdLow(int value)
{
    if (value > mDetectionParameters.cannyThresholdHigh)
    {
        CannyThresholdHighSlider->setValue(value);
    }

    mDetectionParameters.cannyThresholdLow = value;
    CannyThresholdLowLabel->setText(QString::number(value));

}

void ParameterWidget::setCannyThresholdHigh(int value)
{
    if (value < mDetectionParameters.cannyThresholdLow)
    {
        CannyThresholdLowSlider->setValue(value);
    }

    mDetectionParameters.cannyThresholdHigh = value;
    CannyThresholdHighLabel->setText(QString::number(value));
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

void ParameterWidget::setAlphaPrediction(double value)
{
    mDetectionParameters.alphaPrediction = value;
    AlphaPredictionLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setAlphaMiscellaneous(double value)
{
    mDetectionParameters.alphaMiscellaneous = value;
    AlphaMiscellaneousLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setAlphaMomentum(double value)
{
    mDetectionParameters.alphaMomentum = value;
    AlphaMomentumLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdCircumference(double value)
{
    mDetectionParameters.circumferenceChangeThreshold = value;
    ThresholdCircumferenceLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setThresholdAspectRatio(double value)
{
    mDetectionParameters.aspectRatioChangeThreshold = value;
    ThresholdAspectRatioLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setHaarOffset(int value)
{
    mDetectionParameters.pupilOffset = value;
    HaarOffsetLabel->setText(QString::number(value));
}

void ParameterWidget::setGlintSize(int value)
{
    int newValue = 2 * value;
    mDetectionParameters.glintSize = newValue;
    GlintSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCurvatureFactor(double value)
{
    mDetectionParameters.curvatureFactor = value;
    CurvatureFactorLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setCurvatureOffset(double value)
{
    mDetectionParameters.curvatureOffsetMin = value;
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


