#include "headers/parameterwidget.h"

ParameterWidget::ParameterWidget(QWidget *parent) : QWidget(parent)
{
    // Averages and Limits/Thresholds

    QLabel *ParametersTextBox = new QLabel;
    ParametersTextBox->setText("<b>Eye-tracking parameters</b>");
    ParametersTextBox->setAlignment(Qt::AlignCenter);

    // Pupil circumference

    QLabel *PupilCircumferenceMinTextBox = new QLabel;
    PupilCircumferenceMinTextBox->setText("<b>Circumference min:</b>");

    PupilCircumferenceMinLabel  = new QLabel;
    PupilCircumferenceMinSlider = new SliderDouble;
    PupilCircumferenceMinSlider->setPrecision(1);
    PupilCircumferenceMinSlider->setDoubleRange(0, 500);
    PupilCircumferenceMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilCircumferenceMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumferenceMin(double)));

    QLabel *PupilCircumferenceMaxTextBox = new QLabel;
    PupilCircumferenceMaxTextBox->setText("<b>Circumference max:</b>");

    PupilCircumferenceMaxLabel  = new QLabel;
    PupilCircumferenceMaxSlider = new SliderDouble;
    PupilCircumferenceMaxSlider->setPrecision(1);
    PupilCircumferenceMaxSlider->setDoubleRange(0, 500);
    PupilCircumferenceMaxSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilCircumferenceMaxSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumferenceMax(double)));

    // Pupil fraction

    QLabel *PupilAspectRatioMinTextBox = new QLabel;
    PupilAspectRatioMinTextBox->setText("<b>Fraction minimum:</b>");

    PupilAspectRatioMinLabel  = new QLabel;
    PupilAspectRatioMinSlider = new SliderDouble;
    PupilAspectRatioMinSlider->setPrecision(2);
    PupilAspectRatioMinSlider->setDoubleRange(0.0, 1.0);
    PupilAspectRatioMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilAspectRatioMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilAspectRatioMin(double)));

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

    QLabel *PupilHaarOffsetTextBox = new QLabel;
    PupilHaarOffsetTextBox->setText("<b>Pupil Haar-offset:</b>");

    PupilHaarOffsetLabel  = new QLabel;
    PupilHaarOffsetSlider = new QSlider;
    PupilHaarOffsetSlider->setRange(0, 50);
    PupilHaarOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilHaarOffsetSlider, SIGNAL(valueChanged(int)), this, SLOT(setPupilHaarOffset(int)));

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

    QWidget *ParameterWidget = new QWidget;
    QGridLayout *ParameterLayout = new QGridLayout(ParameterWidget);

    ParameterLayout->addWidget(PupilCircumferenceMaxTextBox,     1, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(PupilCircumferenceMinTextBox,     2, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(PupilAspectRatioMinTextBox,       3, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(EdgeIntensityOffsetTextBox,       4, 0, 1, 1, Qt::AlignRight);

    ParameterLayout->addWidget(CannyThresholdHighTextBox,        6, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(CannyThresholdLowTextBox,         7, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(CannyKernelSizeTextBox,           8, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(CannyBlurLevelTextBox,            9, 0, 1, 1, Qt::AlignRight);

    ParameterLayout->addWidget(AlphaPredictionTextBox,          11, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(AlphaAverageTextBox,             12, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(AlphaMomentumTextBox,            13, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(AlphaMiscellaneousTextBox,       14, 0, 1, 1, Qt::AlignRight);

    ParameterLayout->addWidget(ThresholdCircumferenceTextBox,   16, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(ThresholdAspectRatioTextBox,     17, 0, 1, 1, Qt::AlignRight);

    ParameterLayout->addWidget(PupilHaarOffsetTextBox,          19, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(GlintSizeTextBox,                20, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(CurvatureFactorTextBox,          21, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(CurvatureOffsetTextBox,          22, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(EdgeLengthMinimumTextBox,        23, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(EllipseFitNumberMaximumTextBox,  24, 0, 1, 1, Qt::AlignRight);
    ParameterLayout->addWidget(EllipseFitErrorMaximumTextBox,   25, 0, 1, 1, Qt::AlignRight);

    ParameterLayout->addWidget(TitleLimitTextBox,                0, 1, 1, 1, Qt::AlignCenter);
    ParameterLayout->addWidget(PupilCircumferenceMaxSlider,      1, 1);
    ParameterLayout->addWidget(PupilCircumferenceMinSlider,      2, 1);
    ParameterLayout->addWidget(PupilAspectRatioMinSlider,        3, 1);
    ParameterLayout->addWidget(EdgeIntensityOffsetSlider,        4, 1);
    ParameterLayout->addWidget(TitleCannyTextBox,                5, 1, 1, 1, Qt::AlignCenter);
    ParameterLayout->addWidget(CannyThresholdHighSlider,         6, 1);
    ParameterLayout->addWidget(CannyThresholdLowSlider,          7, 1);
    ParameterLayout->addWidget(CannyKernelSizeSlider,            8, 1);
    ParameterLayout->addWidget(CannyBlurLevelSlider,             9, 1);
    ParameterLayout->addWidget(TitleLearnTextBox,               10, 1, 1, 1, Qt::AlignCenter);
    ParameterLayout->addWidget(AlphaPredictionSlider,           11, 1);
    ParameterLayout->addWidget(AlphaAverageSlider,              12, 1);
    ParameterLayout->addWidget(AlphaMomentumSlider,             13, 1);
    ParameterLayout->addWidget(AlphaMiscellaneousSlider,        14, 1);
    ParameterLayout->addWidget(TitleChangeTextBox,              15, 1, 1, 1, Qt::AlignCenter);
    ParameterLayout->addWidget(ThresholdCircumferenceSlider,    16, 1);
    ParameterLayout->addWidget(ThresholdAspectRatioSlider,      17, 1);
    ParameterLayout->addWidget(TitleMiscTextBox,                18, 1, 1, 1, Qt::AlignCenter);
    ParameterLayout->addWidget(PupilHaarOffsetSlider,           19, 1);
    ParameterLayout->addWidget(GlintSizeSlider,                 20, 1);
    ParameterLayout->addWidget(CurvatureFactorSlider,           21, 1);
    ParameterLayout->addWidget(CurvatureOffsetSlider,           22, 1);
    ParameterLayout->addWidget(EdgeLengthMinimumSlider,         23, 1);
    ParameterLayout->addWidget(EllipseFitNumberMaximumSlider,   24, 1);
    ParameterLayout->addWidget(EllipseFitErrorMaximumSlider,    25, 1);

    ParameterLayout->addWidget(PupilCircumferenceMaxLabel,    1, 2);
    ParameterLayout->addWidget(PupilCircumferenceMinLabel,    2, 2);
    ParameterLayout->addWidget(PupilAspectRatioMinLabel,      3, 2);
    ParameterLayout->addWidget(EdgeIntensityOffsetLabel,      4, 2);

    ParameterLayout->addWidget(CannyThresholdHighLabel,       6, 2);
    ParameterLayout->addWidget(CannyThresholdLowLabel,        7, 2);
    ParameterLayout->addWidget(CannyKernelSizeLabel,          8, 2);
    ParameterLayout->addWidget(CannyBlurLevelLabel,           9, 2);

    ParameterLayout->addWidget(AlphaPredictionLabel,          11, 2);
    ParameterLayout->addWidget(AlphaAverageLabel,             12, 2);
    ParameterLayout->addWidget(AlphaMomentumLabel,            13, 2);
    ParameterLayout->addWidget(AlphaMiscellaneousLabel,       14, 2);

    ParameterLayout->addWidget(ThresholdCircumferenceLabel,   16, 2);
    ParameterLayout->addWidget(ThresholdAspectRatioLabel,     17, 2);

    ParameterLayout->addWidget(PupilHaarOffsetLabel,          19, 2);
    ParameterLayout->addWidget(GlintSizeLabel,                20, 2);
    ParameterLayout->addWidget(CurvatureFactorLabel,          21, 2);
    ParameterLayout->addWidget(CurvatureOffsetLabel,          22, 2);
    ParameterLayout->addWidget(EdgeLengthMinimumLabel,        23, 2);
    ParameterLayout->addWidget(EllipseFitNumberMaximumLabel,  24, 2);
    ParameterLayout->addWidget(EllipseFitErrorMaximumLabel,   25, 2);

    QScrollArea *ParameterScrollArea = new QScrollArea();
    ParameterScrollArea->setWidget(ParameterWidget);
    ParameterScrollArea->setWidgetResizable(true);

    QVBoxLayout *MainLayout = new QVBoxLayout;
    MainLayout->addWidget(ParameterScrollArea);

    setLayout(MainLayout);
}

ParameterWidget::~ParameterWidget()
{

}

eyePropertiesParameters ParameterWidget::getStructure()
{
    return mEyePropertiesParameters;
}

void ParameterWidget::setStructure(eyePropertiesParameters mEyePropertiesParametersNew)
{
    mEyePropertiesParameters = mEyePropertiesParametersNew;
    this->reset();
}

void ParameterWidget::reset()
{
    PupilCircumferenceMinSlider->setDoubleValue(mEyePropertiesParameters.circumferenceMin);
    PupilCircumferenceMinLabel ->setText(QString::number(mEyePropertiesParameters.circumferenceMin, 'f', 1));

    PupilCircumferenceMaxSlider->setDoubleValue(mEyePropertiesParameters.circumferenceMax);
    PupilCircumferenceMaxLabel ->setText(QString::number(mEyePropertiesParameters.circumferenceMax, 'f', 1));

    PupilAspectRatioMinSlider->setDoubleValue(mEyePropertiesParameters.aspectRatioMin);
    PupilAspectRatioMinLabel ->setText(QString::number(mEyePropertiesParameters.aspectRatioMin, 'f', 2));

    EdgeIntensityOffsetSlider->setDoubleValue(mEyePropertiesParameters.edgeIntensityOffset);
    EdgeIntensityOffsetLabel ->setText(QString::number(mEyePropertiesParameters.edgeIntensityOffset, 'f', 1));

    CannyThresholdLowSlider->setValue(mEyePropertiesParameters.cannyThresholdLow);
    CannyThresholdLowLabel ->setText(QString::number(mEyePropertiesParameters.cannyThresholdLow));

    CannyThresholdHighSlider->setValue(mEyePropertiesParameters.cannyThresholdHigh);
    CannyThresholdHighLabel ->setText(QString::number(mEyePropertiesParameters.cannyThresholdHigh));

    CannyBlurLevelSlider->setValue(mEyePropertiesParameters.cannyBlurLevel);
    CannyBlurLevelLabel ->setText(QString::number(mEyePropertiesParameters.cannyBlurLevel));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mEyePropertiesParameters.cannyKernelSize));
    CannyKernelSizeLabel ->setText(QString::number(mEyePropertiesParameters.cannyKernelSize));

    AlphaAverageSlider      ->setDoubleValue(mEyePropertiesParameters.alphaAverage);
    AlphaPredictionSlider   ->setDoubleValue(mEyePropertiesParameters.alphaPrediction);
    AlphaMiscellaneousSlider->setDoubleValue(mEyePropertiesParameters.alphaMiscellaneous);
    AlphaMomentumSlider     ->setDoubleValue(mEyePropertiesParameters.alphaMomentum);

    AlphaAverageLabel      ->setText(QString::number(mEyePropertiesParameters.alphaAverage, 'f', 3));
    AlphaPredictionLabel   ->setText(QString::number(mEyePropertiesParameters.alphaPrediction, 'f', 2));
    AlphaMiscellaneousLabel->setText(QString::number(mEyePropertiesParameters.alphaMiscellaneous, 'f', 2));
    AlphaMomentumLabel     ->setText(QString::number(mEyePropertiesParameters.alphaMomentum, 'f', 2));

    ThresholdCircumferenceSlider->setDoubleValue(mEyePropertiesParameters.circumferenceChangeThreshold);
    ThresholdCircumferenceLabel ->setText(QString::number(mEyePropertiesParameters.circumferenceChangeThreshold, 'f', 1));

    ThresholdAspectRatioSlider->setDoubleValue(mEyePropertiesParameters.aspectRatioChangeThreshold);
    ThresholdAspectRatioLabel ->setText(QString::number(mEyePropertiesParameters.aspectRatioChangeThreshold, 'f', 2));

    PupilHaarOffsetSlider->setValue(mEyePropertiesParameters.pupilOffset);
    PupilHaarOffsetLabel ->setText(QString::number(mEyePropertiesParameters.pupilOffset));

    GlintSizeSlider->setValue(round(0.5 * mEyePropertiesParameters.glintSize));
    GlintSizeLabel->setText(QString::number(mEyePropertiesParameters.glintSize));

    CurvatureFactorSlider->setDoubleValue(mEyePropertiesParameters.curvatureFactor);
    CurvatureFactorLabel ->setText(QString::number(mEyePropertiesParameters.curvatureFactor, 'f', 2));

    CurvatureOffsetSlider->setDoubleValue(mEyePropertiesParameters.curvatureOffsetMin);
    CurvatureOffsetLabel ->setText(QString::number(mEyePropertiesParameters.curvatureOffsetMin, 'f', 1));

    EdgeLengthMinimumSlider->setDoubleValue(mEyePropertiesParameters.edgeLengthMinimum);
    EdgeLengthMinimumLabel ->setText(QString::number(mEyePropertiesParameters.edgeLengthMinimum, 'f', 2));

    EllipseFitNumberMaximumSlider->setValue(mEyePropertiesParameters.ellipseFitNumberMaximum);
    EllipseFitNumberMaximumLabel ->setText(QString::number(mEyePropertiesParameters.ellipseFitNumberMaximum));

    EllipseFitErrorMaximumSlider->setDoubleValue(mEyePropertiesParameters.ellipseFitErrorMaximum);
    EllipseFitErrorMaximumLabel ->setText(QString::number(mEyePropertiesParameters.ellipseFitErrorMaximum, 'f', 1));
}

void ParameterWidget::setPupilCircumferenceMin(double value)
{
    if (mEyePropertiesParameters.circumferenceMax < value)
    {
        PupilCircumferenceMaxSlider->setDoubleValue(value);
    }

    mEyePropertiesParameters.circumferenceMin = value;
    PupilCircumferenceMinLabel->setText(QString::number(mEyePropertiesParameters.circumferenceMin, 'f', 1));
}

void ParameterWidget::setPupilCircumferenceMax(double value)
{
    if (mEyePropertiesParameters.circumferenceMin > value)
    {
        PupilCircumferenceMinSlider->setDoubleValue(value);
    }

    mEyePropertiesParameters.circumferenceMax = value;
    PupilCircumferenceMaxLabel->setText(QString::number(mEyePropertiesParameters.circumferenceMax, 'f', 1));
}

void ParameterWidget::setPupilAspectRatioMin(double value)
{
    mEyePropertiesParameters.aspectRatioMin = value;
    PupilAspectRatioMinLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setEdgeIntensityOffset(double value)
{
    mEyePropertiesParameters.edgeIntensityOffset = value;
    EdgeIntensityOffsetLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setEdgeLengthMinimum(double value)
{
    mEyePropertiesParameters.edgeLengthMinimum = value;
    EdgeLengthMinimumLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setCannyThresholdLow(int value)
{
    if (value > mEyePropertiesParameters.cannyThresholdHigh)
    {
        CannyThresholdHighSlider->setValue(value);
    }

    mEyePropertiesParameters.cannyThresholdLow = value;
    CannyThresholdLowLabel->setText(QString::number(value));

}

void ParameterWidget::setCannyThresholdHigh(int value)
{
    if (value < mEyePropertiesParameters.cannyThresholdLow)
    {
        CannyThresholdLowSlider->setValue(value);
    }

    mEyePropertiesParameters.cannyThresholdHigh = value;
    CannyThresholdHighLabel->setText(QString::number(value));
}

void ParameterWidget::setCannyKernelSize(int value)
{
    int newValue = 2 * value - 1;
    mEyePropertiesParameters.cannyKernelSize = newValue;
    CannyKernelSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCannyBlurLevel(int value)
{
    mEyePropertiesParameters.cannyBlurLevel = value;
    CannyBlurLevelLabel->setText(QString::number(value));
}

void ParameterWidget::setAlphaAverage(double value)
{
    mEyePropertiesParameters.alphaAverage = value;
    AlphaAverageLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setAlphaPrediction(double value)
{
    mEyePropertiesParameters.alphaPrediction = value;
    AlphaPredictionLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setAlphaMiscellaneous(double value)
{
    mEyePropertiesParameters.alphaMiscellaneous = value;
    AlphaMiscellaneousLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setAlphaMomentum(double value)
{
    mEyePropertiesParameters.alphaMomentum = value;
    AlphaMomentumLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdCircumference(double value)
{
    mEyePropertiesParameters.circumferenceChangeThreshold = value;
    ThresholdCircumferenceLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setThresholdAspectRatio(double value)
{
    mEyePropertiesParameters.aspectRatioChangeThreshold = value;
    ThresholdAspectRatioLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setPupilHaarOffset(int value)
{
    mEyePropertiesParameters.pupilOffset = value;
    PupilHaarOffsetLabel->setText(QString::number(value));
}

void ParameterWidget::setGlintSize(int value)
{
    int newValue = 2 * value;
    mEyePropertiesParameters.glintSize = newValue;
    GlintSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCurvatureFactor(double value)
{
    mEyePropertiesParameters.curvatureFactor = value;
    CurvatureFactorLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setCurvatureOffset(double value)
{
    mEyePropertiesParameters.curvatureOffsetMin = value;
    CurvatureOffsetLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setEllipseFitNumberMaximum(int value)
{
    mEyePropertiesParameters.ellipseFitNumberMaximum = value;
    EllipseFitNumberMaximumLabel->setText(QString::number(value));
}

void ParameterWidget::setEllipseFitErrorMaximum(double value)
{
    mEyePropertiesParameters.ellipseFitErrorMaximum = value;
    EllipseFitErrorMaximumLabel->setText(QString::number(value, 'f', 1));
}


