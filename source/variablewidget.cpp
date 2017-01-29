#include "headers/variablewidget.h"

VariableWidget::VariableWidget(QWidget *parent) : QWidget(parent)
{
    // Circumference

    QLabel *CircumferenceTextBox = new QLabel;
    CircumferenceTextBox->setText("<b> Circumference:</b>");

    CircumferenceLabel  = new QLabel();
    CircumferenceSlider = new SliderDouble();
    CircumferenceSlider->setPrecision(1);
    CircumferenceSlider->setDoubleRange(0, 500);
    CircumferenceSlider->setOrientation(Qt::Horizontal);

    // Aspect ratio

    QLabel *AspectRatioTextBox = new QLabel;
    AspectRatioTextBox->setText("<b> Aspect ratio:</b>");

    AspectRatioLabel  = new QLabel();
    AspectRatioSlider = new SliderDouble();
    AspectRatioSlider->setPrecision(2);
    AspectRatioSlider->setDoubleRange(0, 1.0);
    AspectRatioSlider->setOrientation(Qt::Horizontal);

    QObject::connect(CircumferenceSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCircumference(double)));
    QObject::connect(AspectRatioSlider,   SIGNAL(doubleValueChanged(double)), this, SLOT(setAspectRatio(double)));

    // Title

    QLabel *TitleWidget = new QLabel;
    TitleWidget->setText("<b>Real-time variables</b>");

    // Set-up layout

    QGridLayout *MainLayout = new QGridLayout;

    MainLayout->addWidget(TitleWidget,          0, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(CircumferenceTextBox, 1, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CircumferenceSlider,  1, 1);
    MainLayout->addWidget(CircumferenceLabel,   1, 2);
    MainLayout->addWidget(AspectRatioTextBox,   2, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AspectRatioSlider,    2, 1);
    MainLayout->addWidget(AspectRatioLabel,     2 ,2);

    MainLayout->setColumnStretch(0,1);
    MainLayout->setColumnStretch(1,3);
    MainLayout->setColumnStretch(2,1);

    setLayout(MainLayout);
}

VariableWidget::~VariableWidget()
{


}

detectionVariables VariableWidget::getStructure()
{
    return mDetectionVariables;
}

void VariableWidget::setWidgets(const detectionVariables& mDetectionVariables)
{
    CircumferenceSlider->setDoubleValue(mDetectionVariables.circumferencePrediction);
    CircumferenceLabel->setText(QString::number(mDetectionVariables.circumferencePrediction, 'f', 1));

    AspectRatioSlider->setDoubleValue(mDetectionVariables.aspectRatioPrediction);
    AspectRatioLabel->setText(QString::number(mDetectionVariables.aspectRatioPrediction, 'f', 2));
}

void VariableWidget::resetStructure(const detectionParameters& mDetectionParameters)
{
    mDetectionVariables.aspectRatioAverage    = 1.0;
    mDetectionVariables.aspectRatioExact      = 0.0;
    mDetectionVariables.aspectRatioMomentum   = 0.0;
    mDetectionVariables.aspectRatioPrediction = 1.0;

    mDetectionVariables.circumferenceAverage    = 0.5 * (mDetectionParameters.circumferenceMax + mDetectionParameters.circumferenceMin);
    mDetectionVariables.circumferenceExact      = 0;
    mDetectionVariables.circumferenceMomentum   = 0;
    mDetectionVariables.circumferencePrediction = mDetectionVariables.circumferenceAverage;

    mDetectionVariables.heightAverage    = mDetectionVariables.circumferencePrediction / M_PI;
    mDetectionVariables.heightPrediction = mDetectionVariables.heightAverage;
    mDetectionVariables.heightMomentum   = 0;

    mDetectionVariables.radiusMomentum   = 0;
    mDetectionVariables.radiusPrediction = 0.5 * mDetectionVariables.circumferencePrediction / M_PI;

    mDetectionVariables.searchRadius = 0.5 * Parameters::eyeAOIWdth;

    mDetectionVariables.thresholdCircumferenceChange = mDetectionParameters.circumferenceMax;
    mDetectionVariables.thresholdAspectRatioChange   = 1.0;

    mDetectionVariables.widthAverage    = mDetectionVariables.circumferencePrediction / M_PI;
    mDetectionVariables.widthPrediction = mDetectionVariables.widthAverage;
    mDetectionVariables.widthMomentum   = 0;

    mDetectionVariables.xPosAbsolute  = 0;
    mDetectionVariables.xPosExact     = 0;
    mDetectionVariables.xPosPredicted = 0.5 * Parameters::eyeAOIWdth;
    mDetectionVariables.xVelocity     = 0;

    mDetectionVariables.yPosAbsolute  = 0;
    mDetectionVariables.yPosExact     = 0;
    mDetectionVariables.yPosPredicted = 0.5 * Parameters::eyeAOIHght;
    mDetectionVariables.yVelocity     = 0;

    mDetectionVariables.priorCertainty = certaintyLowerLimit;

    setWidgets(mDetectionVariables);
}

void VariableWidget::setCircumference(double value)
{
    if (!Parameters::ONLINE_PROCESSING)
    {
        mDetectionVariables.circumferencePrediction = value;
        CircumferenceLabel->setText(QString::number(value, 'f', 1));
    }
}

void VariableWidget::setAspectRatio(double value)
{
    if (!Parameters::ONLINE_PROCESSING)
    {
        mDetectionVariables.aspectRatioPrediction = value;
        AspectRatioLabel->setText(QString::number(value, 'f', 2));
    }
}
