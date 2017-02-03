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
    CircumferenceSlider->setDoubleValue(mDetectionVariables.predictionCircumference);
    CircumferenceLabel->setText(QString::number(mDetectionVariables.predictionCircumference, 'f', 1));

    AspectRatioSlider->setDoubleValue(mDetectionVariables.predictionAspectRatio);
    AspectRatioLabel->setText(QString::number(mDetectionVariables.predictionAspectRatio, 'f', 2));
}

void VariableWidget::resetStructure(const detectionParameters& mDetectionParameters)
{
    mDetectionVariables.averageAspectRatio   = 0.9; // close to perfect circle
    mDetectionVariables.averageCircumference = 0.5 * (mDetectionParameters.circumferenceMax + mDetectionParameters.circumferenceMin);
    mDetectionVariables.averageWidth         = mDetectionVariables.predictionCircumference / M_PI;
    mDetectionVariables.averageHeight        = mDetectionVariables.predictionCircumference / M_PI;
    mDetectionVariables.averageIntensity     = 255;

    mDetectionVariables.exactAspectRatio   = 0;
    mDetectionVariables.exactCircumference = 0;
    mDetectionVariables.exactXPos          = 0;
    mDetectionVariables.exactYPos          = 0;

    mDetectionVariables.predictionAspectRatio   = mDetectionVariables.averageAspectRatio;
    mDetectionVariables.predictionCircumference = mDetectionVariables.averageCircumference;
    mDetectionVariables.predictionWidth         = mDetectionVariables.averageWidth;
    mDetectionVariables.predictionHeight        = mDetectionVariables.averageHeight;
    mDetectionVariables.predictionXPos          = 0.5 * Parameters::eyeAOI.wdth; // centre of image
    mDetectionVariables.predictionYPos          = 0.5 * Parameters::eyeAOI.hght;

    mDetectionVariables.momentumAspectRatio   = 0;
    mDetectionVariables.momentumCircumference = 0;
    mDetectionVariables.momentumWidth         = 0;
    mDetectionVariables.momentumHeight        = 0;
    mDetectionVariables.momentumXPos          = 0;
    mDetectionVariables.momentumYPos          = 0;

    mDetectionVariables.absoluteXPos  = 0;
    mDetectionVariables.absoluteYPos  = 0;

    int imgSize;
    if (Parameters::eyeAOI.wdth > Parameters::eyeAOI.hght) { imgSize = Parameters::eyeAOI.wdth; }
    else                                                   { imgSize = Parameters::eyeAOI.hght; }

    mDetectionVariables.changeThresholdAspectRatio   = 1.0;
    mDetectionVariables.changeThresholdCircumference = mDetectionParameters.circumferenceMax;
    mDetectionVariables.changeThresholdPosition      = imgSize;
    mDetectionVariables.searchRadius                 = imgSize;
    mDetectionVariables.curvatureOffset              = 180;

    mDetectionVariables.certaintyPosition = 0;
    mDetectionVariables.certaintyFeatures = 0;

    setWidgets(mDetectionVariables);
}

void VariableWidget::setCircumference(double value)
{
    if (!Parameters::ONLINE_PROCESSING)
    {
        mDetectionVariables.predictionCircumference = value;
        CircumferenceLabel->setText(QString::number(value, 'f', 1));
    }
}

void VariableWidget::setAspectRatio(double value)
{
    if (!Parameters::ONLINE_PROCESSING)
    {
        mDetectionVariables.predictionAspectRatio = value;
        AspectRatioLabel->setText(QString::number(value, 'f', 2));
    }
}
