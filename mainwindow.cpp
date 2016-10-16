//  Copyright (C) 2016  Terence Brouns

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#include "headers/mainwindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
    // Load data directory name

    {
        QSettings settings("config.ini", QSettings::IniFormat);
        dataDirectory = (settings.value("dataDirectory", "").toString()).toStdString();
    }

    // Get current date

    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(currentDate, 80, "%Y_%m_%d", timeinfo);

    // Initialize default values

    mUEyeOpencvCam.setDeviceInfo(5129, 5445);

    APP_EXIT = false;
    APP_RUNNING = true;
    cameraAOIFractionHghtDefaultLeft = 0.19;
    cameraAOIFractionHghtDefaultRght = 0.22;
    cameraAOIFractionHght = cameraAOIFractionHghtDefaultLeft;
    cameraAOIFractionWdthDefaultLeft = 0.25;
    cameraAOIFractionWdthDefaultRght = 0.31;
    cameraAOIFractionWdth = cameraAOIFractionWdthDefaultLeft;
    cameraAOIFractionXPosDefaultLeft = 0.20;
    cameraAOIFractionXPosDefaultRght = 0.52;
    cameraAOIFractionXPos = cameraAOIFractionXPosDefaultLeft;
    cameraAOIFractionYPosDefaultLeft = 0.41;
    cameraAOIFractionYPosDefaultRght = 0.37;
    cameraAOIFractionYPos = cameraAOIFractionYPosDefaultLeft;
    cameraAOIHghtMin = 4;
    cameraAOIHghtStepSize = 2;
    cameraAOIWdthMin = 32;
    cameraAOIWdthStepSize = 4;
    cameraFrameRateDesired = 250;
    cameraPixelClock = 24;
    cameraSubSamplingFactor = 2;
    camImageHght = 200;
    camImageWdth = 480; // size of image in widget
    PROCESSING_ALL_IMAGES = false;
    PROCESSING_ALL_TRIALS = false;
    editImageIndex = 0;
    editImageTotal = 0;
    experimentIndex = 0;
    EXPERIMENT_TRIAL_RECORDING = false;
    eyeAOIHghtFraction = 0.50;
    eyeAOIHghtMin = 75;
    eyeAOIWdthFraction = 0.50;
    eyeAOIWdthMin = 100;
    eyeImageHght = 200;
    eyeImageWdth = 320;
    FLASH_STANDBY = false;
    frameCount = 0;
    GAIN_AUTO = true;
    GAIN_BOOST = false;
    guiUpdateFrequency = 30;
    relativeTime = 0;
    SAVE_EYE_IMAGE = false;
    startTime = 0;
    subjectIdentifier = "";
    trialIndex = 0;
    trialTimeLength = 1000;

    Parameters::cameraXResolution = 1280;
    Parameters::cameraYResolution = 1024;

    Parameters::CAMERA_RUNNING = false;
    Parameters::REALTIME_PROCESSING = true;

    // Grab parameters from ini file

    LastUsedSettingsFileName = "config_user.ini";

    CamEyeAOIWdthSlider = new SliderDouble;
    CamEyeAOIWdthSlider->setPrecision(2);
    CamEyeAOIWdthSlider->setDoubleRange(0, 1.0);
    CamEyeAOIWdthSlider->setOrientation(Qt::Horizontal);

    CamEyeAOIHghtSlider = new SliderDouble;
    CamEyeAOIHghtSlider->setPrecision(2);
    CamEyeAOIHghtSlider->setDoubleRange(0, 1.0);
    CamEyeAOIHghtSlider->setOrientation(Qt::Vertical);
    CamEyeAOIHghtSlider->setInvertedAppearance(true);

    CamEyeAOIXPosSlider = new SliderDouble;
    CamEyeAOIXPosSlider->setPrecision(2);
    CamEyeAOIXPosSlider->setDoubleRange(0, 1.0);
    CamEyeAOIXPosSlider->setOrientation(Qt::Horizontal);

    CamEyeAOIYPosSlider = new SliderDouble;
    CamEyeAOIYPosSlider->setPrecision(2);
    CamEyeAOIYPosSlider->setDoubleRange(0, 1.0);
    CamEyeAOIYPosSlider->setOrientation(Qt::Vertical);
    CamEyeAOIYPosSlider->setInvertedAppearance(true);

    loadSettings(LastUsedSettingsFileName);

    cameraAOIWdthMax = Parameters::cameraXResolution / (double) cameraSubSamplingFactor; // maximum possible AOI size
    cameraAOIHghtMax = Parameters::cameraYResolution / (double) cameraSubSamplingFactor;

    updateCamAOIx();
    updateCamAOIy();

    // Set pupil parameters

    mEyePropertiesVariables.pupilCircumference = mEyePropertiesParameters.pupilCircfIni;

    mEyePropertiesVariables.pupilFraction = mEyePropertiesParameters.pupilFractIni;

    mEyePropertiesVariables.edgeIntensity = mEyePropertiesParameters.edgeIntensityIni;

    mEyePropertiesVariables.momentumCircumference = 0;
    mEyePropertiesVariables.momentumFraction = 0;
    mEyePropertiesVariables.momentumRadius = 0;

    mEyePropertiesVariables.xPosAbs = 0;
    mEyePropertiesVariables.yPosAbs = 0;

    // Set search area parameters

    mEyePropertiesVariables.xPosPredicted = Parameters::eyeAOIWdth;
    mEyePropertiesVariables.yPosPredicted = Parameters::eyeAOIHght;
    mEyePropertiesVariables.searchRadius = Parameters::eyeAOIWdth;

    mEyePropertiesVariables.xVelocity = 0;
    mEyePropertiesVariables.yVelocity = 0;

    // Drawing parameters

    Parameters::drawFlags.haar = true;
    Parameters::drawFlags.edge = true;
    Parameters::drawFlags.elps = true;

    // Camera feed

    cv::Mat imgCam(camImageWdth, camImageWdth, CV_8UC3, cv::Scalar(150, 150, 150));
    cv::Mat imgEye(eyeImageWdth, eyeImageHght, CV_8UC3, cv::Scalar(150, 150, 150));

    CamQImage = new QImageOpenCV(1);
    CamQImage->setSize(camImageWdth, camImageHght);
    CamQImage->setEyeAOI(Parameters::eyeAOIXPos, Parameters::eyeAOIYPos, Parameters::eyeAOIWdth, Parameters::eyeAOIHght);
    CamQImage->setFlashAOI(Parameters::flashAOIXPos, Parameters::flashAOIYPos, Parameters::flashAOIWdth, Parameters::flashAOIHght);

    CamQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    CamQImage->loadImage(imgCam);
    CamQImage->setImage();
    QObject::connect(CamQImage, SIGNAL(updateImage()), this , SLOT(updateRawImage()));

    EyeQImage  = new QImageOpenCV(2);
    EyeQImage->setSize(eyeImageWdth, eyeImageHght);
    EyeQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    EyeQImage->loadImage(imgEye);
    EyeQImage->setImage();
    QObject::connect(EyeQImage, SIGNAL(imageMouseClick(double, double)), this, SLOT(setPupilPosition(double, double)));

    // Cam AOI sliders

    CamEyeAOIWdthSlider->setDoubleValue(cameraAOIFractionWdth);
    QObject::connect(CamEyeAOIWdthSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIWdth(double)));

    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHght);
    QObject::connect(CamEyeAOIHghtSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIHght(double)));

    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPos);
    QObject::connect(CamEyeAOIXPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIXPos(double)));

    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPos);
    QObject::connect(CamEyeAOIYPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIYPos(double)));

    // Eye AOI sliders

    EyeWdthROISlider = new SliderDouble;
    EyeWdthROISlider->setPrecision(2);
    EyeWdthROISlider->setDoubleRange(0, 1.0);
    EyeWdthROISlider->setDoubleValue(eyeAOIWdthFraction);
    EyeWdthROISlider->setOrientation(Qt::Horizontal);
    QObject::connect(EyeWdthROISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEyeROIWdth(double)));

    EyeHghtROISlider = new SliderDouble;
    EyeHghtROISlider->setPrecision(2);
    EyeHghtROISlider->setDoubleRange(0, 1.0);
    EyeHghtROISlider->setDoubleValue(eyeAOIHghtFraction);
    EyeHghtROISlider->setOrientation(Qt::Vertical);
    QObject::connect(EyeHghtROISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEyeROIHght(double)));

    //

    QPushButton* AOICropButton = new QPushButton("&Crop AOI");
    QObject::connect(AOICropButton, SIGNAL(clicked()), this, SLOT(cropAOI()));

    // Check boxes for draw functions

    QLabel *HaarTextBox = new QLabel;
    HaarTextBox->setText("<b>Box: </b>");

    QLabel *EdgeTextBox = new QLabel;
    EdgeTextBox->setText("<b>Edges: </b>");

    QLabel *ElpsTextBox = new QLabel;
    ElpsTextBox->setText("<b>Ellipse: </b>");

    QCheckBox *HaarCheckBox = new QCheckBox;
    QCheckBox *EdgeCheckBox = new QCheckBox;
    QCheckBox *ElpsCheckBox = new QCheckBox;

    HaarCheckBox->setChecked(true);
    EdgeCheckBox->setChecked(true);
    ElpsCheckBox->setChecked(true);

    QObject::connect(HaarCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setDrawHaar(int)));
    QObject::connect(EdgeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setDrawEdge(int)));
    QObject::connect(ElpsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setDrawElps(int)));

    QHBoxLayout *DrawFunctionsLayout = new QHBoxLayout;
    DrawFunctionsLayout->addWidget(HaarTextBox);
    DrawFunctionsLayout->addWidget(HaarCheckBox);
    DrawFunctionsLayout->addWidget(EdgeTextBox);
    DrawFunctionsLayout->addWidget(EdgeCheckBox);
    DrawFunctionsLayout->addWidget(ElpsTextBox);
    DrawFunctionsLayout->addWidget(ElpsCheckBox);
    DrawFunctionsLayout->addStretch();

    QLabel *RealTimeEyeTrackingTextBox = new QLabel;
    RealTimeEyeTrackingTextBox->setText("<b>Real-time eye tracking: </b>");

    QCheckBox *RealTimeEyeTrackingCheckBox = new QCheckBox;
    RealTimeEyeTrackingCheckBox->setChecked(true);
    QObject::connect(RealTimeEyeTrackingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRealTimeEyeTracking(int)));

    QLabel *ReviewModeTextBox = new QLabel;
    ReviewModeTextBox->setText("<b>Review mode: </b>");

    QCheckBox *ReviewModeCheckBox = new QCheckBox;
    ReviewModeCheckBox->setChecked(false);
    QObject::connect(ReviewModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setReviewMode(int)));

    QHBoxLayout *RealTimeEyeTrackingLayout = new QHBoxLayout;
    RealTimeEyeTrackingLayout->addStretch();
    RealTimeEyeTrackingLayout->addWidget(RealTimeEyeTrackingTextBox);
    RealTimeEyeTrackingLayout->addWidget(RealTimeEyeTrackingCheckBox);
    RealTimeEyeTrackingLayout->addWidget(ReviewModeTextBox);
    RealTimeEyeTrackingLayout->addWidget(ReviewModeCheckBox);
    RealTimeEyeTrackingLayout->addStretch();

    QPushButton *AOILeftEyeButton = new QPushButton("&Left eye");
    QObject::connect(AOILeftEyeButton, SIGNAL(clicked()), this, SLOT(setAOILeftEye()));

    QPushButton *AOIRghtEyeButton = new QPushButton("&Right eye");
    QObject::connect(AOIRghtEyeButton, SIGNAL(clicked()), this, SLOT(setAOIRghtEye()));

    QHBoxLayout *EyeDetectionsLayout = new QHBoxLayout;
    EyeDetectionsLayout->addStretch();
    EyeDetectionsLayout->addWidget(AOILeftEyeButton);
    EyeDetectionsLayout->addWidget(AOIRghtEyeButton);
    EyeDetectionsLayout->addWidget(AOICropButton);
    EyeDetectionsLayout->addStretch();

    QPushButton *ReviewLoadSessionButton = new QPushButton("Load session");
    QObject::connect(ReviewLoadSessionButton, SIGNAL(clicked(bool)), this, SLOT(loadReviewSession()));

    ReviewImageSlider = new QSlider;
    ReviewImageSlider->setRange(0, 0);
    ReviewImageSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ReviewImageSlider, SIGNAL(valueChanged(int)), this, SLOT(setReviewImageFrame(int)));

    ReviewImageFrameTextBox = new QLabel;
    ReviewImageFrameTextBox->setText("<b>0 / 0</b>");
    ReviewImageFrameTextBox->setAlignment(Qt::AlignCenter);
    ReviewImageFrameTextBox->setMinimumWidth(70);

    QPushButton *ReviewPrevImageButton = new QPushButton("<");
    QObject::connect(ReviewPrevImageButton, SIGNAL(clicked(bool)), this, SLOT(prevReviewImage()));

    QPushButton *ReviewNextImageButton = new QPushButton(">");
    QObject::connect(ReviewNextImageButton, SIGNAL(clicked(bool)), this, SLOT(nextReviewImage()));

    QLabel *ReviewPupilDetectionTextBox = new QLabel;
    ReviewPupilDetectionTextBox->setText("<b>Detect pupil: </b> ");

    QPushButton *ReviewPupilDetectionOneButton = new QPushButton("Current");
    QObject::connect(ReviewPupilDetectionOneButton, SIGNAL(clicked(bool)), this, SLOT(detectPupilOneFrame()));

    QPushButton *ReviewPupilDetectionAllFramesButton = new QPushButton("All frames");
    QObject::connect(ReviewPupilDetectionAllFramesButton, SIGNAL(clicked(bool)), this, SLOT(detectPupilAllFrames()));

    QPushButton *ReviewPupilDetectionAllTrialsButton = new QPushButton("All trials");
    QObject::connect(ReviewPupilDetectionAllTrialsButton, SIGNAL(clicked(bool)), this, SLOT(detectPupilAllTrials()));

    QPushButton *SavePupilDataButton = new QPushButton("Save");
    QObject::connect(SavePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(onSavePupilData()));

    EyeTrackingReviewWidget = new QWidget;
    QHBoxLayout *EyeTrackingReviewLayout = new QHBoxLayout(EyeTrackingReviewWidget);
    EyeTrackingReviewLayout->addStretch();
    EyeTrackingReviewLayout->addWidget(ReviewLoadSessionButton);
    EyeTrackingReviewLayout->addWidget(SavePupilDataButton);
    EyeTrackingReviewLayout->addWidget(ReviewPrevImageButton);
    EyeTrackingReviewLayout->addWidget(ReviewImageSlider);
    EyeTrackingReviewLayout->addWidget(ReviewNextImageButton);
    EyeTrackingReviewLayout->addWidget(ReviewImageFrameTextBox);
    EyeTrackingReviewLayout->addWidget(ReviewPupilDetectionTextBox);
    EyeTrackingReviewLayout->addWidget(ReviewPupilDetectionOneButton);
    EyeTrackingReviewLayout->addWidget(ReviewPupilDetectionAllFramesButton);
    EyeTrackingReviewLayout->addWidget(ReviewPupilDetectionAllTrialsButton);
    EyeTrackingReviewLayout->addStretch();

    EyeTrackingReviewWidget->setVisible(false);

    QLabel* ReviewTrialTitle = new QLabel;
    ReviewTrialTitle->setText("<b>Review mode - Trial:</b>");

    ReviewTrialSpinBox = new QSpinBox;
    ReviewTrialSpinBox->setValue(1);
    ReviewTrialSpinBox->setMinimum(1);
    ReviewTrialSpinBox->setAlignment(Qt::AlignRight);

    ReviewTrialSlider = new QSlider;
    ReviewTrialSlider->setValue(1);
    ReviewTrialSlider->setMinimum(1);
    ReviewTrialSlider->setOrientation(Qt::Horizontal);

    QObject::connect(ReviewTrialSlider, SIGNAL(valueChanged(int)),this,SLOT(changeReviewSession(int)));

    QObject::connect(ReviewTrialSlider, SIGNAL(valueChanged(int)),ReviewTrialSpinBox,SLOT(setValue(int)));
    QObject::connect(ReviewTrialSpinBox, SIGNAL(valueChanged(int)),ReviewTrialSlider,SLOT(setValue(int)));

    QGridLayout *ReviewSessionTitleLayout = new QGridLayout;
    ReviewSessionTitleLayout->addWidget(ReviewTrialTitle, 0, 1);
    ReviewSessionTitleLayout->addWidget(ReviewTrialSpinBox, 0, 2);
    ReviewSessionTitleLayout->addWidget(ReviewTrialSlider, 0, 3);
    ReviewSessionTitleLayout->setColumnStretch(0, 1);
    ReviewSessionTitleLayout->setColumnStretch(1, 3);
    ReviewSessionTitleLayout->setColumnStretch(2, 1);
    ReviewSessionTitleLayout->setColumnStretch(3, 3);
    ReviewSessionTitleLayout->setColumnStretch(4, 1);

    QWidget *CameraSettings = new QWidget;
    QGridLayout* CameraOutputLayout = new QGridLayout(CameraSettings);

    CameraOutputLayout->addLayout(ReviewSessionTitleLayout, 0, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(EyeTrackingReviewWidget, 5, 1, 1, 4, Qt::AlignCenter);

    CameraOutputLayout->addWidget(CamEyeAOIXPosSlider, 0, 2);
    CameraOutputLayout->addWidget(CamEyeAOIYPosSlider, 1, 1);
    CameraOutputLayout->addWidget(CamQImage, 1, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(CamEyeAOIWdthSlider, 2, 2);
    CameraOutputLayout->addWidget(CamEyeAOIHghtSlider, 1, 3);
    CameraOutputLayout->addWidget(EyeHghtROISlider, 1, 5);
    CameraOutputLayout->addWidget(EyeQImage, 1, 4, Qt::AlignCenter);
    CameraOutputLayout->addWidget(EyeWdthROISlider, 2, 4);
    CameraOutputLayout->addLayout(RealTimeEyeTrackingLayout, 4, 2);
    CameraOutputLayout->addLayout(EyeDetectionsLayout, 3, 2);
    CameraOutputLayout->addLayout(DrawFunctionsLayout, 3, 4);

    CameraOutputLayout->setColumnStretch(0, 1);
    CameraOutputLayout->setColumnStretch(6, 1);

    // Camera settings

    QLabel *CameraPixelClockTextBox = new QLabel;
    CameraPixelClockTextBox->setText("<b>Pixel clock (MHz): </b>");

    CameraPixelClockSlider = new QSlider;
    CameraPixelClockSlider->setRange(0, 0);
    CameraPixelClockSlider->setValue(0);
    CameraPixelClockSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraPixelClockSlider, SIGNAL(valueChanged(int)), this, SLOT(setCameraPixelClock(int)));

    CameraPixelClockLabel = new QLabel;
    CameraPixelClockLabel->setText(QString::number(0));

    // Frame rate

    QLabel *CameraFrameRateTextBox = new QLabel;
    CameraFrameRateTextBox->setText("<b>Frame-rate (Hz): </b>");

    CameraFrameRateSlider = new SliderDouble;
    CameraFrameRateSlider->setPrecision(1);
    CameraFrameRateSlider->setDoubleRange(0, 0);
    CameraFrameRateSlider->setDoubleValue(0);
    CameraFrameRateSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraFrameRateSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCameraFrameRate(double)));

    CameraFrameRateLabel = new QLabel;
    CameraFrameRateLabel->setText(QString::number(0.0, 'f', 1));

    QLabel *CameraFrameRateDesiredTextBox = new QLabel;
    CameraFrameRateDesiredTextBox->setText("<b>Desired frame-rate (Hz): </b>");

    CameraFrameRateDesiredSpinBox = new QSpinBox;
    CameraFrameRateDesiredSpinBox->setRange(0, 500);
    CameraFrameRateDesiredSpinBox->setValue(cameraFrameRateDesired);
    CameraFrameRateDesiredSpinBox->setAlignment(Qt::AlignRight);

    // Exposure

    QLabel *CameraExposureTextBox = new QLabel;
    CameraExposureTextBox->setText("<b>Exposure time (ms): </b>");

    CameraExposureSlider = new SliderDouble;
    CameraExposureSlider->setPrecision(2);
    CameraExposureSlider->setDoubleRange(0.0, 0.0);
    CameraExposureSlider->setDoubleValue(0.0);
    CameraExposureSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraExposureSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCameraExposure(double)));

    CameraExposureLabel = new QLabel;
    CameraExposureLabel->setText(QString::number(0.0, 'f', 2));

    // Black level

    QLabel *CameraBlackLevelOffsetTextBox = new QLabel;
    CameraBlackLevelOffsetTextBox->setText("<b>Black level correction: </b>");

    CameraBlackLevelOffsetSlider = new QSlider;
    CameraBlackLevelOffsetSlider->setRange(0, 0);
    CameraBlackLevelOffsetSlider->setValue(0);
    CameraBlackLevelOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraBlackLevelOffsetSlider, SIGNAL(valueChanged(int)), this, SLOT(setCameraBlackLevelOffset(int)));

    CameraBlackLevelOffsetLabel = new QLabel;
    CameraBlackLevelOffsetLabel->setText(QString::number(0));

    QLabel *CameraBlackLevelModeTextBox = new QLabel;
    CameraBlackLevelModeTextBox->setText("<b>Auto: </b>");

    QCheckBox *CameraBlackLevelModeCheckBox = new QCheckBox;
    CameraBlackLevelModeCheckBox->setChecked(true);
    QObject::connect(CameraBlackLevelModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setCameraBlackLevelMode(int)));

    //

    QLabel *CameraHardwareGainTextBox = new QLabel;
    CameraHardwareGainTextBox->setText("<b>Gain: </b>");

    CameraHardwareGainSlider = new QSlider;
    CameraHardwareGainSlider->setRange(0, 100);
    CameraHardwareGainSlider->setValue(0);
    CameraHardwareGainSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraHardwareGainSlider, SIGNAL(valueChanged(int)), this, SLOT(setCameraHardwareGain(int)));

    CameraHardwareGainLabel = new QLabel;
    CameraHardwareGainLabel->setText(QString::number(0));

    QLabel *CameraHardwareGainAutoTextBox = new QLabel;
    CameraHardwareGainAutoTextBox->setText("<b>Auto: </b>");

    CameraHardwareGainAutoCheckBox = new QCheckBox;
    CameraHardwareGainAutoCheckBox->setChecked(GAIN_AUTO);
    QObject::connect(CameraHardwareGainAutoCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setCameraAutoGain(int)));

    QLabel *CameraHardwareGainBoostTextBox = new QLabel;
    CameraHardwareGainBoostTextBox->setText("<b>Boost: </b>");

    QCheckBox *CameraHardwareGainBoostCheckBox = new QCheckBox;
    CameraHardwareGainBoostCheckBox->setChecked(GAIN_BOOST);
    QObject::connect(CameraHardwareGainBoostCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setCameraGainBoost(int)));

    //

    QLabel *CameraSubSamplingTextBox = new QLabel;
    CameraSubSamplingTextBox->setText("<b>Sub-sampling: </b>");

    QCheckBox *CameraSubSamplingCheckBox = new QCheckBox;

    if (cameraSubSamplingFactor == 2)
    {
        CameraSubSamplingCheckBox->setChecked(true);
    }
    else
    {
        CameraSubSamplingCheckBox->setChecked(false);
    }

    QObject::connect(CameraSubSamplingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setCameraSubSampling(int)));

    //

    CameraParametersWidget = new QWidget;
    QGridLayout *CameraParametersLayout = new QGridLayout(CameraParametersWidget);
    CameraParametersLayout->addWidget(CameraPixelClockTextBox, 0, 0);
    CameraParametersLayout->addWidget(CameraPixelClockSlider, 0, 1);
    CameraParametersLayout->addWidget(CameraPixelClockLabel, 0, 2);
    CameraParametersLayout->addWidget(CameraFrameRateTextBox, 1, 0);
    CameraParametersLayout->addWidget(CameraFrameRateSlider, 1, 1);
    CameraParametersLayout->addWidget(CameraFrameRateLabel, 1, 2);

    QHBoxLayout *CameraFrameRateDesiredLayout = new QHBoxLayout;
    CameraFrameRateDesiredLayout->addWidget(CameraFrameRateDesiredTextBox);
    CameraFrameRateDesiredLayout->addWidget(CameraFrameRateDesiredSpinBox);
    CameraFrameRateDesiredLayout->addStretch();

    CameraParametersLayout->addLayout(CameraFrameRateDesiredLayout, 2, 1);

    CameraParametersLayout->addWidget(CameraExposureTextBox, 3, 0);
    CameraParametersLayout->addWidget(CameraExposureSlider, 3, 1);
    CameraParametersLayout->addWidget(CameraExposureLabel, 3, 2);
    CameraParametersLayout->addWidget(CameraBlackLevelOffsetTextBox, 4, 0);
    CameraParametersLayout->addWidget(CameraBlackLevelOffsetSlider, 4, 1);
    CameraParametersLayout->addWidget(CameraBlackLevelOffsetLabel, 4, 2);

    QHBoxLayout *CameraBlackLevelModeLayout = new QHBoxLayout;
    CameraBlackLevelModeLayout->addWidget(CameraBlackLevelModeTextBox);
    CameraBlackLevelModeLayout->addWidget(CameraBlackLevelModeCheckBox);
    CameraBlackLevelModeLayout->addStretch();

    CameraParametersLayout->addLayout(CameraBlackLevelModeLayout, 5, 1);

    CameraParametersLayout->addWidget(CameraHardwareGainTextBox, 6, 0);
    CameraParametersLayout->addWidget(CameraHardwareGainSlider, 6, 1);
    CameraParametersLayout->addWidget(CameraHardwareGainLabel, 6, 2);

    QHBoxLayout *CameraHardwareGainOptionsLayout = new QHBoxLayout;
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainAutoTextBox);
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainAutoCheckBox);
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainBoostTextBox);
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainBoostCheckBox);
    CameraHardwareGainOptionsLayout->addStretch();

    CameraParametersLayout->addLayout(CameraHardwareGainOptionsLayout, 7, 1);

    //    CameraParametersLayout->addWidget(CameraSubSamplingTextBox, 7, 0);
    //    CameraParametersLayout->addWidget(CameraSubSamplingCheckBox, 7, 1);

    // Parameter settings layout

    // Real-time parameter tracking

    // Pupil circumference

    QLabel *PupilCircfTextBox = new QLabel;
    PupilCircfTextBox->setText("<b>Pupil circumference:</b>");

    PupilCircfSliderDouble = new SliderDouble();
    PupilCircfSliderDouble->setPrecision(1);
    PupilCircfSliderDouble->setDoubleRange(mEyePropertiesParameters.pupilCircumferenceMin, mEyePropertiesParameters.pupilCircumferenceMax);
    PupilCircfSliderDouble->setDoubleValue(mEyePropertiesParameters.pupilCircfIni);
    PupilCircfSliderDouble->setOrientation(Qt::Horizontal);

    PupilCircfLabel = new QLabel();
    PupilCircfLabel->setText(QString::number(mEyePropertiesParameters.pupilCircfIni, 'f', 1));

    // Pupil fraction

    QLabel *PupilFractTextBox = new QLabel;
    PupilFractTextBox->setText("<b>Pupil fraction:</b>");

    PupilFractSliderDouble = new SliderDouble();
    PupilFractSliderDouble->setPrecision(2);
    PupilFractSliderDouble->setDoubleRange(mEyePropertiesParameters.pupilFractMin, 1.0);
    PupilFractSliderDouble->setDoubleValue(mEyePropertiesParameters.pupilFractIni);
    PupilFractSliderDouble->setOrientation(Qt::Horizontal);

    PupilFractLabel = new QLabel();
    PupilFractLabel->setText(QString::number(mEyePropertiesParameters.pupilFractIni, 'f', 2));

    // Edge intensity

    QLabel *EdgeIntensityTextBox = new QLabel;
    EdgeIntensityTextBox->setText("<b>Edge intensity:</b>");

    EdgeIntensitySliderDouble = new SliderDouble();
    EdgeIntensitySliderDouble->setPrecision(1);
    EdgeIntensitySliderDouble->setDoubleRange(0.0, 255.0);
    EdgeIntensitySliderDouble->setDoubleValue(mEyePropertiesParameters.edgeIntensityIni);
    EdgeIntensitySliderDouble->setOrientation(Qt::Horizontal);

    EdgeIntensityLabel = new QLabel();
    EdgeIntensityLabel->setText(QString::number(mEyePropertiesParameters.edgeIntensityIni, 'f', 1));

    // Review mode

    QObject::connect(PupilCircfSliderDouble, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumference(double)));
    QObject::connect(PupilFractSliderDouble, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilFraction(double)));
    QObject::connect(EdgeIntensitySliderDouble, SIGNAL(doubleValueChanged(double)), this, SLOT(setEdgeIntensity(double)));

    // Averages and Limits/Thresholds

    QLabel *ParametersTextBox = new QLabel;
    ParametersTextBox->setText("<b>Eye-tracking parameters</b>");
    ParametersTextBox->setAlignment(Qt::AlignCenter);

    // Pupil circumference

    QLabel *PupilCircumferenceIniTextBox = new QLabel;
    PupilCircumferenceIniTextBox->setText("<b>Circumference average:</b>");

    PupilCircumferenceIniSlider = new SliderDouble;
    PupilCircumferenceIniSlider->setPrecision(1);
    PupilCircumferenceIniSlider->setDoubleRange(0.1, 300);
    PupilCircumferenceIniSlider->setDoubleValue(mEyePropertiesParameters.pupilCircfIni);
    PupilCircumferenceIniSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilCircumferenceIniSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumferenceIni(double)));

    PupilCircumferenceIniLabel = new QLabel;
    PupilCircumferenceIniLabel->setText(QString::number(mEyePropertiesParameters.pupilCircfIni, 'f', 1));

    QLabel *PupilCircumferenceOffTextBox = new QLabel;
    PupilCircumferenceOffTextBox->setText("<b>Circumference offset:</b>");

    circumferenceOffset = mEyePropertiesParameters.pupilCircumferenceMax - mEyePropertiesParameters.pupilCircfIni;

    PupilCircumferenceOffSlider = new SliderDouble;
    PupilCircumferenceOffSlider->setPrecision(1);
    PupilCircumferenceOffSlider->setDoubleRange(0.1, mEyePropertiesParameters.pupilCircfIni);
    PupilCircumferenceOffSlider->setDoubleValue(circumferenceOffset);
    PupilCircumferenceOffSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilCircumferenceOffSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumferenceOff(double)));

    PupilCircumferenceOffLabel = new QLabel;
    PupilCircumferenceOffLabel->setText(QString::number(circumferenceOffset, 'f', 1));

    // Pupil fraction

    QLabel *PupilFractionIniTextBox = new QLabel;
    PupilFractionIniTextBox->setText("<b>Fraction average:</b>");

    PupilFractionIniSlider = new SliderDouble;
    PupilFractionIniSlider->setPrecision(2);
    PupilFractionIniSlider->setDoubleRange(mEyePropertiesParameters.pupilFractMin, 1.00);
    PupilFractionIniSlider->setDoubleValue(mEyePropertiesParameters.pupilFractIni);
    PupilFractionIniSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilFractionIniSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilFractionIni(double)));

    PupilFractionIniLabel = new QLabel;
    PupilFractionIniLabel->setText(QString::number(mEyePropertiesParameters.pupilFractIni, 'f', 2));

    QLabel *PupilFractionMinTextBox = new QLabel;
    PupilFractionMinTextBox->setText("<b>Fraction minimum:</b>");

    PupilFractionMinSlider = new SliderDouble;
    PupilFractionMinSlider->setPrecision(2);
    PupilFractionMinSlider->setDoubleRange(0.0, mEyePropertiesParameters.pupilFractIni);
    PupilFractionMinSlider->setDoubleValue(mEyePropertiesParameters.pupilFractMin);
    PupilFractionMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilFractionMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilFractionMin(double)));

    PupilFractionMinLabel = new QLabel;
    PupilFractionMinLabel->setText(QString::number(mEyePropertiesParameters.pupilFractMin, 'f', 2));

    // Edge intensity

    QLabel *EdgeIntensityOffsetTextBox = new QLabel;
    EdgeIntensityOffsetTextBox->setText("<b>Edge intensity offset:</b>");

    EdgeIntensityOffsetSlider = new SliderDouble();
    EdgeIntensityOffsetSlider->setPrecision(1);
    EdgeIntensityOffsetSlider->setDoubleRange(0.0, 255.0);
    EdgeIntensityOffsetSlider->setDoubleValue(mEyePropertiesParameters.edgeIntensityOffset);
    EdgeIntensityOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EdgeIntensityOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEdgeIntensityOffset(double)));

    EdgeIntensityOffsetLabel = new QLabel;
    EdgeIntensityOffsetLabel->setText(QString::number(mEyePropertiesParameters.edgeIntensityOffset, 'f', 1));

    // Sliders for canny edge parameters

    QLabel *CannyEdgeTextBox = new QLabel;
    CannyEdgeTextBox->setText("<b>Canny edge detection</b>");
    CannyEdgeTextBox->setAlignment(Qt::AlignCenter);

    QLabel *CannyLowerLimitTextBox = new QLabel;
    CannyLowerLimitTextBox->setText("<b>Lower limit:</b>");

    CannyLowerLimitSlider = new QSlider;
    CannyLowerLimitSlider->setRange(1, mEyePropertiesParameters.cannyUpperLimit);
    CannyLowerLimitSlider->setValue(mEyePropertiesParameters.cannyLowerLimit);
    CannyLowerLimitSlider->setOrientation(Qt::Horizontal);
    CannyLowerLimitSlider->setSingleStep(1);
    QObject::connect(CannyLowerLimitSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyLowerLimit(int)));

    CannyLowerLimitLabel = new QLabel;
    CannyLowerLimitLabel->setText(QString::number(mEyePropertiesParameters.cannyLowerLimit));

    QLabel *CannyUpperLimitTextBox = new QLabel;
    CannyUpperLimitTextBox->setText("<b>Upper limit:</b>");

    CannyUpperLimitSlider = new QSlider;
    CannyUpperLimitSlider->setRange(mEyePropertiesParameters.cannyLowerLimit, 4 * mEyePropertiesParameters.cannyLowerLimit);
    CannyUpperLimitSlider->setValue(mEyePropertiesParameters.cannyUpperLimit);
    CannyUpperLimitSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyUpperLimitSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyUpperLimit(int)));

    CannyUpperLimitLabel = new QLabel;
    CannyUpperLimitLabel->setText(QString::number(mEyePropertiesParameters.cannyUpperLimit));

    QLabel *CannyBlurLevelTextBox = new QLabel;
    CannyBlurLevelTextBox->setText("<b>Blur level:</b>");

    CannyBlurLevelSlider = new QSlider;
    CannyBlurLevelSlider->setRange(0, 10);
    CannyBlurLevelSlider->setValue(ceil(0.5 * mEyePropertiesParameters.cannyBlurLevel));
    CannyBlurLevelSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyBlurLevelSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyBlurLevel(int)));

    CannyBlurLevelLabel = new QLabel;
    CannyBlurLevelLabel->setText(QString::number(ceil(0.5 * mEyePropertiesParameters.cannyBlurLevel)));

    QLabel *CannyKernelSizeTextBox = new QLabel;
    CannyKernelSizeTextBox->setText("<b>Kernel size:</b>");

    CannyKernelSizeSlider = new QSlider;
    CannyKernelSizeSlider->setRange(2, 4);
    CannyKernelSizeSlider->setValue(ceil(0.5 * mEyePropertiesParameters.cannyKernelSize));
    CannyKernelSizeSlider->setOrientation(Qt::Horizontal);
    CannyKernelSizeSlider->setSingleStep(1);
    QObject::connect(CannyKernelSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyKernelSize(int)));

    CannyKernelSizeLabel = new QLabel;
    CannyKernelSizeLabel->setText(QString::number(mEyePropertiesParameters.cannyKernelSize));

    // Learning rate parameters

    QLabel *LearningRatesTextBox = new QLabel;
    LearningRatesTextBox->setText("<b>Learning rates</b>");
    LearningRatesTextBox->setAlignment(Qt::AlignCenter);

    QLabel *AlphaPupilTextBox = new QLabel;
    QLabel *AlphaSearchAreaTextBox = new QLabel;
    QLabel *AlphaGeneralTextBox = new QLabel;
    QLabel *AlphaMomentumTextBox = new QLabel;
    QLabel *AlphaVelocityTextBox = new QLabel;

    AlphaPupilTextBox->setText("<b>Pupil:</b>");
    AlphaSearchAreaTextBox->setText("<b>Search area:</b>");
    AlphaGeneralTextBox->setText("<b>Search radius:</b>");
    AlphaMomentumTextBox->setText("<b>Momentum:</b>");
    AlphaVelocityTextBox->setText("<b>Velocity:</b>");

    AlphaPupilSlider = new SliderDouble;
    AlphaSearchAreaSlider = new SliderDouble;
    AlphaGeneralSlider = new SliderDouble;
    AlphaMomentumSlider = new SliderDouble;
    AlphaVelocitySlider = new SliderDouble;

    AlphaPupilSlider->setPrecision(2);
    AlphaSearchAreaSlider->setPrecision(2);
    AlphaGeneralSlider->setPrecision(2);
    AlphaMomentumSlider->setPrecision(2);
    AlphaVelocitySlider->setPrecision(2);

    AlphaPupilSlider->setDoubleRange(0, 1.0);
    AlphaSearchAreaSlider->setDoubleRange(0, 1.0);
    AlphaGeneralSlider->setDoubleRange(0, 1.0);
    AlphaMomentumSlider->setDoubleRange(0, 1.0);
    AlphaVelocitySlider->setDoubleRange(0, 1.0);

    AlphaPupilSlider->setDoubleValue(mEyePropertiesParameters.alphaPupil);
    AlphaSearchAreaSlider->setDoubleValue(mEyePropertiesParameters.alphaPosition);
    AlphaGeneralSlider->setDoubleValue(mEyePropertiesParameters.alphaGeneral);
    AlphaMomentumSlider->setDoubleValue(mEyePropertiesParameters.alphaMomentum);
    AlphaVelocitySlider->setDoubleValue(mEyePropertiesParameters.alphaVelocity);

    AlphaPupilSlider->setOrientation(Qt::Horizontal);
    AlphaSearchAreaSlider->setOrientation(Qt::Horizontal);
    AlphaGeneralSlider->setOrientation(Qt::Horizontal);
    AlphaMomentumSlider->setOrientation(Qt::Horizontal);
    AlphaVelocitySlider->setOrientation(Qt::Horizontal);

    QObject::connect(AlphaPupilSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaPupil(double)));
    QObject::connect(AlphaSearchAreaSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaSearchArea(double)));
    QObject::connect(AlphaGeneralSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaGeneral(double)));
    QObject::connect(AlphaMomentumSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaMomentum(double)));
    QObject::connect(AlphaVelocitySlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAlphaVelocity(double)));

    AlphaPupilLabel = new QLabel;
    AlphaSearchAreaLabel = new QLabel;
    AlphaGeneralLabel = new QLabel;
    AlphaMomentumLabel = new QLabel;
    AlphaVelocityLabel = new QLabel;

    AlphaPupilLabel->setText(QString::number(mEyePropertiesParameters.alphaPupil, 'f', 2));
    AlphaSearchAreaLabel->setText(QString::number(mEyePropertiesParameters.alphaPosition, 'f', 2));
    AlphaGeneralLabel->setText(QString::number(mEyePropertiesParameters.alphaGeneral, 'f', 2));
    AlphaMomentumLabel->setText(QString::number(mEyePropertiesParameters.alphaMomentum, 'f', 2));
    AlphaVelocityLabel->setText(QString::number(mEyePropertiesParameters.alphaVelocity, 'f', 2));

    // Threshold parameters

    QLabel *ThresholdParametersTextBox = new QLabel;
    ThresholdParametersTextBox->setText("<b>Change thresholds</b>");
    ThresholdParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *ThresholdCircumferenceTextBox = new QLabel;
    ThresholdCircumferenceTextBox->setText("<b>Pupil circumference:</b>");

    QLabel *ThresholdFractionTextBox = new QLabel;
    ThresholdFractionTextBox->setText("<b>Pupil fraction:</b>");

    ThresholdCircumferenceSlider = new SliderDouble;
    ThresholdCircumferenceSlider->setPrecision(1);
    ThresholdCircumferenceSlider->setDoubleRange(0, 50);
    ThresholdCircumferenceSlider->setDoubleValue(mEyePropertiesParameters.thresholdCircumferenceChangeMin);
    ThresholdCircumferenceSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdCircumferenceSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdCircumference(double)));

    ThresholdFractionSlider = new SliderDouble;
    ThresholdFractionSlider->setPrecision(2);
    ThresholdFractionSlider->setDoubleRange(0, 1.0);
    ThresholdFractionSlider->setDoubleValue(mEyePropertiesParameters.thresholdFractionChangeMin);
    ThresholdFractionSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdFractionSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdFraction(double)));

    ThresholdCircumferenceLabel = new QLabel;
    ThresholdCircumferenceLabel->setText(QString::number(mEyePropertiesParameters.thresholdCircumferenceChangeMin, 'f', 1));

    ThresholdFractionLabel = new QLabel;
    ThresholdFractionLabel->setText(QString::number(mEyePropertiesParameters.thresholdFractionChangeMin, 'f', 2));

    // Miscellaneous parameters

    QLabel *MiscParametersTextBox = new QLabel;
    MiscParametersTextBox->setText("<b>Miscellaneous</b>");
    MiscParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *PupilHaarOffsetTextBox = new QLabel;
    PupilHaarOffsetTextBox->setText("<b>Pupil Haar-offset:</b>");

    PupilHaarOffsetSlider = new SliderDouble;
    PupilHaarOffsetSlider->setPrecision(2);
    PupilHaarOffsetSlider->setDoubleRange(0, 1.0);
    PupilHaarOffsetSlider->setDoubleValue(mEyePropertiesParameters.pupilOffset);
    PupilHaarOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilHaarOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilHaarOffset(double)));

    PupilHaarOffsetLabel = new QLabel;
    PupilHaarOffsetLabel->setText(QString::number(mEyePropertiesParameters.pupilOffset, 'f', 2));

    QLabel *GlintRadiusTextBox = new QLabel;
    GlintRadiusTextBox->setText("<b>Glint radius:</b>");

    GlintRadiusSlider = new QSlider;
    GlintRadiusSlider->setRange(0, 10);
    GlintRadiusSlider->setValue(mEyePropertiesParameters.glintRadius);
    GlintRadiusSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GlintRadiusSlider, SIGNAL(valueChanged(int)), this, SLOT(setGlintRadius(int)));

    GlintRadiusLabel = new QLabel;
    GlintRadiusLabel->setText(QString::number(mEyePropertiesParameters.glintRadius));

    QLabel *CurvatureOffsetTextBox = new QLabel;
    CurvatureOffsetTextBox->setText("<b>Curvature offset:</b>");

    CurvatureOffsetSlider = new SliderDouble;
    CurvatureOffsetSlider->setPrecision(1);
    CurvatureOffsetSlider->setDoubleRange(0, 30);
    CurvatureOffsetSlider->setDoubleValue(mEyePropertiesParameters.curvatureOffset);
    CurvatureOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureOffset(double)));

    CurvatureOffsetLabel = new QLabel;
    CurvatureOffsetLabel->setText(QString::number(mEyePropertiesParameters.curvatureOffset, 'f', 1));

    QLabel *EdgeMaximumFitNumberTextBox = new QLabel;
    EdgeMaximumFitNumberTextBox->setText("<b>Edge maximum fit number:</b>");

    EdgeMaximumFitNumberSlider = new QSlider;
    EdgeMaximumFitNumberSlider->setRange(0, 7);
    EdgeMaximumFitNumberSlider->setValue(mEyePropertiesParameters.edgeMaximumFitNumber);
    EdgeMaximumFitNumberSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EdgeMaximumFitNumberSlider, SIGNAL(valueChanged(int)), this, SLOT(setEdgeMaximumFitNumber(int)));

    EdgeMaximumFitNumberLabel = new QLabel;
    EdgeMaximumFitNumberLabel->setText(QString::number(mEyePropertiesParameters.edgeMaximumFitNumber));

    QLabel *EllipseFitErrorMaximumTextBox = new QLabel;
    EllipseFitErrorMaximumTextBox->setText("<b>Ellipse maximum fit error:</b>");

    EllipseFitErrorMaximumSlider = new SliderDouble;
    EllipseFitErrorMaximumSlider->setPrecision(1);
    EllipseFitErrorMaximumSlider->setDoubleRange(0, 80);
    EllipseFitErrorMaximumSlider->setDoubleValue(mEyePropertiesParameters.ellipseFitErrorMaximum);
    EllipseFitErrorMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EllipseFitErrorMaximumSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEllipseFitErrorMaximum(double)));

    EllipseFitErrorMaximumLabel = new QLabel;
    EllipseFitErrorMaximumLabel->setText(QString::number(mEyePropertiesParameters.ellipseFitErrorMaximum, 'f', 1));

    // Experimental tabs

    NameInputLineEdit = new QLineEdit;
    NameInputLineEdit->setAlignment(Qt::AlignRight);
    NameInputLineEdit->setText(subjectIdentifier);

    QLabel *NameInputTextBox = new QLabel;
    NameInputTextBox->setText("<b>Subject identifier:</b>");

    QLabel* DataDirectoryTitleTextBox = new QLabel;
    DataDirectoryTitleTextBox->setText("<b>Data path:</b>");

    DataDirectoryTextBox = new QLabel;
    DataDirectoryTextBox->setText(QString::fromStdString(dataDirectory));

    QPushButton* DataDirectoryButton = new QPushButton("&Browse...");
    QObject::connect(DataDirectoryButton, SIGNAL(clicked()), this, SLOT(selectDirectory()));

    QLabel *DataFilenameTextBox = new QLabel;
    DataFilenameTextBox->setText("<b>File name:</b>");

    DataFilenameLineEdit = new QLineEdit;
    DataFilenameLineEdit->setAlignment(Qt::AlignRight);
    DataFilenameLineEdit->setText(QString::fromStdString(dataFilename));

    QLabel *TrialTimeLengthTextBox = new QLabel;
    TrialTimeLengthTextBox->setText("<b>Trial length (ms):</b>");

    TrialTimeLengthLineEdit = new QLineEdit;
    TrialTimeLengthLineEdit->setAlignment(Qt::AlignRight);
    TrialTimeLengthLineEdit->setText(QString::number(trialTimeLength));

    FlashStandbySlider = new QSlider;
    FlashStandbySlider->setRange(0, 1);
    FlashStandbySlider->setOrientation(Qt::Horizontal);
    QObject::connect(FlashStandbySlider, SIGNAL(valueChanged(int)), this, SLOT(onFlashStandbySlider(int)));

    QLabel* FlashStandbyTextBox = new QLabel;
    FlashStandbyTextBox->setText("<b>Flash standby mode:</b>");

    FlashStandbyLabel = new QLabel;
    FlashStandbyLabel->setText("<font color='red'><b>OFF</b></font>");

    QLabel* FlashThresholdTextBox = new QLabel;
    FlashThresholdTextBox->setText("<b>Flash threshold:</b>");

    SliderDouble* FlashThresholdSlider = new SliderDouble;
    FlashThresholdSlider->setPrecision(3);
    FlashThresholdSlider->setDoubleRange(0.80, 1.00);
    FlashThresholdSlider->setDoubleValue(flashThreshold);
    FlashThresholdSlider->setOrientation(Qt::Horizontal);

    FlashThresholdLabel = new QLabel;
    FlashThresholdLabel->setText(QString::number(flashThreshold, 'f', 3));

    QObject::connect(FlashThresholdSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setFlashThreshold(double)));

    QSpinBox* FlashXPosSpinBox = new QSpinBox;
    QSpinBox* FlashYPosSpinBox = new QSpinBox;
    QSpinBox* FlashWdthSpinBox = new QSpinBox;
    QSpinBox* FlashHghtSpinBox = new QSpinBox;

    FlashXPosSpinBox->setRange(0, Parameters::cameraXResolution);
    FlashYPosSpinBox->setRange(0, Parameters::cameraYResolution);
    FlashWdthSpinBox->setRange(0, Parameters::cameraXResolution);
    FlashHghtSpinBox->setRange(0, Parameters::cameraYResolution);

    FlashXPosSpinBox->setValue(Parameters::flashAOIXPos);
    FlashYPosSpinBox->setValue(Parameters::flashAOIYPos);
    FlashWdthSpinBox->setValue(Parameters::flashAOIWdth);
    FlashHghtSpinBox->setValue(Parameters::flashAOIHght);

    QObject::connect(FlashXPosSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setFlashAOIXPos(int)));
    QObject::connect(FlashYPosSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setFlashAOIYPos(int)));
    QObject::connect(FlashWdthSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setFlashAOIWdth(int)));
    QObject::connect(FlashHghtSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setFlashAOIHght(int)));

    QLabel* FlashCoordsTextBox = new QLabel;
    FlashCoordsTextBox->setText("<b>Flash AOI</b>");

    QLabel* FlashXPosTextBox = new QLabel;
    FlashXPosTextBox->setText("<i>X:</i>");

    QLabel* FlashYPosTextBox = new QLabel;
    FlashYPosTextBox->setText("<i>Y:</i>");

    QLabel* FlashWdthTextBox = new QLabel;
    FlashWdthTextBox->setText("<i>W:</i>");

    QLabel* FlashHghtTextBox = new QLabel;
    FlashHghtTextBox->setText("<i>H:</i>");

    QHBoxLayout* FlashCoordinatesLayout = new QHBoxLayout;
    FlashCoordinatesLayout->addWidget(FlashCoordsTextBox);
    FlashCoordinatesLayout->addWidget(FlashXPosTextBox);
    FlashCoordinatesLayout->addWidget(FlashXPosSpinBox);
    FlashCoordinatesLayout->addWidget(FlashYPosTextBox);
    FlashCoordinatesLayout->addWidget(FlashYPosSpinBox);
    FlashCoordinatesLayout->addWidget(FlashWdthTextBox);
    FlashCoordinatesLayout->addWidget(FlashWdthSpinBox);
    FlashCoordinatesLayout->addWidget(FlashHghtTextBox);
    FlashCoordinatesLayout->addWidget(FlashHghtSpinBox);

    QLabel* TrialIndexTextBox = new QLabel;
    TrialIndexTextBox->setText("<b>Trial index:</b>");

    TrialIndexSpinBox = new QSpinBox;
    TrialIndexSpinBox->setRange(0, 999);
    TrialIndexSpinBox->setValue(trialIndex);
    TrialIndexSpinBox->setAlignment(Qt::AlignRight);
    QObject::connect(TrialIndexSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setTrialIndex(int)));

    QHBoxLayout* TrialIndexSpinBoxLayout = new QHBoxLayout;
    TrialIndexSpinBoxLayout->addStretch();
    TrialIndexSpinBoxLayout->addWidget(TrialIndexSpinBox);

    QWidget* ExperimentTabWidget = new QWidget;
    QGridLayout *ExperimentTabLayout = new QGridLayout(ExperimentTabWidget);
    ExperimentTabLayout->addWidget(NameInputTextBox, 0, 0);
    ExperimentTabLayout->addWidget(NameInputLineEdit, 0, 1);
    ExperimentTabLayout->addWidget(DataDirectoryTitleTextBox, 1, 0);
    ExperimentTabLayout->addWidget(DataDirectoryTextBox, 1, 1);
    ExperimentTabLayout->addWidget(DataDirectoryButton, 1, 2);
    ExperimentTabLayout->addWidget(DataFilenameTextBox, 2, 0);
    ExperimentTabLayout->addWidget(DataFilenameLineEdit, 2, 1);
    ExperimentTabLayout->addWidget(TrialIndexTextBox, 3, 0);
    ExperimentTabLayout->addLayout(TrialIndexSpinBoxLayout, 3, 1);
    ExperimentTabLayout->addWidget(TrialTimeLengthTextBox, 4, 0);
    ExperimentTabLayout->addWidget(TrialTimeLengthLineEdit, 4, 1);
    ExperimentTabLayout->addWidget(FlashStandbyTextBox, 5, 0);
    ExperimentTabLayout->addWidget(FlashStandbySlider, 5, 1);
    ExperimentTabLayout->addWidget(FlashStandbyLabel, 5, 2);
    ExperimentTabLayout->addWidget(FlashThresholdTextBox, 6, 0);
    ExperimentTabLayout->addWidget(FlashThresholdSlider, 6, 1);
    ExperimentTabLayout->addWidget(FlashThresholdLabel, 6, 2);
    ExperimentTabLayout->addLayout(FlashCoordinatesLayout, 7, 0, 1, 3);

    ExperimentTabLayout->setColumnStretch(0, 1);
    ExperimentTabLayout->setColumnStretch(1, 2);
    ExperimentTabLayout->setColumnStretch(2, 1);
    ExperimentTabLayout->setColumnStretch(3, 5);

    // Set-up parameter settings layout

    QWidget* RealTimeVariablesWidget = new QWidget;
    QGridLayout *RealTimeVariablesLayout = new QGridLayout(RealTimeVariablesWidget);
    RealTimeVariablesLayout->addWidget(PupilCircfTextBox, 0, 0);
    RealTimeVariablesLayout->addWidget(PupilCircfSliderDouble, 0, 1);
    RealTimeVariablesLayout->addWidget(PupilCircfLabel, 0, 2);
    RealTimeVariablesLayout->addWidget(PupilFractTextBox, 1, 0);
    RealTimeVariablesLayout->addWidget(PupilFractSliderDouble, 1, 1);
    RealTimeVariablesLayout->addWidget(PupilFractLabel, 1 ,2);
    RealTimeVariablesLayout->addWidget(EdgeIntensityTextBox, 2, 0);
    RealTimeVariablesLayout->addWidget(EdgeIntensitySliderDouble, 2, 1);
    RealTimeVariablesLayout->addWidget(EdgeIntensityLabel, 2, 2);

    QWidget *ParameterLimitsWidget = new QWidget;
    QGridLayout *ParameterLimitsLayout = new QGridLayout(ParameterLimitsWidget);
    ParameterLimitsLayout->addWidget(PupilCircumferenceIniTextBox, 0, 0);
    ParameterLimitsLayout->addWidget(PupilCircumferenceIniSlider, 0, 1);
    ParameterLimitsLayout->addWidget(PupilCircumferenceIniLabel, 0, 2);
    ParameterLimitsLayout->addWidget(PupilCircumferenceOffTextBox, 1, 0);
    ParameterLimitsLayout->addWidget(PupilCircumferenceOffSlider, 1, 1);
    ParameterLimitsLayout->addWidget(PupilCircumferenceOffLabel, 1, 2);
    ParameterLimitsLayout->addWidget(PupilFractionIniTextBox, 2, 0);
    ParameterLimitsLayout->addWidget(PupilFractionIniSlider, 2, 1);
    ParameterLimitsLayout->addWidget(PupilFractionIniLabel, 2, 2);
    ParameterLimitsLayout->addWidget(PupilFractionMinTextBox, 3, 0);
    ParameterLimitsLayout->addWidget(PupilFractionMinSlider, 3, 1);
    ParameterLimitsLayout->addWidget(PupilFractionMinLabel, 3, 2);
    ParameterLimitsLayout->addWidget(EdgeIntensityOffsetTextBox, 4, 0);
    ParameterLimitsLayout->addWidget(EdgeIntensityOffsetSlider, 4, 1);
    ParameterLimitsLayout->addWidget(EdgeIntensityOffsetLabel, 4, 2);

    QWidget *LearningRateWidget = new QWidget;
    QGridLayout *LearningRateLayout = new QGridLayout(LearningRateWidget);
    LearningRateLayout->addWidget(AlphaPupilTextBox, 0, 0);
    LearningRateLayout->addWidget(AlphaPupilSlider, 0, 1);
    LearningRateLayout->addWidget(AlphaPupilLabel, 0, 2);
    LearningRateLayout->addWidget(AlphaMomentumTextBox, 1, 0);
    LearningRateLayout->addWidget(AlphaMomentumSlider, 1, 1);
    LearningRateLayout->addWidget(AlphaMomentumLabel, 1, 2);
    LearningRateLayout->addWidget(AlphaSearchAreaTextBox, 2, 0);
    LearningRateLayout->addWidget(AlphaSearchAreaSlider, 2, 1);
    LearningRateLayout->addWidget(AlphaSearchAreaLabel, 2, 2);
    LearningRateLayout->addWidget(AlphaGeneralTextBox, 3, 0);
    LearningRateLayout->addWidget(AlphaGeneralSlider, 3, 1);
    LearningRateLayout->addWidget(AlphaGeneralLabel, 3, 2);
    LearningRateLayout->addWidget(AlphaVelocityTextBox, 4, 0);
    LearningRateLayout->addWidget(AlphaVelocitySlider, 4, 1);
    LearningRateLayout->addWidget(AlphaVelocityLabel, 4, 2);

    QWidget *ThresholdParametersWidget = new QWidget;
    QGridLayout *ThresholdParametersLayout = new QGridLayout(ThresholdParametersWidget);
    ThresholdParametersLayout->addWidget(ThresholdCircumferenceTextBox, 0, 0);
    ThresholdParametersLayout->addWidget(ThresholdCircumferenceSlider, 0, 1);
    ThresholdParametersLayout->addWidget(ThresholdCircumferenceLabel, 0, 2);
    ThresholdParametersLayout->addWidget(ThresholdFractionTextBox, 1, 0);
    ThresholdParametersLayout->addWidget(ThresholdFractionSlider, 1, 1);
    ThresholdParametersLayout->addWidget(ThresholdFractionLabel, 1, 2);

    QWidget *CannyEdgeWidget = new QWidget;
    QGridLayout *CannyEdgeLayout = new QGridLayout(CannyEdgeWidget);
    CannyEdgeLayout->addWidget(CannyUpperLimitTextBox, 0, 0);
    CannyEdgeLayout->addWidget(CannyUpperLimitSlider, 0, 1);
    CannyEdgeLayout->addWidget(CannyUpperLimitLabel, 0, 2);
    CannyEdgeLayout->addWidget(CannyLowerLimitTextBox, 1, 0);
    CannyEdgeLayout->addWidget(CannyLowerLimitSlider, 1, 1);
    CannyEdgeLayout->addWidget(CannyLowerLimitLabel, 1, 2);
    CannyEdgeLayout->addWidget(CannyKernelSizeTextBox, 2, 0);
    CannyEdgeLayout->addWidget(CannyKernelSizeSlider, 2, 1);
    CannyEdgeLayout->addWidget(CannyKernelSizeLabel, 2, 2);
    CannyEdgeLayout->addWidget(CannyBlurLevelTextBox, 3, 0);
    CannyEdgeLayout->addWidget(CannyBlurLevelSlider, 3, 1);
    CannyEdgeLayout->addWidget(CannyBlurLevelLabel, 3, 2);

    QWidget *MiscellaneousParametersWidget = new QWidget;
    QGridLayout *MiscellaneousParametersLayout = new QGridLayout(MiscellaneousParametersWidget);
    MiscellaneousParametersLayout->addWidget(PupilHaarOffsetTextBox, 0, 0);
    MiscellaneousParametersLayout->addWidget(PupilHaarOffsetSlider, 0, 1);
    MiscellaneousParametersLayout->addWidget(PupilHaarOffsetLabel, 0, 2);
    MiscellaneousParametersLayout->addWidget(GlintRadiusTextBox, 1, 0);
    MiscellaneousParametersLayout->addWidget(GlintRadiusSlider, 1, 1);
    MiscellaneousParametersLayout->addWidget(GlintRadiusLabel, 1, 2);
    MiscellaneousParametersLayout->addWidget(CurvatureOffsetTextBox, 2, 0);
    MiscellaneousParametersLayout->addWidget(CurvatureOffsetSlider, 2, 1);
    MiscellaneousParametersLayout->addWidget(CurvatureOffsetLabel, 2, 2);
    MiscellaneousParametersLayout->addWidget(EdgeMaximumFitNumberTextBox, 3, 0);
    MiscellaneousParametersLayout->addWidget(EdgeMaximumFitNumberSlider, 3, 1);
    MiscellaneousParametersLayout->addWidget(EdgeMaximumFitNumberLabel, 3, 2);
    MiscellaneousParametersLayout->addWidget(EllipseFitErrorMaximumTextBox, 4, 0);
    MiscellaneousParametersLayout->addWidget(EllipseFitErrorMaximumSlider, 4, 1);
    MiscellaneousParametersLayout->addWidget(EllipseFitErrorMaximumLabel, 4, 2);

    EyeTrackingParameterTabWidget = new QTabWidget;
    EyeTrackingParameterTabWidget->addTab(CameraParametersWidget,           tr("Camera"));
    EyeTrackingParameterTabWidget->addTab(RealTimeVariablesWidget,          tr("Variables"));
    EyeTrackingParameterTabWidget->addTab(ParameterLimitsWidget,            tr("Limits"));
    EyeTrackingParameterTabWidget->addTab(CannyEdgeWidget,                  tr("Canny-edge"));
    EyeTrackingParameterTabWidget->addTab(LearningRateWidget,               tr("Learning rates"));
    EyeTrackingParameterTabWidget->addTab(ThresholdParametersWidget,        tr("Thresholds"));
    EyeTrackingParameterTabWidget->addTab(MiscellaneousParametersWidget,    tr("Miscellaneous"));
    EyeTrackingParameterTabWidget->addTab(ExperimentTabWidget,              tr("Experimental"));

    QWidget *EyeTrackingWidget = new QWidget;
    QVBoxLayout *EyeTrackingLayout = new QVBoxLayout(EyeTrackingWidget);
    EyeTrackingLayout->addWidget(CameraSettings);
    EyeTrackingLayout->addWidget(EyeTrackingParameterTabWidget);

    // Tab widget

    QSize screenResolution = QDesktopWidget().availableGeometry(this).size();
    double screenOffset = 0.9;

    EyeTrackingWidget->setMaximumSize(screenResolution * screenOffset);
    EyeTrackingWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    QScrollArea *EyeTrackingWidgetScrollArea = new QScrollArea();
    EyeTrackingWidgetScrollArea->setWidget(EyeTrackingWidget);
    EyeTrackingWidgetScrollArea->setWidgetResizable(true);

    // External widgets

    // Quit button

    QPushButton *QuitButton = new QPushButton("Quit");
    QObject::connect(QuitButton, SIGNAL(clicked()), this, SLOT(onQuitButtonClicked()));

    QGridLayout *ExternalLayout = new QGridLayout;
    ExternalLayout->addWidget(QuitButton, 0, 1);
    ExternalLayout->setColumnStretch(0, 1);
    ExternalLayout->setColumnStretch(1, 0);

    // Main layout

    QVBoxLayout *MainLayout = new QVBoxLayout;
    MainLayout->addWidget(EyeTrackingWidgetScrollArea);
    MainLayout->addLayout(ExternalLayout);

    // Central widget set-up

    QWidget* centralWidget = new QWidget();
    centralWidget->setMaximumSize(screenResolution);
    centralWidget->setLayout(MainLayout);
    setCentralWidget(centralWidget);
    resize(0.7 * screenResolution);

    // Start camera

    std::thread findCameraThread(&MainWindow::findCamera, this);
    findCameraThread.detach();

    UpdateCameraImageTimer = new QTimer;
    QObject::connect(UpdateCameraImageTimer, SIGNAL(timeout()), this, SLOT(updateCameraImage()));
    QObject::connect(this, SIGNAL(startTimer(int)), UpdateCameraImageTimer, SLOT(start(int)));
    QObject::connect(this, SIGNAL(stopTimer()), UpdateCameraImageTimer, SLOT(stop()));
    emit startTimer(round(1000 / guiUpdateFrequency));
}

MainWindow::~MainWindow()
{

}

void MainWindow::pupilTracking()
{
    while(APP_RUNNING && Parameters::CAMERA_RUNNING && Parameters::REALTIME_PROCESSING)
    {
        eyeProperties mEyePropertiesTemp;

        int eyeAOIXPosTemp;
        int eyeAOIYPosTemp;
        int eyeAOIWdthTemp;
        int eyeAOIHghtTemp;

        int flashAOIXPosTemp;
        int flashAOIYPosTemp;
        int flashAOIWdthTemp;
        int flashAOIHghtTemp;

        imageInfo mImageInfo = mUEyeOpencvCam.getFrame(); // get new frame from camera

        cv::Mat imageOriginal = mImageInfo.image;
        absoluteTime = mImageInfo.time; // Get frame timestamp

        double relativeTimeNew = (absoluteTime - startTime) / (double) 10000; // in ms

        if (relativeTimeNew <= (relativeTime + 0.9 * (1000 / cameraFrameRate)))
        {
            continue; // ignore frame if time interval was too short (possible camera error)
        }

        relativeTime = relativeTimeNew;

        std::lock_guard<std::mutex> secondaryMutexLock(Parameters::secondaryMutex);

        { std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);

            imageCamera = imageOriginal.clone();

            mEyePropertiesTemp.v = mEyePropertiesVariables;
            mEyePropertiesTemp.p = mEyePropertiesParameters;

            eyeAOIXPosTemp = Parameters::eyeAOIXPos;
            eyeAOIYPosTemp = Parameters::eyeAOIYPos;
            eyeAOIWdthTemp = Parameters::eyeAOIWdth;
            eyeAOIHghtTemp = Parameters::eyeAOIHght;

            if (FLASH_STANDBY)
            {
                flashAOIXPosTemp = Parameters::flashAOIXPos - Parameters::cameraAOIXPos;
                flashAOIYPosTemp = Parameters::flashAOIYPos - Parameters::cameraAOIYPos;
                flashAOIWdthTemp = Parameters::flashAOIWdth;
                flashAOIHghtTemp = Parameters::flashAOIHght;

                if (flashAOIXPosTemp < 0)
                {
                    if (flashAOIXPosTemp + flashAOIWdthTemp < 0) // Flash AOI not in camera AOI
                    {
                        FlashStandbySlider->setValue(0);
                        continue;
                    }
                    else
                    {
                        flashAOIWdthTemp = flashAOIXPosTemp + flashAOIWdthTemp;
                        flashAOIXPosTemp = 0;
                    }
                }
                else if (flashAOIXPosTemp + flashAOIWdthTemp >= Parameters::cameraAOIWdth)
                {
                    flashAOIWdthTemp = Parameters::cameraAOIWdth - flashAOIXPosTemp;

                    if (flashAOIWdthTemp <= 0)
                    {
                        FlashStandbySlider->setValue(0);
                        continue;
                    }
                }

                if (flashAOIYPosTemp < 0)
                {
                    if (flashAOIYPosTemp + flashAOIHghtTemp < 0) // Flash AOI not in camera AOI
                    {
                        FlashStandbySlider->setValue(0);
                        continue;
                    }
                    else
                    {
                        flashAOIHghtTemp = flashAOIYPosTemp + flashAOIHghtTemp;
                        flashAOIYPosTemp = 0;
                    }
                }
                else if (flashAOIYPosTemp + flashAOIHghtTemp >= Parameters::cameraAOIHght)
                {
                    flashAOIHghtTemp = Parameters::cameraAOIHght - flashAOIYPosTemp;

                    if (flashAOIHghtTemp <= 0)
                    {
                        FlashStandbySlider->setValue(0);
                        continue;
                    }
                }
            }
        }

        // Check limits

        int imgWdth = imageOriginal.cols;
        int imgHght = imageOriginal.rows;

        if (imgHght < (eyeAOIYPosTemp + eyeAOIHghtTemp) || imgWdth < (eyeAOIXPosTemp + eyeAOIWdthTemp))
        {
            continue;
        }

        if (!EXPERIMENT_TRIAL_RECORDING)
        {
            if (FLASH_STANDBY)
            {
                // reset variables

                mEyePropertiesTemp.v.pupilFraction = mEyePropertiesTemp.v.pupilFraction + mEyePropertiesTemp.p.alphaPupil * (mEyePropertiesTemp.p.pupilFractIni - mEyePropertiesTemp.v.pupilFraction);
                mEyePropertiesTemp.v.pupilCircumference = mEyePropertiesTemp.v.pupilCircumference + mEyePropertiesTemp.p.alphaPupil * (mEyePropertiesTemp.p.pupilCircfIni - mEyePropertiesTemp.v.pupilCircumference);

                mEyePropertiesTemp.v.momentumFraction = mEyePropertiesTemp.v.momentumFraction * mEyePropertiesTemp.p.alphaMomentum;
                mEyePropertiesTemp.v.momentumCircumference = mEyePropertiesTemp.v.momentumCircumference * mEyePropertiesTemp.p.alphaMomentum;
                mEyePropertiesTemp.v.momentumRadius = mEyePropertiesTemp.v.momentumRadius * mEyePropertiesTemp.p.alphaMomentum;

                mEyePropertiesTemp.v.edgeIntensity = mEyePropertiesTemp.v.edgeIntensity + mEyePropertiesTemp.p.alphaPupil * (mEyePropertiesTemp.p.edgeIntensityIni - mEyePropertiesTemp.v.edgeIntensity);

                mEyePropertiesTemp.v.xVelocity = mEyePropertiesTemp.v.xVelocity * mEyePropertiesTemp.p.alphaVelocity;
                mEyePropertiesTemp.v.yVelocity = mEyePropertiesTemp.v.yVelocity * mEyePropertiesTemp.p.alphaVelocity;

                mEyePropertiesTemp.v.searchRadius = mEyePropertiesTemp.v.searchRadius * (2 - mEyePropertiesTemp.p.alphaGeneral);
                mEyePropertiesTemp.v.thresholdCircumferenceChange = mEyePropertiesTemp.v.thresholdCircumferenceChange * (2 - mEyePropertiesTemp.p.alphaGeneral);
                mEyePropertiesTemp.v.thresholdFractionChange = mEyePropertiesTemp.v.thresholdFractionChange * (2 - mEyePropertiesTemp.p.alphaGeneral);


                cv::Rect flashRegion(flashAOIXPosTemp, flashAOIYPosTemp, flashAOIWdthTemp, flashAOIHghtTemp);

                if (flashDetection(imageOriginal(flashRegion), flashThreshold))
                {
                    startTime = mImageInfo.time;
                    startTrialRecording();
                }
            }
            else
            {
                cv::Rect eyeRegion(eyeAOIXPosTemp, eyeAOIYPosTemp, eyeAOIWdthTemp, eyeAOIHghtTemp);
                mEyePropertiesTemp = pupilDetector(imageOriginal(eyeRegion), mEyePropertiesTemp); // Pupil tracking algorithm
            }
        }
        else
        {
            if (!SAVE_EYE_IMAGE)
            {
                cv::Rect eyeRegion(eyeAOIXPosTemp, eyeAOIYPosTemp, eyeAOIWdthTemp, eyeAOIHghtTemp);
                mEyePropertiesTemp = pupilDetector(imageOriginal(eyeRegion), mEyePropertiesTemp); // Pupil tracking algorithm

                timeStamps[frameCount] = relativeTime; // save time stamps

                mEyePropertiesTemp.v.xPosAbs = mEyePropertiesTemp.v.xPos + eyeAOIXPosTemp;
                mEyePropertiesTemp.v.yPosAbs = mEyePropertiesTemp.v.yPos + eyeAOIYPosTemp;

                eyeXPositions[frameCount] = mEyePropertiesTemp.v.xPosAbs;
                eyeYPositions[frameCount] = mEyePropertiesTemp.v.yPosAbs;

                eyeDetectionFlags[frameCount] = mEyePropertiesTemp.v.pupilDetected;

                frameCount++;
            }
            else
            {
                timeStamps[frameCount] = relativeTime; // save time stamps
                eyeDetectionFlags[frameCount] = true;

                // Saving screens
                std::stringstream filename;
                filename << dataDirectory << "/"
                         << currentDate << "/"
                         << (NameInputLineEdit->text()).toStdString()
                         << "/trial_" << trialIndex
                         << "/raw/" << frameCount << ".png";

                cv::imwrite(filename.str(), imageOriginal);

                frameCount++;
            }

            if (frameCount >= trialFrameTotal)
            {
                mUEyeOpencvCam.stopRecording();
                saveTrialData();
                trialIndex++;
                TrialIndexSpinBox->setValue(trialIndex);
                EXPERIMENT_TRIAL_RECORDING = false;
                emit startTimer(round(1000 / guiUpdateFrequency));
            }
        }

        // Update structures

        {
            std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);

            mEyePropertiesVariables = mEyePropertiesTemp.v;
            mEyePropertiesMiscellaneous = mEyePropertiesTemp.m;
        }
    }

    if (!APP_RUNNING)
    {
        std::unique_lock<std::mutex> mtxLock(mtx);
        APP_EXIT = true;
        cv.notify_all();
    }
    else if (!Parameters::REALTIME_PROCESSING)
    {
        std::unique_lock<std::mutex> mtxLock(mtx);
        Parameters::CAMERA_RUNNING = false;
        cv.notify_all();
    }
    else if (!Parameters::CAMERA_RUNNING)
    {
        std::thread findCameraThread(&MainWindow::findCamera, this);
        findCameraThread.detach();
    }
}


void MainWindow::updateCameraImage()
{
    if (Parameters::REALTIME_PROCESSING)
    {
        if (!EXPERIMENT_TRIAL_RECORDING)
        {
            if (Parameters::CAMERA_RUNNING)
            {
                eyeProperties mEyePropertiesTemp;

                cv::Mat imageOriginal;

                int eyeAOIXPos;
                int eyeAOIYPos;
                int eyeAOIWdth;
                int eyeAOIHght;

                int flashAOIXPos;
                int flashAOIYPos;
                int flashAOIWdth;
                int flashAOIHght;

                { std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);

                    if (!imageCamera.empty())
                    {
                        imageOriginal = imageCamera.clone();

                        mEyePropertiesTemp.v = mEyePropertiesVariables;
                        mEyePropertiesTemp.m = mEyePropertiesMiscellaneous;

                        eyeAOIXPos = Parameters::eyeAOIXPos;
                        eyeAOIYPos = Parameters::eyeAOIYPos;
                        eyeAOIWdth = Parameters::eyeAOIWdth;
                        eyeAOIHght = Parameters::eyeAOIHght;

                        flashAOIXPos = Parameters::flashAOIXPos;
                        flashAOIYPos = Parameters::flashAOIYPos;
                        flashAOIWdth = Parameters::flashAOIWdth;
                        flashAOIHght = Parameters::flashAOIHght;
                    }
                    else
                    {
                        return;
                    }
                }

                CamQImage->loadImage(imageOriginal);
                CamQImage->setEyeAOI(eyeAOIXPos, eyeAOIYPos, eyeAOIWdth, eyeAOIHght);
                CamQImage->setFlashAOI(flashAOIXPos, flashAOIYPos, flashAOIWdth, flashAOIHght);

                CamQImage->setImage();

                if (eyeAOIWdth >= eyeAOIWdthMin && eyeAOIHght >= eyeAOIHghtMin)
                {
                    if (!mEyePropertiesTemp.m.errorDetected)
                    {
                        cv::Mat imageEye = mEyePropertiesTemp.m.image.clone();
                        drawAll(imageEye, mEyePropertiesTemp);
                        EyeQImage->loadImage(imageEye);
                        EyeQImage->setImage();
                    }
                }
                else
                {
                    EyeQImage->setAOIError();
                }

                // update sliders

                PupilCircfSliderDouble->setDoubleValue(mEyePropertiesTemp.v.pupilCircumference);
                PupilCircfLabel->setText(QString::number(mEyePropertiesTemp.v.pupilCircumference, 'f', 1));

                PupilFractSliderDouble->setDoubleValue(mEyePropertiesTemp.v.pupilFraction);
                PupilFractLabel->setText(QString::number(mEyePropertiesTemp.v.pupilFraction, 'f', 2));

                EdgeIntensitySliderDouble->setDoubleValue(mEyePropertiesTemp.v.edgeIntensity);
                EdgeIntensityLabel->setText(QString::number(mEyePropertiesTemp.v.edgeIntensity, 'f', 1));

                if (CameraHardwareGainAutoCheckBox->checkState())
                {
                    int hardwareGain = mUEyeOpencvCam.getHardwareGain();
                    CameraHardwareGainSlider->setValue(hardwareGain);
                    CameraHardwareGainLabel->setText(QString::number(hardwareGain));
                }

                // increase pixel clock if desired frame-rate has not been reached
                if ((cameraPixelClock < CameraPixelClockSlider->maximum()) && (CameraFrameRateDesiredSpinBox->value() > cameraFrameRate))
                {
                    cameraPixelClock = cameraPixelClock + 1;
                    CameraPixelClockSlider->setValue(cameraPixelClock);
                }

            }
            else if (!Parameters::CAMERA_RUNNING)
            {
                CamQImage->setFindingCamera();
                EyeQImage->setSpinner();
            }
        }
    }
    else
    {
        emit stopTimer();
    }
}

void MainWindow::getCameraParameters()
{
    std::vector<int> pixelClockRange = mUEyeOpencvCam.getPixelClockRange();

    CameraPixelClockSlider->setRange(pixelClockRange[0], pixelClockRange[1]);
    CameraPixelClockSlider->setValue(cameraPixelClock);
    CameraPixelClockLabel->setText(QString::number(cameraPixelClock));

    std::vector<double> frameRateRange = mUEyeOpencvCam.getFrameRateRange();

    cameraFrameRate = frameRateRange[1]; // set to max
    CameraFrameRateSlider->setDoubleRange(frameRateRange[0], cameraFrameRate);
    CameraFrameRateSlider->setDoubleValue(cameraFrameRate);
    CameraFrameRateLabel->setText(QString::number(cameraFrameRate, 'f', 1));

    std::vector<int> blackLevelOffsetRange = mUEyeOpencvCam.getBlackLevelOffsetRange();

    CameraBlackLevelOffsetSlider->setRange(blackLevelOffsetRange[0], blackLevelOffsetRange[1]);
    CameraBlackLevelOffsetLabel->setText(QString::number(blackLevelOffsetRange[1]));
}

void MainWindow::findCamera()
{
    bool CAMERA_START = false;

    while (APP_RUNNING && !Parameters::CAMERA_RUNNING && Parameters::REALTIME_PROCESSING)
    {
        if (mUEyeOpencvCam.findCamera())
        {
            int retInt = mUEyeOpencvCam.initCamera();

            if (retInt == 0)
            {
                continue;
            }
            else if (retInt == 1)
            {
                if (!mUEyeOpencvCam.setColorMode())
                {
                    continue;
                }
            }
            else if (retInt == 2)
            {
                mUEyeOpencvCam.exitCamera();

                if (!mUEyeOpencvCam.initCamera())
                {
                    continue;
                }
            }

            if (mUEyeOpencvCam.setSubSampling(cameraSubSamplingFactor))
            {
                if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOIWdth, Parameters::cameraAOIHght))
                {
                    if (mUEyeOpencvCam.setAOI(Parameters::cameraAOIXPos, Parameters::cameraAOIYPos, Parameters::cameraAOIWdth, Parameters::cameraAOIHght))
                    {
                        CAMERA_START = true;
                        break;
                    }
                }
            }
        }
    }

    if (CAMERA_START)
    {
        if (mUEyeOpencvCam.startVideoCapture())
        {
            mUEyeOpencvCam.setAutoGain(true);

            Parameters::CAMERA_RUNNING = true;
            Parameters::CAMERA_READY = true;

            std::thread pupilTrackingThread(&MainWindow::pupilTracking, this);
            pupilTrackingThread.detach();

            getCameraParameters();
            emit startTimer(round(1000 / guiUpdateFrequency));
        }
    }
    else
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }
}

// General functions

void MainWindow::msWait(int ms)
{
    QTime dieTime = QTime::currentTime().addMSecs(ms);
    while( QTime::currentTime() < dieTime)
    {
        QCoreApplication::processEvents(QEventLoop::AllEvents);
    }
}

// General slots

void MainWindow::onQuitButtonClicked()
{
    if (Parameters::CAMERA_RUNNING)
    {
        APP_RUNNING = false;
        std::unique_lock<std::mutex> lck(mtx);
        while (!APP_EXIT) cv.wait(lck);
    }

    mUEyeOpencvCam.exitCamera();

    saveSettings(LastUsedSettingsFileName);

    qApp->quit();
}

void MainWindow::setReviewMode(int state)
{
    if (!state)
    {
        editImageIndex = 0;

        EyeTrackingParameterTabWidget->setUpdatesEnabled(false);
        EyeTrackingParameterTabWidget->insertTab(0, CameraParametersWidget, tr("Camera"));
        EyeTrackingParameterTabWidget->setUpdatesEnabled(true);

        CamEyeAOIWdthSlider->setVisible(true);
        CamEyeAOIHghtSlider->setVisible(true);
        CamEyeAOIXPosSlider->setVisible(true);
        CamEyeAOIYPosSlider->setVisible(true);

        EyeTrackingReviewWidget->setVisible(false);

        Parameters::REALTIME_PROCESSING = true;
        Parameters::CAMERA_RUNNING = true;
        Parameters::CAMERA_READY = true;

        if (mUEyeOpencvCam.startVideoCapture())
        {
            std::thread pupilTrackingThread(&MainWindow::pupilTracking, this);
            pupilTrackingThread.detach();

            getCameraParameters();
            emit startTimer(round(1000 / guiUpdateFrequency));
        }
        else
        {
            std::thread findCameraThread(&MainWindow::findCamera, this);
            findCameraThread.detach();
        }
    }
    else
    {
        startReviewSession();
    }
}

void MainWindow::setParameterWidgets()
{
    EyeWdthROISlider->setDoubleValue(eyeAOIWdthFraction);
    EyeHghtROISlider->setDoubleValue(eyeAOIHghtFraction);

    PupilCircumferenceIniSlider->setDoubleValue(mEyePropertiesParameters.pupilCircfIni);
    PupilCircumferenceIniLabel->setText(QString::number(mEyePropertiesParameters.pupilCircfIni, 'f', 1));

    circumferenceOffset = mEyePropertiesParameters.pupilCircumferenceMax - mEyePropertiesParameters.pupilCircfIni;

    PupilCircumferenceOffSlider->setDoubleRange(0.1, mEyePropertiesParameters.pupilCircfIni);
    PupilCircumferenceOffSlider->setDoubleValue(circumferenceOffset);
    PupilCircumferenceOffLabel->setText(QString::number(circumferenceOffset, 'f', 1));

    PupilFractionIniSlider->setDoubleRange(mEyePropertiesParameters.pupilFractMin, 1.00);
    PupilFractionIniSlider->setDoubleValue(mEyePropertiesParameters.pupilFractIni);
    PupilFractionIniLabel->setText(QString::number(mEyePropertiesParameters.pupilFractIni, 'f', 2));

    PupilFractionMinSlider->setDoubleRange(0.0, mEyePropertiesParameters.pupilFractIni);
    PupilFractionMinSlider->setDoubleValue(mEyePropertiesParameters.pupilFractMin);
    PupilFractionMinLabel->setText(QString::number(mEyePropertiesParameters.pupilFractMin, 'f', 2));

    EdgeIntensityOffsetSlider->setDoubleValue(mEyePropertiesParameters.edgeIntensityOffset);
    EdgeIntensityOffsetLabel->setText(QString::number(mEyePropertiesParameters.edgeIntensityOffset, 'f', 1));

    CannyLowerLimitSlider->setRange(1, mEyePropertiesParameters.cannyUpperLimit);
    CannyLowerLimitSlider->setValue(mEyePropertiesParameters.cannyLowerLimit);
    CannyLowerLimitLabel->setText(QString::number(mEyePropertiesParameters.cannyLowerLimit));

    CannyUpperLimitSlider->setRange(mEyePropertiesParameters.cannyLowerLimit, 4 * mEyePropertiesParameters.cannyLowerLimit);
    CannyUpperLimitSlider->setValue(mEyePropertiesParameters.cannyUpperLimit);
    CannyUpperLimitLabel->setText(QString::number(mEyePropertiesParameters.cannyUpperLimit));

    CannyBlurLevelSlider->setValue(ceil(0.5 * mEyePropertiesParameters.cannyBlurLevel));
    CannyBlurLevelLabel->setText(QString::number(ceil(0.5 * mEyePropertiesParameters.cannyBlurLevel)));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mEyePropertiesParameters.cannyKernelSize));
    CannyKernelSizeLabel->setText(QString::number(mEyePropertiesParameters.cannyKernelSize));

    AlphaPupilSlider->setDoubleValue(mEyePropertiesParameters.alphaPupil);
    AlphaSearchAreaSlider->setDoubleValue(mEyePropertiesParameters.alphaPosition);
    AlphaGeneralSlider->setDoubleValue(mEyePropertiesParameters.alphaGeneral);
    AlphaMomentumSlider->setDoubleValue(mEyePropertiesParameters.alphaMomentum);
    AlphaVelocitySlider->setDoubleValue(mEyePropertiesParameters.alphaVelocity);

    AlphaPupilLabel->setText(QString::number(mEyePropertiesParameters.alphaPupil, 'f', 2));
    AlphaSearchAreaLabel->setText(QString::number(mEyePropertiesParameters.alphaPosition, 'f', 2));
    AlphaGeneralLabel->setText(QString::number(mEyePropertiesParameters.alphaGeneral, 'f', 2));
    AlphaMomentumLabel->setText(QString::number(mEyePropertiesParameters.alphaMomentum, 'f', 2));
    AlphaVelocityLabel->setText(QString::number(mEyePropertiesParameters.alphaVelocity, 'f', 2));

    ThresholdCircumferenceSlider->setDoubleValue(mEyePropertiesParameters.thresholdCircumferenceChangeMin);
    ThresholdCircumferenceLabel->setText(QString::number(mEyePropertiesParameters.thresholdCircumferenceChangeMin, 'f', 1));

    ThresholdFractionSlider->setDoubleValue(mEyePropertiesParameters.thresholdFractionChangeMin);
    ThresholdFractionLabel->setText(QString::number(mEyePropertiesParameters.thresholdFractionChangeMin, 'f', 2));

    PupilHaarOffsetSlider->setDoubleValue(mEyePropertiesParameters.pupilOffset);
    PupilHaarOffsetLabel->setText(QString::number(mEyePropertiesParameters.pupilOffset, 'f', 2));

    GlintRadiusSlider->setValue(mEyePropertiesParameters.glintRadius);
    GlintRadiusLabel->setText(QString::number(mEyePropertiesParameters.glintRadius));

    CurvatureOffsetSlider->setDoubleValue(mEyePropertiesParameters.curvatureOffset);
    CurvatureOffsetLabel->setText(QString::number(mEyePropertiesParameters.curvatureOffset, 'f', 1));

    EdgeMaximumFitNumberSlider->setValue(mEyePropertiesParameters.edgeMaximumFitNumber);
    EdgeMaximumFitNumberLabel->setText(QString::number(mEyePropertiesParameters.edgeMaximumFitNumber));

    EllipseFitErrorMaximumSlider->setDoubleValue(mEyePropertiesParameters.ellipseFitErrorMaximum);
    EllipseFitErrorMaximumLabel->setText(QString::number(mEyePropertiesParameters.ellipseFitErrorMaximum, 'f', 1));
}

void MainWindow::setVariableWidgets(const eyePropertiesVariables& leftEyeVariables)
{
    PupilCircfSliderDouble->setDoubleValue(leftEyeVariables.pupilCircumference);
    PupilCircfLabel->setText(QString::number(leftEyeVariables.pupilCircumference, 'f', 1));

    PupilFractSliderDouble->setDoubleValue(leftEyeVariables.pupilFraction);
    PupilFractLabel->setText(QString::number(leftEyeVariables.pupilFraction, 'f', 2));

    EdgeIntensitySliderDouble->setDoubleValue(leftEyeVariables.edgeIntensity);
    EdgeIntensityLabel->setText(QString::number(leftEyeVariables.edgeIntensity, 'f', 1));
}

int MainWindow::getCurrentTime()
{
    const boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
    const boost::posix_time::time_duration td = now.time_of_day();
    return td.total_milliseconds(); // return time of day in milliseconds
}

void MainWindow::selectDirectory()
{
    QString directory = QFileDialog::getExistingDirectory(this, tr("Open Directory"), "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    DataDirectoryTextBox->setText(directory);
    dataDirectory = directory.toStdString();
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_F8)
    {
        FlashStandbySlider->setValue(1);
    }
    else if (event->key() == Qt::Key_F5)
    {
        FlashStandbySlider->setValue(0);
    }
}
