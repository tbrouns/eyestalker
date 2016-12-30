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
    // Get current date

    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(currentDate, 80, "%Y_%m_%d", timeinfo);

    // Initialize default values

    mUEyeOpencvCam.setDeviceInfo(5129, 5445);

    APP_EXIT    = false;
    APP_RUNNING = true;
    Parameters::cameraXResolution    = 1280;
    Parameters::cameraYResolution    = 1024;
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
    cameraAOIHghtMin        = 4;
    cameraAOIHghtStepSize   = 2;
    cameraAOIWdthMin        = 32;
    cameraAOIWdthStepSize   = 4;
    cameraFrameRateDesired  = 250;
    cameraPixelClock        = 24;
    cameraSubSamplingFactor = 2;
    camImageHght            = 200;
    camImageWdth            = 480; // size of image in widget
    PROCESSING_ALL_IMAGES   = false;
    PROCESSING_ALL_TRIALS   = false;
    imageIndexOffline       = 0;
    imageTotalOffline       = 0;
    experimentIndex         = 0;
    TRIAL_RECORDING         = false;
    eyeAOIHghtFraction      = 0.50;
    eyeAOIHghtMin           = 75;
    eyeAOIWdthFraction      = 0.50;
    eyeAOIWdthMin           = 100;
    eyeImageHght            = 200;
    eyeImageWdth            = 320;
    FLASH_STANDBY           = false;
    frameCount              = 0;
    GAIN_AUTO               = true;
    GAIN_BOOST              = false;
    guiUpdateFrequency      = 30;
    pupilOffsetMin          =  5;
    pupilOffsetIni          = 15;
    relativeTime            = 0;
    SAVE_EYE_IMAGE          = true;
    startTime               = 0;
    subjectIdentifier       = "";
    trialIndex              = 0;
    trialTimeLength         = 1500;

    Parameters::CAMERA_RUNNING      = false;
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
    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHght);
    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPos);
    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPos);

    QObject::connect(CamEyeAOIWdthSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIWdth(double)));
    QObject::connect(CamEyeAOIHghtSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIHght(double)));
    QObject::connect(CamEyeAOIXPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIXPos(double)));
    QObject::connect(CamEyeAOIYPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCamEyeAOIYPos(double)));

    // Eye AOI sliders

    EyeWdthROISlider = new SliderDouble;
    EyeWdthROISlider->setPrecision(2);
    EyeWdthROISlider->setDoubleRange(0, 1.0);
    EyeWdthROISlider->setOrientation(Qt::Horizontal);
    QObject::connect(EyeWdthROISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEyeROIWdth(double)));

    EyeHghtROISlider = new SliderDouble;
    EyeHghtROISlider->setPrecision(2);
    EyeHghtROISlider->setDoubleRange(0, 1.0);
    EyeHghtROISlider->setOrientation(Qt::Vertical);
    QObject::connect(EyeHghtROISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEyeROIHght(double)));

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
    RealTimeEyeTrackingCheckBox->setChecked(!SAVE_EYE_IMAGE);
    QObject::connect(RealTimeEyeTrackingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRealTimeEyeTracking(int)));

    QLabel *OfflineModeTextBox = new QLabel;
    OfflineModeTextBox->setText("<b>Offline mode: </b>");

    QCheckBox *OfflineModeCheckBox = new QCheckBox;
    OfflineModeCheckBox->setChecked(false);
    QObject::connect(OfflineModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setOfflineMode(int)));

    QHBoxLayout *RealTimeEyeTrackingLayout = new QHBoxLayout;
    RealTimeEyeTrackingLayout->addStretch();
    RealTimeEyeTrackingLayout->addWidget(RealTimeEyeTrackingTextBox);
    RealTimeEyeTrackingLayout->addWidget(RealTimeEyeTrackingCheckBox);
    RealTimeEyeTrackingLayout->addWidget(OfflineModeTextBox);
    RealTimeEyeTrackingLayout->addWidget(OfflineModeCheckBox);
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

    QPushButton *OfflineLoadSessionButton = new QPushButton("Load session");
    QObject::connect(OfflineLoadSessionButton, SIGNAL(clicked(bool)), this, SLOT(loadOfflineSession()));

    OfflineImageSlider = new QSlider;
    OfflineImageSlider->setRange(0, 0);
    OfflineImageSlider->setOrientation(Qt::Horizontal);
    QObject::connect(OfflineImageSlider, SIGNAL(valueChanged(int)), this, SLOT(setOfflineImageFrame(int)));

    OfflineImageFrameTextBox = new QLabel;
    OfflineImageFrameTextBox->setText("<b>0 / 0</b>");
    OfflineImageFrameTextBox->setAlignment(Qt::AlignCenter);
    OfflineImageFrameTextBox->setMinimumWidth(70);

    QPushButton *OfflinePrevImageButton = new QPushButton("<");
    QObject::connect(OfflinePrevImageButton, SIGNAL(clicked(bool)), this, SLOT(prevOfflineImage()));

    QPushButton *OfflineNextImageButton = new QPushButton(">");
    QObject::connect(OfflineNextImageButton, SIGNAL(clicked(bool)), this, SLOT(nextOfflineImage()));

    QLabel *OfflinePupilDetectionTextBox = new QLabel;
    OfflinePupilDetectionTextBox->setText("<b>Detect pupil: </b> ");

    QPushButton *OfflinePupilDetectionOneButton = new QPushButton("Current");
    QObject::connect(OfflinePupilDetectionOneButton, SIGNAL(clicked(bool)), this, SLOT(detectPupilOneFrame()));

    QPushButton *OfflinePupilDetectionAllFramesButton = new QPushButton("All frames");
    QObject::connect(OfflinePupilDetectionAllFramesButton, SIGNAL(clicked(bool)), this, SLOT(detectPupilAllFrames()));

    QPushButton *OfflinePupilDetectionAllTrialsButton = new QPushButton("All trials");
    QObject::connect(OfflinePupilDetectionAllTrialsButton, SIGNAL(clicked(bool)), this, SLOT(detectPupilAllTrials()));

    QPushButton *SavePupilDataButton = new QPushButton("Save");
    QObject::connect(SavePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(offlineSaveExperimentData()));

    QPushButton *CombinePupilDataButton = new QPushButton("Combine");
    QObject::connect(CombinePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(offlineCombineExperimentData()));

    EyeTrackingOfflineWidget = new QWidget;
    QHBoxLayout *EyeTrackingOfflineLayout = new QHBoxLayout(EyeTrackingOfflineWidget);
    EyeTrackingOfflineLayout->addStretch();
    EyeTrackingOfflineLayout->addWidget(OfflineLoadSessionButton);
    EyeTrackingOfflineLayout->addWidget(SavePupilDataButton);
    EyeTrackingOfflineLayout->addWidget(CombinePupilDataButton);
    EyeTrackingOfflineLayout->addWidget(OfflinePrevImageButton);
    EyeTrackingOfflineLayout->addWidget(OfflineImageSlider);
    EyeTrackingOfflineLayout->addWidget(OfflineNextImageButton);
    EyeTrackingOfflineLayout->addWidget(OfflineImageFrameTextBox);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionTextBox);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionOneButton);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionAllFramesButton);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionAllTrialsButton);
    EyeTrackingOfflineLayout->addStretch();

    EyeTrackingOfflineWidget->setVisible(false);

    QLabel* OfflineTrialTitle = new QLabel;
    OfflineTrialTitle->setText("<b>Offline mode - Trial:</b>");

    OfflineTrialSpinBox = new QSpinBox;
    OfflineTrialSpinBox->setValue(1);
    OfflineTrialSpinBox->setMinimum(1);
    OfflineTrialSpinBox->setAlignment(Qt::AlignRight);

    OfflineTrialSlider = new QSlider;
    OfflineTrialSlider->setValue(1);
    OfflineTrialSlider->setMinimum(1);
    OfflineTrialSlider->setOrientation(Qt::Horizontal);

    QObject::connect(OfflineTrialSlider, SIGNAL(valueChanged(int)),this,SLOT(changeOfflineSession(int)));

    QObject::connect(OfflineTrialSlider, SIGNAL(valueChanged(int)),OfflineTrialSpinBox,SLOT(setValue(int)));
    QObject::connect(OfflineTrialSpinBox, SIGNAL(valueChanged(int)),OfflineTrialSlider,SLOT(setValue(int)));

    QGridLayout *OfflineSessionTitleLayout = new QGridLayout;
    OfflineSessionTitleLayout->addWidget(OfflineTrialTitle, 0, 1);
    OfflineSessionTitleLayout->addWidget(OfflineTrialSpinBox, 0, 2);
    OfflineSessionTitleLayout->addWidget(OfflineTrialSlider, 0, 3);
    OfflineSessionTitleLayout->setColumnStretch(0, 1);
    OfflineSessionTitleLayout->setColumnStretch(1, 3);
    OfflineSessionTitleLayout->setColumnStretch(2, 1);
    OfflineSessionTitleLayout->setColumnStretch(3, 3);
    OfflineSessionTitleLayout->setColumnStretch(4, 1);

    QWidget *CameraSettings = new QWidget;
    QGridLayout* CameraOutputLayout = new QGridLayout(CameraSettings);

    CameraOutputLayout->addLayout(OfflineSessionTitleLayout, 0, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(EyeTrackingOfflineWidget, 5, 1, 1, 4, Qt::AlignCenter);

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
    CameraFrameRateDesiredSpinBox->setRange(0, cameraFrameRateUpperLimit);
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

    if (cameraSubSamplingFactor == 2) { CameraSubSamplingCheckBox->setChecked(true);  }
    else                              { CameraSubSamplingCheckBox->setChecked(false); }

    QObject::connect(CameraSubSamplingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setCameraSubSampling(int)));

    CameraParametersWidget = new QWidget;

    QGridLayout *CameraParametersLayout = new QGridLayout(CameraParametersWidget);

    CameraParametersLayout->addWidget(CameraPixelClockTextBox, 0, 0);
    CameraParametersLayout->addWidget(CameraPixelClockSlider,  0, 1);
    CameraParametersLayout->addWidget(CameraPixelClockLabel,   0, 2);

    CameraParametersLayout->addWidget(CameraFrameRateTextBox, 1, 0);
    CameraParametersLayout->addWidget(CameraFrameRateSlider,  1, 1);
    CameraParametersLayout->addWidget(CameraFrameRateLabel,   1, 2);

    QHBoxLayout *CameraFrameRateDesiredLayout = new QHBoxLayout;
    CameraFrameRateDesiredLayout->addWidget(CameraFrameRateDesiredTextBox);
    CameraFrameRateDesiredLayout->addWidget(CameraFrameRateDesiredSpinBox);
    CameraFrameRateDesiredLayout->addStretch();

    CameraParametersLayout->addLayout(CameraFrameRateDesiredLayout, 2, 1);

    CameraParametersLayout->addWidget(CameraExposureTextBox, 3, 0);
    CameraParametersLayout->addWidget(CameraExposureSlider,  3, 1);
    CameraParametersLayout->addWidget(CameraExposureLabel,   3, 2);

    CameraParametersLayout->addWidget(CameraBlackLevelOffsetTextBox, 4, 0);
    CameraParametersLayout->addWidget(CameraBlackLevelOffsetSlider,  4, 1);
    CameraParametersLayout->addWidget(CameraBlackLevelOffsetLabel,   4, 2);

    QHBoxLayout *CameraBlackLevelModeLayout = new QHBoxLayout;
    CameraBlackLevelModeLayout->addWidget(CameraBlackLevelModeTextBox);
    CameraBlackLevelModeLayout->addWidget(CameraBlackLevelModeCheckBox);
    CameraBlackLevelModeLayout->addStretch();

    CameraParametersLayout->addLayout(CameraBlackLevelModeLayout, 5, 1);

    CameraParametersLayout->addWidget(CameraHardwareGainTextBox, 6, 0);
    CameraParametersLayout->addWidget(CameraHardwareGainSlider,  6, 1);
    CameraParametersLayout->addWidget(CameraHardwareGainLabel,   6, 2);

    QHBoxLayout *CameraHardwareGainOptionsLayout = new QHBoxLayout;
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainAutoTextBox);
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainAutoCheckBox);
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainBoostTextBox);
    CameraHardwareGainOptionsLayout->addWidget(CameraHardwareGainBoostCheckBox);
    CameraHardwareGainOptionsLayout->addStretch();

    CameraParametersLayout->addLayout(CameraHardwareGainOptionsLayout, 7, 1);

    CameraParametersLayout->addWidget(CameraSubSamplingTextBox,  8, 0);
    CameraParametersLayout->addWidget(CameraSubSamplingCheckBox, 8, 1);

    // Parameter settings layout

    // Real-time parameter tracking

    // Pupil circumference

    QLabel *PupilCircumferenceTextBox = new QLabel;
    PupilCircumferenceTextBox->setText("<b>Pupil circumference:</b>");

    PupilCircumferenceLabel  = new QLabel();
    PupilCircumferenceSlider = new SliderDouble();
    PupilCircumferenceSlider->setPrecision(1);
    PupilCircumferenceSlider->setDoubleRange(mEyePropertiesParameters.pupilCircumferenceMin, mEyePropertiesParameters.pupilCircumferenceMax);
    PupilCircumferenceSlider->setOrientation(Qt::Horizontal);

    // Pupil aspect ratio

    QLabel *PupilAspectRatioTextBox = new QLabel;
    PupilAspectRatioTextBox->setText("<b>Pupil aspect ratio:</b>");

    PupilAspectRatioLabel  = new QLabel();
    PupilAspectRatioSlider = new SliderDouble();
    PupilAspectRatioSlider->setPrecision(2);
    PupilAspectRatioSlider->setDoubleRange(mEyePropertiesParameters.pupilAspectRatioMin, 1.0);
    PupilAspectRatioSlider->setOrientation(Qt::Horizontal);

    // Edge intensity

    QLabel *EdgeIntensityTextBox = new QLabel;
    EdgeIntensityTextBox->setText("<b>Edge intensity:</b>");

    EdgeIntensityLabel  = new QLabel();
    EdgeIntensitySlider = new SliderDouble();
    EdgeIntensitySlider->setPrecision(1);
    EdgeIntensitySlider->setDoubleRange(0.0, 255.0);
    EdgeIntensitySlider->setOrientation(Qt::Horizontal);

    QObject::connect(PupilCircumferenceSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumference(double)));
    QObject::connect(PupilAspectRatioSlider,   SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilAspectRatio(double)));
    QObject::connect(EdgeIntensitySlider,      SIGNAL(doubleValueChanged(double)), this, SLOT(setEdgeIntensity(double)));

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
    PupilCircumferenceMinSlider->setDoubleRange(0, pupilCircumferenceUpperLimit);
    PupilCircumferenceMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilCircumferenceMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setPupilCircumferenceMin(double)));   

    QLabel *PupilCircumferenceMaxTextBox = new QLabel;
    PupilCircumferenceMaxTextBox->setText("<b>Circumference max:</b>");

    PupilCircumferenceMaxLabel  = new QLabel;
    PupilCircumferenceMaxSlider = new SliderDouble;
    PupilCircumferenceMaxSlider->setPrecision(1);
    PupilCircumferenceMaxSlider->setDoubleRange(0, pupilCircumferenceUpperLimit);
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

    QLabel *CannyLowerLimitTextBox = new QLabel;
    CannyLowerLimitTextBox->setText("<b>Lower limit:</b>");

    CannyLowerLimitLabel  = new QLabel;
    CannyLowerLimitSlider = new QSlider;
    CannyLowerLimitSlider->setRange(1, mEyePropertiesParameters.cannyUpperLimit);
    CannyLowerLimitSlider->setOrientation(Qt::Horizontal);
    CannyLowerLimitSlider->setSingleStep(1);
    QObject::connect(CannyLowerLimitSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyLowerLimit(int)));

    QLabel *CannyUpperLimitTextBox = new QLabel;
    CannyUpperLimitTextBox->setText("<b>Upper limit:</b>");

    CannyUpperLimitLabel  = new QLabel;
    CannyUpperLimitSlider = new QSlider;
    CannyUpperLimitSlider->setRange(mEyePropertiesParameters.cannyLowerLimit, 4 * mEyePropertiesParameters.cannyLowerLimit);
    CannyUpperLimitSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyUpperLimitSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyUpperLimit(int)));

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

    AlphaAverageSlider      ->setPrecision(2);
    AlphaPredictionSlider   ->setPrecision(2);
    AlphaMiscellaneousSlider->setPrecision(2);
    AlphaMomentumSlider     ->setPrecision(2);

    AlphaAverageSlider      ->setDoubleRange(0, 1.0);
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
    ThresholdCircumferenceTextBox->setText("<b>Pupil circumference:</b>");

    QLabel *ThresholdAspectRatioTextBox = new QLabel;
    ThresholdAspectRatioTextBox->setText("<b>Pupil fraction:</b>");

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
    PupilHaarOffsetSlider->setRange(pupilOffsetMin, 30);
    PupilHaarOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(PupilHaarOffsetSlider, SIGNAL(valueChanged(int)), this, SLOT(setPupilHaarOffset(int)));   

    QLabel *GlintRadiusTextBox = new QLabel;
    GlintRadiusTextBox->setText("<b>Glint radius:</b>");

    GlintRadiusLabel  = new QLabel;
    GlintRadiusSlider = new QSlider;
    GlintRadiusSlider->setRange(0, 10);
    GlintRadiusSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GlintRadiusSlider, SIGNAL(valueChanged(int)), this, SLOT(setGlintRadius(int)));   

    QLabel *CurvatureOffsetTextBox = new QLabel;
    CurvatureOffsetTextBox->setText("<b>Curvature offset:</b>");

    CurvatureOffsetLabel  = new QLabel;
    CurvatureOffsetSlider = new SliderDouble;
    CurvatureOffsetSlider->setPrecision(1);
    CurvatureOffsetSlider->setDoubleRange(0, 30);
    CurvatureOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureOffset(double)));  

    QLabel *EdgeMaximumFitNumberTextBox = new QLabel;
    EdgeMaximumFitNumberTextBox->setText("<b>Edge maximum fit number:</b>");

    EdgeMaximumFitNumberLabel  = new QLabel;
    EdgeMaximumFitNumberSlider = new QSlider;
    EdgeMaximumFitNumberSlider->setRange(0, 7);
    EdgeMaximumFitNumberSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EdgeMaximumFitNumberSlider, SIGNAL(valueChanged(int)), this, SLOT(setEdgeMaximumFitNumber(int)));

    QLabel *EllipseFitErrorMaximumTextBox = new QLabel;
    EllipseFitErrorMaximumTextBox->setText("<b>Ellipse maximum fit error:</b>");

    EllipseFitErrorMaximumLabel  = new QLabel;
    EllipseFitErrorMaximumSlider = new SliderDouble;
    EllipseFitErrorMaximumSlider->setPrecision(1);
    EllipseFitErrorMaximumSlider->setDoubleRange(0, 80);
    EllipseFitErrorMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EllipseFitErrorMaximumSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEllipseFitErrorMaximum(double)));  

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

    FlashThresholdSlider = new QSlider;
    FlashThresholdSlider->setRange(0, 255);
    FlashThresholdSlider->setValue(0);
    FlashThresholdSlider->setOrientation(Qt::Horizontal);

    FlashThresholdLabel = new QLabel;
    FlashThresholdLabel->setText(QString::number(flashThreshold));

    QObject::connect(FlashThresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(setFlashThreshold(int)));

    QPushButton* FlashMinIntensityResetButton = new QPushButton("Reset");
    QObject::connect(FlashMinIntensityResetButton, SIGNAL(clicked()), this, SLOT(resetFlashMinIntensity()));

    QLabel* FlashThresholdTextBox = new QLabel;
    FlashThresholdTextBox->setText("<b>Flash threshold:</b>");

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

    QPushButton* StartRecordingButton = new QPushButton("Start");
    QObject::connect(StartRecordingButton, SIGNAL(clicked()), this, SLOT(startRecordingManual()));

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
    ExperimentTabLayout->addWidget(StartRecordingButton, 4, 2);
    ExperimentTabLayout->addWidget(FlashStandbyTextBox, 5, 0);
    ExperimentTabLayout->addWidget(FlashStandbySlider, 5, 1);
    ExperimentTabLayout->addWidget(FlashStandbyLabel, 5, 2);
    ExperimentTabLayout->addWidget(FlashThresholdTextBox, 6, 0);
    ExperimentTabLayout->addWidget(FlashThresholdSlider, 6, 1);
    ExperimentTabLayout->addWidget(FlashThresholdLabel, 6, 2);
    ExperimentTabLayout->addWidget(FlashMinIntensityResetButton, 6, 3, Qt::AlignLeft);
    ExperimentTabLayout->addLayout(FlashCoordinatesLayout, 8, 0, 1, 3);

    ExperimentTabLayout->setColumnStretch(0, 1);
    ExperimentTabLayout->setColumnStretch(1, 2);
    ExperimentTabLayout->setColumnStretch(2, 1);
    ExperimentTabLayout->setColumnStretch(3, 1);
    ExperimentTabLayout->setColumnStretch(4, 5);

    // Set-up parameter settings layout

    QWidget* RealTimeVariablesWidget = new QWidget;
    QGridLayout *RealTimeVariablesLayout = new QGridLayout(RealTimeVariablesWidget);
    RealTimeVariablesLayout->addWidget(PupilCircumferenceTextBox, 0, 0);
    RealTimeVariablesLayout->addWidget(PupilCircumferenceSlider, 0, 1);
    RealTimeVariablesLayout->addWidget(PupilCircumferenceLabel, 0, 2);
    RealTimeVariablesLayout->addWidget(PupilAspectRatioTextBox, 1, 0);
    RealTimeVariablesLayout->addWidget(PupilAspectRatioSlider, 1, 1);
    RealTimeVariablesLayout->addWidget(PupilAspectRatioLabel, 1 ,2);
    RealTimeVariablesLayout->addWidget(EdgeIntensityTextBox, 2, 0);
    RealTimeVariablesLayout->addWidget(EdgeIntensitySlider, 2, 1);
    RealTimeVariablesLayout->addWidget(EdgeIntensityLabel, 2, 2);

    QWidget *ParameterLimitsWidget = new QWidget;
    QGridLayout *ParameterLimitsLayout = new QGridLayout(ParameterLimitsWidget);
    ParameterLimitsLayout->addWidget(PupilCircumferenceMaxTextBox, 0, 0);
    ParameterLimitsLayout->addWidget(PupilCircumferenceMaxSlider, 0, 1);
    ParameterLimitsLayout->addWidget(PupilCircumferenceMaxLabel, 0, 2);
    ParameterLimitsLayout->addWidget(PupilCircumferenceMinTextBox, 1, 0);
    ParameterLimitsLayout->addWidget(PupilCircumferenceMinSlider, 1, 1);
    ParameterLimitsLayout->addWidget(PupilCircumferenceMinLabel, 1, 2);
    ParameterLimitsLayout->addWidget(PupilAspectRatioMinTextBox, 2, 0);
    ParameterLimitsLayout->addWidget(PupilAspectRatioMinSlider, 2, 1);
    ParameterLimitsLayout->addWidget(PupilAspectRatioMinLabel, 2, 2);
    ParameterLimitsLayout->addWidget(EdgeIntensityOffsetTextBox, 3, 0);
    ParameterLimitsLayout->addWidget(EdgeIntensityOffsetSlider, 3, 1);
    ParameterLimitsLayout->addWidget(EdgeIntensityOffsetLabel, 3, 2);

    QWidget *LearningRateWidget = new QWidget;
    QGridLayout *LearningRateLayout = new QGridLayout(LearningRateWidget);
    LearningRateLayout->addWidget(AlphaPredictionTextBox, 0, 0);
    LearningRateLayout->addWidget(AlphaPredictionSlider, 0, 1);
    LearningRateLayout->addWidget(AlphaPredictionLabel, 0, 2);
    LearningRateLayout->addWidget(AlphaAverageTextBox, 1, 0);
    LearningRateLayout->addWidget(AlphaAverageSlider, 1, 1);
    LearningRateLayout->addWidget(AlphaAverageLabel, 1, 2);
    LearningRateLayout->addWidget(AlphaMomentumTextBox, 2, 0);
    LearningRateLayout->addWidget(AlphaMomentumSlider, 2, 1);
    LearningRateLayout->addWidget(AlphaMomentumLabel, 2, 2);
    LearningRateLayout->addWidget(AlphaMiscellaneousTextBox, 3, 0);
    LearningRateLayout->addWidget(AlphaMiscellaneousSlider, 3, 1);
    LearningRateLayout->addWidget(AlphaMiscellaneousLabel, 3, 2);

    QWidget *ThresholdParametersWidget = new QWidget;
    QGridLayout *ThresholdParametersLayout = new QGridLayout(ThresholdParametersWidget);
    ThresholdParametersLayout->addWidget(ThresholdCircumferenceTextBox, 0, 0);
    ThresholdParametersLayout->addWidget(ThresholdCircumferenceSlider,  0, 1);
    ThresholdParametersLayout->addWidget(ThresholdCircumferenceLabel,   0, 2);
    ThresholdParametersLayout->addWidget(ThresholdAspectRatioTextBox,   1, 0);
    ThresholdParametersLayout->addWidget(ThresholdAspectRatioSlider,    1, 1);
    ThresholdParametersLayout->addWidget(ThresholdAspectRatioLabel,     1, 2);

    QWidget *CannyEdgeWidget = new QWidget;
    QGridLayout *CannyEdgeLayout = new QGridLayout(CannyEdgeWidget);
    CannyEdgeLayout->addWidget(CannyUpperLimitTextBox, 0, 0);
    CannyEdgeLayout->addWidget(CannyUpperLimitSlider,  0, 1);
    CannyEdgeLayout->addWidget(CannyUpperLimitLabel,   0, 2);
    CannyEdgeLayout->addWidget(CannyLowerLimitTextBox, 1, 0);
    CannyEdgeLayout->addWidget(CannyLowerLimitSlider,  1, 1);
    CannyEdgeLayout->addWidget(CannyLowerLimitLabel,   1, 2);
    CannyEdgeLayout->addWidget(CannyKernelSizeTextBox, 2, 0);
    CannyEdgeLayout->addWidget(CannyKernelSizeSlider,  2, 1);
    CannyEdgeLayout->addWidget(CannyKernelSizeLabel,   2, 2);
    CannyEdgeLayout->addWidget(CannyBlurLevelTextBox,  3, 0);
    CannyEdgeLayout->addWidget(CannyBlurLevelSlider,   3, 1);
    CannyEdgeLayout->addWidget(CannyBlurLevelLabel,    3, 2);

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
    EyeTrackingParameterTabWidget->addTab(CameraParametersWidget,        tr("Camera"));
    EyeTrackingParameterTabWidget->addTab(RealTimeVariablesWidget,       tr("Variables"));
    EyeTrackingParameterTabWidget->addTab(ParameterLimitsWidget,         tr("Limits"));
    EyeTrackingParameterTabWidget->addTab(CannyEdgeWidget,               tr("Canny-edge"));
    EyeTrackingParameterTabWidget->addTab(LearningRateWidget,            tr("Learning rates"));
    EyeTrackingParameterTabWidget->addTab(ThresholdParametersWidget,     tr("Thresholds"));
    EyeTrackingParameterTabWidget->addTab(MiscellaneousParametersWidget, tr("Miscellaneous"));
    EyeTrackingParameterTabWidget->addTab(ExperimentTabWidget,           tr("Experimental"));

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

    // Set widgets to initial values

    setParameterWidgets();
    resetVariables();

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

        int eyeAOIXPosTemp = 0;
        int eyeAOIYPosTemp = 0;
        int eyeAOIWdthTemp = 0;
        int eyeAOIHghtTemp = 0;

        int flashAOIXPosTemp = 0;
        int flashAOIYPosTemp = 0;
        int flashAOIWdthTemp = 0;
        int flashAOIHghtTemp = 0;

        bool FLASH_REGION_VISIBLE = true;

        imageInfo mImageInfo = mUEyeOpencvCam.getFrame(); // get new frame from camera

        cv::Mat imageOriginal = mImageInfo.image;
        absoluteTime = mImageInfo.time; // Get frame timestamp

        double relativeTimeNew = (absoluteTime - startTime) / (double) 10000; // in ms

        // ignore frame if time interval was too short (possible camera error)
        if (relativeTimeNew <= (relativeTime + 0.9 * (1000 / cameraFrameRate))) { continue; }

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

            if (!TRIAL_RECORDING)
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
                        FLASH_REGION_VISIBLE = false;
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
                        FLASH_REGION_VISIBLE = false;
                    }
                }

                if (flashAOIYPosTemp < 0)
                {
                    if (flashAOIYPosTemp + flashAOIHghtTemp < 0) // Flash AOI not in camera AOI
                    {
                        FlashStandbySlider->setValue(0);
                        FLASH_REGION_VISIBLE = false;
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
                        FLASH_REGION_VISIBLE = false;
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

        if (!TRIAL_RECORDING)
        {
            double avgIntensity = 0;

            if (FLASH_REGION_VISIBLE)
            {
                cv::Rect flashRegion(flashAOIXPosTemp, flashAOIYPosTemp, flashAOIWdthTemp, flashAOIHghtTemp);
                avgIntensity = flashDetection(imageOriginal(flashRegion));
            }

            if (FLASH_STANDBY)
            {
                // reset variables

                mEyePropertiesTemp.v.pupilAspectRatioPrediction = pupilAspectRatioIni;
                mEyePropertiesTemp.v.pupilCircumferencePrediction = mEyePropertiesTemp.v.pupilCircumferenceAverage;
                mEyePropertiesTemp.v.pupilRadiusPrediction = mEyePropertiesTemp.v.pupilCircumferencePrediction / (2 * M_PI);

                mEyePropertiesTemp.v.momentumAspectRatio = 0;
                mEyePropertiesTemp.v.momentumCircumference = 0;
                mEyePropertiesTemp.v.momentumRadius = 0;

                mEyePropertiesTemp.v.edgeIntensityAverage = edgeIntensityIni;
                mEyePropertiesTemp.v.edgeIntensityPrediction = edgeIntensityIni;

                mEyePropertiesTemp.v.xVelocity = 0;
                mEyePropertiesTemp.v.yVelocity = 0;

                mEyePropertiesTemp.v.searchRadius = imgWdth;
                mEyePropertiesTemp.v.thresholdCircumferenceChange = mEyePropertiesTemp.p.pupilCircumferenceMax;
                mEyePropertiesTemp.v.thresholdAspectRatioChange = 1.0;

                if (avgIntensity > flashThreshold)
                {
                    startTime = mImageInfo.time;
                    startTrialRecording();
                }
            }
            else // Default mode
            {
                if (avgIntensity > flashMinIntensity)
                {
                    flashMinIntensity = avgIntensity;
                    FlashThresholdSlider->setMinimum(flashMinIntensity);
                }

                cv::Rect eyeRegion(eyeAOIXPosTemp, eyeAOIYPosTemp, eyeAOIWdthTemp, eyeAOIHghtTemp);
                mEyePropertiesTemp = pupilDetection(imageOriginal(eyeRegion), mEyePropertiesTemp); // Pupil tracking algorithm
            }
        }
        else // Trial recording
        {
            if (!SAVE_EYE_IMAGE)
            {
                cv::Rect eyeRegion(eyeAOIXPosTemp, eyeAOIYPosTemp, eyeAOIWdthTemp, eyeAOIHghtTemp);
                mEyePropertiesTemp = pupilDetection(imageOriginal(eyeRegion), mEyePropertiesTemp); // Pupil tracking algorithm

                timeStamps[frameCount] = relativeTime; // save time stamps

                mEyePropertiesTemp.v.xPosAbsolute = mEyePropertiesTemp.v.xPosExact + eyeAOIXPosTemp;
                mEyePropertiesTemp.v.yPosAbsolute = mEyePropertiesTemp.v.yPosExact + eyeAOIYPosTemp;

                eyeXPositions[frameCount] = mEyePropertiesTemp.v.xPosAbsolute;
                eyeYPositions[frameCount] = mEyePropertiesTemp.v.yPosAbsolute;

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
                         << currentDate   << "/"
                         << (NameInputLineEdit->text()).toStdString()
                         << "/trial_"     << trialIndex
                         << "/raw/"       << frameCount
                         << ".png";

                std::vector<int> compression_params;
                compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
                compression_params.push_back(0);

                cv::Mat imageOriginalGray;
                cv::cvtColor(imageOriginal, imageOriginalGray, cv::COLOR_BGR2GRAY);
                cv::imwrite(filename.str(), imageOriginalGray, compression_params);

                frameCount++;
            }

            if (frameCount >= trialFrameTotal)
            {
                mUEyeOpencvCam.stopRecording();
                saveTrialData();
                trialIndex++;
                TrialIndexSpinBox->setValue(trialIndex);
                TRIAL_RECORDING = false;
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
        if (!TRIAL_RECORDING)
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

                setVariableWidgets(mEyePropertiesTemp.v); // update sliders

                if (CameraHardwareGainAutoCheckBox->checkState())
                {
                    int hardwareGain = mUEyeOpencvCam.getHardwareGain();
                    CameraHardwareGainSlider->setValue(hardwareGain);
                    CameraHardwareGainLabel->setText(QString::number(hardwareGain));
                }

                // increase pixel clock if desired frame-rate has not been reached
                int desiredFrameRate = CameraFrameRateDesiredSpinBox->value();

                if ((cameraPixelClock < CameraPixelClockSlider->maximum()) && (desiredFrameRate > cameraFrameRate))
                {
                    cameraPixelClock = cameraPixelClock + 1;
                    CameraPixelClockSlider->setValue(cameraPixelClock);

                    if (cameraFrameRate > desiredFrameRate)
                    {
                        CameraFrameRateSlider->setDoubleValue(desiredFrameRate);
                    }
                }

                if (desiredFrameRate < cameraFrameRate - 1.0)
                {
                    cameraPixelClock = cameraPixelClock - 1;
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

void MainWindow::setOfflineMode(int state)
{
    if (!state)
    {
        imageIndexOffline = 0;

        EyeTrackingParameterTabWidget->setUpdatesEnabled(false);
        EyeTrackingParameterTabWidget->insertTab(0, CameraParametersWidget, tr("Camera"));
        EyeTrackingParameterTabWidget->setUpdatesEnabled(true);

        CamEyeAOIWdthSlider->setVisible(true);
        CamEyeAOIHghtSlider->setVisible(true);
        CamEyeAOIXPosSlider->setVisible(true);
        CamEyeAOIYPosSlider->setVisible(true);

        EyeTrackingOfflineWidget->setVisible(false);

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
        startOfflineSession();
    }
}

void MainWindow::setParameterWidgets()
{
    EyeWdthROISlider->setDoubleValue(eyeAOIWdthFraction);
    EyeHghtROISlider->setDoubleValue(eyeAOIHghtFraction);

    PupilCircumferenceMinSlider->setDoubleValue(mEyePropertiesParameters.pupilCircumferenceMin);
    PupilCircumferenceMinLabel ->setText(QString::number(mEyePropertiesParameters.pupilCircumferenceMin, 'f', 1));

    PupilCircumferenceMaxSlider->setDoubleValue(mEyePropertiesParameters.pupilCircumferenceMax);
    PupilCircumferenceMaxLabel ->setText(QString::number(mEyePropertiesParameters.pupilCircumferenceMax, 'f', 1));

    PupilAspectRatioMinSlider->setDoubleValue(mEyePropertiesParameters.pupilAspectRatioMin);
    PupilAspectRatioMinLabel ->setText(QString::number(mEyePropertiesParameters.pupilAspectRatioMin, 'f', 2));

    EdgeIntensityOffsetSlider->setDoubleValue(mEyePropertiesParameters.edgeIntensityOffset);
    EdgeIntensityOffsetLabel ->setText(QString::number(mEyePropertiesParameters.edgeIntensityOffset, 'f', 1));

    CannyLowerLimitSlider->setValue(mEyePropertiesParameters.cannyLowerLimit);
    CannyLowerLimitLabel ->setText(QString::number(mEyePropertiesParameters.cannyLowerLimit));

    CannyUpperLimitSlider->setValue(mEyePropertiesParameters.cannyUpperLimit);
    CannyUpperLimitLabel ->setText(QString::number(mEyePropertiesParameters.cannyUpperLimit));

    CannyBlurLevelSlider->setValue(mEyePropertiesParameters.cannyBlurLevel);
    CannyBlurLevelLabel ->setText(QString::number(mEyePropertiesParameters.cannyBlurLevel));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mEyePropertiesParameters.cannyKernelSize));
    CannyKernelSizeLabel ->setText(QString::number(mEyePropertiesParameters.cannyKernelSize));

    AlphaAverageSlider      ->setDoubleValue(mEyePropertiesParameters.alphaAverage);
    AlphaPredictionSlider   ->setDoubleValue(mEyePropertiesParameters.alphaPrediction);
    AlphaMiscellaneousSlider->setDoubleValue(mEyePropertiesParameters.alphaMiscellaneous);
    AlphaMomentumSlider     ->setDoubleValue(mEyePropertiesParameters.alphaMomentum);

    AlphaAverageLabel      ->setText(QString::number(mEyePropertiesParameters.alphaAverage, 'f', 2));
    AlphaPredictionLabel   ->setText(QString::number(mEyePropertiesParameters.alphaPrediction, 'f', 2));
    AlphaMiscellaneousLabel->setText(QString::number(mEyePropertiesParameters.alphaMiscellaneous, 'f', 2));
    AlphaMomentumLabel     ->setText(QString::number(mEyePropertiesParameters.alphaMomentum, 'f', 2));

    ThresholdCircumferenceSlider->setDoubleValue(mEyePropertiesParameters.circumferenceChangeThreshold);
    ThresholdCircumferenceLabel ->setText(QString::number(mEyePropertiesParameters.circumferenceChangeThreshold, 'f', 1));

    ThresholdAspectRatioSlider->setDoubleValue(mEyePropertiesParameters.aspectRatioChangeThreshold);
    ThresholdAspectRatioLabel ->setText(QString::number(mEyePropertiesParameters.aspectRatioChangeThreshold, 'f', 2));

    PupilHaarOffsetSlider->setValue(mEyePropertiesParameters.pupilOffset);
    PupilHaarOffsetLabel ->setText(QString::number(mEyePropertiesParameters.pupilOffset));

    GlintRadiusSlider->setValue(mEyePropertiesParameters.glintRadius);
    GlintRadiusLabel->setText(QString::number(mEyePropertiesParameters.glintRadius));

    CurvatureOffsetSlider->setDoubleValue(mEyePropertiesParameters.curvatureOffsetMin);
    CurvatureOffsetLabel->setText(QString::number(mEyePropertiesParameters.curvatureOffsetMin, 'f', 1));

    EdgeMaximumFitNumberSlider->setValue(mEyePropertiesParameters.edgeMaximumFitNumber);
    EdgeMaximumFitNumberLabel->setText(QString::number(mEyePropertiesParameters.edgeMaximumFitNumber));

    EllipseFitErrorMaximumSlider->setDoubleValue(mEyePropertiesParameters.ellipseFitErrorMaximum);
    EllipseFitErrorMaximumLabel->setText(QString::number(mEyePropertiesParameters.ellipseFitErrorMaximum, 'f', 1));
}

void MainWindow::setVariableWidgets(const eyePropertiesVariables& mEyePropertiesVariables)
{
    PupilCircumferenceSlider->setDoubleValue(mEyePropertiesVariables.pupilCircumferencePrediction);
    PupilCircumferenceLabel->setText(QString::number(mEyePropertiesVariables.pupilCircumferencePrediction, 'f', 1));

    PupilAspectRatioSlider->setDoubleValue(mEyePropertiesVariables.pupilAspectRatioPrediction);
    PupilAspectRatioLabel->setText(QString::number(mEyePropertiesVariables.pupilAspectRatioPrediction, 'f', 2));

    EdgeIntensitySlider->setDoubleValue(mEyePropertiesVariables.edgeIntensityPrediction);
    EdgeIntensityLabel->setText(QString::number(mEyePropertiesVariables.edgeIntensityPrediction, 'f', 1));
}

void MainWindow::resetVariables()
{
    mEyePropertiesVariables.pupilCircumferenceAverage    = 0.5 * (mEyePropertiesParameters.pupilCircumferenceMax + mEyePropertiesParameters.pupilCircumferenceMin);
    mEyePropertiesVariables.pupilCircumferencePrediction = mEyePropertiesVariables.pupilCircumferenceAverage;
    mEyePropertiesVariables.pupilAspectRatioAverage      = pupilAspectRatioIni;
    mEyePropertiesVariables.pupilAspectRatioPrediction   = pupilAspectRatioIni;
    mEyePropertiesVariables.edgeIntensityAverage         = edgeIntensityIni;
    mEyePropertiesVariables.edgeIntensityPrediction      = edgeIntensityIni;

    mEyePropertiesVariables.momentumCircumference = 0;
    mEyePropertiesVariables.momentumAspectRatio   = 0;
    mEyePropertiesVariables.momentumRadius        = 0;

    mEyePropertiesVariables.xPosExact     = 0.5 * Parameters::eyeAOIWdth;
    mEyePropertiesVariables.yPosExact     = 0.5 * Parameters::eyeAOIHght;
    mEyePropertiesVariables.xPosPredicted = mEyePropertiesVariables.xPosExact;
    mEyePropertiesVariables.yPosPredicted = mEyePropertiesVariables.yPosExact;
    mEyePropertiesVariables.xVelocity     = 0;
    mEyePropertiesVariables.yVelocity     = 0;

    mEyePropertiesVariables.searchRadius = 0.5 * Parameters::eyeAOIWdth;

    setVariableWidgets(mEyePropertiesVariables);
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
