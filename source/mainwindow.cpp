//  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

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

#include "mainwindow.h"

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

    Parameters::CAMERA_READY      = false;
    Parameters::CAMERA_RUNNING    = false;
    Parameters::ONLINE_PROCESSING = true;

    TRIAL_RECORDING = false;
    SAVE_EYE_IMAGE  = true;
    FLASH_STANDBY   = false;

    Parameters::ellipseDrawCrossSize    = 5;
    Parameters::ellipseDrawOutlineWidth = 0.032;

    cameraPixelClock = 24;

    frameCount         = 0;
    guiUpdateFrequency = 30;
    relativeTime       = 0;

    startTime         = 0;
    subjectIdentifier = "";
    trialIndex        = 0;

    mParameterWidgetEye  = new ParameterWidget;
    mParameterWidgetBead = new ParameterWidget;

    mVariableWidgetEye  = new VariableWidget;
    mVariableWidgetBead = new VariableWidget;

    // AOI

    Parameters::cameraXResolution = 1280;
    Parameters::cameraYResolution = 1024;

    cameraAOIFractionHghtDefaultLeft = 0.19;
    cameraAOIFractionHghtDefaultRght = 0.22;
    cameraAOIFractionWdthDefaultLeft = 0.25;
    cameraAOIFractionWdthDefaultRght = 0.31;
    cameraAOIFractionXPosDefaultLeft = 0.20;
    cameraAOIFractionXPosDefaultRght = 0.52;
    cameraAOIFractionYPosDefaultLeft = 0.41;
    cameraAOIFractionYPosDefaultRght = 0.37;

    camImageHght = 200;
    camImageWdth = 480; // size of image in widget

    eyeImageHght = 200;
    eyeImageWdth = 320; // size of image in widget

    eyeAOIHghtMin = 75;
    eyeAOIWdthMin = 100;

    cameraAOIHghtMin = 4;
    cameraAOIWdthMin = 32;

    cameraAOIHghtStepSize = 2;
    cameraAOIWdthStepSize = 4;

    // Offline mode

    PROCESSING_ALL_IMAGES = false;
    PROCESSING_ALL_TRIALS = false;
    PROCESSING_ALL_EXPS   = false;

    imageIndexOffline = 0;
    imageTotalOffline = 0;

    // Advanced options

    mAdvancedOptions.CURVATURE_MEASUREMENT = false;

    SAVE_DATA_FIT   = false;
    SAVE_DATA_EDGE  = false;
    SAVE_DATA_EXTRA = false;

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

    CameraHardwareGainAutoCheckBox  = new QCheckBox;
    CameraHardwareGainBoostCheckBox = new QCheckBox;

    loadSettings(LastUsedSettingsFileName);

    Parameters::drawFlags.haar = true;
    Parameters::drawFlags.edge = true;
    Parameters::drawFlags.elps = true;

    // Camera feed

    cv::Mat imgCam(camImageWdth, camImageWdth, CV_8UC3, cv::Scalar(150, 150, 150));
    cv::Mat imgEye(eyeImageWdth, eyeImageHght, CV_8UC3, cv::Scalar(150, 150, 150));

    CamQImage = new QImageOpenCV(1);
    CamQImage->setSize(camImageWdth, camImageHght);
    CamQImage->setAOIEye  (Parameters::eyeAOI);
    CamQImage->setAOIBead (Parameters::beadAOI);
    CamQImage->setAOIFlash(flashAOI);

    CamQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    CamQImage->loadImage(imgCam);
    CamQImage->setImage();
    QObject::connect(CamQImage, SIGNAL(updateImage(int)), this , SLOT(onUpdateImageRaw(int)));

    EyeQImage  = new QImageOpenCV(2);
    EyeQImage->setSize(eyeImageWdth, eyeImageHght);
    EyeQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    EyeQImage->loadImage(imgEye);
    EyeQImage->setImage();
    QObject::connect(EyeQImage, SIGNAL(imageMouseClick(double, double)), this, SLOT(onSetPupilPosition(double, double)));

    // Cam AOI sliders

    CamEyeAOIWdthSlider->setDoubleValue(cameraAOIFractionWdth);
    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHght);
    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPos);
    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPos);

    QObject::connect(CamEyeAOIWdthSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamEyeAOIWdth(double)));
    QObject::connect(CamEyeAOIHghtSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamEyeAOIHght(double)));
    QObject::connect(CamEyeAOIXPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamEyeAOIXPos(double)));
    QObject::connect(CamEyeAOIYPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamEyeAOIYPos(double)));

    // Eye AOI sliders

    EyeWdthAOISlider = new SliderDouble;
    EyeWdthAOISlider->setPrecision(2);
    EyeWdthAOISlider->setDoubleRange(0, 1.0);
    EyeWdthAOISlider->setOrientation(Qt::Horizontal);
    QObject::connect(EyeWdthAOISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetEyeAOIWdth(double)));

    EyeHghtAOISlider = new SliderDouble;
    EyeHghtAOISlider->setPrecision(2);
    EyeHghtAOISlider->setDoubleRange(0, 1.0);
    EyeHghtAOISlider->setOrientation(Qt::Vertical);
    QObject::connect(EyeHghtAOISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetEyeAOIHght(double)));

    EyeWdthAOISlider->setDoubleValue(eyeAOIWdthFraction);
    EyeHghtAOISlider->setDoubleValue(eyeAOIHghtFraction);

    QPushButton* AOICropButton = new QPushButton("&Crop AOI");
    QObject::connect(AOICropButton, SIGNAL(clicked()), this, SLOT(onCropAOI()));

    // Qwt plot

    mQwtPlotWidget = new QwtPlotWidget;
    mQwtPlotWidget->setWidth(350);
    mQwtPlotWidget->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    mQwtPlotWidget->setVisible(false);
    QObject::connect(this, SIGNAL(showPlot()), this, SLOT(onPlotTrialData()));

    /////////////////// Draw checkboxes ///////////////////////

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

    QObject::connect(HaarCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetDrawHaar(int)));
    QObject::connect(EdgeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetDrawEdge(int)));
    QObject::connect(ElpsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetDrawElps(int)));

    QHBoxLayout *DrawFunctionsLayout = new QHBoxLayout;
    DrawFunctionsLayout->addStretch();
    DrawFunctionsLayout->addWidget(HaarTextBox);
    DrawFunctionsLayout->addWidget(HaarCheckBox);
    DrawFunctionsLayout->addWidget(EdgeTextBox);
    DrawFunctionsLayout->addWidget(EdgeCheckBox);
    DrawFunctionsLayout->addWidget(ElpsTextBox);
    DrawFunctionsLayout->addWidget(ElpsCheckBox);
    DrawFunctionsLayout->addStretch();

    /////////////////// Options checkboxes ///////////////////////

    QLabel *RealTimeEyeTrackingTextBox = new QLabel;
    RealTimeEyeTrackingTextBox->setText("Real-time eye tracking:");

    RealTimeEyeTrackingCheckBox = new QCheckBox;
    RealTimeEyeTrackingCheckBox->setChecked(!SAVE_EYE_IMAGE);
    QObject::connect(RealTimeEyeTrackingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetRealTimeTracking(int)));

    QLabel *OfflineModeTextBox = new QLabel;
    OfflineModeTextBox->setText("Offline mode:");

    QCheckBox *OfflineModeCheckBox = new QCheckBox;
    OfflineModeCheckBox->setChecked(false);
    QObject::connect(OfflineModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetOfflineMode(int)));

    QPushButton *ResetParametersPushButton = new QPushButton("Reset parameters");
    QObject::connect(ResetParametersPushButton, SIGNAL(clicked()), this, SLOT(onResetParameters()));

    QLabel *BeadDetectionTextBox = new QLabel;
    BeadDetectionTextBox->setText("Bead detection:");

    BeadDetectionCheckBox = new QCheckBox;
    BeadDetectionCheckBox->setChecked(mParameterWidgetBead->getState());
    QObject::connect(BeadDetectionCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetBeadDetection(int)));

    QHBoxLayout *OptionsLayout = new QHBoxLayout;
    OptionsLayout->addStretch();
    OptionsLayout->addWidget(RealTimeEyeTrackingTextBox);
    OptionsLayout->addWidget(RealTimeEyeTrackingCheckBox);
    OptionsLayout->addWidget(OfflineModeTextBox);
    OptionsLayout->addWidget(OfflineModeCheckBox);
    OptionsLayout->addWidget(BeadDetectionTextBox);
    OptionsLayout->addWidget(BeadDetectionCheckBox);
    OptionsLayout->addWidget(ResetParametersPushButton);
    OptionsLayout->addStretch();

    QPushButton *AOILeftEyeButton = new QPushButton("&Left eye");
    QObject::connect(AOILeftEyeButton, SIGNAL(clicked()), this, SLOT(onSetAOIEyeLeft()));

    QPushButton *AOIRghtEyeButton = new QPushButton("&Right eye");
    QObject::connect(AOIRghtEyeButton, SIGNAL(clicked()), this, SLOT(onSetAOIEyeRght()));

    QHBoxLayout *EyeDetectionsLayout = new QHBoxLayout;
    EyeDetectionsLayout->addStretch();
    EyeDetectionsLayout->addWidget(AOILeftEyeButton);
    EyeDetectionsLayout->addWidget(AOIRghtEyeButton);
    EyeDetectionsLayout->addWidget(AOICropButton);
    EyeDetectionsLayout->addStretch();

    /////////////////// Offline mode ///////////////////////

    QPushButton *OfflineLoadSessionButton = new QPushButton("Load session");
    QObject::connect(OfflineLoadSessionButton, SIGNAL(clicked(bool)), this, SLOT(onLoadSession()));

    OfflineImageSlider = new QSlider;
    OfflineImageSlider->setRange(0, 0);
    OfflineImageSlider->setOrientation(Qt::Horizontal);
    QObject::connect(OfflineImageSlider, SIGNAL(valueChanged(int)), this, SLOT(onSetOfflineImage(int)));

    OfflineImageFrameTextBox = new QLabel;
    OfflineImageFrameTextBox->setText("<b>0 / 0</b>");
    OfflineImageFrameTextBox->setAlignment(Qt::AlignCenter);
    OfflineImageFrameTextBox->setMinimumWidth(70);

    QPushButton *OfflinePrevImageButton = new QPushButton("<");
    QObject::connect(OfflinePrevImageButton, SIGNAL(clicked(bool)), this, SLOT(onImagePrevious()));

    QPushButton *OfflineNextImageButton = new QPushButton(">");
    QObject::connect(OfflineNextImageButton, SIGNAL(clicked(bool)), this, SLOT(onImageNext()));

    QLabel *OfflineDetectionTextBox = new QLabel;
    OfflineDetectionTextBox->setText("<b>Detect pupil in:</b>");

    QPushButton *OfflineDetectOneFrameButton = new QPushButton("One frame");
    QObject::connect(OfflineDetectOneFrameButton, SIGNAL(clicked(bool)), this, SLOT(onDetectCurrentFrame()));

    QPushButton *OfflineDetectAllFramesButton = new QPushButton("All frames");
    QObject::connect(OfflineDetectAllFramesButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllFrames()));

    QPushButton *OfflineDetectAllTrialsButton = new QPushButton("All trials");
    QObject::connect(OfflineDetectAllTrialsButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllTrials()));

    QPushButton *OfflineDetectAllExperimentsButton = new QPushButton("All experiments");
    QObject::connect(OfflineDetectAllExperimentsButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllExperiments()));

    QPushButton *SavePupilDataButton = new QPushButton("Save");
    QObject::connect(SavePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(onSaveTrialData()));

    QPushButton *CombinePupilDataButton = new QPushButton("Combine");
    QObject::connect(CombinePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(onCombineData()));

    OfflineModeMainWidget = new QWidget;
    QHBoxLayout *EyeTrackingOfflineLayout = new QHBoxLayout(OfflineModeMainWidget);
    EyeTrackingOfflineLayout->addStretch();
    EyeTrackingOfflineLayout->addWidget(OfflineLoadSessionButton);
    EyeTrackingOfflineLayout->addWidget(SavePupilDataButton);
    EyeTrackingOfflineLayout->addWidget(CombinePupilDataButton);
    EyeTrackingOfflineLayout->addWidget(OfflinePrevImageButton);
    EyeTrackingOfflineLayout->addWidget(OfflineImageSlider);
    EyeTrackingOfflineLayout->addWidget(OfflineNextImageButton);
    EyeTrackingOfflineLayout->addWidget(OfflineImageFrameTextBox);
    EyeTrackingOfflineLayout->addWidget(OfflineDetectionTextBox);
    EyeTrackingOfflineLayout->addWidget(OfflineDetectOneFrameButton);
    EyeTrackingOfflineLayout->addWidget(OfflineDetectAllFramesButton);
    EyeTrackingOfflineLayout->addWidget(OfflineDetectAllTrialsButton);
    EyeTrackingOfflineLayout->addWidget(OfflineDetectAllExperimentsButton);
    EyeTrackingOfflineLayout->addStretch();

    OfflineModeMainWidget->setVisible(false);

    QLabel* OfflineTrialTitle = new QLabel;
    OfflineTrialTitle->setText("<b>Offline mode - Trial:</b>");

    OfflineTrialSpinBox = new QSpinBox;
    OfflineTrialSpinBox->setValue(0);
    OfflineTrialSpinBox->setAlignment(Qt::AlignRight);

    OfflineTrialSlider = new QSlider;
    OfflineTrialSlider->setValue(0);
    OfflineTrialSlider->setOrientation(Qt::Horizontal);

    QObject::connect(OfflineTrialSlider,  SIGNAL(valueChanged(int)), this, SLOT(onSetTrialOffline(int)));

    QObject::connect( OfflineTrialSlider, SIGNAL(valueChanged(int)), OfflineTrialSpinBox, SLOT(setValue(int)));
    QObject::connect(OfflineTrialSpinBox, SIGNAL(valueChanged(int)),  OfflineTrialSlider, SLOT(setValue(int)));

    OfflineModeHeaderWidget = new QWidget;
    QGridLayout *OfflineSessionTitleLayout = new QGridLayout(OfflineModeHeaderWidget);
    OfflineSessionTitleLayout->addWidget(OfflineTrialTitle, 0, 1);
    OfflineSessionTitleLayout->addWidget(OfflineTrialSpinBox, 0, 2);
    OfflineSessionTitleLayout->addWidget(OfflineTrialSlider, 0, 3);
    OfflineSessionTitleLayout->setColumnStretch(0, 1);
    OfflineSessionTitleLayout->setColumnStretch(1, 3);
    OfflineSessionTitleLayout->setColumnStretch(2, 1);
    OfflineSessionTitleLayout->setColumnStretch(3, 3);
    OfflineSessionTitleLayout->setColumnStretch(4, 1);

    OfflineModeHeaderWidget->setVisible(false);

    QWidget *CameraSettings = new QWidget;
    QGridLayout* CameraOutputLayout = new QGridLayout(CameraSettings);

    CameraOutputLayout->addWidget(OfflineModeHeaderWidget, 0, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(OfflineModeMainWidget, 5, 1, 1, 4, Qt::AlignCenter);

    CameraOutputLayout->addWidget(CamEyeAOIXPosSlider, 0, 2);
    CameraOutputLayout->addWidget(CamEyeAOIYPosSlider, 1, 1);
    CameraOutputLayout->addWidget(CamQImage,           1, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(CamEyeAOIWdthSlider, 2, 2);
    CameraOutputLayout->addWidget(CamEyeAOIHghtSlider, 1, 3);
    CameraOutputLayout->addWidget(EyeHghtAOISlider,    1, 5);
    CameraOutputLayout->addWidget(EyeQImage,           1, 4, Qt::AlignCenter);
    CameraOutputLayout->addWidget(mQwtPlotWidget,      0, 4, 3, 2);
    CameraOutputLayout->addWidget(EyeWdthAOISlider,    2, 4);
    CameraOutputLayout->addLayout(OptionsLayout,       4, 2);
    CameraOutputLayout->addLayout(EyeDetectionsLayout, 3, 2);
    CameraOutputLayout->addLayout(DrawFunctionsLayout, 3, 4, Qt::AlignCenter);

    CameraOutputLayout->setColumnStretch(0, 1);
    CameraOutputLayout->setColumnStretch(6, 1);

    ///////////////////////////////////////////////////////////////
    /////////////////////// CAMERA TAB  ///////////////////////////
    ///////////////////////////////////////////////////////////////

    // Pixel clock

    QLabel *CameraPixelClockTextBox = new QLabel;
    CameraPixelClockTextBox->setText("<b>Pixel clock (MHz): </b>");

    CameraPixelClockSlider = new QSlider;
    CameraPixelClockSlider->setRange(0, 0);
    CameraPixelClockSlider->setValue(0);
    CameraPixelClockSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraPixelClockSlider, SIGNAL(valueChanged(int)), this, SLOT(onSetCameraPixelClock(int)));

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
    QObject::connect(CameraFrameRateSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCameraFrameRate(double)));

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
    QObject::connect(CameraExposureSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCameraExposure(double)));

    CameraExposureLabel = new QLabel;
    CameraExposureLabel->setText(QString::number(0.0, 'f', 2));

    // Black level

    QLabel *CameraBlackLevelOffsetTextBox = new QLabel;
    CameraBlackLevelOffsetTextBox->setText("<b>Black level correction: </b>");

    CameraBlackLevelOffsetSlider = new QSlider;
    CameraBlackLevelOffsetSlider->setRange(0, 0);
    CameraBlackLevelOffsetSlider->setValue(0);
    CameraBlackLevelOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraBlackLevelOffsetSlider, SIGNAL(valueChanged(int)), this, SLOT(onSetCameraBlackLevelOffset(int)));

    CameraBlackLevelOffsetLabel = new QLabel;
    CameraBlackLevelOffsetLabel->setText(QString::number(0));

    QLabel *CameraBlackLevelModeTextBox = new QLabel;
    CameraBlackLevelModeTextBox->setText("<b>Auto: </b>");

    QCheckBox *CameraBlackLevelModeCheckBox = new QCheckBox;
    CameraBlackLevelModeCheckBox->setChecked(true);
    QObject::connect(CameraBlackLevelModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetCameraBlackLevelMode(int)));

    // Gain level

    QLabel *CameraHardwareGainTextBox = new QLabel;
    CameraHardwareGainTextBox->setText("<b>Gain: </b>");

    CameraHardwareGainSlider = new QSlider;
    CameraHardwareGainSlider->setRange(0, 100);
    CameraHardwareGainSlider->setValue(0);
    CameraHardwareGainSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CameraHardwareGainSlider, SIGNAL(valueChanged(int)), this, SLOT(onSetCameraHardwareGain(int)));

    CameraHardwareGainLabel = new QLabel;
    CameraHardwareGainLabel->setText(QString::number(0));

    QLabel *CameraHardwareGainAutoTextBox = new QLabel;
    CameraHardwareGainAutoTextBox->setText("<b>Auto: </b>");

    QObject::connect(CameraHardwareGainAutoCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetCameraAutoGain(int)));

    QLabel *CameraHardwareGainBoostTextBox = new QLabel;
    CameraHardwareGainBoostTextBox->setText("<b>Boost: </b>");

    QObject::connect(CameraHardwareGainBoostCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetCameraGainBoost(int)));

    // Sub-sampling

    QLabel *CameraSubSamplingTextBox = new QLabel;
    CameraSubSamplingTextBox->setText("<b>Sub-sampling: </b>");

    QCheckBox *CameraSubSamplingCheckBox = new QCheckBox;

    if (cameraSubSamplingFactor == 2) { CameraSubSamplingCheckBox->setChecked(true);  }
    else                              { CameraSubSamplingCheckBox->setChecked(false); }

    QObject::connect(CameraSubSamplingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetCameraSubSampling(int)));

    // Set-up layout

    QWidget* CameraParametersWidget = new QWidget;
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

    ///////////////////////////////////////////////////////////////
    ///////////////////// EXPERIMENT TAB  /////////////////////////
    ///////////////////////////////////////////////////////////////

    // Subject input

    NameInputLineEdit = new QLineEdit;
    NameInputLineEdit->setAlignment(Qt::AlignRight);
    NameInputLineEdit->setText(subjectIdentifier);

    QLabel *NameInputTextBox = new QLabel;
    NameInputTextBox->setText("<b>Subject identifier:</b>");

    // Data directory

    QLabel* DataDirectoryTitleTextBox = new QLabel;
    DataDirectoryTitleTextBox->setText("<b>Data path:</b>");

    DataDirectoryTextBox = new QLabel;
    DataDirectoryTextBox->setText(QString::fromStdString(dataDirectory));

    QPushButton* DataDirectoryButton = new QPushButton("&Browse...");
    QObject::connect(DataDirectoryButton, SIGNAL(clicked()), this, SLOT(onDirectorySelect()));

    // Data filename

    QLabel *DataFilenameTextBox = new QLabel;
    DataFilenameTextBox->setText("<b>File name:</b>");

    DataFilenameLineEdit = new QLineEdit;
    DataFilenameLineEdit->setAlignment(Qt::AlignRight);
    DataFilenameLineEdit->setText(QString::fromStdString(dataFilename));

    // Trial time length

    QLabel *TrialTimeLengthTextBox = new QLabel;
    TrialTimeLengthTextBox->setText("<b>Trial length (ms):</b>");

    TrialTimeLengthLineEdit = new QLineEdit;
    TrialTimeLengthLineEdit->setAlignment(Qt::AlignRight);
    TrialTimeLengthLineEdit->setText(QString::number(trialTimeLength));

    // Flash trigger button

    FlashStandbySlider = new QSlider;
    FlashStandbySlider->setRange(0, 1);
    FlashStandbySlider->setOrientation(Qt::Horizontal);
    QObject::connect(FlashStandbySlider, SIGNAL(valueChanged(int)), this, SLOT(onFlashStandbySlider(int)));

    QLabel* FlashStandbyTextBox = new QLabel;
    FlashStandbyTextBox->setText("<b>Flash standby mode:</b>");

    FlashStandbyLabel = new QLabel;
    FlashStandbyLabel->setText("<font color='red'><b>OFF</b></font>");

    // Flash threshold

    FlashThresholdSlider = new QSlider;
    FlashThresholdSlider->setRange(0, 255);
    FlashThresholdSlider->setValue(0);
    FlashThresholdSlider->setOrientation(Qt::Horizontal);

    FlashThresholdLabel = new QLabel;
    FlashThresholdLabel->setText(QString::number(flashThreshold));

    QObject::connect(FlashThresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashThreshold(int)));

    QPushButton* FlashThresholdResetButton = new QPushButton("Reset");
    QObject::connect(FlashThresholdResetButton, SIGNAL(clicked()), this, SLOT(onResetFlashIntensity()));

    QLabel* FlashThresholdTextBox = new QLabel;
    FlashThresholdTextBox->setText("<b>Flash threshold:</b>");

    // Flash AOI coordinates

    QSpinBox* FlashXPosSpinBox = new QSpinBox;
    QSpinBox* FlashYPosSpinBox = new QSpinBox;
    QSpinBox* FlashWdthSpinBox = new QSpinBox;
    QSpinBox* FlashHghtSpinBox = new QSpinBox;

    FlashXPosSpinBox->setRange(0, Parameters::cameraXResolution);
    FlashYPosSpinBox->setRange(0, Parameters::cameraYResolution);
    FlashWdthSpinBox->setRange(0, Parameters::cameraXResolution);
    FlashHghtSpinBox->setRange(0, Parameters::cameraYResolution);

    FlashXPosSpinBox->setValue(flashAOI.xPos);
    FlashYPosSpinBox->setValue(flashAOI.yPos);
    FlashWdthSpinBox->setValue(flashAOI.wdth);
    FlashHghtSpinBox->setValue(flashAOI.hght);

    QObject::connect(FlashXPosSpinBox, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIXPos(int)));
    QObject::connect(FlashYPosSpinBox, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIYPos(int)));
    QObject::connect(FlashWdthSpinBox, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIWdth(int)));
    QObject::connect(FlashHghtSpinBox, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIHght(int)));

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

    // Trial index

    QLabel* TrialIndexTextBox = new QLabel;
    TrialIndexTextBox->setText("<b>Trial index:</b>");

    TrialIndexSpinBox = new QSpinBox;
    TrialIndexSpinBox->setRange(0, 999);
    TrialIndexSpinBox->setValue(trialIndex);
    TrialIndexSpinBox->setAlignment(Qt::AlignRight);
    QObject::connect(TrialIndexSpinBox, SIGNAL(valueChanged(int)), this, SLOT(onSetTrialIndex(int)));

    QHBoxLayout* TrialIndexSpinBoxLayout = new QHBoxLayout;
    TrialIndexSpinBoxLayout->addStretch();
    TrialIndexSpinBoxLayout->addWidget(TrialIndexSpinBox);

    // Manual record button

    QPushButton* StartRecordingButton = new QPushButton("Start");
    QObject::connect(StartRecordingButton, SIGNAL(clicked()), this, SLOT(onStartRecordingManual()));

    // Data save

    QLabel* SaveDataTextBox = new QLabel;
    SaveDataTextBox->setText("<b>Data to be saved:</b>");

    QLabel* SaveDataAspectRatioTextBox = new QLabel;
    SaveDataAspectRatioTextBox->setText("Aspect ratio:");

    QLabel* SaveDataCircumferenceTextBox = new QLabel;
    SaveDataCircumferenceTextBox->setText("Circumference:");

    QLabel* SaveDataPositionTextBox = new QLabel;
    SaveDataPositionTextBox->setText("Position:");

    QCheckBox* SaveDataAspectRatioCheckBox   = new QCheckBox;
    QCheckBox* SaveDataCircumferenceCheckBox = new QCheckBox;
    QCheckBox* SaveDataPositionCheckBox      = new QCheckBox;

    SaveDataAspectRatioCheckBox  ->setChecked(SAVE_ASPECT_RATIO);
    SaveDataCircumferenceCheckBox->setChecked(SAVE_CIRCUMFERENCE);
    SaveDataPositionCheckBox     ->setChecked(SAVE_POSITION);

    QObject::connect(SaveDataAspectRatioCheckBox,   SIGNAL(stateChanged(int)), this, SLOT(onSetSaveDataAspectRatio(int)));
    QObject::connect(SaveDataCircumferenceCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetSaveDataCircumference(int)));
    QObject::connect(SaveDataPositionCheckBox,      SIGNAL(stateChanged(int)), this, SLOT(onSetSaveDataPosition(int)));

    // Set-up layouts

    QWidget* ExperimentTabWidget     = new QWidget;
    QHBoxLayout *ExperimentTabLayout = new QHBoxLayout(ExperimentTabWidget);

    QGridLayout *ExperimentTabLeftLayout = new QGridLayout;
    QGridLayout *ExperimentTabRghtLayout = new QGridLayout;

    ExperimentTabLeftLayout->addWidget(NameInputTextBox,            0, 0);
    ExperimentTabLeftLayout->addWidget(NameInputLineEdit,           0, 1);
    ExperimentTabLeftLayout->addWidget(DataDirectoryTitleTextBox,   1, 0);
    ExperimentTabLeftLayout->addWidget(DataDirectoryTextBox,        1, 1);
    ExperimentTabLeftLayout->addWidget(DataDirectoryButton,         1, 2);
    ExperimentTabLeftLayout->addWidget(DataFilenameTextBox,         2, 0);
    ExperimentTabLeftLayout->addWidget(DataFilenameLineEdit,        2, 1);
    ExperimentTabLeftLayout->addWidget(TrialIndexTextBox,           3, 0);
    ExperimentTabLeftLayout->addLayout(TrialIndexSpinBoxLayout,     3, 1);
    ExperimentTabLeftLayout->addWidget(TrialTimeLengthTextBox,      4, 0);
    ExperimentTabLeftLayout->addWidget(TrialTimeLengthLineEdit,     4, 1);
    ExperimentTabLeftLayout->addWidget(StartRecordingButton,        4, 2);
    ExperimentTabLeftLayout->addWidget(FlashStandbyTextBox,         5, 0);
    ExperimentTabLeftLayout->addWidget(FlashStandbySlider,          5, 1);
    ExperimentTabLeftLayout->addWidget(FlashStandbyLabel,           5, 2);
    ExperimentTabLeftLayout->addWidget(FlashThresholdTextBox,       6, 0);
    ExperimentTabLeftLayout->addWidget(FlashThresholdSlider,        6, 1);
    ExperimentTabLeftLayout->addWidget(FlashThresholdLabel,         6, 2);
    ExperimentTabLeftLayout->addWidget(FlashThresholdResetButton,   6, 3, Qt::AlignLeft);
    ExperimentTabLeftLayout->addLayout(FlashCoordinatesLayout,      8, 0, 1, 3);

    ExperimentTabLeftLayout->setColumnStretch(0, 1);
    ExperimentTabLeftLayout->setColumnStretch(1, 2);
    ExperimentTabLeftLayout->setColumnStretch(2, 1);
    ExperimentTabLeftLayout->setColumnStretch(3, 1);
    ExperimentTabLeftLayout->setColumnStretch(4, 5);

    ExperimentTabRghtLayout->addWidget(SaveDataTextBox,               0, 0);
    ExperimentTabRghtLayout->addWidget(SaveDataPositionTextBox,       1, 0);
    ExperimentTabRghtLayout->addWidget(SaveDataPositionCheckBox,      1, 1);
    ExperimentTabRghtLayout->addWidget(SaveDataCircumferenceTextBox,  2, 0);
    ExperimentTabRghtLayout->addWidget(SaveDataCircumferenceCheckBox, 2, 1);
    ExperimentTabRghtLayout->addWidget(SaveDataAspectRatioTextBox,    3, 0);
    ExperimentTabRghtLayout->addWidget(SaveDataAspectRatioCheckBox,   3, 1);

    ExperimentTabRghtLayout->setColumnStretch(0, 2);
    ExperimentTabRghtLayout->setColumnStretch(1, 1);
    ExperimentTabRghtLayout->setColumnStretch(2, 3);
    ExperimentTabRghtLayout->setRowStretch(4, 3);

    ExperimentTabLayout->addLayout(ExperimentTabLeftLayout);
    ExperimentTabLayout->addLayout(ExperimentTabRghtLayout);

    ///////////////////////////////////////////////////////////////
    //////////////////// DEVELOPMENT TAB  /////////////////////////
    ///////////////////////////////////////////////////////////////

    QLabel *CurvatureMeasurementTextBox = new QLabel;
    CurvatureMeasurementTextBox->setText("<b>Curvature measurement:</b>");

    QCheckBox *CurvatureMeasurementCheckBox = new QCheckBox;
    QObject::connect(CurvatureMeasurementCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetCurvatureMeasurement(int)));

    QLabel *SaveDataEdgeTextBox = new QLabel;
    SaveDataEdgeTextBox->setText("<b>Save edge data:</b>");

    QCheckBox *SaveDataEdgeCheckBox = new QCheckBox;
    QObject::connect(SaveDataEdgeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetSaveDataEdge(int)));

    QLabel *SaveDataFitTextBox = new QLabel;
    SaveDataFitTextBox->setText("<b>Save fit data:</b>");

    QCheckBox *SaveDataFitCheckBox = new QCheckBox;
    QObject::connect(SaveDataFitCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetSaveDataFit(int)));

    QLabel *SaveDataExtraTextBox = new QLabel;
    SaveDataExtraTextBox->setText("<b>Save extra data:</b>");

    QCheckBox *SaveDataExtraCheckBox = new QCheckBox;
    QObject::connect(SaveDataExtraCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetSaveDataExtra(int)));

    QWidget* AdvancedOptionsWidget = new QWidget;
    QGridLayout* AdvancedOptionsLayout = new QGridLayout(AdvancedOptionsWidget);
    AdvancedOptionsLayout->addWidget(CurvatureMeasurementTextBox,  0, 0, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(CurvatureMeasurementCheckBox, 0, 1, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(SaveDataEdgeTextBox,          1, 0, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(SaveDataEdgeCheckBox,         1, 1, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(SaveDataFitTextBox,           2, 0, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(SaveDataFitCheckBox,          2, 1, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(SaveDataExtraTextBox,         3, 0, Qt::AlignRight);
    AdvancedOptionsLayout->addWidget(SaveDataExtraCheckBox,        3, 1, Qt::AlignRight);

    AdvancedOptionsLayout->setRowStretch(4, 1);
    AdvancedOptionsLayout->setColumnStretch(2, 1);

    /////////////////// Tab layout ///////////////////////

    QWidget* EyeTrackingWidget = new QWidget;
    QVBoxLayout *EyeTrackingLayout = new QVBoxLayout(EyeTrackingWidget);
    EyeTrackingLayout->addWidget(mVariableWidgetEye);
    EyeTrackingLayout->addWidget(mParameterWidgetEye);

    QScrollArea *EyeTrackingScrollArea = new QScrollArea();
    EyeTrackingScrollArea->setWidget(EyeTrackingWidget);
    EyeTrackingScrollArea->setWidgetResizable(true);

    QWidget* BeadTrackingWidget = new QWidget;
    QVBoxLayout *BeadTrackingLayout = new QVBoxLayout(BeadTrackingWidget);
    BeadTrackingLayout->addWidget(mVariableWidgetBead);
    BeadTrackingLayout->addWidget(mParameterWidgetBead);

    BeadTrackingScrollArea = new QScrollArea();
    BeadTrackingScrollArea->setWidget(BeadTrackingWidget);
    BeadTrackingScrollArea->setWidgetResizable(true);

    AdvancedScrollArea = new QScrollArea();
    AdvancedScrollArea->setWidget(AdvancedOptionsWidget);
    AdvancedScrollArea->setWidgetResizable(true);

    MainTabWidget = new QTabWidget;
    MainTabWidget->addTab(CameraParametersWidget,  tr("Camera"));
    MainTabWidget->addTab(ExperimentTabWidget,     tr("Experimental"));
    MainTabWidget->addTab(EyeTrackingScrollArea,   tr("Eye-tracking"));
    MainTabWidget->addTab(BeadTrackingScrollArea,  tr("Bead-tracking"));
    MainTabWidget->addTab(AdvancedScrollArea,   tr("Advanced"));

    MainTabWidget->setTabEnabled(3, 0); // disable bead-tracking
    MainTabWidget->setTabEnabled(4, 0); // disable development mode

    MainTabWidget->setStyleSheet("QTabBar::tab { height: 30px; width: 125px; }");

    QWidget *CentralWidget = new QWidget;
    QVBoxLayout *CentralLayout = new QVBoxLayout(CentralWidget);
    CentralLayout->addWidget(CameraSettings);
    CentralLayout->addWidget(MainTabWidget);

    // Tab widget

    QSize screenResolution = QDesktopWidget().availableGeometry(this).size();
    double screenOffset = 0.9;

    CentralWidget->setMaximumSize(screenResolution * screenOffset);
    CentralWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    QScrollArea *CentralWidgetScrollArea = new QScrollArea();
    CentralWidgetScrollArea->setWidget(CentralWidget);
    CentralWidgetScrollArea->setWidgetResizable(true);

    ///////////////////////////////////////////////////////////////
    //////////////////// EXTERNAL WIDGETS  ////////////////////////
    ///////////////////////////////////////////////////////////////

    // Quit button

    QPushButton *QuitButton = new QPushButton("Quit");
    QObject::connect(QuitButton, SIGNAL(clicked()), this, SLOT(onQuitButtonClicked()));

    QGridLayout *ExternalLayout = new QGridLayout;
    ExternalLayout->addWidget(QuitButton, 0, 1);
    ExternalLayout->setColumnStretch(0, 1);
    ExternalLayout->setColumnStretch(1, 0);

    // Menu bar

    QAction *options = new QAction("&Advanced mode", this);
    options->setCheckable(true);
    QObject::connect(options, SIGNAL(triggered(bool)), this, SLOT(onSetAdvancedMode(bool)));

    QAction *about   = new QAction("&About EyeStalker...", this);
    QObject::connect(about, &QAction::triggered, this, &MainWindow::onDialogueOpen);

    QMenu *file;
    file = menuBar()->addMenu("&Options");
    file->addAction(options);
    file = menuBar()->addMenu("&Help");
    file->addAction(about);


    ///////////////////////////////////////////////////////////////
    /////////////////////// MAIN LAYOUT ///////////////////////////
    ///////////////////////////////////////////////////////////////

    // Main layout

    QVBoxLayout *MainLayout = new QVBoxLayout;
    MainLayout->addWidget(CentralWidgetScrollArea);
    MainLayout->addLayout(ExternalLayout);

    // Central widget set-up

    QWidget* centralWidget = new QWidget();
    centralWidget->setMaximumSize(screenResolution);
    centralWidget->setLayout(MainLayout);
    setCentralWidget(centralWidget);
    resize(0.7 * screenResolution);

    ///////////////////////////////////////////////////////////////
    //////////////////////// START-UP /////////////////////////////
    ///////////////////////////////////////////////////////////////

    // Start camera

    std::thread findCameraThread(&MainWindow::findCamera, this);
    findCameraThread.detach();

    UpdateCameraImageTimer = new QTimer;
    QObject::connect(UpdateCameraImageTimer, SIGNAL(timeout()), this, SLOT(onUpdateCameraImage()));
    QObject::connect(this, SIGNAL(startTimer(int)), UpdateCameraImageTimer, SLOT(start(int)));
    QObject::connect(this, SIGNAL(stopTimer()),     UpdateCameraImageTimer, SLOT(stop()));
    emit startTimer(round(1000 / guiUpdateFrequency));
}

MainWindow::~MainWindow()
{

}

void MainWindow::pupilTracking()
{   
    resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
    resetVariablesHard(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);

    while(APP_RUNNING && Parameters::CAMERA_RUNNING && Parameters::ONLINE_PROCESSING)
    {
        detectionVariables mDetectionVariablesTemp;
        detectionParameters mDetectionParametersTemp;

        dataVariables mDataVariablesTemp;
        drawVariables mDrawVariablesTemp;

        AOIProperties AOICameraTemp;
        AOIProperties AOIFlashTemp;
        AOIProperties AOIEyeTemp;

        imageInfo mImageInfo = mUEyeOpencvCam.getFrame(); // get new frame from camera
        cv::Mat imageOriginal = mImageInfo.image;
        absoluteTime = mImageInfo.time; // Get frame timestamp

        double relativeTimeNew = (absoluteTime - startTime) / (double) 10000; // in ms

        // ignore frame if time interval was too short (possible camera error)
        if (relativeTimeNew <= (relativeTime + 0.9 * (1000 / cameraFrameRate))) { continue; }

        relativeTime = relativeTimeNew;

        { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

            imageCamera = imageOriginal.clone();

            mDetectionVariablesTemp  = mDetectionVariablesEye;
            mDetectionParametersTemp = mParameterWidgetEye->getStructure();

            AOIFlashTemp  = flashAOI;
            AOICameraTemp = Parameters::cameraAOI;
            AOIEyeTemp    = Parameters::eyeAOI;
        }

        AOIProperties AOIFlashRelative;
        bool FLASH_AOI_VISIBLE = false;

        if (!TRIAL_RECORDING)
        {
            FLASH_AOI_VISIBLE = checkFlashAOI(AOIFlashRelative, AOIFlashTemp, AOICameraTemp);
        }

        // Check limits

        int imgWdth = imageOriginal.cols;
        int imgHght = imageOriginal.rows;

        if (imgWdth < (AOIEyeTemp.xPos + AOIEyeTemp.wdth)) { continue; }
        if (imgHght < (AOIEyeTemp.yPos + AOIEyeTemp.hght)) { continue; }

        if (!TRIAL_RECORDING)
        {
            double avgIntensity = 0;

            if (FLASH_AOI_VISIBLE)
            {
                cv::Rect flashRegion(AOIFlashRelative.xPos, AOIFlashRelative.yPos, AOIFlashRelative.wdth, AOIFlashRelative.hght);
                avgIntensity = flashDetection(imageOriginal(flashRegion));
            }

            if (FLASH_STANDBY)
            {
                if (avgIntensity > flashThreshold)
                {
                    resetVariablesSoft(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
                    resetVariablesSoft(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);

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

                mDetectionVariablesTemp = eyeStalker(imageOriginal, AOIEyeTemp, mDetectionVariablesTemp, mDetectionParametersTemp, mDataVariablesTemp, mDrawVariablesTemp); // Pupil tracking algorithm
            }
        }
        else // Trial recording
        {
            if (!SAVE_EYE_IMAGE)
            {
                mDetectionVariablesTemp         = eyeStalker(imageOriginal, AOIEyeTemp, mDetectionVariablesTemp, mDetectionParametersTemp, mDataVariablesTemp, mDrawVariablesTemp); // Pupil tracking algorithm
                mDataVariablesTemp.absoluteXPos = mDataVariablesTemp.exactXPos + AOIEyeTemp.xPos + AOICameraTemp.xPos;
                mDataVariablesTemp.absoluteYPos = mDataVariablesTemp.exactYPos + AOIEyeTemp.yPos + AOICameraTemp.yPos;
                mDataVariablesTemp.timestamp    = relativeTime; // save time stamps
                vDataVariables[frameCount]      = mDataVariablesTemp;
                frameCount++;
            }
            else
            {
                vDataVariables[frameCount].timestamp = relativeTime; // save time stamps

                // Saving camera frame

                std::stringstream filename;
                filename << dataDirectory << "/"
                         << currentDate   << "/"
                         << (NameInputLineEdit->text()).toStdString()
                         << "/trial_"
                         << trialIndex
                         << "/raw/"
                         << frameCount
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
                if (!SAVE_EYE_IMAGE) { emit showPlot(); }
                emit startTimer(round(1000 / guiUpdateFrequency));
            }
        }

        // Update structures

        {
            std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

            mDetectionVariablesEye = mDetectionVariablesTemp;
            mDrawVariables         = mDrawVariablesTemp;
            mDataVariables         = mDataVariablesTemp;
        }
    }

    if (!APP_RUNNING)
    {
        std::unique_lock<std::mutex> mtxLock(mtx);
        APP_EXIT = true;
        cv.notify_all();
    }
    else if (!Parameters::ONLINE_PROCESSING)
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


void MainWindow::onUpdateCameraImage()
{
    if (Parameters::ONLINE_PROCESSING)
    {
        if (!TRIAL_RECORDING)
        {
            if (Parameters::CAMERA_RUNNING)
            {
                drawVariables mDrawVariablesTemp;
                dataVariables mDataVariablesTemp;

                cv::Mat imageOriginal;

                AOIProperties   eyeAOITemp;
                AOIProperties  beadAOITemp;
                AOIProperties AOIFlashTemp;

                { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
                    if (!imageCamera.empty())
                    {
                        imageOriginal = imageCamera.clone();

                        mDrawVariablesTemp = mDrawVariables;
                        mDataVariablesTemp = mDataVariables;

                        eyeAOITemp   = Parameters::eyeAOI;
                        beadAOITemp  = Parameters::beadAOI;
                        AOIFlashTemp = flashAOI;

                    } else { return; }
                }

                CamQImage->loadImage(imageOriginal);
                CamQImage->setAOIEye  (  eyeAOITemp);
                CamQImage->setAOIBead ( beadAOITemp);
                CamQImage->setAOIFlash(AOIFlashTemp);
                CamQImage->setImage();

                if (eyeAOITemp.wdth >= eyeAOIWdthMin && eyeAOITemp.hght >= eyeAOIHghtMin)
                {
                    cv::Mat imageProcessed = imageOriginal.clone();
                    drawAll(imageProcessed, mDrawVariablesTemp);
                    EyeQImage->loadImage(imageProcessed);
                    EyeQImage->setImage();
                }
                else { EyeQImage->setAOIError(); }

                mVariableWidgetEye->setWidgets(mDataVariablesTemp); // update sliders

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
    CameraPixelClockLabel ->setText(QString::number(cameraPixelClock));

    std::vector<double> frameRateRange = mUEyeOpencvCam.getFrameRateRange();

    cameraFrameRate = frameRateRange[1]; // set to max
    CameraFrameRateSlider->setDoubleRange(frameRateRange[0], cameraFrameRate);
    CameraFrameRateSlider->setDoubleValue(cameraFrameRate);
    CameraFrameRateLabel ->setText(QString::number(cameraFrameRate, 'f', 1));

    std::vector<int> blackLevelOffsetRange = mUEyeOpencvCam.getBlackLevelOffsetRange();

    CameraBlackLevelOffsetSlider->setRange(blackLevelOffsetRange[0], blackLevelOffsetRange[1]);
    CameraBlackLevelOffsetLabel ->setText(QString::number(blackLevelOffsetRange[1]));
}

void MainWindow::findCamera()
{
    bool CAMERA_START = false;

    while (APP_RUNNING && !Parameters::CAMERA_RUNNING && Parameters::ONLINE_PROCESSING)
    {
        if (mUEyeOpencvCam.findCamera())
        {
            int retInt = mUEyeOpencvCam.initCamera();

            if      (retInt == 0) { continue; }
            else if (retInt == 1)
            {
                if (!mUEyeOpencvCam.setColorMode()) { continue; }
            }
            else if (retInt == 2)
            {
                mUEyeOpencvCam.exitCamera();

                if (!mUEyeOpencvCam.initCamera()) { continue; }
            }

            if (mUEyeOpencvCam.setSubSampling(cameraSubSamplingFactor))
            {
                if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
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

            Parameters::CAMERA_RUNNING  = true;
            Parameters::CAMERA_READY    = true;

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

// Variables

void MainWindow::resetVariablesHard(detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const AOIProperties& mAOI)
{
    // Reset all variables

    mDetectionVariables.averageAspectRatio   = initialAspectRatio; // close to perfect circle
    mDetectionVariables.averageCircumference = 0.5 * (mDetectionParameters.circumferenceMax + mDetectionParameters.circumferenceMin); // calculate first
    mDetectionVariables.averageCurvature     = initialCurvature;
    mDetectionVariables.averageGradient      = 0;
    mDetectionVariables.averageHeight        = mDetectionVariables.averageCircumference / M_PI;
    mDetectionVariables.averageIntensity     = initialIntensity;
    mDetectionVariables.averageWidth         = mDetectionVariables.averageCircumference / M_PI;

    mDetectionVariables.certaintyAverages      = 0;
    mDetectionVariables.certaintyAveragesPrime = 0;

    resetVariablesSoft(mDetectionVariables, mDetectionParameters, mAOI);
}

void MainWindow::resetVariablesSoft(detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const AOIProperties& mAOI)
{
    // Reset everything but averages

    mDetectionVariables.predictedAngle         = 0;
    mDetectionVariables.predictedAspectRatio   = mDetectionVariables.averageAspectRatio;
    mDetectionVariables.predictedCircumference = mDetectionVariables.averageCircumference;
    mDetectionVariables.predictedCurvature     = mDetectionVariables.averageCurvature;
    mDetectionVariables.predictedGradient      = mDetectionVariables.averageGradient;
    mDetectionVariables.predictedHeight        = mDetectionVariables.averageHeight;
    mDetectionVariables.predictedIntensity     = mDetectionVariables.averageIntensity;
    mDetectionVariables.predictedWidth         = mDetectionVariables.averageWidth;
    mDetectionVariables.predictedXPos          = 0.5 * mAOI.wdth; // centre of image
    mDetectionVariables.predictedYPos          = 0.5 * mAOI.hght;

    mDetectionVariables.momentumAspectRatio   = 0;
    mDetectionVariables.momentumCircumference = 0;
    mDetectionVariables.momentumCurvature     = 0;
    mDetectionVariables.momentumGradient      = 0;
    mDetectionVariables.momentumHeight        = 0;
    mDetectionVariables.momentumIntensity     = 0;
    mDetectionVariables.momentumWidth         = 0;
    mDetectionVariables.momentumXPos          = 0;
    mDetectionVariables.momentumYPos          = 0;

    int AOISize;
    if (mAOI.wdth > mAOI.hght) { AOISize = mAOI.wdth; }
    else                       { AOISize = mAOI.hght; }

    mDetectionVariables.thresholdChangeAspectRatio   = 1.0 - mDetectionParameters.aspectRatioMin;
    mDetectionVariables.thresholdChangeCircumference = std::abs(mDetectionParameters.circumferenceMax - mDetectionParameters.circumferenceMin) / mDetectionParameters.circumferenceMax;;
    mDetectionVariables.thresholdChangePosition      = AOISize;
    mDetectionVariables.thresholdScore               = 0;

    mDetectionVariables.offsetCircumference = mDetectionVariables.thresholdChangeCircumference;

    mDetectionVariables.certaintyFeatures      = 0;
    mDetectionVariables.certaintyPosition      = 0;
    mDetectionVariables.certaintyFeaturesPrime = 0;
    mDetectionVariables.certaintyPositionPrime = 0;

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

void MainWindow::onSetOfflineMode(int state)
{
    MainTabWidget->setTabEnabled(0, state);

    if (!state)
    {
        imageIndexOffline = 0;

        CamEyeAOIWdthSlider->setVisible(true);
        CamEyeAOIHghtSlider->setVisible(true);
        CamEyeAOIXPosSlider->setVisible(true);
        CamEyeAOIYPosSlider->setVisible(true);

        OfflineModeHeaderWidget->setVisible(false);
        OfflineModeMainWidget  ->setVisible(false);

        Parameters::ONLINE_PROCESSING = true;
        Parameters::CAMERA_RUNNING    = true;
        Parameters::CAMERA_READY      = true;

        if (mUEyeOpencvCam.startVideoCapture())
        {
            std::thread pupilTrackingThread(&MainWindow::pupilTracking, this);
            pupilTrackingThread.detach();

            getCameraParameters();
        }
        else
        {
            std::thread findCameraThread(&MainWindow::findCamera, this);
            findCameraThread.detach();
        }

        emit startTimer(round(1000 / guiUpdateFrequency));
    }
    else
    {
        onStartOfflineSession();
    }
}

void MainWindow::onSetAdvancedMode(bool state)
{
    MainTabWidget->setTabEnabled(4, state);
}

int MainWindow::getCurrentTime()
{
    const boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
    const boost::posix_time::time_duration td = now.time_of_day();
    return td.total_milliseconds(); // return time of day in milliseconds
}

void MainWindow::onDirectorySelect()
{
    QString directory = QFileDialog::getExistingDirectory(this, tr("Open Directory"), "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    DataDirectoryTextBox->setText(directory);
    dataDirectory = directory.toStdString();
}

void MainWindow::onDialogueOpen()
{
    QString text = "Copyright 2016 - 2017 Terence Brouns, t.s.n.brouns@gmail.com. <br> <br> "
                   "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; <br> "
                   "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. <br><br>"
                   "See the GNU General Public License for more details. <br><br>"
                   "Included 3rd party libraries: <br><br>"
                   "<b>OpenCV:</b> <br><br>"
                   "Intel License Agreement <br>"
                   "For Open Source Computer Vision Library <br>"
                   "Copyright (C) 2000, Intel Corporation, all rights reserved. <br>"
                   "Copyright (C) 2014, Itseez Inc., all rights reserved. <br>"
                   "Third party copyrights are property of their respective owners. <br><br>"
                   "<b>Eigen:</b> <br><br> http://eigen.tuxfamily.org/ <br><br>"
                   "<b>Boost:</b> <br><br>"
                   "Copyright Joe Coder 2004 - 2006. <br>"
                   "Distributed under the Boost Software License, Version 1.0. <br>"
                   "(See accompanying file LICENSE_1_0.txt or copy at <br>"
                   "http://www.boost.org/LICENSE_1_0.txt) <br><br>"
                   "<b>UEye:</b> <br><br> (c) 2016, IDS Imaging Advanced Systems GmbH <br><br>"
                   "<b>libusb:</b> <br><br> libusb is released under version 2.1 of the GNU Lesser General Public License (LGPL).<br>"
                   "http://libusb.info/ <br>";

    ConfirmationWindow mConfirmationWindow(text, false);
    mConfirmationWindow.setWindowTitle("About EyeStalker");
    mConfirmationWindow.exec();
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    if      (event->key() == Qt::Key_F8) { FlashStandbySlider->setValue(1); }
    else if (event->key() == Qt::Key_F5) { FlashStandbySlider->setValue(0); }
}

void MainWindow::onSetPupilPosition(double xPos, double yPos)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    if (xPos > 0 && xPos < Parameters::eyeAOI.wdth && yPos > 0 && yPos < Parameters::eyeAOI.hght)
    {
        mDetectionVariablesEye.predictedXPos = xPos;
        mDetectionVariablesEye.predictedYPos = yPos;
    }
}

bool MainWindow::checkFlashAOI(AOIProperties& AOIFlashRelative, const AOIProperties& AOIFlash, const AOIProperties& AOICamera)
{
    bool FLASH_AOI_VISIBLE = true;

    AOIFlashRelative      = AOIFlash;
    AOIFlashRelative.xPos = AOIFlash.xPos - AOICamera.xPos;
    AOIFlashRelative.yPos = AOIFlash.yPos - AOICamera.yPos;

    if (AOIFlashRelative.xPos < 0)
    {
        if (AOIFlashRelative.xPos + AOIFlashRelative.wdth < 0) // Flash AOI not in camera AOI
        {
            FlashStandbySlider->setValue(0);
            FLASH_AOI_VISIBLE = false;
        }
        else
        {
            AOIFlashRelative.wdth = AOIFlashRelative.xPos + AOIFlashRelative.wdth;
            AOIFlashRelative.xPos = 0;
        }
    }
    else if (AOIFlashRelative.xPos + AOIFlashRelative.wdth >= AOICamera.wdth)
    {
        AOIFlashRelative.wdth = AOICamera.wdth - AOIFlashRelative.xPos;

        if (AOIFlashRelative.wdth <= 0)
        {
            FlashStandbySlider->setValue(0);
            FLASH_AOI_VISIBLE = false;
        }
    }

    if (AOIFlashRelative.yPos < 0)
    {
        if (AOIFlashRelative.yPos + AOIFlashRelative.hght < 0) // Flash AOI not in camera AOI
        {
            FlashStandbySlider->setValue(0);
            FLASH_AOI_VISIBLE = false;
        }
        else
        {
            AOIFlashRelative.hght = AOIFlashRelative.yPos + AOIFlashRelative.hght;
            AOIFlashRelative.yPos = 0;
        }
    }
    else if (AOIFlashRelative.yPos + AOIFlashRelative.hght >= AOICamera.hght)
    {
        AOIFlashRelative.hght = AOICamera.hght - AOIFlashRelative.yPos;

        if (AOIFlashRelative.hght <= 0)
        {
            FlashStandbySlider->setValue(0);
            FLASH_AOI_VISIBLE = false;
        }
    }

    return FLASH_AOI_VISIBLE;
}

///////////////////////////////////////////////////////////////
////////////////// EXPERIMENT FUNCTIONS  //////////////////////
///////////////////////////////////////////////////////////////

void MainWindow::startTrialRecording()
{
    if (Parameters::CAMERA_RUNNING)
    {
        if (!TRIAL_RECORDING)
        {
            if (!boost::filesystem::exists(dataDirectory))
            {
                QString text = "Data directory does not exist. Please choose a valid path.";
                ConfirmationWindow mConfirmationWindow(text, false);
                mConfirmationWindow.setWindowTitle("Warning");
                mConfirmationWindow.exec();
                return;
            }

            emit stopTimer(); // stop showing camera feed

            // start recording

            TRIAL_RECORDING = true;
            mUEyeOpencvCam.startRecording();

            // get start times

            absoluteTime = startTime;
            relativeTime = 0;

            trialStartTime  = getCurrentTime();
            trialFrameTotal = ceil((trialTimeLength * cameraFrameRate) / 1000);

            trialTimeLength = (TrialTimeLengthLineEdit->text()).toInt();

            frameCount = 0;

            // preallocate vector space for data

            vDataVariables.resize(trialFrameTotal);
            vDataVariablesBead.resize(trialFrameTotal);

            vDetectionVariablesEye.resize(trialFrameTotal);
            vDetectionVariablesBead.resize(trialFrameTotal);

            // Reset variables

            if (SAVE_EYE_IMAGE) // create directories if necessary
            {
                std::stringstream directoryName;
                directoryName << dataDirectory << "/" << currentDate;

                if (!boost::filesystem::exists(directoryName.str()))
                {
                    boost::filesystem::create_directory(directoryName.str().c_str());
                }

                if (!NameInputLineEdit->text().isEmpty())
                {
                    directoryName << "-" << (NameInputLineEdit->text()).toStdString();
                }

                if (!boost::filesystem::exists(directoryName.str()))
                {
                    boost::filesystem::create_directory(directoryName.str().c_str());
                }

                directoryName << "/trial_" << trialIndex;
                if (!boost::filesystem::exists(directoryName.str()))
                {
                    boost::filesystem::create_directory(directoryName.str().c_str());
                }

                directoryName << "/raw/";

                if (!boost::filesystem::exists(directoryName.str()))
                {
                    boost::filesystem::create_directory(directoryName.str().c_str());
                }
            }
        }
        else
        {
            QString text = "Camera not running";
            ConfirmationWindow mConfirmationWindow(text, false);
            mConfirmationWindow.setWindowTitle("Warning");
            mConfirmationWindow.exec();
        }
    }
}

void MainWindow::onStartRecordingManual()
{
    if (!TRIAL_RECORDING && !PROCESSING_ALL_EXPS && !PROCESSING_ALL_TRIALS && !PROCESSING_ALL_IMAGES)
    {
        imageInfo mImageInfo = mUEyeOpencvCam.getFrame();
        startTime = mImageInfo.time;
        startTrialRecording();
    }
}

void MainWindow::saveTrialData()
{
    if (!SAVE_EYE_IMAGE)
    {
        dataFilename = (DataFilenameLineEdit->text()).toStdString();

        std::stringstream filename;
        filename << dataDirectory
                 << "/"
                 << dataFilename;

        std::ofstream file;

        if (!boost::filesystem::exists(filename.str()))
        {
            file.open(filename.str(), std::ios::out | std::ios::ate);
        }
        else
        {
            file.open(filename.str(), std::ios_base::app);
            file << "\n";
        }

        std::string delimiter = ";";

        file << std::setw(3) << std::setfill('0') << trialIndex << delimiter; // print with leading zeros
        file << frameCount << delimiter;                                      // data samples
        file << trialStartTime << delimiter;                                  // system clock time

        file << std::fixed;
        file << std::setprecision(3);

        for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].DETECTED  << delimiter; }
        for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].timestamp << delimiter; } // saving ALL timestamps allows for accurate frame-rate calculation during data processing

        if (SAVE_POSITION)
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].absoluteXPos  << delimiter; }
            for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].absoluteYPos  << delimiter; }
        }

        if (SAVE_CIRCUMFERENCE)
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].exactCircumference << delimiter; }
        }

        if (SAVE_ASPECT_RATIO)
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].exactAspectRatio << delimiter; }
        }

        file.close();
    }
    else
    {
        // save time stamps

        std::stringstream filename;
        filename << dataDirectory << "/"
                 << currentDate   << "/"
                 << (NameInputLineEdit->text()).toStdString()
                 << "/trial_"     << trialIndex
                 << "/"
                 << "timestamps.dat";

        std::ofstream file;

        file.open(filename.str(), std::ios::out | std::ios::trunc); // open file and remove any existing data

        std::string delimiter = " "; // space delimiter allows for easier reading when combining data

        file << std::setw(3) << std::setfill('0') << trialIndex << delimiter; // print with leading zeros
        file << trialStartTime << delimiter; // system clock time

        file << std::fixed;
        file << std::setprecision(3);

        for (int i = 0; i < frameCount; i++) { file << vDataVariables[i].timestamp << delimiter; }

        file.close();
    }
}

void MainWindow::onPlotTrialData()
{
    if (SAVE_POSITION)
    {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> t;

        for (int i = 0; i < frameCount; i++)
        {
            if (vDataVariables[i].DETECTED)
            {
                x.push_back(vDataVariables[i].absoluteXPos - Parameters::cameraAOI.xPos);
                y.push_back(Parameters::eyeAOI.hght - (vDataVariables[i].absoluteYPos - Parameters::cameraAOI.yPos));
                t.push_back(0.001 * vDataVariables[i].timestamp);
            }
        }

        //        mQwtPlotWidget->plotTrajectory(x,y);
        mQwtPlotWidget->plotTimeSeries(x,y,t,0.001*trialTimeLength);
    }
}

void MainWindow::onFlashStandbySlider(int val)
{
    if (Parameters::CAMERA_RUNNING && Parameters::ONLINE_PROCESSING && !TRIAL_RECORDING)
    {
        if (val == 0)
        {
            {
                std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
                FLASH_STANDBY = false;
            }

            FlashStandbyLabel->setText("<font color='red'><b>OFF</b></font>");

            if (CameraHardwareGainAutoCheckBox->isChecked())
            {
                mUEyeOpencvCam.setAutoGain(true);
            }

            // display real-time eye-tracking

            EyeQImage       ->setVisible(true);
            EyeWdthAOISlider->setVisible(true);
            EyeHghtAOISlider->setVisible(true);
            mQwtPlotWidget  ->setVisible(false);
        }
        else
        {
            dataFilename = (DataFilenameLineEdit->text()).toStdString();

            std::stringstream filename;
            filename << dataDirectory
                     << "/"
                     << dataFilename
                     << ".dat";

            if (boost::filesystem::exists(filename.str()))
            {
                QString text = "The file <b>" + QString::fromStdString(dataFilename) + ".dat</b> already exists in <b>" + QString::fromStdString(dataDirectory) + "/</b>. Do you wish to add data to the end of this file?";
                ConfirmationWindow mConfirmationWindow(text);
                mConfirmationWindow.setWindowTitle("Please select option");

                if(mConfirmationWindow.exec() == QDialog::Rejected)
                {
                    FlashStandbySlider->setValue(0);
                    return;
                }
            }

            {
                std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
                FLASH_STANDBY = true;
            }

            FlashStandbyLabel->setText("<font color='green'><b>ON</b></font>");

            if (CameraHardwareGainAutoCheckBox->isChecked())
            {
                mUEyeOpencvCam.setAutoGain(false);
            }

            // display plot window

            EyeQImage       ->setVisible(false);
            EyeWdthAOISlider->setVisible(false);
            EyeHghtAOISlider->setVisible(false);
            mQwtPlotWidget  ->setVisible(true);
        }
    }
    else
    {
        FlashStandbySlider->setValue(0);
    }
}

///////////////////////////////////////////////////////////////
///////////////// OFFLINE MODE FUNCTIONS  /////////////////////
///////////////////////////////////////////////////////////////

void MainWindow::onStartOfflineSession()
{
    trialIndexOffline = 0;
    imageIndexOffline = 0;
    Parameters::ONLINE_PROCESSING = false; // turn off pupil tracking

    { // unlock frame grabbing threads
        std::unique_lock<std::mutex> frameCaptureMutexLock(Parameters::frameCaptureMutex);
        Parameters::frameCaptureCV.notify_all(); // unlock getFrame() thread
    }

    { // wait for threads to finish
        std::unique_lock<std::mutex> lck(mtx);
        while (Parameters::CAMERA_RUNNING) cv.wait(lck);
    }

    // hide camera AOI sliders
    CamEyeAOIWdthSlider->setVisible(false);
    CamEyeAOIHghtSlider->setVisible(false);
    CamEyeAOIXPosSlider->setVisible(false);
    CamEyeAOIYPosSlider->setVisible(false);

    OfflineModeMainWidget  ->setVisible(true);
    OfflineModeHeaderWidget->setVisible(true);

    CamQImage->clearImage();
    EyeQImage->clearImage();

    setupOfflineSession();
}

void MainWindow::countNumTrials()
{
    trialTotalOffline = 0;

    while (1)
    {
        std::stringstream folderName;
        folderName << dataDirectoryOffline.toStdString()
                   << "/images/trial_"
                   << trialTotalOffline;
        if (!boost::filesystem::exists(folderName.str())) { break; }
        trialTotalOffline++;
    }

    OfflineTrialSlider ->setMaximum(trialTotalOffline - 1);
    OfflineTrialSpinBox->setMaximum(trialTotalOffline - 1);
}

void MainWindow::countNumImages()
{
    imageTotalOffline = 0;

    while (1)
    {
        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/images/trial_"
                 << trialIndexOffline
                 << "/raw/"
                 << imageTotalOffline
                 << ".png";

        if (!boost::filesystem::exists(filename.str())) { break; }
        else { imageTotalOffline++; }
    }
}

void MainWindow::onLoadSession()
{
    QString dataDirectoryTemp = QFileDialog::getExistingDirectory(this, tr("Select data directory"), dataDirectoryOffline, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    if (!dataDirectoryTemp.isEmpty())
    {
        dataDirectoryOffline = dataDirectoryTemp;
        timeMatrix.clear();
        setupOfflineSession();
    }
}

void MainWindow::setupOfflineSession()
{
    if (!dataDirectoryOffline.isEmpty())
    {
        countNumTrials();

        if (trialTotalOffline > 0)
        {
            if (trialIndexOffline != 0) { OfflineTrialSlider->setValue(0); } // start with first trial
            else { onSetTrialOffline(0); }

            std::stringstream directoryName;
            directoryName << dataDirectoryOffline.toStdString()
                          << "/images/trial_"
                          << trialIndexOffline
                          << "/processed";

            if (!boost::filesystem::exists(directoryName.str()))
            {    boost::filesystem::create_directory(directoryName.str().c_str()); }

            if (timeMatrix.empty()) // Grab time stamps
            {
                std::stringstream directory;
                directory << dataDirectoryOffline.toStdString()
                          << "/images/timestamps.dat";

                std::ifstream data;
                data.open(directory.str().c_str());

                if (data.is_open())
                {
                    std::string str;

                    for (int i = 0; i < trialTotalOffline && std::getline(data, str); i++)
                    {
                        std::vector<double> times;
                        std::istringstream sin(str);
                        double time;
                        while (sin >> time) { times.push_back(time); }
                        timeMatrix.push_back(times);
                    }
                }
            }
        }
    }
}

void MainWindow::onSetTrialOffline(int index)
{
    if (index >= 0) {

        trialIndexOffline = index;

        countNumImages();

        if (imageTotalOffline > 0)
        {
            vDataVariables.resize(imageTotalOffline);
            vDataVariablesBead.resize(imageTotalOffline);

            vDetectionVariablesEye.resize(imageTotalOffline + 1);
            vDetectionVariablesBead.resize(imageTotalOffline + 1);

            if (index == 0)
            {
                resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
                resetVariablesHard(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);
            }
            else
            {
                resetVariablesSoft(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
                resetVariablesSoft(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);
            }

            vDetectionVariablesEye[0]  = mDetectionVariablesEye;
            vDetectionVariablesBead[0] = mDetectionVariablesBead;

            if (imageIndexOffline != 0) { OfflineImageSlider->setValue(0); } // start with first frame
            else { onSetOfflineImage(0); }

            OfflineImageSlider->setMaximum(imageTotalOffline - 1);

            // Create folder

            std::stringstream directoryPath;
            directoryPath << dataDirectoryOffline.toStdString()
                          << "/images/trial_"
                          << trialIndexOffline
                          << "/processed/";

            if (!boost::filesystem::exists(directoryPath.str()))
            {    boost::filesystem::create_directory(directoryPath.str().c_str()); }
        }
    }
}

void MainWindow::onUpdateImageRaw(int imgIndex) // for signal from qimageopencv
{
    if (imgIndex < 0) { imgIndex = imageIndexOffline; }

    std::stringstream imagePath;
    imagePath << dataDirectoryOffline.toStdString()
              << "/images/trial_"
              << trialIndexOffline
              << "/raw/"
              << imgIndex
              << ".png";

    if (boost::filesystem::exists(imagePath.str()))
    {
        cv::Mat eyeImageRaw = cv::imread(imagePath.str(), CV_LOAD_IMAGE_COLOR);
        CamQImage->loadImage(eyeImageRaw);
        { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
            Parameters::cameraAOI.wdth = eyeImageRaw.cols;
            Parameters::cameraAOI.hght = eyeImageRaw.rows;
            updateEyeAOIx();
            updateEyeAOIy();
            CamQImage->setAOIEye (Parameters::eyeAOI);
            CamQImage->setAOIBead(Parameters::beadAOI);
        }
        CamQImage->setImage();
    } else { CamQImage->clearImage(); }
}

void MainWindow::onUpdateImageProcessed(int imgIndex)
{
    std::stringstream fileName;
    fileName << dataDirectoryOffline.toStdString()
             << "/images/trial_"
             << trialIndexOffline
             << "/processed/"
             << imgIndex
             << ".png";

    if (boost::filesystem::exists(fileName.str()))
    {
        cv::Mat eyeImage = cv::imread(fileName.str(), CV_LOAD_IMAGE_COLOR);
        EyeQImage->loadImage(eyeImage);
        EyeQImage->setImage();
    }
    else
    {
        EyeQImage->clearImage();
    }
}

void MainWindow::onImagePrevious()
{
    if (imageIndexOffline > 0)
    {
        imageIndexOffline--;
        OfflineImageSlider->setValue(imageIndexOffline);
    }
}

void MainWindow::onImageNext()
{
    if (imageIndexOffline < imageTotalOffline - 1)
    {
        imageIndexOffline++;
        OfflineImageSlider->setValue(imageIndexOffline);
    }
}

void MainWindow::onSetOfflineImage(int imgIndex)
{
    if (!PROCESSING_ALL_IMAGES) { imageIndexOffline = imgIndex; }

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        mVariableWidgetEye ->setWidgets(vDataVariables    [imgIndex]);
        mVariableWidgetBead->setWidgets(vDataVariablesBead[imgIndex]);
    }

    if (imageTotalOffline > 0)
    {
        std::stringstream ss;
        ss << "<b>" << imgIndex + 1 << " / " << imageTotalOffline << "</b>";
        QString title = QString::fromStdString(ss.str());
        OfflineImageFrameTextBox->setText(title);

        onUpdateImageRaw(imgIndex);
        onUpdateImageProcessed(imgIndex);
    }
    else
    {
        CamQImage->clearImage();
        EyeQImage->clearImage();
        OfflineImageFrameTextBox->setText("<b>0 / 0</b>");
    }
}

void MainWindow::detectCurrentFrame(int imageIndex)
{
    // Grab raw images

    cv::Mat imageRaw;

    std::stringstream imagePathRaw;
    imagePathRaw << dataDirectoryOffline.toStdString()
                 << "/images/trial_"
                 << trialIndexOffline
                 << "/raw/"
                 << imageIndex
                 << ".png";

    if (boost::filesystem::exists(imagePathRaw.str())) { imageRaw = cv::imread(imagePathRaw.str(), CV_LOAD_IMAGE_COLOR); }
    else                                               { return; }

    // Detect pupil

    int cameraAOIXPos;
    int cameraAOIYPos;

    detectionVariables mDetectionVariablesEyeTemp;
    detectionVariables mDetectionVariablesBeadTemp;

    detectionParameters mDetectionParametersEyeTemp;
    detectionParameters mDetectionParametersBeadTemp;

    AOIProperties AOIEyeTemp;
    AOIProperties AOIBeadTemp;

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

        mDetectionVariablesEyeTemp  = vDetectionVariablesEye[imageIndex];
        mDetectionParametersEyeTemp = mParameterWidgetEye->getStructure();

        mDetectionVariablesBeadTemp  = vDetectionVariablesBead[imageIndex];
        mDetectionParametersBeadTemp = mParameterWidgetBead->getStructure();

        cameraAOIXPos = Parameters::cameraAOI.xPos;
        cameraAOIYPos = Parameters::cameraAOI.yPos;

        Parameters::cameraAOI.wdth = imageRaw.cols;
        Parameters::cameraAOI.hght = imageRaw.rows;

        updateEyeAOIx();
        updateEyeAOIy();

        AOIEyeTemp  = Parameters::eyeAOI;
        AOIBeadTemp = Parameters::beadAOI;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    detectionVariables mDetectionVariablesEyeNew = eyeStalker(imageRaw,
                                                              AOIEyeTemp,
                                                              mDetectionVariablesEyeTemp,
                                                              mDetectionParametersEyeTemp,
                                                              mDataVariables,
                                                              mDrawVariables,
                                                              mAdvancedOptions);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

    // Save data

    mDataVariables.duration     = fp_ms.count();
    mDataVariables.absoluteXPos = mDataVariables.exactXPos + cameraAOIXPos;
    mDataVariables.absoluteYPos = mDataVariables.exactYPos + cameraAOIYPos;
    vDataVariables[imageIndex]  = mDataVariables;

    cv::Mat imageProcessed = imageRaw.clone();
    drawAll(imageProcessed, mDrawVariables);

    detectionVariables mDetectionVariablesBeadNew;

    if (mParameterWidgetBead->getState())
    {
        mDetectionVariablesBeadNew      = eyeStalker(imageRaw, AOIBeadTemp, mDetectionVariablesBeadTemp, mDetectionParametersBeadTemp, mDataVariablesBead, mDrawVariablesBead);
        mDataVariablesBead.absoluteXPos = mDataVariablesBead.exactXPos + cameraAOIXPos;
        mDataVariablesBead.absoluteYPos = mDataVariablesBead.exactYPos + cameraAOIYPos;
        vDataVariablesBead[imageIndex]  = mDataVariablesBead;
        drawAll(imageProcessed, mDrawVariablesBead);
    }

    // Save image

    std::stringstream imagePath;
    imagePath << dataDirectoryOffline.toStdString()
              << "/images/trial_"
              << trialIndexOffline
              << "/processed/"
              << imageIndex
              << ".png";

    cv::imwrite(imagePath.str(), imageProcessed);

    // Record variables for next frame(s)

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        mDetectionVariablesEye = mDetectionVariablesEyeNew;
        vDetectionVariablesEye[imageIndex + 1] = mDetectionVariablesEye;
        if (mParameterWidgetBead->getState())
        {
            mDetectionVariablesBead = mDetectionVariablesBeadNew;
            vDetectionVariablesBead[imageIndex + 1] = mDetectionVariablesBead; }
    }
}

void MainWindow::onDetectCurrentFrame()
{
    detectCurrentFrame(imageIndexOffline);
    onUpdateImageProcessed(imageIndexOffline);
}

void MainWindow::onDetectAllFrames()
{
    if (!PROCESSING_ALL_IMAGES)
    {
        PROCESSING_ALL_IMAGES = true;
        OFFLINE_SAVE_DATA     = false;

        std::thread pupilDetectionThread(&MainWindow::detectAllFrames, this);
        pupilDetectionThread.detach();

        while (!OFFLINE_SAVE_DATA && PROCESSING_ALL_IMAGES)
        {
            std::stringstream ss;
            ss << "<b>" << imageIndexOffline + 1 << " / " << imageTotalOffline << "</b>";
            QString title = QString::fromStdString(ss.str());
            OfflineImageFrameTextBox->setText(title);

            OfflineImageSlider->setValue(imageIndexOffline - 1); // Show progress

            msWait(1000 / guiUpdateFrequency);
        }
    }

    {
        std::unique_lock<std::mutex> lck(mtxOffline);
        PROCESSING_ALL_IMAGES = false;
        cvOffline.notify_one(); // notify save-thread that (this) main-thread has exited while loop
    }
    {
        std::unique_lock<std::mutex> lck(mtxOffline);
        while (OFFLINE_SAVE_DATA) cvOffline.wait(lck); // wait for save-thread to complete saving
    }
}

void MainWindow::detectAllFrames()
{
    int initialIndex = imageIndexOffline; // needed for progressbar

    for (imageIndexOffline = initialIndex; imageIndexOffline < imageTotalOffline && PROCESSING_ALL_IMAGES; imageIndexOffline++)
    {
        detectCurrentFrame(imageIndexOffline);
    }

    imageIndexOffline--; // for-loop overshoots value

    if (PROCESSING_ALL_IMAGES)
    {
        OFFLINE_SAVE_DATA = true;

        { std::unique_lock<std::mutex> lck(mtxOffline);
            while (PROCESSING_ALL_IMAGES) cvOffline.wait(lck); } // wait for main-thread to exit while loop

        onSaveTrialData();

        { std::unique_lock<std::mutex> lck(mtxOffline);
            OFFLINE_SAVE_DATA = false;
            cvOffline.notify_one(); } // notify main-thread that saving has been completed

        OfflineImageSlider->setValue(imageIndexOffline);
    }
}

void MainWindow::onDetectAllTrials()
{
    if (!PROCESSING_ALL_TRIALS)
    {
        PROCESSING_ALL_TRIALS = true;

        for (int iTrial = trialIndexOffline; iTrial < trialTotalOffline && PROCESSING_ALL_TRIALS; iTrial++)
        {
            if (iTrial < 5)  { continue; }
            if (iTrial > 15) { break; } // needs to be removed

            OfflineTrialSlider->setValue(iTrial);
            onDetectAllFrames();

        }
    }

    PROCESSING_ALL_IMAGES = false;
    PROCESSING_ALL_TRIALS = false;
}

void MainWindow::onDetectAllExperiments()
{
    if (!PROCESSING_ALL_EXPS)
    {
        QString mainDirectoryTemp = QFileDialog::getExistingDirectory(this, tr("Select data directory"), dataDirectoryOffline, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
        std::string mainDirectory = mainDirectoryTemp.toStdString();

        if (boost::filesystem::exists(mainDirectory))
        {
            PROCESSING_ALL_EXPS = true;

            for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(mainDirectory), {}))
            {
                if (PROCESSING_ALL_EXPS)
                {
                    std::string subDirectory = entry.path().string();
                    dataDirectoryOffline = QString::fromStdString(subDirectory);
                    timeMatrix.clear();
                    setupOfflineSession();
                    onDetectAllTrials();
                } else { break; }
            }
        }
    }

    PROCESSING_ALL_EXPS   = false;
    PROCESSING_ALL_TRIALS = false;
    PROCESSING_ALL_IMAGES = false;
}


void MainWindow::onSaveTrialData()
{
    std::string delimiter = ";";

    { // save pupil data

        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/images/trial_"
                 << trialIndexOffline
                 << "/tracking_data.dat";

        std::ofstream file;
        file.open(filename.str(), std::ios::trunc); // open file and remove any existing data

        if (timeMatrix.size() > 0) { file << std::setw(3) << std::setfill('0') << timeMatrix[trialIndexOffline][0] << ";"; } // trial index
        file << imageTotalOffline << ";";  // data samples
        if (timeMatrix.size() > 0) { file << (int) timeMatrix[trialIndexOffline][1] << ";"; } // system clock time

        file << std::fixed;
        file << std::setprecision(3);

        // write data

        for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariables[i].DETECTED << delimiter; }
        if (timeMatrix.size() > 0) { for (int i = 0; i < imageTotalOffline; i++) { file << timeMatrix[trialIndexOffline][i + 2] << delimiter; } }

        if (SAVE_POSITION)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariables[i].absoluteXPos  << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariables[i].absoluteYPos  << delimiter; }
        }

        if (SAVE_CIRCUMFERENCE)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariables[i].exactCircumference << delimiter; }
        }

        if (SAVE_ASPECT_RATIO)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariables[i].exactAspectRatio << delimiter; }
        }

        if (mParameterWidgetBead->getState())
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesBead[i].DETECTED      << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesBead[i].absoluteXPos  << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesBead[i].absoluteYPos  << delimiter; }
        }

        if (SAVE_DATA_EXTRA)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedXPos          << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedYPos          << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedCircumference << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedAspectRatio   << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedCurvature     << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedIntensity     << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedGradient      << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDetectionVariablesEye[i].predictedAngle         << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariables[i].duration                       << delimiter; }
        }

        file.close();
    }

    if (SAVE_DATA_EDGE)
    {

        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/images/trial_"
                 << trialIndexOffline
                 << "/edge_data.dat";

        std::ofstream file;
        file.open(filename.str());

        for (int i = 0; i < imageTotalOffline; i++)
        {
            int numEdges = vDataVariables[i].edgeData.size();
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].tag          << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].curvature    << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].curvatureMax << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].curvatureMin << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].length       << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].radius       << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].radiusVar    << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].intensity    << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariables[i].edgeData[j].gradient     << delimiter; }
            file << "\n";
        }

        file.close();
    }

    if (SAVE_DATA_FIT)
    {

        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/images/trial_"
                 << trialIndexOffline
                 << "/fit_data.dat";

        std::ofstream file;
        file.open(filename.str());

        for (int i = 0; i < imageTotalOffline; i++)
        {
            int numFits = vDataVariables[i].ellipseData.size();

            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].tag           << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].xPos          << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].yPos          << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].circumference << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].aspectRatio   << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].fitError      << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].edgeLength    << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariables[i].ellipseData[j].angle         << delimiter; }

            file << "\n";
        }

        file.close();
    }
}

void MainWindow::onCombineData()
{
    dataFilename = (DataFilenameLineEdit->text()).toStdString();
    std::stringstream fileNameWriteSS;
    fileNameWriteSS << dataDirectoryOffline.toStdString()
                    << "/"
                    << dataFilename
                    << ".dat";

    std::string fileNameWrite = fileNameWriteSS.str();

    if (boost::filesystem::exists(fileNameWrite))
    {
        QString text = "The file <b>" + QString::fromStdString(dataFilename) + ".dat</b> already exists in <b>" + dataDirectoryOffline + "/</b>. Do you wish to add data to the end of this file?";
        ConfirmationWindow mConfirmationWindow(text);
        mConfirmationWindow.setWindowTitle("Please select option");

        if(mConfirmationWindow.exec() == QDialog::Rejected) { return; }
    }

    for (int iTrial = 0; iTrial < trialTotalOffline; iTrial++)
    {
        OfflineTrialSlider->setValue(iTrial);

        std::stringstream fileNameRead;
        fileNameRead << dataDirectoryOffline.toStdString()
                     << "/images/trial_"
                     << iTrial
                     << "/tracking_data.dat";

        std::ifstream dataFile;
        dataFile.open(fileNameRead.str().c_str());

        std::string line;
        while (dataFile) { if (!std::getline(dataFile, line)) { break; }}

        std::ofstream file;
        file.open(fileNameWrite, std::ios_base::app);
        file << line;
        file << "\n";
        file.close();
    }
}

///////////////////////////////////////////////////////////////
////////////////// PARAMETER FUNCTIONS  ///////////////////////
///////////////////////////////////////////////////////////////

void MainWindow::loadSettings(QString filename)
{
    static const std::vector<double> parametersEye  = {0.005,   // Alpha average
                                                       0.40,    // Alpha features
                                                       0.20,    // Alpha certainty
                                                       0.75,    // Alpha position
                                                       4,       // Canny blur level
                                                       5,       // Canny kernel size
                                                       300.0,   // Canny threshold low
                                                       600.0,   // Canny threshold high
                                                       8,       // Curvature offset
                                                       0.05,    // Ellipse edge fraction
                                                       3,       // Ellipse fit number maximum
                                                       0.60,    // Ellipse fit error maximum
                                                       12,      // Glint size
                                                       305,     // Circumference max
                                                       60,      // Circumference min
                                                       0.4,     // Aspect ratio min
                                                       0.35,    // Circumference offset
                                                       0.10,    // Circumference change threshold
                                                       0.10,    // Aspect ratio change threshold
                                                       12,      // Displacement change threshold
                                                       0.30,    // Score threshold
                                                       0.60,    // Score difference threshold edge
                                                       0.10,    // Score difference threshold fit
                                                       7};      // Edge window length

    static const std::vector<double> parametersBead = {0.005,   // Alpha average
                                                       0.40,    // Alpha features
                                                       0.20,    // Alpha certainty
                                                       0.75,    // Alpha predicted
                                                       4,       // Canny blur level
                                                       5,       // Canny kernel size
                                                       300.0,   // Canny threshold low
                                                       600.0,   // Canny threshold high
                                                       8,       // Curvature offset
                                                       0.05,    // Ellipse edge fraction
                                                       3,       // Ellipse fit number maximum
                                                       0.60,    // Ellipse fit error maximum
                                                       0,       // Glint size
                                                       130,     // Circumference max
                                                       90,      // Circumference min
                                                       0.8,     // Aspect ratio min
                                                       0.10,    // Circumference offset
                                                       0.05,    // Circumference change threshold
                                                       0.05,    // Aspect ratio change threshold
                                                       10,      // Displacement change threshold
                                                       0.30,    // Score threshold
                                                       0.60,    // Score difference threshold
                                                       0.10,    // Score difference threshold fit
                                                       7};      // Edge window length

    QSettings settings(filename, QSettings::IniFormat);

    cameraAOIFractionHght           = settings.value("CamAOIHghtFraction",          cameraAOIFractionHghtDefaultLeft).toDouble();
    cameraAOIFractionWdth           = settings.value("CamAOIWdthFraction",          cameraAOIFractionWdthDefaultLeft).toDouble();
    cameraAOIFractionXPos           = settings.value("CamAOIXPosFraction",          cameraAOIFractionXPosDefaultLeft).toDouble();
    cameraAOIFractionYPos           = settings.value("CamAOIYPosFraction",          cameraAOIFractionYPosDefaultLeft).toDouble();
    cameraFrameRateDesired          = settings.value("CameraFrameRateDesired",      250).toInt();
    cameraSubSamplingFactor         = settings.value("SubSamplingFactor",           1).toInt();
    dataDirectory                   = settings.value("DataDirectory",               "").toString().toStdString();
    dataDirectoryOffline            = settings.value("DataDirectoryOffline",        "").toString();
    dataFilename                    = settings.value("DataFilename",                "experiment_data").toString().toStdString();
    trialIndexOffline               = settings.value("trialIndexOffline",           0).toInt();
    imageTotalOffline               = settings.value("imageTotalOffline",           0).toInt();
    subjectIdentifier               = settings.value("SubjectName",                 "").toString();
    eyeAOIHghtFraction              = settings.value("AOIHghtFraction",             1.0).toDouble();
    eyeAOIWdthFraction              = settings.value("AOIWdthFraction",             1.0).toDouble();
    beadAOIHghtFraction             = settings.value("AOIBeadHghtFraction",         0.6).toDouble();
    beadAOIWdthFraction             = settings.value("AOIBeadWdthFraction",         0.3).toDouble();
    flashThreshold                  = settings.value("FlashThreshold",              230).toInt();
    Parameters::eyeAOIXPosFraction  = settings.value("AOIXPosRelative",             0.0).toDouble();
    Parameters::eyeAOIYPosFraction  = settings.value("AOIYPosRelative",             0.0).toDouble();
    Parameters::beadAOIXPosFraction = settings.value("AOIBeadXPosRelative",         0.2).toDouble();
    Parameters::beadAOIYPosFraction = settings.value("AOIBeadYPosRelative",         0.5).toDouble();
    flashAOI.hght                   = settings.value("FlashAOIHght",                100).toInt();
    flashAOI.wdth                   = settings.value("FlashAOIWdth",                60).toInt();
    flashAOI.xPos                   = settings.value("FlashAOIXPos",                227).toInt();
    flashAOI.yPos                   = settings.value("FlashAOIYPos",                500).toInt();
    SAVE_ASPECT_RATIO               = settings.value("SaveAspectRatio",             true).toBool();
    SAVE_CIRCUMFERENCE              = settings.value("SaveCircumference",           true).toBool();
    SAVE_POSITION                   = settings.value("SavePosition",                true).toBool();
    SAVE_EYE_IMAGE                  = settings.value("SaveEyeImage",                true).toBool();
    subjectIdentifier               = settings.value("SubjectIdentifier",           "").toString();
    trialTimeLength                 = settings.value("TrialTimeLength",             1500).toInt();

    CameraHardwareGainAutoCheckBox ->setChecked(settings.value("GainAuto",   true).toBool());
    CameraHardwareGainBoostCheckBox->setChecked(settings.value("GainBoost", false).toBool());

    cameraAOIWdthMax = Parameters::cameraXResolution / (double) cameraSubSamplingFactor; // maximum possible AOI size
    cameraAOIHghtMax = Parameters::cameraYResolution / (double) cameraSubSamplingFactor;

    updateCamAOIx();
    updateCamAOIy();

    detectionParameters mDetectionParametersEye  = loadParameters(filename, "Eye",  parametersEye);
    detectionParameters mDetectionParametersBead = loadParameters(filename, "Bead", parametersBead);

    mParameterWidgetEye ->setStructure(mDetectionParametersEye);
    mParameterWidgetBead->setStructure(mDetectionParametersBead);

    resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
    resetVariablesHard(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);
}

detectionParameters MainWindow::loadParameters(QString filename, QString prefix, std::vector<double> parameters)
{
    QSettings settings(filename, QSettings::IniFormat);

    detectionParameters mDetectionParameters;
    mDetectionParameters.alphaAverages                      = settings.value(prefix + "AlphaAverage",                    parameters[0]).toDouble();
    mDetectionParameters.alphaFeatures                      = settings.value(prefix + "AlphaFeatures",                   parameters[1]).toDouble();
    mDetectionParameters.alphaCertainty                     = settings.value(prefix + "AlphaCertainty",                  parameters[2]).toDouble();
    mDetectionParameters.alphaPosition                      = settings.value(prefix + "AlphaPosition",                   parameters[3]).toDouble();
    mDetectionParameters.cannyBlurLevel                     = settings.value(prefix + "CannyBlurLevel",                  parameters[4]).toInt();
    mDetectionParameters.cannyKernelSize                    = settings.value(prefix + "CannyKernelSize",                 parameters[5]).toInt();
    mDetectionParameters.cannyThresholdLow                  = settings.value(prefix + "CannyThresholdLow",               parameters[6]).toDouble();
    mDetectionParameters.cannyThresholdHigh                 = settings.value(prefix + "CannyThresholdHigh",              parameters[7]).toDouble();
    mDetectionParameters.curvatureOffset                    = settings.value(prefix + "CurvatureOffset",                 parameters[8]).toDouble();
    mDetectionParameters.fitEdgeFraction                    = settings.value(prefix + "FitEdgeFraction",                 parameters[9]).toDouble();
    mDetectionParameters.fitEdgeMaximum                     = settings.value(prefix + "FitEdgeMaximum",                  parameters[10]).toInt();
    mDetectionParameters.thresholdFitError                  = settings.value(prefix + "ThresholdFitError",               parameters[11]).toDouble();
    mDetectionParameters.glintWdth                          = settings.value(prefix + "GlintSize",                       parameters[12]).toInt();
    mDetectionParameters.circumferenceMax                   = settings.value(prefix + "CircumferenceMax",                parameters[13]).toDouble();
    mDetectionParameters.circumferenceMin                   = settings.value(prefix + "CircumferenceMin",                parameters[14]).toDouble();
    mDetectionParameters.aspectRatioMin                     = settings.value(prefix + "AspectRatioMin",                  parameters[15]).toDouble();
    mDetectionParameters.circumferenceOffset                = settings.value(prefix + "CircumferenceOffset",             parameters[16]).toDouble();
    mDetectionParameters.thresholdChangeCircumference       = settings.value(prefix + "CircumferenceChangeThreshold",    parameters[17]).toDouble();
    mDetectionParameters.thresholdChangeAspectRatio         = settings.value(prefix + "AspectRatioChangeThreshold",      parameters[18]).toDouble();
    mDetectionParameters.thresholdChangePosition            = settings.value(prefix + "DisplacementChangeThreshold",     parameters[19]).toDouble();
    mDetectionParameters.thresholdScore                     = settings.value(prefix + "ScoreThreshold",                  parameters[20]).toDouble();
    mDetectionParameters.thresholdScoreDiffEdge             = settings.value(prefix + "ScoreThresholdDiffEdge",          parameters[21]).toDouble();
    mDetectionParameters.thresholdScoreDiffFit              = settings.value(prefix + "ScoreThresholdDiffFit",           parameters[22]).toDouble();
    mDetectionParameters.windowLengthEdge                   = settings.value(prefix + "WindowLengthEdge",                parameters[23]).toDouble();

    return mDetectionParameters;
}

void MainWindow::saveSettings(QString filename)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue("AOIHghtFraction",                eyeAOIHghtFraction);
    settings.setValue("AOIWdthFraction",                eyeAOIWdthFraction);
    settings.setValue("AOIXPosRelative",                Parameters::eyeAOIXPosFraction);
    settings.setValue("AOIYPosRelative",                Parameters::eyeAOIYPosFraction);
    settings.setValue("AOIBeadHghtFraction",            beadAOIHghtFraction);
    settings.setValue("AOIBeadWdthFraction",            beadAOIWdthFraction);
    settings.setValue("AOIBeadXPosRelative",            Parameters::beadAOIXPosFraction);
    settings.setValue("AOIBeadYPosRelative",            Parameters::beadAOIYPosFraction);
    settings.setValue("DataDirectory",                  QString::fromStdString(dataDirectory));
    settings.setValue("DataDirectoryOffline",           dataDirectoryOffline);
    settings.setValue("GainAuto",                       CameraHardwareGainAutoCheckBox ->checkState());
    settings.setValue("GainBoost",                      CameraHardwareGainBoostCheckBox->checkState());
    settings.setValue("CamAOIHghtFraction",             cameraAOIFractionHght);
    settings.setValue("CamAOIWdthFraction",             cameraAOIFractionWdth);
    settings.setValue("CamAOIXPosFraction",             cameraAOIFractionXPos);
    settings.setValue("CamAOIYPosFraction",             cameraAOIFractionYPos);
    settings.setValue("CameraFrameRateDesired",         cameraFrameRateDesired);
    settings.setValue("DataFilename",                   QString::fromStdString(dataFilename));
    settings.setValue("FlashAOIHght",                   flashAOI.hght);
    settings.setValue("FlashAOIWdth",                   flashAOI.wdth);
    settings.setValue("FlashAOIXPos",                   flashAOI.xPos);
    settings.setValue("FlashAOIYPos",                   flashAOI.yPos);
    settings.setValue("FlashThreshold",                 flashThreshold);
    settings.setValue("SaveAspectRatio",                SAVE_ASPECT_RATIO);
    settings.setValue("SaveCircumference",              SAVE_CIRCUMFERENCE);
    settings.setValue("SavePosition",                   SAVE_POSITION);
    settings.setValue("SaveEyeImage",                   SAVE_EYE_IMAGE);
    settings.setValue("SubjectName",                    subjectIdentifier);
    settings.setValue("SubSamplingFactor",              cameraSubSamplingFactor);
    settings.setValue("TrialTimeLength",                TrialTimeLengthLineEdit->text().toInt());

    detectionParameters mDetectionParametersEye  = mParameterWidgetEye->getStructure();
    saveParameters(filename,  "Eye", mDetectionParametersEye);

    detectionParameters mDetectionParametersBead = mParameterWidgetBead->getStructure();
    saveParameters(filename, "Bead", mDetectionParametersBead);
}

void MainWindow::saveParameters(QString filename, QString prefix, detectionParameters mDetectionParameters)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue(prefix + "AlphaAverage",                   mDetectionParameters.alphaAverages);
    settings.setValue(prefix + "AlphaFeatures",                  mDetectionParameters.alphaFeatures);
    settings.setValue(prefix + "AlphaPosition",                  mDetectionParameters.alphaPosition);
    settings.setValue(prefix + "AlphaCertainty",                 mDetectionParameters.alphaCertainty);
    settings.setValue(prefix + "CannyBlurLevel",                 mDetectionParameters.cannyBlurLevel);
    settings.setValue(prefix + "CannyKernelSize",                mDetectionParameters.cannyKernelSize);
    settings.setValue(prefix + "CannyThresholdLow",              mDetectionParameters.cannyThresholdLow);
    settings.setValue(prefix + "CannyThresholdHigh",             mDetectionParameters.cannyThresholdHigh);
    settings.setValue(prefix + "CircumferenceOffset",            mDetectionParameters.circumferenceOffset);
    settings.setValue(prefix + "CircumferenceMax",               mDetectionParameters.circumferenceMax);
    settings.setValue(prefix + "CircumferenceMin",               mDetectionParameters.circumferenceMin);
    settings.setValue(prefix + "CurvatureOffset",                mDetectionParameters.curvatureOffset);
    settings.setValue(prefix + "FitEdgeMaximum",                 mDetectionParameters.fitEdgeMaximum);
    settings.setValue(prefix + "ThresholdFitError",              mDetectionParameters.thresholdFitError);
    settings.setValue(prefix + "AspectRatioMin",                 mDetectionParameters.aspectRatioMin);
    settings.setValue(prefix + "GlintSize",                      mDetectionParameters.glintWdth);
    settings.setValue(prefix + "CircumferenceChangeThreshold",   mDetectionParameters.thresholdChangeCircumference);
    settings.setValue(prefix + "AspectRatioChangeThreshold",     mDetectionParameters.thresholdChangeAspectRatio);
    settings.setValue(prefix + "DisplacementChangeThreshold",    mDetectionParameters.thresholdChangePosition);
    settings.setValue(prefix + "ScoreThreshold",                 mDetectionParameters.thresholdScore);
    settings.setValue(prefix + "ScoreThresholdDiffEdge",         mDetectionParameters.thresholdScoreDiffEdge);
    settings.setValue(prefix + "ScoreThresholdDiffFit",          mDetectionParameters.thresholdScoreDiffFit);
    settings.setValue(prefix + "WindowLengthEdge",               mDetectionParameters.windowLengthEdge);
}

void MainWindow::onResetParameters()
{
    QString text = "Do you wish to reset all parameters to their default values?";
    ConfirmationWindow mConfirmationWindow(text);
    mConfirmationWindow.setWindowTitle("Please select option");

    if(mConfirmationWindow.exec() == QDialog::Rejected) { return; }

    QString filename = "";
    loadSettings(filename);

    mParameterWidgetEye ->reset();
    mParameterWidgetBead->reset();

    resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
    resetVariablesHard(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);
}

void MainWindow::onSetBeadDetection(int state)
{
    if (!TRIAL_RECORDING && !PROCESSING_ALL_TRIALS && !PROCESSING_ALL_IMAGES)
    {
        MainTabWidget->setTabEnabled(3, state);
        mParameterWidgetBead->setState(state);

        if (!state) { MainTabWidget->setCurrentIndex(0); } // set to camera tab

        CamQImage->showAOIBead(state);
        CamQImage->setImage();
    }
    else
    {
        BeadDetectionCheckBox->setChecked(false);
    }
}

void MainWindow::onSetRealTimeTracking(int state)
{
    if (Parameters::ONLINE_PROCESSING && !TRIAL_RECORDING)
    {
        if (state) { SAVE_EYE_IMAGE = false; }
        else       { SAVE_EYE_IMAGE = true;  }
    }
    else
    {
        RealTimeEyeTrackingCheckBox->setChecked(false);
    }
}

void MainWindow::onResetFlashIntensity()
{
    flashMinIntensity = 0;
    FlashThresholdSlider->setMinimum(0);
    FlashThresholdSlider->setValue(0);
}

void MainWindow::onSetCameraPixelClock(int value)
{
    cameraPixelClock = value;
    mUEyeOpencvCam.setPixelClock(cameraPixelClock);
    CameraPixelClockLabel->setText(QString::number(cameraPixelClock));

    // Set new frame rate

    std::vector<double> frameRateRange = mUEyeOpencvCam.getFrameRateRange();
    CameraFrameRateSlider->setDoubleRange(frameRateRange[0], frameRateRange[1]);
    CameraFrameRateSlider->setDoubleValue(frameRateRange[1]);
    onSetCameraFrameRate(frameRateRange[1]); // set frame-rate to maximum
}

void MainWindow::onSetCameraFrameRate(double value)
{
    cameraFrameRate = mUEyeOpencvCam.setFrameRate(value);
    CameraFrameRateLabel->setText(QString::number(cameraFrameRate, 'f', 1));

    // Set new exposure

    std::vector<double> exposureRange = mUEyeOpencvCam.getExposureRange();
    CameraExposureSlider->setDoubleRange(exposureRange[0], exposureRange[1]);
    CameraExposureSlider->setDoubleValue(exposureRange[1]);
    onSetCameraExposure(exposureRange[1]);
}

void MainWindow::onSetCameraExposure(double value)
{
    mUEyeOpencvCam.setExposure(value);
    CameraExposureLabel->setText(QString::number(value, 'f', 2));
}

void MainWindow::onSetCameraBlackLevelOffset(int value)
{
    double blackLevelOffset = mUEyeOpencvCam.setBlackLevelOffset(value);
    CameraBlackLevelOffsetLabel->setText(QString::number(blackLevelOffset));
}

void MainWindow::onSetCameraBlackLevelMode(int state) { mUEyeOpencvCam.setBlackLevelMode(state); }
void MainWindow::onSetCameraGainBoost     (int state) { mUEyeOpencvCam.setGainBoost(state); }
void MainWindow::onSetCameraAutoGain      (int state) { mUEyeOpencvCam.setAutoGain(state); }

void MainWindow::onSetCameraHardwareGain(int val)
{
    if (!CameraHardwareGainAutoCheckBox->checkState())
    {
        mUEyeOpencvCam.setHardwareGain(val);
    }

    CameraHardwareGainLabel->setText(QString::number(val));
}

void MainWindow::onSetCameraSubSampling(int state)
{
    double subSamplingChange = 0;

    if (state == 1)
    {
        subSamplingChange = 0.5 * cameraSubSamplingFactor;
        cameraSubSamplingFactor = 2;
    }
    else if (state == 0)
    {
        subSamplingChange = cameraSubSamplingFactor;
        cameraSubSamplingFactor = 1;
    }

    cameraAOIWdthMax = cameraAOIWdthMax * subSamplingChange;
    cameraAOIHghtMax = cameraAOIHghtMax * subSamplingChange;

    updateCamAOIx();
    updateCamAOIy();

    if (Parameters::CAMERA_RUNNING)
    {
        if (mUEyeOpencvCam.freeImageMemory())
        {
            if (mUEyeOpencvCam.setSubSampling(cameraSubSamplingFactor))
            {
                if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght);
                    Parameters::CAMERA_READY = true;
                }
            }
        }

        getCameraParameters();
    }
}

void MainWindow::onCropAOI()
{
    int absXPos = Parameters::eyeAOI.xPos + Parameters::cameraAOI.xPos;
    int absYPos = Parameters::eyeAOI.yPos + Parameters::cameraAOI.yPos;

    double fracXPos = absXPos / (double) cameraAOIWdthMax;
    double fracYPos = absYPos / (double) cameraAOIHghtMax;
    double fracWdth = (Parameters::eyeAOI.wdth - cameraAOIWdthMin) / (double) (cameraAOIWdthMax - cameraAOIWdthMin);
    double fracHght = (Parameters::eyeAOI.hght - cameraAOIHghtMin) / (double) (cameraAOIHghtMax - cameraAOIHghtMin);

    Parameters::eyeAOIXPosFraction = 1.0;
    Parameters::eyeAOIYPosFraction = 1.0;
    eyeAOIWdthFraction = 1.0;
    eyeAOIHghtFraction = 1.0;

    CamEyeAOIXPosSlider->setDoubleValue(fracXPos);
    CamEyeAOIYPosSlider->setDoubleValue(fracYPos);
    CamEyeAOIWdthSlider->setDoubleValue(fracWdth);
    CamEyeAOIHghtSlider->setDoubleValue(fracHght);
}

void MainWindow::updateCamAOIx()
{
    Parameters::cameraAOI.wdth = floor(((cameraAOIWdthMax - cameraAOIWdthMin) * cameraAOIFractionWdth + cameraAOIWdthMin) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    Parameters::cameraAOI.xPos = floor((cameraAOIWdthMax * cameraAOIFractionXPos) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    CamEyeAOIXPosSlider->setDoubleMaximum((cameraAOIWdthMax - Parameters::cameraAOI.wdth) / (double) cameraAOIWdthMax);

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        updateEyeAOIx(); }
}

void MainWindow::updateCamAOIy()
{
    Parameters::cameraAOI.hght = floor(((cameraAOIHghtMax - cameraAOIHghtMin) * cameraAOIFractionHght + cameraAOIHghtMin) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    Parameters::cameraAOI.yPos = floor((cameraAOIHghtMax * cameraAOIFractionYPos) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    CamEyeAOIYPosSlider->setDoubleMaximum((cameraAOIHghtMax - Parameters::cameraAOI.hght) / (double) cameraAOIHghtMax);

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        updateEyeAOIy(); }
}

void MainWindow::updateEyeAOIx()
{
    Parameters::eyeAOI.wdth = round(Parameters::cameraAOI.wdth * eyeAOIWdthFraction);
    Parameters::eyeAOI.xPos = round(Parameters::cameraAOI.wdth * Parameters::eyeAOIXPosFraction);

    if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::cameraAOI.wdth)
    {   Parameters::eyeAOI.xPos = Parameters::cameraAOI.wdth - Parameters::eyeAOI.wdth; }

    Parameters::beadAOI.wdth = round(Parameters::cameraAOI.wdth * beadAOIWdthFraction);
    Parameters::beadAOI.xPos = round(Parameters::cameraAOI.wdth * Parameters::beadAOIXPosFraction);

    if (Parameters::beadAOI.xPos + Parameters::beadAOI.wdth > Parameters::cameraAOI.wdth)
    {   Parameters::beadAOI.xPos = Parameters::cameraAOI.wdth - Parameters::beadAOI.wdth; }
}

void MainWindow::updateEyeAOIy()
{
    Parameters::eyeAOI.hght = round(Parameters::cameraAOI.hght * eyeAOIHghtFraction);
    Parameters::eyeAOI.yPos = round(Parameters::cameraAOI.hght * Parameters::eyeAOIYPosFraction);

    if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght > Parameters::cameraAOI.hght)
    {   Parameters::eyeAOI.yPos = Parameters::cameraAOI.hght - Parameters::eyeAOI.hght; }

    Parameters::beadAOI.hght = round(Parameters::cameraAOI.hght * beadAOIHghtFraction);
    Parameters::beadAOI.yPos = round(Parameters::cameraAOI.hght * Parameters::beadAOIYPosFraction);

    if (Parameters::beadAOI.yPos + Parameters::beadAOI.hght > Parameters::cameraAOI.hght)
    {   Parameters::beadAOI.yPos = Parameters::cameraAOI.hght - Parameters::beadAOI.hght; }
}

void MainWindow::onSetCamEyeAOIWdth(double fraction)
{
    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        cameraAOIFractionWdth = fraction;
        Parameters::cameraAOI.wdth = floor(((cameraAOIWdthMax - cameraAOIWdthMin) * cameraAOIFractionWdth + cameraAOIWdthMin) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
        updateEyeAOIx(); }

    if (mUEyeOpencvCam.freeImageMemory())
    {
        if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
        {
            if (Parameters::cameraAOI.xPos + Parameters::cameraAOI.wdth <= cameraAOIWdthMax)
            {
                if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    Parameters::CAMERA_READY = true;
                    getCameraParameters();
                }
            }
        }
    }

    CamEyeAOIXPosSlider->setDoubleMaximum((cameraAOIWdthMax - Parameters::cameraAOI.wdth) / (double) cameraAOIWdthMax);
}

void MainWindow::onSetCamEyeAOIHght(double fraction)
{
    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        cameraAOIFractionHght = fraction;
        Parameters::cameraAOI.hght = floor(((cameraAOIHghtMax - cameraAOIHghtMin) * cameraAOIFractionHght + cameraAOIHghtMin) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
        updateEyeAOIy(); }

    if (mUEyeOpencvCam.freeImageMemory())
    {
        if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
        {
            if (Parameters::cameraAOI.yPos + Parameters::cameraAOI.hght <= cameraAOIHghtMax)
            {
                if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    Parameters::CAMERA_READY = true;
                    getCameraParameters();
                }
            }
        }
    }

    CamEyeAOIYPosSlider->setDoubleMaximum((cameraAOIHghtMax - Parameters::cameraAOI.hght) / (double) cameraAOIHghtMax);
}

void MainWindow::onSetCamEyeAOIXPos(double fraction)
{
    cameraAOIFractionXPos = fraction;

    Parameters::cameraAOI.xPos = floor((cameraAOIWdthMax * cameraAOIFractionXPos) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;

    if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
    {
        Parameters::CAMERA_READY = true;
        getCameraParameters();
    }
}

void MainWindow::onSetCamEyeAOIYPos(double fraction)
{
    cameraAOIFractionYPos = fraction;

    Parameters::cameraAOI.yPos = floor((cameraAOIHghtMax * cameraAOIFractionYPos) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;

    if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
    {
        Parameters::CAMERA_READY = true;
        getCameraParameters();
    }
}

void MainWindow::onSetEyeAOIWdth(double fraction)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    eyeAOIWdthFraction = fraction;

    Parameters::eyeAOI.wdth = round(Parameters::cameraAOI.wdth * eyeAOIWdthFraction);

    if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::cameraAOI.wdth)
    {   Parameters::eyeAOI.xPos = Parameters::cameraAOI.wdth - Parameters::eyeAOI.wdth; }
}

void MainWindow::onSetEyeAOIHght(double fraction)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    eyeAOIHghtFraction = fraction;

    Parameters::eyeAOI.hght = round(Parameters::cameraAOI.hght * eyeAOIHghtFraction);

    if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght > Parameters::cameraAOI.hght)
    {   Parameters::eyeAOI.yPos = Parameters::cameraAOI.hght - Parameters::eyeAOI.hght; }
}

void MainWindow::onSetFlashAOIXPos(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.xPos = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::onSetFlashAOIYPos(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.yPos = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::onSetFlashAOIWdth(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.wdth = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::onSetFlashAOIHght(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.hght = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::onSetFlashThreshold(int val)
{
    flashThreshold = val;
    FlashThresholdLabel->setText(QString::number(flashThreshold));
}

void MainWindow::onSetAOIEyeLeft()
{
    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPosDefaultLeft);
    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPosDefaultLeft);
    CamEyeAOIWdthSlider->setDoubleValue(cameraAOIFractionWdthDefaultLeft);
    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHghtDefaultLeft);
}

void MainWindow::onSetAOIEyeRght()
{
    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPosDefaultRght);
    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPosDefaultRght);
    CamEyeAOIWdthSlider->setDoubleValue(cameraAOIFractionWdthDefaultRght);
    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHghtDefaultRght);
}

void MainWindow::onSetTrialIndex           (int val)   { trialIndex = val; }
void MainWindow::onSetSaveDataAspectRatio  (int state) { SAVE_ASPECT_RATIO  = state; }
void MainWindow::onSetSaveDataCircumference(int state) { SAVE_CIRCUMFERENCE = state; }
void MainWindow::onSetSaveDataPosition     (int state) { SAVE_POSITION      = state; }
void MainWindow::onSetDrawHaar             (int state) { Parameters::drawFlags.haar = state; }
void MainWindow::onSetDrawEdge             (int state) { Parameters::drawFlags.edge = state; }
void MainWindow::onSetDrawElps             (int state) { Parameters::drawFlags.elps = state; }

void MainWindow::onSetCurvatureMeasurement(int state) { mAdvancedOptions.CURVATURE_MEASUREMENT = state; }

void MainWindow::onSetSaveDataEdge (int state) { SAVE_DATA_EDGE  = state; }
void MainWindow::onSetSaveDataFit  (int state) { SAVE_DATA_FIT   = state; }
void MainWindow::onSetSaveDataExtra(int state) { SAVE_DATA_EXTRA = state; }
