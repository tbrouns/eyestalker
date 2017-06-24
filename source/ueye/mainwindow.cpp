//  EyeStalker: robust video-based eye tracking
//  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

//  EyeStalker is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  EyeStalker is distributed in the hope that it will be useful,
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

    Parameters::CAMERA_READY    = false;
    Parameters::CAMERA_RUNNING  = false;
    Parameters::ONLINE_MODE     = true;

    TRIAL_RECORDING = false;
    FLASH_STANDBY   = false;

    SET_FRAME_RATE = true;

    flashThresholdMin = 0;

    Parameters::ellipseDrawCrossSize    = 5;
    Parameters::ellipseDrawOutlineWidth = 0.032;

    cameraPixelClock   = 24;  // MHz
    frameRateOffset    = 1.0; // Hz
    frameCount         = 0;
    guiUpdateFrequency = 30;  // Hz
    relativeTime       = 0;

    startTime         = 0;
    trialIndex        = 0;

    mParameterWidgetEye  = new ParameterWidget;
    mParameterWidgetBead = new ParameterWidget;

    mVariableWidgetEye  = new VariableWidget;
    mVariableWidgetBead = new VariableWidget;

    // AOI

    Parameters::cameraXResolution = 1280;
    Parameters::cameraYResolution = 1024;

    camAOIRatioLeft.hght = 0.19;
    camAOIRatioRght.hght = 0.22;
    camAOIRatioLeft.wdth = 0.25;
    camAOIRatioRght.wdth = 0.31;
    camAOIRatioLeft.xPos = 0.20;
    camAOIRatioRght.xPos = 0.52;
    camAOIRatioLeft.yPos = 0.41;
    camAOIRatioRght.yPos = 0.37;

    camImageHght = 200;
    camImageWdth = 400; // size of image in widget

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

    CamAOIWdthSlider = new SliderDouble;
    CamAOIWdthSlider->setPrecision(2);
    CamAOIWdthSlider->setDoubleRange(0, 1.0);
    CamAOIWdthSlider->setOrientation(Qt::Horizontal);

    CamAOIHghtSlider = new SliderDouble;
    CamAOIHghtSlider->setPrecision(2);
    CamAOIHghtSlider->setDoubleRange(0, 1.0);
    CamAOIHghtSlider->setOrientation(Qt::Vertical);
    CamAOIHghtSlider->setInvertedAppearance(true);

    CamAOIXPosSlider = new SliderDouble;
    CamAOIXPosSlider->setPrecision(2);
    CamAOIXPosSlider->setDoubleRange(0, 1.0);
    CamAOIXPosSlider->setOrientation(Qt::Horizontal);

    CamAOIYPosSlider = new SliderDouble;
    CamAOIYPosSlider->setPrecision(2);
    CamAOIYPosSlider->setDoubleRange(0, 1.0);
    CamAOIYPosSlider->setOrientation(Qt::Vertical);
    CamAOIYPosSlider->setInvertedAppearance(true);

    CameraHardwareGainAutoCheckBox  = new QCheckBox;
    CameraHardwareGainBoostCheckBox = new QCheckBox;

    loadSettings(LastUsedSettingsFileName);

    // Camera feed

    flashAOI    = flashAOILeft;
    camAOITemp  = Parameters::camAOI;

    cv::Mat imgCam(camImageWdth, camImageWdth, CV_8UC3, cv::Scalar(150, 150, 150));

    CamQImage = new QImageOpenCV();
    CamQImage->setSize(camImageWdth, camImageHght);
    CamQImage->setAOIEye  (Parameters::eyeAOI);
    CamQImage->setAOIBead (Parameters::beadAOI);
    CamQImage->setAOIFlash(flashAOI);

    CamQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    CamQImage->loadImage(imgCam);
    CamQImage->setImage();
    QObject::connect(CamQImage, SIGNAL(updateImage(int)), this , SLOT(onUpdateImageRaw(int)));

    // Cam AOI sliders

    CamAOIWdthSlider->setDoubleValue(camAOIRatio.wdth);
    CamAOIHghtSlider->setDoubleValue(camAOIRatio.hght);
    CamAOIXPosSlider->setDoubleValue(camAOIRatio.xPos);
    CamAOIYPosSlider->setDoubleValue(camAOIRatio.yPos);

    QObject::connect(CamAOIWdthSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamAOIWdth(double)));
    QObject::connect(CamAOIHghtSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamAOIHght(double)));
    QObject::connect(CamAOIXPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamAOIXPos(double)));
    QObject::connect(CamAOIYPosSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCamAOIYPos(double)));

    // Eye AOI sliders

    EyeAOIWdthSlider = new SliderDouble;
    EyeAOIWdthSlider->setPrecision(2);
    EyeAOIWdthSlider->setDoubleRange(0, 1.0);
    EyeAOIWdthSlider->setOrientation(Qt::Horizontal);
    QObject::connect(EyeAOIWdthSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetEyeAOIWdth(double)));

    EyeAOIHghtSlider = new SliderDouble;
    EyeAOIHghtSlider->setPrecision(2);
    EyeAOIHghtSlider->setDoubleRange(0, 1.0);
    EyeAOIHghtSlider->setOrientation(Qt::Vertical);
    EyeAOIHghtSlider->setInvertedAppearance(true);
    QObject::connect(EyeAOIHghtSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetEyeAOIHght(double)));

    EyeAOIWdthSlider->setDoubleValue(Parameters::eyeAOIRatio.wdth);
    EyeAOIHghtSlider->setDoubleValue(Parameters::eyeAOIRatio.hght);

    // Bead AOI sliders

    BeadAOIWdthSlider = new SliderDouble;
    BeadAOIWdthSlider->setPrecision(2);
    BeadAOIWdthSlider->setDoubleRange(0, 1.0);
    BeadAOIWdthSlider->setOrientation(Qt::Horizontal);
    QObject::connect(BeadAOIWdthSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetBeadAOIWdth(double)));

    BeadAOIHghtSlider = new SliderDouble;
    BeadAOIHghtSlider->setPrecision(2);
    BeadAOIHghtSlider->setDoubleRange(0, 1.0);
    BeadAOIHghtSlider->setOrientation(Qt::Horizontal);
    QObject::connect(BeadAOIHghtSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetBeadAOIHght(double)));

    BeadAOIWdthSlider->setDoubleValue(Parameters::beadAOIRatio.wdth);
    BeadAOIHghtSlider->setDoubleValue(Parameters::beadAOIRatio.hght);

    //

    QPushButton *AOISetButton = new QPushButton("&Set AOI");
    QObject::connect(AOISetButton,  SIGNAL(clicked()), this, SLOT(onSetCamAOI()));

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

    HaarCheckBox->setChecked(Parameters::drawFlags.haar); // Haar-like features
    EdgeCheckBox->setChecked(Parameters::drawFlags.edge); // Canny edges
    ElpsCheckBox->setChecked(Parameters::drawFlags.elps); // Ellipse fit

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

    QLabel *OnlineProcessingTextBox = new QLabel;
    OnlineProcessingTextBox->setText("Online processing:");

    OnlineProcessingCheckBox = new QCheckBox;
    OnlineProcessingCheckBox->setChecked(!SAVE_EYE_IMAGE);
    QObject::connect(OnlineProcessingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetOnlineProcessing(int)));

    QLabel *OfflineModeTextBox = new QLabel;
    OfflineModeTextBox->setText("Offline mode:");

    QCheckBox *OfflineModeCheckBox = new QCheckBox;
    OfflineModeCheckBox->setChecked(false);
    QObject::connect(OfflineModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetOfflineMode(int)));

    QLabel *BeadDetectionTextBox = new QLabel;
    BeadDetectionTextBox->setText("Bead detection:");

    BeadDetectionCheckBox = new QCheckBox;
    BeadDetectionCheckBox->setChecked(mParameterWidgetBead->getState());
    QObject::connect(BeadDetectionCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetBeadDetection(int)));

    QHBoxLayout *OptionsLayout = new QHBoxLayout;
    OptionsLayout->addStretch();
    OptionsLayout->addWidget(OnlineProcessingTextBox);
    OptionsLayout->addWidget(OnlineProcessingCheckBox);
    OptionsLayout->addWidget(OfflineModeTextBox);
    OptionsLayout->addWidget(OfflineModeCheckBox);
    OptionsLayout->addWidget(BeadDetectionTextBox);
    OptionsLayout->addWidget(BeadDetectionCheckBox);
    OptionsLayout->addStretch();

    QPushButton *AOILeftEyeButton = new QPushButton("&Left eye");
    QObject::connect(AOILeftEyeButton, SIGNAL(clicked()), this, SLOT(onSetAOIEyeLeft()));

    QPushButton *AOIRghtEyeButton = new QPushButton("&Right eye");
    QObject::connect(AOIRghtEyeButton, SIGNAL(clicked()), this, SLOT(onSetAOIEyeRght()));

    AOIEyeOptionsWidget = new QWidget;
    QHBoxLayout* AOIEyeOptionsLayout = new QHBoxLayout(AOIEyeOptionsWidget);
    AOIEyeOptionsLayout->addStretch();
    AOIEyeOptionsLayout->addWidget(AOILeftEyeButton);
    AOIEyeOptionsLayout->addWidget(AOIRghtEyeButton);
    AOIEyeOptionsLayout->addWidget(AOISetButton);
    AOIEyeOptionsLayout->addWidget(AOICropButton);
    AOIEyeOptionsLayout->addStretch();

    QPushButton *ResetPushButton = new QPushButton("Reset parameters");
    QObject::connect(ResetPushButton, SIGNAL(clicked()), this, SLOT(onResetParameters()));

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

    QLabel *OfflineTextBox = new QLabel;
    OfflineTextBox->setText("<b>Detect pupil in:</b>");

    QPushButton *OfflineOneFrameButton  = new QPushButton("One frame");
    QPushButton *OfflineAllFramesButton = new QPushButton("All frames");
    QPushButton *OfflineAllTrialsButton = new QPushButton("All trials");
    QPushButton *OfflineAllExpsButton   = new QPushButton("All experiments");

    QObject::connect(OfflineOneFrameButton,  SIGNAL(clicked(bool)), this, SLOT(onDetectCurrentFrame()));
    QObject::connect(OfflineAllFramesButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllFrames()));
    QObject::connect(OfflineAllTrialsButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllTrials()));
    QObject::connect(OfflineAllExpsButton,   SIGNAL(clicked(bool)), this, SLOT(onDetectAllExperiments()));

    QPushButton *SavePupilDataButton = new QPushButton("Save");
    QObject::connect(SavePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(onSaveTrialData()));

    QPushButton *CombinePupilDataButton = new QPushButton("Combine");
    QObject::connect(CombinePupilDataButton, SIGNAL(clicked(bool)), this, SLOT(onCombineData()));

    OfflineModeWidget = new QWidget;
    QVBoxLayout *EyeTrackingOfflineLayout = new QVBoxLayout(OfflineModeWidget);

    QHBoxLayout *OfflineOptionsLayout = new QHBoxLayout;
    OfflineOptionsLayout->addWidget(OfflineLoadSessionButton);
    OfflineOptionsLayout->addWidget(SavePupilDataButton);
    OfflineOptionsLayout->addWidget(CombinePupilDataButton);
    OfflineOptionsLayout->addWidget(OfflinePrevImageButton);
    OfflineOptionsLayout->addWidget(OfflineImageSlider);
    OfflineOptionsLayout->addWidget(OfflineNextImageButton);
    OfflineOptionsLayout->addWidget(OfflineImageFrameTextBox);

    QHBoxLayout *OfflineDetectionLayout = new QHBoxLayout;
    OfflineDetectionLayout->addWidget(OfflineTextBox);
    OfflineDetectionLayout->addWidget(OfflineOneFrameButton);
    OfflineDetectionLayout->addWidget(OfflineAllFramesButton);
    OfflineDetectionLayout->addWidget(OfflineAllTrialsButton);
    OfflineDetectionLayout->addWidget(OfflineAllExpsButton);

    QLabel* OfflineTrialTitle = new QLabel;
    OfflineTrialTitle->setText("<b>Offline mode - Trial:</b>");

    OfflineTrialSpinBox = new QSpinBox;
    OfflineTrialSpinBox->setValue(0);
    OfflineTrialSpinBox->setAlignment(Qt::AlignRight);

    OfflineTrialSlider = new QSlider;
    OfflineTrialSlider->setValue(0);
    OfflineTrialSlider->setOrientation(Qt::Horizontal);

    QObject::connect(OfflineTrialSlider,  SIGNAL(valueChanged(int)), this,                SLOT(onSetTrialOffline(int)));
    QObject::connect(OfflineTrialSlider,  SIGNAL(valueChanged(int)), OfflineTrialSpinBox, SLOT(setValue(int)));
    QObject::connect(OfflineTrialSpinBox, SIGNAL(valueChanged(int)), OfflineTrialSlider,  SLOT(setValue(int)));

    QGridLayout *OfflineSessionTitleLayout = new QGridLayout;
    OfflineSessionTitleLayout->addWidget(OfflineTrialTitle,     0, 1);
    OfflineSessionTitleLayout->addWidget(OfflineTrialSpinBox,   0, 2);
    OfflineSessionTitleLayout->addWidget(OfflineTrialSlider,    0, 3);
    OfflineSessionTitleLayout->setColumnStretch(0, 1);
    OfflineSessionTitleLayout->setColumnStretch(1, 3);
    OfflineSessionTitleLayout->setColumnStretch(2, 1);
    OfflineSessionTitleLayout->setColumnStretch(3, 3);
    OfflineSessionTitleLayout->setColumnStretch(4, 1);

    EyeTrackingOfflineLayout->addLayout(OfflineSessionTitleLayout);
    EyeTrackingOfflineLayout->addLayout(OfflineOptionsLayout);
    EyeTrackingOfflineLayout->addLayout(OfflineDetectionLayout);

    OfflineModeWidget->setVisible(false);

    QWidget *CameraSettings = new QWidget;
    QGridLayout* CameraOutputLayout = new QGridLayout(CameraSettings);

    CameraOutputLayout->addLayout(OptionsLayout,       0, 2);
    CameraOutputLayout->addWidget(CamAOIXPosSlider,    1, 2);
    CameraOutputLayout->addWidget(CamAOIYPosSlider,    2, 1);
    CameraOutputLayout->addWidget(CamQImage,           2, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(CamAOIWdthSlider,    3, 2);
    CameraOutputLayout->addWidget(EyeAOIWdthSlider,    4, 2);
    CameraOutputLayout->addWidget(CamAOIHghtSlider,    2, 3);
    CameraOutputLayout->addWidget(EyeAOIHghtSlider,    2, 4);
    CameraOutputLayout->addWidget(mQwtPlotWidget,      2, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(OfflineModeWidget,   3, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(AOIEyeOptionsWidget, 5, 2);
    CameraOutputLayout->addLayout(DrawFunctionsLayout, 6, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(ResetPushButton,     7, 2, Qt::AlignCenter);

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
    QObject::connect(CameraFrameRateDesiredSpinBox, SIGNAL(valueChanged(int)), this, SLOT(onSetCameraFrameRateDesired(int)));

    QPushButton* CameraFrameRateButton = new QPushButton("Set");
    QObject::connect(CameraFrameRateButton, SIGNAL(clicked(bool)), this, SLOT(onCalibrateFrameRate()));

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
    CameraFrameRateDesiredLayout->addWidget(CameraFrameRateButton);
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

    //    CameraParametersLayout->addWidget(CameraSubSamplingTextBox,  8, 0);
    //    CameraParametersLayout->addWidget(CameraSubSamplingCheckBox, 8, 1);

    ///////////////////////////////////////////////////////////////
    ///////////////////// EXPERIMENT TAB  /////////////////////////
    ///////////////////////////////////////////////////////////////

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

    // Left

    QSpinBox* FlashXPosSpinBoxLeft = new QSpinBox;
    QSpinBox* FlashYPosSpinBoxLeft = new QSpinBox;
    QSpinBox* FlashWdthSpinBoxLeft = new QSpinBox;
    QSpinBox* FlashHghtSpinBoxLeft = new QSpinBox;

    FlashXPosSpinBoxLeft->setRange(0, Parameters::cameraXResolution);
    FlashYPosSpinBoxLeft->setRange(0, Parameters::cameraYResolution);
    FlashWdthSpinBoxLeft->setRange(0, Parameters::cameraXResolution);
    FlashHghtSpinBoxLeft->setRange(0, Parameters::cameraYResolution);

    FlashXPosSpinBoxLeft->setValue(flashAOILeft.xPos);
    FlashYPosSpinBoxLeft->setValue(flashAOILeft.yPos);
    FlashWdthSpinBoxLeft->setValue(flashAOILeft.wdth);
    FlashHghtSpinBoxLeft->setValue(flashAOILeft.hght);

    QObject::connect(FlashXPosSpinBoxLeft, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIXPosLeft(int)));
    QObject::connect(FlashYPosSpinBoxLeft, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIYPosLeft(int)));
    QObject::connect(FlashWdthSpinBoxLeft, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIWdthLeft(int)));
    QObject::connect(FlashHghtSpinBoxLeft, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIHghtLeft(int)));

    QLabel* FlashCoordsTextBoxLeft = new QLabel;
    FlashCoordsTextBoxLeft->setText("<b>Flash AOI left</b>");

    QLabel* FlashXPosTextBoxLeft = new QLabel;
    FlashXPosTextBoxLeft->setText("<i>X:</i>");

    QLabel* FlashYPosTextBoxLeft = new QLabel;
    FlashYPosTextBoxLeft->setText("<i>Y:</i>");

    QLabel* FlashWdthTextBoxLeft = new QLabel;
    FlashWdthTextBoxLeft->setText("<i>W:</i>");

    QLabel* FlashHghtTextBoxLeft = new QLabel;
    FlashHghtTextBoxLeft->setText("<i>H:</i>");

    QHBoxLayout* FlashCoordinatesLayoutLeft = new QHBoxLayout;
    FlashCoordinatesLayoutLeft->addWidget(FlashXPosTextBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashXPosSpinBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashYPosTextBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashYPosSpinBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashWdthTextBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashWdthSpinBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashHghtTextBoxLeft);
    FlashCoordinatesLayoutLeft->addWidget(FlashHghtSpinBoxLeft);

    // Right

    QSpinBox* FlashXPosSpinBoxRght = new QSpinBox;
    QSpinBox* FlashYPosSpinBoxRght = new QSpinBox;
    QSpinBox* FlashWdthSpinBoxRght = new QSpinBox;
    QSpinBox* FlashHghtSpinBoxRght = new QSpinBox;

    FlashXPosSpinBoxRght->setRange(0, Parameters::cameraXResolution);
    FlashYPosSpinBoxRght->setRange(0, Parameters::cameraYResolution);
    FlashWdthSpinBoxRght->setRange(0, Parameters::cameraXResolution);
    FlashHghtSpinBoxRght->setRange(0, Parameters::cameraYResolution);

    FlashXPosSpinBoxRght->setValue(flashAOIRght.xPos);
    FlashYPosSpinBoxRght->setValue(flashAOIRght.yPos);
    FlashWdthSpinBoxRght->setValue(flashAOIRght.wdth);
    FlashHghtSpinBoxRght->setValue(flashAOIRght.hght);

    QObject::connect(FlashXPosSpinBoxRght, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIXPosRght(int)));
    QObject::connect(FlashYPosSpinBoxRght, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIYPosRght(int)));
    QObject::connect(FlashWdthSpinBoxRght, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIWdthRght(int)));
    QObject::connect(FlashHghtSpinBoxRght, SIGNAL(valueChanged(int)), this, SLOT(onSetFlashAOIHghtRght(int)));

    QLabel* FlashCoordsTextBoxRght = new QLabel;
    FlashCoordsTextBoxRght->setText("<b>Flash AOI right</b>");

    QLabel* FlashXPosTextBoxRght = new QLabel;
    FlashXPosTextBoxRght->setText("<i>X:</i>");

    QLabel* FlashYPosTextBoxRght = new QLabel;
    FlashYPosTextBoxRght->setText("<i>Y:</i>");

    QLabel* FlashWdthTextBoxRght = new QLabel;
    FlashWdthTextBoxRght->setText("<i>W:</i>");

    QLabel* FlashHghtTextBoxRght = new QLabel;
    FlashHghtTextBoxRght->setText("<i>H:</i>");

    QHBoxLayout* FlashCoordinatesLayoutRght = new QHBoxLayout;
    FlashCoordinatesLayoutRght->addWidget(FlashXPosTextBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashXPosSpinBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashYPosTextBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashYPosSpinBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashWdthTextBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashWdthSpinBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashHghtTextBoxRght);
    FlashCoordinatesLayoutRght->addWidget(FlashHghtSpinBoxRght);

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
    SaveDataTextBox->setText("<b>Save data:</b>");

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

    QHBoxLayout* ExperimentSaveLayout = new QHBoxLayout;
    ExperimentSaveLayout->addWidget(SaveDataPositionTextBox);
    ExperimentSaveLayout->addWidget(SaveDataPositionCheckBox);
    ExperimentSaveLayout->addWidget(SaveDataCircumferenceTextBox);
    ExperimentSaveLayout->addWidget(SaveDataCircumferenceCheckBox);
    ExperimentSaveLayout->addWidget(SaveDataAspectRatioTextBox);
    ExperimentSaveLayout->addWidget(SaveDataAspectRatioCheckBox);

    // Set-up layouts

    QWidget* ExperimentTabWidget     = new QWidget;
    QGridLayout *ExperimentTabLayout = new QGridLayout(ExperimentTabWidget);

    ExperimentTabLayout->addWidget(DataDirectoryTitleTextBox,   1, 0);
    ExperimentTabLayout->addWidget(DataDirectoryTextBox,        1, 1);
    ExperimentTabLayout->addWidget(DataDirectoryButton,         1, 2);
    ExperimentTabLayout->addWidget(DataFilenameTextBox,         2, 0);
    ExperimentTabLayout->addWidget(DataFilenameLineEdit,        2, 1);
    ExperimentTabLayout->addWidget(TrialIndexTextBox,           3, 0);
    ExperimentTabLayout->addLayout(TrialIndexSpinBoxLayout,     3, 1);
    ExperimentTabLayout->addWidget(TrialTimeLengthTextBox,      4, 0);
    ExperimentTabLayout->addWidget(TrialTimeLengthLineEdit,     4, 1);
    ExperimentTabLayout->addWidget(StartRecordingButton,        4, 2);
    ExperimentTabLayout->addWidget(FlashStandbyTextBox,         5, 0);
    ExperimentTabLayout->addWidget(FlashStandbySlider,          5, 1);
    ExperimentTabLayout->addWidget(FlashStandbyLabel,           5, 2);
    ExperimentTabLayout->addWidget(FlashThresholdTextBox,       6, 0);
    ExperimentTabLayout->addWidget(FlashThresholdSlider,        6, 1);
    ExperimentTabLayout->addWidget(FlashThresholdLabel,         6, 2);
    ExperimentTabLayout->addWidget(FlashThresholdResetButton,   6, 3, Qt::AlignLeft);
    ExperimentTabLayout->addWidget(FlashCoordsTextBoxLeft,      7, 0);
    ExperimentTabLayout->addLayout(FlashCoordinatesLayoutLeft,  7, 1, 1, 3);
    ExperimentTabLayout->addWidget(FlashCoordsTextBoxRght,      8, 0);
    ExperimentTabLayout->addLayout(FlashCoordinatesLayoutRght,  8, 1, 1, 3);
    ExperimentTabLayout->addWidget(SaveDataTextBox,             9, 0);
    ExperimentTabLayout->addLayout(ExperimentSaveLayout,        9, 1, 1, 3);

    ExperimentTabLayout->setColumnStretch(0, 1);
    ExperimentTabLayout->setColumnStretch(1, 2);
    ExperimentTabLayout->setColumnStretch(2, 1);
    ExperimentTabLayout->setColumnStretch(3, 1);
    ExperimentTabLayout->setColumnStretch(4, 1);
    ExperimentTabLayout->setColumnStretch(5, 5);

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

    QLabel* BeadAOIWdthTextBox = new QLabel;
    BeadAOIWdthTextBox->setText("<b>Bead AOI width:</b>");

    QLabel* BeadAOIHghtTextBox = new QLabel;
    BeadAOIHghtTextBox->setText("<b>Bead AOI height:</b>");

    QLabel* BeadAOITitleTextBox = new QLabel;
    BeadAOITitleTextBox->setText("<b>Bead AOI dimensions</b>");

    QGridLayout* BeadAOILayout = new QGridLayout;
    BeadAOILayout->addWidget(BeadAOITitleTextBox, 0, 1, 1, 1, Qt::AlignCenter);
    BeadAOILayout->addWidget(BeadAOIWdthTextBox,  1, 0, 1, 1, Qt::AlignRight);
    BeadAOILayout->addWidget(BeadAOIWdthSlider,   1, 1);
    BeadAOILayout->addWidget(BeadAOIHghtTextBox,  2, 0, 1, 1, Qt::AlignRight);
    BeadAOILayout->addWidget(BeadAOIHghtSlider,   2, 1);

    BeadAOILayout->setColumnStretch(0,1);
    BeadAOILayout->setColumnStretch(1,3);
    BeadAOILayout->setColumnStretch(2,1);

    QWidget* BeadTrackingWidget = new QWidget;
    QVBoxLayout *BeadTrackingLayout = new QVBoxLayout(BeadTrackingWidget);
    BeadTrackingLayout->addLayout(BeadAOILayout);
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
    MainTabWidget->addTab(EyeTrackingScrollArea,   tr("Eye tracking"));
    MainTabWidget->addTab(BeadTrackingScrollArea,  tr("Bead tracking"));
    MainTabWidget->addTab(AdvancedScrollArea,      tr("Advanced"));

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
    { std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
        resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
    }

    { std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
        resetVariablesHard(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);
    }

    while(APP_RUNNING && Parameters::CAMERA_RUNNING && Parameters::ONLINE_MODE)
    {
        std::lock_guard<std::mutex> AOILock_1(mutexAOI_1);

        detectionVariables mDetectionVariablesEyeTemp;
        detectionParameters mDetectionParametersEyeTemp;

        dataVariables mDataVariablesEyeTemp;
        drawVariables mDrawVariablesEyeTemp;

        detectionVariables mDetectionVariablesBeadTemp;
        detectionParameters mDetectionParametersBeadTemp;

        dataVariables mDataVariablesBeadTemp;
        drawVariables mDrawVariablesBeadTemp;

        AOIProperties AOICameraTemp;
        AOIProperties AOIFlashTemp;
        AOIProperties AOIBeadTemp;
        AOIProperties AOIEyeTemp;

        imageInfo mImageInfo = mUEyeOpencvCam.getFrame(); // get new frame from camera
        cv::Mat imageOriginal = mImageInfo.image;
        absoluteTime = mImageInfo.time; // Get frame timestamp

        double relativeTimeNew = (absoluteTime - startTime) / (double) 10000; // in ms

        // ignore frame if time interval was too short (possible camera error)
        if (relativeTimeNew <= (relativeTime + 0.9 * (1000 / cameraFrameRate))) { continue; }

        relativeTime = relativeTimeNew;

        { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
            imageCamera = imageOriginal.clone();

            mDetectionVariablesEyeTemp  = mDetectionVariablesEye;
            mDetectionParametersEyeTemp = mParameterWidgetEye->getStructure();

            mDetectionVariablesBeadTemp  = mDetectionVariablesBead;
            mDetectionParametersBeadTemp = mParameterWidgetBead->getStructure();

            AOIFlashTemp  = flashAOI;
            AOICameraTemp = Parameters::camAOI;
        }

        { std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
            AOIEyeTemp = Parameters::eyeAOI;
        }

        { std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
            AOIBeadTemp = Parameters::beadAOI;
        }

        AOIProperties AOIFlashRelative;

        bool FLASH_AOI_VISIBLE = false;
        if (!TRIAL_RECORDING) { FLASH_AOI_VISIBLE = checkFlashAOI(AOIFlashRelative, AOIFlashTemp, AOICameraTemp); }

        // Check limits

        int imgWdth = imageOriginal.cols;
        int imgHght = imageOriginal.rows;

        if
                (imgWdth         >= eyeAOIWdthMin &&
                 imgHght         >= eyeAOIHghtMin &&
                 AOIEyeTemp.wdth >= eyeAOIWdthMin &&
                 AOIEyeTemp.hght >= eyeAOIHghtMin)
        {
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
                        resetVariablesSoft(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), AOIEyeTemp);
                        resetVariablesSoft(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), AOIBeadTemp);
                        startTime = mImageInfo.time;
                        startTrialRecording();
                    }
                }
                else // Default mode
                {
                    if (avgIntensity > flashThreshold)
                    {
                        FlashThresholdSlider->setMinimum(floor(avgIntensity));
                        FlashThresholdSlider->setValue(  floor(avgIntensity));
                    }

                    mDetectionVariablesEyeTemp = eyeStalker(imageOriginal, AOIEyeTemp, mDetectionVariablesEyeTemp, mDetectionParametersEyeTemp, mDataVariablesEyeTemp, mDrawVariablesEyeTemp); // Pupil tracking algorithm

                    if (mDetectionParametersBeadTemp.DETECTION_ON) // bead detection
                    {
                        mDetectionVariablesBeadTemp = eyeStalker(imageOriginal, AOIBeadTemp, mDetectionVariablesBeadTemp, mDetectionParametersBeadTemp, mDataVariablesBeadTemp, mDrawVariablesBeadTemp);
                    }
                }
            }
            else // Trial recording
            {
                if (!SAVE_EYE_IMAGE)
                {
                    mDetectionVariablesEyeTemp         = eyeStalker(imageOriginal, AOIEyeTemp, mDetectionVariablesEyeTemp, mDetectionParametersEyeTemp, mDataVariablesEyeTemp, mDrawVariablesEyeTemp); // Pupil tracking algorithm
                    mDataVariablesEyeTemp.absoluteXPos = mDataVariablesEyeTemp.exactXPos + AOIEyeTemp.xPos + AOICameraTemp.xPos;
                    mDataVariablesEyeTemp.absoluteYPos = mDataVariablesEyeTemp.exactYPos + AOIEyeTemp.yPos + AOICameraTemp.yPos;
                    vDataVariablesEye[frameCount]      = mDataVariablesEyeTemp;

                    if (mDetectionParametersBeadTemp.DETECTION_ON) // bead detection
                    {
                        mDetectionVariablesBeadTemp         = eyeStalker(imageOriginal, AOIBeadTemp, mDetectionVariablesBeadTemp, mDetectionParametersBeadTemp, mDataVariablesBeadTemp, mDrawVariablesBeadTemp); // Pupil tracking algorithm
                        mDataVariablesBeadTemp.absoluteXPos = mDataVariablesBeadTemp.exactXPos + AOIBeadTemp.xPos + AOICameraTemp.xPos;
                        mDataVariablesBeadTemp.absoluteYPos = mDataVariablesBeadTemp.exactYPos + AOIBeadTemp.yPos + AOICameraTemp.yPos;
                        vDataVariablesBead[frameCount]      = mDataVariablesBeadTemp;
                    }

                    mDataVariablesEyeTemp.timestamp    = relativeTime; // save time stamps

                    frameCount++;
                }
                else
                {
                    vDataVariablesEye[frameCount].timestamp = relativeTime; // save time stamps

                    // Saving camera frame

                    std::stringstream filename;
                    filename << dataDirectory
                             << "/"
                             << currentDate
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
                    Parameters::frameCaptureCV.notify_all(); // continue regular frame capture
                    if (!SAVE_EYE_IMAGE) { emit showPlot(); }
                    emit startTimer(round(1000 / guiUpdateFrequency));
                }
            }
        }

        // Update structures

        {
            std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

            mDetectionVariablesEye = mDetectionVariablesEyeTemp;
            mDrawVariablesEye      = mDrawVariablesEyeTemp;
            mDataVariablesEye      = mDataVariablesEyeTemp;
            mDrawVariablesBead     = mDrawVariablesBeadTemp;
            mDataVariablesBead     = mDataVariablesBeadTemp;
        }
    }

    if (!APP_RUNNING)
    {
        std::unique_lock<std::mutex> mtxLock(mutexQuit);
        APP_EXIT = true;
        cv.notify_all();
    }
    else if (!Parameters::ONLINE_MODE)
    {
        std::unique_lock<std::mutex> mtxLock(mutexQuit);
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
    if (Parameters::ONLINE_MODE)
    {
        if (!TRIAL_RECORDING)
        {
            if (Parameters::CAMERA_RUNNING)
            {
                std::lock_guard<std::mutex> AOILock_2(mutexAOI_2);

                drawVariables mDrawVariablesEyeTemp;
                dataVariables mDataVariablesEyeTemp;

                drawVariables mDrawVariablesBeadTemp;
                dataVariables mDataVariablesBeadTemp;

                cv::Mat imageOriginal;

                int imgWdth = 0;
                int imgHght = 0;

                AOIProperties AOIEyeTemp;
                AOIProperties AOIBeadTemp;
                AOIProperties AOIFlashTemp;

                bool DRAW_BEAD = false;

                { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
                    if (!imageCamera.empty())
                    {
                        imageOriginal = imageCamera.clone();

                        imgWdth = imageOriginal.cols;
                        imgHght = imageOriginal.rows;

                        mDrawVariablesEyeTemp = mDrawVariablesEye;
                        mDataVariablesEyeTemp = mDataVariablesEye;

                        mDrawVariablesBeadTemp = mDrawVariablesBead;
                        mDataVariablesBeadTemp = mDataVariablesBead;

                        detectionParameters mDetectionParameters = mParameterWidgetBead->getStructure();
                        DRAW_BEAD = mDetectionParameters.DETECTION_ON;

                        AOIFlashTemp       = flashAOI;

                    } else { return; }
                }

                { std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
                    AOIEyeTemp = Parameters::eyeAOI;
                }

                { std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
                    AOIBeadTemp  = Parameters::beadAOI;
                }

                if
                        (imgWdth         >= eyeAOIWdthMin &&
                         imgHght         >= eyeAOIHghtMin &&
                         AOIEyeTemp.wdth >= eyeAOIWdthMin &&
                         AOIEyeTemp.hght >= eyeAOIHghtMin)
                {

                    { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

                        // Increase pixel clock if desired frame-rate has not been reached

                        if (SET_FRAME_RATE)
                        {
                            int desiredFrameRate = CameraFrameRateDesiredSpinBox->value();

                            if (desiredFrameRate > cameraFrameRate + frameRateOffset)
                            {
                                if (cameraPixelClock < CameraPixelClockSlider->maximum())
                                {
                                    cameraPixelClock = cameraPixelClock + 1;
                                    CameraPixelClockSlider->setValue(cameraPixelClock);

                                    CameraFrameRateSlider->setDoubleValue(desiredFrameRate);
                                }
                                else
                                {
                                    CameraFrameRateSlider->setDoubleValue(CameraFrameRateSlider->maximum());
                                    SET_FRAME_RATE = false; // frame-rate is accepted
                                }
                            }

                            if (desiredFrameRate < cameraFrameRate - frameRateOffset) // frame-rate should be within 1 Hz
                            {
                                if (cameraPixelClock > CameraPixelClockSlider->minimum())
                                {
                                    cameraPixelClock = cameraPixelClock - 1;
                                    CameraPixelClockSlider->setValue(cameraPixelClock);
                                }
                                else
                                {
                                    CameraFrameRateSlider->setDoubleValue(desiredFrameRate);
                                }
                            }

                            if (desiredFrameRate >= cameraFrameRate - frameRateOffset && desiredFrameRate <= cameraFrameRate + frameRateOffset)
                            {
                                SET_FRAME_RATE = false; // frame-rate is accepted
                            }
                        }

                        // For flash

                        if (CameraHardwareGainAutoCheckBox->checkState())
                        {
                            int hardwareGain = mUEyeOpencvCam.getHardwareGain();
                            CameraHardwareGainSlider->setValue(hardwareGain);
                            CameraHardwareGainLabel->setText(QString::number(hardwareGain));
                        }
                    }

                    mVariableWidgetEye->setWidgets(mDataVariablesEyeTemp); // update sliders

                    cv::Mat imageProcessed = imageOriginal.clone();
                    drawAll(imageProcessed, mDrawVariablesEyeTemp);  // draw eye features
                    if (DRAW_BEAD) { drawAll(imageProcessed, mDrawVariablesBeadTemp); } // draw bead features
                    CamQImage->loadImage(imageProcessed);
                    CamQImage->setAOIEye  (  AOIEyeTemp);
                    CamQImage->setAOIBead ( AOIBeadTemp);
                    CamQImage->setAOIFlash(AOIFlashTemp);
                    CamQImage->setImage();
                }
                else { CamQImage->setAOIError(); }
            }
            else if (!Parameters::CAMERA_RUNNING)
            {
                CamQImage->setSpinner();
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

    while (APP_RUNNING && !Parameters::CAMERA_RUNNING && Parameters::ONLINE_MODE)
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
                if (mUEyeOpencvCam.allocateMemory(Parameters::camAOI.wdth, Parameters::camAOI.hght))
                {
                    if (mUEyeOpencvCam.setAOI(Parameters::camAOI.xPos, Parameters::camAOI.yPos, Parameters::camAOI.wdth, Parameters::camAOI.hght))
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
    else { std::this_thread::sleep_for(std::chrono::milliseconds(2000)); }
}

// Variables

void MainWindow::resetVariablesHard(detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const AOIProperties& mAOI)
{
    // Reset all variables

    mDetectionVariables.averageAspectRatio   = initialAspectRatio; // close to perfect circle
    mDetectionVariables.averageCircumference = 0.5 * (mDetectionParameters.thresholdCircumferenceMax + mDetectionParameters.thresholdCircumferenceMin); // calculate first
    mDetectionVariables.averageCurvature     = initialCurvature;
    mDetectionVariables.averageGradient      = 0;
    mDetectionVariables.averageHeight        = mDetectionVariables.averageCircumference / M_PI;
    mDetectionVariables.averageIntensity     = initialIntensity;
    mDetectionVariables.averageHaarResponse  = 0;
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
    mDetectionVariables.predictedHaarResponse  = mDetectionVariables.averageHaarResponse;
    mDetectionVariables.predictedHeight        = mDetectionVariables.averageHeight;
    mDetectionVariables.predictedIntensity     = mDetectionVariables.averageIntensity;
    mDetectionVariables.predictedWidth         = mDetectionVariables.averageWidth;
    mDetectionVariables.predictedXPos          = 0.5 * (mAOI.wdth - 1); // centre of image
    mDetectionVariables.predictedYPos          = 0.5 * (mAOI.hght - 1);

    mDetectionVariables.momentumAspectRatio   = 0;
    mDetectionVariables.momentumCircumference = 0;
    mDetectionVariables.momentumCurvature     = 0;
    mDetectionVariables.momentumGradient      = 0;
    mDetectionVariables.momentumHaarResponse  = 0;
    mDetectionVariables.momentumHeight        = 0;
    mDetectionVariables.momentumIntensity     = 0;
    mDetectionVariables.momentumWidth         = 0;
    mDetectionVariables.momentumXPos          = 0;
    mDetectionVariables.momentumYPos          = 0;

    double maxChangeThresholdAspectRatio_1 = mDetectionVariables.predictedAspectRatio - mDetectionParameters.thresholdAspectRatioMin;
    double maxChangeThresholdAspectRatio_2 = 1.0 - mDetectionVariables.predictedAspectRatio;
    double maxChangeThresholdAspectRatio   = std::max(maxChangeThresholdAspectRatio_1, maxChangeThresholdAspectRatio_2);
    double rangeChangeThresholdAspectRatio = maxChangeThresholdAspectRatio   - mDetectionParameters.thresholdChangeAspectRatioUpper;
    mDetectionVariables.thresholdChangeAspectRatioUpper   = rangeChangeThresholdAspectRatio   + mDetectionParameters.thresholdChangeAspectRatioUpper;

    double maxChangeThresholdCircumference_1 = (mDetectionParameters.thresholdCircumferenceMax - mDetectionVariables.predictedCircumference) / mDetectionParameters.thresholdCircumferenceMax;
    double maxChangeThresholdCircumference_2 = (mDetectionVariables.predictedCircumference - mDetectionParameters.thresholdCircumferenceMin) / mDetectionVariables.predictedCircumference;
    double maxChangeThresholdCircumference   = std::max(maxChangeThresholdCircumference_1, maxChangeThresholdCircumference_2);
    double rangeChangeThresholdCircumference = maxChangeThresholdCircumference - mDetectionParameters.thresholdChangeCircumferenceUpper;
    mDetectionVariables.thresholdChangeCircumferenceUpper = rangeChangeThresholdCircumference + mDetectionParameters.thresholdChangeCircumferenceUpper;

    double maxChangeThresholdPositionX = mAOI.wdth - mDetectionVariables.predictedWidth;
    double maxChangeThresholdPositionY = mAOI.hght - mDetectionVariables.predictedHeight;
    double maxChangeThresholdPosition  = std::max(maxChangeThresholdPositionX, maxChangeThresholdPositionY);
    double rangeChangeThresholdPosition = maxChangeThresholdPosition - mDetectionParameters.thresholdChangePositionUpper;
    mDetectionVariables.thresholdChangePositionUpper = rangeChangeThresholdPosition + mDetectionParameters.thresholdChangePositionUpper;

    mDetectionVariables.thresholdScoreEdge = 0;
    mDetectionVariables.thresholdScoreFit  = 0;

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
        std::unique_lock<std::mutex> lck(mutexQuit);
        while (!APP_EXIT) cv.wait(lck);
    }

    mUEyeOpencvCam.exitCamera();

    saveSettings(LastUsedSettingsFileName);

    qApp->quit();
}

void MainWindow::onSetOfflineMode(int state)
{
    MainTabWidget->setTabEnabled(0, state);

    CamAOIWdthSlider->setVisible(!state);
    CamAOIHghtSlider->setVisible(!state);
    CamAOIXPosSlider->setVisible(!state);
    CamAOIYPosSlider->setVisible(!state);
    EyeAOIWdthSlider->setVisible(!state);
    EyeAOIHghtSlider->setVisible(!state);

    AOIEyeOptionsWidget->setVisible(!state);
    OfflineModeWidget  ->setVisible(state);

    Parameters::ONLINE_MODE = !state;
    Parameters::CAMERA_RUNNING    = !state;
    Parameters::CAMERA_READY      = !state;

    trialIndexOffline = 0;
    imageIndexOffline = 0;

    if (!state)
    {
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
        { // unlock frame grabbing threads
            std::unique_lock<std::mutex> frameCaptureMutexLock(Parameters::frameCaptureMutex);
            Parameters::frameCaptureCV.notify_all(); // unlock getFrame() thread
        }

        { // wait for threads to finish
            std::unique_lock<std::mutex> lck(mutexQuit);
            while (Parameters::CAMERA_RUNNING) cv.wait(lck);
        }

        CamQImage->clearImage();
        setupOfflineSession();
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
    QString text = "EyeStalker: robust video-based eye tracking <br> <br>"
                   "Copyright 2016 - 2017 Terence Brouns, t.s.n.brouns@gmail.com <br> <br> "
                   "EyeStalker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; <br> "
                   "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. <br><br>"
                   "See the GNU General Public License for more details. <br><br>"
                   "Included 3rd party libraries: <br><br>"
                   "<b>Qt:</b> <br><br>"
                   "(c) 2017 The Qt Company Ltd. https://www.qt.io/ <br><br>"
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
                   "<b>UEye API:</b> <br><br> (c) 2016, IDS Imaging Advanced Systems GmbH, https://en.ids-imaging.com/<br><br>"
                   "<b>libusb:</b> <br><br> libusb is released under version 2.1 of the GNU Lesser General Public License (LGPL).<br>"
                   "http://libusb.info/ <br><br>"
                   "<b>QDarkStyleSheet:</b> <br><br>"
                   "Copyright (c) 2013-2014, Colin Duquesnoy. Released under the MIT license. <br>";

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
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

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

            vDataVariablesEye.resize(trialFrameTotal);
            vDataVariablesBead.resize(trialFrameTotal);

            vDetectionVariablesEye.resize(trialFrameTotal);
            vDetectionVariablesBead.resize(trialFrameTotal);

            // Reset variables

            if (SAVE_EYE_IMAGE) // create directories if necessary
            {
                std::stringstream directoryName;
                directoryName << dataDirectory << "/" << currentDate;

                if (!boost::filesystem::exists(directoryName.str()))
                {    boost::filesystem::create_directory(directoryName.str().c_str()); }

                directoryName << "/trial_" << trialIndex;
                if (!boost::filesystem::exists(directoryName.str()))
                {    boost::filesystem::create_directory(directoryName.str().c_str()); }

                directoryName << "/raw/";
                if (!boost::filesystem::exists(directoryName.str()))
                {    boost::filesystem::create_directory(directoryName.str().c_str()); }
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
                 << dataFilename
                 << ".dat";

        std::ofstream file;

        if (!boost::filesystem::exists(filename.str()))
        {   file.open(filename.str(), std::ios::out | std::ios::ate); }
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

        for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].DETECTED  << delimiter; }
        for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].timestamp << delimiter; } // saving ALL timestamps allows for accurate frame-rate calculation during data processing

        if (SAVE_POSITION)
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].absoluteXPos  << delimiter; }
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].absoluteYPos  << delimiter; }
        }

        if (SAVE_CIRCUMFERENCE)
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].exactCircumference << delimiter; }
        }

        if (SAVE_ASPECT_RATIO)
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].exactAspectRatio << delimiter; }
        }

        if (mParameterWidgetBead->getState())
        {
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesBead[i].DETECTED      << delimiter; }
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesBead[i].absoluteXPos  << delimiter; }
            for (int i = 0; i < frameCount; i++) { file << vDataVariablesBead[i].absoluteYPos  << delimiter; }
        }

        file.close();
    }
    else
    {
        // save time stamps

        std::stringstream filename;
        filename << dataDirectory << "/"
                 << currentDate   << "/"
                 << "trial_"      << trialIndex
                 << "/"
                 << "timestamps.dat";

        std::ofstream file;

        file.open(filename.str(), std::ios::out | std::ios::trunc); // open file and remove any existing data

        std::string delimiter = " "; // space delimiter allows for easier reading when combining data

        file << std::setw(3) << std::setfill('0') << trialIndex << delimiter; // print with leading zeros
        file << trialStartTime << delimiter; // system clock time

        file << std::fixed;
        file << std::setprecision(3);

        for (int i = 0; i < frameCount; i++) { file << vDataVariablesEye[i].timestamp << delimiter; }

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
            if (vDataVariablesEye[i].DETECTED)
            {
                x.push_back(vDataVariablesEye[i].absoluteXPos - Parameters::camAOI.xPos);
                y.push_back(Parameters::eyeAOI.hght - (vDataVariablesEye[i].absoluteYPos - Parameters::camAOI.yPos));
                t.push_back(0.001 * vDataVariablesEye[i].timestamp);
            }
        }

        //        mQwtPlotWidget->plotTrajectory(x,y);
        mQwtPlotWidget->plotTimeSeries(x,y,t,0.001*trialTimeLength);
    }
}

void MainWindow::onFlashStandbySlider(int state)
{
    if (Parameters::CAMERA_RUNNING && Parameters::ONLINE_MODE && !TRIAL_RECORDING)
    {
        CamQImage       ->setVisible(!state);
        CamAOIXPosSlider->setVisible(!state);
        CamAOIYPosSlider->setVisible(!state);
        CamAOIWdthSlider->setVisible(!state);
        CamAOIHghtSlider->setVisible(!state);
        EyeAOIWdthSlider->setVisible(!state);
        EyeAOIHghtSlider->setVisible(!state);
        mQwtPlotWidget  ->setVisible( state);

        if (CameraHardwareGainAutoCheckBox->isChecked()) { mUEyeOpencvCam.setAutoGain(!state); }

        { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
            FLASH_STANDBY = state;
        }

        if (state == 0) { FlashStandbyLabel->setText("<font color='red'><b>OFF</b></font>"); }
        else
        {
            if (!boost::filesystem::exists(dataDirectory)) // check if data directory is valid
            {
                QString text = "Data directory does not exist. Please choose a valid path.";
                ConfirmationWindow mConfirmationWindow(text, false);
                mConfirmationWindow.setWindowTitle("Warning");
                mConfirmationWindow.exec();
                FlashStandbySlider->setValue(0);
                return;
            }

            dataFilename = (DataFilenameLineEdit->text()).toStdString();

            std::stringstream filename;
            filename << dataDirectory
                     << "/"
                     << dataFilename
                     << ".dat";

            if (boost::filesystem::exists(filename.str())) // check if DAT file already exists
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

            FlashStandbyLabel->setText("<font color='green'><b>ON</b></font>");
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
            vDataVariablesEye.resize( imageTotalOffline);
            vDataVariablesBead.resize(imageTotalOffline);

            vDetectionVariablesEye.resize( imageTotalOffline + 1);
            vDetectionVariablesBead.resize(imageTotalOffline + 1);

            resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
            resetVariablesHard(mDetectionVariablesBead, mParameterWidgetBead->getStructure(), Parameters::beadAOI);

            vDetectionVariablesEye [0] = mDetectionVariablesEye;
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
        { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
            Parameters::camAOI.wdth = eyeImageRaw.cols;
            Parameters::camAOI.hght = eyeImageRaw.rows;
        }
        updateAOIx();
        updateAOIy();
        CamQImage->setAOIEye (Parameters::eyeAOI);
        CamQImage->setAOIBead(Parameters::beadAOI);
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
        CamQImage->loadImage(eyeImage);
        CamQImage->setImage();
    }
    else
    {
        CamQImage->clearImage();
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

    { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
        mVariableWidgetEye ->setWidgets(vDataVariablesEye    [imgIndex]);
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

    detectionVariables mDetectionVariablesEyeTemp;
    detectionVariables mDetectionVariablesBeadTemp;

    detectionParameters mDetectionParametersEyeTemp;
    detectionParameters mDetectionParametersBeadTemp;

    AOIProperties AOIEyeTemp;
    AOIProperties AOIBeadTemp;

    { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

        mDetectionVariablesEyeTemp  = vDetectionVariablesEye[imageIndex];
        mDetectionParametersEyeTemp = mParameterWidgetEye->getStructure();

        mDetectionVariablesBeadTemp  = vDetectionVariablesBead[imageIndex];
        mDetectionParametersBeadTemp = mParameterWidgetBead->getStructure();

        Parameters::camAOI.wdth = imageRaw.cols;
        Parameters::camAOI.hght = imageRaw.rows;

        if (mAdvancedOptions.CURVATURE_MEASUREMENT) { setCurvatureMeasurement(mDetectionParametersEyeTemp, imageRaw.cols); }

        updateAOIx();
        updateAOIy();

        AOIEyeTemp  = Parameters::eyeAOI;
        AOIBeadTemp = Parameters::beadAOI;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    detectionVariables mDetectionVariablesEyeNew = eyeStalker(imageRaw,
                                                              AOIEyeTemp,
                                                              mDetectionVariablesEyeTemp,
                                                              mDetectionParametersEyeTemp,
                                                              mDataVariablesEye,
                                                              mDrawVariablesEye,
                                                              mAdvancedOptions);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

    // Save data

    mDataVariablesEye.duration         = fp_ms.count();
    mDataVariablesEye.absoluteXPos     = mDataVariablesEye.exactXPos;
    mDataVariablesEye.absoluteYPos     = mDataVariablesEye.exactYPos;
    vDataVariablesEye[imageIndex]   = mDataVariablesEye;

    cv::Mat imageProcessed = imageRaw.clone();
    drawAll(imageProcessed, mDrawVariablesEye);

    detectionVariables mDetectionVariablesBeadNew;

    if (mParameterWidgetBead->getState())
    {
        mDetectionVariablesBeadNew      = eyeStalker(imageRaw, AOIBeadTemp, mDetectionVariablesBeadTemp, mDetectionParametersBeadTemp, mDataVariablesBead, mDrawVariablesBead);
        mDataVariablesBead.absoluteXPos = mDataVariablesBead.exactXPos;
        mDataVariablesBead.absoluteYPos = mDataVariablesBead.exactYPos;
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

    { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
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
        std::unique_lock<std::mutex> lck(mutexOffline);
        PROCESSING_ALL_IMAGES = false;
        cvOffline.notify_one(); // notify save-thread that (this) main-thread has exited while loop
    }
    {
        std::unique_lock<std::mutex> lck(mutexOffline);
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

        { std::unique_lock<std::mutex> lck(mutexOffline);
            while (PROCESSING_ALL_IMAGES) cvOffline.wait(lck); } // wait for main-thread to exit while loop

        onSaveTrialData();

        { std::unique_lock<std::mutex> lck(mutexOffline);
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
        file << cameraFrameRate   << ";";  // sampling rate
        if (timeMatrix.size() > 0) { file << (int) timeMatrix[trialIndexOffline][1] << ";"; } // system clock time

        file << std::fixed;
        file << std::setprecision(3);

        // write data

        for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].DETECTED << delimiter; }
        if (timeMatrix.size() > 0) { for (int i = 0; i < imageTotalOffline; i++) { file << timeMatrix[trialIndexOffline][i + 2] << delimiter; } }

        if (SAVE_POSITION)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].absoluteXPos  << delimiter; }
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].absoluteYPos  << delimiter; }
        }

        if (SAVE_CIRCUMFERENCE)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].exactCircumference << delimiter; }
        }

        if (SAVE_ASPECT_RATIO)
        {
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].exactAspectRatio << delimiter; }
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
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].duration                       << delimiter; }
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
            int numEdges = vDataVariablesEye[i].edgeData.size();
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].tag          << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].curvature    << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].curvatureMax << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].curvatureMin << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].length       << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].radius       << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].radiusVar    << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].intensity    << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vDataVariablesEye[i].edgeData[j].gradient     << delimiter; }
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
            int numFits = vDataVariablesEye[i].ellipseData.size();
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].tag           << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].xPos          << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].yPos          << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].circumference << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].aspectRatio   << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].fitError      << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].edgeLength    << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].angle         << delimiter; }
            for (int j = 0; j < numFits; j++) { file << vDataVariablesEye[i].ellipseData[j].edgeScore     << delimiter; }
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
    static const std::vector<double> parametersBead = {0.005,   //  0. Gain average
                                                       0.40,    //  1. Gain features
                                                       0.25,    //  2. Gain certainty
                                                       0.75,    //  3. Gain predicted
                                                       4,       //  4. Canny blur level
                                                       5,       //  5. Canny kernel size
                                                       300.0,   //  6. Canny threshold low
                                                       600.0,   //  7. Canny threshold high
                                                       2.5,     //  8. Curvature offset
                                                       0.05,    //  9. Ellipse edge fraction
                                                       4,       // 10. Maximum number of edges
                                                       0.60,    // 11. Maximum fit error
                                                       0,       // 12. Glint size
                                                       130,     // 13. Circumference max
                                                       90,      // 14. Circumference min
                                                       0.8,     // 15. Aspect ratio min
                                                       0.05,    // 16. Circumference change threshold upper
                                                       0.02,    // 17. Circumference change threshold lower
                                                       0.05,    // 18. Aspect ratio  change threshold upper
                                                       0.02,    // 19. Aspect ratio  change threshold lower
                                                       4,       // 20. Displacement  change threshold upper
                                                       2,       // 21. Displacement  change threshold lower
                                                       0.41,    // 22. Score threshold edge
                                                       0.29,    // 23. Score threshold fit
                                                       0.60,    // 24. Score difference threshold edge
                                                       0.10,    // 25. Score difference threshold fit
                                                       7,       // 26. Edge window length
                                                       6};      // 27. Maximum number of fits

    QSettings settings(filename, QSettings::IniFormat);

    camAOIRatio.hght                = settings.value("CamAOIHghtFraction",          camAOIRatioLeft.hght).toDouble();
    camAOIRatio.wdth                = settings.value("CamAOIWdthFraction",          camAOIRatioLeft.wdth).toDouble();
    camAOIRatio.xPos                = settings.value("CamAOIXPosFraction",          camAOIRatioLeft.xPos).toDouble();
    camAOIRatio.yPos                = settings.value("CamAOIYPosFraction",          camAOIRatioLeft.yPos).toDouble();
    cameraFrameRateDesired          = settings.value("CameraFrameRateDesired",      250).toInt();
    cameraSubSamplingFactor         = settings.value("SubSamplingFactor",             1).toInt();
    dataDirectory                   = settings.value("DataDirectory",                "").toString().toStdString();
    dataDirectoryOffline            = settings.value("DataDirectoryOffline",         "").toString();
    dataFilename                    = settings.value("DataFilename",                "experiment_data").toString().toStdString();
    Parameters::drawFlags.haar      = settings.value("DrawHaar",                  false).toBool();
    Parameters::drawFlags.edge      = settings.value("DrawEdge",                  false).toBool();
    Parameters::drawFlags.elps      = settings.value("DrawElps",                   true).toBool();
    trialIndexOffline               = settings.value("trialIndexOffline",             0).toInt();
    imageTotalOffline               = settings.value("imageTotalOffline",             0).toInt();
    flashThreshold                  = settings.value("FlashThreshold",              230).toInt();
    Parameters::eyeAOIRatio.xPos    = settings.value("AOIXPosRatio",                0.0).toDouble();
    Parameters::eyeAOIRatio.yPos    = settings.value("AOIYPosRatio",                0.0).toDouble();
    Parameters::eyeAOIRatio.hght    = settings.value("AOIHghtRatio",                1.0).toDouble();
    Parameters::eyeAOIRatio.wdth    = settings.value("AOIWdthRatio",                1.0).toDouble();
    Parameters::beadAOIRatio.xPos   = settings.value("AOIBeadXPosRatio",            0.2).toDouble();
    Parameters::beadAOIRatio.yPos   = settings.value("AOIBeadYPosRatio",            0.5).toDouble();
    Parameters::beadAOIRatio.hght   = settings.value("AOIBeadHghtRatio",            0.6).toDouble();
    Parameters::beadAOIRatio.wdth   = settings.value("AOIBeadWdthRatio",            0.3).toDouble();
    flashAOILeft.hght               = settings.value("FlashAOIHghtLeft",            100).toInt();
    flashAOILeft.wdth               = settings.value("FlashAOIWdthLeft",             60).toInt();
    flashAOILeft.xPos               = settings.value("FlashAOIXPosLeft",            227).toInt();
    flashAOILeft.yPos               = settings.value("FlashAOIYPosLeft",            500).toInt();
    flashAOIRght.hght               = settings.value("FlashAOIHghtRght",            100).toInt();
    flashAOIRght.wdth               = settings.value("FlashAOIWdthRght",             60).toInt();
    flashAOIRght.xPos               = settings.value("FlashAOIXPosRght",           1020).toInt();
    flashAOIRght.yPos               = settings.value("FlashAOIYPosRght",            450).toInt();
    SAVE_ASPECT_RATIO               = settings.value("SaveAspectRatio",             true).toBool();
    SAVE_CIRCUMFERENCE              = settings.value("SaveCircumference",           true).toBool();
    SAVE_POSITION                   = settings.value("SavePosition",                true).toBool();
    SAVE_EYE_IMAGE                  = settings.value("SaveEyeImage",                false).toBool();
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
    mDetectionParameters.gainAverages                       = settings.value(prefix + "GainAverage",                    parameters[ 0]).toDouble();
    mDetectionParameters.gainAppearance                     = settings.value(prefix + "GainAppearance",                 parameters[ 1]).toDouble();
    mDetectionParameters.gainCertainty                      = settings.value(prefix + "GainCertainty",                  parameters[ 2]).toDouble();
    mDetectionParameters.gainPosition                       = settings.value(prefix + "GainPosition",                   parameters[ 3]).toDouble();
    mDetectionParameters.cannyBlurLevel                     = settings.value(prefix + "CannyBlurLevel",                 parameters[ 4]).toInt();
    mDetectionParameters.cannyKernelSize                    = settings.value(prefix + "CannyKernelSize",                parameters[ 5]).toInt();
    mDetectionParameters.cannyThresholdLow                  = settings.value(prefix + "CannyThresholdLow",              parameters[ 6]).toDouble();
    mDetectionParameters.cannyThresholdHigh                 = settings.value(prefix + "CannyThresholdHigh",             parameters[ 7]).toDouble();
    mDetectionParameters.curvatureOffset                    = settings.value(prefix + "CurvatureOffset",                parameters[ 8]).toDouble();
    mDetectionParameters.fitEdgeFraction                    = settings.value(prefix + "FitEdgeFraction",                parameters[ 9]).toDouble();
    mDetectionParameters.fitEdgeMaximum                     = settings.value(prefix + "FitEdgeMaximum",                 parameters[10]).toInt();
    mDetectionParameters.thresholdFitError                  = settings.value(prefix + "ThresholdFitError",              parameters[11]).toDouble();
    mDetectionParameters.glintWdth                          = settings.value(prefix + "GlintSize",                      parameters[12]).toInt();
    mDetectionParameters.thresholdCircumferenceMax          = settings.value(prefix + "CircumferenceMax",               parameters[13]).toDouble();
    mDetectionParameters.thresholdCircumferenceMin          = settings.value(prefix + "CircumferenceMin",               parameters[14]).toDouble();
    mDetectionParameters.thresholdAspectRatioMin            = settings.value(prefix + "AspectRatioMin",                 parameters[15]).toDouble();
    mDetectionParameters.thresholdChangeCircumferenceUpper  = settings.value(prefix + "CircumferenceChangeUpper",       parameters[16]).toDouble();
    mDetectionParameters.thresholdChangeCircumferenceLower  = settings.value(prefix + "CircumferenceChangeLower",       parameters[17]).toDouble();
    mDetectionParameters.thresholdChangeAspectRatioUpper    = settings.value(prefix + "AspectRatioChangeUpper",         parameters[18]).toDouble();
    mDetectionParameters.thresholdChangeAspectRatioLower    = settings.value(prefix + "AspectRatioChangeLower",         parameters[19]).toDouble();
    mDetectionParameters.thresholdChangePositionUpper       = settings.value(prefix + "PositionChangeUpper",            parameters[20]).toDouble();
    mDetectionParameters.thresholdChangePositionLower       = settings.value(prefix + "PositionChangeLower",            parameters[21]).toDouble();
    mDetectionParameters.thresholdScoreEdge                 = settings.value(prefix + "ScoreThresholdEdge",             parameters[22]).toDouble();
    mDetectionParameters.thresholdScoreFit                  = settings.value(prefix + "ScoreThresholdFit",              parameters[23]).toDouble();
    mDetectionParameters.thresholdScoreDiffEdge             = settings.value(prefix + "ScoreThresholdDiffEdge",         parameters[24]).toDouble();
    mDetectionParameters.thresholdScoreDiffFit              = settings.value(prefix + "ScoreThresholdDiffFit",          parameters[25]).toDouble();
    mDetectionParameters.windowLengthEdge                   = settings.value(prefix + "WindowLengthEdge",               parameters[26]).toDouble();
    mDetectionParameters.fitMaximum                         = settings.value(prefix + "FitMaximum",                     parameters[27]).toDouble();

    return mDetectionParameters;
}

void MainWindow::saveSettings(QString filename)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue("AOIHghtRatio",           Parameters::eyeAOIRatio.hght);
    settings.setValue("AOIWdthRatio",           Parameters::eyeAOIRatio.wdth);
    settings.setValue("AOIXPosRatio",           Parameters::eyeAOIRatio.xPos);
    settings.setValue("AOIYPosRatio",           Parameters::eyeAOIRatio.yPos);
    settings.setValue("AOIBeadHghtRatio",       Parameters::beadAOIRatio.hght);
    settings.setValue("AOIBeadWdthRatio",       Parameters::beadAOIRatio.wdth);
    settings.setValue("AOIBeadXPosRatio",       Parameters::beadAOIRatio.xPos);
    settings.setValue("AOIBeadYPosRatio",       Parameters::beadAOIRatio.yPos);
    settings.setValue("DataDirectory",          QString::fromStdString(dataDirectory));
    settings.setValue("DataDirectoryOffline",   dataDirectoryOffline);
    settings.setValue("DrawHaar",               Parameters::drawFlags.haar);
    settings.setValue("DrawEdge",               Parameters::drawFlags.edge);
    settings.setValue("DrawElps",               Parameters::drawFlags.elps);
    settings.setValue("GainAuto",               CameraHardwareGainAutoCheckBox ->checkState());
    settings.setValue("GainBoost",              CameraHardwareGainBoostCheckBox->checkState());
    settings.setValue("CamAOIHghtFraction",     camAOIRatio.hght);
    settings.setValue("CamAOIWdthFraction",     camAOIRatio.wdth);
    settings.setValue("CamAOIXPosFraction",     camAOIRatio.xPos);
    settings.setValue("CamAOIYPosFraction",     camAOIRatio.yPos);
    settings.setValue("CameraFrameRateDesired", cameraFrameRateDesired);
    settings.setValue("DataFilename",           QString::fromStdString(dataFilename));
    settings.setValue("FlashAOIHghtLeft",       flashAOILeft.hght);
    settings.setValue("FlashAOIWdthLeft",       flashAOILeft.wdth);
    settings.setValue("FlashAOIXPosLeft",       flashAOILeft.xPos);
    settings.setValue("FlashAOIYPosLeft",       flashAOILeft.yPos);
    settings.setValue("FlashAOIHghtRght",       flashAOIRght.hght);
    settings.setValue("FlashAOIWdthRght",       flashAOIRght.wdth);
    settings.setValue("FlashAOIXPosRght",       flashAOIRght.xPos);
    settings.setValue("FlashAOIYPosRght",       flashAOIRght.yPos);
    settings.setValue("FlashThreshold",         flashThreshold);
    settings.setValue("SaveAspectRatio",        SAVE_ASPECT_RATIO);
    settings.setValue("SaveCircumference",      SAVE_CIRCUMFERENCE);
    settings.setValue("SavePosition",           SAVE_POSITION);
    settings.setValue("SaveEyeImage",           SAVE_EYE_IMAGE);
    settings.setValue("SubSamplingFactor",      cameraSubSamplingFactor);
    settings.setValue("TrialTimeLength",        TrialTimeLengthLineEdit->text().toInt());

    detectionParameters mDetectionParametersEye  = mParameterWidgetEye->getStructure();
    saveParameters(filename,  "Eye", mDetectionParametersEye);

    detectionParameters mDetectionParametersBead = mParameterWidgetBead->getStructure();
    saveParameters(filename, "Bead", mDetectionParametersBead);
}

void MainWindow::saveParameters(QString filename, QString prefix, detectionParameters mDetectionParameters)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue(prefix + "GainAverage",               mDetectionParameters.gainAverages);
    settings.setValue(prefix + "GainAppearance",            mDetectionParameters.gainAppearance);
    settings.setValue(prefix + "GainPosition",              mDetectionParameters.gainPosition);
    settings.setValue(prefix + "GainCertainty",             mDetectionParameters.gainCertainty);
    settings.setValue(prefix + "CannyBlurLevel",            mDetectionParameters.cannyBlurLevel);
    settings.setValue(prefix + "CannyKernelSize",           mDetectionParameters.cannyKernelSize);
    settings.setValue(prefix + "CannyThresholdLow",         mDetectionParameters.cannyThresholdLow);
    settings.setValue(prefix + "CannyThresholdHigh",        mDetectionParameters.cannyThresholdHigh);
    settings.setValue(prefix + "CircumferenceMax",          mDetectionParameters.thresholdCircumferenceMax);
    settings.setValue(prefix + "CircumferenceMin",          mDetectionParameters.thresholdCircumferenceMin);
    settings.setValue(prefix + "CurvatureOffset",           mDetectionParameters.curvatureOffset);
    settings.setValue(prefix + "FitEdgeMaximum",            mDetectionParameters.fitEdgeMaximum);
    settings.setValue(prefix + "FitMaximum",                mDetectionParameters.fitMaximum);
    settings.setValue(prefix + "ThresholdFitError",         mDetectionParameters.thresholdFitError);
    settings.setValue(prefix + "AspectRatioMin",            mDetectionParameters.thresholdAspectRatioMin);
    settings.setValue(prefix + "GlintSize",                 mDetectionParameters.glintWdth);
    settings.setValue(prefix + "CircumferenceChangeUpper",  mDetectionParameters.thresholdChangeCircumferenceUpper);
    settings.setValue(prefix + "CircumferenceChangeLower",  mDetectionParameters.thresholdChangeCircumferenceLower);
    settings.setValue(prefix + "AspectRatioChangeUpper",    mDetectionParameters.thresholdChangeAspectRatioUpper);
    settings.setValue(prefix + "AspectRatioChangeLower",    mDetectionParameters.thresholdChangeAspectRatioLower);
    settings.setValue(prefix + "PositionChangeUpper",       mDetectionParameters.thresholdChangePositionUpper);
    settings.setValue(prefix + "PositionChangeLower",       mDetectionParameters.thresholdChangePositionLower);
    settings.setValue(prefix + "ScoreThresholdEdge",        mDetectionParameters.thresholdScoreEdge);
    settings.setValue(prefix + "ScoreThresholdFit",         mDetectionParameters.thresholdScoreFit);
    settings.setValue(prefix + "ScoreThresholdDiffEdge",    mDetectionParameters.thresholdScoreDiffEdge);
    settings.setValue(prefix + "ScoreThresholdDiffFit",     mDetectionParameters.thresholdScoreDiffFit);
    settings.setValue(prefix + "WindowLengthEdge",          mDetectionParameters.windowLengthEdge);
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

void MainWindow::onSetOnlineProcessing(int state)
{
    if (state) { SAVE_EYE_IMAGE = false; }
    else       { SAVE_EYE_IMAGE = true;  }
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

void MainWindow::onSetCameraFrameRateDesired(int value)
{
    cameraFrameRateDesired = value;
}

void MainWindow::onCalibrateFrameRate()
{
    SET_FRAME_RATE = true;
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
void MainWindow::onSetCameraGainBoost     (int state) { mUEyeOpencvCam.setGainBoost(state);      }
void MainWindow::onSetCameraAutoGain      (int state) { mUEyeOpencvCam.setAutoGain(state);       }

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
                if (mUEyeOpencvCam.allocateMemory(Parameters::camAOI.wdth, Parameters::camAOI.hght))
                {
                    mUEyeOpencvCam.setAOI(Parameters::camAOI.xPos, Parameters::camAOI.yPos, Parameters::camAOI.wdth, Parameters::camAOI.hght);
                    Parameters::CAMERA_READY = true;
                }
            }
        }

        getCameraParameters();
    }
}

void MainWindow::onCropAOI()
{
    double fracXPos;
    double fracYPos;
    double fracWdth;
    double fracHght;

    {
        std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
        std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);

        int absXPos = Parameters::eyeAOI.xPos + Parameters::camAOI.xPos;
        int absYPos = Parameters::eyeAOI.yPos + Parameters::camAOI.yPos;

        fracXPos = absXPos / (double) cameraAOIWdthMax;
        fracYPos = absYPos / (double) cameraAOIHghtMax;
        fracWdth = (Parameters::eyeAOI.wdth - cameraAOIWdthMin) / (double) (cameraAOIWdthMax - cameraAOIWdthMin);
        fracHght = (Parameters::eyeAOI.hght - cameraAOIHghtMin) / (double) (cameraAOIHghtMax - cameraAOIHghtMin);

        Parameters::eyeAOIRatio.xPos = 1.0;
        Parameters::eyeAOIRatio.yPos = 1.0;
        Parameters::eyeAOIRatio.wdth = 1.0;
        Parameters::eyeAOIRatio.hght = 1.0;
    }

    CamAOIXPosSlider->setDoubleValue(fracXPos);
    CamAOIYPosSlider->setDoubleValue(fracYPos);
    CamAOIWdthSlider->setDoubleValue(fracWdth);
    CamAOIHghtSlider->setDoubleValue(fracHght);
    onSetCamAOI();
}

void MainWindow::updateCamAOIx()
{
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
    Parameters::camAOI.wdth = floor(((cameraAOIWdthMax - cameraAOIWdthMin) * camAOIRatio.wdth + cameraAOIWdthMin) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    Parameters::camAOI.xPos = floor((cameraAOIWdthMax * camAOIRatio.xPos) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    CamAOIXPosSlider->setDoubleMaximum((cameraAOIWdthMax - Parameters::camAOI.wdth) / (double) cameraAOIWdthMax);
    updateAOIx();
}

void MainWindow::updateCamAOIy()
{
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
    Parameters::camAOI.hght = floor(((cameraAOIHghtMax - cameraAOIHghtMin) * camAOIRatio.hght + cameraAOIHghtMin) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    Parameters::camAOI.yPos = floor((cameraAOIHghtMax * camAOIRatio.yPos) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    CamAOIYPosSlider->setDoubleMaximum((cameraAOIHghtMax - Parameters::camAOI.hght) / (double) cameraAOIHghtMax);
    updateAOIy();
}

void MainWindow::updateAOIx()
{
    { std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
        Parameters::eyeAOI.wdth = round(Parameters::camAOI.wdth * Parameters::eyeAOIRatio.wdth);
        Parameters::eyeAOI.xPos = round(Parameters::camAOI.wdth * Parameters::eyeAOIRatio.xPos);
        if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::camAOI.wdth)
        {   Parameters::eyeAOI.xPos = Parameters::camAOI.wdth - Parameters::eyeAOI.wdth; }
    }

    { std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
        Parameters::beadAOI.wdth = round(Parameters::camAOI.wdth * Parameters::beadAOIRatio.wdth);
        Parameters::beadAOI.xPos = round(Parameters::camAOI.wdth * Parameters::beadAOIRatio.xPos);
        if (Parameters::beadAOI.xPos + Parameters::beadAOI.wdth > Parameters::camAOI.wdth)
        {   Parameters::beadAOI.xPos = Parameters::camAOI.wdth - Parameters::beadAOI.wdth; }
    }
}

void MainWindow::updateAOIy()
{
    { std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
        Parameters::eyeAOI.hght = round(Parameters::camAOI.hght * Parameters::eyeAOIRatio.hght);
        Parameters::eyeAOI.yPos = round(Parameters::camAOI.hght * Parameters::eyeAOIRatio.yPos);
        if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght > Parameters::camAOI.hght)
        {   Parameters::eyeAOI.yPos = Parameters::camAOI.hght - Parameters::eyeAOI.hght; }
    }

    { std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
        Parameters::beadAOI.hght = round(Parameters::camAOI.hght * Parameters::beadAOIRatio.hght);
        Parameters::beadAOI.yPos = round(Parameters::camAOI.hght * Parameters::beadAOIRatio.yPos);
        if (Parameters::beadAOI.yPos + Parameters::beadAOI.hght > Parameters::camAOI.hght)
        {   Parameters::beadAOI.yPos = Parameters::camAOI.hght - Parameters::beadAOI.hght; }
    }
}

void MainWindow::onSetCamAOI()
{
    std::lock_guard<std::mutex> AOILock_1(mutexAOI_1);
    std::lock_guard<std::mutex> AOILock_2(mutexAOI_2);

    Parameters::camAOI = camAOITemp;

    updateAOIx();
    updateAOIy();

    if (mUEyeOpencvCam.freeImageMemory())
    {
        if (mUEyeOpencvCam.allocateMemory(Parameters::camAOI.wdth, Parameters::camAOI.hght))
        {
            if
                    (Parameters::camAOI.xPos + Parameters::camAOI.wdth <= cameraAOIWdthMax &&
                     Parameters::camAOI.yPos + Parameters::camAOI.hght <= cameraAOIHghtMax)
            {
                if (mUEyeOpencvCam.setAOI(Parameters::camAOI.xPos, Parameters::camAOI.yPos, Parameters::camAOI.wdth, Parameters::camAOI.hght))
                {
                    Parameters::CAMERA_READY = true;
                    getCameraParameters();
                }
            }
        }
    }
}

void MainWindow::onSetCamAOIWdth(double fraction)
{
    camAOIRatio.wdth = fraction;
    camAOITemp.wdth  = floor(((cameraAOIWdthMax - cameraAOIWdthMin) * camAOIRatio.wdth + cameraAOIWdthMin) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    CamAOIXPosSlider->setDoubleMaximum((cameraAOIWdthMax - camAOITemp.wdth) / (double) cameraAOIWdthMax);
}

void MainWindow::onSetCamAOIHght(double fraction)
{
    camAOIRatio.hght = fraction;
    camAOITemp.hght  = floor(((cameraAOIHghtMax - cameraAOIHghtMin) * camAOIRatio.hght + cameraAOIHghtMin) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    CamAOIYPosSlider->setDoubleMaximum((cameraAOIHghtMax - camAOITemp.hght) / (double) cameraAOIHghtMax);
}

void MainWindow::onSetCamAOIXPos(double fraction)
{
    camAOIRatio.xPos = fraction;
    camAOITemp.xPos  = floor((cameraAOIWdthMax * camAOIRatio.xPos) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
}

void MainWindow::onSetCamAOIYPos(double fraction)
{
    camAOIRatio.yPos = fraction;
    camAOITemp.yPos  = floor((cameraAOIHghtMax * camAOIRatio.yPos) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
}

void MainWindow::onSetEyeAOIWdth(double ratio)
{
    std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
    Parameters::eyeAOIRatio.wdth = ratio;
    Parameters::eyeAOI.wdth = round(Parameters::camAOI.wdth * Parameters::eyeAOIRatio.wdth);
    if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::camAOI.wdth)
    {   Parameters::eyeAOI.xPos = Parameters::camAOI.wdth - Parameters::eyeAOI.wdth; }
}

void MainWindow::onSetEyeAOIHght(double ratio)
{
    std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
    Parameters::eyeAOIRatio.hght = ratio;
    Parameters::eyeAOI.hght      = round(Parameters::camAOI.hght * Parameters::eyeAOIRatio.hght);
    if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght > Parameters::camAOI.hght)
    {   Parameters::eyeAOI.yPos = Parameters::camAOI.hght - Parameters::eyeAOI.hght; }
}

void MainWindow::onSetBeadAOIWdth(double ratio)
{
    std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
    Parameters::beadAOIRatio.wdth = ratio;
    Parameters::beadAOI.wdth      = round(Parameters::camAOI.wdth * Parameters::beadAOIRatio.wdth);
    if (Parameters::beadAOI.xPos + Parameters::beadAOI.wdth > Parameters::camAOI.wdth)
    {   Parameters::beadAOI.xPos = Parameters::camAOI.wdth - Parameters::beadAOI.wdth; }
}

void MainWindow::onSetBeadAOIHght(double ratio)
{
    std::lock_guard<std::mutex> AOIBeadLock(Parameters::AOIBeadMutex);
    Parameters::beadAOIRatio.hght = ratio;
    Parameters::beadAOI.hght      = round(Parameters::camAOI.hght * Parameters::beadAOIRatio.hght);
    if (Parameters::beadAOI.yPos + Parameters::beadAOI.hght > Parameters::camAOI.hght)
    {   Parameters::beadAOI.yPos = Parameters::camAOI.hght  - Parameters::beadAOI.hght; }
}

void MainWindow::onSetFlashAOIXPosLeft(int val) { flashAOILeft.xPos = val; }
void MainWindow::onSetFlashAOIYPosLeft(int val) { flashAOILeft.yPos = val; }
void MainWindow::onSetFlashAOIWdthLeft(int val) { flashAOILeft.wdth = val; }
void MainWindow::onSetFlashAOIHghtLeft(int val) { flashAOILeft.hght = val; }
void MainWindow::onSetFlashAOIXPosRght(int val) { flashAOIRght.xPos = val; }
void MainWindow::onSetFlashAOIYPosRght(int val) { flashAOIRght.yPos = val; }
void MainWindow::onSetFlashAOIWdthRght(int val) { flashAOIRght.wdth = val; }
void MainWindow::onSetFlashAOIHghtRght(int val) { flashAOIRght.hght = val; }

void MainWindow::onSetFlashThreshold(int val)
{
    flashThreshold = val;
    FlashThresholdLabel->setText(QString::number(flashThreshold));
}

void MainWindow::onSetAOIEyeLeft()
{
    CamAOIXPosSlider->setDoubleValue(camAOIRatioLeft.xPos);
    CamAOIYPosSlider->setDoubleValue(camAOIRatioLeft.yPos);
    CamAOIWdthSlider->setDoubleValue(camAOIRatioLeft.wdth);
    CamAOIHghtSlider->setDoubleValue(camAOIRatioLeft.hght);
    onSetCamAOI();

    flashAOI = flashAOILeft;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::onSetAOIEyeRght()
{
    CamAOIXPosSlider->setDoubleValue(camAOIRatioRght.xPos);
    CamAOIYPosSlider->setDoubleValue(camAOIRatioRght.yPos);
    CamAOIWdthSlider->setDoubleValue(camAOIRatioRght.wdth);
    CamAOIHghtSlider->setDoubleValue(camAOIRatioRght.hght);
    onSetCamAOI();

    flashAOI = flashAOIRght;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::onSetTrialIndex           (int val)   { trialIndex                 =   val; }
void MainWindow::onSetSaveDataAspectRatio  (int state) { SAVE_ASPECT_RATIO          = state; }
void MainWindow::onSetSaveDataCircumference(int state) { SAVE_CIRCUMFERENCE         = state; }
void MainWindow::onSetSaveDataPosition     (int state) { SAVE_POSITION              = state; }
void MainWindow::onSetDrawHaar             (int state) { Parameters::drawFlags.haar = state; }
void MainWindow::onSetDrawEdge             (int state) { Parameters::drawFlags.edge = state; }
void MainWindow::onSetDrawElps             (int state) { Parameters::drawFlags.elps = state; }

void MainWindow::onSetCurvatureMeasurement(int state) { mAdvancedOptions.CURVATURE_MEASUREMENT = state; }
void MainWindow::setCurvatureMeasurement(detectionParameters& mDetectionParameters, int imgWdth)
{
    mDetectionParameters.thresholdAspectRatioMin    = 0.0;
    mDetectionParameters.thresholdCircumferenceMax  = M_PI * imgWdth;
    mDetectionParameters.thresholdCircumferenceMin  = 1.0; // 1.0 avoids inf
    mDetectionParameters.thresholdScoreEdge         = 0.0;
    mDetectionParameters.thresholdScoreFit          = 0.0;
    mDetectionParameters.thresholdScoreDiffEdge     = 1.0;
    mDetectionParameters.curvatureOffset            = 360;
    mDetectionParameters.glintWdth                  = 0.0;
}

void MainWindow::onSetSaveDataEdge (int state) { SAVE_DATA_EDGE  = state; }
void MainWindow::onSetSaveDataFit  (int state) { SAVE_DATA_FIT   = state; }
void MainWindow::onSetSaveDataExtra(int state) { SAVE_DATA_EXTRA = state; }

double MainWindow::flashDetection(const cv::Mat& imgBGR)
{
    int imgSize = imgBGR.cols * imgBGR.rows;

    if (imgSize > 0)
    {
        unsigned long long intensityTotal = 0;
        cv::Mat img;
        cv::cvtColor(imgBGR, img, cv::COLOR_BGR2GRAY);
        uchar *ptr = img.data;
        for (int iPixel = 0; iPixel < imgSize; iPixel++) { intensityTotal += ptr[iPixel]; }
        return (intensityTotal / (double) imgSize);
    } else { return 0; }
}
