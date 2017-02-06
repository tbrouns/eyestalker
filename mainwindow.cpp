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

    Parameters::cameraXResolution       = 1280;
    Parameters::cameraYResolution       = 1024;
    Parameters::CAMERA_READY            = false;
    Parameters::CAMERA_RUNNING          = false;
    Parameters::ONLINE_PROCESSING       = true;
    Parameters::cannyKernelSize         = 3;
    Parameters::ellipseDrawCrossSize    = 5;
    Parameters::ellipseDrawOutlineWidth = 0.032;

    APP_EXIT    = false;
    APP_RUNNING = true;
    cameraAOIFractionHghtDefaultLeft = 0.19;
    cameraAOIFractionHghtDefaultRght = 0.22;
    cameraAOIFractionWdthDefaultLeft = 0.25;
    cameraAOIFractionWdthDefaultRght = 0.31;
    cameraAOIFractionXPosDefaultLeft = 0.20;
    cameraAOIFractionXPosDefaultRght = 0.52;
    cameraAOIFractionYPosDefaultLeft = 0.41;
    cameraAOIFractionYPosDefaultRght = 0.37;
    cameraAOIHghtMin        = 4;
    cameraAOIHghtStepSize   = 2;
    cameraAOIWdthMin        = 32;
    cameraAOIWdthStepSize   = 4;
    cameraPixelClock        = 24;
    cameraSubSamplingFactor = 2;
    camImageHght            = 200;
    camImageWdth            = 480; // size of image in widget
    PROCESSING_ALL_IMAGES   = false;
    PROCESSING_ALL_TRIALS   = false;
    imageIndexOffline       = 0;
    imageTotalOffline       = 0;
    TRIAL_RECORDING         = false;
    eyeAOIHghtMin           = 75;
    eyeAOIWdthMin           = 100;
    eyeImageHght            = 200;
    eyeImageWdth            = 320;
    FLASH_STANDBY           = false;
    frameCount              = 0;
    guiUpdateFrequency      = 30;
    relativeTime            = 0;
    SAVE_EYE_IMAGE          = true;
    startTime               = 0;
    subjectIdentifier       = "";
    trialIndex              = 0;

    mParameterWidgetEye  = new ParameterWidget;
    mParameterWidgetBead = new ParameterWidget;

    mVariableWidgetEye  = new VariableWidget;
    mVariableWidgetBead = new VariableWidget;

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
    CamQImage->setAOIEye  (Parameters::eyeAOI);
    CamQImage->setAOIBead (Parameters::beadAOI);
    CamQImage->setAOIFlash(flashAOI);

    CamQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    CamQImage->loadImage(imgCam);
    CamQImage->setImage();
    QObject::connect(CamQImage, SIGNAL(updateImage(int)), this , SLOT(updateImageRaw(int)));

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

    EyeWdthAOISlider = new SliderDouble;
    EyeWdthAOISlider->setPrecision(2);
    EyeWdthAOISlider->setDoubleRange(0, 1.0);
    EyeWdthAOISlider->setOrientation(Qt::Horizontal);
    QObject::connect(EyeWdthAOISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEyeAOIWdth(double)));

    EyeHghtAOISlider = new SliderDouble;
    EyeHghtAOISlider->setPrecision(2);
    EyeHghtAOISlider->setDoubleRange(0, 1.0);
    EyeHghtAOISlider->setOrientation(Qt::Vertical);
    QObject::connect(EyeHghtAOISlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setEyeAOIHght(double)));

    EyeWdthAOISlider->setDoubleValue(eyeAOIWdthFraction);
    EyeHghtAOISlider->setDoubleValue(eyeAOIHghtFraction);

    QPushButton* AOICropButton = new QPushButton("&Crop AOI");
    QObject::connect(AOICropButton, SIGNAL(clicked()), this, SLOT(onCropAOI()));

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

    QObject::connect(HaarCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setDrawHaar(int)));
    QObject::connect(EdgeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setDrawEdge(int)));
    QObject::connect(ElpsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setDrawElps(int)));

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

    QCheckBox *RealTimeEyeTrackingCheckBox = new QCheckBox;
    RealTimeEyeTrackingCheckBox->setChecked(!SAVE_EYE_IMAGE);
    QObject::connect(RealTimeEyeTrackingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetRealTime(int)));

    QLabel *OfflineModeTextBox = new QLabel;
    OfflineModeTextBox->setText("Offline mode:");

    QCheckBox *OfflineModeCheckBox = new QCheckBox;
    OfflineModeCheckBox->setChecked(false);
    QObject::connect(OfflineModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onSetOfflineMode(int)));

    QPushButton *ResetParametersPushButton = new QPushButton("Reset parameters");
    QObject::connect(ResetParametersPushButton, SIGNAL(clicked()), this, SLOT(onResetParameters()));

    QLabel *BeadDetectionTextBox = new QLabel;
    BeadDetectionTextBox->setText("Bead detection:");

    QCheckBox *BeadDetectionCheckBox = new QCheckBox;
    BeadDetectionCheckBox->setChecked(mParameterWidgetBead->getState());
    QObject::connect(BeadDetectionCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setBeadDetection(int)));

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

    QLabel *OfflinePupilDetectionTextBox = new QLabel;
    OfflinePupilDetectionTextBox->setText("<b>Detect pupil: </b> ");

    QPushButton *OfflinePupilDetectionOneButton = new QPushButton("Current");
    QObject::connect(OfflinePupilDetectionOneButton, SIGNAL(clicked(bool)), this, SLOT(onDetectCurrentFrame()));

    QPushButton *OfflinePupilDetectionAllFramesButton = new QPushButton("All frames");
    QObject::connect(OfflinePupilDetectionAllFramesButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllFrames()));

    QPushButton *OfflinePupilDetectionAllTrialsButton = new QPushButton("All trials");
    QObject::connect(OfflinePupilDetectionAllTrialsButton, SIGNAL(clicked(bool)), this, SLOT(onDetectAllTrials()));

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
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionTextBox);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionOneButton);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionAllFramesButton);
    EyeTrackingOfflineLayout->addWidget(OfflinePupilDetectionAllTrialsButton);
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
    CameraOutputLayout->addWidget(CamQImage, 1, 2, Qt::AlignCenter);
    CameraOutputLayout->addWidget(CamEyeAOIWdthSlider, 2, 2);
    CameraOutputLayout->addWidget(CamEyeAOIHghtSlider, 1, 3);
    CameraOutputLayout->addWidget(EyeHghtAOISlider, 1, 5);
    CameraOutputLayout->addWidget(EyeQImage, 1, 4, Qt::AlignCenter);
    CameraOutputLayout->addWidget(EyeWdthAOISlider, 2, 4);
    CameraOutputLayout->addLayout(OptionsLayout, 4, 2);
    CameraOutputLayout->addLayout(EyeDetectionsLayout, 3, 2);
    CameraOutputLayout->addLayout(DrawFunctionsLayout, 3, 4, Qt::AlignCenter);

    CameraOutputLayout->setColumnStretch(0, 1);
    CameraOutputLayout->setColumnStretch(6, 1);

    ///////////////////////////////////////////////////////////////
    /////////////////////// TAB WIDGET  ///////////////////////////
    ///////////////////////////////////////////////////////////////

    /////////////////// Camera settings tab ///////////////////////

    // Pixel clock

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

    // Gain level

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

    // Sub-sampling

    QLabel *CameraSubSamplingTextBox = new QLabel;
    CameraSubSamplingTextBox->setText("<b>Sub-sampling: </b>");

    QCheckBox *CameraSubSamplingCheckBox = new QCheckBox;

    if (cameraSubSamplingFactor == 2) { CameraSubSamplingCheckBox->setChecked(true);  }
    else                              { CameraSubSamplingCheckBox->setChecked(false); }

    QObject::connect(CameraSubSamplingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setCameraSubSampling(int)));

    // Set-up layout

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

    /////////////////// Experiment tab ///////////////////////

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

    QObject::connect(FlashThresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(setFlashThreshold(int)));

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
    QObject::connect(StartRecordingButton, SIGNAL(clicked()), this, SLOT(startRecordingManual()));

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

    QObject::connect(SaveDataAspectRatioCheckBox,   SIGNAL(stateChanged(int)), this, SLOT(setSaveDataAspectRatio(int)));
    QObject::connect(SaveDataCircumferenceCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setSaveDataCircumference(int)));
    QObject::connect(SaveDataPositionCheckBox,      SIGNAL(stateChanged(int)), this, SLOT(setSaveDataPosition(int)));

    // Set-up layouts

    QWidget* ExperimentTabWidget = new QWidget;
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

    MainTabWidget = new QTabWidget;
    MainTabWidget->addTab(CameraParametersWidget,  tr("Camera"));
    MainTabWidget->addTab(ExperimentTabWidget,     tr("Experimental"));
    MainTabWidget->addTab(EyeTrackingScrollArea,   tr("Eye-tracking"));
    if (mParameterWidgetBead->getState()) { MainTabWidget->addTab(BeadTrackingScrollArea,  tr("Bead-tracking")); }


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

    QAction *about = new QAction("&About EyeStalker...", this);

    QMenu *file;
    file = menuBar()->addMenu("&Help");
    file->addAction(about);

    QObject::connect(about, &QAction::triggered, this, &MainWindow::onDialogueOpen);

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
    QObject::connect(UpdateCameraImageTimer, SIGNAL(timeout()), this, SLOT(updateCameraImage()));
    QObject::connect(this, SIGNAL(startTimer(int)), UpdateCameraImageTimer, SLOT(start(int)));
    QObject::connect(this, SIGNAL(stopTimer()),     UpdateCameraImageTimer, SLOT(stop()));
    emit startTimer(round(1000 / guiUpdateFrequency));
}

MainWindow::~MainWindow()
{

}

void MainWindow::pupilTracking()
{
    while(APP_RUNNING && Parameters::CAMERA_RUNNING && Parameters::ONLINE_PROCESSING)
    {
        detectionProperties mDetectionPropertiesEyeTemp;

        dataVariables mDataVariablesTemp;
        drawVariables mDrawVariablesTemp;

        AOIProperties flashAOITemp;

        flashAOITemp.xPos = 0;
        flashAOITemp.yPos = 0;
        flashAOITemp.wdth = 0;
        flashAOITemp.hght = 0;

        bool FLASH_REGION_VISIBLE = true;

        imageInfo mImageInfo = mUEyeOpencvCam.getFrame(); // get new frame from camera

        cv::Mat imageOriginal = mImageInfo.image;
        absoluteTime = mImageInfo.time; // Get frame timestamp

        double relativeTimeNew = (absoluteTime - startTime) / (double) 10000; // in ms

        // ignore frame if time interval was too short (possible camera error)
        if (relativeTimeNew <= (relativeTime + 0.9 * (1000 / cameraFrameRate))) { continue; }

        relativeTime = relativeTimeNew;

        int cameraAOIXPos;
        int cameraAOIYPos;
        int cameraAOIWdth;
        int cameraAOIHght;

        { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

            imageCamera = imageOriginal.clone();

            mDetectionPropertiesEyeTemp.v = mDetectionVariablesEye;
            mDetectionPropertiesEyeTemp.p = mParameterWidgetEye->getStructure();

            cameraAOIXPos = Parameters::cameraAOI.xPos;
            cameraAOIYPos = Parameters::cameraAOI.yPos;
            cameraAOIWdth = Parameters::cameraAOI.wdth;
            cameraAOIHght = Parameters::cameraAOI.hght;

            mDetectionPropertiesEyeTemp.p.AOIXPos = Parameters::eyeAOI.xPos;
            mDetectionPropertiesEyeTemp.p.AOIYPos = Parameters::eyeAOI.yPos;
            mDetectionPropertiesEyeTemp.p.AOIWdth = Parameters::eyeAOI.wdth;
            mDetectionPropertiesEyeTemp.p.AOIHght = Parameters::eyeAOI.hght;

            if (!TRIAL_RECORDING)
            {
                flashAOITemp.xPos = flashAOI.xPos - cameraAOIXPos;
                flashAOITemp.yPos = flashAOI.yPos - cameraAOIYPos;
                flashAOITemp.wdth = flashAOI.wdth;
                flashAOITemp.hght = flashAOI.hght;

                if (flashAOITemp.xPos < 0)
                {
                    if (flashAOITemp.xPos + flashAOITemp.wdth < 0) // Flash AOI not in camera AOI
                    {
                        FlashStandbySlider->setValue(0);
                        FLASH_REGION_VISIBLE = false;
                    }
                    else
                    {
                        flashAOITemp.wdth = flashAOITemp.xPos + flashAOITemp.wdth;
                        flashAOITemp.xPos = 0;
                    }
                }
                else if (flashAOITemp.xPos + flashAOITemp.wdth >= cameraAOIWdth)
                {
                    flashAOITemp.wdth = cameraAOIWdth - flashAOITemp.xPos;

                    if (flashAOITemp.wdth <= 0)
                    {
                        FlashStandbySlider->setValue(0);
                        FLASH_REGION_VISIBLE = false;
                    }
                }

                if (flashAOITemp.yPos < 0)
                {
                    if (flashAOITemp.yPos + flashAOITemp.hght < 0) // Flash AOI not in camera AOI
                    {
                        FlashStandbySlider->setValue(0);
                        FLASH_REGION_VISIBLE = false;
                    }
                    else
                    {
                        flashAOITemp.hght = flashAOITemp.yPos + flashAOITemp.hght;
                        flashAOITemp.yPos = 0;
                    }
                }
                else if (flashAOITemp.yPos + flashAOITemp.hght >= cameraAOIHght)
                {
                    flashAOITemp.hght = cameraAOIHght - flashAOITemp.yPos;

                    if (flashAOITemp.hght <= 0)
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

        if (imgWdth < (mDetectionPropertiesEyeTemp.p.AOIXPos + mDetectionPropertiesEyeTemp.p.AOIWdth)) { continue; }
        if (imgHght < (mDetectionPropertiesEyeTemp.p.AOIYPos + mDetectionPropertiesEyeTemp.p.AOIHght)) { continue; }

        if (!TRIAL_RECORDING)
        {
            double avgIntensity = 0;

            if (FLASH_REGION_VISIBLE)
            {
                cv::Rect flashRegion(flashAOITemp.xPos, flashAOITemp.yPos, flashAOITemp.wdth, flashAOITemp.hght);
                avgIntensity = flashDetection(imageOriginal(flashRegion));
            }

            if (FLASH_STANDBY)
            {
                mVariableWidgetEye->resetStructure(mDetectionPropertiesEyeTemp.p);

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

                mDetectionPropertiesEyeTemp = eyeStalker(imageOriginal, mDetectionPropertiesEyeTemp, mDataVariablesTemp, mDrawVariablesTemp); // Pupil tracking algorithm
            }
        }
        else // Trial recording
        {
            if (!SAVE_EYE_IMAGE)
            {

                mDetectionPropertiesEyeTemp = eyeStalker(imageOriginal, mDetectionPropertiesEyeTemp, mDataVariablesTemp, mDrawVariablesTemp); // Pupil tracking algorithm

                mDataVariablesTemp.absoluteXPos = mDataVariablesTemp.exactXPos + cameraAOIXPos;
                mDataVariablesTemp.absoluteYPos = mDataVariablesTemp.exactYPos + cameraAOIYPos;
                mDataVariablesTemp.timestamp    = relativeTime; // save time stamps
                vDataVariables[frameCount]  = mDataVariablesTemp;

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
                emit startTimer(round(1000 / guiUpdateFrequency));
            }
        }

        // Update structures

        {
            std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

            mDetectionVariablesEye = mDetectionPropertiesEyeTemp.v;
            mDrawVariables         = mDrawVariablesTemp;
            mDataVariablesTemp         = mDataVariablesTemp;
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


void MainWindow::updateCameraImage()
{
    if (Parameters::ONLINE_PROCESSING)
    {
        if (!TRIAL_RECORDING)
        {
            if (Parameters::CAMERA_RUNNING)
            {
                detectionProperties mDetectionPropertiesEyeTemp;
                drawVariables mDrawVariablesTemp;
                dataVariables mDataVariablesTemp;

                cv::Mat imageOriginal;

                AOIProperties   eyeAOITemp;
                AOIProperties  beadAOITemp;
                AOIProperties flashAOITemp;

                { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
                    if (!imageCamera.empty())
                    {
                        imageOriginal = imageCamera.clone();

                        mDetectionPropertiesEyeTemp.v = mDetectionVariablesEye;
                        mDrawVariablesTemp = mDrawVariables;
                        mDataVariablesTemp = mDataVariables;

                         eyeAOITemp = Parameters::eyeAOI;
                        beadAOITemp = Parameters::beadAOI;
                       flashAOITemp = flashAOI;

                    } else { return; }
                }

                CamQImage->loadImage(imageOriginal);
                CamQImage->setAOIEye  (  eyeAOITemp);
                CamQImage->setAOIBead ( beadAOITemp);
                CamQImage->setAOIFlash(flashAOITemp);
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
    if (!state)
    {
        imageIndexOffline = 0;

        MainTabWidget->setUpdatesEnabled(false);
        MainTabWidget->insertTab(0, CameraParametersWidget, tr("Camera"));
        MainTabWidget->setUpdatesEnabled(true);

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
        startOfflineSession();
    }
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
    QString text = "Copyright 2016 - 2017 Terence Brouns. <br> <br> "
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
                   "<b>UEye:</b> <br><br> (c) 2016, IDS Imaging Development Systems GmbH <br><br>"
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

void MainWindow::setPupilPosition(double xPos, double yPos)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    if (xPos > 0 && xPos < Parameters::eyeAOI.wdth && yPos > 0 && yPos < Parameters::eyeAOI.hght)
    {
        mDetectionVariablesEye.predictedXPos = xPos;
        mDetectionVariablesEye.predictedYPos = yPos;
    }
}
