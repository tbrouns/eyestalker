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
    // Initialize default values

    Parameters::ellipseDrawCrossSize    = 5;
    Parameters::ellipseDrawOutlineWidth = 0.032;

    guiUpdateFrequency = 30;

    mParameterWidgetEye  = new ParameterWidget;

    mVariableWidgetEye  = new VariableWidget;

    Parameters::eyeAOIRatio.xPos = 0.0;
    Parameters::eyeAOIRatio.yPos = 0.0;
    Parameters::eyeAOIRatio.hght = 1.0;
    Parameters::eyeAOIRatio.wdth = 1.0;

    // AOI

    camImageHght = 200;
    camImageWdth = 400; // size of image in widget

    // Offline mode

    PROCESSING_ALL_IMAGES = false;
    PROCESSING_ALL_TRIALS = false;
    PROCESSING_ALL_EXPS   = false;

    trialIndexOffline = 0;
    imageIndexOffline = 0;
    imageTotalOffline = 0;

    // Advanced options

    mAdvancedOptions.CURVATURE_MEASUREMENT = false;

    SAVE_DATA_FIT   = false;
    SAVE_DATA_EDGE  = false;
    SAVE_DATA_EXTRA = false;

    // Grab parameters from ini file

    LastUsedSettingsFileName = "config_user.ini";

    loadSettings(LastUsedSettingsFileName);

    // Camera feed

    CamQImage = new QImageOpenCV();
    CamQImage->setSize(camImageWdth, camImageHght);
    CamQImage->setAOIEye(Parameters::eyeAOI);
    CamQImage->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    QObject::connect(CamQImage, SIGNAL(updateImage(int)), this , SLOT(onUpdateImageRaw(int)));

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

    QPushButton *ResetPushButton = new QPushButton("Reset parameters");
    QObject::connect(ResetPushButton, SIGNAL(clicked()), this, SLOT(onResetParameters()));

    QLabel *CameraFrameRateTextBox = new QLabel;
    CameraFrameRateTextBox->setText("<b>Frame-rate (Hz):</b>");

    CameraFrameRateSpinBox = new QDoubleSpinBox;
    CameraFrameRateSpinBox->setDecimals(1);
    CameraFrameRateSpinBox->setRange(30, 1000);
    CameraFrameRateSpinBox->setValue(cameraFrameRate);

    QObject::connect(CameraFrameRateSpinBox, SIGNAL(doubleValueChanged(double)), this, SLOT(onSetCameraFrameRate(double)));

    QHBoxLayout* CameraFrameRateLayout = new QHBoxLayout;
    CameraFrameRateLayout->addWidget(CameraFrameRateTextBox);
    CameraFrameRateLayout->addWidget(CameraFrameRateSpinBox);

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

    QWidget *CameraSettings = new QWidget;
    QGridLayout* CameraOutputLayout = new QGridLayout(CameraSettings);

    CameraOutputLayout->addWidget(CamQImage,             0, 1, Qt::AlignCenter);
    CameraOutputLayout->addWidget(OfflineModeWidget,     1, 1, Qt::AlignCenter);
    CameraOutputLayout->addLayout(CameraFrameRateLayout, 2, 1 ,Qt::AlignCenter);
    CameraOutputLayout->addLayout(DrawFunctionsLayout,   3, 1, Qt::AlignCenter);
    CameraOutputLayout->addWidget(ResetPushButton,       4, 1, Qt::AlignCenter);

    CameraOutputLayout->setColumnStretch(0, 1);
    CameraOutputLayout->setColumnStretch(3, 1);

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

    AdvancedScrollArea = new QScrollArea();
    AdvancedScrollArea->setWidget(AdvancedOptionsWidget);
    AdvancedScrollArea->setWidgetResizable(true);

    MainTabWidget = new QTabWidget;
    MainTabWidget->addTab(EyeTrackingScrollArea,   tr("Eye tracking"));
    MainTabWidget->addTab(AdvancedScrollArea,      tr("Advanced"));

    MainTabWidget->setTabEnabled(1, 0); // disable development mode

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
}

MainWindow::~MainWindow()
{

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

    cameraFrameRate = 250;

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
    saveSettings(LastUsedSettingsFileName);
    qApp->quit();
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

void MainWindow::updateAOIx()
{
    { std::lock_guard<std::mutex> AOIEyeLock(Parameters::AOIEyeMutex);
        Parameters::eyeAOI.wdth = round(Parameters::camAOI.wdth * Parameters::eyeAOIRatio.wdth);
        Parameters::eyeAOI.xPos = round(Parameters::camAOI.wdth * Parameters::eyeAOIRatio.xPos);
        if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::camAOI.wdth)
        {   Parameters::eyeAOI.xPos = Parameters::camAOI.wdth - Parameters::eyeAOI.wdth; }
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
        else
        {
            if (imageTotalOffline == 0)
            {
                cv::Mat imageRaw = cv::imread(filename.str(), CV_LOAD_IMAGE_COLOR);
                Parameters::eyeAOI.wdth = imageRaw.cols;
                Parameters::eyeAOI.hght = imageRaw.rows;
            }

            imageTotalOffline++;
        }
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

            vDetectionVariablesEye.resize( imageTotalOffline + 1);

            resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);

            vDetectionVariablesEye [0] = mDetectionVariablesEye;

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
        mVariableWidgetEye ->setWidgets(vDataVariablesEye[imgIndex]);
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
    detectionParameters mDetectionParametersEyeTemp;

    AOIProperties AOIEyeTemp;

    { std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

        mDetectionVariablesEyeTemp  = vDetectionVariablesEye[imageIndex];
        mDetectionParametersEyeTemp = mParameterWidgetEye->getStructure();

        mDetectionParametersEyeTemp.cameraFrameRate = cameraFrameRate;

        Parameters::camAOI.wdth = imageRaw.cols;
        Parameters::camAOI.hght = imageRaw.rows;

        if (mAdvancedOptions.CURVATURE_MEASUREMENT) { setCurvatureMeasurement(mDetectionParametersEyeTemp, imageRaw.cols); }

        updateAOIx();
        updateAOIy();

        AOIEyeTemp  = Parameters::eyeAOI;
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

    mDataVariablesEye.duration      = fp_ms.count();
    mDataVariablesEye.absoluteXPos  = mDataVariablesEye.exactXPos;
    mDataVariablesEye.absoluteYPos  = mDataVariablesEye.exactYPos;
    vDataVariablesEye[imageIndex]   = mDataVariablesEye;

    cv::Mat imageProcessed = imageRaw.clone();
    drawAll(imageProcessed, mDrawVariablesEye);

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
        if (timeMatrix.size() > 0) { file << (int) timeMatrix[trialIndexOffline][1] << ";"; } // system clock time

        file << std::fixed;
        file << std::setprecision(3);

        // write data

        for (int i = 0; i < imageTotalOffline; i++)                              { file << vDataVariablesEye[i].DETECTED        << delimiter; }
        if (timeMatrix.size() > 0) { for (int i = 0; i < imageTotalOffline; i++) { file << timeMatrix[trialIndexOffline][i + 2] << delimiter; } }

        for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].absoluteXPos         << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].absoluteYPos         << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].exactCircumference   << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].exactAspectRatio     << delimiter; }

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
            for (int i = 0; i < imageTotalOffline; i++) { file << vDataVariablesEye[i].duration                    << delimiter; }
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
    std::stringstream fileNameWriteSS;
    fileNameWriteSS << dataDirectoryOffline.toStdString()
                    << "/combined_data.dat";

    std::string fileNameWrite = fileNameWriteSS.str();

    if (boost::filesystem::exists(fileNameWrite))
    {
        QString text = "The file <b>combined_data.dat</b> already exists in <b>" + dataDirectoryOffline + "/</b>. Do you wish to add data to the end of this file?";
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
    QSettings settings(filename, QSettings::IniFormat);

    dataDirectoryOffline            = settings.value("DataDirectoryOffline",       "").toString();
    Parameters::drawFlags.haar      = settings.value("DrawHaar",                false).toBool();
    Parameters::drawFlags.edge      = settings.value("DrawEdge",                false).toBool();
    Parameters::drawFlags.elps      = settings.value("DrawElps",                 true).toBool();

    detectionParameters mDetectionParametersEye  = loadParameters(filename, "Eye",  parametersEye);
    mParameterWidgetEye ->setStructure(mDetectionParametersEye);

    resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye->getStructure(), Parameters::eyeAOI);
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
    cameraFrameRate                                         = settings.value(prefix + "CameraFrameRate",                           250).toDouble();

    return mDetectionParameters;
}

void MainWindow::saveSettings(QString filename)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue("DataDirectoryOffline",   dataDirectoryOffline);
    settings.setValue("DrawHaar",               Parameters::drawFlags.haar);
    settings.setValue("DrawEdge",               Parameters::drawFlags.edge);
    settings.setValue("DrawElps",               Parameters::drawFlags.elps);

    detectionParameters mDetectionParametersEye  = mParameterWidgetEye->getStructure();
    saveParameters(filename,  "Eye", mDetectionParametersEye);
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
    settings.setValue(prefix + "CameraFramerRate",          cameraFrameRate);
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
    resetVariablesHard(mDetectionVariablesEye,  mParameterWidgetEye ->getStructure(), Parameters::eyeAOI);
}

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

void MainWindow::onSetCameraFrameRate(double val) { cameraFrameRate = val; }
