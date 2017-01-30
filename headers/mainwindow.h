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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Files

#include "headers/confirmationwindow.h"
#include "headers/constants.h"
#include "headers/drawfunctions.h"
#include "headers/pupildetection.h"
#include "headers/parameters.h"
#include "headers/parameterwidget.h"
#include "headers/sliderdouble.h"
#include "headers/structures.h"
#include "headers/qimageopencv.h"
#include "headers/ueyeopencv.h"
#include "headers/variablewidget.h"

// Standard Template

#include <chrono>
#include <condition_variable> // std::condition_variable
#include <ctime>        // for current date
#include <fstream>      // std::ofstream
#include <iomanip>
#include <iostream>     // std::ofstream
#include <thread>

// Directory creation

#include <sys/types.h>
#include <sys/stat.h>

// Boost

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

// OpenCV

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

// uEye

#include <ueye.h>

// QT

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QColor>
#include <QDebug>
#include <QDesktopWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QMouseEvent>
#include <QPainter>
#include <QPen>
#include <QPixmap>
#include <QProgressBar>
#include <QPushButton>
#include <QRadioButton>
#include <QRect>
#include <QScrollArea>
#include <QSettings>
#include <QSize>
#include <QSlider>
#include <QSpinBox>
#include <QTabWidget>
#include <QTime>
#include <QTimer>
#include <QVBoxLayout>
#include <QWidget>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:

    AOIProperties flashAOI;
    bool APP_EXIT;
    bool APP_RUNNING;
    bool GAIN_AUTO;
    bool GAIN_BOOST;
    bool TRIAL_RECORDING;
    bool FLASH_STANDBY;
    bool OFFLINE_SAVE_DATA;
    bool PROCESSING_ALL_IMAGES;
    bool PROCESSING_ALL_TRIALS;
    bool SAVE_ASPECT_RATIO;
    bool SAVE_CIRCUMFERENCE;
    bool SAVE_POSITION;
    bool SAVE_EYE_IMAGE;
    char currentDate[80];
    cv::Mat imageCamera;
    cv::Mat offlineRawEyeImage;
    double cameraAOIFractionHghtDefaultLeft;
    double cameraAOIFractionWdthDefaultLeft;
    double cameraAOIFractionXPosDefaultLeft;
    double cameraAOIFractionYPosDefaultLeft;
    double cameraAOIFractionHghtDefaultRght;
    double cameraAOIFractionWdthDefaultRght;
    double cameraAOIFractionXPosDefaultRght;
    double cameraAOIFractionYPosDefaultRght;
    double cameraAOIFractionHght;
    double cameraAOIFractionWdth;
    double cameraAOIFractionXPos;
    double cameraAOIFractionYPos;
    double cameraFrameRate;
    double eyeAOIHghtFraction;
    double eyeAOIWdthFraction;
    double beadAOIHghtFraction;
    double beadAOIWdthFraction;
    double flashMinIntensity;
    double guiUpdateFrequency; // refresh frequency of GUI (in Hz)
    double relativeTime;       // in ms
    std::vector<std::vector<double>> timeMatrix;
    detectionMiscellaneous mDetectionMiscellaneousEye;
    detectionMiscellaneous mDetectionMiscellaneousBead;
    detectionVariables mDetectionVariablesEye;
    detectionVariables mDetectionVariablesBead;
    int cameraAOIHghtMax;
    int cameraAOIHghtMin;
    int cameraAOIHghtStepSize;
    int cameraAOIWdthMax;
    int cameraAOIWdthMin;
    int cameraAOIWdthStepSize;
    int cameraSubSamplingFactor;
    int camImageHght;
    int camImageWdth;
    int cameraFrameRateDesired;
    int cameraPixelClock;
    int trialIndexOffline;
    int trialTotalOffline;
    int imageIndexOffline;
    int imageTotalOffline;
    int eyeAOIHghtMin;
    int eyeAOIWdthMin;
    int eyeImageHght;
    int eyeImageWdth;
    int flashThreshold;
    int frameCount;
    int getCurrentTime();
    int trialFrameTotal;
    int trialIndex;
    int trialStartTime;
    int trialTimeLength;
    ParameterWidget *mParameterWidgetBead;
    ParameterWidget *mParameterWidgetEye;
    VariableWidget *mVariableWidgetEye;
    VariableWidget *mVariableWidgetBead;
    QCheckBox *CameraHardwareGainAutoCheckBox;
    QImageOpenCV *CamQImage;
    QImageOpenCV *EyeQImage;
    QLabel *CameraBlackLevelOffsetLabel;
    QLabel *CameraExposureLabel;
    QLabel *CameraFrameRateLabel;
    QLabel *CameraHardwareGainLabel;
    QLabel *CameraPixelClockLabel;
    QLabel *DataAnalysisTitleTextBox;
    QLabel *DataDirectoryTextBox;
    QLabel *FlashStandbyLabel;
    QLabel *FlashThresholdLabel;
    QLabel *OfflineImageFrameTextBox;
    QLineEdit *DataFilenameLineEdit;
    QLineEdit *NameInputLineEdit;
    QLineEdit *TrialTimeLengthLineEdit;
    QScrollArea *BeadTrackingScrollArea;
    QSlider *CameraBlackLevelOffsetSlider;
    QSlider *CameraHardwareGainSlider;
    QSlider *CameraPixelClockSlider;
    QSlider *ExperimentEyeVideoSlider;
    QSlider *FlashThresholdSlider;
    QSlider *FlashStandbySlider;
    QSlider *OfflineImageSlider;
    QSlider *OfflineTrialSlider;
    QSpinBox *CameraFrameRateDesiredSpinBox;
    QSpinBox *OfflineTrialSpinBox;
    QSpinBox *TrialIndexSpinBox;
    QString dataDirectoryOffline;
    QString editSubjectName;
    QString subjectIdentifier;
    QString LastUsedSettingsFileName;
    QTabWidget *MainTabWidget;
    QTimer *UpdateCameraImageTimer;
    QWidget *CameraParametersWidget;
    QWidget *OfflineModeMainWidget;
    QWidget *OfflineModeHeaderWidget;
    SliderDouble *CameraExposureSlider;
    SliderDouble *CameraFrameRateSlider;
    SliderDouble *CamEyeAOIHghtSlider;
    SliderDouble *CamEyeAOIWdthSlider;
    SliderDouble *CamEyeAOIXPosSlider;
    SliderDouble *CamEyeAOIYPosSlider;
    SliderDouble *EyeHghtROISlider;
    SliderDouble *EyeWdthROISlider;
    std::condition_variable cv;
    std::condition_variable saveConditionVariable;
    std::condition_variable cvOffline;
    std::mutex mtx;
    std::mutex mtxOffline;
    std::mutex saveMutex;
    std::string dataDirectory;
    std::string dataFilename;
    std::vector<double> timeStamps;
    std::vector<detectionMiscellaneous> vDetectionMiscellaneousEye;
    std::vector<detectionVariables> vDetectionVariablesBead;
    std::vector<detectionVariables> vDetectionVariablesEye;
    UEyeOpencvCam mUEyeOpencvCam;
    unsigned long long absoluteTime; // in units of 0.1 microseconds
    unsigned long long startTime;
    void countNumTrials();
    void countNumImages();
    void findCamera();
    void getCameraParameters();
    detectionParameters loadParameters(QString, QString, std::vector<double> parameters);
    void loadSettings(QString);
    void msWait(int ms);
    void pupilTracking();
    void resetVariables();
    void detectAllFrames();
    void detectCurrentFrame(int);
    void saveParameters(QString, QString, detectionParameters);
    void saveSettings(QString);
    void saveTrialData();
    void setFlashStandby(bool);
    void setParameterWidgets();
    void setVariableWidgets(const detectionVariables &mDetectionVariablesEye);
    void setupOfflineSession();
    void startTrialRecording();
    void updateCamAOIx();
    void updateCamAOIy();
    void updateEyeAOIx();
    void updateEyeAOIy();
    void updateOfflineTrial();
    void writeToFile(std::ofstream& file, const std::vector<bool>&, const std::vector<double>&, std::string);

protected:

    void keyPressEvent(QKeyEvent *event);

signals:

    void startTimer(int);
    void stopTimer();

public slots:

private slots:

    void onImageNext();
    void onCombineData();
    void onCropAOI();
    void onDetectAllFrames();
    void onDetectAllTrials();
    void onDetectCurrentFrame();
    void onDialogueOpen();
    void onFlashStandbySlider(int);
    void onImagePrevious();
    void onLoadSession();
    void onQuitButtonClicked();
    void onResetFlashIntensity();
    void onResetParameters();
    void onSaveTrialData();
    void onSetOfflineImage(int);
    void onSetOfflineMode(int);
    void onSetTrialIndexOffline(int);
    void onDirectorySelect();
    void onSetAOIEyeLeft();
    void onSetAOIEyeRght();
    void setBeadDetection(int);
    void setCamEyeAOIHght(double);
    void setCamEyeAOIWdth(double);
    void setCamEyeAOIXPos(double);
    void setCamEyeAOIYPos(double);
    void setCameraAutoGain(int);
    void setCameraBlackLevelMode(int);
    void setCameraBlackLevelOffset(int);
    void setCameraExposure(double);
    void setCameraFrameRate(double);
    void setCameraGainBoost(int);
    void setCameraHardwareGain(int);
    void setCameraPixelClock(int);
    void setCameraSubSampling(int);
    void setDrawEdge(int);
    void setDrawElps(int);
    void setDrawHaar(int);
    void setEyeAOIHght(double);
    void setEyeAOIWdth(double);
    void setFlashAOIHght(int);
    void setFlashAOIWdth(int);
    void setFlashAOIXPos(int);
    void setFlashAOIYPos(int);
    void setFlashThreshold(int);
    void setPupilPosition(double, double);
    void onSetRealTime(int);
    void setSaveDataAspectRatio(int);
    void setSaveDataCircumference(int);
    void setSaveDataPosition(int);
    void onSetTrialIndex(int);
    void startOfflineSession();
    void startRecordingManual();
    void updateCameraImage();
    void updateImageRaw(int);
    void updateImageProcessed(int);

};

#endif // MAINWINDOW_H
