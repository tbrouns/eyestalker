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

// Files

#include "headers/confirmationwindow.h"
#include "headers/constants.h"
#include "headers/drawfunctions.h"
#include "headers/eyetracking.h"
#include "headers/parameters.h"
#include "headers/sliderdouble.h"
#include "headers/structures.h"
#include "headers/qimageopencv.h"
#include <headers/ueyeopencv.h>

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

    bool APP_EXIT;
    bool APP_RUNNING;
    bool GAIN_AUTO;
    bool GAIN_BOOST;
    bool EXPERIMENT_TRIAL_RECORDING;
    bool FLASH_STANDBY;
    bool PROCESSING_ALL_IMAGES;
    bool PROCESSING_ALL_TRIALS;
    bool SAVE_EYE_IMAGE;
    char currentDate[80];
    cv::Mat imageCamera;
    cv::Mat reviewRawEyeImage;
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
    double circumferenceOffset;
    double eyeAOIHghtFraction;
    double eyeAOIWdthFraction;
    double flashThreshold;
    double guiUpdateFrequency; // refresh frequency of GUI (in Hz)
    double relativeTime; // in ms
    std::vector<double> timeVector;
    eyePropertiesMiscellaneous mEyePropertiesMiscellaneous;
    eyePropertiesParameters    mEyePropertiesParameters;
    eyePropertiesVariables     mEyePropertiesVariables;
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
    int editDataIndex;
    int editDataTotal;
    int editImageIndex;
    int editImageTotal;
    int experimentIndex;
    int eyeAOIHghtMin;
    int eyeAOIWdthMin;
    int eyeImageHght;
    int eyeImageWdth;
    int frameCount;
    int getCurrentTime();
    int trialFrameTotal;
    int trialIndex;
    int trialStartTime;
    int trialTimeLength;
    QCheckBox *CameraHardwareGainAutoCheckBox;
    QImageOpenCV *CamQImage;
    QImageOpenCV *EyeQImage;
    QLabel *AlphaGeneralLabel;
    QLabel *AlphaMomentumLabel;
    QLabel *AlphaPupilLabel;
    QLabel *AlphaAverageLabel;
    QLabel *AlphaPredictionLabel;
    QLabel *AlphaVelocityLabel;
    QLabel *CameraBlackLevelOffsetLabel;
    QLabel *CameraExposureLabel;
    QLabel *CameraFrameRateLabel;
    QLabel *CameraHardwareGainLabel;
    QLabel *CameraPixelClockLabel;
    QLabel *CannyBlurLevelLabel;
    QLabel *CannyKernelSizeLabel;
    QLabel *CannyLowerLimitLabel;
    QLabel *CannyUpperLimitLabel;
    QLabel *CurvatureOffsetLabel;
    QLabel *DataAnalysisTitleTextBox;
    QLabel* DataDirectoryTextBox;
    QLabel *EdgeIntensityLabel;
    QLabel *EdgeIntensityOffsetLabel;
    QLabel *EdgeMaximumFitNumberLabel;
    QLabel *EllipseFitErrorMaximumLabel;
    QLabel *FlashStandbyLabel;
    QLabel *FlashThresholdLabel;
    QLabel *GlintRadiusLabel;
    QLabel *PupilCircfLabel;
    QLabel *PupilCircumferenceIniLabel;
    QLabel *PupilCircumferenceOffLabel;
    QLabel *PupilFractionIniLabel;
    QLabel *PupilFractionMinLabel;
    QLabel *PupilFractLabel;
    QLabel *PupilHaarOffsetLabel;
    QLabel *ReviewImageFrameTextBox;
    QLabel *ThresholdCircumferenceLabel;
    QLabel *ThresholdFractionLabel;
    QLineEdit *DataFilenameLineEdit;
    QLineEdit *NameInputLineEdit;
    QLineEdit *TrialTimeLengthLineEdit;
    QSlider *CameraBlackLevelOffsetSlider;
    QSlider *CameraHardwareGainSlider;
    QSlider *CameraPixelClockSlider;
    QSlider *CannyBlurLevelSlider;
    QSlider *CannyKernelSizeSlider;
    QSlider *CannyLowerLimitSlider;
    QSlider *CannyUpperLimitSlider;
    QSlider *EdgeMaximumFitNumberSlider;
    QSlider *ExperimentEyeVideoSlider;
    QSlider *FlashStandbySlider;
    QSlider *GlintRadiusSlider;
    QSlider *ReviewImageSlider;
    QSlider *ReviewTrialSlider;
    QSpinBox *CameraFrameRateDesiredSpinBox;
    QSpinBox *ReviewTrialSpinBox;
    QSpinBox *TrialIndexSpinBox;
    QString editDataDirectory;
    QString editSubjectName;
    QString subjectIdentifier;
    QString LastUsedSettingsFileName;
    QTabWidget *EyeTrackingParameterTabWidget;
    QTimer *UpdateCameraImageTimer;
    QWidget *CameraParametersWidget;
    QWidget *EyeTrackingReviewWidget;
    SliderDouble *AlphaGeneralSlider;
    SliderDouble *AlphaMomentumSlider;
    SliderDouble *AlphaPupilSlider;
    SliderDouble *AlphaAverageSlider;
    SliderDouble *AlphaPredictionSlider;
    SliderDouble *AlphaVelocitySlider;
    SliderDouble *CameraExposureSlider;
    SliderDouble *CameraFrameRateSlider;
    SliderDouble *CamEyeAOIHghtSlider;
    SliderDouble *CamEyeAOIWdthSlider;
    SliderDouble *CamEyeAOIXPosSlider;
    SliderDouble *CamEyeAOIYPosSlider;
    SliderDouble *CurvatureOffsetSlider;
    SliderDouble *EdgeIntensityOffsetSlider;
    SliderDouble *EdgeIntensitySliderDouble;
    SliderDouble *EllipseFitErrorMaximumSlider;
    SliderDouble *EyeHghtROISlider;
    SliderDouble *EyeWdthROISlider;
    SliderDouble *PupilCircfSliderDouble;
    SliderDouble *PupilCircumferenceIniSlider;
    SliderDouble *PupilCircumferenceOffSlider;
    SliderDouble *PupilFractionIniSlider;
    SliderDouble *PupilFractionMinSlider;
    SliderDouble *PupilFractSliderDouble;
    SliderDouble *PupilHaarOffsetSlider;
    SliderDouble *ThresholdCircumferenceSlider;
    SliderDouble *ThresholdFractionSlider;
    std::condition_variable cv;
    std::condition_variable saveConditionVariable;
    std::mutex mtx;
    std::mutex reviewMutex;
    std::mutex saveMutex;
    std::string dataDirectory;
    std::string dataFilename;
    std::vector<bool> eyeDetectionFlags;
    std::vector<double> eyeXPositions;
    std::vector<double> eyeYPositions;
    std::vector<double> timeStamps;
    std::vector<eyePropertiesVariables> vEyePropertiesVariables;
    UEyeOpencvCam mUEyeOpencvCam;
    unsigned long long absoluteTime; // in units of 0.1 microseconds
    unsigned long long startTime;
    void findCamera();
    void getCameraParameters();
    void loadSettings(QString);
    void msWait(int ms);
    void pupilTracking();
    void reviewPupilDetectionAllFrames();
    void reviewPupilDetectionOneFrame();
    void reviewSaveExperimentData();
    void saveSettings(QString);
    void saveTrialData();
    void setFlashStandby(bool);
    void setParameterWidgets();
    void setVariableWidgets(const eyePropertiesVariables &);
    void startTrialRecording();
    void updateCamAOIx();
    void updateCamAOIy();
    void updateEyeAOIx();
    void updateEyeAOIy();
    void updateReviewImages(int);
    void updateReviewSession();
    void writeToFile(std::ofstream& file, const std::vector<bool>&, const std::vector<double>&);

protected:

    void keyPressEvent(QKeyEvent *event);

signals:

    void startTimer(int);
    void stopTimer();

public slots:

private slots:

    void changeReviewSession(int index);
    void cropAOI();
    void detectPupilAllFrames();
    void detectPupilAllTrials();
    void detectPupilOneFrame();
    void loadReviewSession();
    void nextReviewImage();
    void onQuitButtonClicked();
    void onSavePupilData();
    void onFlashStandbySlider(int);
    void prevReviewImage();
    void selectDirectory();
    void setAlphaGeneral(double value);
    void setAlphaMomentum(double value);
    void setAlphaPupil(double value);
    void setAlphaAverage(double value);
    void setAlphaPrediction(double value);
    void setAlphaVelocity(double value);
    void setCameraAutoGain(int state);
    void setCameraBlackLevelMode(int state);
    void setCameraBlackLevelOffset(int value);
    void setCameraExposure(double value);
    void setCameraFrameRate(double value);
    void setCameraGainBoost(int state);
    void setCameraHardwareGain(int value);
    void setCameraPixelClock(int value);
    void setCameraSubSampling(int state);
    void setCamEyeAOIHght(double);
    void setCamEyeAOIWdth(double);
    void setCamEyeAOIXPos(double);
    void setCamEyeAOIYPos(double);
    void setCannyBlurLevel(int value);
    void setCannyKernelSize(int value);
    void setCannyLowerLimit(int value);
    void setCannyUpperLimit(int value);
    void setCurvatureOffset(double value);
    void setDrawEdge(int state);
    void setDrawElps(int state);
    void setDrawHaar(int state);
    void setEdgeIntensity(double);
    void setEdgeIntensityOffset(double value);
    void setEdgeMaximumFitNumber(int value);
    void setEllipseFitErrorMaximum(double value);
    void setEyeROIHght(double);
    void setEyeROIWdth(double);
    void setFlashAOIXPos(int);
    void setFlashAOIYPos(int);
    void setFlashAOIWdth(int);
    void setFlashAOIHght(int);
    void setFlashThreshold(double);
    void setGlintRadius(int value);
    void setAOILeftEye();
    void setPupilCircumference(double);
    void setPupilCircumferenceIni(double value);
    void setPupilCircumferenceOff(double value);
    void setPupilFraction(double);
    void setPupilFractionIni(double);
    void setPupilFractionMin(double);
    void setPupilHaarOffset(double value);
    void setPupilPosition(double xPos, double yPos);
    void setRealTimeEyeTracking(int state);
    void setReviewMode(int state);
    void setReviewImageFrame(int);
    void setAOIRghtEye();
    void setThresholdCircumference(double);
    void setThresholdFraction(double);
    void setTrialIndex(int);
    void startReviewSession();
    void updateCameraImage();
    void updateRawImage();

};

#endif // MAINWINDOW_H
