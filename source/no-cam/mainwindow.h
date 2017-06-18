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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Files

#include "../confirmationwindow.h"
#include "../constants.h"
#include "../drawfunctions.h"
#include "../eyestalker.h"
#include "../parameters.h"
#include "../parameterwidget.h"
#include "../sliderdouble.h"
#include "../structures.h"
#include "../qimageopencv.h"
#include "../variablewidget.h"

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
#include <boost/range/iterator_range.hpp>

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

    bool OFFLINE_SAVE_DATA;
    bool PROCESSING_ALL_IMAGES;
    bool PROCESSING_ALL_TRIALS;
    bool PROCESSING_ALL_EXPS;

    double guiUpdateFrequency; // refresh frequency of GUI (in Hz)
    std::vector<std::vector<double>> timeMatrix;
    int camImageHght;
    int camImageWdth;
    int trialIndexOffline;
    int trialTotalOffline;
    int imageIndexOffline;
    int imageTotalOffline;
    int getCurrentTime();

    QImageOpenCV *CamQImage;

    QLabel      *OfflineImageFrameTextBox;
    QSlider     *OfflineImageSlider;
    QSlider     *OfflineTrialSlider;
    QSpinBox    *OfflineTrialSpinBox;

    QString dataDirectoryOffline;

    QDoubleSpinBox *CameraFrameRateSpinBox;
//    SliderDouble  *CameraFrameRateSpinBox;
    double         cameraFrameRate;

    QString LastUsedSettingsFileName;
    QTabWidget *MainTabWidget;

    std::condition_variable cvOffline;
    std::mutex mutexOffline;

    // Options menu

    // Interface

    QScrollArea *AdvancedScrollArea;

    // AOIs

    void updateAOIx();
    void updateAOIy();

    // Variables and parameters

    void resetVariablesHard(detectionVariables&, const detectionParameters&, const AOIProperties&);
    void resetVariablesSoft(detectionVariables&, const detectionParameters&, const AOIProperties&);

    detectionVariables  mDetectionVariablesEye;
    drawVariables       mDrawVariablesEye;
    dataVariables       mDataVariablesEye;

    std::vector<detectionVariables> vDetectionVariablesEye;
    std::vector<dataVariables>      vDataVariablesEye;

    ParameterWidget *mParameterWidgetEye;
    VariableWidget  *mVariableWidgetEye;

    // Saving and loading settings

    detectionParameters loadParameters(QString, QString, std::vector<double> parameters);
    void saveParameters(QString, QString, detectionParameters);
    void loadSettings(QString);
    void saveSettings(QString);

    // Experimental

    bool SAVE_ASPECT_RATIO;
    bool SAVE_CIRCUMFERENCE;
    bool SAVE_POSITION;
    bool SAVE_EYE_IMAGE;

    // Offline interface

    QWidget *OfflineModeWidget;

    // Offline functions and variables

    void countNumTrials();
    void countNumImages();

    void setupOfflineSession();
    void updateOfflineTrial();

    void detectCurrentFrame(int);
    void detectAllFrames();

    void setCurvatureMeasurement(detectionParameters&, int);

    // Advanced

    developmentOptions mAdvancedOptions;

    bool SAVE_DATA_EDGE;
    bool SAVE_DATA_FIT;
    bool SAVE_DATA_EXTRA;

    // General

    void msWait(int ms);

protected:

signals:

public slots:

private slots:

    void onCombineData              ();
    void onDetectAllExperiments     ();
    void onDetectAllFrames          ();
    void onDetectAllTrials          ();
    void onDetectCurrentFrame       ();
    void onDialogueOpen             ();
    void onImageNext                ();
    void onImagePrevious            ();
    void onLoadSession              ();
    void onQuitButtonClicked        ();
    void onResetParameters          ();
    void onSaveTrialData            ();
    void onSetAdvancedMode          (bool);
    void onSetCameraFrameRate       (double);
    void onSetDrawEdge              (int);
    void onSetDrawElps              (int);
    void onSetDrawHaar              (int);
    void onSetOfflineImage          (int);
    void onSetTrialOffline          (int);
    void onUpdateImageProcessed     (int);
    void onUpdateImageRaw           (int);

    // Advanced

    void onSetCurvatureMeasurement  (int);

    void onSetSaveDataEdge (int);
    void onSetSaveDataFit  (int);
    void onSetSaveDataExtra(int);

};

#endif // MAINWINDOW_H
