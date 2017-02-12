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

void MainWindow::startRecordingManual()
{
    if (!TRIAL_RECORDING && !PROCESSING_ALL_TRIALS && !PROCESSING_ALL_IMAGES)
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

void MainWindow::plotTrialData()
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

