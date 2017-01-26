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
        if (!boost::filesystem::exists(dataDirectory))
        {
            QString text = "Data directory does not exist. Please choose a valid path.";
            ConfirmationWindow mConfirmationWindow(text, false);
            mConfirmationWindow.setWindowTitle("Warning");
            mConfirmationWindow.exec();
        }

        emit stopTimer();

        TRIAL_RECORDING = true;

        mUEyeOpencvCam.startRecording();

        absoluteTime = startTime;
        relativeTime = 0;

        trialStartTime  = getCurrentTime();
        trialFrameTotal = ceil((trialTimeLength * cameraFrameRate) / 1000);

        trialTimeLength = (TrialTimeLengthLineEdit->text()).toInt();

        frameCount = 0;

        vDetectionVariablesEye.resize(trialFrameTotal);
        vDetectionVariablesBead.resize(trialFrameTotal);

        timeStamps.assign(trialFrameTotal, 0);

        if (SAVE_EYE_IMAGE)
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

void MainWindow::writeToFile(std::ofstream& file, const std::vector<bool>& vFlags, const std::vector<double>& v, std::string delimiter)
{
    int paddingTotal = 0;

    for (int i = 0, vSize = v.size(); i < vSize; i++) // write data
    {
        if (vFlags[i]) { file << v[i] << delimiter; }
        else           { paddingTotal++; }
    }

    for (int i = 0; i < paddingTotal; i++) // write padding of zeroes
    {
        file << 0 << delimiter;
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

        for (int i = 0; i < frameCount; i++) { file << vDetectionVariablesEye[i].pupilDetected << delimiter; }
        for (int i = 0; i < frameCount; i++) { file << timeStamps[i] << delimiter; } // saving ALL timestamps allows for accurate frame-rate calculation during data processing

        if (SAVE_POSITION)
        {
            for (int i = 0; i < frameCount; i++) { file << vDetectionVariablesEye[i].xPosAbsolute  << delimiter; }
            for (int i = 0; i < frameCount; i++) { file << vDetectionVariablesEye[i].yPosAbsolute  << delimiter; }
        }

        if (SAVE_CIRCUMFERENCE)
        {
            for (int i = 0; i < frameCount; i++) { file << vDetectionVariablesEye[i].circumferenceExact << delimiter; }
        }

        if (SAVE_ASPECT_RATIO)
        {
            for (int i = 0; i < frameCount; i++) { file << vDetectionVariablesEye[i].aspectRatioExact << delimiter; }
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

        for (int i = 0; i < frameCount; i++) { file << timeStamps[i] << delimiter; }
        
        file.close();
    }
}
