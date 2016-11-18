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
        emit stopTimer();

        EXPERIMENT_TRIAL_RECORDING = true;

        mUEyeOpencvCam.startRecording();

        absoluteTime = startTime;
        relativeTime = 0;

        trialStartTime = getCurrentTime();
        trialFrameTotal = ceil((trialTimeLength * cameraFrameRate) / 1000);

        trialTimeLength = (TrialTimeLengthLineEdit->text()).toInt();

        frameCount = 0;

        timeStamps.assign(trialFrameTotal, 0);
        eyeXPositions.assign(trialFrameTotal, 0);
        eyeYPositions.assign(trialFrameTotal, 0);
        eyeDetectionFlags.assign(trialFrameTotal, 0);
    }

    if (SAVE_EYE_IMAGE)
    {
        std::stringstream directoryName;
        directoryName << dataDirectory << "/" << currentDate;

        if (!boost::filesystem::exists(directoryName.str()))
        {
            boost::filesystem::create_directory(directoryName.str().c_str());
        }

        directoryName << "/" << (NameInputLineEdit->text()).toStdString();

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

void MainWindow::writeToFile(std::ofstream& file, const std::vector<bool>& vFlags, const std::vector<double>& v, std::string delimiter)
{
    int paddingTotal = 0;

    for (int i = 0, vSize = v.size(); i < vSize; i++) // write data
    {
        if (vFlags[i]) // only detected pupil data
        {
            file << v[i] << delimiter;
        }
        else
        {
            paddingTotal++;
        }
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
                 << "/" << dataFilename;

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

        file << std::setw(3) << std::setfill('0') << trialIndex << " "; // print with leading zeros
        file << trialStartTime << " "; // time of day in milliseconds

        file << std::fixed;
        file << std::setprecision(3);

        writeToFile(file, eyeDetectionFlags, timeStamps, ";");
        writeToFile(file, eyeDetectionFlags, eyeXPositions, ";");
        writeToFile(file, eyeDetectionFlags, eyeYPositions, ";");

        file.close();
    }
    else
    {
        // save info in .ini file

        std::stringstream saveFileNameSS;
        saveFileNameSS << dataDirectory
                       << "/" << currentDate
                       << "/" << (NameInputLineEdit->text()).toStdString()
                       << "/trial_" << trialIndex;

        editDataDirectory = QString::fromStdString(saveFileNameSS.str());

        saveFileNameSS << "/info.ini";

        editDataIndex = trialIndex;

        editImageTotal = frameCount;

        saveSettings(QString::fromStdString(saveFileNameSS.str()));

        // save time stamps

        std::stringstream filename;
        filename << dataDirectory << "/"
                 << currentDate << "/"
                 << "timestamps.dat";

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

        file << std::setw(3) << std::setfill('0') << trialIndex << " "; // print with leading zeros
        file << trialStartTime << " "; // time of day in milliseconds

        file << std::fixed;
        file << std::setprecision(3);

        writeToFile(file, eyeDetectionFlags, timeStamps, " ");

        file.close();
    }
}
