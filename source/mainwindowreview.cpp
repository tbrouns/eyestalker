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

void MainWindow::startReviewSession()
{
    editImageIndex = 0;

    Parameters::REALTIME_PROCESSING = false; // turn off pupil tracking

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

    EyeTrackingReviewWidget->setVisible(true);

    CamQImage->clearImage();
    EyeQImage->clearImage();

    // hide tab with camera parameters
    EyeTrackingParameterTabWidget->removeTab(0);

    updateReviewImages(0);
}

void MainWindow::loadReviewSession()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Please select an info.ini file"), QString::fromStdString(dataDirectory), tr("Ini file (*.ini)"));

    if (!fileName.isEmpty())
    {
        int prevSubSampling = cameraSubSamplingFactor;

        loadSettings(fileName);

        int nextSubSampling = cameraSubSamplingFactor;
        cameraSubSamplingFactor = prevSubSampling;
        setCameraSubSampling(nextSubSampling - 1);

        editDataDirectory = QFileInfo(fileName).absolutePath();

        updateReviewSession();
        updateReviewImages(0);
    }
}

void MainWindow::prevReviewSession()
{
    int index = editDataIndex - 1;

    if (index >= 0)
    {
        QString infoFilePath = editDataDirectory;

        int numOfDigits;

        if (editDataIndex != 0)
        {
            numOfDigits = floor(log10(editDataIndex)) + 1;
        }
        else
        {
            numOfDigits = 1;
        }

        infoFilePath.replace(infoFilePath.size() - numOfDigits, numOfDigits, QString::number(index)); // replace last character(s)

        QString infoFileDataDirectory = infoFilePath;
        infoFilePath.append("/info.ini");

        if (boost::filesystem::exists(infoFilePath.toStdString()))
        {
            loadSettings(infoFilePath);
            editDataDirectory = infoFileDataDirectory;

            updateReviewSession();
            updateReviewImages(0);
        }
    }
}

void MainWindow::nextReviewSession()
{
    int index = editDataIndex + 1;

    QString infoFilePath = editDataDirectory;

    int numOfDigits;

    if (editDataIndex != 0)
    {
        numOfDigits = floor(log10(editDataIndex)) + 1;
    }
    else
    {
        numOfDigits = 1;
    }

    infoFilePath.replace(infoFilePath.size() - numOfDigits, numOfDigits, QString::number(index)); // replace last character(s)

    QString infoFileDataDirectory = infoFilePath;
    infoFilePath.append("/info.ini");

    if (boost::filesystem::exists(infoFilePath.toStdString()))
    {
        loadSettings(infoFilePath);
        editDataDirectory = infoFileDataDirectory;

        updateReviewSession();
        updateReviewImages(0);
    }
}

void MainWindow::updateReviewSession()
{
    if (editImageTotal > 0)
    {
        setParameterWidgets();

        std::stringstream ss;
        ss << "<b>SubjectID: </b>" << "<i>" << editSubjectName.toStdString() << "</i><b> - Trial " << editDataIndex + 1 << "</b>";
        QString title = QString::fromStdString(ss.str());
        ReviewSessionTitleTextBox->setText(title);

        vEyePropertiesVariables.resize(editImageTotal + 1);
        vEyePropertiesVariables[0] = mEyePropertiesVariables;

        editImageIndex = 0;

        ReviewImageSlider->setValue(0);
        ReviewImageSlider->setMaximum(editImageTotal - 1);

        // Create folder

        std::stringstream directoryName;
        directoryName << editDataDirectory.toStdString() << "/images/processed/";

        if (!boost::filesystem::exists(directoryName.str()))
        {
            mkdir(directoryName.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
    }
}

void MainWindow::setPupilPosition(double xPos, double yPos)
{
    std::lock_guard<std::mutex> secondaryMutexLock(Parameters::secondaryMutex);

    if (xPos > 0 && xPos < Parameters::eyeAOIWdth && yPos > 0 && yPos < Parameters::eyeAOIHght)
    {
        int pupilHaarWdth = round(mEyePropertiesVariables.pupilCircumference / M_PI);
        int pupilHaarWdthOffset = pupilHaarWdth + round(pupilHaarWdth * mEyePropertiesParameters.pupilOffset * 2);

        mEyePropertiesVariables.searchRadius = ceil(0.5 * pupilHaarWdthOffset);
        mEyePropertiesVariables.searchXPos = xPos;
        mEyePropertiesVariables.searchYPos = yPos;

        if (editImageIndex == 0)
        {
            mEyePropertiesVariables.xPos = xPos;
            mEyePropertiesVariables.yPos = yPos;
        }
    }
}

void MainWindow::updateReviewImages(int imgIndex)
{
    if (editImageTotal > 0)
    {
        std::stringstream ss;
        ss << "<b>" << imgIndex + 1 << " / " << editImageTotal << "</b>";
        QString title = QString::fromStdString(ss.str());
        ReviewImageFrameTextBox->setText(title);

        // Raw

        std::stringstream fileNameRaw;
        fileNameRaw << editDataDirectory.toStdString()
                    << "/images/raw/" << editSubjectName.toStdString()
                    << "_" << editDataIndex << "_" << imgIndex << ".png";

        if (boost::filesystem::exists(fileNameRaw.str()))
        {
            cv::Mat rawEyeImage = cv::imread(fileNameRaw.str(), CV_LOAD_IMAGE_COLOR);

            CamQImage->loadImage(rawEyeImage);
            CamQImage->setEyeAOI(Parameters::eyeAOIXPos, Parameters::eyeAOIYPos, Parameters::eyeAOIWdth, Parameters::eyeAOIHght);
            CamQImage->setImage();
        }
        else
        {
            CamQImage->clearImage();
        }

        std::stringstream fileName;
        fileName << editDataDirectory.toStdString()
                 << "/images/processed/" << editSubjectName.toStdString()
                 << "_" << editDataIndex << "_" << imgIndex << ".png";

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
    else
    {
        CamQImage->clearImage();
        EyeQImage->clearImage();

        ReviewImageFrameTextBox->setText("<b>0 / 0</b>");
    }
}

void MainWindow::updateRawImage()
{
    std::stringstream fileNameRaw;
    fileNameRaw << editDataDirectory.toStdString()
                << "/images/raw/" << editSubjectName.toStdString()
                << "_" << editDataIndex << "_" << editImageIndex << ".png";

    if (boost::filesystem::exists(fileNameRaw.str()))
    {
        cv::Mat rawEyeImage = cv::imread(fileNameRaw.str(), CV_LOAD_IMAGE_COLOR);

        CamQImage->loadImage(rawEyeImage);
        CamQImage->setEyeAOI(Parameters::eyeAOIXPos, Parameters::eyeAOIYPos, Parameters::eyeAOIWdth, Parameters::eyeAOIHght);
        CamQImage->setImage();
    }
    else
    {
        CamQImage->clearImage();
    }
}

void MainWindow::prevReviewImage()
{
    if (editImageIndex > 0)
    {
        editImageIndex--;
        ReviewImageSlider->setValue(editImageIndex);
    }
}

void MainWindow::nextReviewImage()
{
    if (editImageIndex < editImageTotal - 1)
    {
        editImageIndex++;
        ReviewImageSlider->setValue(editImageIndex);
    }
}

void MainWindow::setReviewImageFrame(int frame)
{
    if (!DETECTION_IN_PROGRESS)
    {
        editImageIndex = frame;
        setVariableWidgets(vEyePropertiesVariables[editImageIndex]);
        updateReviewImages(editImageIndex);
    }
}

void MainWindow::reviewPupilDetectionOneFrame()
{
    cv::Mat rawEyeImage;

    std::stringstream fileNameRaw;
    fileNameRaw << editDataDirectory.toStdString()
                << "/images/raw/" << editSubjectName.toStdString()
                << "_" << editDataIndex << "_" << editImageIndex << ".png";

    if (boost::filesystem::exists(fileNameRaw.str()))
    {
        rawEyeImage = cv::imread(fileNameRaw.str(), CV_LOAD_IMAGE_COLOR);
    }
    else
    {
        return;
    }

    eyeProperties mEyePropertiesTemp;

    int eyeAOIXPosTemp;
    int eyeAOIYPosTemp;

    int eyeAOIWdthTemp;
    int eyeAOIHghtTemp;

    {
        std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);

        mEyePropertiesTemp.v = mEyePropertiesVariables;
        mEyePropertiesTemp.p = mEyePropertiesParameters;

        eyeAOIXPosTemp = Parameters::eyeAOIXPos;
        eyeAOIYPosTemp = Parameters::eyeAOIYPos;

        eyeAOIWdthTemp = Parameters::eyeAOIWdth;
        eyeAOIHghtTemp = Parameters::eyeAOIHght;
    }

    // Grab raw images

    // Pupil tracking algorithm

    cv::Rect eyeRegion(eyeAOIXPosTemp, eyeAOIYPosTemp, eyeAOIWdthTemp, eyeAOIHghtTemp);
    mEyePropertiesTemp = pupilDetector(rawEyeImage(eyeRegion), mEyePropertiesTemp);

    // Save processed images

    cv::Mat imageEye = mEyePropertiesTemp.m.image.clone();
    drawAll(imageEye, mEyePropertiesTemp);

    std::stringstream filename;
    filename << editDataDirectory.toStdString()
             << "/images/processed/" << editSubjectName.toStdString()
             << "_" << editDataIndex << "_" << editImageIndex << ".png";

    cv::imwrite(filename.str(), imageEye);

    // Get pupil coords in screen coords

    mEyePropertiesTemp.v.xPosAbs = mEyePropertiesTemp.v.xPos + eyeAOIXPosTemp;
    mEyePropertiesTemp.v.yPosAbs = mEyePropertiesTemp.v.yPos + eyeAOIYPosTemp;

    // Record pupil positions

    {
        std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);

        vEyePropertiesVariables[editImageIndex + 1] = mEyePropertiesTemp.v;
        mEyePropertiesMiscellaneous = mEyePropertiesTemp.m;
    }
}

void MainWindow::detectPupilOneFrame()
{
    reviewPupilDetectionOneFrame();
    updateReviewImages(editImageIndex);
}

void MainWindow::detectPupilAllFrames()
{
    if (!DETECTION_IN_PROGRESS)
    {
        DETECTION_IN_PROGRESS = true;

        std::thread pupilDetectionThread(&MainWindow::reviewPupilDetectionAllFrames, this);
        pupilDetectionThread.detach();

        while (DETECTION_IN_PROGRESS)
        {
            ReviewImageSlider->setValue(editImageIndex);

            std::stringstream ss;
            ss << "<b>" << editImageIndex + 1 << " / " << editImageTotal << "</b>";
            QString title = QString::fromStdString(ss.str());
            ReviewImageFrameTextBox->setText(title);

            msWait(1000 / guiUpdateFrequency);
        }
    }
    else
    {
        DETECTION_IN_PROGRESS = false;
    }
}

void MainWindow::reviewPupilDetectionAllFrames()
{
    int initialIndex = editImageIndex; // needed for progressbar

    for (int i = initialIndex; i < editImageTotal && DETECTION_IN_PROGRESS; i++)
    {
        editImageIndex = i;
        mEyePropertiesVariables = vEyePropertiesVariables[editImageIndex];
        reviewPupilDetectionOneFrame();
    }

    if (DETECTION_IN_PROGRESS)
    {
        ReviewImageSlider->setValue(editImageIndex);

        reviewSaveExperimentData();

        DETECTION_IN_PROGRESS = false;
    }
}

void MainWindow::reviewPupilDetectionAllSessions()
{
    int i = editDataIndex;
    int stringSize = editDataDirectory.size();

    while(true)
    {
        QString editDataDirectoryNew = editDataDirectory;
        editDataDirectoryNew.replace(stringSize - 1, 1, QString::number(i));

        if (boost::filesystem::exists(editDataDirectoryNew.toStdString()))
        {
            editDataIndex = i;
            editDataDirectory = editDataDirectoryNew;

            std::thread pupilDetectionThread(&MainWindow::reviewPupilDetectionAllFrames, this);
            pupilDetectionThread.detach();

            // add lock

            i++;
        }
        else
        {
            break;
        }
    }
}


void MainWindow::onSavePupilData()
{
    reviewSaveExperimentData();
}

void MainWindow::reviewSaveExperimentData()
{
    std::vector<double> timeStamps;

    {
        // Grab time stamps

        std::stringstream filename;
        filename << editDataDirectory.toStdString()
                 << "/" << editSubjectName.toStdString()
                 << "_tracking_data_timestamps_" << editDataIndex << ".dat";

        std::ifstream data(filename.str().c_str(), std::ios::in);

        if (data.is_open())
        {
            std::string line;
            while (std::getline(data, line))
            {
                std::istringstream iss(line);
                double time;

                if (!(iss >> time))
                {
                    break; // error
                }

                timeStamps.push_back(time);
            }
        }
    }

    // save data

    std::stringstream filename;
    filename << editDataDirectory.toStdString()
             << "/" << editSubjectName.toStdString()
             << "_tracking_data_raw_" << editDataIndex << ".dat";

    std::ofstream file;
    file.open(filename.str(), std::ios::out | std::ios::ate);

    for (int i = 0; i < editImageTotal; i++)
    {
        file << vEyePropertiesVariables[i + 1].xPosAbs
                << " "
                << vEyePropertiesVariables[i + 1].yPosAbs
                << " "
                << vEyePropertiesVariables[i + 1].pupilDetected
                << " "
                << timeStamps[i]
                   << "\n";
    }

    file.close();
}




