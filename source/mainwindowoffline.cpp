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

void MainWindow::startOfflineSession()
{
    imageIndexOffline = 0;

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

    EyeTrackingOfflineWidget->setVisible(true);

    CamQImage->clearImage();
    EyeQImage->clearImage();

    // hide tab with camera parameters
    EyeTrackingParameterTabWidget->removeTab(0);

    updateOfflineImages(0);
}

void MainWindow::countNumTrials()
{
    trialTotalOffline = 0;

    while (true)
    {
        std::stringstream folderName;
        folderName << dataDirectoryOffline.toStdString()
                   << "/trial_"
                   << trialTotalOffline;
        if (!boost::filesystem::exists(folderName.str())) { break; }
        trialTotalOffline++;
    }

    OfflineTrialSlider ->setMaximum(trialTotalOffline);
    OfflineTrialSpinBox->setMaximum(trialTotalOffline);
}

void MainWindow::countNumImages()
{
    imageTotalOffline = 0;

    while (true)
    {
        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/trial_"
                 << trialIndexOffline
                 << "/raw/"
                 << imageTotalOffline
                 << ".png";

        if (!boost::filesystem::exists(filename.str())) { break; }
        else { imageTotalOffline++; }
    }
}


void MainWindow::loadOfflineSession()
{
    dataDirectoryOffline = QFileDialog::getExistingDirectory(this, tr("Select data directory"), "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    trialIndexOffline = 0; // start with first trial

    countNumTrials();

    if (trialTotalOffline > 0)
    {
        updateOfflineSession();
        updateOfflineImages(0);

        std::stringstream directoryName;
        directoryName << dataDirectoryOffline.toStdString()
                      << "/trial_"
                      << trialIndexOffline
                      << "/processed";

        if (!boost::filesystem::exists(directoryName.str()))
        {
            boost::filesystem::create_directory(directoryName.str().c_str());
        }

        if (timeMatrix.empty()) // Grab time stamps
        {
            std::stringstream directory;
            directory << dataDirectoryOffline.toStdString()
                      << "/timestamps.dat";

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

void MainWindow::changeOfflineSession(int index)
{
    trialIndexOffline = index;
    updateOfflineSession();
    updateOfflineImages(0);
}

void MainWindow::updateOfflineSession()
{
    countNumImages();

    if (imageTotalOffline > 0)
    {
        vEyePropertiesMiscellaneous.resize(imageTotalOffline);
        vEyePropertiesVariables.resize(imageTotalOffline + 1);
        vEyePropertiesVariables[0] = mEyePropertiesVariables;

        imageIndexOffline = 0;
        OfflineImageSlider->setValue(0);
        OfflineImageSlider->setMaximum(imageTotalOffline - 1);

        // Create folder

        std::stringstream directoryName;
        directoryName << dataDirectoryOffline.toStdString()
                      << "/trial_"
                      << trialIndexOffline
                      << "/processed/";

        if (!boost::filesystem::exists(directoryName.str()))
        {
            boost::filesystem::create_directory(directoryName.str().c_str());
        }
    }
}

void MainWindow::updateOfflineImages(int imgIndex)
{
    if (imageTotalOffline > 0)
    {
        std::stringstream ss;
        ss << "<b>" << imgIndex + 1 << " / " << imageTotalOffline << "</b>";
        QString title = QString::fromStdString(ss.str());
        OfflineImageFrameTextBox->setText(title);

        // Raw

        std::stringstream fileNameRaw;
        fileNameRaw << dataDirectoryOffline.toStdString()
                    << "/trial_"
                    << trialIndexOffline
                    << "/raw/"
                    << imgIndex
                    << ".png";

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
        fileName << dataDirectoryOffline.toStdString()
                 << "/trial_"
                 << trialIndexOffline
                 << "/processed/"
                 << imgIndex
                 << ".png";

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
        OfflineImageFrameTextBox->setText("<b>0 / 0</b>");
    }
}

void MainWindow::updateRawImage() // for signal from qimageopencv
{
    std::stringstream fileNameRaw;
    fileNameRaw << dataDirectoryOffline.toStdString()
                << "/trial_"
                << trialIndexOffline
                << "/raw/"
                << imageIndexOffline
                << ".png";

    if (boost::filesystem::exists(fileNameRaw.str()))
    {
        cv::Mat rawEyeImage = cv::imread(fileNameRaw.str(), CV_LOAD_IMAGE_COLOR);
        CamQImage->loadImage(rawEyeImage);
        CamQImage->setEyeAOI(Parameters::eyeAOIXPos, Parameters::eyeAOIYPos, Parameters::eyeAOIWdth, Parameters::eyeAOIHght);
        CamQImage->setImage();
    } else { CamQImage->clearImage(); }
}

void MainWindow::prevOfflineImage()
{
    if (imageIndexOffline > 0)
    {
        imageIndexOffline--;
        OfflineImageSlider->setValue(imageIndexOffline);
    }
}

void MainWindow::nextOfflineImage()
{
    if (imageIndexOffline < imageTotalOffline - 1)
    {
        imageIndexOffline++;
        OfflineImageSlider->setValue(imageIndexOffline);
    }
}

void MainWindow::setOfflineImageFrame(int imageIndex)
{
    if (!PROCESSING_ALL_IMAGES)
    {
        imageIndexOffline = imageIndex;
        mEyePropertiesVariables = vEyePropertiesVariables[imageIndexOffline];
        setVariableWidgets(mEyePropertiesVariables);
        updateOfflineImages(imageIndexOffline);
    }
}

void MainWindow::offlinePupilDetectionOneFrame()
{
    // Grab raw images

    cv::Mat rawEyeImage;

    std::stringstream fileNameRaw;
    fileNameRaw << dataDirectoryOffline.toStdString()
                << "/trial_"
                << trialIndexOffline
                << "/raw/"
                << imageIndexOffline
                << ".png";

    if (boost::filesystem::exists(fileNameRaw.str())) { rawEyeImage = cv::imread(fileNameRaw.str(), CV_LOAD_IMAGE_COLOR); }
    else                                              { return; }

    // Detect pupil

    eyeProperties mEyePropertiesTemp;

    int eyeAOIXPosTemp;
    int eyeAOIYPosTemp;
    int eyeAOIWdthTemp;
    int eyeAOIHghtTemp;

    { // Mutex lock
        std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);

        mEyePropertiesTemp.v = mEyePropertiesVariables;
        mEyePropertiesTemp.p = mEyePropertiesParameters;

        eyeAOIXPosTemp = Parameters::eyeAOIXPos;
        eyeAOIYPosTemp = Parameters::eyeAOIYPos;
        eyeAOIWdthTemp = Parameters::eyeAOIWdth;
        eyeAOIHghtTemp = Parameters::eyeAOIHght;
    }

    cv::Rect eyeRegion(eyeAOIXPosTemp, eyeAOIYPosTemp, eyeAOIWdthTemp, eyeAOIHghtTemp);
    eyeProperties mEyePropertiesNew = pupilDetection(rawEyeImage(eyeRegion), mEyePropertiesTemp);

    // Save processed images

    cv::Mat imageEye = mEyePropertiesNew.m.image.clone();
    drawAll(imageEye, mEyePropertiesNew);

    std::stringstream filename;
    filename << dataDirectoryOffline.toStdString()
             << "/trial_"
             << trialIndexOffline
             << "/processed/"
             << imageIndexOffline
             << ".png";

    cv::imwrite(filename.str(), imageEye);

    // Get pupil coords in screen coords

    mEyePropertiesNew.v.xPosAbsolute = mEyePropertiesNew.v.xPosExact + eyeAOIXPosTemp;
    mEyePropertiesNew.v.yPosAbsolute = mEyePropertiesNew.v.yPosExact + eyeAOIYPosTemp;

    // Record pupil positions

    { // Mutex lock
        std::lock_guard<std::mutex> primaryMutexLock(Parameters::primaryMutex);
        vEyePropertiesVariables[imageIndexOffline + 1] = mEyePropertiesNew.v;
        vEyePropertiesMiscellaneous[imageIndexOffline] = mEyePropertiesNew.m;
        mEyePropertiesMiscellaneous = mEyePropertiesNew.m;
    }
}

void MainWindow::detectPupilOneFrame()
{
    offlinePupilDetectionOneFrame();
    updateOfflineImages(imageIndexOffline);
}

void MainWindow::detectPupilAllFrames()
{
    if (!PROCESSING_ALL_IMAGES)
    {
        PROCESSING_ALL_IMAGES = true;

        std::thread pupilDetectionThread(&MainWindow::offlinePupilDetectionAllFrames, this);
        pupilDetectionThread.detach();

        while (PROCESSING_ALL_IMAGES)
        {
            OfflineImageSlider->setValue(imageIndexOffline);

            std::stringstream ss;
            ss << "<b>" << imageIndexOffline + 1 << " / " << imageTotalOffline << "</b>";
            QString title = QString::fromStdString(ss.str());
            OfflineImageFrameTextBox->setText(title);

            msWait(1000 / guiUpdateFrequency);
        }
    }
    else
    {
        PROCESSING_ALL_IMAGES = false;
    }
}

void MainWindow::offlinePupilDetectionAllFrames()
{
    int initialIndex = imageIndexOffline; // needed for progressbar

    for (imageIndexOffline = initialIndex; imageIndexOffline < imageTotalOffline && PROCESSING_ALL_IMAGES; imageIndexOffline++)
    {
        mEyePropertiesVariables = vEyePropertiesVariables[imageIndexOffline];
        offlinePupilDetectionOneFrame();
    }

    if (PROCESSING_ALL_IMAGES)
    {
        OfflineImageSlider->setValue(imageIndexOffline);
        offlineSaveExperimentData();
        PROCESSING_ALL_IMAGES = false;
    }
}

void MainWindow::detectPupilAllTrials()
{
    if (!PROCESSING_ALL_TRIALS)
    {
        PROCESSING_ALL_TRIALS = true;

        for (int i = trialIndexOffline; i < trialTotalOffline && PROCESSING_ALL_TRIALS; i++)
        {
            OfflineTrialSlider->setValue(i + 1);
            detectPupilAllFrames();
        }
    }

    PROCESSING_ALL_IMAGES = false;
    PROCESSING_ALL_TRIALS = false;
}

void MainWindow::offlineSaveExperimentData()
{
    { // save pupil data

        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/tracking_data.dat";

        std::ofstream file;
        file.open(filename.str());

        file << std::setw(3) << std::setfill('0') << timeMatrix[trialIndexOffline][0] << ";"; // print with leading zeros
        file << (int) timeMatrix[trialIndexOffline][1] << ";"; // time of day in milliseconds

        file << std::fixed;
        file << std::setprecision(3);

        std::string delimiter = ";";

        // write data

        for (int i = 0; i < imageTotalOffline; i++) { file << timeMatrix[trialIndexOffline][i + 2]             << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].xPosAbsolute  << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].yPosAbsolute  << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].xPosExact     << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].yPosExact     << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].pupilDetected << delimiter; }

        // additional data

        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].pupilCircumferenceExact << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].pupilAspectRatioExact   << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i + 1].pupilDetected           << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i].edgeCurvaturePrediction     << delimiter; }
        for (int i = 0; i < imageTotalOffline; i++) { file << vEyePropertiesVariables[i].edgeIntensityPrediction     << delimiter; }

        file.close();
    }

    { // save edge data

        std::stringstream filename;
        filename << dataDirectoryOffline.toStdString()
                 << "/edge_data.dat";

        std::ofstream file;
        file.open(filename.str());

        for (int i = 0; i < imageTotalOffline; i++)
        {
            int numEdges = vEyePropertiesMiscellaneous[i].edgePropertiesAll.size();

            std::string delimiter = ";";

            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].flag         << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].curvatureMax << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].curvatureMin << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].curvatureAvg << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].length       << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].size         << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].distance     << delimiter; }
            for (int j = 0; j < numEdges; j++) { file << vEyePropertiesMiscellaneous[i].edgePropertiesAll[j].intensity    << delimiter; }
            file << "\n";
        }

        file.close();
    }

    { // save pupil image

        std::stringstream directoryName;
        directoryName << dataDirectoryOffline.toStdString() << "/pupil";

        if (!boost::filesystem::exists(directoryName.str()))
        {
            boost::filesystem::create_directory(directoryName.str().c_str());
        }

        for (int i = 0; i < imageTotalOffline; i++)
        {
            if (vEyePropertiesVariables[i + 1].pupilDetected)
            {
                std::stringstream filename;
                filename << dataDirectoryOffline.toStdString()
                         << "/trial_"
                         << trialIndexOffline
                         << "/pupil/"
                         << i
                         << ".png";

                std::vector<int> compression_params;
                compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
                compression_params.push_back(0);

                cv::imwrite(filename.str(), vEyePropertiesMiscellaneous[i].imagePupil, compression_params);
            }
        }
    }
}

void MainWindow::offlineCombineExperimentData()
{
    for (int iTrial = 0; iTrial < trialTotalOffline; iTrial++)
    {
        OfflineTrialSlider->setValue(iTrial + 1);

        std::stringstream fileNameRead;
        fileNameRead << dataDirectoryOffline.toStdString()
                     << "/tracking_data.dat";

        std::ifstream dataFile;
        dataFile.open(fileNameRead.str().c_str());

        std::string line;
        while (dataFile)
        {
            if (!std::getline(dataFile, line)) { break; }
        }

        dataFilename = (DataFilenameLineEdit->text()).toStdString();

        std::stringstream fileNameWrite;
        fileNameWrite << dataDirectory
                      << "/" << dataFilename;

        std::ofstream file;
        if (!boost::filesystem::exists(fileNameWrite.str()))
        {
            file.open(fileNameWrite.str(), std::ios::out | std::ios::ate);
        }
        else
        {
            file.open(fileNameWrite.str(), std::ios_base::app);
            file << "\n";
        }

        file << line;

        file.close();
    }
}

void MainWindow::setPupilPosition(double xPos, double yPos)
{
    std::lock_guard<std::mutex> secondaryMutexLock(Parameters::secondaryMutex);

    if (xPos > 0 && xPos < Parameters::eyeAOIWdth && yPos > 0 && yPos < Parameters::eyeAOIHght)
    {
        int pupilHaarWdth = round(mEyePropertiesVariables.pupilCircumferencePrediction / M_PI);
        int pupilHaarWdthOffset = pupilHaarWdth + round(pupilHaarWdth * mEyePropertiesParameters.pupilOffset * 2);

        mEyePropertiesVariables.searchRadius  = ceil(0.5 * pupilHaarWdthOffset);
        mEyePropertiesVariables.xPosPredicted = xPos;
        mEyePropertiesVariables.yPosPredicted = yPos;

        if (imageIndexOffline == 0)
        {
            mEyePropertiesVariables.xPosExact = xPos;
            mEyePropertiesVariables.yPosExact = yPos;
        }
    }
}
