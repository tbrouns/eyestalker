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

#include "mainwindow.h"

void MainWindow::loadSettings(QString filename)
{
    QSettings settings(filename, QSettings::IniFormat);

    cameraAOIFractionHght           = settings.value("CamAOIHghtFraction",          cameraAOIFractionHghtDefaultLeft).toDouble();
    cameraAOIFractionWdth           = settings.value("CamAOIWdthFraction",          cameraAOIFractionWdthDefaultLeft).toDouble();
    cameraAOIFractionXPos           = settings.value("CamAOIXPosFraction",          cameraAOIFractionXPosDefaultLeft).toDouble();
    cameraAOIFractionYPos           = settings.value("CamAOIYPosFraction",          cameraAOIFractionYPosDefaultLeft).toDouble();
    cameraFrameRateDesired          = settings.value("CameraFrameRateDesired",      250).toInt();
    cameraSubSamplingFactor         = settings.value("SubSamplingFactor",           1).toInt();
    dataDirectory                   = settings.value("DataDirectory",               "").toString().toStdString();
    dataDirectoryOffline            = settings.value("DataDirectoryOffline",        "").toString();
    dataFilename                    = settings.value("DataFilename",                "experiment_data").toString().toStdString();
    trialIndexOffline               = settings.value("trialIndexOffline",           0).toInt();
    imageTotalOffline               = settings.value("imageTotalOffline",           0).toInt();
    subjectIdentifier               = settings.value("SubjectName",                 "").toString();
    eyeAOIHghtFraction              = settings.value("AOIHghtFraction",             1.0).toDouble();
    eyeAOIWdthFraction              = settings.value("AOIWdthFraction",             1.0).toDouble();
    beadAOIHghtFraction             = settings.value("AOIBeadHghtFraction",         0.6).toDouble();
    beadAOIWdthFraction             = settings.value("AOIBeadWdthFraction",         0.3).toDouble();
    flashThreshold                  = settings.value("FlashThreshold",              230).toInt();
    GAIN_AUTO                       = settings.value("GainAuto",                    true).toBool();
    GAIN_BOOST                      = settings.value("GainBoost",                   false).toBool();
    Parameters::eyeAOIXPosFraction  = settings.value("AOIXPosRelative",             0.0).toDouble();
    Parameters::eyeAOIYPosFraction  = settings.value("AOIYPosRelative",             0.0).toDouble();
    Parameters::beadAOIXPosFraction = settings.value("AOIBeadXPosRelative",         0.2).toDouble();
    Parameters::beadAOIYPosFraction = settings.value("AOIBeadYPosRelative",         0.5).toDouble();
    flashAOI.hght                   = settings.value("FlashAOIHght",                100).toInt();
    flashAOI.wdth                   = settings.value("FlashAOIWdth",                60).toInt();
    flashAOI.xPos                   = settings.value("FlashAOIXPos",                227).toInt();
    flashAOI.yPos                   = settings.value("FlashAOIYPos",                500).toInt();
    SAVE_ASPECT_RATIO               = settings.value("SaveAspectRatio",             true).toBool();
    SAVE_CIRCUMFERENCE              = settings.value("SaveCircumference",           true).toBool();
    SAVE_POSITION                   = settings.value("SavePosition",                true).toBool();
    SAVE_EYE_IMAGE                  = settings.value("SaveEyeImage",                true).toBool();
    subjectIdentifier               = settings.value("SubjectIdentifier",           "").toString();
    trialTimeLength                 = settings.value("TrialTimeLength",             1500).toInt();

    cameraAOIWdthMax = Parameters::cameraXResolution / (double) cameraSubSamplingFactor; // maximum possible AOI size
    cameraAOIHghtMax = Parameters::cameraYResolution / (double) cameraSubSamplingFactor;

    updateCamAOIx();
    updateCamAOIy();

    detectionParameters mDetectionParametersEye  = loadParameters(filename, "Eye",  parametersEye);
    mParameterWidgetEye ->setStructure(mDetectionParametersEye);
    mVariableWidgetEye->resetStructure(mDetectionParametersEye, Parameters::eyeAOI);

    detectionParameters mDetectionParametersBead = loadParameters(filename, "Bead", parametersBead);
    mParameterWidgetBead ->setStructure(mDetectionParametersBead);
    mVariableWidgetBead->resetStructure(mDetectionParametersBead, Parameters::beadAOI);
}

detectionParameters MainWindow::loadParameters(QString filename, QString prefix, std::vector<double> parameters)
{
    QSettings settings(filename, QSettings::IniFormat);

    detectionParameters mDetectionParameters;
    mDetectionParameters.alphaAverages                       = settings.value(prefix + "AlphaAverage",                    parameters[0]).toDouble();
    mDetectionParameters.alphaFeatures                      = settings.value(prefix + "AlphaFeatures",                   parameters[1]).toDouble();
    mDetectionParameters.alphaCertainty                     = settings.value(prefix + "AlphaCertainty",                  parameters[2]).toDouble();
    mDetectionParameters.alphaPosition                      = settings.value(prefix + "AlphaPosition",                   parameters[3]).toDouble();
    mDetectionParameters.cannyBlurLevel                     = settings.value(prefix + "CannyBlurLevel",                  parameters[4]).toInt();
    mDetectionParameters.cannyKernelSize                    = settings.value(prefix + "CannyKernelSize",                 parameters[5]).toInt();
    mDetectionParameters.cannyThresholdLow                  = settings.value(prefix + "CannyThresholdLow",               parameters[6]).toDouble();
    mDetectionParameters.cannyThresholdHigh                 = settings.value(prefix + "CannyThresholdHigh",              parameters[7]).toDouble();
    mDetectionParameters.curvatureOffset                    = settings.value(prefix + "CurvatureOffset",                 parameters[8]).toDouble();
    mDetectionParameters.edgeLengthFraction                 = settings.value(prefix + "EdgeLengthFraction",              parameters[9]).toDouble();
    mDetectionParameters.ellipseFitNumberMaximum            = settings.value(prefix + "EllipseFitNumberMaximum",         parameters[10]).toInt();
    mDetectionParameters.ellipseFitErrorMaximum             = settings.value(prefix + "EllipseFitErrorMaximum",          parameters[11]).toDouble();
    mDetectionParameters.glintWdth                          = settings.value(prefix + "GlintSize",                       parameters[12]).toInt();
    mDetectionParameters.circumferenceMax                   = settings.value(prefix + "CircumferenceMax",                parameters[13]).toDouble();
    mDetectionParameters.circumferenceMin                   = settings.value(prefix + "CircumferenceMin",                parameters[14]).toDouble();
    mDetectionParameters.aspectRatioMin                     = settings.value(prefix + "AspectRatioMin",                  parameters[15]).toDouble();
    mDetectionParameters.circumferenceOffset                = settings.value(prefix + "CircumferenceOffset",             parameters[16]).toDouble();
    mDetectionParameters.changeThresholdCircumference       = settings.value(prefix + "CircumferenceChangeThreshold",    parameters[17]).toDouble();
    mDetectionParameters.changeThresholdAspectRatio         = settings.value(prefix + "AspectRatioChangeThreshold",      parameters[18]).toDouble();
    mDetectionParameters.changeThresholdPosition            = settings.value(prefix + "DisplacementChangeThreshold",     parameters[19]).toDouble();
    mDetectionParameters.scoreThreshold                     = settings.value(prefix + "ScoreThreshold",                  parameters[20]).toDouble();

    return mDetectionParameters;
}

void MainWindow::saveSettings(QString filename)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue("AOIHghtFraction",                eyeAOIHghtFraction);
    settings.setValue("AOIWdthFraction",                eyeAOIWdthFraction);
    settings.setValue("AOIXPosRelative",                Parameters::eyeAOIXPosFraction);
    settings.setValue("AOIYPosRelative",                Parameters::eyeAOIYPosFraction);
    settings.setValue("AOIBeadHghtFraction",            beadAOIHghtFraction);
    settings.setValue("AOIBeadWdthFraction",            beadAOIWdthFraction);
    settings.setValue("AOIBeadXPosRelative",            Parameters::beadAOIXPosFraction);
    settings.setValue("AOIBeadYPosRelative",            Parameters::beadAOIYPosFraction);
    settings.setValue("DataDirectory",                  QString::fromStdString(dataDirectory));
    settings.setValue("DataDirectoryOffline",           dataDirectoryOffline);
    settings.setValue("GainAuto",                       GAIN_AUTO);
    settings.setValue("GainBoost",                      GAIN_BOOST);
    settings.setValue("CamAOIHghtFraction",             cameraAOIFractionHght);
    settings.setValue("CamAOIWdthFraction",             cameraAOIFractionWdth);
    settings.setValue("CamAOIXPosFraction",             cameraAOIFractionXPos);
    settings.setValue("CamAOIYPosFraction",             cameraAOIFractionYPos);
    settings.setValue("CameraFrameRateDesired",         cameraFrameRateDesired);
    settings.setValue("DataFilename",                   QString::fromStdString(dataFilename));
    settings.setValue("FlashAOIHght",                   flashAOI.hght);
    settings.setValue("FlashAOIWdth",                   flashAOI.wdth);
    settings.setValue("FlashAOIXPos",                   flashAOI.xPos);
    settings.setValue("FlashAOIYPos",                   flashAOI.yPos);
    settings.setValue("FlashThreshold",                 flashThreshold);
    settings.setValue("SaveAspectRatio",                SAVE_ASPECT_RATIO);
    settings.setValue("SaveCircumference",              SAVE_CIRCUMFERENCE);
    settings.setValue("SavePosition",                   SAVE_POSITION);
    settings.setValue("SaveEyeImage",                   SAVE_EYE_IMAGE);
    settings.setValue("SubjectName",                    subjectIdentifier);
    settings.setValue("SubSamplingFactor",              cameraSubSamplingFactor);
    settings.setValue("TrialTimeLength",                TrialTimeLengthLineEdit->text().toInt());

    detectionParameters mDetectionParametersEye  = mParameterWidgetEye->getStructure();
    saveParameters(filename,  "Eye", mDetectionParametersEye);

    detectionParameters mDetectionParametersBead = mParameterWidgetBead->getStructure();
    saveParameters(filename, "Bead", mDetectionParametersBead);
}

void MainWindow::saveParameters(QString filename, QString prefix, detectionParameters mDetectionParameters)
{
    QSettings settings(filename, QSettings::IniFormat);

    settings.setValue(prefix + "AlphaAverage",                   mDetectionParameters.alphaAverages);
    settings.setValue(prefix + "AlphaFeatures",                  mDetectionParameters.alphaFeatures);
    settings.setValue(prefix + "AlphaPosition",                  mDetectionParameters.alphaPosition);
    settings.setValue(prefix + "AlphaCertainty",                 mDetectionParameters.alphaCertainty);
    settings.setValue(prefix + "CannyBlurLevel",                 mDetectionParameters.cannyBlurLevel);
    settings.setValue(prefix + "CannyKernelSize",                mDetectionParameters.cannyKernelSize);
    settings.setValue(prefix + "CannyThresholdLow",              mDetectionParameters.cannyThresholdLow);
    settings.setValue(prefix + "CannyThresholdHigh",             mDetectionParameters.cannyThresholdHigh);
    settings.setValue(prefix + "CircumferenceOffset",            mDetectionParameters.circumferenceOffset);
    settings.setValue(prefix + "CircumferenceMax",               mDetectionParameters.circumferenceMax);
    settings.setValue(prefix + "CircumferenceMin",               mDetectionParameters.circumferenceMin);
    settings.setValue(prefix + "CurvatureOffset",                mDetectionParameters.curvatureOffset);
    settings.setValue(prefix + "EllipseFitNumberMaximum",        mDetectionParameters.ellipseFitNumberMaximum);
    settings.setValue(prefix + "EllipseFitErrorMaximum",         mDetectionParameters.ellipseFitErrorMaximum);
    settings.setValue(prefix + "AspectRatioMin",                 mDetectionParameters.aspectRatioMin);
    settings.setValue(prefix + "GlintSize",                      mDetectionParameters.glintWdth);
    settings.setValue(prefix + "CircumferenceChangeThreshold",   mDetectionParameters.changeThresholdCircumference);
    settings.setValue(prefix + "AspectRatioChangeThreshold",     mDetectionParameters.changeThresholdAspectRatio);
    settings.setValue(prefix + "DisplacementChangeThreshold",    mDetectionParameters.changeThresholdPosition);
    settings.setValue(prefix + "ScoreThreshold",                 mDetectionParameters.scoreThreshold);
}

void MainWindow::onResetParameters()
{
    QString text = "Do you wish to reset all parameters to their default values?";
    ConfirmationWindow mConfirmationWindow(text);
    mConfirmationWindow.setWindowTitle("Please select option");

    if(mConfirmationWindow.exec() == QDialog::Rejected) { return; }

    QString filename = "";
    loadSettings(filename);
    mParameterWidgetEye->reset();
    mVariableWidgetEye->resetStructure(mParameterWidgetEye->getStructure(), Parameters::eyeAOI);
}

void MainWindow::setBeadDetection(int state)
{
    if (!TRIAL_RECORDING && !PROCESSING_ALL_TRIALS && !PROCESSING_ALL_IMAGES)
    {
        int index = 2 + (int) Parameters::ONLINE_PROCESSING;

        if (state)
        {
            MainTabWidget->setUpdatesEnabled(false);
            MainTabWidget->insertTab(index, BeadTrackingScrollArea, tr("Bead-tracking"));
            MainTabWidget->setUpdatesEnabled(true);
        }
        else
        {
            MainTabWidget->removeTab(index);
        }

        mParameterWidgetBead->setState(state);

        CamQImage->showAOIBead(state);
        CamQImage->setImage();
    }
    else
    {
        BeadDetectionCheckBox->setChecked(false);
    }
}

void MainWindow::onSetRealTime(int state)
{
    if (Parameters::ONLINE_PROCESSING && !TRIAL_RECORDING)
    {
        if (state) { SAVE_EYE_IMAGE = false; }
        else       { SAVE_EYE_IMAGE = true;  }
    }
    else
    {
        RealTimeEyeTrackingCheckBox->setChecked(false);
    }
}




void MainWindow::onResetFlashIntensity()
{
    flashMinIntensity = 0;
    FlashThresholdSlider->setMinimum(0);
    FlashThresholdSlider->setValue(0);
}

void MainWindow::setCameraPixelClock(int value)
{
    cameraPixelClock = value;
    mUEyeOpencvCam.setPixelClock(cameraPixelClock);
    CameraPixelClockLabel->setText(QString::number(cameraPixelClock));

    // Set new frame rate

    std::vector<double> frameRateRange = mUEyeOpencvCam.getFrameRateRange();
    CameraFrameRateSlider->setDoubleRange(frameRateRange[0], frameRateRange[1]);
    CameraFrameRateSlider->setDoubleValue(frameRateRange[1]);
    setCameraFrameRate(frameRateRange[1]); // set frame-rate to maximum
}

void MainWindow::setCameraFrameRate(double value)
{
    cameraFrameRate = mUEyeOpencvCam.setFrameRate(value);
    CameraFrameRateLabel->setText(QString::number(cameraFrameRate, 'f', 1));

    // Set new exposure

    std::vector<double> exposureRange = mUEyeOpencvCam.getExposureRange();
    CameraExposureSlider->setDoubleRange(exposureRange[0], exposureRange[1]);
    CameraExposureSlider->setDoubleValue(exposureRange[1]);
    setCameraExposure(exposureRange[1]);
}

void MainWindow::setCameraExposure(double value)
{
    mUEyeOpencvCam.setExposure(value);
    CameraExposureLabel->setText(QString::number(value, 'f', 2));
}

void MainWindow::setCameraBlackLevelOffset(int value)
{
    double blackLevelOffset = mUEyeOpencvCam.setBlackLevelOffset(value);
    CameraBlackLevelOffsetLabel->setText(QString::number(blackLevelOffset));
}

void MainWindow::setCameraBlackLevelMode(int state)
{
    if (state) { mUEyeOpencvCam.setBlackLevelMode(true);  }
    else       { mUEyeOpencvCam.setBlackLevelMode(false); }
}

void MainWindow::setCameraGainBoost(int state)
{
    if (state) { mUEyeOpencvCam.setGainBoost(true);  }
    else       { mUEyeOpencvCam.setGainBoost(false); }
}

void MainWindow::setCameraAutoGain(int state)
{
    GAIN_AUTO = true;
    if (state) { mUEyeOpencvCam.setAutoGain(GAIN_AUTO); }
    else       { mUEyeOpencvCam.setAutoGain(GAIN_AUTO); }
}

void MainWindow::setCameraHardwareGain(int val)
{
    if (!CameraHardwareGainAutoCheckBox->checkState())
    {
        mUEyeOpencvCam.setHardwareGain(val);
    }

    CameraHardwareGainLabel->setText(QString::number(val));
}

void MainWindow::setCameraSubSampling(int state)
{
    double subSamplingChange = 0;

    if (state == 1)
    {
        subSamplingChange = 0.5 * cameraSubSamplingFactor;
        cameraSubSamplingFactor = 2;
    }
    else if (state == 0)
    {
        subSamplingChange = cameraSubSamplingFactor;
        cameraSubSamplingFactor = 1;
    }

    cameraAOIWdthMax = cameraAOIWdthMax * subSamplingChange;
    cameraAOIHghtMax = cameraAOIHghtMax * subSamplingChange;

    updateCamAOIx();
    updateCamAOIy();

    if (Parameters::CAMERA_RUNNING)
    {
        if (mUEyeOpencvCam.freeImageMemory())
        {
            if (mUEyeOpencvCam.setSubSampling(cameraSubSamplingFactor))
            {
                if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght);
                    Parameters::CAMERA_READY = true;
                }
            }
        }

        getCameraParameters();
    }
}

void MainWindow::onCropAOI()
{
    int absXPos = Parameters::eyeAOI.xPos + Parameters::cameraAOI.xPos;
    int absYPos = Parameters::eyeAOI.yPos + Parameters::cameraAOI.yPos;

    double fracXPos = absXPos / (double) cameraAOIWdthMax;
    double fracYPos = absYPos / (double) cameraAOIHghtMax;
    double fracWdth = (Parameters::eyeAOI.wdth - cameraAOIWdthMin) / (double) (cameraAOIWdthMax - cameraAOIWdthMin);
    double fracHght = (Parameters::eyeAOI.hght - cameraAOIHghtMin) / (double) (cameraAOIHghtMax - cameraAOIHghtMin);

    Parameters::eyeAOIXPosFraction = 1.0;
    Parameters::eyeAOIYPosFraction = 1.0;
    eyeAOIWdthFraction = 1.0;
    eyeAOIHghtFraction = 1.0;

    CamEyeAOIXPosSlider->setDoubleValue(fracXPos);
    CamEyeAOIYPosSlider->setDoubleValue(fracYPos);
    CamEyeAOIWdthSlider->setDoubleValue(fracWdth);
    CamEyeAOIHghtSlider->setDoubleValue(fracHght);
}

void MainWindow::updateCamAOIx()
{
    Parameters::cameraAOI.wdth = floor(((cameraAOIWdthMax - cameraAOIWdthMin) * cameraAOIFractionWdth + cameraAOIWdthMin) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    Parameters::cameraAOI.xPos = floor((cameraAOIWdthMax * cameraAOIFractionXPos) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
    CamEyeAOIXPosSlider->setDoubleMaximum((cameraAOIWdthMax - Parameters::cameraAOI.wdth) / (double) cameraAOIWdthMax);

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        updateEyeAOIx(); }
}

void MainWindow::updateCamAOIy()
{
    Parameters::cameraAOI.hght = floor(((cameraAOIHghtMax - cameraAOIHghtMin) * cameraAOIFractionHght + cameraAOIHghtMin) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    Parameters::cameraAOI.yPos = floor((cameraAOIHghtMax * cameraAOIFractionYPos) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
    CamEyeAOIYPosSlider->setDoubleMaximum((cameraAOIHghtMax - Parameters::cameraAOI.hght) / (double) cameraAOIHghtMax);

    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        updateEyeAOIy(); }
}

void MainWindow::updateEyeAOIx()
{
    Parameters::eyeAOI.wdth = round(Parameters::cameraAOI.wdth * eyeAOIWdthFraction);
    Parameters::eyeAOI.xPos = round(Parameters::cameraAOI.wdth * Parameters::eyeAOIXPosFraction);

    if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::cameraAOI.wdth)
    {   Parameters::eyeAOI.xPos = Parameters::cameraAOI.wdth - Parameters::eyeAOI.wdth; }

    Parameters::beadAOI.wdth = round(Parameters::cameraAOI.wdth * beadAOIWdthFraction);
    Parameters::beadAOI.xPos = round(Parameters::cameraAOI.wdth * Parameters::beadAOIXPosFraction);

    if (Parameters::beadAOI.xPos + Parameters::beadAOI.wdth > Parameters::cameraAOI.wdth)
    {   Parameters::beadAOI.xPos = Parameters::cameraAOI.wdth - Parameters::beadAOI.wdth; }
}

void MainWindow::updateEyeAOIy()
{
    Parameters::eyeAOI.hght = round(Parameters::cameraAOI.hght * eyeAOIHghtFraction);
    Parameters::eyeAOI.yPos = round(Parameters::cameraAOI.hght * Parameters::eyeAOIYPosFraction);

    if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght > Parameters::cameraAOI.hght)
    {   Parameters::eyeAOI.yPos = Parameters::cameraAOI.hght - Parameters::eyeAOI.hght; }

    Parameters::beadAOI.hght = round(Parameters::cameraAOI.hght * beadAOIHghtFraction);
    Parameters::beadAOI.yPos = round(Parameters::cameraAOI.hght * Parameters::beadAOIYPosFraction);

    if (Parameters::beadAOI.yPos + Parameters::beadAOI.hght > Parameters::cameraAOI.hght)
    {   Parameters::beadAOI.yPos = Parameters::cameraAOI.hght - Parameters::beadAOI.hght; }
}

void MainWindow::setCamEyeAOIWdth(double fraction)
{
    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        cameraAOIFractionWdth = fraction;
        Parameters::cameraAOI.wdth = floor(((cameraAOIWdthMax - cameraAOIWdthMin) * cameraAOIFractionWdth + cameraAOIWdthMin) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;
        updateEyeAOIx(); }

    if (mUEyeOpencvCam.freeImageMemory())
    {
        if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
        {
            if (Parameters::cameraAOI.xPos + Parameters::cameraAOI.wdth <= cameraAOIWdthMax)
            {
                if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    Parameters::CAMERA_READY = true;
                    getCameraParameters();
                }
            }
        }
    }

    CamEyeAOIXPosSlider->setDoubleMaximum((cameraAOIWdthMax - Parameters::cameraAOI.wdth) / (double) cameraAOIWdthMax);
}

void MainWindow::setCamEyeAOIHght(double fraction)
{
    { std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);
        cameraAOIFractionHght = fraction;
        Parameters::cameraAOI.hght = floor(((cameraAOIHghtMax - cameraAOIHghtMin) * cameraAOIFractionHght + cameraAOIHghtMin) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;
        updateEyeAOIy(); }

    if (mUEyeOpencvCam.freeImageMemory())
    {
        if (mUEyeOpencvCam.allocateMemory(Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
        {
            if (Parameters::cameraAOI.yPos + Parameters::cameraAOI.hght <= cameraAOIHghtMax)
            {
                if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
                {
                    Parameters::CAMERA_READY = true;
                    getCameraParameters();
                }
            }
        }
    }

    CamEyeAOIYPosSlider->setDoubleMaximum((cameraAOIHghtMax - Parameters::cameraAOI.hght) / (double) cameraAOIHghtMax);
}

void MainWindow::setCamEyeAOIXPos(double fraction)
{
    cameraAOIFractionXPos = fraction;

    Parameters::cameraAOI.xPos = floor((cameraAOIWdthMax * cameraAOIFractionXPos) / (double) cameraAOIWdthStepSize) * cameraAOIWdthStepSize;

    if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
    {
        Parameters::CAMERA_READY = true;
        getCameraParameters();
    }
}

void MainWindow::setCamEyeAOIYPos(double fraction)
{
    cameraAOIFractionYPos = fraction;

    Parameters::cameraAOI.yPos = floor((cameraAOIHghtMax * cameraAOIFractionYPos) / (double) cameraAOIHghtStepSize) * cameraAOIHghtStepSize;

    if (mUEyeOpencvCam.setAOI(Parameters::cameraAOI.xPos, Parameters::cameraAOI.yPos, Parameters::cameraAOI.wdth, Parameters::cameraAOI.hght))
    {
        Parameters::CAMERA_READY = true;
        getCameraParameters();
    }
}

void MainWindow::setEyeAOIWdth(double fraction)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    eyeAOIWdthFraction = fraction;

    Parameters::eyeAOI.wdth = round(Parameters::cameraAOI.wdth * eyeAOIWdthFraction);

    if (Parameters::eyeAOI.xPos + Parameters::eyeAOI.wdth > Parameters::cameraAOI.wdth)
    {   Parameters::eyeAOI.xPos = Parameters::cameraAOI.wdth - Parameters::eyeAOI.wdth; }
}

void MainWindow::setEyeAOIHght(double fraction)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    eyeAOIHghtFraction = fraction;

    Parameters::eyeAOI.hght = round(Parameters::cameraAOI.hght * eyeAOIHghtFraction);

    if (Parameters::eyeAOI.yPos + Parameters::eyeAOI.hght > Parameters::cameraAOI.hght)
    {   Parameters::eyeAOI.yPos = Parameters::cameraAOI.hght - Parameters::eyeAOI.hght; }
}

void MainWindow::setFlashAOIXPos(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.xPos = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::setFlashAOIYPos(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.yPos = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::setFlashAOIWdth(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.wdth = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::setFlashAOIHght(int val)
{
    std::lock_guard<std::mutex> mainMutexLock(Parameters::mainMutex);

    flashAOI.hght = val;
    CamQImage->setAOIFlash(flashAOI);
}

void MainWindow::setFlashThreshold(int val)
{
    flashThreshold = val;
    FlashThresholdLabel->setText(QString::number(flashThreshold));
}

void MainWindow::onSetTrialIndex(int val)
{
    trialIndex = val;
}

void MainWindow::setSaveDataAspectRatio(int state)
{
    if (state) { SAVE_ASPECT_RATIO = true;  }
    else       { SAVE_ASPECT_RATIO = false; }
}

void MainWindow::setSaveDataCircumference(int state)
{
    if (state) { SAVE_CIRCUMFERENCE = true;  }
    else       { SAVE_CIRCUMFERENCE = false; }
}

void MainWindow::setSaveDataPosition(int state)
{
    if (state) { SAVE_POSITION = true;  }
    else       { SAVE_POSITION = false; }
}

void MainWindow::onSetAOIEyeLeft()
{
    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPosDefaultLeft);
    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPosDefaultLeft);
    CamEyeAOIWdthSlider->setDoubleValue(cameraAOIFractionWdthDefaultLeft);
    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHghtDefaultLeft);
}

void MainWindow::onSetAOIEyeRght()
{
    CamEyeAOIXPosSlider->setDoubleValue(cameraAOIFractionXPosDefaultRght);
    CamEyeAOIYPosSlider->setDoubleValue(cameraAOIFractionYPosDefaultRght);
    CamEyeAOIWdthSlider->setDoubleValue(cameraAOIFractionWdthDefaultRght);
    CamEyeAOIHghtSlider->setDoubleValue(cameraAOIFractionHghtDefaultRght);
}

void MainWindow::setDrawHaar(int state)
{
    if (state) { Parameters::drawFlags.haar = true;  }
    else       { Parameters::drawFlags.haar = false; }
}

void MainWindow::setDrawEdge(int state)
{
    if (state) { Parameters::drawFlags.edge = true;  }
    else       { Parameters::drawFlags.edge = false; }
}

void MainWindow::setDrawElps(int state)
{
    if (state) { Parameters::drawFlags.elps = true;  }
    else       { Parameters::drawFlags.elps = false; }
}
