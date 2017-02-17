//  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

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

#ifndef UEYEOPENCV_H
#define UEYEOPENCV_H

// Files

#include "constants.h"
#include "parameters.h"
#include "structures.h"

// Libraries

// Standard Template

#include <chrono>
#include <iostream>
#include <stdio.h>
#include <thread>

// OpenCV

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

// UEye

#include <ueye.h>

// LibUSB

#ifdef __linux__
    #include <libusb-1.0/libusb.h>
#endif

class UEyeOpencvCam
{

public:

    bool allocateMemory(int wdth, int hght);
    bool findCamera();
    bool freeImageMemory();
    bool setAOI(int xAOI, int yAOI, int wAOI, int hAOI);
    bool setColorMode();
    bool setSubSampling(int);
    bool startVideoCapture();
    double getExposure();
    double setFrameRate(double FPS);
    imageInfo getFrame(); // returns cv::Mat of camera frame
    int getHardwareGain();
    int initCamera();
    int setBlackLevelOffset(int nOffset);
    std::vector<double> getExposureRange();
    std::vector<double> getFrameRateRange();
    std::vector<int>    getBlackLevelOffsetRange();
    std::vector<int>    getPixelClockRange();
    UEyeOpencvCam(); // default constructor
    void exitCamera();
    void setAutoGain(bool FLAG);
    void setBlackLevelMode(bool FLAG);
    void setDeviceInfo(int, int);
    void setExposure(double pExp);
    void setGainBoost(bool FLAG);
    void setHardwareGain(int nMaster);
    void setPixelClock(int pixelClock);
    void startRecording();
    void stopRecording();
    void threadFrameCapture();
    ~UEyeOpencvCam(); // deconstructor

private:

    bool DEVICE_INITIALIZED;
    bool EVENT_ENABLED;
    bool TRIAL_RECORDING;
    bool THREAD_ACTIVE;
    char* ppcImgMem;
    HIDS hCam;
    int frameCount;
    int frameIndex;
    int height;
    int numberOfImageBuffers;
    int pid;
    int idVendor;
    int idProduct;
    int width;
    std::condition_variable exitCV;
    std::condition_variable fullBufferCV;
    std::mutex exitMutex;
    std::vector<imageInfo> vImageInfo;

};

#endif // UEYEOPENCV
