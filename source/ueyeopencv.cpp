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

#include "ueyeopencv.h"

UEyeOpencvCam::UEyeOpencvCam()
{
    numberOfImageBuffers = 1000;
    vImageInfo.resize(numberOfImageBuffers);

    frameIndex = 0;
    hCam = 0;

    DEVICE_INITIALIZED = false;
    EVENT_ENABLED = true;
    TRIAL_RECORDING = false;
    THREAD_ACTIVE = false;
}

void UEyeOpencvCam::setDeviceInfo(int idV, int idP)
{
    idVendor = idV;
    idProduct = idP;
}

bool UEyeOpencvCam::findCamera()
{
    bool returnVal = false;

#ifdef __linux__

    THREAD_ACTIVE = true;

    libusb_context *context = NULL;
    libusb_device **list = NULL;

    libusb_init(&context); // initializing a library session
    ssize_t count = libusb_get_device_list(context, &list); // get list of all usb devices

    for (ssize_t iDevice = 0; iDevice < count; iDevice++)
    {
        libusb_device *device = list[iDevice];
        libusb_device_descriptor desc;
        libusb_get_device_descriptor(device, &desc);

        if (desc.idVendor == idVendor && desc.idProduct == idProduct)
        {
            returnVal = true;
            break;
        }
    }

    libusb_free_device_list(list, 1);
    libusb_exit(context); //close the session

    std::unique_lock<std::mutex> lck(exitMutex);
    THREAD_ACTIVE = false;
    exitCV.notify_one();

#else

    returnVal = true;

#endif

    return returnVal;
}

int UEyeOpencvCam::initCamera()
{
    if (DEVICE_INITIALIZED) { return 2; }
    else
    {
        int retInt = is_InitCamera(&hCam, 0);
        if (retInt != IS_SUCCESS) { return 0; }
    }

    DEVICE_INITIALIZED = true;
    return 1;
}

bool UEyeOpencvCam::setSubSampling(int subSamplingFactor)
{
    // Set subsampling

    if (subSamplingFactor == 2)
    {
        if (is_SetSubSampling(hCam, IS_SUBSAMPLING_2X_VERTICAL | IS_SUBSAMPLING_2X_HORIZONTAL) != IS_SUCCESS) { return false; }
    }
    else if (subSamplingFactor == 1)
    {
        if (is_SetSubSampling(hCam, IS_SUBSAMPLING_DISABLE) != IS_SUCCESS) { return false; }
    }
    else { return false; }

    return true; // SUCCESS
}

bool UEyeOpencvCam::freeImageMemory()
{
    std::lock_guard<std::mutex> lock(Parameters::frameCaptureMutex);

    if (is_FreeImageMem(hCam, ppcImgMem, pid) != IS_SUCCESS) { return false; }
    else
    {
        vImageInfo.clear();
        frameCount = 0;
        Parameters::CAMERA_READY = false;
    }

    return true;
}

bool UEyeOpencvCam::setAOI(int xAOI, int yAOI, int wAOI, int hAOI)
{
    std::lock_guard<std::mutex> lock(Parameters::frameCaptureMutex);

    // Crop to area of interest

    IS_RECT rectAOI;
    rectAOI.s32X = xAOI;
    rectAOI.s32Y = yAOI;
    rectAOI.s32Width  = wAOI;
    rectAOI.s32Height = hAOI;

    int retInt = is_AOI(hCam, IS_AOI_IMAGE_SET_AOI, (void*)&rectAOI, sizeof(rectAOI));

    if (retInt != IS_SUCCESS)
    {
        return false;
    }

    return true;
}

bool UEyeOpencvCam::allocateMemory(int wdth, int hght)
{
    std::lock_guard<std::mutex> lock(Parameters::frameCaptureMutex);

    width  = wdth;
    height = hght;

    // Allocate and set memory

    if (is_AllocImageMem(hCam, width, height, 24, &ppcImgMem, &pid) != IS_SUCCESS) { return false; }
    if (is_SetImageMem(hCam, ppcImgMem, pid) != IS_SUCCESS) { return false; }

    return true;
}

bool UEyeOpencvCam::setColorMode()
{
    int retInt = is_SetColorMode(hCam, IS_CM_BGR8_PACKED);

    if (retInt != IS_SUCCESS) { return false; }

    return true;
}

bool UEyeOpencvCam::startVideoCapture()
{
    if (is_CaptureVideo(hCam, IS_DONT_WAIT) != IS_SUCCESS)
    {
        Parameters::CAMERA_RUNNING = false;
        return false;
    }

    std::thread frameCaptureThread(&UEyeOpencvCam::threadFrameCapture, this);
    frameCaptureThread.detach();

    return true;
}

void UEyeOpencvCam::exitCamera()
{
    is_StopLiveVideo(hCam, IS_FORCE_VIDEO_STOP);

    DEVICE_INITIALIZED = false;
    EVENT_ENABLED      = false;

    std::unique_lock<std::mutex> lck(exitMutex);
    while (THREAD_ACTIVE) exitCV.wait(lck);

    is_ExitCamera(hCam); // also releases image memory
}

UEyeOpencvCam::~UEyeOpencvCam() {}

void UEyeOpencvCam::threadFrameCapture()
{
    THREAD_ACTIVE = true;

#ifdef _WIN32
    HANDLE hEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
    is_InitEvent(hCam,hEvent,IS_SET_EVENT_FRAME);
#endif

    is_EnableEvent(hCam, IS_SET_EVENT_FRAME);

    while(Parameters::CAMERA_RUNNING && Parameters::ONLINE_PROCESSING)
    {

#ifdef __linux__
        if (is_WaitEvent(hCam, IS_SET_EVENT_FRAME, 1000) == IS_SUCCESS)
#else
        if (WaitForSingleObject(hEvent,1000) == WAIT_OBJECT_0)
#endif

        {
            std::unique_lock<std::mutex> lck(Parameters::frameCaptureMutex);

            if (Parameters::CAMERA_READY)
            {
                VOID* pMem;

                if (is_GetImageMem(hCam, &pMem) == IS_SUCCESS)
                {
                    UEYEIMAGEINFO mUEYEIMAGEINFO;
                    is_GetImageInfo(hCam, pid, &mUEYEIMAGEINFO, sizeof(mUEYEIMAGEINFO));
                    unsigned long long timeStamp = mUEYEIMAGEINFO.u64TimestampDevice;

                    cv::Mat img = cv::Mat(height, width, CV_8UC3);
                    memcpy(img.ptr(), pMem, width * height * 3);

                    if (TRIAL_RECORDING)
                    {
                        vImageInfo[frameIndex].image = img;
                        vImageInfo[frameIndex].time = timeStamp;

                        frameIndex = (frameIndex + 1) % numberOfImageBuffers;
                        frameCount++;

                        if (frameCount >= numberOfImageBuffers)
                        {
                            while (frameCount >= numberOfImageBuffers && Parameters::ONLINE_PROCESSING && TRIAL_RECORDING) Parameters::frameCaptureCV.wait(lck); // wait if image buffer is full
                            frameCount = frameCount % numberOfImageBuffers;
                        }
                    }
                    else
                    {
                        vImageInfo[0].time = timeStamp;
                        vImageInfo[0].image = img;
                    }
                }
            }

            Parameters::frameCaptureCV.notify_one();  // notify waiting thread that new image has arrived
        }
        else
        {
            Parameters::CAMERA_RUNNING = false;
        }
    }

    is_DisableEvent(hCam, IS_SET_EVENT_FRAME);

#ifdef _WIN32
    is_ExitEvent(hCam,IS_SET_EVENT_FRAME);
    CloseHandle(hEvent);
#endif

    is_StopLiveVideo(hCam, IS_FORCE_VIDEO_STOP);

    std::unique_lock<std::mutex> lck(exitMutex);
    THREAD_ACTIVE = false;
    exitCV.notify_one();
}

imageInfo UEyeOpencvCam::getFrame()
{
    imageInfo mImageInfoNew;

    if (Parameters::ONLINE_PROCESSING && Parameters::CAMERA_RUNNING)
    {
        std::unique_lock<std::mutex> lck(Parameters::frameCaptureMutex);

        if (TRIAL_RECORDING)
        {
            if (frameCount <= 0)
            {
                while (frameCount <= 0 && Parameters::ONLINE_PROCESSING && TRIAL_RECORDING) Parameters::frameCaptureCV.wait(lck); // wait for new images to arrive
            }

            int index = frameIndex - frameCount;

            if (index < 0)
            {
                index = index + numberOfImageBuffers;
            }

            if (index >= 0) // condition needed in case of force exit from function
            {
                mImageInfoNew = vImageInfo[index];
                frameCount--;
            }
        }
        else
        {
            mImageInfoNew = vImageInfo[0];
        }
    }

    return mImageInfoNew;
}

void UEyeOpencvCam::startRecording()
{
    frameCount = 0;
    frameIndex = 0;
    TRIAL_RECORDING = true;
}

void UEyeOpencvCam::stopRecording()
{
    TRIAL_RECORDING = false;
}

// Pixel clock

std::vector<int> UEyeOpencvCam::getPixelClockRange()
{
    UINT nRange[3];
    is_PixelClock(hCam, IS_PIXELCLOCK_CMD_GET_RANGE, (void*)nRange, sizeof(nRange));
    std::vector<int> newRange(std::begin(nRange), std::end(nRange));
    return newRange;
}

void UEyeOpencvCam::setPixelClock(int pixelClock)
{
    is_PixelClock(hCam, IS_PIXELCLOCK_CMD_SET, (void*)&pixelClock, sizeof(pixelClock));
}

// Frame rate

std::vector<double> UEyeOpencvCam::getFrameRateRange()
{
    double minFrameTime;
    double maxFrameTime;
    double incFrameTime;

    is_GetFrameTimeRange(hCam, &minFrameTime, &maxFrameTime, &incFrameTime);

    std::vector<double> nRange(3);
    nRange[0] = 1 / maxFrameTime;
    nRange[1] = 1 / minFrameTime;
    nRange[2] = 1 / incFrameTime;

    return nRange;
}

double UEyeOpencvCam::setFrameRate(double FPS)
{
    double newFPS;
    is_SetFrameRate(hCam, FPS, &newFPS);
    return newFPS;
}

// Exposure

std::vector<double> UEyeOpencvCam::getExposureRange()
{
    double nRange[3];

    is_Exposure(hCam, IS_EXPOSURE_CMD_GET_EXPOSURE_RANGE, (void*)nRange, sizeof(nRange));

    std::vector<double> newRange(std::begin(nRange), std::end(nRange));

    return newRange;
}

void UEyeOpencvCam::setExposure(double pExp)
{
    is_Exposure(hCam, IS_EXPOSURE_CMD_SET_EXPOSURE, &pExp, sizeof(pExp));
}

// Black level

std::vector<int> UEyeOpencvCam::getBlackLevelOffsetRange()
{
    IS_RANGE_S32 nRange;
    is_Blacklevel(hCam, IS_BLACKLEVEL_CMD_GET_OFFSET_RANGE, (void*)&nRange, sizeof(nRange));

    std::vector<int> newRange(3);
    newRange[0] = nRange.s32Min;
    newRange[1] = nRange.s32Max;
    newRange[2] = nRange.s32Inc;
    return newRange;
}

void UEyeOpencvCam::setBlackLevelMode(bool FLAG)
{
    if (FLAG)
    {
        INT nMode = IS_AUTO_BLACKLEVEL_ON;
        is_Blacklevel(hCam, IS_BLACKLEVEL_CMD_SET_MODE, (void*)&nMode , sizeof(nMode));
    }
    else
    {
        INT nMode = IS_AUTO_BLACKLEVEL_OFF;
        is_Blacklevel(hCam, IS_BLACKLEVEL_CMD_SET_MODE, (void*)&nMode , sizeof(nMode));
    }
}

int UEyeOpencvCam::setBlackLevelOffset(int nOffset)
{
    is_Blacklevel(hCam, IS_BLACKLEVEL_CMD_SET_OFFSET, (void*)&nOffset, sizeof(nOffset));
    is_Blacklevel(hCam, IS_BLACKLEVEL_CMD_GET_OFFSET, (void*)&nOffset, sizeof(nOffset));
    return nOffset;
}

// Hardware gain

void UEyeOpencvCam::setAutoGain(bool FLAG)
{
    double dEnable;

    if (FLAG) { dEnable = 1; }
    else      { dEnable = 0; }

    is_SetAutoParameter(hCam, IS_SET_ENABLE_AUTO_GAIN, &dEnable, 0);
}

void UEyeOpencvCam::setGainBoost(bool FLAG)
{
    if (FLAG) { is_SetGainBoost (hCam, IS_SET_GAINBOOST_ON);  }
    else      { is_SetGainBoost (hCam, IS_SET_GAINBOOST_OFF); }
}

void UEyeOpencvCam::setHardwareGain(int nMaster)
{
    is_SetHardwareGain(hCam, nMaster, IS_IGNORE_PARAMETER, IS_IGNORE_PARAMETER, IS_IGNORE_PARAMETER);
}

int UEyeOpencvCam::getHardwareGain()
{
    return (is_SetHardwareGain(hCam, IS_GET_MASTER_GAIN, IS_IGNORE_PARAMETER, IS_IGNORE_PARAMETER, IS_IGNORE_PARAMETER));
}
