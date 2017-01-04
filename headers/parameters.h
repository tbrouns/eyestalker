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

#ifndef PARAMETERS_H
#define PARAMETERS_H

// Files

#include "constants.h"
#include "structures.h"

// Libraries

// Standard Template

#include <cmath>
#include <condition_variable>
#include <mutex>              // std::mutex, std::unique_lock

// QT

#include <QColor>

class Parameters
{

public:

    static bool CAMERA_READY;
    static bool CAMERA_RUNNING;
    static bool ONLINE_PROCESSING;
    static double ellipseDrawOutlineWidth;
    static double eyeAOIXPosFraction;
    static double eyeAOIYPosFraction;
    static drawBooleans drawFlags;
    static int cameraAOIHght;
    static int cameraAOIWdth;
    static int cameraAOIXPos;
    static int cameraAOIYPos;
    static int cameraXResolution;
    static int cameraYResolution;
    static int cannyKernelSize;
    static int ellipseDrawCrossSize;
    static int eyeAOIHght;
    static int eyeAOIWdth;
    static int eyeAOIXPos;
    static int eyeAOIYPos;
    static int flashAOIXPos;
    static int flashAOIYPos;
    static int flashAOIWdth;
    static int flashAOIHght;
    static std::condition_variable frameCaptureCV;
    static std::mutex frameCaptureMutex;
    static std::mutex primaryMutex;
    static std::mutex secondaryMutex;

};

#endif // PARAMETERS_H
