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
    static AOIProperties eyeAOI;

    static double beadAOIXPosFraction;
    static double beadAOIYPosFraction;
    static AOIProperties beadAOI;

    static drawBooleans drawFlags;

    static int cameraXResolution;
    static int cameraYResolution;
    static AOIProperties cameraAOI;

    static int cannyKernelSize;
    static int ellipseDrawCrossSize;

    static std::condition_variable frameCaptureCV;

    static std::mutex frameCaptureMutex;
    static std::mutex mainMutex;

};

#endif // PARAMETERS_H
