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

#include "headers/parameters.h"

bool Parameters::CAMERA_READY;
bool Parameters::CAMERA_RUNNING;
bool Parameters::ONLINE_PROCESSING;

double Parameters::eyeAOIXPosFraction;
double Parameters::eyeAOIYPosFraction;
int    Parameters::eyeAOIHght;
int    Parameters::eyeAOIWdth;
int    Parameters::eyeAOIXPos;
int    Parameters::eyeAOIYPos;

double Parameters::beadAOIXPosFraction;
double Parameters::beadAOIYPosFraction;
int    Parameters::beadAOIHght;
int    Parameters::beadAOIWdth;
int    Parameters::beadAOIXPos;
int    Parameters::beadAOIYPos;

int Parameters::cameraAOIHght;
int Parameters::cameraAOIWdth;
int Parameters::cameraAOIXPos;
int Parameters::cameraAOIYPos;
int Parameters::cameraXResolution;
int Parameters::cameraYResolution;
int Parameters::cannyKernelSize;

double Parameters::ellipseDrawOutlineWidth;
int    Parameters::ellipseDrawCrossSize;

drawBooleans Parameters::drawFlags;

std::condition_variable Parameters::frameCaptureCV;

std::mutex Parameters::frameCaptureMutex;
std::mutex Parameters::mainMutex;
