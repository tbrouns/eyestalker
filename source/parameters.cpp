//  EyeStalker: robust video-based eye tracking
//  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

//  EyeStalker is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  EyeStalker is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#include "parameters.h"

bool Parameters::CAMERA_READY;
bool Parameters::CAMERA_RUNNING;
bool Parameters::ONLINE_PROCESSING;

AOIProperties Parameters::eyeAOI;
AOIPropertiesDouble Parameters::eyeAOIRatio;

AOIProperties Parameters::beadAOI;
AOIPropertiesDouble Parameters::beadAOIRatio;

int Parameters::cameraXResolution;
int Parameters::cameraYResolution;
AOIProperties Parameters::camAOI;

double Parameters::ellipseDrawOutlineWidth;
int    Parameters::ellipseDrawCrossSize;

drawBooleans Parameters::drawFlags;

std::condition_variable Parameters::frameCaptureCV;

std::mutex Parameters::frameCaptureMutex;

std::mutex Parameters::AOICamMutex;
std::mutex Parameters::AOIEyeMutex;
std::mutex Parameters::AOIBeadMutex;
