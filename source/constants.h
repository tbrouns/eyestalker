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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef __linux__
#else
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

// Parameter initial values

const double initialAspectRatio = 0.9;
const double initialCurvature   =  30;
const double initialIntensity   =  60;

// Other

const double fitMinRange      = 0.30;
const double fitMinEdgeLength = 0.30;
const double scoreStepSize    = 0.20;
const double aspectRatioSlope = 154;

const double frameRateOffset = 1.0; // frequency that camera frame-rate should be within (in Hz)

// Functions

const double windowLengthFraction = 0.04;
const int    windowLengthMin      = 5;

const double certaintyAsymptoteX = 0.50;
const double certaintyAsymptoteY = 0.99;

const double certaintyReduction = 0.90;
const double certaintyOffset    = 0.10;
const double certaintyThreshold = 0.75;
const double certaintyLatency   = 10.0; // turn into adjustable parameter.

#endif
