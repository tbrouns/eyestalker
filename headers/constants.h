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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef __linux__
#else
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif


// Parameter initial values

const double aspectRatioMinIni = 0.4;
const double aspectRatioIni    = 0.9;

const double fitErrorFraction = 0.05;

// Parameter limits

const int cameraFrameRateUpperLimit = 500;

// Other

const bool DETECT_BEAD      = false;

const double pupilHaarReductionFactor = 0.75;
const double pupilImageFactor = 1.15;

// Edge processing

const int lengthWindowLength    = 5;
const int curvatureWindowLength = 5; // window length for curvature
const int minimumEdgeLength     = 5;
const int edgeIntensitiesPositionOffset = 2;

// Edge classification

const double scoreFactorRadius        = 0.67;
const double scoreFactorCurvature     = 0.05;
const double scoreFactorCircumference = 0.81;
const double scoreFactorIntensity     = 0.55;
const double scoreFactorGradient      = 0.00; // needs to be measured still

const double scoreThreshold = 0.3324;

#endif
