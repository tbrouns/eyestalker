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

const double initialAspectRatio = 0.9;
const double initialCurvature   =  30;
const double initialIntensity   =  60;

const double fitErrorFraction   = 0.05;

// Other

// Edge processing

const int lengthWindowLength    = 5;
const int curvatureWindowLength = 5; // window length for curvature
const int minimumEdgeLength     = 5;

// Edge classification

const double scoreFactorCircumference = 0.50;
const double scoreFactorIntensity     = 0.75;
const double scoreFactorGradient      = 0.60;
const double scoreFactorCurvature     = 0.80;
const double scoreFactorRadius        = 0.65;
const double scoreFactorRadiusVar     = 0.65;

// Functions

const double certaintyAsymptoteX = 0.50;
const double certaintyAsymptoteY = 0.99;

const double certaintyLatency = 10.0; // turn into adjustable parameter.

const double certaintyOffset = 0.1;

#endif
