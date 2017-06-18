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

#include <vector>

// Parameter initial values

static const std::vector<double> parametersEye  = {0.005,   //  0. Gain average
                                                   0.40,    //  1. Gain features
                                                   0.25,    //  2. Gain certainty
                                                   0.75,    //  3. Gain position
                                                   4,       //  4. Canny blur level
                                                   5,       //  5. Canny kernel size
                                                   300.0,   //  6. Canny threshold low
                                                   600.0,   //  7. Canny threshold high
                                                   2.5,     //  8. Curvature offset
                                                   0.05,    //  9. Ellipse edge fraction
                                                   4,       // 10. Maximum number of edges
                                                   0.60,    // 11. Maximum fit error
                                                   12,      // 12. Glint size
                                                   290,     // 13. Circumference max
                                                   60,      // 14. Circumference min
                                                   0.4,     // 15. Aspect ratio min
                                                   0.12,    // 16. Circumference change threshold upper
                                                   0.03,    // 17. Circumference change threshold lower
                                                   0.09,    // 18. Aspect ratio  change threshold upper
                                                   0.03,    // 19. Aspect ratio  change threshold lower
                                                   6,       // 20. Displacement  change threshold upper
                                                   3,       // 21. Displacement  change threshold lower
                                                   0.30,    // 22. Score threshold edge
                                                   0.10,    // 23. Score threshold fit
                                                   0.60,    // 24. Score difference threshold edge
                                                   0.10,    // 25. Score difference threshold fit
                                                   7,       // 26. Edge window length
                                                   6};      // 27. Maximum number of fits

const double initialAspectRatio  = 0.9;
const double initialCurvature    =  30;
const double initialIntensity    =  60;
const double initialHaarResponse =   0;
const double initialGradient     =   0;

// Should probably be adjustable parameters

const double fitMinRange      = 0.30;
const double fitMinEdgeLength = 0.30;
const double scoreStepSize    = 0.20;
const double aspectRatioSlope = 154;

// Can be static constants

const double windowLengthFraction = 0.04;
const int    windowLengthMin      = 5;

const double certaintyAsymptoteX = 0.50;
const double certaintyAsymptoteY = 0.99;

const double certaintyReduction = 0.90;
const double certaintyOffset    = 0.10;
const double certaintyThreshold = 0.75;
const double certaintyLatency   = 10.0; // turn into adjustable parameter.

#endif
