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

#ifndef ARRAYS
#define ARRAYS

#include <vector>

const std::vector<double> parametersEye  = {0.005,  // Alpha average
                                             0.20,  // Alpha features
                                             0.10,  // Alpha certainty
                                             0.75,  // Alpha position
                                                4,  // Canny blur level
                                                5,  // Canny kernel size
                                            300.0,  // Canny threshold low
                                            600.0,  // Canny threshold high
                                               15,  // Curvature offset
                                             0.40,  // Edge length fraction
                                                3,  // Ellipse fit number maximum
                                              120,  // Ellipse fit error maximum
                                               12,  // Glint size
                                              320,  // Circumference max
                                               60,  // Circumference min
                                              0.4,  // Aspect ratio min
                                             1.20,  // Circumference offset
                                             1.15,  // Circumference change threshold
                                             1.15,  // Aspect ratio change threshold
                                               10,  // Displacement change threshold
                                             0.39,  // Score threshold
                                             0.50}; // Score difference threshold


const std::vector<double> parametersBead = {0.005,  // Alpha average
                                             0.75,  // Alpha miscellaneous
                                             0.10,  // Alpha certainty
                                             0.75,  // Alpha predicted
                                                4,  // Canny blur level
                                                3,  // Canny kernel size
                                            300.0,  // Canny threshold low
                                            600.0,  // Canny threshold high
                                               15,  // Curvature offset
                                             0.40,  // Edge length fraction
                                                3,  // Ellipse fit number maximum
                                               80,  // Ellipse fit error maximum
                                                0,  // Glint size
                                              130,  // Circumference max
                                               90,  // Circumference min
                                              0.8,  // Aspect ratio min
                                             1.10,  // Circumference offset
                                             1.10,  // Circumference change threshold
                                             1.10,  // Aspect ratio change threshold
                                               10,  // Displacement change threshold
                                             0.39,  // Score threshold
                                             0.50}; // Score difference threshold

const std::vector<double> parametersEdgeRadius        = { 0.83431,  0.99374,  0.062139,           0.22416,   0.92573,  0.11619};
const std::vector<double> parametersEdgeCircumference = {  0.9201,   1.0419,   0.41725,           0.59665,     0.542,  0.31648};
const std::vector<double> parametersEdgeCurvature     = { -1.1992, -0.69392,    1.4081, 168571135634303.5, -130.9298,  23.1155};
const std::vector<double> parametersEdgeIntensity     = { 0.92512,   1.3436,    9.6657,          0.069804,   -2.3099,  23.4447};
const std::vector<double> parametersEdgeGradient      = {       0, -16.7695,    2.4122,            1.0705,   -1.6043,   6.7995};
const std::vector<double> parametersEdgeRadiusVar     = {0.036046, 0.002672, 0.0017742,           0.98751, 0.0075968, 0.053661};

const std::vector<double> parametersFitAspectRatio   = {0.73173, 0.0005102, 0.0093987, 0.23184, 0.0005102, 0.033331};
const std::vector<double> parametersFitCircumference = {0.81525,   0.99994, 0.0051646, 0.17624,    1.0089, 0.034909};
const std::vector<double> parametersFitDisplacement  = {0.37568,  0.040816,   0.66938, 0.59611,  0.040816,   3.1728};
const std::vector<double> parametersFitLength        = {0.12018,   0.76041,  0.033629, 0.97603,   0.99961,  0.26897};
const std::vector<double> parametersFitError         = {0.83881,   0.20388,  0.075325, 0.28425,   0.41703,  0.24919};

const std::vector<double> arrayCircumferences = {60,65.918,71.837,77.755,83.673,89.592,95.51,101.43,107.35,113.27,119.18,125.1,131.02,136.94,142.86,148.78,154.69,160.61,166.53,172.45,178.37,184.29,190.2,196.12,202.04,207.96,213.88,219.8,225.71,231.63,237.55,243.47,249.39,255.31,261.22,267.14,273.06,278.98,284.9,290.82,296.73,302.65,308.57,314.49,320.41,326.33,332.24,338.16,344.08,350};
const std::vector<double> arrayAspectRatios   = {0.4,0.41224,0.42449,0.43673,0.44898,0.46122,0.47347,0.48571,0.49796,0.5102,0.52245,0.53469,0.54694,0.55918,0.57143,0.58367,0.59592,0.60816,0.62041,0.63265,0.6449,0.65714,0.66939,0.68163,0.69388,0.70612,0.71837,0.73061,0.74286,0.7551,0.76735,0.77959,0.79184,0.80408,0.81633,0.82857,0.84082,0.85306,0.86531,0.87755,0.8898,0.90204,0.91429,0.92653,0.93878,0.95102,0.96327,0.97551,0.98776,1};


const std::vector<std::vector<double>> arrayCurvatureMax(3);

#endif // ARRAYS

