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

const std::vector<int> arrayCurvatureWindowLengths = { 5, 8, 12, 16};

#endif // ARRAYS

