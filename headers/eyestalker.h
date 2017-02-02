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

#ifndef PUPILDETECTION
#define PUPILDETECTION

// Files

#include "arrays.h"
#include "constants.h"
#include "structures.h"

// Libraries

// OpenCV library

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>



// Eigen library

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Standard Template

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring> // used for saving images
#include <fstream>
#include <iostream>
#include <numeric> // used for 'accumulate'
#include <stdio.h>
#include <string>
#include <vector>

/// general functions

double calculateMean(const std::vector<double>&);

/// eye-tracking functions

detectionProperties pupilDetection(const cv::Mat&, detectionProperties);
double flashDetection(const cv::Mat&);

#endif // PUPILDETECTION
