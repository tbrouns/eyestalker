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

#ifndef STRUCTURES_H
#define STRUCTURES_H

// Libraries

// Standard Template

#include <vector>

// QT

#include <QColor>

// OpenCV

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

struct eyePropertiesParameters
{
    double alphaMiscellaneous;
    double alphaMomentum;
    double alphaAverage;
    double alphaPrediction;
    double curvatureOffset;
    double edgeIntensityOffset;
    double ellipseFitErrorMaximum;
    double pupilCircumferenceMax;
    double pupilCircumferenceMin;
    double pupilFractMin;
    double pupilOffset;
    double thresholdCircumferenceChangeMin;
    double thresholdFractionChangeMin;
    int cannyBlurLevel;
    int cannyKernelSize;
    int cannyLowerLimit;
    int cannyUpperLimit;
    int edgeMaximumFitNumber;
    int glintRadius;
};

struct eyePropertiesVariables
{
    bool pupilDetected;
    double edgeIntensityPrediction;
    double edgeIntensityAverage;
    double momentumCircumference;
    double momentumFraction;
    double momentumRadius;
    double pupilCircumferencePrediction;
    double pupilCircumferenceAverage;
    double pupilCircumferenceExact;
    double pupilFractionPrediction;
    double pupilFractionAverage;
    double pupilFractionExact;
    double pupilRadiusPrediction;
    double searchRadius;
    double thresholdCircumferenceChange;
    double thresholdFractionChange;
    double xPosExact;
    double xPosAbsolute;
    double xPosPredicted;
    double xVelocity;
    double yPosExact;
    double yPosAbsolute;
    double yPosPredicted;
    double yVelocity;
};

struct eyePropertiesMiscellaneous
{
    bool errorDetected;
    cv::Mat image;
    int glintHaarWidth;
    int glintHaarXPos;
    int glintHaarYPos;
    int offsetPupilHaarWdth;
    int offsetPupilHaarXPos;
    int offsetPupilHaarYPos;
    int pupilHaarWdth;
    int pupilHaarXPos;
    int pupilHaarYPos;
    std::vector<char> cannyEdges;
    std::vector<double> ellipseCoefficients;
    std::vector<std::vector<int>> edgeIndices;
    std::vector<std::vector<int>> edgeIndicesNew;
};

struct eyeProperties
{
    eyePropertiesParameters     p;
    eyePropertiesVariables      v;
    eyePropertiesMiscellaneous  m;
};

struct haarProperties
{
    int x_pos;
    int y_pos;
};

struct edgeProperties
{
    std::vector<double> curvatures;
    std::vector<double> intensities;
    std::vector<int> lengths;
    std::vector<int> sizes;
    std::vector<std::vector<int>> indices;
};

struct ellipseProperties
{
    bool pupilDetected;
    double circumference;
    double fraction;
    double intensity;
    double radius;
    double xPos;
    double yPos;
    std::vector<double> coefficients;
};

struct imageInfo
{
    unsigned long long time;
    cv::Mat image;
};

struct drawBooleans
{
    bool haar;
    bool edge;
    bool elps;
};

#endif // STRUCTURES_H
