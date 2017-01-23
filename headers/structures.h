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

struct haarProperties
{
    int xPos;
    int yPos;
};

struct edgeProperties
{
    double curvatureMax;
    double curvatureMin;
    double curvatureAvg;
    double intensity;
    double distance;
    int index;
    int length;
    int size;
    int flag;
    std::vector<int> pointIndices;
};

struct ellipseProperties
{
    bool pupilDetected;
    double fitError;
    double circumference;
    double aspectRatio;
    double intensity;
    double radius;
    double xPos;
    double yPos;
    double width;
    double height;
    int edgeLength;
    std::vector<double> coefficients;
    std::vector<int> edgeIndices;

};

struct eyePropertiesParameters
{
    double alphaMiscellaneous;
    double alphaMomentum;
    double alphaAverage;
    double alphaPrediction;
    double curvatureFactor;
    double curvatureOffsetMin;
    double edgeIntensityOffset;
    double ellipseFitErrorMaximum;
    double circumferenceMax;
    double circumferenceMin;
    double aspectRatioMin;
    double circumferenceChangeThreshold;
    double aspectRatioChangeThreshold;
    int cannyBlurLevel;
    int cannyKernelSize;
    int cannyThresholdLow;
    int cannyThresholdHigh;
    int edgeMaximumFitNumber;
    int glintSize;
    int pupilOffset;
};

struct eyePropertiesVariables
{
    bool pupilDetected;
    double aspectRatioAverage;
    double aspectRatioExact;
    double aspectRatioMomentum;
    double aspectRatioPrediction;
    double circumferenceAverage;
    double circumferenceExact;
    double circumferenceMomentum;
    double circumferencePrediction;
    double curvatureOffset;
    double edgeCurvaturePrediction;
    double edgeIntensityAverage;
    double edgeIntensityPrediction;
    double heightAverage;
    double heightMomentum;
    double heightPrediction;
    double radiusMomentum;
    double radiusPrediction;
    double searchRadius;
    double thresholdAspectRatioChange;
    double thresholdCircumferenceChange;
    double priorCertainty;
    double widthAverage;
    double widthMomentum;
    double widthPrediction;
    double xPosAbsolute;
    double xPosExact;
    double xPosPredicted;
    double xVelocity;
    double yPosAbsolute;
    double yPosExact;
    double yPosPredicted;
    double yVelocity;
};

struct eyePropertiesMiscellaneous
{
    bool errorDetected;
    cv::Mat image;
    cv::Mat imagePupil;
    int glintSize;
    int glintXPos;
    int glintYPos;
    int offsetPupilHaarWdth;
    int offsetPupilHaarHght;
    int offsetPupilHaarXPos;
    int offsetPupilHaarYPos;
    int pupilHaarHght;
    int pupilHaarWdth;
    int pupilHaarXPos;
    int pupilHaarYPos;
    std::vector<int> cannyEdges;
    std::vector<double> ellipseCoefficients;
    std::vector<edgeProperties> edgePropertiesAll;
};

struct eyeProperties
{
    eyePropertiesParameters     p;
    eyePropertiesVariables      v;
    eyePropertiesMiscellaneous  m;
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
