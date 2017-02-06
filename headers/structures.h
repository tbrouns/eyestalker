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

struct AOIProperties
{
    int xPos;
    int yPos;
    int wdth;
    int hght;
};

struct edgeProperties
{
    double curvatureMax;
    double curvatureMin;
    double curvature;
    double intensity;
    double gradient;
    double radius;
    double radiusVar;
    int index;
    int length;
    int size;
    int tag;
    std::vector<int> pointIndices;
    std::vector<int> gradients;
    std::vector<int> intensities;
    std::vector<double> radii;
    std::vector<double> curvatures;
    std::vector<double> xnormals;
    std::vector<double> ynormals;
};

struct ellipseProperties
{
    bool DETECTED;
    double fitError;
    double intensity;
    double gradient;
    double circumference;
    double aspectRatio;
    double radius;
    double xPos;
    double yPos;
    double width;
    double height;
    int edgeLength;
    int tag;
    std::vector<double> coefficients;
    std::vector<int> edgeIndices;

};

struct detectionParameters
{
    detectionParameters(): DETECTION_ON(false) { }

    double alphaFeatures;
    double alphaAverage;
    double alphaPosition;
    double alphaCertainty;
    double curvatureOffset;
    double edgeLengthFraction;
    double ellipseFitErrorMaximum;
    double circumferenceMax;
    double circumferenceMin;
    double circumferenceOffset;
    double aspectRatioMin;
    double changeThresholdCircumference;
    double changeThresholdAspectRatio;
    double changeThresholdPosition;
    double scoreThreshold;
    double scoreThresholdPoints;
    int AOIXPos;
    int AOIYPos;
    int AOIWdth;
    int AOIHght;
    int cannyBlurLevel;
    int cannyKernelSize;
    double cannyThresholdLow;
    double cannyThresholdHigh;
    int ellipseFitNumberMaximum;
    int glintWdth;
    bool DETECTION_ON;
};

struct detectionVariables
{
    double averageAspectRatio;
    double averageCircumference;
    double averageHeight;
    double averageWidth;
    double averageIntensity;
    double averageGradient;
    double certaintyAverage;
    double certaintyFeatures;
    double certaintyPosition;
    double changeThresholdAspectRatio;
    double changeThresholdCircumference;
    double changeThresholdPosition;
    double offsetCircumference;
    double momentumAspectRatio;
    double momentumCircumference;
    double momentumHeight;
    double momentumWidth;
    double momentumXPos;
    double momentumYPos;
    double predictedAspectRatio;
    double predictedCircumference;
    double predictedHeight;
    double predictedWidth;
    double predictedXPos;
    double predictedYPos;
    double predictedXPosRelative;
    double predictedYPosRelative;
    double predictedCurvature;
};

struct dataVariables
{
    bool DETECTED;
    double absoluteXPos;
    double absoluteYPos;
    std::vector<edgeProperties> edgeData;
    std::vector<ellipseProperties> ellipseData;
    double duration;
    double exactAspectRatio;
    double exactCircumference;
    double exactXPos;
    double exactYPos;
    double timestamp;
};

struct drawVariables
{
    bool DETECTED;
    int exactXPos;
    int exactYPos;
    int predictedXPos;
    int predictedYPos;
    AOIProperties glintAOI;
    AOIProperties innerAOI;
    AOIProperties outerAOI;
    std::vector<int> cannyEdgeIndices;
    std::vector<double> ellipseCoefficients;
    std::vector<edgeProperties> edgePropertiesAll;
};

struct detectionProperties
{
    detectionParameters     p;
    detectionVariables      v;
};


struct vertexProperties
{
    int index;
    int pointIndex;
    std::vector<int> connectedBranches;
    std::vector<int> connectedPoints;
};

struct branchProperties
{
    int index;
    std::vector<int> pointIndices;
    std::vector<int> connectedVertices;
    int length;
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
