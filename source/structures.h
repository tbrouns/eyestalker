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

struct AOIPropertiesDouble
{
    double xPos;
    double yPos;
    double wdth;
    double hght;
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
    double length;
    double score;
    int index;
    int tag;
    std::vector<int> edgeIndices;
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
    double angle;
    double curvature;
    double fitError;
    double intensity;
    double gradient;
    double circumference;
    double aspectRatio;
    double xPos;
    double yPos;
    double width;
    double height;
    double edgeLength;
    double edgeScore;
    int tag = 0;
    std::vector<double> coefficients;
    std::vector<int> edgeIndices;

};

struct detectionParameters
{
    detectionParameters(): DETECTION_ON(false) { }

    bool DETECTION_ON;

    double alphaAverages;
    double alphaCertainty;
    double alphaFeatures;
    double alphaPosition;

    double cannyThresholdHigh;
    double cannyThresholdLow;
    int    cannyBlurLevel;
    int    cannyKernelSize;

    int    glintWdth;

    double curvatureOffset;

    double thresholdAspectRatioMin;
    double thresholdCircumferenceMax;
    double thresholdCircumferenceMin;
    double thresholdChangeAspectRatio;
    double thresholdChangeCircumference;
    double thresholdChangePosition;
    double thresholdFitError;
    double thresholdScoreEdge;
    double thresholdScoreFit;
    double thresholdScoreDiffEdge;
    double thresholdScoreDiffFit;

    int    fitMaximum;
    int    fitEdgeMaximum;
    double fitEdgeFraction;

    int windowLengthEdge;
};

struct detectionVariables
{
    double averageAspectRatio;
    double averageCircumference;
    double averageCurvature;
    double averageHeight;
    double averageWidth;
    double averageIntensity;
    double averageGradient;
    double averageHaarResponse;
    double certaintyAverages;
    double certaintyFeatures;
    double certaintyPosition;
    double certaintyAveragesPrime;
    double certaintyFeaturesPrime;
    double certaintyPositionPrime;
    double momentumAspectRatio;
    double momentumCircumference;
    double momentumHeight;
    double momentumWidth;
    double momentumXPos;
    double momentumYPos;
    double momentumGradient;
    double momentumIntensity;
    double momentumCurvature;
    double momentumHaarResponse;
    double predictedAngle;
    double predictedAspectRatio;
    double predictedCircumference;
    double predictedCurvature;
    double predictedIntensity;
    double predictedGradient;
    double predictedHaarResponse;
    double predictedHeight;
    double predictedWidth;
    double predictedXPos;
    double predictedYPos;
    double predictedXPosRelative;
    double predictedYPosRelative;
    double thresholdAspectRatioMax;
    double thresholdAspectRatioMin;
    double thresholdCircumferenceMin;
    double thresholdCircumferenceMax;
    double thresholdChangeAspectRatio;
    double thresholdChangeCircumference;
    double thresholdChangePosition;
    double thresholdScoreEdge;
    double thresholdScoreFit;
    int windowLengthEdge;
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

struct developmentOptions
{
    developmentOptions(): CURVATURE_MEASUREMENT(false) { }

    bool CURVATURE_MEASUREMENT;
};

struct drawVariables
{
    bool DETECTED;
    bool PROCESSED;
    int exactXPos;
    int exactYPos;
    int predictedXPos;
    int predictedYPos;
    AOIProperties glintAOI;
    AOIProperties haarAOI;
    AOIProperties cannyAOI;
    std::vector<int> cannyEdgeIndices;
    std::vector<double> ellipseCoefficients;
    std::vector<edgeProperties> edgeData;
};

struct vertexProperties
{
    int tag;
    std::vector<int> pointIndices;
    std::vector<int> connectedArcs;
    std::vector<int> connectedPoints;
};

struct arcProperties
{
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
