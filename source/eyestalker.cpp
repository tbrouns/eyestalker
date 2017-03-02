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

#include "eyestalker.h"

// General

double calculateMean(const std::vector<double>& v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return (sum / v.size());
}

double calculateMeanInt(const std::vector<int>& v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return (sum / v.size());
}

double calculateVariance(const std::vector<double>& v)
{
    int size    = v.size();
    double mean = calculateMean(v);
    double temp = 0;

    for (int i = 0; i < size; i++)
    {
        double val = v[i];
        temp += (val - mean) * (val - mean);
    }

    return (temp / size);
}

double ramanujansApprox(double a, double b) // ramanujans 2nd approximation
{
    double h = pow((a - b), 2) / pow((a + b), 2);
    return (M_PI * (a + b) * (1 + (3 * h) / (10 + sqrt(4 - 3 * h))));
}


double ceil2( double value )
{
    if (value < 0.0) { return floor(value); }
    else             { return  ceil(value); }
}

inline int findUpperBound(const std::vector<double>& v, const double& val)
{
    int index = 0;
    int vSize = v.size();
    auto itr = std::upper_bound(v.begin(), v.end(), val);
    if      (itr == v.begin()) { index = 0; }
    else if (itr == v.end())   { index = vSize - 1;}
    else                       { index = itr - v.begin(); }
    return index;
}

// Scores

inline double calculateCertainty(double x, double midpoint)
{
    double a = certaintyAsymptoteX * midpoint;
    double k = - log((1 / certaintyAsymptoteY) - 1) / a;
    return (1 - 2 / (1 + exp(-k * (x - midpoint)))); // Logistic
}

inline double calculateGaussian2(double x, std::vector<double> p)
{
    return (p[0] * exp(-pow((x - p[1]) / p[2], 2)) + p[3] * exp(-pow((x - p[4]) / p[5], 2))); // Two-term Gaussian
}

double calculateScoreTotal(const detectionVariables& mDetectionVariables, std::vector<double>& inputVector, bool USE_CERTAINTY)
{   
    static const double scoreFactorCircumference = 0.25;
    static const double scoreFactorIntensity     = 0.66;
    static const double scoreFactorCurvature     = 0.42;
    static const double scoreFactorRadius        = 0.65;
    static const double scoreFactorRadiusVar     = 0.21;
    static const double scoreFactorGradient      = 0.24;
    static const double scoreFactorLength        = 0.49;

    static const std::vector<double> parametersEdgeRadius        = { 0.83431,  0.99374,  0.062139,           0.22416,   0.92573,  0.11619};
    static const std::vector<double> parametersEdgeCircumference = {  0.9201,   1.0419,   0.41725,           0.59665,     0.542,  0.31648};
    static const std::vector<double> parametersEdgeCurvature     = { -1.1992, -0.69392,    1.4081, 168571135634303.5, -130.9298,  23.1155};
    static const std::vector<double> parametersEdgeIntensity     = { 0.92512,   1.3436,    9.6657,          0.069804,   -2.3099,  23.4447};
    static const std::vector<double> parametersEdgeGradient      = {       0, -16.7695,    2.4122,            1.0705,   -1.6043,   6.7995};
    static const std::vector<double> parametersEdgeRadiusVar     = {0.036046, 0.002672, 0.0017742,           0.98751, 0.0075968, 0.053661};

    for (int i = 0, vSize = inputVector.size(); i < vSize; i++) // check for NaNs or Infs
    {
        double val = inputVector[i];
        if (!std::isfinite(val)) { inputVector[i] = 0; }
        else if (val < 0)        { inputVector[i] = std::abs(val); }
    }

    // Do score calculation

    double ratioRadius        = inputVector[0];
    double ratioCircumference = inputVector[1];
    double errorCurvature     = inputVector[2];
    double errorIntensity     = inputVector[3];
    double errorGradient      = inputVector[4];
    double varianceRadius     = inputVector[5];
    double edgeLength         = inputVector[6];

    errorIntensity = std::abs(errorIntensity);
    errorGradient  = std::abs(errorGradient);
    errorCurvature = std::abs(errorCurvature);

    double certaintyFactorPosition = mDetectionVariables.certaintyPosition;
    double certaintyFactorFeatures = mDetectionVariables.certaintyFeatures;
    double certaintyFactorAverages = mDetectionVariables.certaintyAverages;

    if (certaintyFactorAverages > certaintyFactorFeatures + certaintyOffset) { certaintyFactorFeatures = certaintyFactorAverages; }

    double certaintyFactorGradient = certaintyFactorPosition * certaintyFactorFeatures;

    double certaintyFactorIntensity = 1.0;
    double certaintyFactorCurvature = 1.0;

    if (USE_CERTAINTY) // false when comparing two edges; true when comparing edge characteristics with predictions
    {
        certaintyFactorIntensity = certaintyFactorFeatures;
        certaintyFactorCurvature = certaintyFactorFeatures;
    }

    double factorRadius        = certaintyFactorPosition  * scoreFactorRadius;
    double factorCircumference = certaintyFactorFeatures  * scoreFactorCircumference;
    double factorGradient      = certaintyFactorGradient  * scoreFactorGradient;
    double factorIntensity     = certaintyFactorIntensity * scoreFactorIntensity;

    // Importance of curvature and radial variance should be dependent on length of edge

    double factorLength = scoreFactorLength * (edgeLength / mDetectionVariables.predictedCircumference) + (1 - scoreFactorLength);

    double factorRadiusVar = certaintyFactorPosition  * scoreFactorRadiusVar * factorLength;
    double factorCurvature = certaintyFactorCurvature * scoreFactorCurvature * factorLength;

    // Calculate scores

    double scoreRadius        = factorRadius        * calculateGaussian2(ratioRadius,        parametersEdgeRadius       );
    double scoreCircumference = factorCircumference * calculateGaussian2(ratioCircumference, parametersEdgeCircumference);
    double scoreCurvature     = factorCurvature     * calculateGaussian2(errorCurvature,     parametersEdgeCurvature    );
    double scoreIntensity     = factorIntensity     * calculateGaussian2(errorIntensity,     parametersEdgeIntensity    );
    double scoreGradient      = factorGradient      * calculateGaussian2(errorGradient,      parametersEdgeGradient     );
    double scoreRadiusVar     = factorRadiusVar     * calculateGaussian2(varianceRadius,     parametersEdgeRadiusVar    );

    const double norm =  factorRadius + factorRadiusVar + factorCurvature + factorCircumference + factorIntensity + factorGradient;

    double scoreTotal = 0;
    if (norm > 0) { scoreTotal = (scoreRadius + scoreRadiusVar + scoreGradient + scoreCurvature + scoreCircumference + scoreIntensity) / norm; }
    return scoreTotal;
}

// Detection

std::vector<unsigned int> calculateIntImg(const cv::Mat& img, AOIProperties searchAOI)
{
    int imgWidth = img.cols;

    int startX = searchAOI.xPos;
    int startY = searchAOI.yPos;
    int width  = searchAOI.wdth;
    int height = searchAOI.hght;

    uchar *ptr = img.data;

    std::vector<unsigned int> integralImage(width * height); // unsigned due to large positive values

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int i = width * y + x; // integral image coordinates
            int j = imgWidth * (y + startY) + (x + startX); // image coordinates

            double val = ptr[j];

            if (x == 0 && y == 0) { integralImage[i] = val; }                            // first point
            else if (y == 0)      { integralImage[i] = val + integralImage[i - 1]; }     // first row
            else if (x == 0)      { integralImage[i] = val + integralImage[i - width]; } // first column
            else                  { integralImage[i] = val + integralImage[i - 1] + integralImage[i - width] - integralImage[i - width - 1]; }
        }
    }

    return integralImage;
}

AOIProperties detectGlint(const cv::Mat& img, AOIProperties searchAOI, AOIProperties glintAOI)
{
    int imgWidth = img.cols;

    int numDirections  = 4;
    int glintThreshold = 200;

    int gradientWindowLength = glintAOI.wdth;
    int glintRadius = round(0.5 * glintAOI.wdth);

    uchar *ptr = img.data;

    std::vector<double> imageGradient(searchAOI.wdth * searchAOI.hght, 0.0);

    std::vector<int> dZ(numDirections);
    dZ[0] = -searchAOI.wdth - 1;
    dZ[1] = -searchAOI.wdth + 1;
    dZ[2] =  searchAOI.wdth + 1;
    dZ[3] =  searchAOI.wdth - 1;

    for (int y = gradientWindowLength; y < searchAOI.hght - gradientWindowLength; y++)
    {
        for (int x = gradientWindowLength; x < searchAOI.wdth - gradientWindowLength; x++)
        {
            int i = searchAOI.wdth * y + x; // gradient coordinates
            int j = imgWidth * (y + searchAOI.yPos) + (x + searchAOI.xPos); // image coordinates

            int centreIntensity = ptr[j];

            if (centreIntensity > glintThreshold)
            {
                double surroundSum = 0;
                for (int m = 0; m < numDirections; m++) { surroundSum += ptr[j + gradientWindowLength * dZ[m]]; }
                imageGradient[i] = centreIntensity / surroundSum;
            }
        }
    }

    int glintIndex = std::distance(imageGradient.begin(), std::max_element(imageGradient.begin(), imageGradient.end()));

    int x = glintIndex % searchAOI.wdth;
    int y = (glintIndex - x) / searchAOI.wdth;

    glintAOI.xPos = x - glintRadius;
    glintAOI.yPos = y - glintRadius;

    return glintAOI;
}

AOIProperties detectPupilApprox(const std::vector<unsigned int>& I, AOIProperties searchAOI, AOIProperties haarAOI, AOIProperties glintAOI, std::vector<double>& innerIntensities, std::vector<double>& outerIntensitiesLeft, std::vector<double>& outerIntensitiesRght, int imgWidth)
{
    double intensityInnerMin = std::numeric_limits<double>::max(); // set to maximum double value;

    haarAOI.xPos = 0;
    haarAOI.yPos = 0;

    int innerArea = (haarAOI.wdth - 1) * (haarAOI.hght - 1);

    int wdth = searchAOI.wdth - haarAOI.wdth;
    int hght = searchAOI.hght - haarAOI.hght;

    for (int iRow = 0; iRow < hght; iRow++)
    {
        for (int iCol = 0; iCol < wdth; iCol++)
        {
            // vertices of inner square

            int xTopLeft = iCol;
            int yTopLeft = iRow;

            int xBtmRght = xTopLeft + haarAOI.wdth - 1;
            int yBtmRght = yTopLeft + haarAOI.hght - 1;

            int iTopLeft = searchAOI.wdth * yTopLeft + xTopLeft;
            int iTopRght = searchAOI.wdth * yTopLeft + xBtmRght;
            int iBtmLeft = searchAOI.wdth * yBtmRght + xTopLeft;
            int iBtmRght = searchAOI.wdth * yBtmRght + xBtmRght;

            // calculate glint intensity

            double glintIntensity = 0.0;
            double glintArea = 0.0;

            bool GLINT_OVERLAP = false; // flag for glint overlap

            int xTopLeftGlint = glintAOI.xPos;
            int yTopLeftGlint = glintAOI.yPos;
            int xBtmRghtGlint = glintAOI.xPos + glintAOI.wdth - 1;
            int yBtmRghtGlint = glintAOI.yPos + glintAOI.hght - 1;

            std::vector<int> X(2);
            X[0] = xTopLeftGlint;
            X[1] = xBtmRghtGlint;

            std::vector<int> Y(2);
            Y[0] = yTopLeftGlint;
            Y[1] = yBtmRghtGlint;

            // check if glint overlaps with Haar detector

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    if (X[i] > xTopLeft && X[i] < xBtmRght && Y[j] > yTopLeft && Y[j] < yBtmRght)
                    { GLINT_OVERLAP = true; }
                }
            }

            if (GLINT_OVERLAP) // if yes, check how much it overlaps
            {
                // Check overlap with edges

                if      (xTopLeftGlint < xTopLeft) { xTopLeftGlint = xTopLeft; } // left edge
                else if (xBtmRghtGlint > xBtmRght) { xBtmRghtGlint = xBtmRght; } // right edge
                if      (yTopLeftGlint < yTopLeft) { yTopLeftGlint = yTopLeft; } // top edge
                else if (yBtmRghtGlint > yBtmRght) { yBtmRghtGlint = yBtmRght; } // bottom edge

                // coordinates of corners of glint square

                int iTopLeftGlint = searchAOI.wdth * yTopLeftGlint + xTopLeftGlint;
                int iTopRghtGlint = searchAOI.wdth * yTopLeftGlint + xBtmRghtGlint;
                int iBtmLeftGlint = searchAOI.wdth * yBtmRghtGlint + xTopLeftGlint;
                int iBtmRghtGlint = searchAOI.wdth * yBtmRghtGlint + xBtmRghtGlint;

                // calculate area and intensity of glint

                glintIntensity = I[iBtmRghtGlint] - I[iBtmLeftGlint] - I[iTopRghtGlint] + I[iTopLeftGlint];

                int glintWdth = xBtmRghtGlint - xTopLeftGlint;
                int glintHght = yBtmRghtGlint - yTopLeftGlint;
                glintArea = glintWdth * glintHght;
            }

            double intensityInner = I[iBtmRght] - I[iBtmLeft] - I[iTopRght] + I[iTopLeft];
            intensityInner = intensityInner - glintIntensity; // adjust for glint
            intensityInner = intensityInner / (innerArea - glintArea);

            int xTopLeftOuter = round(xTopLeft - 0.5 * haarAOI.wdth);
            int xBtmRghtOuter = round(xBtmRght + 0.5 * haarAOI.wdth);

            if (xTopLeftOuter <               0) { xTopLeftOuter =                  0; }
            if (xBtmRghtOuter >= searchAOI.wdth) { xBtmRghtOuter = searchAOI.wdth - 1; }

            int iTopLeftOuter = searchAOI.wdth * yTopLeft + xTopLeftOuter;
            int iTopRghtOuter = searchAOI.wdth * yTopLeft + xBtmRghtOuter;
            int iBtmLeftOuter = searchAOI.wdth * yBtmRght + xTopLeftOuter;
            int iBtmRghtOuter = searchAOI.wdth * yBtmRght + xBtmRghtOuter;

            double intensityOuterLeft = I[iBtmLeft]      - I[iBtmLeftOuter] - I[iTopLeft]      + I[iTopLeftOuter];
            double intensityOuterRght = I[iBtmRghtOuter] - I[iBtmRght]      - I[iTopRghtOuter] + I[iTopRght];

            int outerWdthLeft = xTopLeft - xTopLeftOuter;
            int outerWdthRght = xBtmRghtOuter - xBtmRght;

            if (outerWdthLeft > 0)
            {
                int outerAreaLeft  = outerWdthLeft * (haarAOI.hght - 1);
                intensityOuterLeft = intensityOuterLeft / outerAreaLeft;
            }
            else { intensityOuterLeft = 0; }

            if (outerWdthRght > 0)
            {
                int outerAreaRght  = outerWdthRght * (haarAOI.hght - 1);
                intensityOuterRght = intensityOuterRght / outerAreaRght;
            }
            else { intensityOuterRght = 0; }

            if (intensityInner < intensityInnerMin)
            {
                haarAOI.xPos = xTopLeft;
                haarAOI.yPos = yTopLeft;
                intensityInnerMin = intensityInner;
            }

            // This section needs to be removed

            int x = iCol;
            int y = iRow;

            x = x + searchAOI.xPos;
            y = y + searchAOI.yPos;

            int i = imgWidth * y + x;

            innerIntensities[i]     = intensityInner;
            outerIntensitiesLeft[i] = intensityOuterLeft;
            outerIntensitiesRght[i] = intensityOuterRght;
        }
    }

    return haarAOI;
}

std::vector<int> sharpenEdges(std::vector<int>& binaryImageVectorRaw, AOIProperties mAOI)
{
    std::vector<int> binaryImageVector = binaryImageVectorRaw;

    std::vector<int> dX = {  0,  1, -1,  1,  0, -1,  1, -1};
    std::vector<int> dY = { -1,  1,  0, -1,  1, -1,  0,  1};

    for (int yCentre = 0; yCentre < mAOI.hght; yCentre++)
    {
        for (int xCentre = 0; xCentre < mAOI.wdth; xCentre++)
        {
            int iCentre = mAOI.wdth * yCentre + xCentre;

            if (binaryImageVector[iCentre] == 1)
            {
                for (int m = 0; m < 4; m++)
                {
                    int numFilledPixels = 0;

                    // check combination of two neighbouring pixels in 4-connected environment

                    for (int n = 0; n < 2; n++) // loop through two connected neighbouring pixels
                    {
                        int q = 2 * (m + n) % 8;

                        int xNeighbour = xCentre + dX[q];
                        int yNeighbour = yCentre + dY[q];

                        if (xNeighbour < 0 || xNeighbour >= mAOI.wdth ||yNeighbour < 0 || yNeighbour >= mAOI.hght)
                        {
                            continue; // neighbour is out-of-bounds
                        }

                        int iNeighbour = mAOI.wdth * yNeighbour + xNeighbour;

                        if (binaryImageVector[iNeighbour] == 1) // check if neighbour is filled
                        {
                            numFilledPixels++;
                        }
                    }

                    if (numFilledPixels == 2) // if two neighbouring pixels in 4-connected environment are filled ...
                    {
                        int q = 2 * m + 1;

                        int xOpposite = xCentre + dX[q];
                        int yOpposite = yCentre + dY[q];
                        int iOpposite = mAOI.wdth * yOpposite + xOpposite;

                        if
                                (xOpposite < 0 || xOpposite >= mAOI.wdth || // ... AND opposite pixel is out-of-bounds
                                 yOpposite < 0 || yOpposite >= mAOI.hght ||
                                 binaryImageVector[iOpposite] ==  0      || // ... OR unfilled ...
                                 binaryImageVector[iOpposite] == -1)
                        {
                            binaryImageVector[iCentre] = -1; // ... THEN remove pixel from edge
                        }
                    }
                }
            }
        }
    }

    return binaryImageVector;
}

std::vector<int> cannyConversion(const cv::Mat& img, AOIProperties mAOI)
{
    int AOIArea = mAOI.wdth * mAOI.hght;
    uchar *ptr_img = img.data;
    std::vector<int> binaryImageVectorRaw(AOIArea);

    for (int i = 0; i < AOIArea; i++)
    {
        if (ptr_img[i] == 255)  { binaryImageVectorRaw[i] = 1; }
        else                    { binaryImageVectorRaw[i] = 0; }
    }

    return binaryImageVectorRaw;
}

std::vector<int> getEdgeIndices(const std::vector<int>& binaryImageVector, int tag)
{
    int AOIArea = binaryImageVector.size();
    std::vector<int> cannyEdgeIndices;
    for (int iEdgePoint = 0; iEdgePoint < AOIArea; iEdgePoint++)
    { if (binaryImageVector[iEdgePoint] == tag) { cannyEdgeIndices.push_back(iEdgePoint); }}
    return cannyEdgeIndices;
}

std::vector<int> findEdges(const detectionVariables& mDetectionVariables, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    double pupilRadius = mDetectionVariables.thresholdChangePosition + mDetectionVariables.predictedCircumference * (1 + mDetectionVariables.thresholdChangeCircumference) / (2 * M_PI);

    // Find a starting edge point using Starburst-like algorithm

    std::vector<int> startIndices;

    for (int m = 0; m < 8; m++)
    {
        bool STOP_SEARCH = false;
        bool EDGE_FOUND  = false;

        int x = round(mDetectionVariables.predictedXPosRelative + dX[m]);
        int y = round(mDetectionVariables.predictedYPosRelative + dY[m]);

        double R = 0;

        while (!STOP_SEARCH)
        {
            int dx = dX[m];
            int dy = dY[m];

            for (int n = 0; n < 2; n++)
            {
                x = x + dx * (1 - n);
                y = y + dy * (n);

                if (x < 0 || x >= mAOI.wdth || y < 0 || y >= mAOI.hght)
                {
                    STOP_SEARCH = true;
                    break;
                }

                int centreIndex = y * mAOI.wdth + x;

                int tag = cannyEdgeVector[centreIndex];

                if (tag == 1)
                {
                    EDGE_FOUND = true;
                    startIndices.push_back(centreIndex);
                }

                if (EDGE_FOUND && R > pupilRadius)
                {
                    STOP_SEARCH = true;
                    break;
                }
            }

            R += sqrt((double) (dx * dx + dy * dy));
        }
    }

    // remove duplicates

    std::sort(startIndices.begin(), startIndices.end());
    startIndices.erase(std::unique(startIndices.begin(), startIndices.end()), startIndices.end());

    return startIndices;
}

void calculateEdgeDirections(const std::vector<int>& edgeIndices, std::vector<double>& edgeXTangents, std::vector<double>& edgeYTangents, AOIProperties mAOI)
{
    // scanned neighbours

    std::vector<int> dZ(8);
    dZ[0] = -1;
    dZ[1] = -mAOI.wdth - 1;
    dZ[2] = -mAOI.wdth;
    dZ[3] = -mAOI.wdth + 1;
    dZ[4] =  1;
    dZ[5] =  mAOI.wdth + 1;
    dZ[6] =  mAOI.wdth;
    dZ[7] =  mAOI.wdth - 1;

    int edgeLength = edgeIndices.size();

    // Calculate directions of edge points

    std::vector<double> xOrientation = { -1.0, -sqrt(0.5),  0.0,  sqrt(0.5), 1.0, sqrt(0.5), 0.0, -sqrt(0.5)};
    std::vector<double> yOrientation = {  0.0, -sqrt(0.5), -1.0, -sqrt(0.5), 0.0, sqrt(0.5), 1.0,  sqrt(0.5)};

    for (int iEdgePoint = 0; iEdgePoint < edgeLength - 1; iEdgePoint++)
    {
        int centreIndex    = edgeIndices[iEdgePoint];
        int neighbourIndex = edgeIndices[iEdgePoint + 1];

        for (int m = 0; m < 8; m++)
        {
            int adjacentIndex = centreIndex + dZ[m];

            if (neighbourIndex == adjacentIndex)
            {
                edgeXTangents[iEdgePoint] = xOrientation[m];
                edgeYTangents[iEdgePoint] = yOrientation[m];
                break;
            }
        }
    }

    // add last indices and give them directions of second-to-last

    edgeXTangents[edgeLength - 1] = edgeXTangents[edgeLength - 2];
    edgeYTangents[edgeLength - 1] = edgeYTangents[edgeLength - 2];
}

inline int calculateDirection(double x, double y)
{
    double thresholdLow  = 0.4142; // ~ tan(  M_PI/8)
    double thresholdHigh = 2.4142; // ~ tan(3*M_PI/8)

    int dir = 0;

    if (x == 0)
    {
        if      (y > 0) { dir = 2; }
        else if (y < 0) { dir = 6; }
        else            { dir = 0; } // arbitrary choice
    }
    else
    {
        double ratio = std::abs(y / x);

        if (x > 0) // Right half
        {
            if (y >= 0) // Top-right quadrant
            {
                if      (ratio < thresholdLow)  { dir = 0; }
                else if (ratio > thresholdHigh) { dir = 6; }
                else                            { dir = 7; }
            }
            else // Bottom-right quadrant
            {
                if      (ratio < thresholdLow)  { dir = 0; }
                else if (ratio > thresholdHigh) { dir = 2; }
                else                            { dir = 1; }
            }
        }
        else // Left half
        {
            if (y >= 0) // Top-left quadrant
            {
                if      (ratio < thresholdLow)  { dir = 4; }
                else if (ratio > thresholdHigh) { dir = 6; }
                else                            { dir = 5; }
            }
            else // Bottom-left quadrant
            {
                if      (ratio < thresholdLow)  { dir = 4; }
                else if (ratio > thresholdHigh) { dir = 2; }
                else                            { dir = 3; }
            }
        }
    }

    return dir;
}

int connectEdges(const detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1};
    std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1};

    std::vector<int> edgePoints = {startIndex};
    int edgePointNew = {startIndex};

    for (int iEdgePoint = 0; iEdgePoint < mDetectionParameters.windowLengthEdge; iEdgePoint++) // move back through edge
    {
        int centreIndex = edgePointNew;
        int centreXPos  =  centreIndex % mAOI.wdth;
        int centreYPos  = (centreIndex - centreXPos) / mAOI.wdth;

        for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
        {
            int neighbourXPos = centreXPos + dX[m];
            int neighbourYPos = centreYPos + dY[m];

            if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

            int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

            if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point is filled ...
            {
                std::vector<int>::iterator itr = find(edgePoints.begin(), edgePoints.end(), neighbourIndex); // check if index has already been stored
                if (itr == edgePoints.end())
                {
                    edgePointNew = neighbourIndex; // edge point to-be-checked
                    edgePoints.push_back(edgePointNew);
                }
            }
        }
    }

    int edgeLength = edgePoints.size();

    if (edgeLength >= mDetectionParameters.windowLengthEdge)
    {
        std::reverse(edgePoints.begin(), edgePoints.end());

        std::vector<double> edgeXTangents(edgeLength);
        std::vector<double> edgeYTangents(edgeLength);

        calculateEdgeDirections(edgePoints, edgeXTangents, edgeYTangents, mAOI);

        double xTangent = calculateMean(edgeXTangents);
        double yTangent = calculateMean(edgeYTangents);

        int dir_1 = calculateDirection(xTangent, yTangent);
        int dir_2 = 2 * dir_1;

        // Check if edge can be connected
        \
        std::vector<int> dX2 = {  2,  2,  2,  1,  0, -1, -2, -2, -2, -2, -2, -1,  0,  1,  2,  2 };
        std::vector<int> dY2 = {  0,  1,  2,  2,  2,  2,  2,  1,  0, -1, -2, -2, -2, -2, -2, -1 };

        int centreXPos  =  startIndex % mAOI.wdth;
        int centreYPos  = (startIndex - centreXPos) / mAOI.wdth;

        int neighbourXPos  = centreXPos + dX[dir_1];
        int neighbourYPos  = centreYPos + dY[dir_1];
        int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

        if (cannyEdgeVector[neighbourIndex] == 0)
        {
            for (int dR = -1; dR <= 1; dR++)
            {
                int k = (dir_2 + dR) % 16;
                if (k < 0) { k = k + 16; }
                int edgePointXPos = centreXPos + dX2[k];
                int edgePointYPos = centreYPos + dY2[k];
                if (edgePointXPos < 0 || edgePointXPos >= mAOI.wdth || edgePointYPos < 0 || edgePointYPos >= mAOI.hght) { continue; }
                int edgePointIndex = mAOI.wdth * edgePointYPos + edgePointXPos;

                int pointValue = cannyEdgeVector[edgePointIndex];

                if (pointValue == 1 || pointValue == 2) { return neighbourIndex; } // make the connection
            }
        }
    }

    return startIndex;
}

int findEdgePoints(const detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    cannyEdgeVector[startIndex] = 2; // tag pixel

    std::vector<int> edgePointsOld = {startIndex};

    int centreIndex   = 0;
    int numEdgePoints = 0;

    do
    {
        std::vector<int> edgePointsNew;

        numEdgePoints = edgePointsOld.size();

        for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++) // loop through all newly added unchecked edge points
        {
            centreIndex = edgePointsOld[iEdgePoint]; // index of current edge point

            int centreXPos = centreIndex % mAOI.wdth;
            int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;

            int nConnections = 0;

            for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];

                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                int neighbourTag = cannyEdgeVector[neighbourIndex];

                if (neighbourTag == 1) // if neighbouring point is filled ...
                {
                    cannyEdgeVector[neighbourIndex] = 2; // ... then tag it
                    edgePointsNew.push_back(neighbourIndex); // edge points to-be-checked
                    nConnections++;
                }
                else if (neighbourTag == 2) { nConnections++; }
            }

            if (nConnections == 1) // start or end of edge
            {
                int edgePointNew = connectEdges(mDetectionParameters, cannyEdgeVector, mAOI, centreIndex); // connect possible edge terminals

                if (centreIndex != edgePointNew)
                {
                    cannyEdgeVector[edgePointNew] = 2; // tag newly added point
                    edgePointsNew.push_back(edgePointNew);
                }
            }
        }

        edgePointsOld = edgePointsNew;
        numEdgePoints  = edgePointsOld.size();
        edgePointsNew.clear();

    } while (numEdgePoints > 0);

    return centreIndex; // return last point
}

std::vector<vertexProperties> findEdgeVertices(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    cannyEdgeVector[startIndex] = 3; // tag pixel

    std::vector<int> edgePointsOld = {startIndex};

    int numEdgePoints = 0;

    std::vector<vertexProperties> verticesAll;

    do
    {
        std::vector<int> edgePointsNew;

        numEdgePoints = edgePointsOld.size();

        for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++) // loop through all newly added unchecked edge points
        {
            int centreIndex = edgePointsOld[iEdgePoint]; // index of current edge point

            int centreXPos = centreIndex % mAOI.wdth;
            int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;

            std::vector<int> connections(8, 0);
            std::vector<int> edgePointsAll;
            int nConnections = 0;

            for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];

                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                int neighbourTag = cannyEdgeVector[neighbourIndex];

                if (neighbourTag == 2 || neighbourTag == 3 || neighbourTag == 4)
                {
                    edgePointsAll.push_back(neighbourIndex);
                    connections[m] = 1;
                    nConnections++;

                    if (neighbourTag == 2) // if neighbouring point is filled ...
                    {
                        cannyEdgeVector[neighbourIndex] = 3; // ... then tag it
                        edgePointsNew.push_back(neighbourIndex); // edge points to-be-checked
                    }
                }
            }

            // Check if a vertex was found

            bool VERTEX_FOUND = false;

            if      (nConnections  > 2 || nConnections == 1){ VERTEX_FOUND = true; }  // branch or terminal vertex
            else if (nConnections == 2) // possible vertex terminal
            {
                for (int m = 0; m < 8; m++)
                {
                    int n = m - 1;
                    if (n < 0) { n = n + 8; }
                    if (connections[m] == connections[n])
                    {
                        VERTEX_FOUND = true;
                        break;
                    }
                }
            }

            if (VERTEX_FOUND)
            {
                vertexProperties vertexNew;
                vertexNew.pointIndex      = centreIndex;
                vertexNew.connectedPoints = edgePointsAll;
                verticesAll.push_back(vertexNew);
                cannyEdgeVector[centreIndex] = 5;
            }
        }

        edgePointsOld = edgePointsNew;
        numEdgePoints  = edgePointsOld.size();
        edgePointsNew.clear();

    } while (numEdgePoints > 0);

    return verticesAll; // return last point
}

std::vector<branchProperties> findEdgeBranches(std::vector<vertexProperties>& verticesAll, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    std::vector<branchProperties> branchesAll;

    int numVertices = verticesAll.size();

    std::vector<int> vertexIndices(numVertices);

    for (int iVertex = 0; iVertex < numVertices; iVertex++)
    {
        vertexIndices[iVertex] = verticesAll[iVertex].pointIndex;
    }

    std::vector<std::vector<int>> vertexNeighbours(numVertices);
    std::vector<int> vertexNeighboursMatrix(numVertices * numVertices, 0);

    for (int iVertex = 0; iVertex < numVertices; iVertex++)
    {
        std::vector<int> connectedPoints = verticesAll[iVertex].connectedPoints;
        int numBranches = connectedPoints.size();

        for (int iBranch = 0; iBranch < numBranches; iBranch++)
        {
            int centreIndex = connectedPoints[iBranch];

            if (centreIndex == 5) // Check if vertex edge point neighbours another vertex edge point
            {
                int jVertex = std::distance(vertexIndices.begin(), find(vertexIndices.begin(), vertexIndices.end(), centreIndex));
                vertexNeighbours[iVertex].push_back(jVertex);
                vertexNeighboursMatrix[iVertex * numVertices + jVertex] = 1;
                continue;
            }

            // Create new branch and connect it with vertex

            branchProperties branchNew;
            branchNew.connectedVertices.resize(2); // each branch is connected to two vertices
            branchNew.connectedVertices[0] = iVertex;

            // Scan through branch until new vertex is found

            bool VERTEX_FOUND = false;

            while (!VERTEX_FOUND)
            {
                int centreXPos = centreIndex % mAOI.wdth;
                int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;

                for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                {
                    int neighbourXPos = centreXPos + dX[m];
                    int neighbourYPos = centreYPos + dY[m];

                    if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                    int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                    int neighbourTag = cannyEdgeVector[neighbourIndex];

                    if (neighbourTag == 5) // check if vertex has been found
                    {
                        int jVertex = std::distance(vertexIndices.begin(), find(vertexIndices.begin(), vertexIndices.end(), centreIndex));
                        branchNew.connectedVertices[1] = jVertex;
                        verticesAll[iVertex].connectedVertices.push_back(jVertex);
                        verticesAll[jVertex].connectedVertices.push_back(iVertex);
                        VERTEX_FOUND = true;
                        break;
                    }
                    else if (neighbourTag == 3) // branch continues
                    {
                        branchNew.pointIndices.push_back(neighbourIndex);
                        cannyEdgeVector[neighbourIndex] = 4;
                        centreIndex = neighbourIndex;
                        break;
                    }
                }
            }

            branchesAll.push_back(branchNew);
        }
    }

    // Adjacent vertices should all be connected to single hub vertex

    for (int iVertex = 0; iVertex < numVertices - 2; iVertex++) // redundent to check last two vertices
    {
        // Find all adjacent vertices to current vertex

        std::vector<int> vertexNeighboursNew = {iVertex}; // include self

        for (int jVertex = iVertex + 1; jVertex < numVertices; jVertex++) // undirected graph, so only check one half of matrix
        {
            int kVertex = iVertex * numVertices + jVertex;

            if (vertexNeighboursMatrix[kVertex] == 1)
            {
                vertexNeighboursNew.push_back(jVertex);
            }
        }

        // Find vertex that has greatest number of connected vertices

        int numVerticesNew = vertexNeighboursNew.size();

        if (numVerticesNew > 2)
        {
            std::vector<int> numNeighbours(numVerticesNew);

            for (int jVertex = 0; jVertex < numVerticesNew; jVertex++)
            {
                numNeighbours[jVertex] = vertexNeighbours[jVertex].size();
            }

            vertexNeighboursNew.erase(std::max_element(numNeighbours.begin(), numNeighbours.end())); // ignore hub vertex

            // Remove interconnections from connection matrix

            for (int jVertex = 0; jVertex < numVerticesNew - 1; jVertex++)
            {
                int xVertex = vertexNeighboursNew[jVertex];

                for (int kVertex = 0; kVertex < numVerticesNew - 1; kVertex++)
                {
                    int yVertex = vertexNeighboursNew[jVertex];

                    int vertexIndex = numVertices * yVertex + xVertex;

                    vertexNeighboursMatrix[vertexIndex] == 0;
                }
            }
        }
    }

    // Connect vertices with hub vertices

    for (int iVertex = 0; iVertex < numVertices - 1; iVertex++)
    {
        for (int jVertex = iVertex + 1; jVertex < numVertices; jVertex++)
        {
            int kVertex = iVertex * numVertices + jVertex;

            if (vertexNeighboursMatrix[kVertex] == 1) // create branch
            {
                branchProperties branchNew;
                branchNew.vertexNeighboursNew.resize(2);
                branchNew.vertexNeighboursNew[0] = iVertex;
                branchNew.vertexNeighboursNew[1] = jVertex;
                branchesAll.push_back(branchNew);
            }
        }
    }

    return branchesAll;
}


//void constructGraphTree(std::vector<vertexProperties>& vVertexProperties, std::vector<branchProperties>& vBranchProperties, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int pointIndexStart)
//{
//    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
//    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

//    vertexProperties vertexStart; // first vertex
//    vertexStart.pointIndex = pointIndexStart;
//    vertexStart.index      = 0;

//    { // Find all branches connected to start vertex
//        int centreXPos =  vertexStart.pointIndex % mAOI.wdth;
//        int centreYPos = (vertexStart.pointIndex - centreXPos) / mAOI.wdth;

//        for (int m = 0; m < 8; m++) // loop through 8-connected environment of the vertex
//        {
//            int neighbourXPos = centreXPos + dX[m];
//            int neighbourYPos = centreYPos + dY[m];

//            if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

//            int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

//            if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point is filled ...
//            {
//                cannyEdgeVector[neighbourIndex] = 3; // ... then tag it
//                vertexStart.connectedPoints.push_back(neighbourIndex);
//            }
//        }
//    }

//    std::vector<vertexProperties> verticesOld = {vertexStart}; // store first vertex

//    int branchNumber = 0; // counters
//    int vertexNumber = 0;

//    int numVertices;

//    do
//    {
//        std::vector<vertexProperties> verticesNew; // store new found vertices

//        numVertices = verticesOld.size();

//        for (int iVertex = 0; iVertex < numVertices; iVertex++)
//        {
//            vertexProperties vertexCurrent = verticesOld[iVertex]; // run through all found vertices
//            cannyEdgeVector[vertexCurrent.pointIndex] = 3; // tag vertex

//            // Loop through all found branches (if any) to find all edge points belonging to each branch

//            int numEdges = vertexCurrent.connectedPoints.size();

//            for (int iBranch = 0; iBranch < numEdges || branchNumber == 0; iBranch++)
//            {
//                branchProperties branchNew;
//                if (branchNumber == 0) { branchNew.pointIndices.push_back(vertexCurrent.pointIndex); } // first vertex should be included in first branch

//                if (numEdges > 0)
//                {
//                    int edgePointIndexOld = vertexCurrent.connectedPoints[iBranch];
//                    branchNew.pointIndices.push_back(edgePointIndexOld);

//                    bool VERTEX_FOUND = false;

//                    while(!VERTEX_FOUND) // keep going until vertex is encountered
//                    {
//                        std::vector<int> edgePointsNew;

//                        int centreIndex = edgePointIndexOld;
//                        int centreXPos  = centreIndex % mAOI.wdth;
//                        int centreYPos  = (centreIndex - centreXPos) / mAOI.wdth;

//                        for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
//                        {
//                            int neighbourXPos = centreXPos + dX[m];
//                            int neighbourYPos = centreYPos + dY[m];

//                            if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

//                            int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

//                            if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point was tagged previously ...
//                            {
//                                cannyEdgeVector[neighbourIndex] = 3; // ... give it a new tag
//                                edgePointsNew.push_back(neighbourIndex);
//                            }
//                        }

//                        int numEdgePoints = edgePointsNew.size();

//                        if (numEdgePoints == 1) // edge continues
//                        {
//                            edgePointIndexOld = edgePointsNew[0]; // update index
//                            branchNew.pointIndices.push_back(edgePointIndexOld); // check new index
//                            edgePointsNew.clear();
//                        }
//                        else // new vertex found
//                        {
//                            vertexProperties vertexNew;
//                            vertexNew.pointIndex = centreIndex;
//                            vertexNew.connectedBranches.push_back(branchNumber); // add current branch to new vertex connections
//                            if (edgePointsNew.size() > 0) { vertexNew.connectedPoints = edgePointsNew; }
//                            verticesNew.push_back(vertexNew); // new vertices to be checked
//                            VERTEX_FOUND = true;
//                        }
//                    }
//                }

//                vertexCurrent.connectedBranches.push_back(branchNumber); // add new branches

//                int branchLength = branchNew.pointIndices.size();

//                if (branchLength > 0)
//                {
//                    branchNew.length = branchLength;
//                    branchNew.index  = branchNumber;
//                    vBranchProperties.push_back(branchNew); // record all branches
//                }

//                branchNumber++;
//            }

//            vertexCurrent.index = vertexNumber;
//            vVertexProperties.push_back(vertexCurrent); // record all vertices
//            vertexNumber++;
//        }

//        numVertices = verticesNew.size();
//        verticesOld = verticesNew;

//    } while (numVertices > 0);

//    // Find vertices that branches are connected to

//    numVertices = vVertexProperties.size();

//    for (int iVertex = 0; iVertex < numVertices; iVertex++)
//    {
//        std::vector<int> vertexConnections = vVertexProperties[iVertex].connectedBranches;

//        for (int iEdge = 0, numConnections = vertexConnections.size(); iEdge < numConnections; iEdge++)
//        {
//            vBranchProperties[vertexConnections[iEdge]].connectedVertices.push_back(iVertex);
//        }
//    }
//}

std::vector<int> maxPath(const std::vector<std::vector<int>>& pathAll, const std::vector<branchProperties>& vBranchProperties)
{
    int numPaths = pathAll.size();
    std::vector<int> pathLengths(numPaths);

    for (int iPath = 0; iPath < numPaths; iPath++)
    {
        int lengthTotal = 0;
        std::vector<int> pathCurrent = pathAll[iPath];
        int numPathBranches = pathCurrent.size();

        for (int iBranch = 0; iBranch < numPathBranches; iBranch++)
        {
            lengthTotal += vBranchProperties[pathCurrent[iBranch]].length;
        }

        pathLengths[iPath] = lengthTotal;
    }

    int pathIndexMax = std::distance(pathLengths.begin(), std::max_element(pathLengths.begin(), pathLengths.end()));
    return pathAll[pathIndexMax];
}

void findPath(std::vector<vertexProperties>& vVertexPropertiesAll, std::vector<std::vector<int>>& pathBranchesAll, std::vector<int>& pathBranchesRoot, std::vector<int>& branchesChecked, int vertexIndex)
{
    vertexProperties mVertex = vVertexPropertiesAll[vertexIndex];
    int numBranches = mVertex.connectedBranches.size();

    for (int iBranch = 0; iBranch < numBranches; iBranch++)
    {
        int branchIndexNew = mVertex.connectedBranches[iBranch];

        std::vector<int> pathBranchesRootNew = pathBranchesRoot;
        pathBranchesRootNew.push_back(branchIndexNew);

        std::vector<int> branchesCheckedNew = branchesChecked;
        branchesCheckedNew[branchIndexNew]  = 1;

        std::vector<std::vector<int>> pathBranchesNew = findPath(pathBranchesAll, pathBranchesRootNew, vVertexPropertiesAll, vertexIndexNew, branchesCheckedNew);

        pathBranchesAll.insert(std::end(pathBranchesAll), std::begin(pathBranchesNew), std::end(pathBranchesNew));
    }
}

void findLongestPath(std::vector<vertexProperties>& vVertexPropertiesAll, std::vector<branchProperties>& vBranchPropertiesAll, std::vector<int>& allPoints, std::vector<int>& pathPoints)
{
    int numBranchesAll = vBranchPropertiesAll.size();
    int numVerticesAll = vVertexPropertiesAll.size();

    // Count number of branches

    // Perform depth-first search for each vertex

    for (int iVertex = 0; iVertex < numVerticesAll; iVertex++) // loop through all vertices
    {
        std::vector<std::vector<int>> pathBranchesAll;
        std::vector<int> pathBranchesRoot;
        std::vector<int> branchesChecked(numBranchesAll, 0);

        std::vector<std::vector<int>> paths = findPath(vVertexPropertiesAll, pathBranchesAll, pathBranchesRoot, branchesChecked, iVertex);

        vertexProperties vertexStart = vVertexPropertiesAll[iVertex]; // starting vertex
        int numBranches = vertexStart.connectedBranches.size();

        for (int iBranch = 0; iBranch < numBranches; iBranch++)
        {
            int vertexNewIndex = vertexStart.connectedVertices[iBranch];

            std::vector<int> pathRootNew = pathRoot;
            pathRootNew.push_back(vertexNewIndex);

            vertexProperties vertexNew = vVertexPropertiesAll[vertexNewIndex]; // starting vertex
            int numBranchesNew = vertexNew.connectedBranches.size();

            for (int jBranch = 0; jBranch < numBranchesNew; jBranch++)
            {
                int vertexNewNewIndex = vertexNew.connectedVertices[iBranch];

                pathRoot.push_back(vertexNewNewIndex);

                vertexProperties vertexNewNew = vVertexPropertiesAll[vertexNewIndex]; // starting vertex
                int numBranchesNewNew = vertexNew.connectedBranches.size();

                {  // etc


                } //

                paths.push_back(pathRootNewNew);
            }

             paths.push_back(pathRootNew);
        }

        paths.push_back(pathRoot);
    }















        std::vector<std::vector<int>> pathAll(numVerticesAll);

        for (int iVertex = 0; iVertex < numVerticesAll; iVertex++) // loop through all vertices
        {
            vertexProperties vertexStart = vVertexPropertiesAll[iVertex]; // starting vertex

            std::vector<std::vector<int>> pathAllCurrent; // record all paths for current starting vertex

            std::vector<bool> branchFlagsAll(numBranchesAll, true);  // record if branches should be included again. Don't include terminals more than once

            while(true) // check all paths for current starting vertex
            {
                vertexProperties vertexNew = vertexStart;

                std::vector<bool> branchFlags(numBranchesAll, true); // no branch should be included twice

                std::vector<int> pathNew; // new path for current vertex

                bool PATH_CONTINUES = true;

                while(PATH_CONTINUES) // while path continues
                {
                    int numConnections = vertexNew.connectedBranches.size();

                    int branchNumber;
                    PATH_CONTINUES = false;

                    for (int iEdge = 0; iEdge < numConnections; iEdge++) // find first new branch that should still be included in path
                    {
                        branchNumber = vertexNew.connectedBranches[iEdge];

                        if (branchFlagsAll[branchNumber] && branchFlags[branchNumber]) // branch should still be included
                        {
                            branchFlags[branchNumber] = false;
                            pathNew.push_back(branchNumber); // record branch indices for current path
                            PATH_CONTINUES = true;
                            break;
                        }
                    }

                    if (PATH_CONTINUES) // Get new vertex
                    {
                        std::vector<int> connectedVertices = vBranchPropertiesAll[branchNumber].connectedVertices; // get the vertices the branch is connected to
                        for (int jVertex = 0, vSize = connectedVertices.size(); jVertex < vSize; jVertex++) // find next vertex
                        {
                            int vertexIndex = connectedVertices[jVertex];

                            if (vertexIndex != vertexNew.index) // don't include current vertex
                            {
                                vertexNew = vVertexPropertiesAll[vertexIndex];
                                break;
                            }
                        }
                    }
                }

                if (pathNew.size() > 0)
                {
                    branchFlagsAll[pathNew.back()] = false; // don't include last included branch again. We assume no alternative paths
                    pathAllCurrent.push_back(pathNew);
                }
                else { break; } // no new path found
            }

            pathAll[iVertex] = maxPath(pathAllCurrent, vBranchPropertiesAll); // get longest path for current starting vertex
        }



























        // Find all paths in edge collection

        std::vector<std::vector<int>> pathAll(numVerticesAll);

        for (int iVertex = 0; iVertex < numVerticesAll; iVertex++) // loop through all vertices
        {
            vertexProperties vertexStart = vVertexPropertiesAll[iVertex]; // starting vertex

            //        if (vertexStart.connectedBranches.size() > 1) { continue; } // only start with terminal vertices

            std::vector<std::vector<int>> pathAllCurrent; // record all paths for current starting vertex

            std::vector<bool> branchFlagsAll(numBranchesAll, true);  // record if branches should be included again. Don't include terminals more than once

            while(true) // check all paths for current starting vertex
            {
                vertexProperties vertexNew = vertexStart;

                std::vector<bool> branchFlags(numBranchesAll, true); // no branch should be included twice

                std::vector<int> pathNew; // new path for current vertex

                bool PATH_CONTINUES = true;

                while(PATH_CONTINUES) // while path continues
                {
                    int numConnections = vertexNew.connectedBranches.size();

                    int branchNumber;
                    PATH_CONTINUES = false;

                    for (int iEdge = 0; iEdge < numConnections; iEdge++) // find first new branch that should still be included in path
                    {
                        branchNumber = vertexNew.connectedBranches[iEdge];

                        if (branchFlagsAll[branchNumber] && branchFlags[branchNumber]) // branch should still be included
                        {
                            branchFlags[branchNumber] = false;
                            pathNew.push_back(branchNumber); // record branch indices for current path
                            PATH_CONTINUES = true;
                            break;
                        }
                    }

                    if (PATH_CONTINUES) // Get new vertex
                    {
                        std::vector<int> connectedVertices = vBranchPropertiesAll[branchNumber].connectedVertices; // get the vertices the branch is connected to
                        for (int jVertex = 0, vSize = connectedVertices.size(); jVertex < vSize; jVertex++) // find next vertex
                        {
                            int vertexIndex = connectedVertices[jVertex];

                            if (vertexIndex != vertexNew.index) // don't include current vertex
                            {
                                vertexNew = vVertexPropertiesAll[vertexIndex];
                                break;
                            }
                        }
                    }
                }

                if (pathNew.size() > 0)
                {
                    branchFlagsAll[pathNew.back()] = false; // don't include last included branch again. We assume no alternative paths
                    pathAllCurrent.push_back(pathNew);
                }
                else { break; } // no new path found
            }

            pathAll[iVertex] = maxPath(pathAllCurrent, vBranchPropertiesAll); // get longest path for current starting vertex
        }

        std::vector<int> pathBranchIndices = maxPath(pathAll, vBranchPropertiesAll); // get longest path overall
        int numBranchesPath = pathBranchIndices.size();

        // Grab branches for path

        std::vector<branchProperties> vBranchPropertiesPath(numBranchesPath);

        for (int iBranch = 0; iBranch < numBranchesPath; iBranch++)
        {
            vBranchPropertiesPath[iBranch] = vBranchPropertiesAll[pathBranchIndices[iBranch]];
        }

        // Reverse branches that are not properly aligned in vector

        for (int iBranch = 0; iBranch < numBranchesPath - 1; iBranch++)
        {
            branchProperties mBranchProperties_1 = vBranchPropertiesPath[iBranch];
            branchProperties mBranchProperties_2 = vBranchPropertiesPath[iBranch + 1];

            int numConnections_1 = mBranchProperties_1.connectedVertices.size();
            if (numConnections_1 == 1) { continue; } // alignment is probably fine

            bool BRANCH_CONNECTED = false;

            for (int iVertex = 0; iVertex < numConnections_1 && !BRANCH_CONNECTED; iVertex++)
            {
                int vertex_1 = mBranchProperties_1.connectedVertices[iVertex];
                int numConnections_2 = mBranchProperties_2.connectedVertices.size();

                for (int jVertex = 0; jVertex < numConnections_2; jVertex++)
                {
                    int vertex_2 = mBranchProperties_2.connectedVertices[jVertex];

                    if (vertex_1 == vertex_2) // last vertex of current branch should be equal to first vertex of next branch
                    {
                        if (iVertex == 0) { std::reverse(vBranchPropertiesPath[iBranch].pointIndices.begin(), vBranchPropertiesPath[iBranch].pointIndices.end()); }
                        BRANCH_CONNECTED = true; // successfully connected
                        break;
                    }
                }
            }
        }

        // All path points

        for (int iBranch = 0; iBranch < numBranchesPath; iBranch++)
        {
            branchProperties mBranchProperties = vBranchPropertiesPath[iBranch];
            pathPoints.insert(std::end(pathPoints), std::begin(mBranchProperties.pointIndices), std::end(mBranchProperties.pointIndices));
        }

        // All edge points

        for (int iBranch = 0; iBranch < numBranchesAll; iBranch++)
        {
            allPoints.insert(std::end(allPoints), std::begin(vBranchPropertiesAll[iBranch].pointIndices), std::end(vBranchPropertiesAll[iBranch].pointIndices));
        }
    }

    std::vector<edgeProperties> edgeSelection(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
    {
        std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all edges

        int numEdges = 0;

        do
        {
            std::vector<int> startIndicesRaw = findEdges(mDetectionVariables, cannyEdgeVector, mAOI);

            numEdges = startIndicesRaw.size();
            std::vector<int> startIndices(numEdges);

            for (int iEdge = 0; iEdge < numEdges; iEdge++)
            {
                startIndices[iEdge] = findEdgePoints(mDetectionParameters, cannyEdgeVector, mAOI, startIndicesRaw[iEdge]); // tag all edge points of found edges
            }

            for (int iEdge = 0; iEdge < numEdges; iEdge++)
            {
                // Find all vertices and all connected branches (i.e. obtain graph tree)



                // Find longest path

                std::vector<int> pathIndices;
                std::vector<int>  allIndices;

                findLongestPath(allIndices, pathIndices, cannyEdgeVector, mAOI, startIndices[iEdge]);

                // Give points in longest path a new tag

                for (int iEdgePoint = 0, edgeSize = pathIndices.size(); iEdgePoint < edgeSize; iEdgePoint++)
                { cannyEdgeVector[pathIndices[iEdgePoint]] = 4; }

                // Remove tag from points that have been tagged before, but not included in path

                for (int iEdgePoint = 0, edgeLength = allIndices.size(); iEdgePoint < edgeLength; iEdgePoint++)
                {
                    int edgePointIndex  = allIndices[iEdgePoint];
                    int edgePointTag    = cannyEdgeVector[edgePointIndex];
                    if (edgePointTag == 2 || edgePointTag == 3)
                    {
                        cannyEdgeVector[edgePointIndex] = 1;
                        std::vector<int>::iterator itr = find(startIndices.begin(), startIndices.end(), edgePointIndex); // check if index has already been stored
                        if (itr != startIndices.end()) { startIndices.erase(itr); }
                    }
                }

                edgeProperties mEdgeProperties;
                mEdgeProperties.pointIndices = pathIndices;

                vEdgePropertiesAll.push_back(mEdgeProperties);
            }
        } while (numEdges > 1);

        return vEdgePropertiesAll;
    }

    void calculateCurvatureLimits(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, double& curvatureUpperLimit, double& curvatureLowerLimit)
    {
        // Calculate curvature limits

        std::vector<double> circumferences(2);
        std::vector<double>   aspectRatios(2);

        double certaintyFactorFeatures = 0.5 * (mDetectionVariables.certaintyFeatures + 1);
        double certaintyFactorAverages = 0.5 * (mDetectionVariables.certaintyAverages + 1);

        if (certaintyFactorAverages > certaintyFactorFeatures + certaintyOffset)
        {
            circumferences[0] = mDetectionVariables.averageCircumference + (mDetectionVariables.averageCircumference * mDetectionVariables.offsetCircumference);
            circumferences[1] = mDetectionVariables.averageCircumference - (mDetectionVariables.averageCircumference * mDetectionVariables.offsetCircumference);
        }
        else
        {
            circumferences[0] = mDetectionVariables.predictedCircumference + (mDetectionVariables.predictedCircumference * mDetectionVariables.thresholdChangeCircumference);
            circumferences[1] = mDetectionVariables.predictedCircumference - (mDetectionVariables.predictedCircumference * mDetectionVariables.thresholdChangeCircumference);
        }

        aspectRatios[0] = mDetectionVariables.predictedAspectRatio + mDetectionVariables.thresholdChangeAspectRatio;
        aspectRatios[1] = mDetectionVariables.predictedAspectRatio - mDetectionVariables.thresholdChangeAspectRatio;

        // Calculate limits

        std::vector<double> curvaturesMax(4);
        std::vector<double> curvaturesMin(4);

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                curvaturesMax[2 * j + i] = getCurvatureUpperLimit(circumferences[i], aspectRatios[j], mDetectionParameters.windowLengthEdge);
                curvaturesMin[2 * j + i] = getCurvatureLowerLimit(circumferences[i], aspectRatios[j], mDetectionParameters.windowLengthEdge);
            }
        }

        curvatureUpperLimit = *std::max_element(std::begin(curvaturesMax), std::end(curvaturesMax)) + mDetectionParameters.curvatureOffset;
        curvatureLowerLimit = *std::min_element(std::begin(curvaturesMin), std::end(curvaturesMin)) - mDetectionParameters.curvatureOffset;
    }

    std::vector<double> calculateCurvatures(const detectionParameters& mDetectionParameters, std::vector<double>& xNormals, std::vector<double>& yNormals, const std::vector<double>& xTangentsAll, const std::vector<double>& yTangentsAll)
    {
        int edgeSize = xTangentsAll.size();

        std::vector<double> curvatures(edgeSize, 0.0);

        for (int iEdgePoint = mDetectionParameters.windowLengthEdge; iEdgePoint < edgeSize - mDetectionParameters.windowLengthEdge; iEdgePoint++)
        {
            // calculate window tangents

            // first window

            std::vector<double> tangentsX_1(xTangentsAll.begin() + iEdgePoint - mDetectionParameters.windowLengthEdge, xTangentsAll.begin() + iEdgePoint);
            std::vector<double> tangentsY_1(yTangentsAll.begin() + iEdgePoint - mDetectionParameters.windowLengthEdge, yTangentsAll.begin() + iEdgePoint);

            double tangentXMean_1 = calculateMean(tangentsX_1);
            double tangentYMean_1 = calculateMean(tangentsY_1);

            // second window

            std::vector<double> tangentsX_2(xTangentsAll.begin() + iEdgePoint + 1, xTangentsAll.begin() + iEdgePoint + 1 + mDetectionParameters.windowLengthEdge);
            std::vector<double> tangentsY_2(yTangentsAll.begin() + iEdgePoint + 1, yTangentsAll.begin() + iEdgePoint + 1 + mDetectionParameters.windowLengthEdge);

            double tangentXMean_2 = calculateMean(tangentsX_2);
            double tangentYMean_2 = calculateMean(tangentsY_2);

            // calculate vector difference

            double vectorAngle = atan2(tangentYMean_2, tangentXMean_2) - atan2(tangentYMean_1, tangentXMean_1);

            if      (vectorAngle >  M_PI) { vectorAngle = vectorAngle - 2 * M_PI; }
            else if (vectorAngle < -M_PI) { vectorAngle = vectorAngle + 2 * M_PI; }

            curvatures[iEdgePoint] = 180 * vectorAngle / M_PI; // in degrees

            xNormals[iEdgePoint] = tangentXMean_2 - tangentXMean_1;
            yNormals[iEdgePoint] = tangentYMean_2 - tangentYMean_1;
        }

        return curvatures;
    }

    std::vector<edgeProperties> edgeSegmentationCurvature(const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties, const double curvatureUpperLimit, const double curvatureLowerLimit)
    {
        std::vector<edgeProperties> vEdgePropertiesAll;
        std::vector<edgeProperties> vEdgePropertiesNew;

        edgeProperties mEdgePropertiesNew = mEdgeProperties;

        do
        {
            int edgeSize = mEdgePropertiesNew.curvatures.size();

            // check majority sign

            int curvatureSign = 1;

            int numPos = 0;
            int numNeg = 0;

            for (int iEdgePoint = mDetectionParameters.windowLengthEdge; iEdgePoint < edgeSize - mDetectionParameters.windowLengthEdge; iEdgePoint++)
            {
                double curvature = mEdgePropertiesNew.curvatures[iEdgePoint];
                if      (curvature > 0) { numPos++; }
                else if (curvature < 0) { numNeg++; }
            }

            if (numNeg > numPos) {curvatureSign = -1;} // if majority sign is negative, then swap sign of threshold

            // find breakpoints based on curvature thresholding

            std::vector<int> breakPoints; // position of breakpoints
            breakPoints.push_back(-1); // add first point (+ 1 is added later)

            for (int iEdgePoint = mDetectionParameters.windowLengthEdge; iEdgePoint < edgeSize - mDetectionParameters.windowLengthEdge; iEdgePoint++)
            {
                double curvature = mEdgePropertiesNew.curvatures[iEdgePoint];

                if (std::abs(curvature) >= curvatureUpperLimit || curvatureSign * curvature <= curvatureLowerLimit)
                {
                    breakPoints.push_back(iEdgePoint);
                    break;
                }
            }

            // always add last point

            breakPoints.push_back(edgeSize - 1);

            // make sure length between breakpoints is long enough

            int numBreakPoints = breakPoints.size();

            std::vector<int> breakPointsNew; // position of breakpoints
            breakPointsNew.push_back(-1); // add first point (+ 1 is added later)

            for (int iBreakPoint = 1; iBreakPoint < numBreakPoints - 1; iBreakPoint++)
            {
                int iEdgePointLeft = breakPoints[iBreakPoint - 1];
                int iEdgePointCntr = breakPoints[iBreakPoint];
                int iEdgePointRght = breakPoints[iBreakPoint + 1];

                int diffLeft = iEdgePointCntr - iEdgePointLeft;
                int diffRght = iEdgePointRght - iEdgePointCntr;

                if (diffLeft > mDetectionParameters.windowLengthEdge || diffRght > mDetectionParameters.windowLengthEdge)
                {
                    breakPointsNew.push_back(iEdgePointCntr);
                }
            }

            breakPointsNew.push_back(edgeSize - 1); // add last point again

            // cut edge at breakpoints

            vEdgePropertiesNew.clear();

            int numBreakPointsNew = breakPointsNew.size();

            for (int iBreakPoint = 0; iBreakPoint < numBreakPointsNew - 1; iBreakPoint++)
            {
                int iStartBreakPoint = breakPointsNew[iBreakPoint] + 1;
                int iEndBreakPoint   = breakPointsNew[iBreakPoint  + 1];

                edgeProperties mEdgePropertiesTemp;
                mEdgePropertiesTemp.pointIndices = std::vector<int>    (mEdgePropertiesNew.pointIndices.begin() + iStartBreakPoint, mEdgePropertiesNew.pointIndices.begin() + iEndBreakPoint);
                mEdgePropertiesTemp.curvatures   = std::vector<double> (mEdgePropertiesNew.curvatures.begin()   + iStartBreakPoint, mEdgePropertiesNew.curvatures.begin()   + iEndBreakPoint);
                mEdgePropertiesTemp.xnormals     = std::vector<double> (mEdgePropertiesNew.xnormals.begin()     + iStartBreakPoint, mEdgePropertiesNew.xnormals.begin()     + iEndBreakPoint);
                mEdgePropertiesTemp.ynormals     = std::vector<double> (mEdgePropertiesNew.ynormals.begin()     + iStartBreakPoint, mEdgePropertiesNew.ynormals.begin()     + iEndBreakPoint);

                vEdgePropertiesNew.push_back(mEdgePropertiesTemp);
            }

            vEdgePropertiesAll.push_back(vEdgePropertiesNew[0]); // record first edge
            mEdgePropertiesNew = vEdgePropertiesNew.back();      // check last edge

        } while (vEdgePropertiesNew.size() > 1);

        return vEdgePropertiesAll;
    }

    void calculateCurvatureStats(const detectionParameters& mDetectionParameters, edgeProperties& mEdgeProperties)
    {
        // Calculate min, max and mean curvature

        int edgeSize = mEdgeProperties.curvatures.size();

        double curvatureAvg = 0;
        double curvatureMax = 0;
        double curvatureMin = std::numeric_limits<double>::max();

        if (edgeSize > 2 * mDetectionParameters.windowLengthEdge)
        {
            for (int iEdgePoint = mDetectionParameters.windowLengthEdge; iEdgePoint < edgeSize - mDetectionParameters.windowLengthEdge; iEdgePoint++)
            {
                double curvature = mEdgeProperties.curvatures[iEdgePoint];

                if (curvature < curvatureMin) { curvatureMin = curvature; }
                if (curvature > curvatureMax) { curvatureMax = curvature; }

                curvatureAvg += std::abs(curvature);
            }

            curvatureAvg = curvatureAvg / edgeSize;
        }
        else
        {
            curvatureAvg = 360;
            curvatureMax = 360;
            curvatureMin = 360;
        }

        mEdgeProperties.curvature    = curvatureAvg;
        mEdgeProperties.curvatureMax = curvatureMax;
        mEdgeProperties.curvatureMin = curvatureMin;
    }

    std::vector<edgeProperties> edgeSegmentationLength(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties)
    {
        // This functions cuts edge terminals to make the edge shorter, if the edge is significantly longer than predicted

        int edgeSize = mEdgeProperties.curvatures.size();

        // find breakpoints based on length thresholding

        std::vector<int> breakPoints; // position of breakpoints
        breakPoints.push_back(-1); // add first point (+ 1 is added later)

        double lengthDifference = mEdgeProperties.length - mDetectionVariables.predictedCircumference;

        if (lengthDifference > mDetectionParameters.fitEdgeFraction * mDetectionVariables.predictedCircumference) // difference should be large enough
        {
            if (edgeSize > 2 * lengthDifference + mDetectionParameters.windowLengthEdge) // should be enough space between breakpoints
            {
                breakPoints.push_back(lengthDifference);
                breakPoints.push_back(edgeSize - lengthDifference - 1);
            }
        }

        // add last point

        breakPoints.push_back(edgeSize - 1);

        // Do segmentation

        int numBreakPoints = breakPoints.size();
        std::vector<edgeProperties> vEdgeProperties(numBreakPoints - 1);

        if (numBreakPoints == 4)
        {
            // cut edge at breakpoints

            for (int iBreakPoint = 0; iBreakPoint < numBreakPoints - 1; iBreakPoint++)
            {
                int iStartBreakPoint = breakPoints[iBreakPoint] + 1;
                int iEndBreakPoint   = breakPoints[iBreakPoint  + 1];

                edgeProperties mEdgePropertiesNew;

                mEdgePropertiesNew.pointIndices = std::vector<int>   (mEdgeProperties.pointIndices.begin() + iStartBreakPoint, mEdgeProperties.pointIndices.begin() + iEndBreakPoint);
                mEdgePropertiesNew.intensities  = std::vector<int>   (mEdgeProperties.intensities.begin()  + iStartBreakPoint, mEdgeProperties.intensities.begin()  + iEndBreakPoint);
                mEdgePropertiesNew.gradients    = std::vector<int>   (mEdgeProperties.gradients.begin()    + iStartBreakPoint, mEdgeProperties.gradients.begin()    + iEndBreakPoint);
                mEdgePropertiesNew.radii        = std::vector<double>(mEdgeProperties.radii.begin()        + iStartBreakPoint, mEdgeProperties.radii.begin()        + iEndBreakPoint);
                mEdgePropertiesNew.curvatures   = std::vector<double>(mEdgeProperties.curvatures.begin()   + iStartBreakPoint, mEdgeProperties.curvatures.begin()   + iEndBreakPoint);

                mEdgePropertiesNew.intensity = calculateMeanInt(mEdgePropertiesNew.intensities);
                mEdgePropertiesNew.gradient  = calculateMeanInt(mEdgePropertiesNew.gradients);
                mEdgePropertiesNew.radius    = calculateMean(mEdgePropertiesNew.radii);
                mEdgePropertiesNew.length    = mEdgePropertiesNew.pointIndices.size();

                calculateCurvatureStats(mDetectionParameters, mEdgePropertiesNew);

                vEdgeProperties[iBreakPoint] = mEdgePropertiesNew;
            }

            // only remove one of the two edge terminals - re-attach other one
            // compare edge characteristics of terminals with middle section
            // 0 = start terminal
            // 1 = main section
            // 2 = end terminal

            std::vector<double> inputVector_1(7); // for start terminal
            inputVector_1[0] = vEdgeProperties[0].radius / vEdgeProperties[1].radius;
            inputVector_1[1] = 0; // don't use circumference
            inputVector_1[2] = vEdgeProperties[0].curvature - vEdgeProperties[1].curvature;
            inputVector_1[3] = vEdgeProperties[0].intensity - vEdgeProperties[1].intensity;
            inputVector_1[4] = vEdgeProperties[0].gradient  - vEdgeProperties[1].gradient;
            inputVector_1[5] = 0; // don't use radial variance
            inputVector_1[6] = vEdgeProperties[0].length;

            std::vector<double> inputVector_2(7); // for end terminal
            inputVector_2[0] = vEdgeProperties[2].radius / vEdgeProperties[1].radius;
            inputVector_2[1] = 0; // don't use circumference
            inputVector_2[2] = vEdgeProperties[2].curvature - vEdgeProperties[1].curvature;
            inputVector_2[3] = vEdgeProperties[2].intensity - vEdgeProperties[1].intensity;
            inputVector_2[4] = vEdgeProperties[2].gradient  - vEdgeProperties[1].gradient;
            inputVector_2[5] = 0; // don't use radial variance
            inputVector_2[6] = vEdgeProperties[2].length;

            double scoreTotal_1 = calculateScoreTotal(mDetectionVariables, inputVector_1, false);
            double scoreTotal_2 = calculateScoreTotal(mDetectionVariables, inputVector_2, false);

            int indexStart = 0;
            int indexEnd   = 0;

            if      (scoreTotal_1 > scoreTotal_2) { indexStart = 0; indexEnd = 1; }
            else if (scoreTotal_1 < scoreTotal_2) { indexStart = 1; indexEnd = 2; }

            if (indexStart != indexEnd) // if scores are equal, cut both terminals
            {
                // concatenate vectors

                edgeProperties mEdgePropertiesNew;

                mEdgePropertiesNew.pointIndices.reserve(vEdgeProperties[indexStart].pointIndices.size() + vEdgeProperties[indexEnd].pointIndices.size());
                mEdgePropertiesNew.pointIndices.insert(mEdgePropertiesNew.pointIndices.end(), vEdgeProperties[indexStart].pointIndices.begin(), vEdgeProperties[indexStart].pointIndices.end() );
                mEdgePropertiesNew.pointIndices.insert(mEdgePropertiesNew.pointIndices.end(), vEdgeProperties[indexEnd].pointIndices.begin(),   vEdgeProperties[indexEnd].pointIndices.end() );

                mEdgePropertiesNew.intensities.reserve(vEdgeProperties[indexStart].intensities.size() + vEdgeProperties[indexEnd].intensities.size());
                mEdgePropertiesNew.intensities.insert(mEdgePropertiesNew.intensities.end(), vEdgeProperties[indexStart].intensities.begin(), vEdgeProperties[indexStart].intensities.end() );
                mEdgePropertiesNew.intensities.insert(mEdgePropertiesNew.intensities.end(), vEdgeProperties[indexEnd].intensities.begin(),   vEdgeProperties[indexEnd].intensities.end() );

                mEdgePropertiesNew.gradients.reserve(vEdgeProperties[indexStart].gradients.size() + vEdgeProperties[indexEnd].gradients.size());
                mEdgePropertiesNew.gradients.insert(mEdgePropertiesNew.gradients.end(), vEdgeProperties[indexStart].gradients.begin(), vEdgeProperties[indexStart].gradients.end() );
                mEdgePropertiesNew.gradients.insert(mEdgePropertiesNew.gradients.end(), vEdgeProperties[indexEnd].gradients.begin(),   vEdgeProperties[indexEnd].gradients.end() );

                mEdgePropertiesNew.radii.reserve(vEdgeProperties[indexStart].radii.size() + vEdgeProperties[indexEnd].radii.size());
                mEdgePropertiesNew.radii.insert(mEdgePropertiesNew.radii.end(), vEdgeProperties[indexStart].radii.begin(), vEdgeProperties[indexStart].radii.end() );
                mEdgePropertiesNew.radii.insert(mEdgePropertiesNew.radii.end(), vEdgeProperties[indexEnd].radii.begin(),   vEdgeProperties[indexEnd].radii.end() );

                mEdgePropertiesNew.curvatures.reserve(vEdgeProperties[indexStart].curvatures.size() + vEdgeProperties[indexEnd].curvatures.size());
                mEdgePropertiesNew.curvatures.insert(mEdgePropertiesNew.curvatures.end(), vEdgeProperties[indexStart].curvatures.begin(), vEdgeProperties[indexStart].curvatures.end() );
                mEdgePropertiesNew.curvatures.insert(mEdgePropertiesNew.curvatures.end(), vEdgeProperties[indexEnd].curvatures.begin(),   vEdgeProperties[indexEnd].curvatures.end() );

                std::vector<edgeProperties> vEdgePropertiesNew(2);
                vEdgePropertiesNew[0] = mEdgePropertiesNew;
                vEdgePropertiesNew[1] = vEdgeProperties[(indexStart + 2) % 3];

                vEdgeProperties = vEdgePropertiesNew; // update return vector
            }
        }
        else // do nothing
        {
            vEdgeProperties[0] = mEdgeProperties;
        }

        return vEdgeProperties;
    }

    std::vector<int> calculateRadialGradients(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const cv::Mat& img, const std::vector<int> edgeIndices)
    {
        std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1 };
        std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1 };

        double pupilXCentre = mDetectionVariables.predictedXPosRelative;
        double pupilYCentre = mDetectionVariables.predictedYPosRelative;

        uchar *ptr = img.data;
        int width  = img.cols;
        int height = img.rows;

        int kernelRadius = (mDetectionParameters.cannyKernelSize - 1) / 2;

        std::vector<int> gradientVector;

        for (int iEdgePoint = 0, edgeLength = edgeIndices.size(); iEdgePoint < edgeLength; iEdgePoint++)
        {
            int centreIndex = edgeIndices[iEdgePoint];
            int centreXPos  =  centreIndex % width;
            int centreYPos  = (centreIndex - centreXPos) / width;

            int dir = calculateDirection(centreXPos - pupilXCentre, pupilYCentre - centreYPos);

            int dx = dX[dir];
            int dy = dY[dir];

            int xPos = centreXPos + dx * kernelRadius;
            int xNeg = centreXPos - dx * kernelRadius;

            int yPos = centreYPos + dy * kernelRadius;
            int yNeg = centreYPos - dy * kernelRadius;

            if (xPos < 0 || xPos >= width || yPos < 0 || yPos >= height || xNeg < 0 || xNeg >= width || yNeg < 0 || yNeg >= height)
            {
                gradientVector.push_back(0); // no gradient information available
                continue;
            }

            int neighbourIndexPositive = width * yPos + xPos;
            int neighbourIndexNegative = width * yNeg + xNeg;

            double positiveIntensity = ptr[neighbourIndexPositive];
            double negativeIntensity = ptr[neighbourIndexNegative];

            gradientVector.push_back(positiveIntensity - negativeIntensity);
        }

        return gradientVector;
    }

    std::vector<int> findEdgeIntensities(const cv::Mat& img, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties, AOIProperties mAOI)
    {
        int edgeSize = mEdgeProperties.pointIndices.size();

        int positionOffset = floor(0.5 * mDetectionParameters.cannyKernelSize);

        // calculate pixel intensities within inner curve of edge points

        std::vector<int> edgeIntensities(edgeSize);

        uchar *ptr_img = img.data;

        for (int iEdgePoint = 0; iEdgePoint < edgeSize; iEdgePoint++)
        {
            int edgePointIndex = mEdgeProperties.pointIndices[iEdgePoint];
            int edgePointXPos  = edgePointIndex % mAOI.wdth;
            int edgePointYPos  = (edgePointIndex - edgePointXPos) / mAOI.wdth;

            int offsetXPos = edgePointXPos + positionOffset * ceil2(mEdgeProperties.xnormals[iEdgePoint]);
            int offsetYPos = edgePointYPos + positionOffset * ceil2(mEdgeProperties.ynormals[iEdgePoint]);

            if (offsetXPos < 0 || offsetXPos > mAOI.wdth || offsetYPos < 0 || offsetYPos > mAOI.hght)
            {       edgeIntensities[iEdgePoint] = (int) ptr_img[edgePointXPos + edgePointYPos * mAOI.wdth]; }
            else {  edgeIntensities[iEdgePoint] = (int) ptr_img[   offsetXPos +    offsetYPos * mAOI.wdth]; }
        }

        return edgeIntensities;
    }

    std::vector<edgeProperties> edgeTerminalFilter(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties, double thresholdScoreDiffEdge)
    {
        // start from end and move through edge with window
        // if window average is significantly below central average --> remove

        std::vector<edgeProperties> vEdgePropertiesAll; // all recorded edges
        std::vector<edgeProperties> vEdgePropertiesOld(1); // edges to be checked
        vEdgePropertiesOld[0] = mEdgeProperties;

        do
        {
            std::vector<edgeProperties> vEdgePropertiesNew;

            for (int iEdge = 0, numEdges = vEdgePropertiesOld.size(); iEdge < numEdges; iEdge++)
            {
                edgeProperties mEdgePropertiesOld = vEdgePropertiesOld[iEdge];

                int edgeSize = mEdgePropertiesOld.pointIndices.size();

                if (edgeSize <= 3 * mDetectionParameters.windowLengthEdge) // don't segment short edges
                {
                    vEdgePropertiesAll.push_back(mEdgePropertiesOld);
                    continue;
                }

                int iEdgePoint  = 0;
                bool BREAK_LOOP = false;

                std::vector<int> breakPointsAll;
                std::vector<double> scoreDifference;
                int iBreakPoint = 0;

                do // run through all edge points
                {
                    iEdgePoint = iEdgePoint + mDetectionParameters.windowLengthEdge;

                    if (iEdgePoint >= edgeSize - mDetectionParameters.windowLengthEdge - 1)
                    {
                        BREAK_LOOP = true;
                        iEdgePoint = edgeSize - mDetectionParameters.windowLengthEdge - 1;
                    }

                    std::vector<int> breakPoints = {0, iEdgePoint, iEdgePoint + 1, edgeSize - 1};

                    std::vector<edgeProperties> vEdgePropertiesTemp(2);

                    for (int i = 0; i < 2; i++) // split edge up into two
                    {
                        edgeProperties mEdgePropertiesNew;

                        mEdgePropertiesNew.pointIndices = std::vector<int>   (mEdgePropertiesOld.pointIndices.begin() + breakPoints[2 * i], mEdgePropertiesOld.pointIndices.begin() + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.intensities  = std::vector<int>   (mEdgePropertiesOld.intensities.begin()  + breakPoints[2 * i], mEdgePropertiesOld.intensities.begin()  + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.gradients    = std::vector<int>   (mEdgePropertiesOld.gradients.begin()    + breakPoints[2 * i], mEdgePropertiesOld.gradients.begin()    + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.radii        = std::vector<double>(mEdgePropertiesOld.radii.begin()        + breakPoints[2 * i], mEdgePropertiesOld.radii.begin()        + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.curvatures   = std::vector<double>(mEdgePropertiesOld.curvatures.begin()   + breakPoints[2 * i], mEdgePropertiesOld.curvatures.begin()   + breakPoints[2 * i + 1]);

                        mEdgePropertiesNew.radius    = calculateMean      (mEdgePropertiesNew.radii);
                        mEdgePropertiesNew.intensity = calculateMeanInt   (mEdgePropertiesNew.intensities);
                        mEdgePropertiesNew.gradient  = calculateMeanInt   (mEdgePropertiesNew.gradients);
                        mEdgePropertiesNew.length    = mEdgePropertiesNew.pointIndices.size();
                        calculateCurvatureStats(mDetectionParameters, mEdgePropertiesNew); // get average curvature

                        vEdgePropertiesTemp[i] = mEdgePropertiesNew;
                    }

                    // Calculate score difference

                    std::vector<double> inputVector_1(7);
                    inputVector_1[0] = vEdgePropertiesTemp[0].radius    / (mDetectionVariables.predictedCircumference / (2 * M_PI));
                    inputVector_1[1] = 0; // don't use circumference difference
                    inputVector_1[2] = vEdgePropertiesTemp[0].curvature - mDetectionVariables.predictedCurvature;
                    inputVector_1[3] = vEdgePropertiesTemp[0].intensity - mDetectionVariables.predictedIntensity;
                    inputVector_1[4] = vEdgePropertiesTemp[0].gradient  - mDetectionVariables.predictedGradient;
                    inputVector_1[5] = 0; // don't use radial variance
                    inputVector_1[6] = vEdgePropertiesTemp[0].length;

                    std::vector<double> inputVector_2(7);
                    inputVector_2[0] = vEdgePropertiesTemp[1].radius    / (mDetectionVariables.predictedCircumference / (2 * M_PI));
                    inputVector_2[1] = 0; // don't use circumference difference
                    inputVector_2[2] = vEdgePropertiesTemp[1].curvature - mDetectionVariables.predictedCurvature;
                    inputVector_2[3] = vEdgePropertiesTemp[1].intensity - mDetectionVariables.predictedIntensity;
                    inputVector_2[4] = vEdgePropertiesTemp[1].gradient  - mDetectionVariables.predictedGradient;
                    inputVector_2[5] = 0; // don't use radial variance
                    inputVector_2[6] = vEdgePropertiesTemp[1].length;

                    double score_1 = calculateScoreTotal(mDetectionVariables, inputVector_1, true);
                    double score_2 = calculateScoreTotal(mDetectionVariables, inputVector_2, true);

                    scoreDifference.push_back(std::abs(score_1 - score_2));
                    breakPointsAll.push_back(iEdgePoint);
                    iBreakPoint++;
                } while (!BREAK_LOOP);

                // Find point of maximum score difference

                auto itr = std::max_element(scoreDifference.begin(), scoreDifference.end());
                double scoreDifferenceMax = *itr;
                int indexMax = std::distance(scoreDifference.begin(), itr);

                std::vector<edgeProperties> vEdgePropertiesTemp;

                if (scoreDifferenceMax > thresholdScoreDiffEdge)
                {
                    std::vector<int> breakPoints = {0, breakPointsAll[indexMax], breakPointsAll[indexMax] + 1, edgeSize - 1};

                    for (int i = 0; i < 2; i++) // split edge up
                    {
                        edgeProperties mEdgePropertiesNew;

                        mEdgePropertiesNew.pointIndices = std::vector<int>   (mEdgePropertiesOld.pointIndices.begin() + breakPoints[2 * i], mEdgePropertiesOld.pointIndices.begin() + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.intensities  = std::vector<int>   (mEdgePropertiesOld.intensities.begin()  + breakPoints[2 * i], mEdgePropertiesOld.intensities.begin()  + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.gradients    = std::vector<int>   (mEdgePropertiesOld.gradients.begin()    + breakPoints[2 * i], mEdgePropertiesOld.gradients.begin()    + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.radii        = std::vector<double>(mEdgePropertiesOld.radii.begin()        + breakPoints[2 * i], mEdgePropertiesOld.radii.begin()        + breakPoints[2 * i + 1]);
                        mEdgePropertiesNew.curvatures   = std::vector<double>(mEdgePropertiesOld.curvatures.begin()   + breakPoints[2 * i], mEdgePropertiesOld.curvatures.begin()   + breakPoints[2 * i + 1]);

                        vEdgePropertiesTemp.push_back(mEdgePropertiesNew);
                    }
                }
                else { vEdgePropertiesTemp.push_back(mEdgePropertiesOld); }

                if      (vEdgePropertiesTemp.size()  > 1) { vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end()); } // edges to be checked
                else if (vEdgePropertiesTemp.size() == 1) { vEdgePropertiesAll.push_back(vEdgePropertiesTemp[0]); } // accepted edges
            }

            vEdgePropertiesOld = vEdgePropertiesNew;

        } while (vEdgePropertiesOld.size() > 1);

        return vEdgePropertiesAll;
    }

    std::vector<double> calculateEdgeRadii(const edgeProperties& mEdgeProperties, AOIProperties mAOI, double xCentre, double yCentre)
    {
        // calculate distance between each edge point and expected pupil centre

        int edgeLength = mEdgeProperties.pointIndices.size();

        std::vector<double> edgePointRadii(edgeLength);

        for (int iEdgePoint = 0; iEdgePoint < edgeLength; iEdgePoint++)
        {
            int edgePointIndex = mEdgeProperties.pointIndices[iEdgePoint];
            int edgePointXPos  =  edgePointIndex % mAOI.wdth;
            int edgePointYPos  = (edgePointIndex - edgePointXPos) / mAOI.wdth;

            double dX = xCentre - edgePointXPos;
            double dY = yCentre - edgePointYPos;

            edgePointRadii[iEdgePoint] = sqrt(dX * dX + dY * dY);
        }

        return edgePointRadii;
    }

    void restoreEdgePoints(edgeProperties& mEdgeProperties, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
    {
        // Add additional adjacent indices that were removed by morphological operation

        std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
        std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

        int edgeSize = mEdgeProperties.pointIndices.size();

        for (int iEdgePoint = 0; iEdgePoint < edgeSize; iEdgePoint++)
        {
            int centreIndex = mEdgeProperties.pointIndices[iEdgePoint];

            int centreXPos = centreIndex % mAOI.wdth;
            int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;

            for (int m = 0; m < 8; m++) // loop through 8-connected environment
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];

                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                if (cannyEdgeVector[neighbourIndex] == -1) // if neighbouring point was canny edge point that was removed by morphological operation then ...
                {
                    cannyEdgeVector[neighbourIndex] = 5; // ... tag it and ...
                    mEdgeProperties.pointIndices.push_back(neighbourIndex); // ... add it to the (partial) edge
                }
            }
        }
    }

    std::vector<edgeProperties> removeShortEdges(const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgeProperties)
    {
        std::vector<edgeProperties> vEdgePropertiesNew;
        for (int iEdge = 0, numEdges = vEdgeProperties.size(); iEdge < numEdges; iEdge++) // ignore short edges
        {
            edgeProperties mEdgeProperties = vEdgeProperties[iEdge];
            int edgeSize = mEdgeProperties.pointIndices.size();
            if (edgeSize >= mDetectionParameters.windowLengthEdge) { vEdgePropertiesNew.push_back(mEdgeProperties); }
        }

        return vEdgePropertiesNew;
    }

    std::vector<int> edgeClassification(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgePropertiesAll)
    {
        int numEdgesMax = mDetectionParameters.fitEdgeMaximum;
        int numEdges    = vEdgePropertiesAll.size();
        if (numEdgesMax > numEdges) { numEdgesMax = numEdges; }

        // Classify edges based on score

        std::vector<double> totalScores(numEdges);
        std::vector<int> pupilEdges;

        if (mDetectionVariables.certaintyAverages > certaintyThreshold || mDetectionVariables.certaintyFeatures > certaintyThreshold || mDetectionVariables.certaintyPosition > certaintyThreshold) // need to have enough certainty
        {
            for (int iEdge = 0; iEdge < numEdges; iEdge++)
            {
                std::vector<double> inputVector(7);

                inputVector[0] = vEdgePropertiesAll[iEdge].radius    / (mDetectionVariables.predictedCircumference / (2 * M_PI));
                inputVector[1] = vEdgePropertiesAll[iEdge].length    / mDetectionVariables.predictedCircumference;
                inputVector[2] = vEdgePropertiesAll[iEdge].curvature - mDetectionVariables.predictedCurvature;
                inputVector[3] = vEdgePropertiesAll[iEdge].intensity - mDetectionVariables.predictedIntensity;
                inputVector[4] = vEdgePropertiesAll[iEdge].gradient  - mDetectionVariables.predictedGradient;
                inputVector[5] = vEdgePropertiesAll[iEdge].radiusVar;
                inputVector[6] = vEdgePropertiesAll[iEdge].length;

                totalScores[iEdge] = calculateScoreTotal(mDetectionVariables, inputVector, true);
            }

            // Only pick edges above threshold

            for (int iEdge = 0; iEdge < numEdges; iEdge++)
            {
                if (totalScores[iEdge] >= mDetectionVariables.thresholdScore)
                {
                    pupilEdges.push_back(iEdge);
                }
            }
        }
        else // take darkest edges
        {
            for (int iEdge = 0; iEdge < numEdges; iEdge++)
            {
                totalScores[iEdge] = 1 / (vEdgePropertiesAll[iEdge].intensity + 1); // avoids inf
                pupilEdges.push_back(iEdge);
            }
        }

        int numEdgesNew = pupilEdges.size();

        if (numEdgesNew > numEdgesMax) // Grab edges with highest score if maximum is exceeded
        {
            std::vector<double> totalScoresUnsorted;

            for (int iEdge = 0; iEdge < numEdgesNew; iEdge++)
            {
                totalScoresUnsorted.push_back(totalScores[pupilEdges[iEdge]]);
            }

            std::vector<int> acceptedEdges(numEdgesMax);

            std::vector<double> totalScoresSorted = totalScoresUnsorted;
            std::sort   (totalScoresSorted.begin(), totalScoresSorted.end());
            std::reverse(totalScoresSorted.begin(), totalScoresSorted.end());

            for (int iEdge = 0; iEdge < numEdgesMax; iEdge++) // THRESHOLD: do not handle more than a fixed number of edges
            {
                for (int jEdge = 0; jEdge < numEdgesNew; jEdge++)
                {
                    if (totalScoresSorted[iEdge] == totalScoresUnsorted[jEdge])
                    {
                        acceptedEdges[iEdge] = pupilEdges[jEdge];
                        totalScoresUnsorted[jEdge] = -1;
                        break;
                    }
                }
            }

            return acceptedEdges;
        }
        else
        {
            return pupilEdges;
        }
    }

    std::vector<double> EllipseRotationTransformation(const std::vector<double>& c)
    {
        double A = c[0];
        double B = c[1];
        double C = c[2];
        double D = c[3];
        double E = c[4];
        double F = c[5];

        double alpha = 0.5 * (atan2(B, (A - C))); // rotation angle

        double AA =  A * cos(alpha) * cos(alpha) + B * cos(alpha) * sin(alpha) + C * sin(alpha) * sin(alpha);
        double CC =  A * sin(alpha) * sin(alpha) - B * cos(alpha) * sin(alpha) + C * cos(alpha) * cos(alpha);
        double DD =  D * cos(alpha) + E * sin(alpha);
        double EE = -D * sin(alpha) + E * cos(alpha);
        double FF =  F;

        // semi axes

        double a = sqrt((-4 * FF * AA * CC + CC * DD * DD + AA * EE * EE)/(4 * AA * CC * CC));
        double b = sqrt((-4 * FF * AA * CC + CC * DD * DD + AA * EE * EE)/(4 * AA * AA * CC));

        double semiMajor = 0;
        double semiMinor = 0;

        if (a >= b) {
            semiMajor = a;
            semiMinor = b;
        }
        else
        {
            semiMajor = b;
            semiMinor = a;
        }

        // coordinates of centre point

        double x = -(DD / (2 * AA)) * cos(alpha) + (EE / (2 * CC)) * sin(alpha);
        double y = -(DD / (2 * AA)) * sin(alpha) - (EE / (2 * CC)) * cos(alpha);

        // width and height

        double w = 2 * sqrt(pow(semiMajor * cos(alpha), 2) + pow(semiMinor * sin(alpha), 2));
        double h = 2 * sqrt(pow(semiMajor * sin(alpha), 2) + pow(semiMinor * cos(alpha), 2));

        std::vector<double> v(7);

        v[0] = semiMajor;
        v[1] = semiMinor;
        v[2] = x;
        v[3] = y;
        v[4] = w;
        v[5] = h;
        v[6] = alpha;

        return v;
    }

    ellipseProperties fitEllipse(std::vector<int> edgeIndices, int edgeSetSize, const AOIProperties& mAOI)
    {
        Eigen::MatrixXd ConstraintMatrix(6, 6); // constraint matrix
        ConstraintMatrix <<  0,  0,  2,  0,  0,  0,
                0, -1,  0,  0,  0,  0,
                2,  0,  0,  0,  0,  0,
                0,  0,  0,  0,  0,  0,
                0,  0,  0,  0,  0,  0,
                0,  0,  0,  0,  0,  0;

        // least squares ellipse fitting

        Eigen::MatrixXd DesignMatrix(edgeSetSize, 6); // design matrix

        for (int iEdgePoint = 0; iEdgePoint < edgeSetSize; iEdgePoint++)
        {
            int edgePointIndex = edgeIndices[iEdgePoint];

            double edgePointX = edgePointIndex % mAOI.wdth;
            double edgePointY = (edgePointIndex - edgePointX) / mAOI.wdth;

            DesignMatrix(iEdgePoint, 0) = edgePointX * edgePointX;
            DesignMatrix(iEdgePoint, 1) = edgePointX * edgePointY;
            DesignMatrix(iEdgePoint, 2) = edgePointY * edgePointY;
            DesignMatrix(iEdgePoint, 3) = edgePointX;
            DesignMatrix(iEdgePoint, 4) = edgePointY;
            DesignMatrix(iEdgePoint, 5) = 1;
        }

        Eigen::MatrixXd ScatterMatrix(6, 6); // scatter matrix
        ScatterMatrix = DesignMatrix.transpose() * DesignMatrix;

        // solving eigensystem

        Eigen::MatrixXd EigenSystem(6, 6);
        EigenSystem = ScatterMatrix.inverse() * ConstraintMatrix;

        Eigen::EigenSolver<Eigen::MatrixXd> EigenSolver(EigenSystem);

        Eigen::VectorXd EigenValues  = EigenSolver.eigenvalues().real();
        Eigen::MatrixXd EigenVectors = EigenSolver.eigenvectors().real();

        double minEigenValue   = std::numeric_limits<double>::max(); // set to maximum
        int minEigenValueIndex = 0;

        for (int iEigenValue = 0; iEigenValue < 6; iEigenValue++)
        {
            if (EigenValues(iEigenValue) < minEigenValue && EigenValues(iEigenValue) > 0.00000000001)
            {
                minEigenValueIndex = iEigenValue;
                minEigenValue = EigenValues(iEigenValue);
            }
        }

        Eigen::VectorXd eigenVector = EigenVectors.col(minEigenValueIndex);
        double normalizationFactor = eigenVector.transpose() * ConstraintMatrix * eigenVector;

        std::vector<double> ellipseFitCoefficients(6); // ellipse parameters

        for (int iCoefs = 0; iCoefs < 6; iCoefs++)
        { ellipseFitCoefficients[iCoefs] = (1 / sqrt(normalizationFactor)) * eigenVector(iCoefs); }

        // calculate size, shape and position of ellipse

        std::vector<double> ellipseParameters = EllipseRotationTransformation(ellipseFitCoefficients);

        double semiMajor = ellipseParameters[0];
        double semiMinor = ellipseParameters[1];

        ellipseProperties mEllipseProperties;
        mEllipseProperties.DETECTED      = true;
        mEllipseProperties.circumference = ramanujansApprox(semiMajor,semiMinor);
        mEllipseProperties.aspectRatio   = semiMinor / semiMajor;
        mEllipseProperties.xPos          = ellipseParameters[2] + mAOI.xPos;
        mEllipseProperties.yPos          = ellipseParameters[3] + mAOI.yPos;
        mEllipseProperties.width         = ellipseParameters[4];
        mEllipseProperties.height        = ellipseParameters[5];
        mEllipseProperties.angle         = ellipseParameters[6];
        mEllipseProperties.coefficients  = ellipseFitCoefficients;

        for (int iParameter = 0; iParameter < 6; iParameter++)
        {
            if (std::isnan(ellipseParameters[iParameter])) { mEllipseProperties.DETECTED = false; } // ERROR
        }

        return mEllipseProperties;
    }

    std::vector<ellipseProperties> getEllipseFits(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI)
    {
        int numEdgesTotal = vEdgePropertiesAll.size(); // total number of edges

        std::vector<ellipseProperties> vEllipsePropertiesAll; // vector to record information for each accepted ellipse fit

        for (int combiNumEdges = numEdgesTotal; combiNumEdges >= 1; combiNumEdges--) // loop through all possible edge set sizes
        {
            std::vector<bool> edgeCombination(numEdgesTotal);
            std::fill(edgeCombination.begin() + numEdgesTotal - combiNumEdges, edgeCombination.end(), true);

            do // loop through all possible edge combinations for the current set size
            {
                std::vector<int> combiEdgeIndices (combiNumEdges);
                std::vector<int> combiEdgeLengths (combiNumEdges);
                std::vector<int> combiEdgeSizes   (combiNumEdges);

                std::vector<std::vector<int>> combiEdgePoints(combiNumEdges);

                for (int iEdge = 0, jEdge = 0; iEdge < numEdgesTotal; ++iEdge)
                {
                    if (edgeCombination[iEdge])
                    {
                        combiEdgeIndices     [jEdge] = vEdgePropertiesAll[iEdge].index;
                        combiEdgeLengths     [jEdge] = vEdgePropertiesAll[iEdge].length;
                        combiEdgeSizes       [jEdge] = vEdgePropertiesAll[iEdge].size;
                        combiEdgePoints[jEdge] = vEdgePropertiesAll[iEdge].pointIndices;
                        jEdge++;
                    }
                }

                // Find properties of edge collection

                // Calculate total length

                int edgeSetLength = std::accumulate(combiEdgeLengths.begin(), combiEdgeLengths.end(), 0);
                int edgeSetSize   = std::accumulate(combiEdgeSizes.begin(),   combiEdgeSizes.end(),   0);

                // Concatenate index vectors

                std::vector<int> edgeIndices; // vector containing all indices for fit
                edgeIndices.reserve(edgeSetSize); // preallocate memory

                for (int iEdge = 0; iEdge < combiNumEdges; iEdge++)
                { edgeIndices.insert(edgeIndices.end(), combiEdgePoints[iEdge].begin(), combiEdgePoints[iEdge].end()); }

                // Calculate range

                std::vector<int> edgeXPositions(edgeSetSize);
                std::vector<int> edgeYPositions(edgeSetSize);

                for (int iEdgePoint = 0; iEdgePoint < edgeSetSize; iEdgePoint++)
                {
                    int edgePointIndex = edgeIndices[iEdgePoint];
                    edgeXPositions[iEdgePoint] = edgePointIndex % mAOI.wdth;
                    edgeYPositions[iEdgePoint] = (edgePointIndex - edgeXPositions[iEdgePoint]) / mAOI.wdth;
                }

                int XPosMin = *std::min_element(std::begin(edgeXPositions), std::end(edgeXPositions));
                int XPosMax = *std::max_element(std::begin(edgeXPositions), std::end(edgeXPositions));
                int YPosMin = *std::min_element(std::begin(edgeYPositions), std::end(edgeYPositions));
                int YPosMax = *std::max_element(std::begin(edgeYPositions), std::end(edgeYPositions));

                double combiWdth = 0.5 * (XPosMax - XPosMin);
                double combiHght = 0.5 * (YPosMax - YPosMin);

                double semiMajor, semiMinor;
                if (combiWdth > combiHght) { semiMajor = combiWdth; semiMinor = combiHght; }
                else                       { semiMajor = combiHght; semiMinor = combiWdth; }

                double circumferenceApprox = ramanujansApprox(semiMajor,semiMinor);

                double circumferenceUpperLimit = mDetectionVariables.averageCircumference + (mDetectionVariables.averageCircumference * mDetectionVariables.offsetCircumference);
                double circumferenceLowerLimit = mDetectionVariables.averageCircumference - (mDetectionVariables.averageCircumference * mDetectionVariables.offsetCircumference);

                if (circumferenceUpperLimit > mDetectionParameters.circumferenceMax) { circumferenceUpperLimit = mDetectionParameters.circumferenceMax; }
                if (circumferenceLowerLimit < mDetectionParameters.circumferenceMin) { circumferenceLowerLimit = mDetectionParameters.circumferenceMin; }

                if (circumferenceApprox > circumferenceUpperLimit) { continue; } // no large ellipse
                if (circumferenceApprox < circumferenceLowerLimit) { continue; } // no small ellipse

                // Fit ellipse

                ellipseProperties mEllipseProperties = fitEllipse(edgeIndices, edgeSetSize, mAOI);

                if (!mEllipseProperties.DETECTED) { continue; } // error

                if (edgeSetLength < fitEdgeFractionMin * mEllipseProperties.circumference) { continue; } // minimum number of edge points required

                // Absolute displacement filter

                //            double dX = mEllipseProperties.xPos - mDetectionVariables.predictedXPos;
                //            double dY = mEllipseProperties.yPos - mDetectionVariables.predictedYPos;
                //            double dR = sqrt(dX * dX + dY * dY);
                //            if (dR > mDetectionVariables.thresholdChangePosition) { continue; } // no large ellipse displacements

                // Absolute size and shape filter

                if (mEllipseProperties.circumference > circumferenceUpperLimit) { continue; } // no large ellipse
                if (mEllipseProperties.circumference < circumferenceLowerLimit) { continue; } // no small ellipse
                if (mEllipseProperties.aspectRatio   < mDetectionParameters.aspectRatioMin)   { continue; } // no extreme deviations from circular shape

                // Size and shape combination filter

                circumferenceUpperLimit = aspectRatioSlope * mEllipseProperties.aspectRatio + mDetectionParameters.circumferenceMax - aspectRatioSlope;
                if (mEllipseProperties.circumference > circumferenceUpperLimit) { continue; } // no large ellipse

                // Relative change in size and shape filter

                double differenceAspectRatio   = std::abs(mEllipseProperties.aspectRatio   - mDetectionVariables.predictedAspectRatio); // absolute difference
                double differenceCircumference = std::abs(mEllipseProperties.circumference - mDetectionVariables.predictedCircumference) / std::max(mEllipseProperties.circumference, mDetectionVariables.predictedCircumference); // relative difference

                if (differenceAspectRatio   > mDetectionVariables.thresholdChangeAspectRatio  ) { continue; } // no large ellipse shape changes
                if (differenceCircumference > mDetectionVariables.thresholdChangeCircumference) { continue; } // no large ellipse size changes

                // Calculate error between fit and every edge point

                double A = mEllipseProperties.coefficients[0];
                double B = mEllipseProperties.coefficients[1];
                double C = mEllipseProperties.coefficients[2];
                double D = mEllipseProperties.coefficients[3];
                double E = mEllipseProperties.coefficients[4];
                double F = mEllipseProperties.coefficients[5];

                std::vector<double> fitErrors(edgeSetSize);

                for (int iEdgePoint = 0; iEdgePoint < edgeSetSize; iEdgePoint++)
                {
                    int edgePointIndex = edgeIndices[iEdgePoint];

                    double x = edgePointIndex % mAOI.wdth;
                    double y = (edgePointIndex - x) / mAOI.wdth;

                    fitErrors[iEdgePoint] = std::abs(A * x * x + B * x * y + C * y * y + D * x + E * y + F);
                }

                std::vector<double> fitErrorsSorted = fitErrors;
                std::sort   (fitErrorsSorted.begin(), fitErrorsSorted.end());
                std::reverse(fitErrorsSorted.begin(), fitErrorsSorted.end());
                std::vector<double> fitErrorsMax(fitErrorsSorted.begin(), fitErrorsSorted.begin() + round(mDetectionParameters.fitEdgeFraction * edgeSetLength));
                double fitErrorAbsolute = calculateMean(fitErrorsMax); // absolute error
                double fitErrorRelative = (fitErrorAbsolute + 0.5588) / mEllipseProperties.circumference; // relative error

                if (fitErrorRelative > mDetectionParameters.thresholdFitError) { continue; } // no large fit errors

                // Save parameters of accepted fit

                mEllipseProperties.fitError    = fitErrorRelative;
                mEllipseProperties.edgeIndices = combiEdgeIndices;
                mEllipseProperties.edgeLength  = edgeSetLength;
                vEllipsePropertiesAll.push_back(mEllipseProperties);
            }
            while (std::next_permutation(edgeCombination.begin(), edgeCombination.end()));
        }

        return vEllipsePropertiesAll;
    }

    std::vector<int> ellipseFitFilter(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, std::vector<ellipseProperties> vEllipseProperties)
    {
        static const double scoreFactorAspectRatio   = 0.815;
        static const double scoreFactorCircumference = 0.415;
        static const double scoreFactorLength        = 1.325;
        static const double scoreFactorError         = 0.239;
        static const double scoreFactorAngle         = 1.014;
        static const double scoreFactorAngleFunction = 0.787;

        static const std::vector<double> parametersAspectRatio   = {0.70922, 0.0005102, 0.0087285, 0.26935, 0.0005102, 0.028891};
        static const std::vector<double> parametersCircumference = {0.74258, 1.0003,    0.0051637, 0.23681, 1.0061,    0.030051};
        static const std::vector<double> parametersLength        = {0.32519, 0.9949,    0.14673,   0.67869, 0.9949,    0.32326 };
        static const std::vector<double> parametersError         = {0.28342, 0.20204,   0.024372,  0.68972, 0.20204,   0.15905 };
        static const std::vector<double> parametersAngle         = {0.31997, 0.0010204, 0.01714,   0.66619, 0.0010204, 0.080493};

        int numFits = vEllipseProperties.size();

        std::vector<double> scoreFits(numFits);

        double certaintyFactorFeatures = 0.5 * (mDetectionVariables.certaintyFeatures + 1);
        double certaintyFactorAverages = 0.5 * (mDetectionVariables.certaintyAverages + 1);

        if (certaintyFactorAverages > certaintyFactorFeatures + certaintyOffset) { certaintyFactorFeatures = certaintyFactorAverages; }

        for (int iFit = 0; iFit < numFits; iFit++)
        {
            ellipseProperties mEllipseProperties = vEllipseProperties[iFit];

            double errorAspectRatio   = std::abs(mEllipseProperties.aspectRatio - mDetectionVariables.predictedAspectRatio);
            double errorAngle         = std::abs(mEllipseProperties.angle       - mDetectionVariables.predictedAngle);
            double ratioCircumference = mEllipseProperties.circumference / mDetectionVariables.predictedCircumference;
            double ratioLength        = mEllipseProperties.edgeLength    / mDetectionVariables.predictedCircumference;
            double fitError           = mEllipseProperties.fitError;

            double factorAngleFunction = scoreFactorAngleFunction * (1 - mEllipseProperties.aspectRatio) + (1 - scoreFactorAngleFunction);

            double factorAngle         = certaintyFactorFeatures * scoreFactorAngle         * factorAngleFunction;
            double factorAspectRatio   = certaintyFactorFeatures * scoreFactorAspectRatio;
            double factorCircumference = certaintyFactorFeatures * scoreFactorCircumference;
            double factorLength        = certaintyFactorFeatures * scoreFactorLength;
            double factorFitError      =                           scoreFactorError;

            // Calculate scores

            double scoreAngle         = factorAngle         * calculateGaussian2(errorAngle,         parametersAngle         );
            double scoreAspectRatio   = factorAspectRatio   * calculateGaussian2(errorAspectRatio,   parametersAspectRatio   );
            double scoreCircumference = factorCircumference * calculateGaussian2(ratioCircumference, parametersCircumference );
            double scoreLength        = factorLength        * calculateGaussian2(ratioLength,        parametersLength        );
            double scoreFitError      = factorFitError      * calculateGaussian2(fitError,           parametersError         );

            double norm =  factorAngle + factorAspectRatio + factorCircumference + factorLength + factorFitError;

            double scoreTotal = 0;
            if (norm > 0) { scoreTotal = (scoreAngle + scoreAspectRatio + scoreCircumference + scoreLength + scoreFitError) / norm; }

            scoreFits[iFit] = scoreTotal;
        }

        // Find all accepted fits

        std::vector<int> acceptedIndices(1);

        if (mDetectionVariables.certaintyPosition > certaintyThreshold || mDetectionVariables.certaintyFeatures > certaintyThreshold || mDetectionVariables.certaintyAverages > certaintyThreshold)
        {
            std::vector<double> scoreFitsSorted = scoreFits;

            std::sort   (scoreFitsSorted.begin(), scoreFitsSorted.end());
            std::reverse(scoreFitsSorted.begin(), scoreFitsSorted.end());

            double scoreFitMax = scoreFitsSorted[0];
            std::vector<double> acceptedScores = {scoreFitMax};

            for (int iFit = 1; iFit < numFits; iFit++)
            {
                double scoreFit = scoreFitsSorted[iFit];

                if (scoreFit + mDetectionParameters.thresholdScoreDiffFit >= scoreFitMax)
                {
                    acceptedScores.push_back(scoreFit);
                }
            }

            int numFitsAccepted = acceptedScores.size();

            acceptedIndices.resize(numFitsAccepted);

            for (int iFit = 0; iFit < numFitsAccepted; iFit++)
            {
                double scoreFit = acceptedScores[iFit];

                for (int jFit = 0; jFit < numFits; jFit++)
                {
                    if (scoreFit == scoreFits[jFit])
                    {
                        acceptedIndices[iFit] = jFit;
                        scoreFits[jFit] = -std::numeric_limits<double>::max();
                    }
                }
            }
        }
        else
        {
            acceptedIndices[0] = std::distance(scoreFits.begin(), std::max_element(scoreFits.begin(), scoreFits.end()));
        }

        return acceptedIndices;
    }


    void ellipseFitAverage(ellipseProperties& mEllipseProperties, const std::vector<ellipseProperties>& vEllipseProperties)
    {
        int numFits = vEllipseProperties.size();

        std::vector<double> aspectRatios    (numFits);
        std::vector<double> circumferences  (numFits);
        std::vector<double> widths          (numFits);
        std::vector<double> heights         (numFits);
        std::vector<double> angles          (numFits);
        std::vector<double> xPositions      (numFits);
        std::vector<double> yPositions      (numFits);

        for (int iFit = 0; iFit < numFits; iFit++)
        {
            ellipseProperties mEllipsePropertiesTemp = vEllipseProperties[iFit];
            aspectRatios    [iFit] = mEllipsePropertiesTemp.aspectRatio;
            circumferences  [iFit] = mEllipsePropertiesTemp.circumference;
            widths          [iFit] = mEllipsePropertiesTemp.width;
            heights         [iFit] = mEllipsePropertiesTemp.height;
            angles          [iFit] = mEllipsePropertiesTemp.angle;
            xPositions      [iFit] = mEllipsePropertiesTemp.xPos;
            yPositions      [iFit] = mEllipsePropertiesTemp.yPos;
        }

        mEllipseProperties.aspectRatio      = calculateMean(aspectRatios);
        mEllipseProperties.circumference    = calculateMean(circumferences);
        mEllipseProperties.width            = calculateMean(widths);
        mEllipseProperties.height           = calculateMean(heights);
        mEllipseProperties.angle            = calculateMean(angles);
        mEllipseProperties.xPos             = calculateMean(xPositions);
        mEllipseProperties.yPos             = calculateMean(yPositions);
    }

    void setCurvatureMeasurement(detectionParameters& mDetectionParameters, int imgWdth)
    {
        // For development

        mDetectionParameters.windowLengthEdge = mDetectionParameters.windowLengthEdge;

        mDetectionParameters.aspectRatioMin   = 0.0;
        mDetectionParameters.circumferenceMax = M_PI * imgWdth;
        mDetectionParameters.circumferenceMin = 1.0; // 1.0 avoids inf

        mDetectionParameters.curvatureOffset = 360;

        mDetectionParameters.thresholdScore     = 0.0;
        mDetectionParameters.thresholdScoreDiffEdge = 1.0;

        mDetectionParameters.glintWdth = 0.0;
    }

    void checkVariableLimits(detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters)
    {
        if (mDetectionVariables.thresholdScore > mDetectionParameters.thresholdScore)
        {   mDetectionVariables.thresholdScore = mDetectionParameters.thresholdScore; }

        if (mDetectionVariables.thresholdChangeAspectRatio < mDetectionParameters.thresholdChangeAspectRatio)
        {   mDetectionVariables.thresholdChangeAspectRatio = mDetectionParameters.thresholdChangeAspectRatio; }

        if (mDetectionVariables.thresholdChangeCircumference < mDetectionParameters.thresholdChangeCircumference)
        {   mDetectionVariables.thresholdChangeCircumference = mDetectionParameters.thresholdChangeCircumference; }

        if (mDetectionVariables.thresholdChangePosition < mDetectionParameters.thresholdChangePosition)
        {   mDetectionVariables.thresholdChangePosition = mDetectionParameters.thresholdChangePosition; }

        if (mDetectionVariables.offsetCircumference < mDetectionParameters.circumferenceOffset)
        {   mDetectionVariables.offsetCircumference = mDetectionParameters.circumferenceOffset; }

        if      (mDetectionVariables.certaintyPositionPrime > 1.0) { mDetectionVariables.certaintyPositionPrime = 1.0; }
        else if (mDetectionVariables.certaintyPositionPrime < 0.0) { mDetectionVariables.certaintyPositionPrime = 0.0; }

        if      (mDetectionVariables.certaintyFeaturesPrime > 1.0) { mDetectionVariables.certaintyFeaturesPrime = 1.0; }
        else if (mDetectionVariables.certaintyFeaturesPrime < 0.0) { mDetectionVariables.certaintyFeaturesPrime = 0.0; }

        if      (mDetectionVariables.certaintyAveragesPrime > 1.0) { mDetectionVariables.certaintyAveragesPrime = 1.0; }
        else if (mDetectionVariables.certaintyAveragesPrime < 0.0) { mDetectionVariables.certaintyAveragesPrime = 0.0; }

    }

    detectionVariables eyeStalker(const cv::Mat& imageOriginalBGR, const AOIProperties& mAOI, detectionVariables& mDetectionVariables, detectionParameters& mDetectionParameters, dataVariables& mDataVariables, drawVariables& mDrawVariables, const developmentOptions& mAdvancedOptions)
    {
        int imageWdth = imageOriginalBGR.cols;
        int imageHght = imageOriginalBGR.rows;

        checkVariableLimits(mDetectionVariables, mDetectionParameters); // keep variables within limits

        if (mAdvancedOptions.CURVATURE_MEASUREMENT) { setCurvatureMeasurement(mDetectionParameters, imageWdth); }

        // Define search area

        double AOIOffset = mDetectionVariables.thresholdChangePosition + (mDetectionVariables.predictedCircumference * mDetectionVariables.thresholdChangeCircumference) / (2 * M_PI);

        AOIProperties searchAOI;
        //    searchAOI.xPos = round(mDetectionVariables.predictedXPos - AOIOffset - 0.5 * mDetectionVariables.predictedWidth);  // needs to be uncommented
        //    searchAOI.yPos = round(mDetectionVariables.predictedYPos - AOIOffset - 0.5 * mDetectionVariables.predictedHeight);
        //    int searchEndX = round(mDetectionVariables.predictedXPos + AOIOffset + 0.5 * mDetectionVariables.predictedWidth);
        //    int searchEndY = round(mDetectionVariables.predictedYPos + AOIOffset + 0.5 * mDetectionVariables.predictedHeight);

        searchAOI.xPos = 0; // needs to be removed
        searchAOI.yPos = 0;
        int searchEndX = imageWdth;
        int searchEndY = imageHght;

        if (searchAOI.xPos < mAOI.xPos)
        {   searchAOI.xPos = mAOI.xPos; }
        if (searchAOI.yPos < mAOI.yPos)
        {   searchAOI.yPos = mAOI.yPos; }
        if (searchEndX > mAOI.xPos + mAOI.wdth - 1)
        {   searchEndX = mAOI.xPos + mAOI.wdth - 1; }
        if (searchEndY > mAOI.yPos + mAOI.hght - 1)
        {   searchEndY = mAOI.yPos + mAOI.hght - 1; }
        searchAOI.wdth = searchEndX - searchAOI.xPos + 1;
        searchAOI.hght = searchEndY - searchAOI.yPos + 1;

        AOIProperties haarAOI;
        haarAOI.wdth = round(mDetectionVariables.predictedWidth);
        haarAOI.hght = round(mDetectionVariables.predictedHeight);

        if (haarAOI.wdth > searchAOI.wdth) { haarAOI.wdth = searchAOI.wdth; }
        if (haarAOI.hght > searchAOI.hght) { haarAOI.hght = searchAOI.hght; }

        // Convert to grayscale

        cv::Mat imageOriginalGray;
        cv::cvtColor(imageOriginalBGR, imageOriginalGray, cv::COLOR_BGR2GRAY);

        ////////////////////////////////////////////////////////////////////
        /////////////////////// INITIAL DETECTION  /////////////////////////
        ////////////////////////////////////////////////////////////////////

        AOIProperties glintAOI;

        double sizeFactorUp   = 2;
        double sizeFactorDown = 1 / sizeFactorUp;

        if (searchAOI.wdth > haarAOI.wdth || searchAOI.hght > haarAOI.hght)
        {
            // Down sample image

            int imgWdthResized = round(imageWdth * sizeFactorDown);
            int imgHghtResized = round(imageHght * sizeFactorDown);

            cv::Size size(imgWdthResized, imgHghtResized);
            cv::Mat imageResized;
            cv::resize(imageOriginalGray, imageResized, size);

            AOIProperties searchAOIResized;
            searchAOIResized.xPos = sizeFactorDown * searchAOI.xPos;
            searchAOIResized.yPos = sizeFactorDown * searchAOI.yPos;
            searchAOIResized.wdth = sizeFactorDown * searchAOI.wdth;
            searchAOIResized.hght = sizeFactorDown * searchAOI.hght;

            AOIProperties haarAOIResized;
            haarAOIResized.xPos = sizeFactorDown * haarAOI.xPos;
            haarAOIResized.yPos = sizeFactorDown * haarAOI.yPos;
            haarAOIResized.wdth = sizeFactorDown * haarAOI.wdth;
            haarAOIResized.hght = sizeFactorDown * haarAOI.hght;

            AOIProperties glintAOIResized;
            glintAOIResized.wdth = sizeFactorDown * mDetectionParameters.glintWdth;
            glintAOIResized.hght = glintAOIResized.wdth;

            // Haar-like feature detection

            std::vector<unsigned int> integralImage = calculateIntImg(imageResized, searchAOIResized);

            glintAOIResized = detectGlint(imageResized, searchAOIResized, glintAOIResized);
            glintAOIResized.xPos = searchAOIResized.xPos + glintAOIResized.xPos;
            glintAOIResized.yPos = searchAOIResized.yPos + glintAOIResized.yPos;

            mDataVariables.intensityInner.resize    (imgWdthResized * imgHghtResized);
            mDataVariables.intensityOuterLeft.resize(imgWdthResized * imgHghtResized);
            mDataVariables.intensityOuterRght.resize(imgWdthResized * imgHghtResized);

            haarAOIResized = detectPupilApprox(integralImage, searchAOIResized, haarAOIResized, glintAOIResized, mDataVariables.intensityInner, mDataVariables.intensityOuterLeft, mDataVariables.intensityOuterRght, imgWdthResized);

            // Upsample to original size

            glintAOI.xPos = sizeFactorUp * glintAOIResized.xPos;
            glintAOI.yPos = sizeFactorUp * glintAOIResized.yPos;
            glintAOI.wdth = sizeFactorUp * glintAOIResized.wdth;
            glintAOI.hght = sizeFactorUp * glintAOIResized.hght;

            searchAOI.xPos = sizeFactorUp * searchAOIResized.xPos;
            searchAOI.yPos = sizeFactorUp * searchAOIResized.yPos;
            searchAOI.wdth = sizeFactorUp * searchAOIResized.wdth;
            searchAOI.hght = sizeFactorUp * searchAOIResized.hght;

            haarAOI.xPos = sizeFactorUp * haarAOIResized.xPos;
            haarAOI.yPos = sizeFactorUp * haarAOIResized.yPos;
            haarAOI.wdth = sizeFactorUp * haarAOIResized.wdth;
            haarAOI.hght = sizeFactorUp * haarAOIResized.hght;

            // Offset Haar AOI

            haarAOI.xPos = searchAOI.xPos + haarAOI.xPos;
            haarAOI.yPos = searchAOI.yPos + haarAOI.yPos;
        }
        else // search not required
        {
            haarAOI.xPos  = searchAOI.xPos;
            haarAOI.yPos  = searchAOI.yPos;
            haarAOI.flag  = false;
            glintAOI.flag = false;
        }

        // Create new AOI for Canny edge deteciton

        double AOIXPos   = haarAOI.xPos + 0.5 * haarAOI.wdth; // centre of haar-like detector
        double AOIYPos   = haarAOI.yPos + 0.5 * haarAOI.hght;
        double AOIWdth   = mDetectionVariables.predictedWidth;
        double AOIHght   = mDetectionVariables.predictedHeight;

        if (AOIOffset > searchAOI.wdth || AOIOffset > searchAOI.hght)
        {
            if (searchAOI.wdth > searchAOI.hght) { AOIOffset = searchAOI.wdth; }
            else                                 { AOIOffset = searchAOI.hght; }
        }

        AOIProperties cannyAOI;
        cannyAOI.xPos = round(AOIXPos - 0.5 * AOIWdth - AOIOffset);
        cannyAOI.yPos = round(AOIYPos - 0.5 * AOIHght - AOIOffset);
        cannyAOI.wdth = round(AOIWdth + 2 * AOIOffset);
        cannyAOI.hght = round(AOIHght + 2 * AOIOffset);

        // Check limits

        if (cannyAOI.xPos < mAOI.xPos) { cannyAOI.xPos = mAOI.xPos; }
        if (cannyAOI.yPos < mAOI.yPos) { cannyAOI.yPos = mAOI.yPos; }
        if (cannyAOI.xPos + cannyAOI.wdth >= mAOI.xPos + mAOI.wdth) { cannyAOI.wdth = mAOI.xPos + mAOI.wdth - cannyAOI.xPos; }
        if (cannyAOI.yPos + cannyAOI.hght >= mAOI.yPos + mAOI.hght) { cannyAOI.hght = mAOI.yPos + mAOI.hght - cannyAOI.yPos; }

        // Crop image to new AOI

        cv::Rect outerRect(cannyAOI.xPos, cannyAOI.yPos, cannyAOI.wdth, cannyAOI.hght);
        cv::Mat imageAOIBGR = imageOriginalBGR(outerRect);

        // Convert back to grayscale

        cv::Mat imageAOIGray;
        cv::cvtColor(imageAOIBGR, imageAOIGray, cv::COLOR_BGR2GRAY);

        ///////////////////////////////////////////////////////////////////////
        /////////////////////// CANNY EDGE DETECTION  /////////////////////////
        ///////////////////////////////////////////////////////////////////////

        cv::Mat imageAOIGrayBlurred;
        int cannyBlurLevel = 2 * mDetectionParameters.cannyBlurLevel - 1; // should be odd
        if (cannyBlurLevel > 0) { cv::GaussianBlur(imageAOIGray, imageAOIGrayBlurred, cv::Size(cannyBlurLevel, cannyBlurLevel), 0, 0);
        } else                  { imageAOIGrayBlurred = imageAOIGray; }

        cv::Mat imageCannyEdges;
        cv::Canny(imageAOIGrayBlurred, imageCannyEdges, mDetectionParameters.cannyThresholdHigh, mDetectionParameters.cannyThresholdLow, 5);
        std::vector<int> cannyEdges          = cannyConversion(imageCannyEdges, cannyAOI); // convert to binary vector
        std::vector<int> cannyEdgesSharpened = sharpenEdges(cannyEdges, cannyAOI); // morphological operation
        std::vector<int> edgeIndices         = getEdgeIndices(cannyEdgesSharpened, 1); // used for drawing function

        /////////////////////////////////////////////////////////////////////////////
        //////////////////////////// EDGE SELECTION   ///////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        { // calculate predicted positions relative to Canny edge AOI
            mDetectionVariables.predictedXPosRelative = (1 - mDetectionVariables.certaintyPosition) * AOIXPos + mDetectionVariables.certaintyPosition * mDetectionVariables.predictedXPos - cannyAOI.xPos;
            mDetectionVariables.predictedYPosRelative = (1 - mDetectionVariables.certaintyPosition) * AOIYPos + mDetectionVariables.certaintyPosition * mDetectionVariables.predictedYPos - cannyAOI.yPos;
        }

        std::vector<edgeProperties> vEdgePropertiesAll = edgeSelection(mDetectionVariables, mDetectionParameters, cannyEdgesSharpened, cannyAOI);
        vEdgePropertiesAll = removeShortEdges(mDetectionParameters, vEdgePropertiesAll);

        /////////////////////////////////////////////////////////////////////////////
        //////////////////////////// EDGE SEGMENTATION   ////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        // Curvature calculation

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            int edgeSize = mEdgeProperties.pointIndices.size();

            std::vector<double> edgeXTangents(edgeSize, 0.0);
            std::vector<double> edgeYTangents(edgeSize, 0.0);

            calculateEdgeDirections(mEdgeProperties.pointIndices, edgeXTangents, edgeYTangents, cannyAOI);

            std::vector<double> edgeXNormals(edgeSize, 0.0);
            std::vector<double> edgeYNormals(edgeSize, 0.0);

            mEdgeProperties.curvatures = calculateCurvatures(mDetectionParameters, edgeXNormals, edgeYNormals, edgeXTangents, edgeYTangents);
            mEdgeProperties.xnormals   = edgeXNormals;
            mEdgeProperties.ynormals   = edgeYNormals;

            vEdgePropertiesAll[iEdge] = mEdgeProperties;
        }

        // Calculate curvature limits

        double curvatureUpperLimit;
        double curvatureLowerLimit;

        calculateCurvatureLimits(mDetectionVariables, mDetectionParameters, curvatureUpperLimit, curvatureLowerLimit);

        // Curvature segmentation

        { std::vector<edgeProperties> vEdgePropertiesNew;

            for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
            {
                std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationCurvature(mDetectionParameters, vEdgePropertiesAll[iEdge], curvatureUpperLimit, curvatureLowerLimit);
                vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
            }

            vEdgePropertiesAll = vEdgePropertiesNew;
        }

        vEdgePropertiesAll = removeShortEdges(mDetectionParameters, vEdgePropertiesAll);

        // Calculate additional edge properties

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            mEdgeProperties.length      = mEdgeProperties.pointIndices.size();
            mEdgeProperties.radii       = calculateEdgeRadii(mEdgeProperties, cannyAOI, mDetectionVariables.predictedXPosRelative, mDetectionVariables.predictedYPosRelative);
            mEdgeProperties.gradients   = calculateRadialGradients(mDetectionVariables, mDetectionParameters, imageAOIGray, mEdgeProperties.pointIndices);
            mEdgeProperties.intensities = findEdgeIntensities(imageAOIGray, mDetectionParameters, mEdgeProperties, cannyAOI);

            vEdgePropertiesAll[iEdge] = mEdgeProperties;
        }

        // Length segmentation

        if (!(mAdvancedOptions.CURVATURE_MEASUREMENT))
        {
            std::vector<edgeProperties> vEdgePropertiesNew;

            for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
            {
                std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationLength(mDetectionVariables, mDetectionParameters, vEdgePropertiesAll[iEdge]);
                vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
            }

            vEdgePropertiesAll = vEdgePropertiesNew;
        }

        vEdgePropertiesAll = removeShortEdges(mDetectionParameters, vEdgePropertiesAll);

        ////////////////////////////////////////////////////////////////////////
        //////////////////////// EDGE TERMINAL FILTER   ////////////////////////
        ////////////////////////////////////////////////////////////////////////

        if (mDetectionVariables.certaintyPosition > certaintyThreshold || mDetectionVariables.certaintyFeatures > certaintyThreshold || mDetectionVariables.certaintyAverages > certaintyThreshold)
        {
            std::vector<edgeProperties> vEdgePropertiesNew;

            for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
            {
                edgeProperties mEdgeProperties    = vEdgePropertiesAll[iEdge];

                std::vector<edgeProperties> vEdgePropertiesTemp = edgeTerminalFilter(mDetectionVariables, mDetectionParameters, mEdgeProperties, mDetectionParameters.thresholdScoreDiffEdge);
                vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
            }

            vEdgePropertiesAll = vEdgePropertiesNew;
        }

        vEdgePropertiesAll = removeShortEdges(mDetectionParameters, vEdgePropertiesAll);

        ///////////////////////////////////////////////////////////////////////////
        ////////////////////////// EDGE CLASSIFICATION  ///////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        // Calculate some edge properties

        { std::vector<edgeProperties> vEdgePropertiesNew;

            for (int iEdge = 0, jEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
            {
                edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

                // re-evaluate edge properties

                mEdgeProperties.length    = mEdgeProperties.pointIndices.size();

                mEdgeProperties.intensity = calculateMeanInt(mEdgeProperties.intensities);
                mEdgeProperties.gradient  = calculateMeanInt(mEdgeProperties.gradients);

                mEdgeProperties.radius    = calculateMean(mEdgeProperties.radii);
                mEdgeProperties.radiusVar = calculateVariance(mEdgeProperties.radii) / mEdgeProperties.radius; // relative variance

                calculateCurvatureStats(mDetectionParameters, mEdgeProperties);

                mEdgeProperties.index     = iEdge;
                mEdgeProperties.tag       = 0;

                restoreEdgePoints(mEdgeProperties, cannyEdgesSharpened, cannyAOI); // Restore some points
                mEdgeProperties.size      = mEdgeProperties.pointIndices.size();

                vEdgePropertiesNew.push_back(mEdgeProperties);

                jEdge++;
            }

            vEdgePropertiesAll = vEdgePropertiesNew;
        }

        // Do edge classification

        std::vector<int> acceptedEdges = edgeClassification(mDetectionVariables, mDetectionParameters, vEdgePropertiesAll);

        std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, numEdges = acceptedEdges.size(); iEdge < numEdges; iEdge++) // grab accepted edges
        {
            int jEdge = acceptedEdges[iEdge];
            vEdgePropertiesAll[jEdge].tag = 1; // new tag
            vEdgePropertiesNew.push_back(vEdgePropertiesAll[jEdge]);
        }

        //////////////////////////////////////////////////////////////////
        /////////////////////// ELLIPSE FITTING  /////////////////////////
        //////////////////////////////////////////////////////////////////

        std::vector<ellipseProperties> vEllipsePropertiesAll = getEllipseFits(mDetectionVariables, mDetectionParameters, vEdgePropertiesNew, cannyAOI); // ellipse fitting

        ellipseProperties mEllipseProperties; // properties of accepted fit

        int numFits = vEllipsePropertiesAll.size();

        if (numFits > 0)
        {
            std::vector<int> acceptedFitIndices;

            if (numFits > 1) { acceptedFitIndices = ellipseFitFilter(mDetectionVariables, mDetectionParameters, vEllipsePropertiesAll); } // grab best fit
            else             { acceptedFitIndices = {0}; }

            int numFitsAccepted = acceptedFitIndices.size();

            std::vector<ellipseProperties> vEllipsePropertiesNew(numFitsAccepted);

            // Tag fits

            for (int iFit = 0; iFit < numFitsAccepted; iFit++)
            {
                ellipseProperties mEllipsePropertiesTemp = vEllipsePropertiesAll[acceptedFitIndices[iFit]];
                vEllipsePropertiesNew[iFit] = mEllipsePropertiesTemp;
                vEllipsePropertiesAll[acceptedFitIndices[iFit]].tag = 1;

                // Tag edges

                for (int iEdge = 0, numEdges = mEllipsePropertiesTemp.edgeIndices.size(); iEdge < numEdges; iEdge++)
                {
                    int jEdge = mEllipsePropertiesTemp.edgeIndices[iEdge];
                    vEdgePropertiesAll[jEdge].tag = 2;
                }
            }

            mEllipseProperties = vEllipsePropertiesNew[0]; // set equal to highest score fit
            mEllipseProperties.DETECTED = true;

            // Get average properties if multiple fits are accepted

            if (numFitsAccepted > 1) { ellipseFitAverage(mEllipseProperties, vEllipsePropertiesNew); }

            // Calculate average properties of fitted edges

            std::vector<double> curvatures;
            std::vector<int> intensities;
            std::vector<int> gradients;
            std::vector<int> numEdgePoints_1;
            std::vector<int> numEdgePoints_2;

            for (int iEdge = 0, numEdgesAll = vEdgePropertiesAll.size(); iEdge < numEdgesAll; iEdge++)
            {
                edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

                if (mEdgeProperties.tag == 2)
                {
                    intensities.push_back(mEdgeProperties.intensity);
                    gradients.push_back(mEdgeProperties.gradient);
                    numEdgePoints_1.push_back(mEdgeProperties.size);

                    if (mEdgeProperties.curvature < 180) // ignore short edges for which no curvature information is available
                    {
                        numEdgePoints_2.push_back(mEdgeProperties.size);
                        curvatures.push_back(mEdgeProperties.curvature);
                    }
                }
            }

            // Calculate intensity and gradient means

            double intensityMean = 0;
            double gradientMean  = 0;
            double numEdgePointsTotal_1 = std::accumulate(numEdgePoints_1.begin(), numEdgePoints_1.end(), 0.0);

            if (numEdgePointsTotal_1 > 0)
            {
                for (int iEdge = 0, numEdges = intensities.size(); iEdge < numEdges; iEdge++)
                {
                    double weight  = numEdgePoints_1[iEdge] / numEdgePointsTotal_1;
                    intensityMean += weight * intensities[iEdge];
                    gradientMean  += weight * gradients[iEdge];
                }
            }
            else
            {
                intensityMean = mDetectionVariables.predictedIntensity;
                gradientMean  = mDetectionVariables.predictedGradient;
            }

            mEllipseProperties.intensity = intensityMean;
            mEllipseProperties.gradient  =  gradientMean;

            // Calculate curvature mean

            double curvatureMean = 0;
            double numEdgePointsTotal_2 = std::accumulate(numEdgePoints_2.begin(), numEdgePoints_2.end(), 0.0);

            if (numEdgePointsTotal_2 > 0)
            {
                for (int iEdge = 0, numEdges = curvatures.size(); iEdge < numEdges; iEdge++)
                {
                    double weight = numEdgePoints_2[iEdge] / numEdgePointsTotal_2;
                    curvatureMean += weight * curvatures[iEdge];
                }
            }
            else { curvatureMean = mDetectionVariables.predictedCurvature; }

            mEllipseProperties.curvature = curvatureMean;
        }
        else
        {
            mEllipseProperties.tag      = -1;
            mEllipseProperties.DETECTED = false;
            vEllipsePropertiesAll.push_back(mEllipseProperties);
        }


        /////////////////////////////////////////////////////////////////
        /////////////////////// SAVING DATA  ////////////////////////////
        /////////////////////////////////////////////////////////////////

        detectionVariables mDetectionVariablesNew = mDetectionVariables; // properties for next frame

        mDataVariables.edgeData    = vEdgePropertiesAll;     // edge data
        mDataVariables.ellipseData = vEllipsePropertiesAll;  // ellipse data

        // Save parameters

        mDataVariables.DETECTED = mEllipseProperties.DETECTED;

        // Calculate new threshold limits. Thresholds are harsher with higher certainties (= lower certainty factors)

        int AOISize;
        if (mAOI.wdth > mAOI.hght) { AOISize = mAOI.wdth; }
        else                       { AOISize = mAOI.hght; }

        double maxChangeThresholdAspectRatio   = 1.0 - mDetectionParameters.aspectRatioMin;
        double maxChangeThresholdCircumference = std::abs(mDetectionParameters.circumferenceMax - mDetectionParameters.circumferenceMin) / mDetectionParameters.circumferenceMax;
        double maxChangeThresholdPosition      = AOISize;

        double rangeChangeThresholdAspectRatio   = maxChangeThresholdAspectRatio   - mDetectionParameters.thresholdChangeAspectRatio;
        double rangeChangeThresholdCircumference = maxChangeThresholdCircumference - mDetectionParameters.thresholdChangeCircumference;
        double rangeChangeThresholdPosition      = maxChangeThresholdPosition      - mDetectionParameters.thresholdChangePosition;
        double rangeOffsetCircumference          = maxChangeThresholdCircumference - mDetectionParameters.circumferenceOffset;

        mDetectionVariablesNew.thresholdChangeAspectRatio   = rangeChangeThresholdAspectRatio   * (1 - mDetectionVariables.certaintyFeatures) + mDetectionParameters.thresholdChangeAspectRatio;
        mDetectionVariablesNew.thresholdChangeCircumference = rangeChangeThresholdCircumference * (1 - mDetectionVariables.certaintyFeatures) + mDetectionParameters.thresholdChangeCircumference;
        mDetectionVariablesNew.thresholdChangePosition      = rangeChangeThresholdPosition      * (1 - mDetectionVariables.certaintyPosition) + mDetectionParameters.thresholdChangePosition;
        mDetectionVariablesNew.offsetCircumference          = rangeOffsetCircumference          * (1 - mDetectionVariables.certaintyAverages) + mDetectionParameters.circumferenceOffset;

        std::vector<double> certaintyFactors = {mDetectionVariables.certaintyAverages, mDetectionVariables.certaintyFeatures, mDetectionVariables.certaintyAverages};
        double certaintyFactorMin = *std::min_element(certaintyFactors.begin(), certaintyFactors.end());
        mDetectionVariablesNew.thresholdScore = (1 - certaintyFactorMin) * mDetectionParameters.thresholdScore;

        if (!mEllipseProperties.DETECTED) // pupil not detected
        {
            // Running averages

            double meanAspectRatio   = initialAspectRatio;
            double meanCircumference = 0.5 * (mDetectionParameters.circumferenceMax + mDetectionParameters.circumferenceMin);
            double meanWidth         = meanCircumference / M_PI;
            double meanHeight        = meanCircumference / M_PI;
            double meanIntensity     = initialIntensity;
            double meanGradient      = 0;

            double curvatureMax = getCurvatureUpperLimit(meanCircumference, meanAspectRatio, mDetectionParameters.windowLengthEdge);
            double curvatureMin = getCurvatureLowerLimit(meanCircumference, meanAspectRatio, mDetectionParameters.windowLengthEdge);

            double meanCurvature = (0.5 * (curvatureMax + curvatureMin));

            // Momentum terms should decay to zero. Assume that pupil was stationary (i.e. error = 0)

            mDetectionVariablesNew.momentumAspectRatio   = mDetectionVariables.momentumAspectRatio   * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumCircumference = mDetectionVariables.momentumCircumference * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumWidth         = mDetectionVariables.momentumWidth         * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumHeight        = mDetectionVariables.momentumHeight        * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumCurvature     = mDetectionVariables.momentumCurvature     * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumIntensity     = mDetectionVariables.momentumIntensity     * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumGradient      = mDetectionVariables.momentumGradient      * (1 - mDetectionParameters.alphaFeatures);
            mDetectionVariablesNew.momentumXPos          = mDetectionVariables.momentumXPos          * (1 - mDetectionParameters.alphaPosition);
            mDetectionVariablesNew.momentumYPos          = mDetectionVariables.momentumYPos          * (1 - mDetectionParameters.alphaPosition);

            // Averages should decay to initial values

            mDetectionVariablesNew.averageAspectRatio   = mDetectionVariables.averageAspectRatio   + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanAspectRatio   - mDetectionVariables.averageAspectRatio);
            mDetectionVariablesNew.averageCircumference = mDetectionVariables.averageCircumference + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanCircumference - mDetectionVariables.averageCircumference);
            mDetectionVariablesNew.averageWidth         = mDetectionVariables.averageWidth         + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanWidth         - mDetectionVariables.averageWidth);
            mDetectionVariablesNew.averageHeight        = mDetectionVariables.averageHeight        + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanHeight        - mDetectionVariables.averageHeight);
            mDetectionVariablesNew.averageCurvature     = mDetectionVariables.averageCurvature     + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanCurvature     - mDetectionVariables.averageCurvature);
            mDetectionVariablesNew.averageIntensity     = mDetectionVariables.averageIntensity     + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanIntensity     - mDetectionVariables.averageIntensity);
            mDetectionVariablesNew.averageGradient      = mDetectionVariables.averageGradient      + mDetectionParameters.alphaAverages * (1 - mDetectionVariables.certaintyAverages) * (meanGradient      - mDetectionVariables.averageGradient);

            // Feature predictions should decay to average terms. Certainty term gives some latency.

            mDetectionVariablesNew.predictedAspectRatio   = mDetectionVariables.predictedAspectRatio   + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageAspectRatio   - mDetectionVariables.predictedAspectRatio)   + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumAspectRatio;
            mDetectionVariablesNew.predictedCircumference = mDetectionVariables.predictedCircumference + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageCircumference - mDetectionVariables.predictedCircumference) + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumCircumference;
            mDetectionVariablesNew.predictedWidth         = mDetectionVariables.predictedWidth         + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageWidth         - mDetectionVariables.predictedWidth)         + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumWidth;
            mDetectionVariablesNew.predictedHeight        = mDetectionVariables.predictedHeight        + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageHeight        - mDetectionVariables.predictedHeight)        + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumHeight;
            mDetectionVariablesNew.predictedCurvature     = mDetectionVariables.predictedCurvature     + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageCurvature     - mDetectionVariables.predictedCurvature)     + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumCurvature;
            mDetectionVariablesNew.predictedIntensity     = mDetectionVariables.predictedIntensity     + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageIntensity     - mDetectionVariables.predictedIntensity)     + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumIntensity;
            mDetectionVariablesNew.predictedGradient      = mDetectionVariables.predictedGradient      + mDetectionParameters.alphaFeatures * (1 - mDetectionVariables.certaintyFeatures) * (mDetectionVariables.averageGradient      - mDetectionVariables.predictedGradient)      + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumGradient;

            mDetectionVariablesNew.predictedAngle = mDetectionVariables.predictedAngle * (1 - mDetectionParameters.alphaFeatures * mDetectionVariables.certaintyFeatures); // reduce to zero

            // Position predictions should decay to approximate detection position. Certainty term gives some latency.

            mDetectionVariablesNew.predictedXPos = mDetectionVariables.predictedXPos + mDetectionParameters.alphaPosition * (1 - mDetectionVariables.certaintyPosition) * (haarAOI.xPos + 0.5 * haarAOI.wdth - mDetectionVariables.predictedXPos) + mDetectionVariables.certaintyPosition * mDetectionVariables.momentumXPos;
            mDetectionVariablesNew.predictedYPos = mDetectionVariables.predictedYPos + mDetectionParameters.alphaPosition * (1 - mDetectionVariables.certaintyPosition) * (haarAOI.yPos + 0.5 * haarAOI.hght - mDetectionVariables.predictedYPos) + mDetectionVariables.certaintyPosition * mDetectionVariables.momentumYPos;

            // Certainty decays to minimum value

            mDetectionVariablesNew.certaintyPositionPrime = mDetectionVariables.certaintyPositionPrime - mDetectionParameters.alphaCertainty * mDetectionParameters.alphaPosition;
            mDetectionVariablesNew.certaintyFeaturesPrime = mDetectionVariables.certaintyFeaturesPrime - mDetectionParameters.alphaCertainty * mDetectionParameters.alphaFeatures;
            mDetectionVariablesNew.certaintyAveragesPrime = mDetectionVariables.certaintyAveragesPrime - mDetectionParameters.alphaCertainty * mDetectionParameters.alphaAverages;
        }
        else // pupil detected
        {
            // Exact values

            mDataVariables.exactAspectRatio   = mEllipseProperties.aspectRatio;
            mDataVariables.exactCircumference = mEllipseProperties.circumference;

            mDataVariables.exactXPos = mEllipseProperties.xPos;
            mDataVariables.exactYPos = mEllipseProperties.yPos;

            // Delta error

            double errorAspectRatio   = mEllipseProperties.aspectRatio   - mDetectionVariables.predictedAspectRatio;
            double errorCircumference = mEllipseProperties.circumference - mDetectionVariables.predictedCircumference;
            double errorWidth         = mEllipseProperties.width         - mDetectionVariables.predictedWidth;
            double errorHeight        = mEllipseProperties.height        - mDetectionVariables.predictedHeight;
            double errorCurvature     = mEllipseProperties.curvature     - mDetectionVariables.predictedCurvature;
            double errorIntensity     = mEllipseProperties.intensity     - mDetectionVariables.predictedIntensity;
            double errorGradient      = mEllipseProperties.gradient      - mDetectionVariables.predictedGradient;
            double errorAngle         = mEllipseProperties.angle         - mDetectionVariables.predictedAngle;
            double errorXPosition     = mDataVariables.exactXPos         - mDetectionVariables.predictedXPos;
            double errorYPosition     = mDataVariables.exactYPos         - mDetectionVariables.predictedYPos;

            // Momentum. Rate of change approximation.

            mDetectionVariablesNew.momentumAspectRatio   =  mDetectionVariables.momentumAspectRatio   + mDetectionParameters.alphaFeatures * (errorAspectRatio   - mDetectionVariables.momentumAspectRatio);
            mDetectionVariablesNew.momentumCircumference =  mDetectionVariables.momentumCircumference + mDetectionParameters.alphaFeatures * (errorCircumference - mDetectionVariables.momentumCircumference);
            mDetectionVariablesNew.momentumWidth         =  mDetectionVariables.momentumWidth         + mDetectionParameters.alphaFeatures * (errorWidth         - mDetectionVariables.momentumWidth);
            mDetectionVariablesNew.momentumHeight        =  mDetectionVariables.momentumHeight        + mDetectionParameters.alphaFeatures * (errorHeight        - mDetectionVariables.momentumHeight);
            mDetectionVariablesNew.momentumGradient      =  mDetectionVariables.momentumGradient      + mDetectionParameters.alphaFeatures * (errorGradient      - mDetectionVariables.momentumGradient);
            mDetectionVariablesNew.momentumIntensity     =  mDetectionVariables.momentumIntensity     + mDetectionParameters.alphaFeatures * (errorIntensity     - mDetectionVariables.momentumIntensity);
            mDetectionVariablesNew.momentumCurvature     =  mDetectionVariables.momentumCurvature     + mDetectionParameters.alphaFeatures * (errorCurvature     - mDetectionVariables.momentumCurvature);
            mDetectionVariablesNew.momentumXPos          =  mDetectionVariables.momentumXPos          + mDetectionParameters.alphaPosition * (errorXPosition     - mDetectionVariables.momentumXPos);
            mDetectionVariablesNew.momentumYPos          =  mDetectionVariables.momentumYPos          + mDetectionParameters.alphaPosition * (errorYPosition     - mDetectionVariables.momentumYPos);

            // Update predictions. Only add momentum when certainty is high

            mDetectionVariablesNew.predictedAspectRatio   = mDetectionVariables.predictedAspectRatio   + mDetectionParameters.alphaFeatures * errorAspectRatio   + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumAspectRatio;
            mDetectionVariablesNew.predictedCircumference = mDetectionVariables.predictedCircumference + mDetectionParameters.alphaFeatures * errorCircumference + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumCircumference;
            mDetectionVariablesNew.predictedWidth         = mDetectionVariables.predictedWidth         + mDetectionParameters.alphaFeatures * errorWidth         + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumWidth;
            mDetectionVariablesNew.predictedHeight        = mDetectionVariables.predictedHeight        + mDetectionParameters.alphaFeatures * errorHeight        + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumHeight;
            mDetectionVariablesNew.predictedCurvature     = mDetectionVariables.predictedCurvature     + mDetectionParameters.alphaFeatures * errorCurvature     + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumCurvature;
            mDetectionVariablesNew.predictedIntensity     = mDetectionVariables.predictedIntensity     + mDetectionParameters.alphaFeatures * errorIntensity     + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumIntensity;
            mDetectionVariablesNew.predictedGradient      = mDetectionVariables.predictedGradient      + mDetectionParameters.alphaFeatures * errorGradient      + mDetectionVariables.certaintyFeatures * mDetectionVariables.momentumGradient;
            mDetectionVariablesNew.predictedAngle         = mDetectionVariables.predictedAngle         + mDetectionParameters.alphaFeatures * errorAngle;

            mDetectionVariablesNew.predictedXPos          = mDetectionVariables.predictedXPos          + mDetectionParameters.alphaPosition * errorXPosition     + mDetectionVariables.certaintyPosition * mDetectionVariables.momentumXPos;
            mDetectionVariablesNew.predictedYPos          = mDetectionVariables.predictedYPos          + mDetectionParameters.alphaPosition * errorYPosition     + mDetectionVariables.certaintyPosition * mDetectionVariables.momentumYPos;

            // Averages

            mDetectionVariablesNew.averageAspectRatio   = mDetectionVariables.averageAspectRatio   + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedAspectRatio   - mDetectionVariables.averageAspectRatio);
            mDetectionVariablesNew.averageCircumference = mDetectionVariables.averageCircumference + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedCircumference - mDetectionVariables.averageCircumference);
            mDetectionVariablesNew.averageWidth         = mDetectionVariables.averageWidth         + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedWidth         - mDetectionVariables.averageWidth);
            mDetectionVariablesNew.averageHeight        = mDetectionVariables.averageHeight        + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedHeight        - mDetectionVariables.averageHeight);
            mDetectionVariablesNew.averageCurvature     = mDetectionVariables.averageCurvature     + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedCurvature     - mDetectionVariables.averageCurvature);
            mDetectionVariablesNew.averageIntensity     = mDetectionVariables.averageIntensity     + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedIntensity     - mDetectionVariables.averageIntensity);
            mDetectionVariablesNew.averageGradient      = mDetectionVariables.averageGradient      + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedGradient      - mDetectionVariables.averageGradient);

            // Determine certainty of current measurement

            double displacement                = sqrt(errorXPosition * errorXPosition + errorYPosition * errorYPosition);
            double relativeChangeCircumference = std::abs(errorCircumference) / std::max(mEllipseProperties.circumference, mDetectionVariables.predictedCircumference);

            double certaintyPosition      = calculateCertainty(displacement,                mDetectionParameters.thresholdChangePosition);
            double certaintyAspectRatio   = calculateCertainty(std::abs(errorAspectRatio),  mDetectionParameters.thresholdChangeAspectRatio);
            double certaintyCircumference = calculateCertainty(relativeChangeCircumference, mDetectionParameters.thresholdChangeCircumference);
            double certaintyFeatures      = 0.5 * (certaintyAspectRatio + certaintyCircumference);

            // Increase or reduce certainty with a fraction of a constant step size (product of the two learning rate parameters)

            mDetectionVariablesNew.certaintyPositionPrime = mDetectionVariables.certaintyPositionPrime + certaintyPosition * mDetectionParameters.alphaCertainty * mDetectionParameters.alphaPosition;
            mDetectionVariablesNew.certaintyFeaturesPrime = mDetectionVariables.certaintyFeaturesPrime + certaintyFeatures * mDetectionParameters.alphaCertainty * mDetectionParameters.alphaFeatures;
            mDetectionVariablesNew.certaintyAveragesPrime = mDetectionVariables.certaintyAveragesPrime + certaintyFeatures * mDetectionParameters.alphaCertainty * mDetectionParameters.alphaAverages;
        }

        // Logistic functions are used to add response latency to changes in certainty.
        // Multiple detections are required until a reduction in certainty leads to changes in parameter limits or variable predictions

        mDetectionVariablesNew.certaintyPosition = 1 / (1 + exp(-certaintyLatency * (mDetectionVariablesNew.certaintyPositionPrime - 0.5)));
        mDetectionVariablesNew.certaintyFeatures = 1 / (1 + exp(-certaintyLatency * (mDetectionVariablesNew.certaintyFeaturesPrime - 0.5)));
        mDetectionVariablesNew.certaintyAverages = 1 / (1 + exp(-certaintyLatency * (mDetectionVariablesNew.certaintyAveragesPrime - 0.5)));

        // For drawing

        mDrawVariables.DETECTED = mEllipseProperties.DETECTED;
        mDrawVariables.haarAOI  = haarAOI;
        mDrawVariables.glintAOI = glintAOI;
        mDrawVariables.cannyAOI = cannyAOI;

        mDrawVariables.exactXPos = round(mDataVariables.exactXPos);
        mDrawVariables.exactYPos = round(mDataVariables.exactYPos);

        mDrawVariables.predictedXPos = round(mDetectionVariables.predictedXPos);
        mDrawVariables.predictedYPos = round(mDetectionVariables.predictedYPos);

        mDrawVariables.cannyEdgeIndices    = edgeIndices;
        mDrawVariables.edgeData            = vEdgePropertiesAll;
        mDrawVariables.ellipseCoefficients = mEllipseProperties.coefficients;

        return mDetectionVariablesNew; // use these variables for next frame
    }

    double flashDetection(const cv::Mat& imgBGR)
    {
        int imgSize = imgBGR.cols * imgBGR.rows;

        if (imgSize > 0)
        {
            unsigned long long intensityTotal = 0;
            cv::Mat img;
            cv::cvtColor(imgBGR, img, cv::COLOR_BGR2GRAY);
            uchar *ptr = img.data;
            for (int iPixel = 0; iPixel < imgSize; iPixel++) { intensityTotal += ptr[iPixel]; }
            return (intensityTotal / (double) imgSize);
        } else { return 0; }
    }

    // Look-up tables


    double getCurvatureUpperLimit(double circumference, double aspectRatio, int windowLength)
    {
        static const std::vector<double> arrayCircumferences = {60,65.918,71.837,77.755,83.673,89.592,95.51,101.43,107.35,113.27,119.18,125.1,131.02,136.94,142.86,148.78,154.69,160.61,166.53,172.45,178.37,184.29,190.2,196.12,202.04,207.96,213.88,219.8,225.71,231.63,237.55,243.47,249.39,255.31,261.22,267.14,273.06,278.98,284.9,290.82,296.73,302.65,308.57,314.49,320.41,326.33,332.24,338.16,344.08,350};
        static const std::vector<double> arrayAspectRatios   = {0.4,0.41224,0.42449,0.43673,0.44898,0.46122,0.47347,0.48571,0.49796,0.5102,0.52245,0.53469,0.54694,0.55918,0.57143,0.58367,0.59592,0.60816,0.62041,0.63265,0.6449,0.65714,0.66939,0.68163,0.69388,0.70612,0.71837,0.73061,0.74286,0.7551,0.76735,0.77959,0.79184,0.80408,0.81633,0.82857,0.84082,0.85306,0.86531,0.87755,0.8898,0.90204,0.91429,0.92653,0.93878,0.95102,0.96327,0.97551,0.98776,1};

        int x = findUpperBound(arrayCircumferences, circumference);
        int y = findUpperBound(arrayAspectRatios,   aspectRatio);

        int arrayWidth = arrayCircumferences.size();

        switch(windowLength)
        {

        case 5:
        {
            static const std::vector<double> arrayCurvatureMax =    {97.478,95.166,92.838,90.504,88.175,85.862,83.575,81.325,79.123,76.98 ,74.906,72.908,70.995,69.173,67.445,65.816,64.285,62.85 ,61.509,60.256,59.084,57.988,56.959,55.99 ,55.076,54.211,53.39 ,52.61 ,51.869,51.163,50.492,49.855,49.251,48.679,48.138,47.629,47.149,46.699,46.277,45.882,45.513,45.169,44.848,44.55 ,44.253,43.49 ,40.985,40.312,40.091,39.902,
                                                                     96.009,93.664,91.305,88.941,86.583,84.241,81.926,79.647,77.416,75.241,73.134,71.102,69.153,67.295,65.534,63.874,62.319,60.869,59.525,58.283,57.139,56.087,55.119,54.228,53.405,52.64 ,51.928,51.26 ,50.63 ,50.035,49.47 ,48.934,48.423,47.937,47.475,47.036,46.62 ,46.226,45.855,45.505,45.176,44.867,44.578,44.303,43.924,41.978,40.418,40.138,39.947,39.772,
                                                                     94.566,92.197,89.816,87.434,85.06 ,82.704,80.376,78.085,75.841,73.654,71.531,69.481,67.512,65.631,63.844,62.156,60.571,59.093,57.722,56.461,55.306,54.255,53.302,52.442,51.666,50.965,50.33 ,49.752,49.221,48.73 ,48.273,47.842,47.435,47.047,46.678,46.324,45.986,45.663,45.355,45.062,44.783,44.518,44.267,44.003,43.091,40.674,40.156,39.962,39.792,39.633,
                                                                     93.137,90.751,88.357,85.965,83.585,81.226,78.897,76.606,74.364,72.179,70.057,68.007,66.035,64.147,62.349,60.646,59.042,57.54 ,56.143,54.853,53.67 ,52.594,51.623,50.752,49.978,49.293,48.69 ,48.16 ,47.692,47.278,46.908,46.573,46.266,45.98 ,45.711,45.454,45.208,44.97 ,44.74 ,44.517,44.302,44.094,43.888,43.513,41.379,40.138,39.92 ,39.762,39.615,39.476,
                                                                     91.715,89.318,86.917,84.523,82.144,79.789,77.468,75.189,72.959,70.787,68.679,66.642,64.682,62.803,61.011,59.309,57.701,56.19 ,54.778,53.468,52.259,51.154,50.151,49.25 ,48.448,47.743,47.129,46.6  ,46.148,45.766,45.444,45.171,44.939,44.739,44.561,44.401,44.25 ,44.107,43.967,43.83 ,43.693,43.555,43.378,42.329,40.129,39.773,39.639,39.517,39.4  ,39.288,
                                                                     90.297,87.893,85.491,83.099,80.727,78.384,76.077,73.816,71.606,69.456,67.372,65.358,63.421,61.564,59.79 ,58.104,56.508,55.003,53.591,52.275,51.054,49.929,48.9  ,47.969,47.132,46.391,45.742,45.182,44.707,44.313,43.991,43.734,43.534,43.381,43.266,43.179,43.112,43.058,43.01 ,42.964,42.916,42.856,42.538,40.343,39.457,39.351,39.277,39.202,39.125,39.049,
                                                                     88.882,86.475,84.076,81.691,79.331,77.004,74.717,72.479,70.295,68.174,66.119,64.136,62.229,60.402,58.658,56.998,55.425,53.939,52.542,51.234,50.015,48.886,47.845,46.894,46.032,45.258,44.571,43.971,43.455,43.02 ,42.665,42.383,42.169,42.016,41.917,41.863,41.844,41.851,41.877,41.912,41.949,41.925,40.781,38.944,38.8  ,38.8  ,38.794,38.78 ,38.758,38.73 ,
                                                                     87.468,85.063,82.67 ,80.297,77.952,75.645,73.382,71.171,69.019,66.93 ,64.909,62.962,61.092,59.301,57.592,55.966,54.424,52.967,51.594,50.307,49.103,47.983,46.945,45.99 ,45.116,44.322,43.608,42.974,42.418,41.939,41.536,41.207,40.95 ,40.761,40.635,40.567,40.549,40.574,40.633,40.717,40.805,40.544,38.414,37.911,37.986,38.074,38.149,38.211,38.259,38.294,
                                                                     86.058,83.658,81.274,78.915,76.59 ,74.306,72.071,69.891,67.772,65.719,63.737,61.829,59.999,58.248,56.578,54.991,53.486,52.064,50.723,49.464,48.284,47.184,46.161,45.214,44.342,43.543,42.817,42.162,41.578,41.064,40.619,40.244,39.936,39.694,39.517,39.403,39.347,39.345,39.391,39.475,39.507,38.199,36.758,36.82 ,36.989,37.157,37.315,37.459,37.588,37.7  ,
                                                                     84.651,82.259,79.888,77.548,75.245,72.987,70.782,68.635,66.553,64.538,62.597,60.731,58.943,57.235,55.608,54.062,52.598,51.214,49.91 ,48.684,47.535,46.461,45.461,44.532,43.672,42.881,42.155,41.494,40.897,40.362,39.89 ,39.478,39.127,38.837,38.606,38.435,38.323,38.267,38.265,38.294,37.874,35.743,35.467,35.644,35.855,36.075,36.298,36.516,36.723,36.916,
                                                                     83.249,80.868,78.514,76.194,73.916,71.688,69.515,67.404,65.36 ,63.386,61.487,59.664,57.921,56.257,54.674,53.172,51.75 ,50.407,49.142,47.953,46.838,45.796,44.823,43.918,43.079,42.303,41.588,40.932,40.334,39.792,39.306,38.873,38.494,38.167,37.894,37.674,37.506,37.391,37.325,37.184,35.498,34.26 ,34.316,34.483,34.686,34.916,35.165,35.424,35.687,35.946,
                                                                     81.854,79.487,77.152,74.855,72.605,70.408,68.271,66.198,64.193,62.262,60.406,58.628,56.93 ,55.312,53.774,52.316,50.937,49.637,48.412,47.261,46.183,45.174,44.233,43.356,42.541,41.786,41.088,40.446,39.856,39.317,38.828,38.387,37.994,37.646,37.345,37.089,36.879,36.714,36.567,35.795,33.568,33.275,33.333,33.447,33.603,33.798,34.026,34.281,34.557,34.846,
                                                                     80.467,78.117,75.803,73.533,71.313,69.15 ,67.049,65.015,63.053,61.165,59.354,57.622,55.969,54.397,52.905,51.492,50.157,48.898,47.715,46.603,45.562,44.588,43.679,42.832,42.045,41.314,40.638,40.013,39.438,38.909,38.426,37.987,37.59 ,37.233,36.917,36.64 ,36.402,36.198,35.851,33.769,32.638,32.55 ,32.553,32.601,32.692,32.824,32.997,33.207,33.451,33.723,
                                                                     79.089,76.76 ,74.47 ,72.229,70.041,67.913,65.851,63.858,61.939,60.095,58.33 ,56.643,55.037,53.511,52.064,50.696,49.405,48.189,47.047,45.975,44.971,44.032,43.157,42.341,41.582,40.877,40.224,39.62 ,39.063,38.549,38.078,37.646,37.253,36.897,36.576,36.29 ,36.036,35.774,34.658,32.448,32.107,32.015,31.967,31.956,31.983,32.048,32.152,32.296,32.477,32.696,
                                                                     77.723,75.416,73.154,70.943,68.789,66.699,64.677,62.726,60.851,59.052,57.333,55.693,54.133,52.654,51.252,49.929,48.681,47.508,46.406,45.373,44.406,43.503,42.66 ,41.875,41.146,40.468,39.839,39.257,38.72 ,38.224,37.767,37.347,36.963,36.613,36.294,36.006,35.739,35.243,32.9  ,31.921,31.742,31.625,31.539,31.483,31.459,31.466,31.508,31.584,31.696,31.845,
                                                                     76.369,74.088,71.854,69.676,67.559,65.508,63.528,61.62 ,59.79 ,58.037,56.364,54.771,53.258,51.824,50.468,49.189,47.984,46.852,45.79 ,44.795,43.865,42.997,42.187,41.433,40.732,40.081,39.478,38.919,38.402,37.924,37.484,37.079,36.706,36.365,36.054,35.768,35.449,34.026,31.974,31.637,31.473,31.339,31.228,31.143,31.082,31.047,31.039,31.06 ,31.112,31.195,
                                                                     75.03 ,72.776,70.574,68.431,66.352,64.342,62.404,60.541,58.755,57.049,55.423,53.876,52.409,51.021,49.71 ,48.474,47.312,46.221,45.198,44.241,43.347,42.512,41.735,41.011,40.339,39.715,39.136,38.6  ,38.104,37.645,37.222,36.832,36.474,36.144,35.842,35.552,34.922,32.448,31.641,31.432,31.267,31.123,30.999,30.895,30.811,30.748,30.707,30.688,30.693,30.723,
                                                                     73.707,71.482,69.314,67.208,65.168,63.2  ,61.305,59.488,57.748,56.088,54.509,53.009,51.588,50.245,48.978,47.786,46.665,45.614,44.63 ,43.709,42.85 ,42.049,41.303,40.609,39.964,39.366,38.811,38.298,37.823,37.383,36.978,36.604,36.259,35.942,35.648,35.29 ,33.576,31.765,31.455,31.268,31.102,30.955,30.825,30.711,30.614,30.534,30.471,30.426,30.399,30.392,
                                                                     72.4  ,70.208,68.075,66.007,64.009,62.083,60.233,58.461,56.768,55.155,53.622,52.169,50.793,49.495,48.272,47.122,46.042,45.03 ,44.083,43.199,42.374,41.605,40.889,40.224,39.606,39.034,38.503,38.011,37.556,37.136,36.748,36.39 ,36.06 ,35.755,35.455,34.673,32.173,31.522,31.308,31.128,30.965,30.818,30.686,30.569,30.465,30.376,30.3  ,30.238,30.191,30.158,
                                                                     71.112,68.954,66.859,64.831,62.874,60.993,59.188,57.462,55.816,54.249,52.763,51.355,50.025,48.771,47.591,46.482,45.442,44.469,43.558,42.709,41.917,41.179,40.493,39.856,39.265,38.716,38.208,37.738,37.303,36.901,36.53 ,36.188,35.872,35.576,35.176,33.19 ,31.653,31.372,31.179,31.005,30.846,30.702,30.571,30.452,30.346,30.252,30.169,30.097,30.036,29.987,
                                                                     69.843,67.721,65.665,63.678,61.765,59.928,58.17 ,56.49 ,54.891,53.371,51.931,50.568,49.283,48.072,46.934,45.866,44.865,43.929,43.054,42.239,41.479,40.772,40.115,39.504,38.938,38.413,37.928,37.478,37.062,36.678,36.323,35.996,35.693,35.385,34.415,31.983,31.455,31.244,31.061,30.892,30.738,30.597,30.468,30.35 ,30.243,30.145,30.057,29.977,29.906,29.842,
                                                                     68.596,66.511,64.495,62.551,60.682,58.891,57.178,55.545,53.993,52.519,51.125,49.808,48.566,47.398,46.301,45.273,44.31 ,43.41 ,42.571,41.788,41.059,40.381,39.752,39.168,38.626,38.124,37.659,37.229,36.832,36.465,36.126,35.813,35.518,35.061,32.841,31.576,31.317,31.124,30.947,30.784,30.634,30.495,30.367,30.249,30.139,30.037,29.941,29.852,29.767,29.687,
                                                                     67.37 ,65.325,63.35 ,61.45 ,59.626,57.88 ,56.214,54.628,53.122,51.695,50.346,49.073,47.875,46.748,45.692,44.702,43.777,42.912,42.106,41.355,40.657,40.007,39.404,38.845,38.327,37.847,37.403,36.992,36.612,36.261,35.937,35.636,35.316,34.118,31.834,31.398,31.19 ,31.003,30.831,30.671,30.522,30.383,30.252,30.129,30.011,29.899,29.79 ,29.685,29.582,29.481,
                                                                     66.167,64.162,62.231,60.374,58.596,56.897,55.277,53.738,52.278,50.897,49.593,48.364,47.208,46.122,45.105,44.153,43.264,42.434,41.66 ,40.94 ,40.271,39.649,39.071,38.536,38.04 ,37.581,37.156,36.763,36.399,36.063,35.752,35.455,34.917,32.521,31.498,31.252,31.054,30.87 ,30.698,30.535,30.382,30.235,30.094,29.958,29.826,29.698,29.572,29.45 ,29.332,29.218,
                                                                     64.989,63.025,61.137,59.326,57.593,55.941,54.368,52.875,51.461,50.125,48.865,47.679,46.565,45.519,44.541,43.626,42.771,41.974,41.232,40.542,39.9  ,39.304,38.751,38.239,37.764,37.324,36.916,36.539,36.19 ,35.866,35.564,35.222,33.757,31.682,31.305,31.088,30.889,30.7  ,30.52 ,30.347,30.181,30.021,29.866,29.716,29.572,29.434,29.303,29.181,29.069,28.966,
                                                                     63.836,61.915,60.07 ,58.304,56.618,55.012,53.486,52.039,50.671,49.379,48.163,47.019,45.945,44.939,43.997,43.118,42.297,41.532,40.82 ,40.158,39.542,38.971,38.441,37.949,37.493,37.069,36.677,36.312,35.974,35.658,35.349,34.689,32.182,31.348,31.093,30.873,30.664,30.463,30.271,30.087,29.911,29.744,29.586,29.439,29.303,29.178,29.065,28.963,28.871,28.79 ,
                                                                     62.71 ,60.831,59.031,57.31 ,55.67 ,54.111,52.631,51.23 ,49.906,48.659,47.484,46.382,45.347,44.379,43.473,42.628,41.839,41.105,40.421,39.785,39.193,38.644,38.133,37.658,37.217,36.806,36.423,36.065,35.73 ,35.412,35.022,33.256,31.417,31.059,30.815,30.586,30.369,30.164,29.972,29.793,29.628,29.476,29.337,29.212,29.099,28.998,28.907,28.825,28.752,28.687,
                                                                     61.614,59.776,58.02 ,56.344,54.75 ,53.236,51.802,50.446,49.167,47.962,46.829,45.766,44.769,43.837,42.966,42.152,41.394,40.686,40.027,39.414,38.842,38.309,37.813,37.349,36.915,36.508,36.127,35.767,35.427,35.085,34.247,31.697,30.999,30.728,30.49 ,30.269,30.065,29.878,29.707,29.552,29.411,29.284,29.169,29.066,28.972,28.888,28.813,28.744,28.682,28.626,
                                                                     60.551,58.752,57.038,55.406,53.857,52.388,50.999,49.687,48.45 ,47.286,46.193,45.167,44.207,43.308,42.467,41.682,40.948,40.264,39.624,39.026,38.467,37.944,37.452,36.991,36.557,36.148,35.763,35.4  ,35.055,34.601,32.549,30.979,30.652,30.416,30.202,30.007,29.829,29.669,29.523,29.391,29.271,29.163,29.064,28.975,28.894,28.82 ,28.753,28.692,28.636,28.585,
                                                                     59.528,57.763,56.087,54.497,52.99 ,51.565,50.218,48.948,47.751,46.626,45.569,44.577,43.648,42.777,41.962,41.198,40.482,39.812,39.183,38.592,38.038,37.516,37.027,36.567,36.136,35.733,35.358,35.01 ,34.658,33.638,31.177,30.625,30.388,30.183,29.996,29.826,29.671,29.53 ,29.402,29.285,29.179,29.082,28.993,28.912,28.838,28.77 ,28.709,28.652,28.6  ,28.553,
                                                                     58.555,56.812,55.168,53.615,52.147,50.761,49.453,48.22 ,47.06 ,45.968,44.941,43.977,43.071,42.22 ,41.421,40.67 ,39.964,39.302,38.679,38.096,37.549,37.039,36.565,36.124,35.717,35.342,34.998,34.675,34.192,31.941,30.676,30.407,30.204,30.021,29.853,29.7  ,29.56 ,29.432,29.315,29.208,29.11 ,29.02 ,28.938,28.862,28.793,28.73 ,28.672,28.619,28.57 ,28.525,
                                                                     57.646,55.909,54.283,52.756,51.319,49.965,48.689,47.486,46.354,45.286,44.281,43.335,42.445,41.608,40.821,40.083,39.391,38.745,38.144,37.585,37.068,36.59 ,36.151,35.746,35.375,35.035,34.722,34.391,33.172,30.895,30.462,30.253,30.069,29.9  ,29.745,29.604,29.474,29.354,29.245,29.144,29.052,28.968,28.89 ,28.819,28.754,28.694,28.639,28.588,28.542,28.5  ,
                                                                     56.825,55.06 ,53.429,51.91 ,50.489,49.155,47.899,46.717,45.603,44.554,43.566,42.639,41.768,40.954,40.194,39.487,38.831,38.223,37.661,37.143,36.667,36.228,35.825,35.454,35.114,34.801,34.504,33.962,31.562,30.557,30.319,30.13 ,29.958,29.8  ,29.654,29.521,29.398,29.286,29.182,29.087,29    ,28.92 ,28.847,28.779,28.717,28.661,28.608,28.56 ,28.516,28.476,
                                                                     56.119,54.277,52.602,51.063,49.637,48.308,47.065,45.9  ,44.808,43.786,42.831,41.94 ,41.11 ,40.34 ,39.626,38.965,38.354,37.79 ,37.269,36.79 ,36.348,35.941,35.566,35.22 ,34.903,34.608,34.273,32.804,30.756,30.398,30.2  ,30.023,29.86 ,29.71 ,29.573,29.446,29.33 ,29.223,29.125,29.035,28.952,28.876,28.806,28.742,28.683,28.629,28.58 ,28.534,28.492,28.454,
                                                                     55.564,53.586,51.82 ,50.229,48.782,47.454,46.229,45.095,44.043,43.067,42.16 ,41.32 ,40.541,39.819,39.151,38.533,37.962,37.434,36.947,36.497,36.081,35.698,35.345,35.019,34.718,34.427,33.78 ,31.293,30.498,30.276,30.092,29.924,29.769,29.627,29.496,29.376,29.266,29.164,29.071,28.985,28.907,28.834,28.768,28.707,28.651,28.6  ,28.552,28.509,28.469,28.433,
                                                                     55.215,53.073,51.178,49.507,48.024,46.694,45.49 ,44.39 ,43.381,42.451,41.591,40.796,40.061,39.379,38.749,38.165,37.625,37.126,36.664,36.238,35.844,35.48 ,35.145,34.836,34.547,34.189,32.446,30.666,30.358,30.165,29.991,29.831,29.684,29.549,29.425,29.31 ,29.205,29.109,29.02 ,28.938,28.863,28.795,28.732,28.674,28.62 ,28.571,28.527,28.485,28.448,28.413,
                                                                     55.133,52.864,50.826,49.037,47.477,46.111,44.899,43.813,42.828,41.929,41.104,40.343,39.64 ,38.991,38.39 ,37.834,37.32 ,36.845,36.405,35.999,35.624,35.278,34.959,34.664,34.371,33.577,31.085,30.453,30.242,30.061,29.896,29.744,29.604,29.475,29.356,29.248,29.148,29.056,28.971,28.894,28.823,28.757,28.697,28.642,28.591,28.545,28.502,28.463,28.427,28.394,
                                                                     55.263,52.98 ,50.831,48.895,47.202,45.739,44.47 ,43.358,42.37 ,41.482,40.675,39.937,39.259,38.635,38.059,37.527,37.036,36.582,36.162,35.775,35.418,35.088,34.783,34.498,34.101,32.094,30.595,30.323,30.134,29.963,29.805,29.66 ,29.527,29.404,29.292,29.188,29.093,29.006,28.925,28.851,28.784,28.721,28.664,28.612,28.564,28.519,28.479,28.441,28.407,28.376,
                                                                     55.412,53.262,51.116,49.075,47.227,45.612,44.227,43.039,42.01 ,41.105,40.299,39.572,38.911,38.306,37.751,37.24 ,36.768,36.334,35.933,35.563,35.221,34.906,34.615,34.316,33.332,30.921,30.413,30.21 ,30.032,29.869,29.719,29.581,29.454,29.337,29.23 ,29.131,29.041,28.958,28.881,28.811,28.747,28.687,28.633,28.583,28.537,28.495,28.457,28.421,28.389,28.359,
                                                                     55.407,53.492,51.497,49.48 ,47.539,45.767,44.219,42.899,41.779,40.82 ,39.988,39.254,38.597,38.004,37.464,36.97 ,36.516,36.099,35.715,35.361,35.035,34.734,34.45 ,33.996,31.763,30.537,30.289,30.104,29.935,29.78 ,29.637,29.505,29.385,29.274,29.171,29.078,28.991,28.912,28.84 ,28.773,28.711,28.655,28.603,28.556,28.512,28.472,28.435,28.402,28.371,28.342,
                                                                     55.193,53.528,51.772,49.925,48.033,46.188,44.487,42.996,41.731,40.668,39.77 ,39    ,38.329,37.735,37.202,36.72 ,36.28 ,35.878,35.509,35.17 ,34.858,34.569,34.259,33.044,30.792,30.376,30.179,30.004,29.843,29.695,29.559,29.434,29.319,29.213,29.116,29.026,28.944,28.869,28.8  ,28.736,28.678,28.624,28.575,28.53 ,28.488,28.45 ,28.415,28.383,28.354,28.327,
                                                                     54.804,53.346,51.832,50.228,48.524,46.754,45.002,43.367,41.927,40.707,39.691,38.844,38.129,37.513,36.974,36.494,36.063,35.672,35.316,34.989,34.69 ,34.405,33.868,31.47 ,30.487,30.257,30.075,29.908,29.755,29.614,29.484,29.365,29.256,29.155,29.063,28.978,28.9  ,28.828,28.762,28.702,28.646,28.595,28.548,28.505,28.466,28.43 ,28.396,28.366,28.338,28.312,
                                                                     54.301,52.996,51.677,50.305,48.844,47.277,45.627,43.966,42.394,40.998,39.815,38.838,38.033,37.363,36.795,36.304,35.871,35.485,35.137,34.821,34.53 ,34.199,32.72 ,30.692,30.341,30.149,29.976,29.817,29.671,29.537,29.413,29.3  ,29.196,29.1  ,29.012,28.931,28.857,28.789,28.726,28.668,28.616,28.567,28.523,28.482,28.444,28.41 ,28.378,28.349,28.322,28.298,
                                                                     53.736,52.541,51.362,50.172,48.933,47.607,46.171,44.637,43.065,41.555,40.198,39.045,38.097,37.326,36.693,36.165,35.715,35.323,34.977,34.666,34.37 ,33.716,31.228,30.445,30.227,30.047,29.882,29.731,29.591,29.464,29.346,29.238,29.139,29.048,28.964,28.887,28.816,28.751,28.692,28.637,28.587,28.541,28.498,28.459,28.424,28.391,28.361,28.333,28.308,28.285,
                                                                     53.144,52.029,50.951,49.89 ,48.82 ,47.704,46.505,45.194,43.772,42.288,40.836,39.514,38.384,37.458,36.712,36.109,35.615,35.199,34.842,34.527,34.152,32.387,30.619,30.313,30.122,29.95 ,29.793,29.648,29.516,29.394,29.282,29.179,29.085,28.998,28.918,28.845,28.778,28.716,28.659,28.607,28.559,28.515,28.475,28.438,28.404,28.373,28.344,28.318,28.294,28.272,
                                                                     52.548,51.494,50.49 ,49.521,48.569,47.61 ,46.608,45.527,44.333,43.018,41.62 ,40.226,38.936,37.823,36.911,36.183,35.603,35.135,34.748,34.4  ,33.564,31.056,30.422,30.207,30.025,29.86 ,29.709,29.571,29.444,29.328,29.221,29.123,29.033,28.951,28.875,28.805,28.741,28.682,28.628,28.579,28.533,28.491,28.453,28.418,28.386,28.356,28.329,28.304,28.281,28.26 ,
                                                                     51.959,50.956,50.009,49.107,48.239,47.388,46.53 ,45.633,44.66 ,43.576,42.364,41.051,39.713,38.453,37.353,36.449,35.732,35.167,34.715,34.211,32.127,30.606,30.313,30.113,29.936,29.775,29.63 ,29.497,29.376,29.265,29.163,29.07 ,28.984,28.906,28.833,28.767,28.706,28.65 ,28.599,28.552,28.508,28.469,28.432,28.399,28.368,28.34 ,28.314,28.29 ,28.268,28.248,
                                                                     51.384,50.426,49.525,48.675,47.869,47.093,46.334,45.569,44.768,43.895,42.912,41.798,40.567,39.287,38.056,36.965,36.063,35.349,34.762,33.579,31.058,30.483,30.237,30.034,29.857,29.699,29.557,29.429,29.312,29.206,29.109,29.019,28.938,28.863,28.794,28.731,28.673,28.62 ,28.571,28.526,28.485,28.447,28.412,28.38 ,28.351,28.324,28.3  ,28.277,28.257,28.238,
                                                                     50.827,49.909,49.049,48.243,47.483,46.763,46.073,45.397,44.716,44.002,43.22 ,42.331,41.31 ,40.16 ,38.936,37.734,36.651,35.74 ,34.821,32.26 ,30.841,30.462,30.196,29.977,29.793,29.633,29.492,29.366,29.253,29.151,29.057,28.972,28.894,28.822,28.757,28.696,28.641,28.591,28.544,28.501,28.462,28.426,28.393,28.363,28.335,28.31 ,28.286,28.265,28.245,28.227,
                                                                     50.289,49.409,48.586,47.817,47.097,46.42 ,45.78 ,45.166,44.566,43.961,43.325,42.626,41.824,40.89 ,39.818,38.652,37.478,36.355,34.459,31.719,30.959,30.532,30.209,29.956,29.752,29.583,29.438,29.312,29.2  ,29.1  ,29.009,28.927,28.852,28.784,28.721,28.664,28.611,28.563,28.519,28.478,28.441,28.407,28.375,28.347,28.32 ,28.296,28.274,28.253,28.234,28.217};

            return arrayCurvatureMax[y * arrayWidth + x];
        }
        case 6:
        {
            static const std::vector<double> arrayCurvatureMax =    {101.5 ,99.257,96.879,94.387,91.803,89.152,86.463,83.767,81.094,78.475,75.936,73.499,71.184,69.003,66.965,65.073,63.327,61.723,60.255,58.915,57.693,56.579,55.564,54.637,53.789,53.011,52.296,51.635,51.024,50.455,49.924,49.427,48.96 ,48.521,48.106,47.715,47.345,46.995,46.665,46.353,46.06 ,45.784,45.525,45.282,45.055,44.842,44.643,44.456,44.28 ,44.112,
                                                                     100.44,98.148,95.733,93.209,90.599,87.93 ,85.231,82.534,79.868,77.262,74.744,72.333,70.048,67.9  ,65.895,64.038,62.326,60.755,59.318,58.006,56.81 ,55.72 ,54.726,53.818,52.988,52.225,51.524,50.876,50.276,49.717,49.196,48.709,48.252,47.823,47.419,47.038,46.679,46.341,46.023,45.724,45.444,45.181,44.936,44.706,44.492,44.291,44.102,43.924,43.753,43.588,
                                                                     99.347,97.019,94.568,92.014,89.381,86.696,83.99 ,81.293,78.635,76.046,73.549,71.166,68.911,66.796,64.827,63.004,61.326,59.787,58.38 ,57.097,55.927,54.86 ,53.888,52.999,52.185,51.439,50.751,50.116,49.528,48.981,48.471,47.995,47.549,47.13 ,46.737,46.368,46.022,45.697,45.392,45.107,44.84 ,44.591,44.359,44.142,43.939,43.748,43.568,43.395,43.227,43.062,
                                                                     98.238,95.871,93.385,90.803,88.148,85.45 ,82.738,80.044,77.397,74.825,72.352,69.997,67.775,65.694,63.759,61.971,60.326,58.82 ,57.443,56.188,55.043,54    ,53.049,52.179,51.383,50.652,49.979,49.357,48.782,48.247,47.749,47.285,46.85 ,46.444,46.064,45.708,45.375,45.064,44.773,44.502,44.249,44.014,43.794,43.589,43.395,43.212,43.037,42.866,42.698,42.53 ,
                                                                     97.109,94.704,92.185,89.576,86.902,84.193,81.478,78.788,76.154,73.602,71.154,68.829,66.639,64.592,62.692,60.939,59.328,57.854,56.507,55.28 ,54.161,53.141,52.211,51.36 ,50.582,49.867,49.208,48.601,48.039,47.517,47.032,46.581,46.16 ,45.767,45.401,45.059,44.74 ,44.443,44.167,43.91 ,43.671,43.448,43.24 ,43.044,42.858,42.68 ,42.506,42.335,42.162,41.987,
                                                                     95.961,93.518,90.969,88.335,85.644,82.925,80.209,77.527,74.908,72.377,69.956,67.661,65.504,63.493,61.628,59.91 ,58.333,56.89 ,55.573,54.373,53.28 ,52.283,51.374,50.543,49.782,49.084,48.442,47.849,47.301,46.794,46.323,45.885,45.479,45.101,44.749,44.422,44.118,43.836,43.575,43.332,43.105,42.894,42.694,42.505,42.324,42.147,41.972,41.795,41.615,41.428,
                                                                     94.793,92.316,89.737,87.08 ,84.374,81.648,78.934,76.261,73.659,71.151,68.757,66.494,64.372,62.396,60.567,58.883,57.34 ,55.929,54.642,53.47 ,52.402,51.428,50.541,49.73 ,48.987,48.306,47.68 ,47.103,46.571,46.079,45.623,45.201,44.809,44.447,44.111,43.799,43.511,43.244,42.996,42.765,42.55 ,42.347,42.155,41.97 ,41.789,41.61 ,41.43 ,41.245,41.053,40.852,
                                                                     93.608,91.096,88.49 ,85.812,83.093,80.363,77.653,74.992,72.408,69.924,67.56 ,65.33 ,63.242,61.302,59.509,57.861,56.351,54.972,53.715,52.57 ,51.527,50.578,49.712,48.921,48.198,47.535,46.926,46.366,45.85 ,45.374,44.934,44.529,44.154,43.807,43.487,43.191,42.918,42.665,42.43 ,42.21 ,42.004,41.808,41.619,41.434,41.25 ,41.065,40.876,40.679,40.473,40.255,
                                                                     92.405,89.861,87.229,84.533,81.803,79.071,76.366,73.719,71.156,68.698,66.365,64.169,62.117,60.213,58.456,56.843,55.367,54.019,52.792,51.675,50.658,49.733,48.889,48.119,47.416,46.771,46.181,45.639,45.14 ,44.682,44.26 ,43.871,43.513,43.183,42.879,42.599,42.34 ,42.099,41.875,41.664,41.464,41.27 ,41.081,40.893,40.703,40.509,40.307,40.095,39.872,39.635,
                                                                     91.186,88.611,85.955,83.243,80.505,77.773,75.077,72.445,69.904,67.475,65.173,63.012,60.996,59.129,57.409,55.831,54.389,53.073,51.876,50.787,49.796,48.895,48.074,47.326,46.643,46.019,45.448,44.924,44.444,44.004,43.6  ,43.229,42.888,42.575,42.287,42.021,41.775,41.545,41.329,41.124,40.926,40.732,40.539,40.345,40.145,39.937,39.72 ,39.491,39.248,38.992,
                                                                     89.951,87.348,84.669,81.943,79.2  ,76.47 ,73.785,71.171,68.655,66.255,63.986,61.86 ,59.882,58.052,56.369,54.826,53.418,52.135,50.968,49.907,48.943,48.067,47.27 ,46.544,45.883,45.28 ,44.729,44.225,43.764,43.343,42.958,42.605,42.281,41.984,41.711,41.458,41.222,41.001,40.79 ,40.587,40.388,40.19 ,39.989,39.784,39.571,39.348,39.113,38.865,38.602,38.325,
                                                                     88.702,86.072,83.374,80.636,77.889,75.165,72.492,69.899,67.408,65.039,62.805,60.716,58.775,56.983,55.337,53.831,52.457,51.206,50.07 ,49.038,48.101,47.251,46.478,45.776,45.137,44.555,44.026,43.543,43.102,42.701,42.334,41.999,41.692,41.409,41.149,40.906,40.679,40.462,40.253,40.049,39.845,39.639,39.428,39.209,38.98 ,38.739,38.485,38.217,37.934,37.636,
                                                                     87.439,84.785,82.069,79.322,76.575,73.858,71.201,68.63 ,66.167,63.83 ,61.632,59.58 ,57.678,55.924,54.316,52.846,51.506,50.289,49.184,48.182,47.273,46.449,45.702,45.024,44.408,43.849,43.341,42.88 ,42.46 ,42.078,41.729,41.411,41.119,40.85 ,40.6  ,40.365,40.142,39.927,39.717,39.507,39.295,39.078,38.853,38.617,38.371,38.111,37.837,37.548,37.245,36.927,
                                                                     86.166,83.488,80.758,78.004,75.259,72.552,69.912,67.365,64.932,62.629,60.468,58.454,56.591,54.877,53.306,51.873,50.57 ,49.386,48.313,47.341,46.461,45.664,44.943,44.29 ,43.699,43.163,42.678,42.238,41.838,41.475,41.144,40.841,40.562,40.304,40.061,39.831,39.609,39.392,39.176,38.958,38.734,38.503,38.261,38.008,37.742,37.463,37.169,36.86 ,36.537,36.2  ,
                                                                     84.882,82.184,79.441,76.683,73.942,71.248,68.628,66.107,63.706,61.438,59.315,57.341,55.518,53.844,52.312,50.916,49.649,48.499,47.458,46.517,45.667,44.899,44.205,43.578,43.011,42.499,42.037,41.618,41.238,40.893,40.578,40.288,40.02 ,39.769,39.531,39.302,39.077,38.854,38.629,38.399,38.161,37.913,37.654,37.382,37.096,36.796,36.483,36.154,35.813,35.458,
                                                                     83.59 ,80.874,78.121,75.362,72.628,69.948,67.351,64.859,62.491,60.26 ,58.176,56.243,54.461,52.827,51.335,49.977,48.746,47.631,46.624,45.715,44.895,44.156,43.489,42.889,42.348,41.86 ,41.419,41.021,40.66 ,40.331,40.029,39.75 ,39.49 ,39.243,39.005,38.773,38.543,38.31 ,38.073,37.828,37.574,37.308,37.03 ,36.739,36.433,36.114,35.781,35.435,35.076,34.705,
                                                                     82.293,79.561,76.801,74.043,71.318,68.656,66.084,63.622,61.289,59.097,57.053,55.162,53.422,51.829,50.377,49.059,47.865,46.786,45.813,44.937,44.148,43.438,42.8  ,42.226,41.71 ,41.246,40.827,40.448,40.103,39.788,39.497,39.226,38.969,38.723,38.482,38.243,38.003,37.759,37.507,37.246,36.974,36.689,36.392,36.081,35.757,35.419,35.068,34.705,34.33 ,33.944,
                                                                     80.991,78.247,75.482,72.728,70.016,67.374,64.828,62.399,60.103,57.951,55.95 ,54.101,52.404,50.854,49.443,48.164,47.009,45.967,45.029,44.185,43.428,42.749,42.139,41.592,41.1  ,40.658,40.259,39.897,39.567,39.263,38.979,38.712,38.456,38.206,37.959,37.71 ,37.458,37.199,36.931,36.652,36.361,36.058,35.741,35.412,35.07 ,34.715,34.348,33.969,33.58 ,33.181,
                                                                     79.689,76.934,74.168,71.421,68.724,66.105,63.589,61.194,58.937,56.827,54.869,53.064,51.411,49.904,48.535,47.298,46.181,45.176,44.274,43.464,42.739,42.09 ,41.508,40.987,40.519,40.098,39.716,39.369,39.05 ,38.753,38.474,38.207,37.947,37.69 ,37.433,37.173,36.906,36.63 ,36.344,36.047,35.738,35.416,35.082,34.735,34.376,34.005,33.623,33.231,32.83 ,32.42 ,
                                                                     78.389,75.627,72.861,70.124,67.446,64.852,62.368,60.01 ,57.794,55.727,53.814,52.055,50.446,48.983,47.658,46.462,45.385,44.418,43.552,42.777,42.083,41.464,40.909,40.413,39.966,39.563,39.197,38.862,38.55 ,38.257,37.978,37.707,37.441,37.174,36.905,36.63 ,36.347,36.054,35.751,35.435,35.108,34.769,34.417,34.054,33.68 ,33.295,32.9  ,32.497,32.084,31.665,
                                                                     77.094,74.327,71.566,68.843,66.185,63.62 ,61.169,58.851,56.677,54.655,52.788,51.076,49.514,48.096,46.815,45.661,44.625,43.696,42.866,42.125,41.463,40.872,40.343,39.869,39.442,39.055,38.701,38.373,38.066,37.773,37.491,37.213,36.937,36.658,36.375,36.084,35.784,35.473,35.152,34.819,34.475,34.119,33.752,33.374,32.987,32.59 ,32.184,31.77 ,31.349,30.921,
                                                                     75.809,73.04 ,70.287,67.579,64.945,62.411,59.998,57.721,55.591,53.616,51.797,50.133,48.618,47.247,46.01 ,44.899,43.904,43.014,42.219,41.511,40.88 ,40.316,39.811,39.357,38.947,38.572,38.226,37.902,37.596,37.3  ,37.01 ,36.723,36.435,36.142,35.843,35.535,35.218,34.891,34.553,34.203,33.843,33.472,33.092,32.701,32.302,31.894,31.479,31.057,30.629,30.195,
                                                                     74.536,71.77 ,69.027,66.338,63.732,61.231,58.857,56.624,54.541,52.615,50.845,49.23 ,47.764,46.439,45.248,44.18 ,43.225,42.373,41.614,40.938,40.335,39.797,39.313,38.876,38.478,38.112,37.771,37.448,37.138,36.836,36.537,36.238,35.936,35.628,35.313,34.989,34.655,34.311,33.957,33.593,33.218,32.834,32.441,32.04 ,31.631,31.215,30.792,30.363,29.929,29.49 ,
                                                                     73.281,70.52 ,67.792,65.125,62.549,60.085,57.753,55.566,53.532,51.655,49.936,48.371,46.954,45.678,44.532,43.507,42.593,41.778,41.053,40.407,39.83 ,39.314,38.848,38.425,38.036,37.675,37.335,37.009,36.693,36.381,36.071,35.759,35.442,35.119,34.787,34.447,34.098,33.739,33.37 ,32.992,32.606,32.21 ,31.807,31.397,30.98 ,30.556,30.128,29.694,29.255,28.813,
                                                                     72.049,69.297,66.586,63.945,61.402,58.978,56.69 ,54.551,52.568,50.743,49.076,47.562,46.195,44.966,43.866,42.884,42.009,41.23 ,40.537,39.918,39.365,38.868,38.416,38.003,37.62 ,37.26 ,36.917,36.586,36.261,35.938,35.615,35.289,34.957,34.618,34.271,33.916,33.552,33.18 ,32.798,32.409,32.012,31.607,31.195,30.778,30.354,29.925,29.492,29.055,28.613,28.169,
                                                                     70.845,68.105,65.415,62.804,60.297,57.915,55.674,53.585,51.654,49.883,48.269,46.807,45.491,44.31 ,43.254,42.313,41.476,40.731,40.067,39.473,38.94 ,38.457,38.017,37.609,37.228,36.866,36.518,36.178,35.843,35.508,35.171,34.831,34.484,34.131,33.77 ,33.401,33.024,32.64 ,32.248,31.848,31.442,31.03 ,30.611,30.188,29.76 ,29.327,28.891,28.451,28.009,27.564,
                                                                     69.674,66.95 ,64.285,61.706,59.239,56.903,54.711,52.675,50.798,49.081,47.52 ,46.111,44.844,43.71 ,42.699,41.798,40.996,40.281,39.643,39.071,38.554,38.083,37.649,37.244,36.862,36.495,36.139,35.789,35.442,35.094,34.744,34.389,34.029,33.662,33.288,32.907,32.52 ,32.125,31.724,31.316,30.903,30.485,30.062,29.634,29.203,28.768,28.33 ,27.889,27.447,27.002,
                                                                     68.542,65.838,63.201,60.659,58.235,55.946,53.807,51.824,50.003,48.341,46.835,45.478,44.26 ,43.172,42.202,41.338,40.569,39.882,39.266,38.711,38.207,37.743,37.313,36.907,36.52 ,36.147,35.782,35.421,35.062,34.701,34.338,33.97 ,33.597,33.218,32.833,32.442,32.045,31.642,31.233,30.82 ,30.401,29.979,29.552,29.122,28.689,28.253,27.815,27.374,26.932,26.489,
                                                                     67.455,64.775,62.171,59.668,57.29 ,55.052,52.966,51.04 ,49.274,47.668,46.216,44.91 ,43.74 ,42.696,41.765,40.936,40.196,39.533,38.936,38.394,37.898,37.439,37.008,36.599,36.206,35.824,35.449,35.077,34.706,34.333,33.957,33.577,33.193,32.804,32.409,32.01 ,31.605,31.195,30.781,30.363,29.941,29.516,29.087,28.656,28.222,27.787,27.349,26.91 ,26.47 ,26.029,
                                                                     66.419,63.767,61.198,58.739,56.409,54.224,52.194,50.325,48.617,47.066,45.667,44.411,43.287,42.284,41.39 ,40.591,39.876,39.233,38.65 ,38.118,37.627,37.169,36.736,36.322,35.922,35.53 ,35.144,34.761,34.378,33.994,33.607,33.217,32.824,32.425,32.023,31.617,31.206,30.792,30.374,29.953,29.529,29.102,28.673,28.241,27.808,27.374,26.938,26.501,26.064,25.626,
                                                                     65.437,62.817,60.289,57.876,55.598,53.468,51.496,49.684,48.032,46.537,45.19 ,43.981,42.901,41.936,41.074,40.302,39.608,38.981,38.41 ,37.884,37.395,36.935,36.498,36.077,35.668,35.267,34.871,34.477,34.084,33.69 ,33.293,32.895,32.493,32.088,31.68 ,31.268,30.853,30.436,30.015,29.592,29.167,28.74 ,28.311,27.881,27.45 ,27.017,26.584,26.151,25.716,25.282,
                                                                     64.511,61.929,59.445,57.082,54.858,52.786,50.872,49.118,47.523,46.081,44.784,43.621,42.58 ,41.65 ,40.816,40.067,39.391,38.776,38.212,37.69 ,37.2  ,36.737,36.294,35.866,35.448,35.038,34.632,34.229,33.827,33.424,33.02 ,32.614,32.205,31.795,31.382,30.967,30.549,30.13 ,29.708,29.285,28.86 ,28.434,28.007,27.578,27.149,26.72 ,26.29 ,25.86 ,25.43 ,25    ,
                                                                     63.641,61.1  ,58.665,56.356,54.19 ,52.176,50.322,48.627,47.087,45.697,44.448,43.327,42.323,41.423,40.614,39.884,39.222,38.616,38.056,37.535,37.043,36.576,36.126,35.69 ,35.265,34.846,34.431,34.02 ,33.609,33.199,32.789,32.377,31.964,31.55 ,31.134,30.716,30.297,29.877,29.455,29.033,28.609,28.184,27.759,27.333,26.907,26.481,26.055,25.629,25.203,24.778,
                                                                     62.817,60.325,57.944,55.693,53.587,51.635,49.841,48.204,46.72 ,45.38 ,44.175,43.093,42.122,41.25 ,40.462,39.748,39.096,38.497,37.939,37.417,36.923,36.45 ,35.994,35.551,35.118,34.692,34.27 ,33.851,33.434,33.019,32.603,32.187,31.771,31.354,30.936,30.517,30.098,29.678,29.257,28.835,28.413,27.991,27.568,27.146,26.723,26.301,25.879,25.457,25.036,24.616,
                                                                     62.027,59.59 ,57.269,55.081,53.039,51.151,49.418,47.84 ,46.409,45.118,43.956,42.911,41.971,41.122,40.353,39.652,39.009,38.414,37.858,37.334,36.836,36.359,35.897,35.448,35.008,34.575,34.147,33.723,33.301,32.881,32.462,32.043,31.625,31.206,30.788,30.369,29.95 ,29.53 ,29.111,28.691,28.271,27.852,27.432,27.013,26.594,26.176,25.758,25.341,24.925,24.509,
                                                                     61.257,58.88 ,56.622,54.501,52.527,50.704,49.035,47.516,46.139,44.896,43.776,42.765,41.853,41.027,40.275,39.586,38.951,38.36 ,37.805,37.28 ,36.779,36.298,35.832,35.378,34.933,34.495,34.063,33.634,33.209,32.786,32.364,31.944,31.524,31.105,30.686,30.267,29.849,29.431,29.014,28.596,28.179,27.762,27.346,26.93 ,26.515,26.101,25.687,25.275,24.863,24.452,
                                                                     60.507,58.188,55.993,53.938,52.031,50.275,48.669,47.209,45.885,44.69 ,43.61 ,42.635,41.751,40.947,40.212,39.535,38.908,38.322,37.77 ,37.246,36.744,36.26 ,35.791,35.334,34.886,34.445,34.01 ,33.579,33.151,32.727,32.304,31.883,31.463,31.044,30.626,30.208,29.792,29.375,28.96 ,28.545,28.13 ,27.717,27.304,26.891,26.48 ,26.069,25.659,25.251,24.843,24.437,
                                                                     59.806,57.535,55.394,53.397,51.549,49.853,48.304,46.898,45.625,44.475,43.434,42.493,41.637,40.857,40.14 ,39.478,38.862,38.283,37.736,37.215,36.716,36.233,35.764,35.307,34.858,34.416,33.98 ,33.549,33.121,32.696,32.274,31.853,31.434,31.016,30.599,30.183,29.768,29.354,28.941,28.528,28.117,27.706,27.296,26.887,26.479,26.072,25.666,25.261,24.857,24.454,
                                                                     59.199,56.963,54.863,52.909,51.106,49.454,47.95 ,46.587,45.354,44.24 ,43.233,42.321,41.491,40.733,40.036,39.39 ,38.787,38.22 ,37.682,37.169,36.675,36.197,35.733,35.278,34.833,34.394,33.96 ,33.531,33.105,32.682,32.261,31.843,31.425,31.01 ,30.595,30.182,29.769,29.357,28.947,28.537,28.128,27.72 ,27.313,26.907,26.502,26.098,25.695,25.293,24.893,24.494,
                                                                     58.714,56.506,54.436,52.513,50.741,49.119,47.643,46.305,45.095,44.003,43.016,42.122,41.31 ,40.567,39.885,39.253,38.663,38.108,37.582,37.08 ,36.596,36.128,35.672,35.226,34.788,34.356,33.928,33.505,33.084,32.666,32.25 ,31.836,31.422,31.01 ,30.599,30.189,29.78 ,29.371,28.964,28.557,28.151,27.746,27.342,26.939,26.537,26.136,25.736,25.337,24.94 ,24.543,
                                                                     58.344,56.163,54.12 ,52.223,50.474,48.872,47.412,46.088,44.888,43.803,42.821,41.931,41.121,40.382,39.703,39.076,38.492,37.944,37.427,36.934,36.46 ,36.003,35.558,35.124,34.697,34.275,33.859,33.445,33.034,32.625,32.217,31.81 ,31.403,30.998,30.592,30.188,29.783,29.379,28.976,28.573,28.171,27.769,27.369,26.969,26.57 ,26.172,25.775,25.379,24.984,24.59 ,
                                                                     58.061,55.908,53.893,52.02 ,50.293,48.708,47.262,45.945,44.749,43.664,42.678,41.782,40.965,40.217,39.529,38.894,38.304,37.752,37.233,36.74 ,36.27 ,35.818,35.381,34.956,34.539,34.13 ,33.725,33.324,32.926,32.529,32.133,31.738,31.342,30.946,30.55 ,30.154,29.757,29.361,28.964,28.567,28.17 ,27.773,27.377,26.981,26.585,26.191,25.797,25.404,25.012,24.621,
                                                                     57.834,55.709,53.72 ,51.872,50.166,48.599,47.166,45.857,44.665,43.579,42.589,41.685,40.857,40.096,39.394,38.744,38.139,37.573,37.041,36.537,36.058,35.6  ,35.16 ,34.734,34.321,33.917,33.52 ,33.129,32.741,32.357,31.974,31.592,31.21 ,30.827,30.444,30.06 ,29.675,29.289,28.901,28.513,28.124,27.735,27.345,26.955,26.566,26.176,25.786,25.397,25.009,24.621,
                                                                     57.633,55.536,53.572,51.748,50.063,48.513,47.093,45.794,44.607,43.523,42.531,41.621,40.785,40.013,39.298,38.634,38.012,37.429,36.879,36.359,35.863,35.39 ,34.936,34.498,34.075,33.665,33.265,32.873,32.488,32.109,31.733,31.361,30.99 ,30.62 ,30.25 ,29.88 ,29.509,29.136,28.762,28.386,28.009,27.631,27.251,26.87 ,26.489,26.107,25.724,25.341,24.958,24.575,
                                                                     57.435,55.362,53.423,51.62 ,49.955,48.422,47.015,45.726,44.546,43.466,42.474,41.563,40.722,39.943,39.22 ,38.544,37.91 ,37.313,36.748,36.211,35.698,35.208,34.736,34.281,33.841,33.415,33.001,32.597,32.203,31.816,31.437,31.063,30.694,30.329,29.966,29.605,29.244,28.884,28.523,28.162,27.799,27.434,27.068,26.7  ,26.33 ,25.959,25.587,25.214,24.839,24.464,
                                                                     57.219,55.168,53.25 ,51.468,49.821,48.304,46.91 ,45.632,44.46 ,43.384,42.396,41.485,40.642,39.86 ,39.131,38.449,37.807,37.2  ,36.624,36.075,35.549,35.044,34.557,34.086,33.63 ,33.186,32.755,32.334,31.923,31.521,31.127,30.741,30.362,29.989,29.621,29.258,28.899,28.542,28.188,27.835,27.482,27.13 ,26.777,26.423,26.067,25.71 ,25.352,24.992,24.63 ,24.267,
                                                                     56.968,54.937,53.039,51.276,49.645,48.143,46.763,45.495,44.331,43.262,42.277,41.368,40.527,39.744,39.013,38.326,37.68 ,37.067,36.484,35.927,35.392,34.877,34.379,33.896,33.427,32.969,32.523,32.086,31.659,31.24 ,30.83 ,30.427,30.031,29.642,29.26 ,28.884,28.514,28.15 ,27.79 ,27.434,27.082,26.733,26.385,26.039,25.694,25.349,25.004,24.657,24.31 ,23.961,
                                                                     56.669,54.656,52.776,51.031,49.416,47.929,46.561,45.304,44.149,43.087,42.107,41.202,40.362,39.58 ,38.848,38.161,37.512,36.896,36.308,35.746,35.206,34.684,34.179,33.688,33.21 ,32.743,32.286,31.838,31.398,30.966,30.541,30.124,29.712,29.307,28.909,28.517,28.13 ,27.75 ,27.376,27.007,26.644,26.286,25.932,25.583,25.237,24.895,24.555,24.216,23.879,23.541,
                                                                     56.312,54.316,52.453,50.724,49.126,47.653,46.297,45.052,43.906,42.851,41.878,40.977,40.141,39.361,38.63 ,37.942,37.292,36.674,36.085,35.52 ,34.976,34.45 ,33.94 ,33.444,32.96 ,32.487,32.023,31.567,31.119,30.678,30.244,29.815,29.393,28.975,28.564,28.157,27.756,27.36 ,26.97 ,26.585,26.206,25.833,25.465,25.103,24.746,24.395,24.048,23.705,23.367,23.032,
                                                                     55.889,53.909,52.064,50.351,48.768,47.31 ,45.968,44.734,43.598,42.552,41.586,40.69 ,39.858,39.081,38.352,37.666,37.016,36.398,35.807,35.24 ,34.694,34.166,33.653,33.154,32.666,32.189,31.72 ,31.26 ,30.807,30.36 ,29.919,29.484,29.054,28.628,28.208,27.791,27.38 ,26.973,26.57 ,26.172,25.779,25.391,25.009,24.631,24.259,23.892,23.532,23.176,22.826,22.482};

            return arrayCurvatureMax[y * arrayWidth + x];
        }
        case 7:
        {
            static const std::vector<double> arrayCurvatureMax =    {127.51,122.28,117.29,112.55,108.08,103.89,99.963,96.309,92.914,89.766,86.851,84.151,81.65 ,79.331,77.178,75.175,73.307,71.559,69.919,68.374,66.911,65.52 ,64.188,62.905,61.658,60.439,59.237,58.049,56.874,55.717,54.588,53.497,52.454,51.461,50.518,49.619,48.756,47.92 ,47.102,46.293,45.485,44.671,43.846,43.004,42.144,41.265,40.369,39.461,38.547,37.638,
                                                                     125.1 ,119.88,114.89,110.16,105.7 ,101.51,97.605,93.969,90.594,87.467,84.574,81.897,79.421,77.128,75.002,73.027,71.188,69.471,67.862,66.348,64.918,63.56 ,62.261,61.01 ,59.796,58.608,57.439,56.284,55.145,54.027,52.94 ,51.895,50.899,49.955,49.06 ,48.207,47.388,46.594,45.815,45.042,44.269,43.487,42.692,41.88 ,41.048,40.197,39.329,38.45 ,37.567,36.689,
                                                                     122.65,117.43,112.45,107.73,103.29,99.124,95.239,91.627,88.277,85.178,82.313,79.666,77.221,74.96 ,72.866,70.925,69.12 ,67.438,65.865,64.387,62.994,61.671,60.408,59.193,58.014,56.861,55.726,54.607,53.506,52.428,51.386,50.388,49.441,48.545,47.697,46.889,46.113,45.359,44.617,43.88 ,43.139,42.388,41.622,40.837,40.031,39.206,38.365,37.513,36.66 ,35.813,
                                                                     120.17,114.96,109.99,105.29,100.87,96.734,92.878,89.296,85.978,82.912,80.082,77.471,75.062,72.838,70.782,68.879,67.114,65.471,63.938,62.5  ,61.146,59.863,58.639,57.461,56.319,55.202,54.103,53.021,51.959,50.924,49.927,48.976,48.077,47.229,46.427,45.663,44.928,44.212,43.506,42.802,42.091,41.369,40.63 ,39.87 ,39.089,38.289,37.473,36.647,35.821,35.004,
                                                                     117.67,112.47,107.53,102.86,98.463,94.357,90.535,86.99 ,83.711,80.684,77.894,75.323,72.956,70.774,68.761,66.901,65.179,63.579,62.088,60.694,59.382,58.14 ,56.956,55.817,54.713,53.633,52.572,51.528,50.505,49.513,48.562,47.659,46.808,46.006,45.249,44.527,43.831,43.151,42.479,41.805,41.123,40.427,39.712,38.976,38.218,37.44 ,36.647,35.846,35.046,34.258,
                                                                     115.16,109.99,105.07,100.43,96.075,92.007,88.226,84.723,81.487,78.505,75.76 ,73.236,70.914,68.779,66.812,64.999,63.323,61.769,60.324,58.974,57.706,56.507,55.364,54.265,53.199,52.157,51.133,50.127,49.145,48.197,47.291,46.434,45.63 ,44.874,44.159,43.477,42.818,42.172,41.53 ,40.885,40.229,39.557,38.864,38.149,37.412,36.655,35.884,35.106,34.331,33.569,
                                                                     112.66,107.52,102.64,98.034,93.72 ,89.697,85.962,82.507,79.32 ,76.387,73.692,71.218,68.946,66.861,64.944,63.179,61.552,60.046,58.648,57.345,56.122,54.966,53.865,52.806,51.778,50.773,49.786,48.818,47.877,46.972,46.111,45.3  ,44.54 ,43.827,43.153,42.509,41.884,41.269,40.656,40.037,39.405,38.754,38.082,37.386,36.667,35.929,35.178,34.422,33.671,32.935,
                                                                     110.18,105.08,100.24,95.678,91.412,87.439,83.755,80.353,77.219,74.34 ,71.698,69.277,67.059,65.026,63.161,61.448,59.871,58.415,57.066,55.809,54.631,53.519,52.459,51.439,50.449,49.48 ,48.53 ,47.6  ,46.699,45.836,45.019,44.252,43.536,42.863,42.227,41.618,41.025,40.439,39.852,39.256,38.645,38.014,37.359,36.68 ,35.979,35.259,34.526,33.791,33.062,32.35 ,
                                                                     107.74,102.68,97.884,93.376,89.163,85.244,81.617,78.271,75.194,72.372,69.787,67.422,65.259,63.281,61.47 ,59.809,58.284,56.878,55.577,54.367,53.234,52.164,51.145,50.164,49.21 ,48.276,47.362,46.469,45.607,44.786,44.011,43.287,42.611,41.978,41.377,40.8  ,40.236,39.677,39.113,38.538,37.946,37.331,36.693,36.029,35.343,34.638,33.924,33.207,32.5  ,31.812,
                                                                     105.35,100.33,95.591,91.139,86.982,83.123,79.555,76.269,73.252,70.489,67.963,65.657,63.551,61.628,59.872,58.265,56.791,55.435,54.182,53.018,51.929,50.901,49.921,48.977,48.058,47.159,46.279,45.423,44.599,43.817,43.084,42.4  ,41.763,41.165,40.597,40.05 ,39.513,38.977,38.434,37.878,37.302,36.702,36.077,35.427,34.754,34.065,33.366,32.668,31.98 ,31.315,
                                                                     103.01,98.052,93.368,88.975,84.88 ,81.082,77.576,74.353,71.399,68.697,66.232,63.984,61.937,60.071,58.369,56.815,55.392,54.086,52.881,51.762,50.716,49.728,48.785,47.876,46.99 ,46.124,45.277,44.456,43.67 ,42.927,42.232,41.587,40.986,40.421,39.884,39.363,38.85 ,38.335,37.81 ,37.27 ,36.708,36.121,35.508,34.87 ,34.209,33.533,32.85 ,32.168,31.501,30.856,
                                                                     100.75,95.848,91.223,86.893,82.861,79.128,75.688,72.529,69.638,66.999,64.595,62.408,60.418,58.608,56.961,55.46 ,54.088,52.83 ,51.67 ,50.595,49.59 ,48.64 ,47.733,46.857,46.003,45.168,44.353,43.566,42.815,42.109,41.452,40.843,40.275,39.741,39.231,38.734,38.242,37.746,37.238,36.711,36.161,35.585,34.982,34.354,33.704,33.04 ,32.371,31.706,31.057,30.433,
                                                                     98.568,93.725,89.164,84.898,80.933,77.266,73.892,70.799,67.973,65.397,63.054,60.926,58.994,57.241,55.647,54.197,52.875,51.663,50.548,49.515,48.549,47.635,46.761,45.916,45.092,44.286,43.502,42.747,42.031,41.36 ,40.738,40.162,39.625,39.119,38.634,38.159,37.686,37.206,36.711,36.196,35.656,35.089,34.495,33.875,33.235,32.583,31.927,31.278,30.646,30.042,
                                                                     96.47 ,91.69 ,87.195,82.996,79.098,75.499,72.192,69.165,66.403,63.891,61.609,59.54 ,57.665,55.966,54.425,53.025,51.75 ,50.584,49.511,48.517,47.588,46.708,45.865,45.048,44.252,43.474,42.719,41.995,41.312,40.675,40.086,39.541,39.032,38.552,38.089,37.633,37.176,36.71 ,36.226,35.721,35.189,34.63 ,34.042,33.431,32.8  ,32.157,31.514,30.88 ,30.266,29.681,
                                                                     94.463,89.748,85.319,81.189,77.36 ,73.829,70.589,67.627,64.93 ,62.479,60.258,58.247,56.427,54.781,53.291,51.939,50.71 ,49.587,48.555,47.598,46.702,45.854,45.039,44.25 ,43.479,42.727,42    ,41.306,40.654,40.049,39.49 ,38.974,38.491,38.034,37.59 ,37.152,36.709,36.254,35.78 ,35.282,34.757,34.203,33.622,33.016,32.393,31.761,31.13 ,30.51 ,29.913,29.346,
                                                                     92.55 ,87.901,83.54 ,79.478,75.717,72.255,69.081,66.185,63.551,61.161,58.999,57.044,55.278,53.684,52.242,50.937,49.751,48.668,47.674,46.752,45.888,45.068,44.28 ,43.515,42.768,42.041,41.34 ,40.675,40.052,39.477,38.947,38.457,37.998,37.561,37.135,36.71 ,36.279,35.834,35.368,34.876,34.356,33.807,33.23 ,32.63 ,32.014,31.391,30.772,30.166,29.585,29.036,
                                                                     90.732,86.15 ,81.858,77.864,74.171,70.776,67.668,64.836,62.263,59.934,57.828,55.928,54.214,52.669,51.274,50.013,48.868,47.824,46.865,45.975,45.141,44.347,43.583,42.84 ,42.115,41.411,40.736,40.097,39.503,38.955,38.452,37.985,37.548,37.129,36.718,36.306,35.885,35.447,34.987,34.499,33.982,33.436,32.863,32.268,31.659,31.045,30.437,29.845,29.28 ,28.749,
                                                                     89.011,84.496,80.271,76.345,72.72 ,69.39 ,66.347,63.577,61.065,58.793,56.742,54.895,53.231,51.733,50.382,49.162,48.056,47.048,46.122,45.262,44.454,43.684,42.942,42.219,41.515,40.833,40.181,39.569,39.001,38.479,38    ,37.556,37.137,36.734,36.336,35.934,35.521,35.089,34.633,34.148,33.634,33.09 ,32.52 ,31.929,31.326,30.721,30.124,29.546,28.996,28.482,
                                                                     87.385,82.936,78.778,74.919,71.36 ,68.095,65.114,62.405,59.951,57.735,55.737,53.94 ,52.324,50.871,49.562,48.381,47.312,46.337,45.441,44.608,43.824,43.076,42.353,41.649,40.964,40.303,39.674,39.086,38.543,38.045,37.588,37.164,36.762,36.372,35.985,35.592,35.185,34.757,34.304,33.821,33.308,32.765,32.197,31.61 ,31.014,30.417,29.83 ,29.265,28.73 ,28.233,
                                                                     85.854,81.47 ,77.377,73.583,70.088,66.886,63.966,61.316,58.918,56.755,54.809,53.059,51.488,50.077,48.809,47.665,46.629,45.685,44.817,44.009,43.247,42.518,41.812,41.125,40.458,39.816,39.209,38.644,38.125,37.65 ,37.213,36.806,36.418,36.041,35.663,35.277,34.874,34.449,33.997,33.514,33.001,32.459,31.893,31.31 ,30.719,30.13 ,29.555,29.003,28.483,28.001,
                                                                     84.413,80.093,76.064,72.334,68.901,65.76 ,62.899,60.305,57.961,55.849,53.951,52.247,50.719,49.349,48.117,47.008,46.003,45.088,44.245,43.459,42.717,42.005,41.314,40.643,39.992,39.37 ,38.783,38.24 ,37.743,37.288,36.87 ,36.478,36.104,35.736,35.366,34.985,34.585,34.161,33.709,33.226,32.712,32.17 ,31.605,31.025,30.44 ,29.86 ,29.295,28.756,28.251,27.785,
                                                                     83.06 ,78.802,74.836,71.167,67.794,64.712,61.907,59.367,57.075,55.012,53.16 ,51.5  ,50.012,48.679,47.483,46.405,45.43 ,44.541,43.721,42.955,42.23 ,41.533,40.857,40.199,39.565,38.96 ,38.393,37.871,37.394,36.958,36.556,36.179,35.815,35.456,35.091,34.713,34.315,33.891,33.438,32.953,32.438,31.896,31.332,30.756,30.176,29.604,29.05 ,28.524,28.033,27.582,
                                                                     81.79 ,77.593,73.687,70.078,66.763,63.737,60.987,58.499,56.256,54.24 ,52.431,50.812,49.363,48.065,46.901,45.853,44.905,44.039,43.24 ,42.492,41.782,41.098,40.435,39.791,39.171,38.583,38.036,37.533,37.075,36.656,36.269,35.904,35.549,35.196,34.836,34.46 ,34.063,33.638,33.182,32.695,32.178,31.635,31.072,30.499,29.925,29.361,28.818,28.305,27.828,27.393,
                                                                     80.6  ,76.461,72.613,69.061,65.803,62.831,60.133,57.694,55.498,53.526,51.759,50.179,48.766,47.501,46.367,45.347,44.423,43.579,42.799,42.067,41.37 ,40.698,40.045,39.414,38.808,38.237,37.707,37.223,36.783,36.38 ,36.006,35.65 ,35.304,34.956,34.598,34.223,33.825,33.398,32.939,32.45 ,31.93 ,31.386,30.824,30.254,29.686,29.131,28.599,28.098,27.635,27.215,
                                                                     79.484,75.401,71.61 ,68.113,64.909,61.988,59.34 ,56.949,54.797,52.868,51.14 ,49.596,48.217,46.983,45.878,44.883,43.982,43.157,42.393,41.675,40.989,40.327,39.685,39.065,38.473,37.918,37.406,36.939,36.514,36.126,35.763,35.417,35.076,34.732,34.376,34    ,33.6  ,33.17 ,32.708,32.214,31.692,31.147,30.586,30.02 ,29.458,28.911,28.39 ,27.902,27.454,27.047,
                                                                     78.438,74.41 ,70.672,67.229,64.076,61.205,58.604,56.258,54.149,52.259,50.569,49.059,47.712,46.507,45.428,44.457,43.576,42.769,42.019,41.313,40.637,39.984,39.352,38.742,38.164,37.624,37.128,36.677,36.268,35.892,35.539,35.2  ,34.864,34.522,34.166,33.789,33.386,32.952,32.486,31.989,31.464,30.917,30.358,29.794,29.239,28.701,28.191,27.716,27.281,26.889,
                                                                     77.457,73.481,69.795,66.403,63.299,60.476,57.921,55.618,53.549,51.697,50.042,48.565,47.247,46.069,45.014,44.064,43.202,42.41 ,41.674,40.978,40.311,39.666,39.042,38.443,37.877,37.352,36.872,36.436,36.04 ,35.675,35.331,34.998,34.666,34.325,33.968,33.589,33.181,32.743,32.272,31.771,31.243,30.695,30.137,29.578,29.029,28.5  ,28.002,27.539,27.118,26.739,
                                                                     76.536,72.61 ,68.975,65.631,62.575,59.797,57.285,55.023,52.993,51.177,49.555,48.108,46.818,45.666,44.633,43.703,42.857,42.08 ,41.354,40.667,40.007,39.369,38.753,38.165,37.611,37.101,36.635,36.214,35.83 ,35.475,35.138,34.809,34.479,34.138,33.78 ,33.397,32.985,32.541,32.065,31.559,31.029,30.48 ,29.923,29.368,28.826,28.307,27.82 ,27.37 ,26.962,26.597,
                                                                     75.671,71.794,68.206,64.909,61.898,59.164,56.693,54.47 ,52.477,50.694,49.104,47.686,46.422,45.293,44.281,43.369,42.538,41.773,41.057,40.377,39.723,39.092,38.483,37.905,37.364,36.867,36.416,36.007,35.634,35.287,34.957,34.631,34.302,33.961,33.599,33.212,32.794,32.345,31.864,31.353,30.82 ,30.27 ,29.715,29.165,28.631,28.121,27.645,27.208,26.813,26.461,
                                                                     74.857,71.026,67.484,64.232,61.265,58.573,56.141,53.955,51.997,50.247,48.686,47.295,46.055,44.948,43.956,43.059,42.242,41.487,40.779,40.105,39.457,38.831,38.231,37.662,37.133,36.649,36.211,35.814,35.451,35.112,34.786,34.463,34.134,33.79 ,33.425,33.032,32.609,32.153,31.666,31.152,30.615,30.066,29.513,28.968,28.441,27.941,27.476,27.052,26.67 ,26.331,
                                                                     74.09 ,70.303,66.806,63.597,60.671,58.019,55.626,53.475,51.549,49.83 ,48.297,46.932,45.715,44.628,43.653,42.772,41.966,41.22 ,40.519,39.85 ,39.207,38.586,37.993,37.434,36.917,36.446,36.02 ,35.634,35.28 ,34.946,34.624,34.302,33.972,33.625,33.255,32.857,32.427,31.965,31.472,30.953,30.414,29.865,29.314,28.775,28.256,27.766,27.313,26.901,26.531,26.205,
                                                                     73.365,69.621,66.166,62.999,60.114,57.499,55.142,53.025,51.131,49.441,47.934,46.593,45.398,44.33 ,43.371,42.503,41.708,40.97 ,40.274,39.609,38.969,38.354,37.768,37.219,36.714,36.255,35.84 ,35.464,35.117,34.789,34.469,34.147,33.815,33.464,33.089,32.684,32.248,31.779,31.281,30.757,30.216,29.666,29.119,28.586,28.075,27.596,27.154,26.754,26.397,26.082,
                                                                     72.679,68.977,65.563,62.435,59.588,57.011,54.688,52.604,50.74 ,49.077,47.595,46.277,45.102,44.051,43.106,42.25 ,41.464,40.733,40.041,39.38 ,38.744,38.133,37.555,37.016,36.521,36.074,35.671,35.303,34.963,34.638,34.319,33.997,33.661,33.306,32.925,32.513,32.07 ,31.594,31.09 ,30.562,30.019,29.47 ,28.927,28.399,27.897,27.428,26.998,26.61 ,26.265,25.962,
                                                                     72.028,68.366,64.991,61.902,59.092,56.549,54.26 ,52.207,50.371,48.734,47.277,45.98 ,44.823,43.788,42.857,42.011,41.233,40.507,39.819,39.161,38.528,37.923,37.351,36.822,36.339,35.903,35.509,35.149,34.814,34.492,34.174,33.849,33.51 ,33.149,32.761,32.343,31.891,31.409,30.899,30.367,29.822,29.275,28.736,28.215,27.721,27.263,26.844,26.468,26.135,25.844,
                                                                     71.408,67.784,64.448,61.396,58.621,56.113,53.855,51.831,50.023,48.411,46.977,45.699,44.56 ,43.54 ,42.621,41.785,41.013,40.291,39.606,38.95 ,38.32 ,37.72 ,37.156,36.636,36.164,35.738,35.353,35    ,34.669,34.349,34.03 ,33.703,33.359,32.993,32.597,32.171,31.712,31.222,30.707,30.171,29.626,29.08 ,28.545,28.031,27.547,27.099,26.692,26.328,26.006,25.725,
                                                                     70.816,67.23 ,63.93 ,60.914,58.174,55.698,53.471,51.475,49.693,48.105,46.692,45.433,44.31 ,43.304,42.396,41.567,40.801,40.083,39.399,38.745,38.118,37.523,36.967,36.457,35.995,35.579,35.202,34.855,34.527,34.208,33.888,33.557,33.208,32.834,32.431,31.997,31.53 ,31.033,30.512,29.973,29.427,28.884,28.354,27.847,27.372,26.935,26.539,26.186,25.876,25.606,
                                                                     70.249,66.699,63.434,60.453,57.746,55.301,53.104,51.136,49.378,47.813,46.42 ,45.179,44.071,43.077,42.179,41.357,40.596,39.88 ,39.198,38.545,37.921,37.332,36.784,36.283,35.831,35.424,35.054,34.712,34.386,34.067,33.744,33.409,33.054,32.673,32.262,31.818,31.344,30.84 ,30.314,29.772,29.226,28.685,28.16 ,27.661,27.196,26.769,26.385,26.043,25.744,25.484,
                                                                     69.704,66.188,62.959,60.01 ,57.336,54.921,52.752,50.81 ,49.077,47.533,46.159,44.934,43.841,42.858,41.968,41.153,40.395,39.68 ,38.999,38.348,37.727,37.143,36.604,36.112,35.669,35.27 ,34.907,34.568,34.244,33.924,33.598,33.257,32.895,32.507,32.087,31.635,31.152,30.642,30.11 ,29.567,29.02 ,28.483,27.964,27.473,27.017,26.601,26.228,25.897,25.608,25.357,
                                                                     69.177,65.696,62.499,59.584,56.94 ,54.555,52.412,50.496,48.785,47.262,45.906,44.697,43.616,42.644,41.762,40.951,40.196,39.482,38.802,38.152,37.535,36.957,36.425,35.943,35.509,35.117,34.759,34.423,34.1  ,33.777,33.447,33.101,32.731,32.334,31.905,31.445,30.954,30.437,29.901,29.355,28.81 ,28.275,27.763,27.28 ,26.834,26.428,26.066,25.746,25.466,25.225,
                                                                     68.667,65.218,62.054,59.17 ,56.556,54.199,52.083,50.191,48.502,46.998,45.659,44.465,43.396,42.434,41.558,40.751,39.997,39.284,38.604,37.955,37.342,36.77 ,36.247,35.773,35.348,34.962,34.608,34.275,33.95 ,33.625,33.29 ,32.937,32.56 ,32.154,31.716,31.246,30.747,30.223,29.683,29.136,28.592,28.062,27.555,27.081,26.645,26.25 ,25.897,25.587,25.317,25.085,
                                                                     68.17 ,64.753,61.62 ,58.767,56.182,53.852,51.761,49.892,48.224,46.739,45.416,44.236,43.178,42.223,41.353,40.549,39.796,39.083,38.403,37.756,37.147,36.582,36.066,35.601,35.183,34.804,34.453,34.12 ,33.795,33.466,33.125,32.765,32.379,31.963,31.516,31.037,30.53 ,30    ,29.456,28.908,28.365,27.84 ,27.34 ,26.874,26.448,26.063,25.721,25.421,25.16 ,24.935,
                                                                     67.683,64.297,61.195,58.371,55.814,53.511,51.445,49.598,47.95 ,46.483,45.175,44.007,42.959,42.012,41.146,40.344,39.592,38.878,38.198,37.553,36.948,36.389,35.882,35.425,35.014,34.639,34.291,33.958,33.63 ,33.297,32.95 ,32.581,32.186,31.761,31.304,30.816,30.301,29.765,29.218,28.669,28.129,27.608,27.115,26.658,26.241,25.866,25.534,25.243,24.991,24.773,
                                                                     67.203,63.848,60.776,57.981,55.451,53.173,51.131,49.305,47.677,46.226,44.933,43.776,42.738,41.797,40.934,40.134,39.381,38.666,37.986,37.343,36.743,36.191,35.692,35.242,34.837,34.467,34.12 ,33.786,33.455,33.116,32.762,32.385,31.98 ,31.545,31.077,30.58 ,30.058,29.517,28.966,28.417,27.88 ,27.364,26.879,26.43 ,26.023,25.658,25.335,25.053,24.809,24.599,
                                                                     66.729,63.403,60.36 ,57.593,55.09 ,52.837,50.817,49.011,47.401,45.967,44.687,43.541,42.511,41.575,40.715,39.915,39.162,38.446,37.766,37.126,36.53 ,35.985,35.493,35.051,34.651,34.283,33.938,33.602,33.266,32.922,32.559,32.173,31.759,31.313,30.835,30.329,29.799,29.253,28.7  ,28.151,27.617,27.107,26.629,26.189,25.791,25.436,25.122,24.849,24.612,24.408,
                                                                     66.256,62.96 ,59.945,57.205,54.727,52.498,50.5  ,48.715,47.122,45.702,44.435,43.299,42.276,41.344,40.487,39.686,38.931,38.215,37.535,36.898,36.307,35.769,35.285,34.848,34.453,34.088,33.742,33.403,33.063,32.711,32.34 ,31.945,31.519,31.063,30.575,30.06 ,29.523,28.972,28.417,27.869,27.338,26.834,26.364,25.933,25.544,25.198,24.894,24.628,24.399,24.201,
                                                                     65.783,62.515,59.528,56.814,54.361,52.155,50.178,48.412,46.836,45.43 ,44.175,43.048,42.031,41.103,40.246,39.445,38.689,37.971,37.292,36.657,36.073,35.542,35.063,34.633,34.241,33.877,33.53 ,33.188,32.842,32.483,32.103,31.697,31.261,30.793,30.296,29.772,29.228,28.673,28.117,27.57 ,27.043,26.545,26.083,25.661,25.282,24.944,24.648,24.391,24.168,23.976,
                                                                     65.306,62.066,59.105,56.417,53.989,51.805,49.849,48.1  ,46.54 ,45.148,43.903,42.785,41.773,40.848,39.991,39.189,38.431,37.712,37.034,36.403,35.824,35.3  ,34.828,34.402,34.014,33.65 ,33.301,32.955,32.602,32.235,31.845,31.428,30.981,30.503,29.995,29.463,28.913,28.354,27.797,27.253,26.731,26.239,25.785,25.372,25.001,24.673,24.386,24.135,23.919,23.733,
                                                                     64.823,61.61 ,58.676,56.013,53.608,51.446,49.509,47.778,46.233,44.854,43.619,42.508,41.5  ,40.577,39.72 ,38.916,38.156,37.436,36.76 ,36.133,35.56 ,35.042,34.576,34.155,33.768,33.404,33.052,32.701,32.342,31.965,31.566,31.138,30.68 ,30.191,29.674,29.133,28.577,28.016,27.459,26.917,26.4  ,25.916,25.47 ,25.066,24.704,24.385,24.105,23.862,23.652,23.471,
                                                                     64.332,61.145,58.236,55.598,53.216,51.075,49.157,47.443,45.913,44.545,43.319,42.215,41.211,40.289,39.431,38.625,37.863,37.143,36.469,35.846,35.279,34.767,34.307,33.89 ,33.504,33.139,32.784,32.427,32.06 ,31.674,31.264,30.826,30.356,29.857,29.33 ,28.782,28.221,27.657,27.101,26.563,26.051,25.575,25.138,24.743,24.391,24.08 ,23.808,23.572,23.369,23.193,
                                                                     63.829,60.668,57.785,55.17 ,52.811,50.69 ,48.791,47.093,45.576,44.22 ,43.003,41.904,40.903,39.982,39.123,38.314,37.55 ,36.831,36.159,35.541,34.98 ,34.475,34.02 ,33.606,33.221,32.854,32.495,32.132,31.757,31.362,30.941,30.492,30.011,29.501,28.966,28.411,27.846,27.281,26.726,26.192,25.687,25.219,24.791,24.405,24.062,23.76 ,23.497,23.268,23.07 ,22.9};

            return arrayCurvatureMax[y * arrayWidth + x];
        }
        case 8:
        {
            static const std::vector<double> arrayCurvatureMax =    {132.25,127.87,123.73,119.79,115.99,112.31,108.74,105.27,101.9 ,98.625,95.456,92.387,89.419,86.554,83.79 ,81.127,78.565,76.102,73.737,71.47 ,69.297,67.218,65.229,63.329,61.515,59.785,58.136,56.565,55.071,53.651,52.304,51.036,49.87 ,48.875,48.212,47.997,47.884,47.453,46.765,45.992,45.221,44.479,43.774,43.106,42.474,41.876,41.311,40.776,40.268,39.786,
                                                                     129.23,124.76,120.59,116.66,112.9 ,109.28,105.78,102.38,99.092,95.905,92.821,89.838,86.957,84.178,81.501,78.924,76.447,74.069,71.788,69.602,67.509,65.507,63.595,61.769,60.027,58.367,56.785,55.28 ,53.848,52.489,51.203,49.998,48.909,48.031,47.545,47.453,47.288,46.779,46.074,45.321,44.58 ,43.871,43.198,42.561,41.959,41.389,40.849,40.337,39.85 ,39.384,
                                                                     126.36,121.78,117.56,113.62,109.9 ,106.33,102.89,99.572,96.363,93.259,90.26 ,87.364,84.57 ,81.877,79.286,76.795,74.403,72.108,69.908,67.802,65.788,63.863,62.024,60.271,58.599,57.006,55.49 ,54.048,52.678,51.378,50.151,49.011,48.007,47.265,46.959,46.932,46.684,46.116,45.408,44.678,43.968,43.291,42.65 ,42.042,41.467,40.923,40.406,39.914,39.443,38.989,
                                                                     123.65,118.95,114.66,110.69,106.98,103.46,100.08,96.838,93.709,90.688,87.774,84.964,82.257,79.651,77.145,74.739,72.43 ,70.216,68.097,66.07 ,64.132,62.282,60.516,58.833,57.229,55.702,54.249,52.868,51.557,50.316,49.148,48.076,47.168,46.581,46.439,46.415,46.075,45.467,44.767,44.063,43.384,42.738,42.126,41.547,40.997,40.476,39.979,39.503,39.044,38.594,
                                                                     121.1 ,116.29,111.89,107.88,104.16,100.67,97.356,94.182,91.131,88.193,85.364,82.64 ,80.018,77.498,75.077,72.754,70.527,68.394,66.354,64.404,62.541,60.763,59.068,57.453,55.915,54.452,53.061,51.74 ,50.486,49.301,48.193,47.193,46.397,45.981,45.965,45.89 ,45.467,44.838,44.153,43.476,42.827,42.211,41.626,41.072,40.546,40.044,39.563,39.098,38.642,38.186,
                                                                     118.69,113.8 ,109.29,105.2 ,101.45,97.984,94.716,91.607,88.632,85.775,83.03 ,80.39 ,77.853,75.417,73.08 ,70.839,68.693,66.64 ,64.677,62.801,61.012,59.305,57.679,56.13 ,54.656,53.255,51.924,50.66 ,49.462,48.333,47.285,46.366,45.699,45.461,45.514,45.354,44.866,44.23 ,43.565,42.915,42.295,41.706,41.147,40.616,40.11 ,39.624,39.153,38.69 ,38.226,37.749,
                                                                     116.37,111.47,106.86,102.67,98.875,95.402,92.168,89.117,86.212,83.434,80.771,78.215,75.762,73.409,71.154,68.994,66.928,64.952,63.064,61.262,59.544,57.906,56.347,54.862,53.451,52.109,50.836,49.627,48.484,47.41 ,46.426,45.598,45.081,45.007,45.065,44.81 ,44.277,43.644,43.002,42.38 ,41.787,41.223,40.687,40.176,39.685,39.208,38.738,38.266,37.779,37.262,
                                                                     114.09,109.27,104.61,100.31,96.44 ,92.938,89.72 ,86.715,83.875,81.172,78.588,76.114,73.743,71.473,69.298,67.218,65.229,63.328,61.514,59.784,58.135,56.565,55.07 ,53.648,52.297,51.013,49.795,48.64 ,47.55 ,46.532,45.616,44.895,44.545,44.597,44.605,44.263,43.704,43.082,42.463,41.867,41.3  ,40.759,40.242,39.745,39.263,38.786,38.306,37.809,37.28 ,36.705,
                                                                     111.83,107.17,102.52,98.137,94.166,90.607,87.382,84.408,81.624,78.989,76.482,74.087,71.797,69.606,67.511,65.508,63.595,61.769,60.026,58.366,56.784,55.279,53.847,52.485,51.192,49.965,48.8  ,47.698,46.66 ,45.699,44.859,44.262,44.085,44.21 ,44.133,43.72 ,43.149,42.542,41.947,41.376,40.831,40.309,39.807,39.318,38.834,38.345,37.838,37.296,36.708,36.065,
                                                                     109.54,105.1 ,100.56,96.137,92.068,88.426,85.167,82.204,79.462,76.889,74.452,72.134,69.921,67.809,65.791,63.864,62.025,60.27 ,58.598,57.005,55.489,54.047,52.676,51.373,50.136,48.962,47.85 ,46.798,45.812,44.91 ,44.157,43.704,43.69 ,43.823,43.648,43.185,42.613,42.025,41.452,40.903,40.376,39.868,39.373,38.882,38.384,37.866,37.312,36.709,36.053,35.353,
                                                                     107.21,103.04,98.675,94.295,90.151,86.41 ,83.088,80.112,77.395,74.873,72.501,70.254,68.117,66.08 ,64.137,62.284,60.517,58.833,57.228,55.701,54.248,52.867,51.555,50.309,49.126,48.004,46.942,45.94 ,45.006,44.168,43.516,43.225,43.34 ,43.425,43.158,42.662,42.096,41.527,40.975,40.443,39.93 ,39.428,38.929,38.422,37.893,37.326,36.709,36.039,35.33 ,34.611,
                                                                     104.85,100.95,96.829,92.572,88.409,84.573,81.163,78.144,75.431,72.945,70.629,68.448,66.382,64.418,62.548,60.767,59.07 ,57.453,55.915,54.452,53.06 ,51.738,50.483,49.291,48.16 ,47.089,46.076,45.123,44.241,43.475,42.941,42.82 ,43.011,43.01 ,42.666,42.153,41.597,41.045,40.51 ,39.991,39.483,38.976,38.46 ,37.919,37.338,36.707,36.025,35.306,34.584,33.899,
                                                                     102.47,98.831,94.978,90.921,86.817,82.917,79.406,76.316,73.58 ,71.111,68.839,66.717,64.718,62.823,61.023,59.311,57.681,56.131,54.657,53.255,51.923,50.658,49.457,48.318,47.238,46.215,45.249,44.344,43.519,42.834,42.439,42.479,42.681,42.579,42.178,41.658,41.114,40.577,40.053,39.538,39.023,38.497,37.944,37.35 ,36.704,36.009,35.282,34.558,33.877,33.27 ,
                                                                     100.09,96.682,93.1  ,89.295,85.334,81.431,77.827,74.642,71.854,69.38 ,67.136,65.063,63.123,61.293,59.56 ,57.914,56.351,54.864,53.452,52.109,50.835,49.625,48.477,47.389,46.357,45.381,44.461,43.605,42.839,42.25 ,42.012,42.182,42.337,42.137,41.697,41.176,40.642,40.114,39.592,39.069,38.533,37.968,37.36 ,36.7  ,35.993,35.259,34.534,33.857,33.258,32.748,
                                                                     97.733,94.519,91.19 ,87.658,83.911,80.085,76.427,73.137,70.268,67.76 ,65.525,63.487,61.599,59.829,58.159,56.577,55.076,53.651,52.298,51.013,49.794,48.637,47.54 ,46.501,45.516,44.586,43.711,42.903,42.205,41.729,41.656,41.905,41.973,41.689,41.222,40.702,40.174,39.647,39.115,38.568,37.991,37.369,36.695,35.976,35.235,34.51 ,33.839,33.249,32.747,32.328,
                                                                     95.409,92.36 ,89.255,85.99 ,82.5  ,78.835,75.189,71.81 ,68.837,66.265,64.013,61.995,60.148,58.43 ,56.818,55.296,53.856,52.49 ,51.194,49.965,48.799,47.694,46.646,45.653,44.714,43.827,42.997,42.241,41.619,41.276,41.362,41.625,41.586,41.236,40.752,40.23 ,39.7  ,39.16 ,38.603,38.013,37.376,36.689,35.958,35.212,34.487,33.823,33.241,32.749,32.336,31.992,
                                                                     93.142,90.222,87.308,84.289,81.071,77.629,74.078,70.654,67.574,64.91 ,62.613,60.592,58.771,57.098,55.538,54.073,52.689,51.38 ,50.139,48.963,47.849,46.792,45.792,44.844,43.948,43.104,42.318,41.618,41.086,40.893,41.11 ,41.325,41.18 ,40.779,40.28 ,39.75 ,39.204,38.636,38.033,37.383,36.681,35.94 ,35.189,34.466,33.809,33.236,32.752,32.347,32.007,31.718,
                                                                     90.959,88.123,85.368,82.563,79.606,76.425,73.047,69.65 ,66.484,63.71 ,61.336,59.287,57.475,55.833,54.319,52.904,51.574,50.318,49.131,48.006,46.941,45.931,44.976,44.072,43.217,42.415,41.675,41.037,40.611,40.579,40.874,40.997,40.756,40.313,39.796,39.246,38.669,38.052,37.388,36.673,35.922,35.167,34.447,33.796,33.233,32.757,32.359,32.024,31.738,31.489,
                                                                     88.897,86.081,83.45 ,80.828,78.107,75.193,72.047,68.76 ,65.56 ,62.675,60.198,58.09 ,56.264,54.64 ,53.161,51.792,50.511,49.305,48.168,47.092,46.073,45.11 ,44.197,43.335,42.521,41.759,41.066,40.499,40.199,40.321,40.629,40.636,40.311,39.83 ,39.284,38.699,38.07 ,37.391,36.663,35.903,35.145,34.429,33.786,33.231,32.764,32.372,32.042,31.758,31.511,31.292,
                                                                     87.018,84.121,81.571,79.1  ,76.584,73.922,71.035,67.933,64.774,61.808,59.212,57.016,55.149,53.523,52.067,50.735,49.498,48.339,47.248,46.219,45.246,44.325,43.455,42.632,41.856,41.135,40.493,40.009,39.851,40.099,40.354,40.238,39.839,39.314,38.726,38.086,37.393,36.653,35.885,35.125,34.413,33.778,33.232,32.773,32.387,32.06 ,31.78 ,31.534,31.315,31.117,
                                                                     85.416,82.276,79.748,77.394,75.051,72.614,69.988,67.122,64.087,61.096,58.389,56.079,54.141,52.489,51.04 ,49.735,48.535,47.419,46.372,45.387,44.456,43.577,42.745,41.96 ,41.222,40.541,39.955,39.568,39.562,39.883,40.033,39.797,39.327,38.746,38.099,37.394,36.641,35.866,35.105,34.398,33.771,33.235,32.783,32.403,32.08 ,31.802,31.557,31.338,31.139,30.958,
                                                                     84.221,80.603,78.002,75.726,73.524,71.28 ,68.896,66.287,63.446,60.511,57.728,55.293,53.253,51.547,50.085,48.794,47.624,46.544,45.538,44.593,43.703,42.862,42.068,41.319,40.617,39.977,39.453,39.181,39.317,39.641,39.653,39.302,38.755,38.107,37.392,36.629,35.848,35.087,34.385,33.767,33.239,32.795,32.42 ,32.101,31.824,31.58 ,31.361,31.163,30.98 ,30.812,
                                                                     83.534,79.192,76.365,74.111,72.019,69.937,67.765,65.408,62.806,60.007,57.212,54.665,52.502,50.711,49.211,47.916,46.765,45.715,44.744,43.837,42.984,42.18 ,41.421,40.706,40.038,39.441,38.988,38.843,39.091,39.344,39.202,38.736,38.106,37.386,36.615,35.83 ,35.07 ,34.375,33.765,33.246,32.809,32.439,32.123,31.848,31.604,31.385,31.186,31.003,30.834,30.676,
                                                                     83.306,78.169,74.89 ,72.569,70.549,68.599,66.607,64.481,62.131,59.534,56.808,54.192,51.898,49.994,48.427,47.108,45.961,44.933,43.992,43.118,42.299,41.529,40.802,40.117,39.481,38.928,38.558,38.544,38.846,38.966,38.663,38.086,37.375,36.6  ,35.811,35.055,34.366,33.765,33.254,32.824,32.459,32.146,31.872,31.629,31.41 ,31.21 ,31.026,30.855,30.697,30.549,
                                                                     83.252,77.644,73.66 ,71.127,69.129,67.28 ,65.439,63.511,61.404,59.046,56.468,53.853,51.449,49.411,47.746,46.377,45.217,44.198,43.279,42.434,41.645,40.905,40.206,39.548,38.943,38.437,38.158,38.263,38.541,38.485,38.027,37.35 ,36.581,35.793,35.04 ,34.359,33.767,33.265,32.84 ,32.48 ,32.169,31.897,31.654,31.435,31.234,31.049,30.878,30.718,30.569,30.429,
                                                                     83.007,77.592,72.797,69.836,67.775,65.993,64.273,62.514,60.624,58.512,56.142,53.609,51.145,48.974,47.184,45.735,44.538,43.515,42.609,41.784,41.021,40.306,39.63 ,38.995,38.416,37.959,37.778,37.963,38.138,37.89 ,37.298,36.553,35.773,35.027,34.354,33.772,33.277,32.859,32.502,32.194,31.922,31.68 ,31.46 ,31.259,31.073,30.9  ,30.739,30.589,30.448,30.316,
                                                                     82.386,77.754,72.42 ,68.774,66.514,64.751,63.125,61.504,59.8  ,57.919,55.786,53.413,50.963,48.684,46.753,45.194,43.934,42.886,41.98 ,41.169,40.423,39.727,39.068,38.447,37.889,37.484,37.398,37.603,37.614,37.19 ,36.507,35.748,35.014,34.351,33.778,33.291,32.878,32.525,32.219,31.948,31.706,31.486,31.283,31.097,30.923,30.761,30.609,30.467,30.334,30.208,
                                                                     81.42 ,77.759,72.526,68.059,65.392,63.57 ,62.005,60.495,58.947,57.271,55.376,53.214,50.861,48.532,46.465,44.769,43.416,42.319,41.396,40.585,39.848,39.161,38.509,37.894,37.351,36.997,36.991,37.147,36.974,36.422,35.713,34.999,34.349,33.785,33.306,32.899,32.549,32.245,31.975,31.732,31.512,31.309,31.121,30.946,30.782,30.63 ,30.486,30.352,30.226,30.107,
                                                                     80.234,77.391,72.882,67.813,64.483,62.474,60.925,59.5  ,58.079,56.576,54.902,52.974,50.788,48.487,46.318,44.474,42.996,41.821,40.858,40.032,39.29 ,38.6  ,37.941,37.318,36.783,36.486,36.533,36.587,36.258,35.651,34.977,34.346,33.794,33.323,32.921,32.575,32.272,32.002,31.76 ,31.538,31.334,31.145,30.969,30.804,30.65 ,30.506,30.371,30.243,30.124,30.011,
                                                                     78.941,76.663,73.118,68.057,63.901,61.506,59.9  ,58.53 ,57.21 ,55.851,54.367,52.669,50.697,48.506,46.297,44.315,42.688,41.402,40.371,39.507,38.741,38.029,37.346,36.701,36.175,35.943,36.016,35.954,35.533,34.938,34.34 ,33.803,33.341,32.944,32.601,32.299,32.03 ,31.787,31.565,31.36 ,31.17 ,30.992,30.827,30.671,30.526,30.389,30.261,30.14 ,30.027,29.92 ,
                                                                     77.614,75.693,72.987,68.586,63.77 ,60.735,58.953,57.596,56.352,55.109,53.786,52.295,50.552,48.536,46.367,44.287,42.501,41.071,39.937,39.005,38.189,37.432,36.705,36.03 ,35.527,35.38 ,35.461,35.307,34.861,34.323,33.81 ,33.359,32.968,32.628,32.328,32.059,31.815,31.592,31.386,31.194,31.016,30.849,30.693,30.546,30.408,30.279,30.157,30.042,29.935,29.833,
                                                                     76.296,74.599,72.481,69.031,64.132,60.27 ,58.122,56.711,55.516,54.363,53.171,51.858,50.334,48.53 ,46.48 ,44.367,42.437,40.835,39.559,38.519,37.62 ,36.79 ,36.004,35.308,34.867,34.833,34.912,34.705,34.28 ,33.811,33.376,32.993,32.656,32.356,32.088,31.844,31.62 ,31.412,31.22 ,31.04 ,30.872,30.714,30.566,30.427,30.296,30.174,30.058,29.95 ,29.847,29.751,
                                                                     75.009,73.459,71.713,69.121,64.81 ,60.238,57.473,55.895,54.709,53.625,52.537,51.371,50.04 ,48.455,46.58 ,44.513,42.479,40.695,39.234,38.042,37.02 ,36.094,35.248,34.567,34.249,34.341,34.409,34.18 ,33.794,33.389,33.017,32.684,32.386,32.117,31.873,31.647,31.439,31.245,31.064,30.895,30.736,30.586,30.446,30.315,30.191,30.074,29.965,29.862,29.765,29.674,
                                                                     73.764,72.318,70.802,68.822,65.444,60.697,57.11 ,55.185,53.946,52.903,51.896,50.846,49.676,48.294,46.621,44.668,42.594,40.638,38.959,37.566,36.386,35.354,34.475,33.87 ,33.731,33.932,33.974,33.739,33.392,33.038,32.712,32.416,32.147,31.902,31.676,31.466,31.271,31.088,30.918,30.757,30.607,30.466,30.333,30.208,30.09 ,29.98 ,29.876,29.778,29.686,29.6  ,
                                                                     72.565,71.201,69.832,68.239,65.737,61.502,57.16 ,54.641,53.244,52.206,51.256,50.296,49.253,48.041,46.564,44.769,42.723,40.633,38.723,37.096,35.739,34.614,33.75 ,33.287,33.346,33.606,33.607,33.369,33.054,32.739,32.446,32.178,31.932,31.704,31.493,31.297,31.113,30.941,30.779,30.628,30.485,30.351,30.225,30.107,29.995,29.89 ,29.792,29.699,29.612,29.53 ,
                                                                     71.415,70.119,68.852,67.495,65.629,62.302,57.696,54.363,52.635,51.542,50.623,49.727,48.778,47.698,46.388,44.758,42.795,40.633,38.52 ,36.659,35.131,33.945,33.146,32.865,33.092,33.345,33.294,33.054,32.763,32.476,32.209,31.962,31.734,31.521,31.323,31.138,30.964,30.802,30.649,30.505,30.37 ,30.243,30.123,30.011,29.905,29.806,29.712,29.624,29.541,29.463,
                                                                     70.312,69.077,67.887,66.678,65.215,62.778,58.603,54.475,52.175,50.926,50.001,49.143,48.256,47.264,46.078,44.594,42.745,40.589,38.35 ,36.305,34.64 ,33.423,32.715,32.616,32.939,33.129,33.024,32.78 ,32.507,32.241,31.994,31.763,31.549,31.35 ,31.163,30.988,30.824,30.67 ,30.525,30.388,30.26 ,30.14 ,30.026,29.92 ,29.819,29.725,29.636,29.553,29.474,29.4  ,
                                                                     69.256,68.078,66.95 ,65.84 ,64.621,62.843,59.542,55.064,51.955,50.384,49.392,48.54 ,47.682,46.737,45.625,44.251,42.532,40.47 ,38.227,36.096,34.338,33.1  ,32.482,32.525,32.85 ,32.942,32.786,32.538,32.276,32.027,31.795,31.579,31.377,31.189,31.012,30.847,30.691,30.545,30.407,30.278,30.156,30.042,29.935,29.833,29.738,29.649,29.564,29.485,29.41 ,29.34 ,
                                                                     68.245,67.12 ,66.047,65.009,63.937,62.578,60.177,56.037,52.095,49.963,48.802,47.912,47.046,46.108,45.029,43.738,42.164,40.28 ,38.17 ,36.066,34.264,33.003,32.451,32.562,32.801,32.778,32.574,32.32 ,32.067,31.83 ,31.61 ,31.406,31.215,31.037,30.87 ,30.713,30.565,30.426,30.296,30.173,30.058,29.949,29.847,29.751,29.661,29.576,29.496,29.421,29.35 ,29.284,
                                                                     67.278,66.203,65.18 ,64.198,63.22 ,62.109,60.385,57.074,52.685,49.739,48.243,47.247,46.335,45.375,44.308,43.094,41.689,40.049,38.177,36.206,34.416,33.135,32.614,32.695,32.779,32.636,32.386,32.123,31.875,31.648,31.438,31.244,31.063,30.893,30.735,30.586,30.446,30.314,30.19 ,30.074,29.965,29.862,29.765,29.674,29.588,29.507,29.431,29.36 ,29.293,29.23 ,
                                                                     66.353,65.324,64.346,63.411,62.497,61.529,60.233,57.814,53.653,49.818,47.744,46.54 ,45.545,44.551,43.506,42.392,41.18 ,39.81 ,38.227,36.469,34.76 ,33.489,32.959,32.904,32.785,32.52 ,32.223,31.947,31.7  ,31.478,31.276,31.091,30.918,30.758,30.607,30.465,30.333,30.208,30.09 ,29.98 ,29.876,29.778,29.686,29.6  ,29.518,29.442,29.37 ,29.302,29.238,29.179,
                                                                     65.467,64.482,63.545,62.648,61.778,60.893,59.837,58.096,54.688,50.289,47.372,45.809,44.694,43.68 ,42.69 ,41.705,40.693,39.582,38.289,36.789,35.243,34.04 ,33.463,33.183,32.832,32.439,32.087,31.791,31.54 ,31.321,31.126,30.947,30.782,30.629,30.486,30.351,30.225,30.107,29.995,29.89 ,29.792,29.699,29.612,29.53 ,29.452,29.38 ,29.312,29.247,29.187,29.13 ,
                                                                     64.617,63.673,62.771,61.906,61.065,60.221,59.288,57.96 ,55.406,51.099,47.239,45.113,43.836,42.829,41.929,41.087,40.256,39.367,38.334,37.11 ,35.806,34.748,34.096,33.536,32.936,32.405,31.985,31.659,31.397,31.177,30.985,30.812,30.654,30.508,30.371,30.243,30.123,30.011,29.905,29.806,29.712,29.624,29.541,29.463,29.39 ,29.321,29.256,29.195,29.138,29.084,
                                                                     63.801,62.891,62.019,61.175,60.349,59.518,58.629,57.51 ,55.612,51.973,47.477,44.576,43.056,42.062,41.263,40.554,39.875,39.162,38.349,37.396,36.393,35.551,34.82 ,33.97 ,33.119,32.431,31.925,31.555,31.274,31.047,30.855,30.686,30.533,30.393,30.263,30.141,30.027,29.92 ,29.819,29.725,29.636,29.553,29.474,29.4  ,29.331,29.265,29.204,29.146,29.092,29.04 ,
                                                                     63.011,62.131,61.279,60.445,59.616,58.772,57.877,56.838,55.344,52.564,48.125,44.379,42.453,41.424,40.706,40.103,39.539,38.963,38.331,37.636,36.966,36.377,35.587,34.489,33.398,32.535,31.918,31.486,31.173,30.932,30.737,30.569,30.42 ,30.285,30.16 ,30.044,29.935,29.834,29.738,29.649,29.564,29.485,29.41 ,29.34 ,29.274,29.212,29.154,29.099,29.048,28.999,
                                                                     62.243,61.381,60.537,59.698,58.846,57.965,57.034,56.013,54.752,52.717,49.003,44.691,42.13 ,40.94 ,40.249,39.718,39.24 ,38.769,38.288,37.834,37.499,37.154,36.345,35.082,33.784,32.732,31.978,31.46 ,31.099,30.836,30.631,30.461,30.315,30.183,30.063,29.952,29.849,29.752,29.661,29.576,29.496,29.421,29.35 ,29.284,29.221,29.162,29.107,29.055,29.006,28.96 ,
                                                                     61.484,60.629,59.777,58.912,58.018,57.086,56.114,55.103,54.001,52.512,49.788,45.552,42.199,40.634,39.885,39.388,38.969,38.58 ,38.232,38.007,37.979,37.827,37.044,35.72 ,34.277,33.038,32.12 ,31.489,31.061,30.762,30.54 ,30.364,30.217,30.088,29.972,29.866,29.767,29.675,29.588,29.507,29.431,29.36 ,29.293,29.23 ,29.17 ,29.114,29.062,29.013,28.966,28.922,
                                                                     60.721,59.856,58.976,58.066,57.12 ,56.142,55.155,54.185,53.225,52.117,50.264,46.742,42.766,40.557,39.615,39.104,38.722,38.402,38.178,38.172,38.396,38.367,37.647,36.362,34.856,33.459,32.362,31.588,31.068,30.716,30.467,30.28 ,30.129,30    ,29.887,29.784,29.689,29.601,29.519,29.442,29.37 ,29.302,29.239,29.179,29.122,29.069,29.019,28.973,28.928,28.887,
                                                                     59.934,59.039,58.113,57.15 ,56.161,55.17 ,54.216,53.33 ,52.508,51.661,50.421,47.873,43.846,40.8  ,39.463,38.867,38.499,38.238,38.139,38.342,38.742,38.769,38.133,36.962,35.486,33.986,32.714,31.773,31.133,30.707,30.417,30.21 ,30.05 ,29.92 ,29.807,29.707,29.616,29.532,29.454,29.38 ,29.312,29.247,29.187,29.13 ,29.077,29.026,28.979,28.935,28.893,28.854,
                                                                     59.1  ,58.158,57.18 ,56.178,55.182,54.23 ,53.357,52.578,51.883,51.212,50.356,48.663,45.237,41.475,39.488,38.691,38.303,38.098,38.131,38.523,39.011,39.044,38.499,37.485,36.118,34.595,33.18 ,32.06 ,31.27 ,30.744,30.396,30.158,29.984,29.847,29.734,29.636,29.547,29.467,29.392,29.322,29.257,29.196,29.138,29.084,29.033,28.986,28.941,28.899,28.859,28.822};
        }
        case 9:
        {
            static const std::vector<double> arrayCurvatureMax =    {133.28,131.02,127.8 ,123.46,118.17,112.57,107.41,103.09,99.58 ,96.647,94.062,91.678,89.425,87.281,85.244,83.318,81.503,79.794,78.18 ,76.647,75.183,73.774,72.408,71.078,69.776,68.497,67.238,65.997,64.773,63.565,62.375,61.202,60.047,58.913,57.798,56.705,55.635,54.586,53.561,52.558,51.576,50.615,49.672,48.745,47.83 ,46.922,46.017,45.107,44.186,43.247,
                                                                     132.1 ,129.75,126.45,122.02,116.66,111.02,105.85,101.53,98.015,95.062,92.446,90.024,87.729,85.54 ,83.458,81.488,79.633,77.886,76.238,74.678,73.19 ,71.762,70.383,69.045,67.74 ,66.462,65.21 ,63.98 ,62.771,61.584,60.417,59.272,58.148,57.048,55.97 ,54.917,53.888,52.882,51.901,50.943,50.006,49.09 ,48.19 ,47.304,46.428,45.556,44.682,43.8  ,42.901,41.981,
                                                                     130.82,128.4 ,125   ,120.49,115.06,109.39,104.21,99.891,96.367,93.395,90.752,88.296,85.963,83.734,81.612,79.604,77.713,75.935,74.261,72.679,71.175,69.736,68.351,67.012,65.711,64.443,63.205,61.993,60.807,59.646,58.509,57.397,56.311,55.249,54.213,53.203,52.218,51.259,50.323,49.411,48.52 ,47.648,46.791,45.945,45.105,44.266,43.421,42.562,41.683,40.777,
                                                                     129.45,126.95,123.47,118.86,113.37,107.67,102.49,98.169,94.641,91.654,88.986,86.501,84.134,81.871,79.715,77.675,75.755,73.953,72.26 ,70.663,69.15 ,67.708,66.326,64.993,63.704,62.453,61.235,60.049,58.892,57.763,56.662,55.588,54.542,53.524,52.533,51.568,50.63 ,49.718,48.83 ,47.964,47.118,46.289,45.473,44.665,43.859,43.05 ,42.229,41.389,40.525,39.632,
                                                                     128   ,125.41,121.84,117.15,111.59,105.86,100.68,96.372,92.843,89.844,87.157,84.647,82.253,79.961,77.778,75.713,73.772,71.953,70.247,68.644,67.13 ,65.692,64.319,63.001,61.731,60.503,59.313,58.157,57.035,55.944,54.883,53.852,52.851,51.878,50.934,50.018,49.128,48.263,47.422,46.603,45.801,45.013,44.235,43.462,42.687,41.902,41.102,40.278,39.426,38.541,
                                                                     126.46,123.79,120.13,115.35,109.74,103.98,98.807,94.507,90.98 ,87.975,85.274,82.745,80.33 ,78.017,75.814,73.731,71.776,69.947,68.237,66.635,65.127,63.701,62.344,61.048,59.803,58.605,57.448,56.329,55.246,54.196,53.18 ,52.195,51.24 ,50.316,49.421,48.553,47.712,46.895,46.101,45.325,44.566,43.817,43.075,42.332,41.583,40.82 ,40.035,39.223,38.379,37.502,
                                                                     124.83,122.08,118.34,113.48,107.8 ,102.04,96.868,92.582,89.063,86.057,83.348,80.807,78.378,76.051,73.835,71.742,69.78 ,67.95 ,66.243,64.649,63.155,61.747,60.414,59.145,57.932,56.768,55.65 ,54.572,53.532,52.528,51.558,50.621,49.716,48.841,47.995,47.177,46.384,45.614,44.864,44.131,43.411,42.698,41.987,41.272,40.545,39.798,39.025,38.221,37.382,36.51 ,
                                                                     123.12,120.3 ,116.47,111.53,105.81,100.03,94.875,90.609,87.104,84.103,81.393,78.847,76.411,74.078,71.856,69.761,67.8  ,65.975,64.279,62.701,61.227,59.844,58.54 ,57.304,56.127,55.004,53.927,52.894,51.9  ,50.944,50.024,49.136,48.281,47.456,46.659,45.889,45.142,44.418,43.71 ,43.017,42.333,41.653,40.97 ,40.277,39.566,38.832,38.066,37.266,36.432,35.566,
                                                                     121.34,118.44,114.53,109.52,103.76,97.974,92.842,88.602,85.117,82.128,79.423,76.879,74.445,72.112,69.893,67.803,65.851,64.039,62.36 ,60.803,59.356,58.003,56.734,55.536,54.4  ,53.32 ,52.289,51.303,50.358,49.451,48.58 ,47.743,46.938,46.162,45.413,44.69 ,43.989,43.306,42.639,41.983,41.331,40.679,40.018,39.343,38.645,37.918,37.156,36.358,35.526,34.667,
                                                                     119.5 ,116.53,112.55,107.47,101.67,95.891,90.786,86.579,83.12 ,80.149,77.458,74.924,72.498,70.174,67.965,65.887,63.95 ,62.157,60.501,58.972,57.556,56.239,55.008,53.852,52.76 ,51.726,50.743,49.806,48.91 ,48.054,47.233,46.446,45.69 ,44.962,44.26 ,43.582,42.923,42.28 ,41.65 ,41.026,40.402,39.773,39.131,38.468,37.777,37.053,36.293,35.495,34.666,33.816,
                                                                     117.63,114.59,110.54,105.4 ,99.574,93.801,88.731,84.563,81.137,78.192,75.52 ,73.004,70.594,68.285,66.093,64.033,62.118,60.35 ,58.723,57.226,55.846,54.568,53.379,52.267,51.221,50.235,49.301,48.413,47.568,46.761,45.99 ,45.251,44.543,43.861,43.204,42.567,41.948,41.342,40.743,40.148,39.548,38.937,38.308,37.652,36.965,36.24 ,35.478,34.68 ,33.856,33.018,
                                                                     115.75,112.65,108.53,103.34,97.493,91.738,86.71 ,82.588,79.201,76.287,73.643,71.15 ,68.762,66.476,64.306,62.271,60.382,58.644,57.049,55.589,54.247,53.011,51.866,50.8  ,49.802,48.864,47.978,47.14 ,46.344,45.586,44.863,44.171,43.508,42.87 ,42.254,41.656,41.072,40.498,39.928,39.355,38.774,38.176,37.554,36.902,36.213,35.485,34.719,33.921,33.104,32.282,
                                                                     113.9 ,110.75,106.57,101.34,95.477,89.747,84.768,80.697,77.356,74.481,71.869,69.406,67.046,64.787,62.645,60.638,58.78 ,57.075,55.516,54.094,52.793,51.6  ,50.5  ,49.48 ,48.529,47.639,46.802,46.011,45.263,44.551,43.874,43.226,42.605,42.007,41.428,40.865,40.313,39.766,39.218,38.664,38.095,37.504,36.885,36.23 ,35.535,34.801,34.032,33.236,32.429,31.626,
                                                                     112.16,108.95,104.72,99.449,93.581,87.885,82.961,78.948,75.657,72.826,70.253,67.824,65.498,63.27 ,61.16 ,59.185,57.361,55.69 ,54.169,52.786,51.527,50.377,49.321,48.347,47.442,46.598,45.807,45.062,44.359,43.691,43.055,42.448,41.865,41.303,40.758,40.224,39.698,39.173,38.643,38.1  ,37.538,36.949,36.326,35.663,34.961,34.219,33.446,32.654,31.86 ,31.082,
                                                                     110.55,107.3 ,103.02,97.719,91.854,86.199,81.337,77.388,74.154,71.372,68.842,66.454,64.165,61.974,59.898,57.958,56.17 ,54.536,53.053,51.71 ,50.493,49.386,48.374,47.444,46.584,45.784,45.037,44.335,43.673,43.046,42.449,41.878,41.329,40.799,40.281,39.773,39.267,38.758,38.239,37.702,37.14 ,36.547,35.915,35.242,34.527,33.778,33.003,32.217,31.44 ,30.69 ,
                                                                     109.06,105.78,101.47,96.142,90.29 ,84.686,79.895,76.018,72.848,70.122,67.643,65.3  ,63.054,60.904,58.868,56.967,55.217,53.623,52.181,50.88 ,49.705,48.641,47.673,46.787,45.971,45.214,44.509,43.848,43.225,42.636,42.075,41.537,41.02 ,40.517,40.024,39.535,39.046,38.548,38.034,37.498,36.932,36.329,35.684,34.998,34.272,33.515,32.741,31.966,31.21 ,30.491,
                                                                     107.61,104.3 ,99.965,94.631,88.806,83.264,78.554,74.758,71.661,68.999,66.577,64.289,62.092,59.99 ,57.999,56.142,54.436,52.886,51.488,50.231,49.102,48.084,47.161,46.32 ,45.548,44.836,44.173,43.554,42.97 ,42.418,41.892,41.387,40.899,40.422,39.952,39.481,39.005,38.515,38.005,37.466,36.892,36.278,35.62 ,34.92 ,34.185,33.425,32.658,31.9  ,31.172,30.489,
                                                                     106.08,102.76,98.407,93.074,87.286,81.817,77.197,73.489,70.473,67.882,65.526,63.297,61.157,59.108,57.168,55.36 ,53.702,52.199,50.849,49.64 ,48.559,47.589,46.714,45.92 ,45.195,44.528,43.909,43.332,42.789,42.275,41.785,41.313,40.854,40.404,39.955,39.502,39.038,38.556,38.047,37.505,36.923,36.298,35.629,34.92 ,34.181,33.427,32.674,31.943,31.25 ,30.609,
                                                                     104.43,101.1 ,96.736,91.413,85.67 ,80.281,75.757,72.143,69.211,66.696,64.409,62.244,60.164,58.171,56.284,54.527,52.919,51.465,50.163,49.003,47.97 ,47.048,46.222,45.475,44.797,44.174,43.6  ,43.064,42.561,42.085,41.63 ,41.19 ,40.76 ,40.334,39.906,39.469,39.016,38.539,38.03 ,37.484,36.894,36.259,35.582,34.87 ,34.135,33.394,32.665,31.968,31.318,30.725,
                                                                     102.67,99.338,94.973,89.666,83.973,78.668,74.243,70.725,67.88 ,65.442,63.224,61.124,59.104,57.167,55.333,53.627,52.067,50.661,49.406,48.292,47.306,46.43 ,45.649,44.947,44.312,43.732,43.198,42.702,42.236,41.794,41.371,40.959,40.554,40.15 ,39.739,39.314,38.867,38.392,37.88 ,37.327,36.728,36.084,35.401,34.688,33.963,33.242,32.544,31.887,31.283,30.741,
                                                                     100.86,97.532,93.168,87.88 ,82.242,77.025,72.701,69.28 ,66.521,64.16 ,62.013,59.976,58.016,56.133,54.35 ,52.693,51.179,49.817,48.606,47.536,46.593,45.759,45.02 ,44.359,43.764,43.222,42.725,42.263,41.83 ,41.418,41.022,40.634,40.25 ,39.861,39.461,39.042,38.597,38.118,37.599,37.034,36.423,35.77 ,35.082,34.374,33.662,32.966,32.304,31.69 ,31.134,30.641,
                                                                     99.051,95.728,91.37 ,86.104,80.525,75.398,71.176,67.852,65.179,62.896,60.817,58.843,56.94 ,55.111,53.377,51.765,50.295,48.975,47.805,46.775,45.871,45.076,44.374,43.75 ,43.191,42.683,42.219,41.787,41.382,40.995,40.621,40.253,39.883,39.505,39.111,38.693,38.243,37.755,37.223,36.644,36.02 ,35.358,34.668,33.967,33.274,32.607,31.982,31.412,30.904,30.459,
                                                                     97.277,93.962,89.611,84.372,78.854,73.82 ,69.701,66.474,63.887,61.678,59.667,57.755,55.908,54.13 ,52.444,50.875,49.445,48.165,47.032,46.038,45.17 ,44.411,43.743,43.153,42.625,42.148,41.711,41.306,40.925,40.559,40.203,39.848,39.488,39.115,38.722,38.299,37.84 ,37.339,36.791,36.196,35.559,34.889,34.2  ,33.511,32.841,32.207,31.622,31.097,30.635,30.235,
                                                                     95.56 ,92.255,87.914,82.705,77.251,72.311,68.294,65.163,62.661,60.527,58.582,56.73 ,54.937,53.208,51.567,50.039,48.648,47.403,46.305,45.345,44.509,43.781,43.145,42.584,42.085,41.635,41.223,40.841,40.479,40.13 ,39.787,39.442,39.088,38.715,38.317,37.885,37.412,36.893,36.327,35.716,35.067,34.393,33.711,33.039,32.396,31.798,31.256,30.775,30.357,29.999,
                                                                     93.913,90.62 ,86.292,81.114,75.727,70.881,66.966,63.93 ,61.511,59.45 ,57.57 ,55.776,54.036,52.354,50.754,49.265,47.909,46.698,45.631,44.702,43.896,43.197,42.588,42.054,41.58 ,41.153,40.763,40.4  ,40.054,39.719,39.385,39.045,38.691,38.315,37.907,37.461,36.971,36.434,35.849,35.223,34.565,33.892,33.221,32.572,31.961,31.402,30.902,30.465,30.089,29.769,
                                                                     92.346,89.066,84.75 ,79.606,74.286,69.535,65.722,62.779,60.441,58.451,56.633,54.895,53.205,51.568,50.009,48.556,47.233,46.051,45.014,44.112,43.333,42.66 ,42.076,41.565,41.113,40.707,40.335,39.987,39.654,39.328,39    ,38.661,38.304,37.919,37.498,37.035,36.525,35.967,35.364,34.725,34.063,33.396,32.742,32.121,31.545,31.027,30.57 ,30.175,29.838,29.555,
                                                                     90.863,87.595,83.294,78.183,72.931,68.274,64.561,61.71 ,59.451,57.529,55.772,54.088,52.446,50.851,49.329,47.91 ,46.617,45.463,44.452,43.575,42.82 ,42.17 ,41.608,41.118,40.685,40.295,39.938,39.603,39.279,38.959,38.632,38.29 ,37.925,37.527,37.09 ,36.607,36.075,35.497,34.878,34.23 ,33.569,32.914,32.283,31.693,31.157,30.68 ,30.265,29.91 ,29.61 ,29.36 ,
                                                                     89.466,86.21 ,81.922,76.845,71.662,67.099,63.484,60.721,58.539,56.683,54.984,53.351,51.754,50.2  ,48.713,47.324,46.059,44.93 ,43.942,43.088,42.353,41.723,41.181,40.708,40.291,39.916,39.571,39.245,38.927,38.608,38.279,37.931,37.553,37.139,36.682,36.177,35.623,35.026,34.394,33.742,33.088,32.451,31.848,31.294,30.797,30.362,29.988,29.671,29.405,29.184,
                                                                     88.154,84.91 ,80.635,75.592,70.477,66.006,62.487,59.811,57.702,55.91 ,54.266,52.682,51.128,49.611,48.157,46.796,45.556,44.45 ,43.482,42.647,41.931,41.319,40.792,40.335,39.931,39.567,39.231,38.91 ,38.594,38.274,37.939,37.579,37.186,36.753,36.273,35.744,35.17 ,34.556,33.915,33.265,32.624,32.01 ,31.44 ,30.924,30.468,30.074,29.738,29.456,29.222,29.027,
                                                                     86.928,83.694,79.43 ,74.421,69.373,64.993,61.567,58.974,56.937,55.205,53.613,52.076,50.562,49.08 ,47.656,46.322,45.103,44.017,43.068,42.25 ,41.55 ,40.952,40.439,39.994,39.601,39.245,38.914,38.595,38.278,37.952,37.607,37.233,36.821,36.365,35.861,35.31 ,34.715,34.088,33.445,32.803,32.181,31.596,31.061,30.585,30.17 ,29.815,29.515,29.266,29.059,28.889,
                                                                     85.784,82.559,78.304,73.329,68.348,64.056,60.721,58.208,56.239,54.564,53.023,51.529,50.053,48.604,47.207,45.896,44.698,43.63 ,42.696,41.892,41.205,40.619,40.117,39.682,39.296,38.946,38.617,38.297,37.975,37.64 ,37.281,36.888,36.455,35.974,35.446,34.872,34.261,33.627,32.986,32.358,31.761,31.209,30.712,30.276,29.901,29.583,29.317,29.097,28.915,28.766,
                                                                     84.721,81.501,77.255,72.312,67.396,63.191,59.944,57.508,55.604,53.984,52.49 ,51.037,49.596,48.177,46.806,45.517,44.336,43.282,42.362,41.57 ,40.893,40.318,39.824,39.395,39.014,38.666,38.336,38.012,37.682,37.333,36.957,36.543,36.085,35.579,35.027,34.433,33.81 ,33.174,32.542,31.934,31.366,30.85 ,30.393,29.996,29.659,29.376,29.141,28.947,28.788,28.659,
                                                                     83.734,80.517,76.279,71.366,66.514,62.394,59.232,56.87 ,55.027,53.46 ,52.01 ,50.595,49.188,47.797,46.449,45.178,44.013,42.972,42.063,41.28 ,40.612,40.043,39.556,39.131,38.752,38.402,38.069,37.737,37.394,37.029,36.632,36.194,35.71 ,35.179,34.604,33.994,33.364,32.732,32.116,31.533,30.999,30.52 ,30.103,29.745,29.443,29.192,28.985,28.816,28.677,28.564,
                                                                     82.82 ,79.604,75.37 ,70.487,65.698,61.661,58.58 ,56.289,54.505,52.987,51.579,50.201,48.824,47.459,46.132,44.878,43.726,42.696,41.795,41.019,40.357,39.793,39.309,38.885,38.505,38.152,37.811,37.467,37.109,36.724,36.303,35.839,35.328,34.772,34.178,33.557,32.926,32.304,31.71 ,31.157,30.659,30.219,29.841,29.52 ,29.252,29.03 ,28.848,28.7  ,28.579,28.481,
                                                                     81.975,78.756,74.525,69.671,64.944,60.986,57.985,55.762,54.034,52.562,51.194,49.849,48.501,47.16 ,45.852,44.612,43.472,42.45 ,41.555,40.784,40.126,39.564,39.08 ,38.656,38.271,37.911,37.559,37.201,36.823,36.416,35.969,35.477,34.94 ,34.361,33.751,33.125,32.5  ,31.895,31.326,30.808,30.347,29.947,29.606,29.319,29.082,28.887,28.728,28.599,28.493,28.407,
                                                                     81.194,77.97 ,73.74 ,68.913,64.246,60.367,57.441,55.283,53.608,52.18 ,50.849,49.536,48.215,46.895,45.604,44.378,43.246,42.231,41.34 ,40.572,39.916,39.354,38.868,38.439,38.047,37.677,37.31 ,36.934,36.534,36.101,35.626,35.106,34.544,33.947,33.327,32.701,32.087,31.505,30.968,30.486,30.063,29.701,29.395,29.141,28.932,28.761,28.622,28.51 ,28.418,28.343,
                                                                     80.472,77.241,73.011,68.21 ,63.602,59.798,56.946,54.85 ,53.225,51.839,50.542,49.259,47.963,46.662,45.387,44.171,43.047,42.036,41.148,40.38 ,39.723,39.158,38.668,38.231,37.83 ,37.446,37.062,36.664,36.239,35.777,35.274,34.728,34.143,33.531,32.907,32.288,31.693,31.138,30.635,30.19 ,29.806,29.481,29.209,28.984,28.8  ,28.651,28.53 ,28.432,28.352,28.287,
                                                                     79.807,76.566,72.332,67.556,63.006,59.276,56.494,54.458,52.881,51.534,50.27 ,49.014,47.74 ,46.458,45.196,43.989,42.871,41.863,40.975,40.206,39.545,38.975,38.477,38.031,37.616,37.215,36.81 ,36.388,35.935,35.445,34.913,34.341,33.738,33.117,32.495,31.89 ,31.319,30.795,30.329,29.922,29.575,29.284,29.043,28.845,28.684,28.554,28.449,28.363,28.294,28.238,
                                                                     79.192,75.938,71.701,66.949,62.455,58.796,56.083,54.103,52.571,51.261,50.028,48.798,47.545,46.28 ,45.029,43.83 ,42.715,41.708,40.818,40.046,39.379,38.802,38.293,37.834,37.403,36.982,36.553,36.103,35.621,35.101,34.542,33.948,33.332,32.708,32.094,31.509,30.967,30.478,30.049,29.68 ,29.368,29.11 ,28.897,28.723,28.583,28.469,28.378,28.304,28.244,28.195,
                                                                     78.625,75.356,71.113,66.385,61.946,58.356,55.709,53.783,52.294,51.019,49.815,48.608,47.375,46.124,44.884,43.691,42.578,41.57 ,40.676,39.898,39.223,38.635,38.114,37.638,37.188,36.743,36.288,35.808,35.296,34.747,34.163,33.552,32.928,32.308,31.709,31.149,30.639,30.186,29.795,29.462,29.184,28.955,28.768,28.616,28.494,28.395,28.316,28.251,28.199,28.157,
                                                                     78.1  ,74.814,70.565,65.859,61.474,57.953,55.368,53.494,52.046,50.803,49.626,48.442,47.227,45.989,44.758,43.569,42.457,41.445,40.546,39.759,39.074,38.473,37.936,37.441,36.968,36.497,36.012,35.502,34.959,34.382,33.777,33.154,32.529,31.919,31.342,30.811,30.336,29.921,29.566,29.268,29.021,28.819,28.655,28.523,28.416,28.33 ,28.261,28.206,28.161,28.124,
                                                                     77.616,74.31 ,70.052,65.369,61.037,57.581,55.058,53.233,51.824,50.612,49.46 ,48.297,47.098,45.873,44.649,43.463,42.349,41.332,40.426,39.628,38.93 ,38.313,37.757,37.24 ,36.74 ,36.24 ,35.725,35.183,34.611,34.009,33.387,32.758,32.138,31.546,30.994,30.496,30.058,29.68 ,29.361,29.095,28.877,28.7  ,28.556,28.441,28.348,28.274,28.214,28.166,28.127,28.095,
                                                                     77.166,73.84 ,69.573,64.911,60.631,57.24 ,54.775,52.998,51.626,50.443,49.315,48.17 ,46.987,45.772,44.554,43.369,42.252,41.229,40.312,39.502,38.788,38.152,37.574,37.031,36.503,35.972,35.424,34.851,34.251,33.63 ,32.997,32.368,31.761,31.19 ,30.67 ,30.207,29.805,29.463,29.178,28.942,28.75 ,28.595,28.47 ,28.369,28.289,28.225,28.173,28.131,28.097,28.069,
                                                                     76.75 ,73.4  ,69.123,64.483,60.255,56.926,54.518,52.786,51.448,50.293,49.187,48.06 ,46.89 ,45.685,44.471,43.286,42.164,41.133,40.204,39.378,38.646,37.988,37.385,36.814,36.255,35.691,35.11 ,34.508,33.884,33.246,32.609,31.988,31.399,30.856,30.368,29.942,29.577,29.27 ,29.016,28.808,28.639,28.503,28.394,28.307,28.237,28.182,28.137,28.1  ,28.071,28.047,
                                                                     76.362,72.988,68.7  ,64.081,59.904,56.636,54.283,52.594,51.29 ,50.161,49.075,47.965,46.807,45.609,44.399,43.212,42.084,41.042,40.098,39.255,38.501,37.819,37.187,36.585,35.993,35.396,34.783,34.153,33.51 ,32.863,32.228,31.62 ,31.055,30.543,30.091,29.701,29.372,29.097,28.872,28.689,28.542,28.423,28.329,28.253,28.192,28.144,28.105,28.074,28.048,28.027,
                                                                     76.001,72.601,68.302,63.704,59.577,56.369,54.069,52.421,51.149,50.043,48.977,47.882,46.736,45.544,44.335,43.145,42.009,40.954,39.993,39.129,38.351,37.641,36.979,36.344,35.717,35.086,34.443,33.79 ,33.132,32.483,31.857,31.269,30.732,30.254,29.838,29.485,29.189,28.945,28.746,28.585,28.457,28.353,28.271,28.206,28.153,28.111,28.078,28.05 ,28.028,28.01 ,
                                                                     75.664,72.236,67.926,63.349,59.272,56.122,53.873,52.265,51.022,49.94 ,48.892,47.811,46.673,45.487,44.278,43.083,41.937,40.866,39.887,38.999,38.194,37.453,36.758,36.088,35.426,34.763,34.093,33.421,32.755,32.11 ,31.499,30.937,30.431,29.989,29.61 ,29.29 ,29.026,28.81 ,28.635,28.495,28.382,28.293,28.221,28.164,28.119,28.083,28.054,28.03 ,28.011,27.995,
                                                                     75.348,71.891,67.569,63.014,58.986,55.893,53.694,52.124,50.909,49.848,48.816,47.748,46.619,45.436,44.226,43.024,41.866,40.778,39.777,38.863,38.027,37.254,36.523,35.817,35.121,34.427,33.734,33.049,32.382,31.748,31.158,30.625,30.154,29.748,29.404,29.117,28.882,28.692,28.538,28.416,28.318,28.24 ,28.178,28.129,28.089,28.058,28.033,28.012,27.996,27.982,
                                                                     75.051,71.564,67.231,62.697,58.718,55.681,53.529,51.995,50.807,49.767,48.75 ,47.693,46.571,45.391,44.177,42.967,41.794,40.687,39.66 ,38.717,37.849,37.04 ,36.273,35.53 ,34.801,34.08 ,33.37 ,32.679,32.018,31.4  ,30.837,30.336,29.9  ,29.53 ,29.219,28.964,28.756,28.588,28.454,28.347,28.261,28.194,28.14 ,28.098,28.063,28.036,28.014,27.997,27.982,27.97 ,
                                                                     74.771,71.253,66.909,62.397,58.466,55.483,53.378,51.879,50.716,49.694,48.691,47.644,46.528,45.348,44.13 ,42.909,41.719,40.59 ,39.536,38.561,37.658,36.812,36.006,35.228,34.469,33.726,33.004,32.313,31.665,31.069,30.536,30.07 ,29.67 ,29.334,29.056,28.828,28.645,28.497,28.38 ,28.287,28.213,28.154,28.107,28.07 ,28.041,28.017,27.998,27.983,27.97 ,27.96};
        }
        case 10:
        {
            static const std::vector<double> arrayCurvatureMax =    {144.17,138.45,133.04,128   ,123.37,119.14,115.27,111.73,108.48,105.48,102.68,100.07,97.602,95.269,93.049,90.926,88.885,86.915,85.005,83.146,81.329,79.546,77.79 ,76.057,74.342,72.641,70.955,69.284,67.632,66.006,64.413,62.862,61.363,59.922,58.543,57.231,55.985,54.805,53.687,52.629,51.627,50.678,49.778,48.924,48.113,47.344,46.612,45.918,45.258,44.631,
                                                                     142.54,136.78,131.33,126.25,121.58,117.3 ,113.39,109.81,106.52,103.48,100.65,97.999,95.507,93.151,90.913,88.779,86.734,84.768,82.869,81.028,79.235,77.482,75.761,74.066,72.392,70.735,69.095,67.47 ,65.867,64.289,62.745,61.244,59.793,58.401,57.072,55.808,54.611,53.478,52.408,51.397,50.441,49.538,48.682,47.872,47.105,46.377,45.687,45.032,44.411,43.822,
                                                                     140.87,135.07,129.58,124.46,119.75,115.44,111.49,107.88,104.55,101.47,98.607,95.929,93.41 ,91.032,88.777,86.632,84.583,82.62 ,80.732,78.908,77.14 ,75.418,73.734,72.08 ,70.452,68.843,67.252,65.68 ,64.128,62.604,61.113,59.665,58.268,56.928,55.651,54.439,53.292,52.209,51.187,50.224,49.315,48.457,47.646,46.88 ,46.155,45.468,44.818,44.203,43.619,43.066,
                                                                     139.15,133.31,127.79,122.64,117.89,113.55,109.57,105.92,102.56,99.457,96.565,93.861,91.32 ,88.923,86.654,84.499,82.446,80.486,78.608,76.801,75.057,73.366,71.719,70.109,68.527,66.969,65.432,63.915,62.42 ,60.953,59.519,58.128,56.787,55.503,54.28 ,53.121,52.027,50.995,50.023,49.108,48.246,47.434,46.667,45.944,45.261,44.615,44.004,43.426,42.879,42.361,
                                                                     137.39,131.51,125.95,120.77,116   ,111.63,107.62,103.95,100.57,97.441,94.528,91.804,89.246,86.834,84.553,82.391,80.337,78.38 ,76.511,74.72 ,72.999,71.339,69.728,68.16 ,66.626,65.121,63.639,62.179,60.744,59.336,57.963,56.632,55.35 ,54.124,52.958,51.855,50.814,49.835,48.913,48.047,47.232,46.466,45.743,45.062,44.42 ,43.813,43.241,42.699,42.187,41.703,
                                                                     135.58,129.67,124.08,118.88,114.09,109.69,105.67,101.98,98.576,95.43 ,92.502,89.764,87.194,84.773,82.485,80.319,78.265,76.313,74.454,72.679,70.98 ,69.348,67.772,66.245,64.758,63.304,61.877,60.476,59.101,57.756,56.445,55.176,53.956,52.79 ,51.683,50.637,49.652,48.725,47.855,47.038,46.271,45.55 ,44.871,44.232,43.63 ,43.062,42.526,42.02 ,41.542,41.091,
                                                                     133.72,127.79,122.18,116.96,112.15,107.74,103.7 ,99.994,96.584,93.427,90.491,87.747,85.172,82.747,80.458,78.293,76.242,74.297,72.449,70.69 ,69.012,67.406,65.863,64.374,62.931,61.526,60.154,58.811,57.496,56.213,54.966,53.76 ,52.603,51.499,50.453,49.465,48.536,47.664,46.846,46.079,45.359,44.683,44.048,43.45 ,42.888,42.357,41.857,41.385,40.94 ,40.52 ,
                                                                     131.83,125.88,120.25,115.02,110.19,105.77,101.72,98.014,94.599,91.44 ,88.502,85.758,83.184,80.762,78.477,76.319,74.276,72.341,70.506,68.764,67.106,65.525,64.012,62.558,61.156,59.797,58.476,57.188,55.932,54.71 ,53.525,52.384,51.29 ,50.249,49.263,48.335,47.464,46.646,45.881,45.164,44.492,43.862,43.27 ,42.713,42.189,41.695,41.23 ,40.792,40.378,39.988,
                                                                     129.91,123.94,118.3 ,113.05,108.22,103.8 ,99.75 ,96.041,92.628,89.472,86.54 ,83.802,81.236,78.823,76.549,74.402,72.372,70.451,68.632,66.908,65.271,63.713,62.228,60.807,59.442,58.124,56.85 ,55.613,54.412,53.249,52.125,51.046,50.015,49.036,48.113,47.244,46.431,45.669,44.957,44.291,43.667,43.082,42.532,42.016,41.531,41.073,40.642,40.236,39.853,39.491,
                                                                     127.96,121.98,116.33,111.08,106.25,101.83,97.783,94.08 ,90.676,87.53 ,84.61 ,81.885,79.333,76.936,74.678,72.547,70.535,68.633,66.833,65.128,63.513,61.98 ,60.521,59.129,57.797,56.517,55.284,54.093,52.943,51.834,50.767,49.747,48.777,47.86 ,46.997,46.189,45.433,44.728,44.069,43.454,42.878,42.338,41.832,41.356,40.908,40.487,40.089,39.715,39.361,39.028,
                                                                     125.99,120   ,114.35,109.1 ,104.28,99.861,95.829,92.138,88.748,85.619,82.716,80.011,77.479,75.103,72.867,70.758,68.768,66.889,65.111,63.43 ,61.838,60.329,58.895,57.531,56.229,54.982,53.785,52.635,51.53 ,50.469,49.455,48.49 ,47.578,46.719,45.916,45.166,44.468,43.818,43.213,42.648,42.121,41.627,41.164,40.728,40.318,39.932,39.568,39.224,38.9  ,38.593,
                                                                     124.01,118.02,112.37,107.13,102.31,97.91 ,93.893,90.221,86.851,83.744,80.865,78.184,75.678,73.328,71.119,69.037,67.074,65.221,63.47 ,61.815,60.248,58.764,57.356,56.018,54.743,53.525,52.361,51.246,50.179,49.161,48.193,47.278,46.418,45.614,44.866,44.173,43.53 ,42.935,42.383,41.87 ,41.391,40.943,40.523,40.128,39.756,39.405,39.074,38.761,38.465,38.185,
                                                                     122.01,116.03,110.39,105.16,100.36,95.979,91.984,88.335,84.991,81.91 ,79.059,76.408,73.933,71.614,69.436,67.386,65.454,63.631,61.909,60.283,58.744,57.287,55.905,54.592,53.343,52.152,51.015,49.93 ,48.896,47.914,46.986,46.114,45.3  ,44.545,43.849,43.208,42.619,42.077,41.577,41.114,40.683,40.281,39.904,39.55 ,39.216,38.901,38.602,38.32 ,38.053,37.799,
                                                                     120.02,114.05,108.42,103.21,98.434,94.075,90.106,86.486,83.171,80.123,77.305,74.687,72.247,69.963,67.82 ,65.805,63.907,62.118,60.429,58.834,57.325,55.896,54.541,53.254,52.029,50.863,49.751,48.692,47.685,46.733,45.838,45.002,44.228,43.516,42.866,42.273,41.733,41.241,40.791,40.377,39.994,39.637,39.304,38.99 ,38.695,38.415,38.15 ,37.899,37.66 ,37.433,
                                                                     118.04,112.08,106.48,101.29,96.533,92.204,88.266,84.678,81.398,78.385,75.604,73.025,70.622,68.377,66.272,64.295,62.435,60.682,59.028,57.467,55.989,54.59 ,53.264,52.003,50.803,49.66 ,48.571,47.534,46.551,45.623,44.754,43.948,43.206,42.531,41.919,41.368,40.873,40.427,40.023,39.656,39.319,39.007,38.717,38.444,38.187,37.944,37.712,37.492,37.282,37.082,
                                                                     116.07,110.13,104.55,99.39 ,94.668,90.371,86.468,82.918,79.676,76.702,73.961,71.422,69.061,66.857,64.793,62.857,61.036,59.322,57.705,56.179,54.735,53.368,52.07 ,50.836,49.661,48.541,47.473,46.457,45.493,44.586,43.738,42.955,42.239,41.592,41.013,40.497,40.04 ,39.634,39.273,38.949,38.656,38.387,38.139,37.907,37.689,37.482,37.284,37.096,36.915,36.742,
                                                                     114.13,108.21,102.66,97.53 ,92.842,88.583,84.719,81.209,78.008,75.077,72.379,69.883,67.565,65.404,63.383,61.489,59.71 ,58.036,56.458,54.968,53.559,52.224,50.956,49.75 ,48.601,47.503,46.456,45.459,44.513,43.622,42.792,42.027,41.331,40.706,40.151,39.664,39.237,38.866,38.541,38.256,38.003,37.775,37.568,37.376,37.196,37.025,36.863,36.706,36.556,36.41 ,
                                                                     112.21,106.32,100.8 ,95.711,91.062,86.843,83.022,79.554,76.398,73.512,70.858,68.408,66.135,64.019,62.043,60.192,58.455,56.822,55.283,53.831,52.458,51.156,49.92 ,48.742,47.618,46.544,45.517,44.538,43.608,42.733,41.917,41.166,40.484,39.875,39.339,38.872,38.47 ,38.126,37.831,37.58 ,37.362,37.171,37.001,36.847,36.705,36.571,36.443,36.319,36.199,36.082,
                                                                     110.33,104.47,98.992,93.939,89.332,85.157,81.38 ,77.959,74.848,72.008,69.402,66.998,64.772,62.701,60.77 ,58.963,57.269,55.678,54.179,52.765,51.428,50.16 ,48.955,47.806,46.708,45.657,44.651,43.69 ,42.776,41.915,41.111,40.371,39.7  ,39.103,38.579,38.126,37.741,37.417,37.147,36.921,36.734,36.575,36.439,36.321,36.214,36.115,36.021,35.93 ,35.842,35.754,
                                                                     108.49,102.67,97.229,92.219,87.657,83.527,79.797,76.423,73.361,70.569,68.01 ,65.654,63.474,61.45 ,59.565,57.802,56.151,54.601,53.143,51.767,50.466,49.232,48.058,46.939,45.867,44.84 ,43.855,42.912,42.013,41.165,40.372,39.641,38.979,38.389,37.873,37.43 ,37.056,36.746,36.492,36.286,36.122,35.99 ,35.883,35.796,35.721,35.656,35.595,35.537,35.48 ,35.422,
                                                                     106.7 ,100.92,95.518,90.553,86.038,81.957,78.276,74.951,71.938,69.194,66.684,64.376,62.243,60.266,58.425,56.707,55.099,53.59 ,52.171,50.833,49.568,48.368,47.226,46.135,45.09 ,44.086,43.122,42.197,41.315,40.479,39.697,38.975,38.32 ,37.735,37.224,36.786,36.418,36.115,35.871,35.679,35.53 ,35.418,35.335,35.274,35.228,35.193,35.164,35.137,35.111,35.083,
                                                                     104.95,99.217,93.863,88.946,84.48 ,80.448,76.817,73.542,70.578,67.884,65.423,63.162,61.077,59.145,57.35 ,55.675,54.11 ,52.642,51.262,49.961,48.731,47.564,46.453,45.391,44.372,43.392,42.449,41.543,40.676,39.854,39.082,38.368,37.719,37.139,36.631,36.195,35.829,35.529,35.289,35.104,34.965,34.866,34.8  ,34.759,34.737,34.728,34.727,34.729,34.733,34.734,
                                                                     103.27,97.574,92.268,87.4  ,82.984,79.003,75.422,72.197,69.284,66.639,64.226,62.013,59.974,58.088,56.337,54.705,53.181,51.753,50.411,49.147,47.951,46.817,45.736,44.702,43.709,42.753,41.831,40.943,40.092,39.283,38.522,37.817,37.173,36.597,36.091,35.655,35.289,34.989,34.749,34.565,34.431,34.339,34.283,34.256,34.251,34.263,34.285,34.314,34.345,34.375,
                                                                     101.64,95.991,90.733,85.916,81.551,77.621,74.091,70.917,68.053,65.457,63.092,60.926,58.933,57.092,55.384,53.794,52.31 ,50.921,49.616,48.387,47.225,46.123,45.071,44.065,43.097,42.163,41.262,40.393,39.558,38.763,38.013,37.316,36.679,36.107,35.602,35.167,34.799,34.495,34.253,34.067,33.931,33.841,33.788,33.769,33.775,33.802,33.843,33.894,33.949,34.004,
                                                                     100.07,94.47 ,89.262,84.495,80.181,76.303,72.824,69.7  ,66.886,64.338,62.021,59.901,57.953,56.155,54.489,52.94 ,51.495,50.143,48.874,47.679,46.55 ,45.477,44.455,43.475,42.532,41.621,40.74 ,39.889,39.07 ,38.288,37.55 ,36.863,36.232,35.664,35.162,34.726,34.355,34.049,33.802,33.61 ,33.47 ,33.376,33.322,33.304,33.315,33.351,33.405,33.471,33.546,33.623,
                                                                     98.56 ,93.011,87.854,83.138,78.876,75.049,71.621,68.547,65.781,63.281,61.009,58.934,57.03 ,55.274,53.649,52.139,50.732,49.416,48.182,47.02 ,45.922,44.879,43.883,42.929,42.01 ,41.121,40.26 ,39.427,38.624,37.856,37.129,36.451,35.828,35.265,34.765,34.329,33.957,33.646,33.393,33.195,33.048,32.947,32.888,32.866,32.877,32.915,32.975,33.052,33.14 ,33.235,
                                                                     97.116,91.616,86.509,81.845,77.635,73.858,70.48 ,67.455,64.737,62.283,60.056,58.025,56.163,54.448,52.862,51.39 ,50.019,48.738,47.537,46.406,45.338,44.323,43.354,42.424,41.527,40.659,39.818,39.002,38.215,37.461,36.746,36.078,35.462,34.905,34.408,33.973,33.599,33.285,33.027,32.822,32.666,32.557,32.489,32.46 ,32.465,32.5  ,32.56 ,32.641,32.737,32.845,
                                                                     95.735,90.285,85.229,80.616,76.456,72.729,69.4  ,66.423,63.751,61.342,59.16 ,57.171,55.349,53.673,52.125,50.689,49.353,48.105,46.936,45.835,44.795,43.806,42.863,41.957,41.082,40.234,39.411,38.612,37.841,37.1  ,36.397,35.739,35.131,34.58 ,34.087,33.654,33.279,32.962,32.699,32.488,32.324,32.204,32.127,32.087,32.083,32.11 ,32.166,32.245,32.343,32.457,
                                                                     94.418,89.017,84.011,79.448,75.338,71.661,68.379,65.449,62.822,60.457,58.316,56.368,54.586,52.948,51.436,50.035,48.732,47.515,46.376,45.304,44.29 ,43.328,42.408,41.524,40.67 ,39.842,39.037,38.254,37.497,36.77 ,36.079,35.431,34.832,34.286,33.798,33.367,32.994,32.674,32.408,32.19 ,32.018,31.889,31.8  ,31.749,31.734,31.75 ,31.797,31.869,31.964,32.077,
                                                                     93.164,87.812,82.855,78.341,74.28 ,70.651,67.416,64.53 ,61.947,59.625,57.525,55.615,53.871,52.269,50.791,49.423,48.152,46.966,45.855,44.809,43.822,42.883,41.986,41.123,40.289,39.48 ,38.691,37.925,37.182,36.468,35.788,35.15 ,34.559,34.021,33.538,33.11 ,32.738,32.418,32.148,31.925,31.746,31.608,31.509,31.447,31.418,31.422,31.456,31.518,31.605,31.713,
                                                                     91.973,86.667,81.759,77.293,73.28 ,69.697,66.508,63.666,61.125,58.843,56.782,54.91 ,53.201,51.634,50.189,48.852,47.611,46.453,45.369,44.35 ,43.386,42.471,41.595,40.753,39.937,39.145,38.373,37.621,36.892,36.19 ,35.522,34.894,34.312,33.781,33.303,32.88 ,32.509,32.189,31.918,31.691,31.506,31.36 ,31.251,31.177,31.136,31.127,31.147,31.196,31.271,31.369,
                                                                     90.841,85.583,80.721,76.303,72.335,68.798,65.652,62.852,60.352,58.109,56.085,54.249,52.574,51.04 ,49.626,48.32 ,47.107,45.976,44.918,43.923,42.982,42.088,41.233,40.409,39.612,38.836,38.079,37.342,36.626,35.936,35.278,34.66 ,34.087,33.563,33.091,32.672,32.304,31.986,31.713,31.483,31.294,31.141,31.024,30.939,30.886,30.863,30.87 ,30.904,30.964,31.05 ,
                                                                     89.769,84.556,79.74 ,75.367,71.444,67.95 ,64.846,62.087,59.626,57.42 ,55.432,53.63 ,51.988,50.484,49.101,47.822,46.636,45.531,44.497,43.525,42.606,41.732,40.896,40.091,39.31 ,38.55 ,37.808,37.084,36.38 ,35.702,35.055,34.446,33.881,33.365,32.899,32.485,32.12 ,31.803,31.531,31.3  ,31.106,30.949,30.824,30.73 ,30.666,30.631,30.623,30.642,30.688,30.758,
                                                                     88.753,83.584,78.813,74.484,70.604,67.152,64.089,61.368,58.944,56.774,54.82 ,53.05 ,51.439,49.965,48.61 ,47.358,46.198,45.117,44.106,43.155,42.257,41.402,40.584,39.796,39.031,38.285,37.557,36.845,36.154,35.486,34.85 ,34.25 ,33.694,33.184,32.725,32.316,31.955,31.64 ,31.369,31.137,30.941,30.779,30.649,30.547,30.474,30.427,30.406,30.411,30.441,30.495,
                                                                     87.792,82.665,77.937,73.651,69.813,66.401,63.376,60.693,58.304,56.168,54.246,52.508,50.926,49.48 ,48.152,46.926,45.789,44.731,43.741,42.811,41.932,41.096,40.294,39.522,38.772,38.04 ,37.324,36.625,35.945,35.288,34.661,34.07 ,33.522,33.02 ,32.567,32.162,31.805,31.494,31.224,30.992,30.795,30.63 ,30.495,30.387,30.305,30.249,30.216,30.208,30.223,30.262,
                                                                     86.883,81.798,77.11 ,72.865,69.067,65.694,62.706,60.058,57.704,55.6  ,53.709,52    ,50.446,49.027,47.724,46.521,45.408,44.371,43.402,42.491,41.63 ,40.81 ,40.025,39.267,38.531,37.813,37.109,36.421,35.752,35.105,34.487,33.905,33.364,32.869,32.422,32.023,31.671,31.362,31.094,30.863,30.665,30.499,30.36 ,30.247,30.159,30.094,30.051,30.031,30.033,30.056,
                                                                     86.024,80.978,76.331,72.124,68.365,65.028,62.076,59.462,57.14 ,55.067,53.205,51.524,49.997,48.603,47.324,46.144,45.052,44.036,43.086,42.193,41.349,40.545,39.775,39.031,38.308,37.602,36.91 ,36.233,35.573,34.936,34.327,33.753,33.219,32.731,32.29 ,31.896,31.548,31.243,30.977,30.748,30.551,30.383,30.242,30.125,30.031,29.959,29.908,29.878,29.867,29.877,
                                                                     85.212,80.205,75.595,71.426,67.703,64.402,61.484,58.902,56.611,54.567,52.733,51.079,49.577,48.207,46.95 ,45.792,44.72 ,43.723,42.791,41.915,41.087,40.298,39.542,38.812,38.101,37.406,36.725,36.058,35.408,34.779,34.178,33.612,33.086,32.604,32.169,31.781,31.437,31.135,30.872,30.645,30.448,30.28 ,30.137,30.018,29.92 ,29.843,29.785,29.745,29.725,29.723,
                                                                     84.446,79.474,74.901,70.768,67.08 ,63.813,60.927,58.376,56.114,54.098,52.291,50.661,49.183,47.835,46.6  ,45.462,44.409,43.43 ,42.515,41.656,40.843,40.068,39.325,38.607,37.908,37.224,36.553,35.895,35.254,34.634,34.041,33.482,32.963,32.488,32.058,31.675,31.336,31.038,30.778,30.552,30.357,30.189,30.045,29.924,29.823,29.742,29.678,29.631,29.602,29.589,
                                                                     83.721,78.784,74.246,70.148,66.493,63.258,60.403,57.882,55.647,53.658,51.876,50.27 ,48.814,47.488,46.273,45.154,44.119,43.157,42.258,41.414,40.615,39.854,39.123,38.417,37.728,37.055,36.393,35.744,35.111,34.499,33.914,33.362,32.85 ,32.38 ,31.957,31.578,31.243,30.949,30.692,30.468,30.275,30.108,29.964,29.842,29.739,29.654,29.585,29.533,29.496,29.475,
                                                                     83.036,78.133,73.628,69.562,65.94 ,62.736,59.91 ,57.416,55.209,53.244,51.486,49.903,48.468,47.162,45.966,44.865,43.847,42.902,42.018,41.188,40.402,39.654,38.935,38.239,37.561,36.897,36.244,35.604,34.979,34.374,33.796,33.251,32.744,32.281,31.863,31.489,31.158,30.868,30.614,30.393,30.201,30.035,29.892,29.769,29.664,29.577,29.505,29.448,29.405,29.377,
                                                                     82.389,77.517,73.044,69.01 ,65.418,62.243,59.445,56.978,54.796,52.855,51.12 ,49.558,48.144,46.857,45.679,44.595,43.593,42.663,41.793,40.976,40.203,39.467,38.759,38.073,37.405,36.75 ,36.105,35.473,34.855,34.258,33.686,33.147,32.647,32.189,31.776,31.407,31.08 ,30.794,30.543,30.324,30.134,29.97 ,29.827,29.704,29.599,29.51 ,29.435,29.375,29.327,29.293,
                                                                     81.777,76.935,72.492,68.488,64.925,61.778,59.007,56.566,54.407,52.489,50.775,49.233,47.839,46.57 ,45.409,44.341,43.355,42.439,41.583,40.778,40.017,39.292,38.594,37.918,37.259,36.612,35.976,35.351,34.74 ,34.149,33.584,33.051,32.556,32.104,31.696,31.331,31.009,30.726,30.478,30.262,30.074,29.911,29.769,29.647,29.541,29.451,29.375,29.311,29.26 ,29.221,
                                                                     81.197,76.384,71.97 ,67.995,64.46 ,61.34 ,58.594,56.177,54.041,52.145,50.451,48.928,47.552,46.3  ,45.155,44.103,43.131,42.228,41.385,40.593,39.843,39.128,38.44 ,37.774,37.123,36.483,35.855,35.237,34.633,34.048,33.489,32.961,32.472,32.025,31.621,31.261,30.943,30.663,30.418,30.205,30.019,29.858,29.717,29.595,29.49 ,29.399,29.321,29.256,29.202,29.159,
                                                                     80.648,75.862,71.476,67.528,64.019,60.925,58.204,55.809,53.695,51.82 ,50.145,48.641,47.282,46.046,44.917,43.879,42.921,42.031,41.2  ,40.419,39.68 ,38.974,38.296,37.638,36.995,36.363,35.741,35.13 ,34.532,33.953,33.399,32.878,32.394,31.951,31.552,31.196,30.881,30.605,30.363,30.153,29.969,29.809,29.67 ,29.549,29.444,29.353,29.275,29.208,29.152,29.106,
                                                                     80.127,75.367,71.008,67.086,63.603,60.533,57.835,55.462,53.369,51.513,49.857,48.37 ,47.027,45.807,44.692,43.668,42.723,41.845,41.026,40.255,39.526,38.83 ,38.16 ,37.51 ,36.875,36.25 ,35.635,35.03 ,34.438,33.864,33.316,32.8  ,32.32 ,31.882,31.488,31.136,30.825,30.552,30.313,30.105,29.924,29.766,29.628,29.508,29.403,29.312,29.234,29.166,29.108,29.06 ,
                                                                     79.632,74.898,70.563,66.666,63.208,60.161,57.485,55.134,53.06 ,51.223,49.584,48.114,46.787,45.582,44.481,43.47 ,42.537,41.67 ,40.862,40.101,39.382,38.694,38.033,37.39 ,36.762,36.144,35.535,34.936,34.349,33.781,33.238,32.726,32.252,31.818,31.428,31.079,30.772,30.502,30.266,30.061,29.882,29.725,29.589,29.47 ,29.367,29.276,29.197,29.129,29.07 ,29.021,
                                                                     79.162,74.451,70.141,66.268,62.833,59.808,57.154,54.822,52.768,50.948,49.326,47.872,46.56 ,45.368,44.281,43.282,42.361,41.506,40.707,39.956,39.245,38.567,37.913,37.278,36.656,36.044,35.441,34.847,34.266,33.703,33.165,32.658,32.188,31.758,31.371,31.027,30.723,30.456,30.223,30.02 ,29.843,29.689,29.554,29.436,29.333,29.243,29.165,29.096,29.037,28.986,
                                                                     78.715,74.027,69.74 ,65.89 ,62.476,59.473,56.839,54.527,52.491,50.688,49.082,47.643,46.345,45.167,44.092,43.105,42.195,41.35 ,40.561,39.819,39.117,38.446,37.8  ,37.171,36.556,35.95 ,35.352,34.764,34.188,33.629,33.096,32.593,32.127,31.702,31.319,30.978,30.677,30.413,30.183,29.982,29.807,29.655,29.522,29.405,29.303,29.214,29.136,29.067,29.008,28.956,
                                                                     78.288,73.622,69.357,65.529,62.137,59.155,56.54 ,54.246,52.227,50.441,48.85 ,47.426,46.141,44.976,43.912,42.937,42.037,41.202,40.423,39.69 ,38.996,38.332,37.693,37.071,36.462,35.861,35.269,34.685,34.114,33.56 ,33.031,32.532,32.07 ,31.649,31.27 ,30.932,30.634,30.373,30.145,29.947,29.774,29.623,29.492,29.377,29.276,29.187,29.11 ,29.041,28.982,28.93};
        }
        case 11:
        {
            static const std::vector<double> arrayCurvatureMax =    {142.71,139.66,136.61,133.56,130.52,127.47,124.41,121.32,118.14,114.81,111.28,107.64,104.14,100.97,98.181,95.698,93.451,91.389,89.479,87.692,86.002,84.38 ,82.796,81.224,79.639,78.025,76.37 ,74.671,72.93 ,71.153,69.35 ,67.533,65.715,63.907,62.121,60.365,58.65 ,56.981,55.364,53.803,52.301,50.861,49.483,48.168,46.916,45.726,44.596,43.526,42.514,41.558,
                                                                     141.56,138.42,135.28,132.14,129.01,125.86,122.7 ,119.48,116.15,112.62,108.92,105.23,101.8 ,98.742,96.035,93.599,91.376,89.328,87.43 ,85.656,83.98 ,82.375,80.811,79.261,77.702,76.116,74.493,72.829,71.126,69.392,67.635,65.867,64.1  ,62.346,60.615,58.917,57.259,55.647,54.088,52.585,51.14 ,49.756,48.433,47.171,45.971,44.831,43.75 ,42.727,41.76 ,40.846,
                                                                     140.38,137.15,133.92,130.69,127.46,124.22,120.95,117.58,114.06,110.34,106.52,102.85,99.522,96.569,93.924,91.522,89.318,87.284,85.399,83.64 ,81.983,80.398,78.856,77.332,75.801,74.246,72.657,71.03 ,69.369,67.678,65.968,64.251,62.536,60.836,59.161,57.519,55.918,54.364,52.862,51.415,50.026,48.696,47.427,46.217,45.067,43.976,42.942,41.964,41.04 ,40.168,
                                                                     139.16,135.84,132.52,129.2 ,125.88,122.53,119.13,115.6 ,111.89,108   ,104.13,100.53,97.315,94.443,91.844,89.467,87.279,85.259,83.39 ,81.649,80.012,78.451,76.936,75.439,73.939,72.417,70.864,69.277,67.658,66.013,64.351,62.684,61.022,59.377,57.757,56.171,54.626,53.129,51.683,50.292,48.957,47.681,46.463,45.304,44.203,43.159,42.171,41.237,40.354,39.523,
                                                                     137.91,134.5 ,131.09,127.68,124.26,120.8 ,117.25,113.54,109.63,105.63,101.78,98.288,95.167,92.356,89.791,87.433,85.26 ,83.257,81.405,79.686,78.073,76.538,75.051,73.586,72.118,70.632,69.117,67.571,65.996,64.397,62.784,61.168,59.559,57.968,56.403,54.873,53.384,51.943,50.552,49.215,47.934,46.709,45.542,44.432,43.379,42.38 ,41.436,40.543,39.701,38.908,
                                                                     136.63,133.13,129.63,126.12,122.6 ,119.01,115.29,111.38,107.31,103.26,99.496,96.111,93.069,90.303,87.763,85.422,83.265,81.279,79.449,77.754,76.168,74.662,73.206,71.774,70.341,68.892,67.418,65.914,64.383,62.832,61.268,59.704,58.147,56.61 ,55.1  ,53.625,52.191,50.804,49.467,48.183,46.954,45.78 ,44.662,43.599,42.591,41.637,40.735,39.883,39.079,38.323,
                                                                     135.31,131.72,128.13,124.53,120.89,117.15,113.25,109.15,104.96,100.95,97.285,93.995,91.014,88.28 ,85.759,83.434,81.295,79.329,77.523,75.855,74.299,72.825,71.403,70.005,68.61 ,67.2  ,65.766,64.306,62.82 ,61.317,59.803,58.29 ,56.786,55.302,53.846,52.425,51.046,49.712,48.428,47.195,46.016,44.891,43.821,42.804,41.84 ,40.929,40.067,39.254,38.488,37.766,
                                                                     133.96,130.29,126.6 ,122.89,119.12,115.22,111.12,106.85,102.63,98.7  ,95.145,91.933,88.995,86.285,83.782,81.473,79.353,77.412,75.632,73.994,72.47 ,71.03 ,69.643,68.283,66.926,65.556,64.165,62.748,61.309,59.853,58.39 ,56.927,55.476,54.045,52.642,51.274,49.947,48.666,47.433,46.251,45.12 ,44.043,43.018,42.046,41.125,40.254,39.431,38.655,37.925,37.237,
                                                                     132.57,128.82,125.04,121.21,117.29,113.2 ,108.91,104.55,100.36,96.53 ,93.071,89.916,87.01 ,84.319,81.832,79.542,77.444,75.528,73.779,72.173,70.683,69.279,67.93 ,66.608,65.291,63.962,62.614,61.242,59.849,58.442,57.028,55.616,54.216,52.837,51.486,50.171,48.895,47.665,46.482,45.348,44.265,43.234,42.253,41.324,40.443,39.611,38.826,38.086,37.39 ,36.735,
                                                                     131.16,127.31,123.43,119.48,115.39,111.1 ,106.66,102.26,98.159,94.436,91.052,87.939,85.056,82.383,79.913,77.643,75.57 ,73.683,71.966,70.395,68.942,67.575,66.265,64.983,63.706,62.42 ,61.114,59.787,58.441,57.082,55.717,54.355,53.006,51.678,50.379,49.114,47.889,46.707,45.573,44.486,43.449,42.462,41.524,40.635,39.795,39    ,38.251,37.546,36.882,36.257,
                                                                     129.71,125.77,121.78,117.68,113.41,108.93,104.39,100.03,96.039,92.409,89.081,85.999,83.135,80.478,78.027,75.78 ,73.734,71.879,70.196,68.662,67.247,65.92 ,64.649,63.408,62.173,60.929,59.667,58.385,57.085,55.773,54.457,53.145,51.846,50.568,49.318,48.103,46.927,45.793,44.705,43.664,42.671,41.727,40.83 ,39.981,39.178,38.419,37.705,37.032,36.399,35.804,
                                                                     128.23,124.19,120.07,115.81,111.34,106.72,102.15,97.885,93.996,90.441,87.154,84.095,81.247,78.607,76.177,73.956,71.94 ,70.118,68.472,66.976,65.601,64.315,63.085,61.885,60.692,59.491,58.273,57.036,55.782,54.516,53.248,51.984,50.734,49.505,48.304,47.136,46.008,44.92 ,43.878,42.881,41.93 ,41.027,40.169,39.358,38.591,37.868,37.186,36.544,35.941,35.375,
                                                                     126.71,122.56,118.3 ,113.86,109.21,104.5 ,99.979,95.82 ,92.023,88.523,85.266,82.225,79.394,76.773,74.366,72.173,70.19 ,68.404,66.796,65.341,64.007,62.762,61.574,60.415,59.264,58.106,56.931,55.738,54.53 ,53.311,52.089,50.873,49.67 ,48.488,47.335,46.214,45.131,44.088,43.089,42.135,41.225,40.361,39.541,38.766,38.034,37.343,36.693,36.081,35.507,34.967,
                                                                     125.15,120.88,116.46,111.83,107.04,102.31,97.886,93.835,90.111,86.652,83.417,80.393,77.579,74.979,72.597,70.435,68.486,66.738,65.171,63.757,62.465,61.262,60.116,58.999,57.89 ,56.774,55.643,54.494,53.33 ,52.156,50.98 ,49.81 ,48.653,47.518,46.41 ,45.334,44.295,43.296,42.339,41.425,40.554,39.728,38.945,38.204,37.505,36.846,36.226,35.643,35.095,34.581,
                                                                     123.55,119.14,114.54,109.73,104.86,100.19,95.879,91.921,88.253,84.823,81.607,78.599,75.803,73.227,70.874,68.744,66.832,65.124,63.598,62.226,60.977,59.816,58.713,57.638,56.571,55.497,54.407,53.301,52.181,51.051,49.92 ,48.795,47.683,46.592,45.528,44.496,43.5  ,42.542,41.625,40.75 ,39.917,39.127,38.378,37.671,37.003,36.374,35.783,35.227,34.705,34.215,
                                                                     121.89,117.32,112.54,107.6 ,102.73,98.158,93.954,90.07 ,86.443,83.036,79.836,76.846,74.071,71.52 ,69.198,67.103,65.23 ,63.563,62.079,60.751,59.545,58.427,57.365,56.332,55.306,54.273,53.225,52.161,51.083,49.996,48.908,47.826,46.758,45.71 ,44.689,43.699,42.743,41.825,40.947,40.109,39.312,38.556,37.841,37.165,36.528,35.927,35.363,34.833,34.335,33.869,
                                                                     120.16,115.43,110.48,105.47,100.67,96.211,92.102,88.274,84.678,81.29 ,78.108,75.136,72.385,69.862,67.572,65.515,63.681,62.056,60.616,59.332,58.169,57.093,56.073,55.081,54.095,53.102,52.095,51.072,50.035,48.99 ,47.944,46.904,45.878,44.871,43.891,42.941,42.024,41.144,40.303,39.501,38.738,38.015,37.331,36.686,36.077,35.504,34.965,34.46 ,33.985,33.541,
                                                                     118.36,113.45,108.39,103.39,98.689,94.347,90.314,86.529,82.957,79.587,76.422,73.472,70.746,68.254,66    ,63.981,62.189,60.607,59.211,57.97 ,56.85 ,55.816,54.837,53.885,52.939,51.986,51.018,50.034,49.037,48.033,47.027,46.027,45.041,44.074,43.133,42.221,41.342,40.499,39.692,38.924,38.194,37.503,36.849,36.232,35.65 ,35.103,34.589,34.107,33.655,33.231,
                                                                     116.48,111.42,106.3 ,101.38,96.8  ,92.556,88.582,84.83 ,81.28 ,77.928,74.783,71.856,69.159,66.7  ,64.482,62.503,60.754,59.216,57.864,56.667,55.589,54.597,53.658,52.745,51.837,50.922,49.992,49.047,48.088,47.122,46.155,45.194,44.246,43.318,42.414,41.538,40.695,39.887,39.114,38.378,37.679,37.018,36.392,35.802,35.247,34.724,34.234,33.773,33.342,32.938,
                                                                     114.52,109.36,104.26,99.453,94.994,90.83 ,86.902,83.178,79.647,76.315,73.192,70.292,67.625,65.2  ,63.021,61.084,59.378,57.884,56.576,55.422,54.387,53.435,52.535,51.66 ,50.79 ,49.911,49.018,48.109,47.188,46.258,45.328,44.404,43.493,42.601,41.732,40.892,40.083,39.307,38.566,37.861,37.192,36.559,35.96 ,35.396,34.865,34.366,33.897,33.458,33.046,32.661,
                                                                     112.5 ,107.3 ,102.29,97.618,93.262,89.161,85.27 ,81.57 ,78.06 ,74.75 ,71.652,68.78 ,66.146,63.759,61.619,59.724,58.062,56.613,55.349,54.238,53.244,52.332,51.469,50.631,49.796,48.953,48.095,47.221,46.334,45.44 ,44.545,43.656,42.78 ,41.922,41.088,40.281,39.504,38.759,38.049,37.373,36.732,36.125,35.553,35.013,34.505,34.028,33.58 ,33.161,32.768,32.4  ,
                                                                     110.46,105.3 ,100.41,95.864,91.593,87.544,83.685,80.008,76.521,73.235,70.164,67.324,64.725,62.376,60.278,58.426,56.808,55.403,54.182,53.113,52.159,51.285,50.46 ,49.656,48.856,48.046,47.221,46.38 ,45.527,44.666,43.805,42.95 ,42.106,41.281,40.479,39.703,38.957,38.242,37.56 ,36.912,36.298,35.716,35.168,34.651,34.165,33.709,33.281,32.88 ,32.505,32.154,
                                                                     108.43,103.36,98.623,94.183,89.982,85.975,82.145,78.494,75.031,71.772,68.732,65.925,63.363,61.054,58.998,57.19 ,55.616,54.255,53.077,52.049,51.134,50.297,49.506,48.736,47.968,47.19 ,46.396,45.587,44.765,43.936,43.107,42.283,41.471,40.677,39.905,39.158,38.441,37.755,37.1  ,36.478,35.888,35.33 ,34.805,34.31 ,33.845,33.408,32.999,32.615,32.257,31.921,
                                                                     106.44,101.52,96.913,92.565,88.421,84.452,80.652,77.027,73.592,70.363,67.356,64.585,62.062,59.794,57.781,56.017,54.487,53.169,52.033,51.045,50.167,49.366,48.608,47.87 ,47.132,46.384,45.62 ,44.841,44.048,43.249,42.449,41.655,40.872,40.107,39.364,38.646,37.956,37.295,36.666,36.068,35.502,34.967,34.463,33.989,33.543,33.125,32.733,32.366,32.023,31.702,
                                                                     104.54,99.763,95.273,91.002,86.908,82.975,79.206,75.611,72.207,69.01 ,66.038,63.305,60.823,58.597,56.627,54.907,53.422,52.146,51.051,50.101,49.26 ,48.491,47.765,47.057,46.348,45.628,44.891,44.14 ,43.375,42.604,41.831,41.065,40.31 ,39.572,38.855,38.164,37.499,36.863,36.258,35.683,35.139,34.626,34.142,33.686,33.259,32.858,32.483,32.131,31.803,31.496,
                                                                     102.71,98.084,93.694,89.489,85.441,81.544,77.808,74.245,70.875,67.715,64.78 ,62.087,59.647,57.464,55.538,53.862,52.42 ,51.186,50.13 ,49.217,48.41 ,47.674,46.977,46.296,45.614,44.92 ,44.209,43.483,42.744,41.999,41.252,40.512,39.783,39.07 ,38.379,37.711,37.07 ,36.458,35.875,35.322,34.799,34.305,33.84 ,33.402,32.992,32.608,32.248,31.911,31.596,31.303,
                                                                     100.98,96.473,92.169,88.023,84.018,80.16 ,76.459,72.933,69.6  ,66.477,63.583,60.932,58.534,56.395,54.514,52.882,51.482,50.289,49.271,48.393,47.619,46.912,46.243,45.588,44.93 ,44.26 ,43.573,42.87 ,42.155,41.434,40.711,39.995,39.289,38.6  ,37.932,37.288,36.669,36.078,35.516,34.983,34.479,34.004,33.556,33.136,32.742,32.372,32.027,31.704,31.402,31.121,
                                                                     99.315,94.922,90.692,86.6  ,82.64 ,78.822,75.161,71.674,68.381,65.3  ,62.448,59.84 ,57.486,55.392,53.555,51.966,50.607,49.454,48.473,47.629,46.884,46.205,45.562,44.93 ,44.295,43.647,42.981,42.3  ,41.607,40.907,40.207,39.512,38.829,38.162,37.515,36.892,36.294,35.723,35.181,34.666,34.18 ,33.722,33.291,32.886,32.507,32.152,31.82 ,31.509,31.22 ,30.95 ,
                                                                     97.72 ,93.421,89.259,85.22 ,81.307,77.533,73.914,70.47 ,67.22 ,64.183,61.375,58.812,56.503,54.453,52.66 ,51.114,49.797,48.681,47.735,46.923,46.207,45.554,44.933,44.323,43.708,43.08 ,42.434,41.772,41.098,40.418,39.738,39.064,38.401,37.754,37.127,36.523,35.944,35.392,34.867,34.37 ,33.901,33.458,33.043,32.653,32.287,31.945,31.625,31.327,31.048,30.789,
                                                                     96.182,91.966,87.868,83.884,80.02 ,76.293,72.72 ,69.322,66.118,63.126,60.365,57.848,55.585,53.58 ,51.831,50.327,49.049,47.97 ,47.058,46.275,45.586,44.956,44.356,43.766,43.169,42.557,41.929,41.285,40.629,39.966,39.304,38.648,38.003,37.375,36.766,36.18 ,35.618,35.083,34.575,34.094,33.64 ,33.212,32.811,32.434,32.081,31.751,31.443,31.156,30.887,30.637,
                                                                     94.692,90.554,86.519,82.589,78.778,75.102,71.579,68.23 ,65.075,62.132,59.418,56.948,54.731,52.772,51.066,49.604,48.364,47.321,46.439,45.684,45.019,44.411,43.83 ,43.256,42.675,42.079,41.466,40.837,40.196,39.55 ,38.904,38.264,37.636,37.024,36.432,35.862,35.316,34.796,34.303,33.837,33.397,32.983,32.594,32.23 ,31.889,31.57 ,31.272,30.994,30.735,30.494,
                                                                     93.246,89.181,85.209,81.339,77.583,73.961,70.492,67.195,64.091,61.198,58.534,56.112,53.942,52.028,50.366,48.944,47.741,46.731,45.88 ,45.15 ,44.507,43.918,43.353,42.795,42.227,41.644,41.044,40.428,39.801,39.169,38.537,37.911,37.298,36.701,36.123,35.567,35.036,34.531,34.051,33.598,33.17 ,32.769,32.391,32.038,31.708,31.399,31.11 ,30.842,30.591,30.357,
                                                                     91.84 ,87.847,83.94 ,80.131,76.435,72.871,69.458,66.216,63.166,60.325,57.712,55.339,53.217,51.348,49.729,48.346,47.18 ,46.201,45.378,44.672,44.049,43.476,42.926,42.379,41.823,41.251,40.662,40.057,39.441,38.821,38.201,37.588,36.988,36.403,35.838,35.296,34.777,34.284,33.817,33.375,32.959,32.568,32.201,31.858,31.537,31.237,30.956,30.695,30.452,30.225,
                                                                     90.471,86.551,82.712,78.968,75.335,71.832,68.478,65.294,62.3  ,59.513,56.952,54.63 ,52.555,50.732,49.154,47.81 ,46.678,45.73 ,44.932,44.248,43.643,43.084,42.546,42.009,41.462,40.899,40.318,39.722,39.116,38.505,37.896,37.294,36.704,36.131,35.577,35.046,34.538,34.056,33.599,33.168,32.761,32.38 ,32.022,31.687,31.373,31.08 ,30.807,30.552,30.314,30.092,
                                                                     89.139,85.293,81.524,77.849,74.282,70.843,67.552,64.428,61.492,58.761,56.254,53.982,51.956,50.177,48.641,47.334,46.236,45.316,44.542,43.877,43.287,42.741,42.212,41.684,41.143,40.586,40.012,39.423,38.824,38.221,37.619,37.026,36.445,35.881,35.337,34.815,34.316,33.843,33.395,32.972,32.574,32.2  ,31.849,31.52 ,31.213,30.926,30.657,30.406,30.172,29.953,
                                                                     87.844,84.073,80.378,76.774,73.277,69.905,66.679,63.618,60.742,58.069,55.615,53.395,51.417,49.683,48.188,46.918,45.851,44.958,44.206,43.559,42.982,42.445,41.924,41.4  ,40.865,40.312,39.741,39.157,38.562,37.965,37.37 ,36.783,36.209,35.652,35.115,34.6  ,34.109,33.643,33.201,32.785,32.392,32.024,31.678,31.354,31.05 ,30.766,30.5  ,30.25 ,30.017,29.797,
                                                                     86.584,82.892,79.273,75.743,72.319,69.017,65.859,62.863,60.048,57.434,55.036,52.868,50.938,49.249,47.794,46.559,45.522,44.654,43.922,43.29 ,42.724,42.194,41.678,41.158,40.624,40.073,39.504,38.922,38.33 ,37.735,37.144,36.561,35.992,35.439,34.907,34.397,33.911,33.449,33.012,32.599,32.21 ,31.844,31.5  ,31.178,30.875,30.59 ,30.323,30.072,29.836,29.613,
                                                                     85.36 ,81.749,78.209,74.757,71.407,68.179,65.09 ,62.161,59.41 ,56.855,54.514,52.399,50.517,48.872,47.456,46.255,45.247,44.403,43.689,43.069,42.512,41.987,41.473,40.953,40.419,39.867,39.297,38.714,38.122,37.528,36.938,36.356,35.788,35.238,34.707,34.199,33.715,33.254,32.818,32.405,32.016,31.65 ,31.305,30.98 ,30.674,30.386,30.114,29.858,29.616,29.387,
                                                                     84.172,80.644,77.187,73.814,70.542,67.389,64.372,61.511,58.824,56.331,54.047,51.984,50.152,48.55 ,47.172,46.005,45.024,44.201,43.502,42.893,42.342,41.819,41.305,40.783,40.245,39.689,39.116,38.529,37.934,37.337,36.744,36.16 ,35.591,35.038,34.506,33.996,33.51 ,33.047,32.607,32.192,31.799,31.427,31.077,30.746,30.434,30.139,29.86 ,29.596,29.347,29.11 ,
                                                                     83.021,79.578,76.204,72.914,69.722,66.645,63.702,60.91 ,58.29 ,55.859,53.632,51.623,49.838,48.279,46.939,45.803,44.848,44.044,43.359,42.758,42.209,41.686,41.168,40.64 ,40.096,39.533,38.952,38.357,37.755,37.152,36.552,35.963,35.387,34.829,34.291,33.775,33.282,32.812,32.366,31.943,31.543,31.164,30.806,30.467,30.147,29.845,29.559,29.289,29.034,28.793,
                                                                     81.905,78.549,75.262,72.055,68.944,65.946,63.077,60.357,57.804,55.435,53.267,51.31 ,49.573,48.056,46.752,45.646,44.714,43.927,43.252,42.656,42.107,41.58 ,41.054,40.517,39.962,39.387,38.794,38.188,37.575,36.96 ,36.35 ,35.75 ,35.164,34.596,34.048,33.522,33.019,32.539,32.084,31.652,31.243,30.856,30.491,30.147,29.822,29.516,29.229,28.959,28.705,28.467,
                                                                     80.824,77.557,74.356,71.235,68.207,65.288,62.495,59.847,57.361,55.055,52.944,51.04 ,49.35 ,47.874,46.604,45.526,44.614,43.841,43.173,42.578,42.025,41.489,40.952,40.401,39.831,39.24 ,38.631,38.008,37.379,36.748,36.123,35.508,34.909,34.327,33.767,33.23 ,32.716,32.228,31.765,31.326,30.912,30.522,30.156,29.812,29.49 ,29.188,28.907,28.644,28.4  ,28.173,
                                                                     79.775,76.598,73.486,70.45 ,67.505,64.666,61.95 ,59.374,56.955,54.711,52.657,50.805,49.16 ,47.723,46.485,45.431,44.537,43.774,43.109,42.51 ,41.949,41.401,40.847,40.278,39.688,39.076,38.446,37.804,37.155,36.506,35.864,35.234,34.621,34.027,33.456,32.911,32.391,31.899,31.434,30.996,30.584,30.199,29.839,29.503,29.19 ,28.9  ,28.63 ,28.38 ,28.149,27.936,
                                                                     78.756,75.669,72.645,69.696,66.834,64.075,61.434,58.929,56.577,54.394,52.396,50.593,48.991,47.59 ,46.381,45.348,44.467,43.71 ,43.045,42.439,41.865,41.299,40.726,40.135,39.522,38.888,38.237,37.574,36.906,36.241,35.585,34.943,34.32 ,33.72 ,33.146,32.599,32.08 ,31.591,31.131,30.7  ,30.298,29.923,29.574,29.25 ,28.951,28.673,28.417,28.181,27.963,27.763,
                                                                     77.76 ,74.763,71.827,68.963,66.183,63.502,60.936,58.5  ,56.213,54.089,52.144,50.388,48.827,47.459,46.276,45.261,44.391,43.637,42.968,42.352,41.763,41.18 ,40.586,39.974,39.34 ,38.686,38.016,37.337,36.657,35.981,35.317,34.67 ,34.046,33.446,32.874,32.332,31.821,31.34 ,30.889,30.469,30.077,29.714,29.377,29.065,28.777,28.511,28.266,28.041,27.833,27.642,
                                                                     76.78 ,73.87 ,71.02 ,68.239,65.54 ,62.935,60.441,58.073,55.848,53.782,51.888,50.177,48.654,47.318,46.159,45.162,44.301,43.55 ,42.878,42.254,41.654,41.056,40.446,39.818,39.169,38.501,37.82 ,37.132,36.445,35.766,35.1  ,34.454,33.832,33.237,32.671,32.136,31.632,31.16 ,30.719,30.307,29.925,29.571,29.243,28.94 ,28.66 ,28.402,28.165,27.947,27.746,27.562,
                                                                     75.801,72.977,70.21 ,67.51 ,64.889,62.359,59.935,57.633,55.469,53.459,51.616,49.95 ,48.467,47.164,46.033,45.056,44.21 ,43.467,42.797,42.172,41.566,40.96 ,40.342,39.705,39.048,38.373,37.687,36.995,36.305,35.625,34.96 ,34.316,33.697,33.105,32.544,32.014,31.515,31.048,30.613,30.207,29.83 ,29.481,29.158,28.86 ,28.586,28.333,28.1  ,27.885,27.689,27.508,
                                                                     74.811,72.07 ,69.384,66.764,64.219,61.762,59.409,57.174,55.074,53.123,51.335,49.72 ,48.282,47.019,45.921,44.973,44.149,43.422,42.762,42.143,41.539,40.932,40.311,39.67 ,39.007,38.328,37.637,36.941,36.249,35.566,34.9  ,34.255,33.635,33.044,32.483,31.954,31.457,30.991,30.557,30.153,29.778,29.431,29.111,28.815,28.542,28.29 ,28.059,27.847,27.652,27.473,
                                                                     73.799,71.141,68.537,65.997,63.531,61.152,58.874,56.712,54.681,52.797,51.072,49.516,48.132,46.918,45.863,44.951,44.156,43.451,42.807,42.198,41.598,40.992,40.369,39.723,39.055,38.369,37.671,36.969,36.27 ,35.581,34.909,34.259,33.635,33.04 ,32.476,31.943,31.444,30.976,30.54 ,30.135,29.759,29.411,29.089,28.793,28.519,28.268,28.037,27.825,27.63 ,27.452,
                                                                     72.773,70.201,67.684,65.23 ,62.85 ,60.556,58.361,56.28 ,54.328,52.52 ,50.867,49.377,48.055,46.896,45.889,45.016,44.253,43.572,42.945,42.344,41.749,41.142,40.513,39.86 ,39.183,38.486,37.777,37.064,36.354,35.655,34.974,34.315,33.682,33.08 ,32.508,31.97 ,31.464,30.992,30.551,30.142,29.762,29.411,29.087,28.789,28.514,28.261,28.028,27.815,27.62 ,27.441};
        }
        }

        return 0;
    }

    double getCurvatureLowerLimit(double circumference, double aspectRatio, int windowLength)
    {
        static const std::vector<double> arrayCircumferences = {60,65.918,71.837,77.755,83.673,89.592,95.51,101.43,107.35,113.27,119.18,125.1,131.02,136.94,142.86,148.78,154.69,160.61,166.53,172.45,178.37,184.29,190.2,196.12,202.04,207.96,213.88,219.8,225.71,231.63,237.55,243.47,249.39,255.31,261.22,267.14,273.06,278.98,284.9,290.82,296.73,302.65,308.57,314.49,320.41,326.33,332.24,338.16,344.08,350};
        static const std::vector<double> arrayAspectRatios   = {0.4,0.41224,0.42449,0.43673,0.44898,0.46122,0.47347,0.48571,0.49796,0.5102,0.52245,0.53469,0.54694,0.55918,0.57143,0.58367,0.59592,0.60816,0.62041,0.63265,0.6449,0.65714,0.66939,0.68163,0.69388,0.70612,0.71837,0.73061,0.74286,0.7551,0.76735,0.77959,0.79184,0.80408,0.81633,0.82857,0.84082,0.85306,0.86531,0.87755,0.8898,0.90204,0.91429,0.92653,0.93878,0.95102,0.96327,0.97551,0.98776,1};

        int x = findUpperBound(arrayCircumferences, circumference);
        int y = findUpperBound(arrayAspectRatios,   aspectRatio);

        int arrayWidth = arrayCircumferences.size();

        switch(windowLength)
        {

        case 5:
        {
            static const std::vector<double> arrayCurvatureMin   = {-0.35553 ,-1.221  ,-2.8381  ,-5.0028 ,-6.8972 ,-8.0299  ,-8.5603 ,-8.7801  ,-8.8665   ,-8.8998  ,-8.9125  ,-8.9173  ,-8.9191   ,-8.9198   ,-8.9201   ,-8.9202   ,-8.9202    ,-8.9202  ,-8.9201   ,-8.9201 ,-8.9201  ,-8.9201 ,-8.92    ,-8.92   ,-8.92    ,-8.9199 ,-8.92    ,-8.9201  ,-8.9206   ,-8.9215 ,-8.9237  ,-8.9283  ,-8.9377,-8.955 ,-8.9824,-9.0171,-9.0504,-9.0746,-9.089 ,-9.0963,-9.0997,-9.101 ,-9.1013,-9.101 ,-9.1004,-9.0996,-9.0986,-9.0975,-9.0962,-9.0946,
                                                                    -0.070614,-0.5865 ,-1.6947  ,-3.5674 ,-5.7376 ,-7.3802  ,-8.2671 ,-8.6606  ,-8.8198   ,-8.8818  ,-8.9056  ,-8.9147  ,-8.9181   ,-8.9194   ,-8.9199   ,-8.9201   ,-8.9201    ,-8.9201  ,-8.9201   ,-8.9201 ,-8.9201  ,-8.92   ,-8.92    ,-8.9199 ,-8.9199  ,-8.9199 ,-8.9199  ,-8.92    ,-8.9204   ,-8.9212 ,-8.9231  ,-8.9273  ,-8.9358,-8.9517,-8.9777,-9.0118,-9.0458,-9.0715,-9.0871,-9.0952,-9.0989,-9.1003,-9.1006,-9.1003,-9.0996,-9.0987,-9.0977,-9.0964,-9.0949,-9.0932,
                                                                    0.089106 ,-0.20364,-0.89012 ,-2.2684 ,-4.3377 ,-6.4003  ,-7.7665 ,-8.4448  ,-8.7335   ,-8.8484  ,-8.8928  ,-8.9098  ,-8.9163   ,-8.9187   ,-8.9196   ,-8.92     ,-8.9201    ,-8.9201  ,-8.9201   ,-8.9201 ,-8.92    ,-8.92   ,-8.9199  ,-8.9199 ,-8.9198  ,-8.9198 ,-8.9198  ,-8.9199  ,-8.9202   ,-8.9209 ,-8.9226  ,-8.9263  ,-8.934 ,-8.9486,-8.9731,-9.0064,-9.041 ,-9.0681,-9.0849,-9.0938,-9.098 ,-9.0996,-9.0999,-9.0995,-9.0988,-9.0978,-9.0966,-9.0951,-9.0935,-9.0916,
                                                                    0.1766   ,0.015126,-0.38356 ,-1.2804 ,-2.9348 ,-5.1071  ,-6.9697 ,-8.0667  ,-8.576    ,-8.7863  ,-8.8689  ,-8.9006  ,-8.9128   ,-8.9174   ,-8.9191   ,-8.9198   ,-8.92      ,-8.92    ,-8.9201   ,-8.92   ,-8.92    ,-8.9199 ,-8.9199  ,-8.9198 ,-8.9198  ,-8.9197 ,-8.9197  ,-8.9198  ,-8.92     ,-8.9206 ,-8.9221  ,-8.9254  ,-8.9324,-8.9458,-8.9688,-9.0011,-9.036 ,-9.0645,-9.0826,-9.0924,-9.097 ,-9.0987,-9.0991,-9.0987,-9.0979,-9.0967,-9.0954,-9.0938,-9.0919,-9.0898,
                                                                    0.22393  ,0.13626 ,-0.086579,-0.62366,-1.7678 ,-3.6719  ,-5.8337 ,-7.439   ,-8.2948   ,-8.6721  ,-8.8243  ,-8.8835  ,-8.9062   ,-8.9149   ,-8.9181   ,-8.9194   ,-8.9198    ,-8.92    ,-8.92     ,-8.92   ,-8.9199  ,-8.9199 ,-8.9198  ,-8.9198 ,-8.9197  ,-8.9196 ,-8.9196  ,-8.9196  ,-8.9198   ,-8.9203 ,-8.9217  ,-8.9246  ,-8.9308,-8.9431,-8.9645,-8.9958,-9.0309,-9.0606,-9.0801,-9.0908,-9.0958,-9.0978,-9.0982,-9.0977,-9.0968,-9.0955,-9.094 ,-9.0922,-9.0902,-9.0878,
                                                                    0.24936  ,0.20217 ,0.08029  ,-0.22534,-0.93843,-2.355   ,-4.4447 ,-6.4845  ,-7.8126   ,-8.4653  ,-8.7418  ,-8.8516  ,-8.894    ,-8.9102   ,-8.9163   ,-8.9187   ,-8.9195    ,-8.9198  ,-8.9199   ,-8.9199 ,-8.9199  ,-8.9198 ,-8.9198  ,-8.9197 ,-8.9196  ,-8.9195 ,-8.9195  ,-8.9195  ,-8.9196   ,-8.9201 ,-8.9212  ,-8.9238  ,-8.9294,-8.9405,-8.9605,-8.9905,-9.0256,-9.0564,-9.0773,-9.089 ,-9.0946,-9.0967,-9.0972,-9.0967,-9.0956,-9.0942,-9.0925,-9.0905,-9.0882,-9.0855,
                                                                    0.26297  ,0.23769 ,0.17182  ,0.00297 ,-0.41272,-1.3416  ,-3.0327 ,-5.2106  ,-7.0403   ,-8.102   ,-8.5911  ,-8.7922  ,-8.8711   ,-8.9014   ,-8.913    ,-8.9174   ,-8.919     ,-8.9196  ,-8.9198   ,-8.9199 ,-8.9198  ,-8.9198 ,-8.9197  ,-8.9196 ,-8.9195  ,-8.9194 ,-8.9194  ,-8.9193  ,-8.9194   ,-8.9198 ,-8.9207  ,-8.923   ,-8.928 ,-8.9381,-8.9566,-8.9853,-9.0201,-9.052 ,-9.0743,-9.087 ,-9.0932,-9.0956,-9.096 ,-9.0955,-9.0943,-9.0927,-9.0908,-9.0886,-9.086 ,-9.083 ,
                                                                    0.27025  ,0.25673 ,0.22137  ,0.12961 ,-0.10323,-0.66223 ,-1.8428 ,-3.7771  ,-5.9284   ,-7.4959  ,-8.3213  ,-8.683   ,-8.8285   ,-8.8851   ,-8.9067   ,-8.915    ,-8.9181    ,-8.9192  ,-8.9196   ,-8.9198 ,-8.9198  ,-8.9197 ,-8.9196  ,-8.9195 ,-8.9194  ,-8.9193 ,-8.9192  ,-8.9192  ,-8.9192   ,-8.9195 ,-8.9203  ,-8.9223  ,-8.9267,-8.9359,-8.9529,-8.9801,-9.0145,-9.0473,-9.0711,-9.0848,-9.0916,-9.0943,-9.0948,-9.0941,-9.0928,-9.0911,-9.0889,-9.0864,-9.0835,-9.0802,
                                                                    0.27413  ,0.26692 ,0.248    ,0.19859 ,0.071086,-0.24795 ,-0.9884 ,-2.4433  ,-4.5516   ,-6.5669  ,-7.8571  ,-8.4849  ,-8.7497   ,-8.8546   ,-8.895    ,-8.9105   ,-8.9164    ,-8.9186  ,-8.9194   ,-8.9196 ,-8.9197  ,-8.9196 ,-8.9195  ,-8.9194 ,-8.9193  ,-8.9192 ,-8.9191  ,-8.919   ,-8.919    ,-8.9192 ,-8.9198  ,-8.9216  ,-8.9255,-8.9337,-8.9494,-8.9751,-9.0089,-9.0424,-9.0675,-9.0824,-9.0898,-9.0928,-9.0933,-9.0926,-9.0912,-9.0892,-9.0868,-9.084 ,-9.0807,-9.077 ,
                                                                    0.27621  ,0.27236 ,0.26226  ,0.23578 ,0.16683 ,-0.009712,-0.44303,-1.4046  ,-3.132    ,-5.3131  ,-7.1089  ,-8.136   ,-8.6054   ,-8.7979   ,-8.8732   ,-8.9021   ,-8.9131    ,-8.9173  ,-8.9188   ,-8.9194 ,-8.9195  ,-8.9195 ,-8.9194  ,-8.9193 ,-8.9192  ,-8.919  ,-8.9189  ,-8.9188  ,-8.9187   ,-8.9188 ,-8.9194  ,-8.9208  ,-8.9243,-8.9317,-8.9461,-8.9702,-9.0031,-9.0372,-9.0636,-9.0798,-9.0879,-9.0911,-9.0918,-9.091 ,-9.0893,-9.0871,-9.0844,-9.0813,-9.0777,-9.0735,
                                                                    0.27733  ,0.27526 ,0.26988  ,0.25573 ,0.2187  ,0.12268  ,-0.12058,-0.70223 ,-1.9196   ,-3.8829  ,-6.0216  ,-7.5509  ,-8.3468   ,-8.6934   ,-8.8325   ,-8.8865   ,-8.9072    ,-8.915   ,-8.9179   ,-8.919  ,-8.9193  ,-8.9194 ,-8.9193  ,-8.9192 ,-8.919   ,-8.9189 ,-8.9187  ,-8.9185  ,-8.9184   ,-8.9185 ,-8.9189  ,-8.9201  ,-8.9232,-8.9298,-8.9429,-8.9655,-8.9974,-9.0317,-9.0595,-9.0769,-9.0857,-9.0893,-9.09  ,-9.0891,-9.0872,-9.0848,-9.0818,-9.0782,-9.0742,-9.0695,
                                                                    0.27794  ,0.27682 ,0.27395  ,0.2664  ,0.24659 ,0.19486  ,0.061489,-0.27147 ,-1.0401   ,-2.5333  ,-4.6582  ,-6.6474  ,-7.9      ,-8.5037   ,-8.7571   ,-8.8573   ,-8.896     ,-8.9107  ,-8.9162   ,-8.9183 ,-8.919   ,-8.9192 ,-8.9192  ,-8.919  ,-8.9189  ,-8.9187 ,-8.9185  ,-8.9183  ,-8.9181   ,-8.9181 ,-8.9184  ,-8.9194  ,-8.9221,-8.9279,-8.9398,-8.9609,-8.9916,-9.0261,-9.055 ,-9.0736,-9.0833,-9.0873,-9.088 ,-9.087 ,-9.0849,-9.0821,-9.0788,-9.0749,-9.0703,-9.0651,
                                                                    0.27829  ,0.27766 ,0.27612  ,0.2721  ,0.26153 ,0.23381  ,0.16164 ,-0.022927,-0.47452  ,-1.4695  ,-3.2324  ,-5.4146  ,-7.1755   ,-8.1685   ,-8.6191   ,-8.8031   ,-8.8751    ,-8.9027  ,-8.9131   ,-8.917  ,-8.9185  ,-8.9189 ,-8.919   ,-8.9189 ,-8.9187  ,-8.9185 ,-8.9183  ,-8.918   ,-8.9178   ,-8.9177 ,-8.9179  ,-8.9187  ,-8.9209,-8.9261,-8.9369,-8.9563,-8.9858,-9.0201,-9.0502,-9.0701,-9.0807,-9.085 ,-9.0858,-9.0846,-9.0823,-9.0792,-9.0755,-9.0711,-9.066 ,-9.0602,
                                                                    0.27853  ,0.27811 ,0.27729  ,0.27515 ,0.26952 ,0.25471  ,0.21594 ,0.11546  ,-0.13865  ,-0.74367 ,-1.9984  ,-3.9891  ,-6.113    ,-7.6041   ,-8.3711   ,-8.7032   ,-8.8362    ,-8.8877  ,-8.9074   ,-8.9148 ,-8.9175  ,-8.9185 ,-8.9187  ,-8.9187 ,-8.9185  ,-8.9183 ,-8.918   ,-8.9177  ,-8.9174   ,-8.9173 ,-8.9173  ,-8.9179  ,-8.9198,-8.9244,-8.934 ,-8.952 ,-8.98  ,-9.014 ,-9.045 ,-9.0662,-9.0777,-9.0824,-9.0833,-9.082 ,-9.0794,-9.076 ,-9.0717,-9.0668,-9.0612,-9.0547,
                                                                    0.27873  ,0.27837 ,0.27792  ,0.27679 ,0.27379 ,0.26589  ,0.24515 ,0.191    ,0.051503  ,-0.29593 ,-1.0934  ,-2.6249  ,-4.7643   ,-6.7259   ,-7.9412   ,-8.5215   ,-8.7641    ,-8.8598  ,-8.8966   ,-8.9106 ,-8.9158  ,-8.9177 ,-8.9183  ,-8.9184 ,-8.9182  ,-8.918  ,-8.9177  ,-8.9174  ,-8.9171   ,-8.9168 ,-8.9167  ,-8.9172  ,-8.9187,-8.9227,-8.9313,-8.9477,-8.9742,-9.0077,-9.0395,-9.062 ,-9.0744,-9.0796,-9.0805,-9.0791,-9.0762,-9.0723,-9.0676,-9.0621,-9.0558,-9.0485,
                                                                    0.27899  ,0.27853 ,0.27826  ,0.27766 ,0.27608 ,0.27188  ,0.26081 ,0.23179  ,0.15626   ,-0.036671,-0.5072  ,-1.5362  ,-3.3338   ,-5.5148   ,-7.2401   ,-8.1997   ,-8.632     ,-8.808   ,-8.8766   ,-8.9029 ,-8.9128  ,-8.9165 ,-8.9177  ,-8.918  ,-8.9179  ,-8.9177 ,-8.9174  ,-8.917   ,-8.9166   ,-8.9163 ,-8.9161  ,-8.9163  ,-8.9176,-8.9209,-8.9286,-8.9435,-8.9684,-9.0012,-9.0336,-9.0573,-9.0707,-9.0764,-9.0774,-9.0758,-9.0725,-9.0682,-9.0629,-9.0568,-9.0497,-9.0416,
                                                                    0.27938  ,0.27866 ,0.27846  ,0.27814 ,0.2773  ,0.27508  ,0.26919 ,0.25369  ,0.21311   ,0.10798  ,-0.15742 ,-0.78656 ,-2.0788   ,-4.0955   ,-6.2027   ,-7.6554   ,-8.3942    ,-8.7124  ,-8.8394   ,-8.8885 ,-8.9072  ,-8.9142 ,-8.9167  ,-8.9175 ,-8.9176  ,-8.9173 ,-8.917   ,-8.9166  ,-8.9161   ,-8.9157 ,-8.9154  ,-8.9154  ,-8.9164,-8.9192,-8.9259,-8.9395,-8.9627,-8.9945,-9.0273,-9.0523,-9.0667,-9.0729,-9.074 ,-9.0721,-9.0685,-9.0636,-9.0577,-9.0509,-9.0429,-9.0339,
                                                                    0.28005  ,0.27879 ,0.27858  ,0.27841 ,0.27797 ,0.2768   ,0.27367 ,0.26542  ,0.24371   ,0.18704  ,0.04115  ,-0.32131 ,-1.1484   ,-2.7179   ,-4.8697   ,-6.8024   ,-7.9807    ,-8.5384  ,-8.7705   ,-8.8618 ,-8.8969  ,-8.9101 ,-8.915   ,-8.9166 ,-8.917   ,-8.9169 ,-8.9166  ,-8.9161  ,-8.9156   ,-8.9151 ,-8.9146  ,-8.9145  ,-8.9151,-8.9174,-8.9233,-8.9354,-8.957 ,-8.9877,-9.0207,-9.0467,-9.0622,-9.0689,-9.0701,-9.068 ,-9.0639,-9.0585,-9.0519,-9.0442,-9.0354,-9.0252,
                                                                    0.28123  ,0.27898 ,0.27867  ,0.27856 ,0.27834 ,0.27773  ,0.27608 ,0.2717   ,0.26013   ,0.22976  ,0.15072  ,-0.050914,-0.54103  ,-1.6046   ,-3.4361   ,-5.6136   ,-7.3026    ,-8.2295  ,-8.644    ,-8.8122 ,-8.8777  ,-8.9027 ,-8.912   ,-8.9153 ,-8.9163  ,-8.9164 ,-8.9161  ,-8.9156  ,-8.915    ,-8.9144 ,-8.9138  ,-8.9135  ,-8.9138,-8.9156,-8.9206,-8.9315,-8.9514,-8.9807,-9.0136,-9.0407,-9.0573,-9.0645,-9.0658,-9.0634,-9.0589,-9.0528,-9.0454,-9.0368,-9.0269,-9.0156,
                                                                    0.28336  ,0.27928 ,0.27875  ,0.27866 ,0.27855 ,0.27824  ,0.27738 ,0.27507  ,0.26892   ,0.25271  ,0.21025  ,0.10028  ,-0.17686  ,-0.83084  ,-2.161    ,-4.2019   ,-6.2904    ,-7.7047  ,-8.4161   ,-8.7208 ,-8.8421  ,-8.8888 ,-8.9065  ,-8.913  ,-8.9152  ,-8.9157 ,-8.9155  ,-8.915   ,-8.9143   ,-8.9136 ,-8.9129  ,-8.9124  ,-8.9124,-8.9138,-8.9179,-8.9275,-8.9457,-8.9736,-9.0062,-9.0342,-9.0518,-9.0596,-9.0609,-9.0583,-9.0532,-9.0464,-9.0382,-9.0286,-9.0175,-9.0048,
                                                                    0.2872   ,0.27981 ,0.27885  ,0.27873 ,0.27868 ,0.27853  ,0.2781  ,0.27689  ,0.27364   ,0.26502  ,0.24232  ,0.18302  ,0.030486  ,-0.34755  ,-1.205    ,-2.8122   ,-4.9744    ,-6.8766  ,-8.0185   ,-8.5541 ,-8.7761  ,-8.8632 ,-8.8965  ,-8.9089 ,-8.9133  ,-8.9146 ,-8.9147  ,-8.9142  ,-8.9135   ,-8.9127 ,-8.9119  ,-8.9111  ,-8.9109,-8.9118,-8.9152,-8.9235,-8.94  ,-8.9664,-8.9985,-9.0271,-9.0457,-9.0542,-9.0555,-9.0525,-9.0468,-9.0392,-9.0301,-9.0193,-9.007 ,-8.9928,
                                                                    0.29412  ,0.28075 ,0.279    ,0.27879 ,0.27876 ,0.2787   ,0.2785  ,0.27789  ,0.27619   ,0.27163  ,0.25955  ,0.22778  ,0.1451    ,-0.065587 ,-0.57595  ,-1.6747   ,-3.5389    ,-5.7107  ,-7.3628   ,-8.2576 ,-8.655   ,-8.8157 ,-8.878   ,-8.9016 ,-8.9102  ,-8.9131 ,-8.9137  ,-8.9134  ,-8.9127   ,-8.9117 ,-8.9107  ,-8.9098  ,-8.9093,-8.9097,-8.9124,-8.9195,-8.9343,-8.9589,-8.9903,-9.0194,-9.039 ,-9.0481,-9.0495,-9.0461,-9.0398,-9.0313,-9.021 ,-9.009 ,-8.9952,-8.9794,
                                                                    0.30664  ,0.28242 ,0.27925  ,0.27886 ,0.27883 ,0.27882  ,0.27874 ,0.27845  ,0.27758   ,0.27519  ,0.26879  ,0.25186  ,0.20745   ,0.092451  ,-0.19688  ,-0.87642  ,-2.2446    ,-4.3081  ,-6.376    ,-7.7518 ,-8.4364  ,-8.728  ,-8.8438  ,-8.8882 ,-8.9047  ,-8.9106 ,-8.9123  ,-8.9123  ,-8.9116   ,-8.9106 ,-8.9095  ,-8.9083  ,-8.9075,-8.9075,-8.9094,-8.9154,-8.9286,-8.9513,-8.9817,-9.0112,-9.0317,-9.0413,-9.0428,-9.039 ,-9.0318,-9.0223,-9.0109,-8.9975,-8.982 ,-8.9644,
                                                                    0.32923  ,0.28544 ,0.27968  ,0.27896 ,0.27889 ,0.27891  ,0.2789  ,0.27878  ,0.27836   ,0.27713  ,0.27378  ,0.2648   ,0.24108   ,0.17907   ,0.019629  ,-0.37455  ,-1.263     ,-2.9075  ,-5.0778   ,-6.9483 ,-8.0541  ,-8.5683 ,-8.7805  ,-8.8634 ,-8.8949  ,-8.9064 ,-8.9101  ,-8.9109  ,-8.9104   ,-8.9094 ,-8.9081  ,-8.9067  ,-8.9056,-8.9051,-8.9063,-8.9112,-8.9227,-8.9436,-8.9727,-9.0022,-9.0235,-9.0338,-9.0353,-9.031 ,-9.023 ,-9.0124,-8.9995,-8.9846,-8.9674,-8.9477,
                                                                    0.36998  ,0.29087 ,0.28043  ,0.27909 ,0.27896 ,0.27898  ,0.27901 ,0.27899  ,0.27881   ,0.27822  ,0.27649  ,0.27177  ,0.25918   ,0.22599   ,0.13955   ,-0.080537 ,-0.6118    ,-1.7463  ,-3.642    ,-5.8058 ,-7.4204  ,-8.2836 ,-8.6644  ,-8.8178 ,-8.8769  ,-8.899  ,-8.9068  ,-8.9089  ,-8.9089   ,-8.9079 ,-8.9065  ,-8.9049  ,-8.9034,-8.9025,-8.903 ,-8.9068,-8.9166,-8.9356,-8.9632,-8.9926,-9.0146,-9.0254,-9.0269,-9.022 ,-9.0131,-9.0012,-8.9869,-8.9702,-8.951 ,-8.929 ,
                                                                    0.44322  ,0.30068 ,0.28176  ,0.27931 ,0.27904 ,0.27905  ,0.27911 ,0.27914  ,0.2791    ,0.27885  ,0.278    ,0.27556  ,0.26893   ,0.25128   ,0.20488   ,0.084673  ,-0.21728   ,-0.92307 ,-2.3293   ,-4.4135 ,-6.4588  ,-7.7962 ,-8.4548  ,-8.7337 ,-8.8439  ,-8.8858 ,-8.9011  ,-8.906   ,-8.9069   ,-8.9062 ,-8.9047  ,-8.9029  ,-8.9011,-8.8996,-8.8995,-8.9022,-8.9104,-8.9273,-8.9533,-8.9823,-9.0048,-9.0161,-9.0175,-9.012 ,-9.002 ,-8.9888,-8.9728,-8.9541,-8.9327,-8.9082,
                                                                    0.57416  ,0.31839 ,0.28414  ,0.27967 ,0.27914 ,0.27913  ,0.2792  ,0.27927  ,0.27931   ,0.27924  ,0.27886  ,0.27766  ,0.27423   ,0.26492   ,0.24019   ,0.17539   ,0.0088091  ,-0.40204 ,-1.3221   ,-3.0033 ,-5.1794  ,-7.017  ,-8.0871  ,-8.5805 ,-8.783   ,-8.8617 ,-8.8911  ,-8.9014  ,-8.9042   ,-8.904  ,-8.9026  ,-8.9006  ,-8.8985,-8.8966,-8.8957,-8.8973,-8.9039,-8.9188,-8.9429,-8.9711,-8.994 ,-9.0058,-9.0071,-9.0009,-8.9897,-8.9749,-8.957 ,-8.9362,-8.9123,-8.885 ,
                                                                    0.80595  ,0.35034 ,0.28842  ,0.28028 ,0.27927 ,0.27922  ,0.27929 ,0.27938  ,0.27947   ,0.27951  ,0.2794   ,0.27886  ,0.27715   ,0.27231   ,0.25925   ,0.22463   ,0.13432    ,-0.09547 ,-0.64824  ,-1.8188 ,-3.7448  ,-5.8983 ,-7.4747  ,-8.307  ,-8.6715  ,-8.8176 ,-8.8734  ,-8.8937  ,-8.9002   ,-8.9013 ,-8.9001  ,-8.898   ,-8.8956,-8.8931,-8.8915,-8.8921,-8.8972,-8.9099,-8.932 ,-8.9592,-8.9822,-8.9943,-8.9954,-8.9884,-8.9759,-8.9593,-8.9394,-8.9162,-8.8896,-8.8591,
                                                                    1.209    ,0.40785 ,0.29612  ,0.28135 ,0.27948 ,0.27932  ,0.27939 ,0.2795   ,0.27963   ,0.27973  ,0.27977  ,0.2796   ,0.27881   ,0.27638   ,0.26958   ,0.25125   ,0.20284    ,0.077292 ,-0.23768  ,-0.97036,-2.4147  ,-4.5175 ,-6.5384  ,-7.8373 ,-8.4702  ,-8.7366 ,-8.8411  ,-8.8803  ,-8.8939   ,-8.8975 ,-8.8971  ,-8.8951  ,-8.8923,-8.8894,-8.887 ,-8.8865,-8.89  ,-8.9006,-8.9205,-8.9463,-8.9692,-8.9815,-8.9824,-8.9745,-8.9605,-8.942 ,-8.9198,-8.8939,-8.8643,-8.8304,
                                                                    1.889    ,0.51093 ,0.31002  ,0.28325 ,0.2798  ,0.27944  ,0.2795  ,0.27963  ,0.27978   ,0.27993  ,0.28007  ,0.2801   ,0.27983   ,0.2787    ,0.27526   ,0.26568   ,0.23999    ,0.17238  ,-0.0015315,-0.42953,-1.3818  ,-3.0989 ,-5.2784  ,-7.0818 ,-8.1164  ,-8.5894 ,-8.7823  ,-8.8565  ,-8.8835   ,-8.892  ,-8.8933  ,-8.8916  ,-8.8886,-8.8852,-8.8821,-8.8805,-8.8824,-8.8908,-8.9083,-8.9325,-8.955 ,-8.9674,-8.9679,-8.959 ,-8.9433,-8.9227,-8.8979,-8.8691,-8.836 ,-8.7983,
                                                                    2.979    ,0.69425 ,0.33508  ,0.28664 ,0.28033 ,0.2796   ,0.27963 ,0.27977  ,0.27994   ,0.28013  ,0.28033  ,0.28048  ,0.2805    ,0.28008   ,0.27846   ,0.2736    ,0.26014    ,0.22415  ,0.12993   ,-0.10982,-0.68464 ,-1.8915 ,-3.8462  ,-5.987  ,-7.5247  ,-8.3264 ,-8.6749  ,-8.8134  ,-8.8654   ,-8.8834 ,-8.888   ,-8.8873  ,-8.8843,-8.8805,-8.8766,-8.8739,-8.8742,-8.8805,-8.8954,-8.9177,-8.9394,-8.9517,-8.9518,-8.9416,-8.9242,-8.9012,-8.8736,-8.8415,-8.8047,-8.7627,
                                                                    4.5912   ,1.0158  ,0.38022  ,0.29271 ,0.28122 ,0.27982  ,0.27977 ,0.27992  ,0.28012   ,0.28034  ,0.28058  ,0.28083  ,0.28102   ,0.281     ,0.28037   ,0.27805   ,0.27118    ,0.25227  ,0.20191   ,0.070955,-0.25735 ,-1.0174 ,-2.4995  ,-4.6187 ,-6.6132  ,-7.8735 ,-8.4811  ,-8.735   ,-8.8334   ,-8.8692 ,-8.8804  ,-8.8819  ,-8.8794,-8.8752,-8.8706,-8.8667,-8.8654,-8.8695,-8.8817,-8.9017,-8.9224,-8.9344,-8.9338,-8.9223,-8.9028,-8.8772,-8.8464,-8.8107,-8.7697,-8.7231,
                                                                    6.7126   ,1.5661  ,0.46127  ,0.30364 ,0.28276 ,0.28014  ,0.27994 ,0.28009  ,0.28031   ,0.28057  ,0.28085  ,0.28116  ,0.28146   ,0.28169   ,0.28162   ,0.28068   ,0.27736    ,0.26766  ,0.24112   ,0.17077 ,-0.010568,-0.45609,-1.4409  ,-3.193  ,-5.3732  ,-7.1411 ,-8.1404  ,-8.5929  ,-8.7758   ,-8.8449 ,-8.8687  ,-8.8746  ,-8.8733,-8.8691,-8.8639,-8.8589,-8.8559,-8.8576,-8.8671,-8.8846,-8.9038,-8.9152,-8.9138,-8.9008,-8.879 ,-8.8505,-8.8162,-8.7764,-8.7309,-8.6791,
                                                                    9.1125   ,2.4699  ,0.60596  ,0.32332 ,0.28547 ,0.28062  ,0.28015 ,0.28029  ,0.28053   ,0.28082  ,0.28114  ,0.2815   ,0.28189   ,0.28227   ,0.28253   ,0.28237   ,0.28099    ,0.27626  ,0.26258   ,0.22536 ,0.12731  ,-0.12253,-0.7198  ,-1.9629 ,-3.9446  ,-6.0701 ,-7.5683  ,-8.3395  ,-8.6716   ,-8.8021 ,-8.8495  ,-8.8641  ,-8.8657,-8.862 ,-8.8563,-8.8501,-8.8455,-8.8449,-8.8515,-8.8661,-8.8835,-8.8939,-8.8916,-8.8768,-8.8525,-8.8207,-8.7826,-8.7384,-8.6878,-8.6303,
                                                                    11.404   ,3.8577  ,0.8614   ,0.35876 ,0.29029 ,0.2814   ,0.28042 ,0.28051  ,0.28077   ,0.28109  ,0.28146  ,0.28187  ,0.28234   ,0.28283   ,0.2833    ,0.2836    ,0.2833     ,0.28127  ,0.27453   ,0.25526 ,0.20314  ,0.066854,-0.27495 ,-1.0628 ,-2.582   ,-4.7151 ,-6.6809  ,-7.9023  ,-8.4845   ,-8.7252 ,-8.8168  ,-8.8481  ,-8.8556,-8.8535,-8.8476,-8.8405,-8.8341,-8.8311,-8.8347,-8.8461,-8.8612,-8.8704,-8.8668,-8.8501,-8.8229,-8.7876,-8.7452,-8.6961,-8.6399,-8.5763,
                                                                    13.265   ,5.783   ,1.3038   ,0.42247 ,0.29892 ,0.2827   ,0.28077 ,0.28076  ,0.28104   ,0.2814   ,0.28181  ,0.28228  ,0.28281   ,0.2834    ,0.28403   ,0.28461   ,0.28493    ,0.28441  ,0.28145   ,0.27189 ,0.24477  ,0.17191 ,-0.016775,-0.48   ,-1.4975  ,-3.2831 ,-5.4611  ,-7.1917  ,-8.1554   ,-8.5871 ,-8.7593  ,-8.8221  ,-8.8411,-8.8428,-8.8375,-8.8296,-8.8216,-8.8162,-8.8165,-8.8246,-8.837 ,-8.8445,-8.8393,-8.8204,-8.7901,-8.7507,-8.7036,-8.6491,-8.5869,-8.5164,
                                                                    14.586   ,8.1096  ,2.0447   ,0.53647 ,0.31442 ,0.28492  ,0.28126 ,0.28105  ,0.28134   ,0.28174  ,0.2822   ,0.28273  ,0.28333   ,0.28401   ,0.28476   ,0.28554   ,0.28626    ,0.28661  ,0.28574   ,0.28145 ,0.26788  ,0.22978 ,0.12816  ,-0.13166,-0.75149 ,-2.0306 ,-4.0368  ,-6.1441  ,-7.6016   ,-8.3418 ,-8.6569  ,-8.778   ,-8.8192,-8.8287,-8.8255,-8.8172,-8.8077,-8.7998,-8.7968,-8.8012,-8.8104,-8.8158,-8.8088,-8.7874,-8.7535,-8.7098,-8.6575,-8.5971,-8.5282,-8.4502,
                                                                    15.438   ,10.494  ,3.2187   ,0.73878 ,0.34226 ,0.28878  ,0.28198 ,0.28139  ,0.28166   ,0.2821   ,0.28262  ,0.28321  ,0.28388   ,0.28465   ,0.28552   ,0.28647   ,0.28746    ,0.28834  ,0.28867   ,0.2873  ,0.28109  ,0.26188 ,0.20842  ,0.067141,-0.28803 ,-1.1035 ,-2.6589  ,-4.8027  ,-6.7373   ,-7.9187 ,-8.4747  ,-8.7011  ,-8.784 ,-8.8089,-8.8105,-8.8029,-8.7921,-8.7818,-8.7755,-8.7759,-8.7815,-8.7842,-8.775 ,-8.7507,-8.713 ,-8.6644,-8.6065,-8.5395,-8.4633,-8.3771,
                                                                    15.953   ,12.56   ,4.925    ,1.0923  ,0.39226 ,0.29557  ,0.28304 ,0.28173  ,0.28195   ,0.28243  ,0.283    ,0.28366  ,0.28442   ,0.28528   ,0.28626   ,0.28736   ,0.28857    ,0.28981  ,0.29088   ,0.29116 ,0.28901  ,0.28007 ,0.25291  ,0.17806 ,-0.017548,-0.49826,-1.548   ,-3.3652  ,-5.5374   ,-7.2286 ,-8.1556  ,-8.5653  ,-8.7248,-8.7791,-8.7911,-8.7861,-8.7747,-8.7621,-8.7523,-8.7485,-8.75  ,-8.7496,-8.7377,-8.7101,-8.6682,-8.6143,-8.5501,-8.476 ,-8.3918,-8.2968,
                                                                    16.253   ,14.103  ,7.117    ,1.6937  ,0.48167 ,0.30746  ,0.28456 ,0.28191  ,0.28199   ,0.28249  ,0.28312  ,0.28384  ,0.28468   ,0.28563   ,0.28672   ,0.28796   ,0.28935    ,0.29087  ,0.29242   ,0.29369 ,0.29383  ,0.29053 ,0.27769  ,0.23936 ,0.13483  ,-0.13448,-0.77651 ,-2.0906  ,-4.1184   ,-6.2038 ,-7.6187  ,-8.3266  ,-8.6228,-8.7321,-8.7644,-8.7659,-8.7553,-8.7408,-8.7276,-8.7193,-8.7164,-8.7123,-8.6972,-8.6659,-8.6193,-8.5596,-8.4885,-8.4067,-8.3139,-8.2093,
                                                                    16.422   ,15.133  ,9.5259   ,2.6711  ,0.64008 ,0.32791  ,0.2862  ,0.28107  ,0.28083   ,0.28131  ,0.28197  ,0.28273  ,0.28362   ,0.28464   ,0.28581   ,0.28715   ,0.28868    ,0.2904   ,0.29227   ,0.29416 ,0.29562  ,0.29545 ,0.29042  ,0.27198 ,0.21796  ,0.072424,-0.29552 ,-1.1381  ,-2.7277   ,-4.8785 ,-6.7784  ,-7.918   ,-8.446 ,-8.6555,-8.7269,-8.7421,-8.7349,-8.7195,-8.703 ,-8.6903,-8.6824,-8.6741,-8.6554,-8.6201,-8.5683,-8.5023,-8.424 ,-8.3338,-8.2317,-8.117 ,
                                                                    16.513   ,15.766  ,11.756   ,4.146   ,0.91562 ,0.36125  ,0.28552 ,0.27577  ,0.27483   ,0.27512  ,0.27567  ,0.27634  ,0.27715   ,0.2781    ,0.27921   ,0.28052   ,0.28204    ,0.28379  ,0.28577   ,0.28793 ,0.29007  ,0.29156 ,0.29071  ,0.28299 ,0.25647  ,0.18037 ,-0.021491,-0.51913 ,-1.6003   ,-3.4463 ,-5.6086  ,-7.2575  ,-8.146 ,-8.5313,-8.6753,-8.7177,-8.7192,-8.7049,-8.6861,-8.669 ,-8.6562,-8.6433,-8.6206,-8.5811,-8.524 ,-8.4514,-8.3654,-8.2665,-8.1547,-8.0293,
                                                                    16.547   ,16.122  ,13.504   ,6.1356  ,1.3788  ,0.40944  ,0.27227 ,0.25387  ,0.25135   ,0.251    ,0.25102  ,0.25119  ,0.25149   ,0.25195   ,0.25258   ,0.25342   ,0.25448    ,0.25581  ,0.25742   ,0.25931 ,0.26139  ,0.26338 ,0.2644   ,0.26208 ,0.25005  ,0.21168 ,0.10447  ,-0.17415 ,-0.83827  ,-2.1865 ,-4.2326  ,-6.2925  ,-7.6628,-8.3367,-8.6116,-8.7063,-8.7266,-8.7175,-8.698 ,-8.6774,-8.6597,-8.6422,-8.6156,-8.5718,-8.5093,-8.4302,-8.3364,-8.2286,-8.1069,-7.9706,
                                                                    16.517   ,16.278  ,14.674   ,8.4381  ,2.1209  ,0.4706   ,0.22358 ,0.18945  ,0.18398   ,0.18239  ,0.1814   ,0.18061  ,0.17998   ,0.17952   ,0.17927   ,0.17925   ,0.1795     ,0.18005  ,0.18096   ,0.18223 ,0.18386  ,0.18573 ,0.18742  ,0.18771 ,0.18324  ,0.16514 ,0.11024  ,-0.039787,-0.42076  ,-1.2906 ,-2.9131  ,-5.0672  ,-6.9295,-8.0255,-8.5234,-8.7134,-8.77  ,-8.7726,-8.7545,-8.7307,-8.7079,-8.6848,-8.6531,-8.604 ,-8.5352,-8.4483,-8.3453,-8.2274,-8.0944,-7.9458,
                                                                    16.43    ,16.296  ,15.355   ,10.688  ,3.278   ,0.58944  ,0.15011 ,0.089394 ,0.080447  ,0.078478 ,0.077572 ,0.076967 ,0.076584  ,0.076438  ,0.076557  ,0.076976  ,0.077734   ,0.078877 ,0.08045   ,0.082499,0.085053 ,0.088089,0.09143  ,0.094474,0.095495 ,0.089797,0.064948 ,-0.011164,-0.21788  ,-0.73099,-1.8415  ,-3.7123  ,-5.8589,-7.4608,-8.3067,-8.663 ,-8.7862,-8.8115,-8.7979,-8.7709,-8.7414,-8.7108,-8.6717,-8.6152,-8.5378,-8.4409,-8.3268,-8.1964,-8.0499,-7.8868,
                                                                    16.36    ,16.286  ,15.749   ,12.596  ,5.0031  ,0.90469  ,0.13342 ,0.026034 ,0.011898  ,0.010289 ,0.01049  ,0.011108 ,0.011982  ,0.013117  ,0.014544  ,0.016302  ,0.018434   ,0.02099  ,0.024025  ,0.027597,0.031763 ,0.036552,0.041923 ,0.047613,0.052768 ,0.054952,0.047578 ,0.013302 ,-0.092196 ,-0.37554,-1.0569  ,-2.4333  ,-4.4895,-6.5147,-7.8326,-8.4645,-8.7097,-8.7808,-8.7807,-8.7536,-8.7184,-8.68  ,-8.6331,-8.5685,-8.4821,-8.3748,-8.2489,-8.1056,-7.9451,-7.7672,
                                                                    16.335   ,16.295  ,15.992   ,14.018  ,7.2305  ,1.5348   ,0.20652 ,0.01417  ,-0.010766 ,-0.013163,-0.012421,-0.011078,-0.0094426,-0.007526 ,-0.0052959,-0.0027094,0.00028303 ,0.0037377,0.007718  ,0.012294,0.017539 ,0.023521,0.030274 ,0.037724,0.045488 ,0.052352,0.054887 ,0.043864 ,-0.0048487,-0.15202,-0.5385  ,-1.4275  ,-3.0734,-5.2187,-7.0316,-8.0733,-8.5301,-8.6891,-8.7197,-8.6993,-8.6608,-8.6154,-8.5606,-8.4879,-8.3923,-8.2745,-8.1365,-7.98  ,-7.8054,-7.6125,
                                                                    16.33    ,16.308  ,16.139   ,14.964  ,9.6314  ,2.5854   ,0.37616 ,0.031919 ,-0.013518 ,-0.018338,-0.017609,-0.015944,-0.013921 ,-0.011584 ,-0.0089013,-0.0058244,-0.0022981 ,0.0017409,0.0063645 ,0.011653,0.017696 ,0.024585,0.032402 ,0.04118 ,0.050799 ,0.060703,0.069168 ,0.071369 ,0.054396  ,-0.01521,-0.22002 ,-0.74163 ,-1.8737,-3.7591,-5.8782,-7.4225,-8.2158,-8.5326,-8.6245,-8.6229,-8.5855,-8.5344,-8.4717,-8.3905,-8.2853,-8.1563,-8.0059,-7.8358,-7.6465,-7.4384,
                                                                    16.329   ,16.318  ,16.224   ,15.547  ,11.8    ,4.1587   ,0.68463 ,0.074394 ,-0.0083728,-0.017967,-0.017658,-0.015816,-0.013512 ,-0.010849 ,-0.0077985,-0.0043092,-0.00031858,0.0042443,0.00946   ,0.01542 ,0.022227 ,0.029992,0.03883  ,0.04884 ,0.060049 ,0.072265,0.084686 ,0.094866 ,0.096022  ,0.070133,-0.028967,-0.31203 ,-1.0053,-2.4031,-4.4597,-6.4402,-7.6966,-8.2779,-8.4843,-8.5236,-8.4957,-8.4414,-8.3713,-8.2812,-8.1659,-8.0253,-7.8618,-7.6774,-7.4731,-7.2492,
                                                                    16.33    ,16.324  ,16.273   ,15.889  ,13.476  ,6.2501   ,1.216   ,0.15286  ,0.0028859 ,-0.015532,-0.01619 ,-0.01425 ,-0.011676 ,-0.0086871,-0.0052632,-0.0013483,0.0031268  ,0.0082414,0.014086  ,0.020762,0.028387 ,0.037089,0.047009 ,0.05829 ,0.071047 ,0.085295,0.10074  ,0.11621  ,0.12825   ,0.12735 ,0.088402 ,-0.051844,-0.4393,-1.3423,-3.004 ,-5.1306,-6.8845,-7.8633,-8.2705,-8.3902,-8.3874,-8.3356,-8.2596,-8.1606,-8.0348,-7.882 ,-7.7049,-7.5058,-7.286 ,-7.0462};

            return arrayCurvatureMin[y * arrayWidth + x];
        }
        case 6:
        {
            static const std::vector<double> arrayCurvatureMin   = {0.37825,0.064607,-0.21663 ,-0.59588  ,-1.2464  ,-2.336   ,-3.8139   ,-5.275   ,-6.3229  ,-6.9141   ,-7.2047  ,-7.3397  ,-7.4025  ,-7.434   ,-7.4524   ,-7.4657   ,-7.4775  ,-7.4892  ,-7.5015  ,-7.5148  ,-7.5294  ,-7.5457  ,-7.5639  ,-7.5842  ,-7.6066  ,-7.6308  ,-7.6564  ,-7.6827  ,-7.7092  ,-7.7354   ,-7.7608  ,-7.7849  ,-7.8076  ,-7.8284  ,-7.8473  ,-7.8642 ,-7.8791  ,-7.8921  ,-7.9032 ,-7.9125 ,-7.9201 ,-7.926 ,-7.9305,-7.9334,-7.9349,-7.9348,-7.9333,-7.9301,-7.9252,-7.9183,
                                                                    0.50177,0.17209 ,-0.074161,-0.33967  ,-0.75702 ,-1.5038  ,-2.7129   ,-4.2332  ,-5.6077  ,-6.523    ,-7.016   ,-7.253   ,-7.3625  ,-7.4141  ,-7.4408   ,-7.4571   ,-7.4696  ,-7.481   ,-7.4925  ,-7.5047  ,-7.5183  ,-7.5336  ,-7.551   ,-7.5707  ,-7.5927  ,-7.6165  ,-7.6419  ,-7.6682  ,-7.6949  ,-7.7215   ,-7.7474  ,-7.7723  ,-7.7957  ,-7.8174  ,-7.8372  ,-7.8551 ,-7.8709  ,-7.8847  ,-7.8965 ,-7.9064 ,-7.9146 ,-7.921 ,-7.9257,-7.9288,-7.9304,-7.9303,-7.9285,-7.925 ,-7.9195,-7.9119,
                                                                    0.63277,0.26811 ,0.025729 ,-0.18256  ,-0.45512 ,-0.93155 ,-1.7902   ,-3.1113  ,-4.6391  ,-5.904    ,-6.6915  ,-7.0994  ,-7.2921  ,-7.3811  ,-7.4237   ,-7.4464   ,-7.461   ,-7.4726  ,-7.4836  ,-7.495   ,-7.5077  ,-7.5221  ,-7.5387  ,-7.5577  ,-7.5791  ,-7.6025  ,-7.6275  ,-7.6537  ,-7.6805  ,-7.7073   ,-7.7337  ,-7.7592  ,-7.7834  ,-7.8059  ,-7.8267  ,-7.8454 ,-7.8621  ,-7.8767  ,-7.8893 ,-7.8999 ,-7.9086 ,-7.9154,-7.9205,-7.9238,-7.9254,-7.9253,-7.9233,-7.9193,-7.9132,-7.9048,
                                                                    0.77808,0.36417 ,0.1056   ,-0.081212 ,-0.27396 ,-0.57268 ,-1.127    ,-2.1081  ,-3.5239  ,-5.0228   ,-6.1632  ,-6.832   ,-7.1673  ,-7.3238  ,-7.3961   ,-7.4313   ,-7.4507  ,-7.4638  ,-7.4748  ,-7.4856  ,-7.4975  ,-7.5111  ,-7.5269  ,-7.5452  ,-7.5659  ,-7.5888  ,-7.6134  ,-7.6394  ,-7.6661  ,-7.6931   ,-7.7198  ,-7.7458  ,-7.7706  ,-7.794   ,-7.8155  ,-7.8352 ,-7.8528  ,-7.8682  ,-7.8816 ,-7.8929 ,-7.9021 ,-7.9095,-7.9149,-7.9184,-7.92  ,-7.9197,-7.9175,-7.9131,-7.9063,-7.897 ,
                                                                    0.94283,0.46756 ,0.17887  ,-0.0086943,-0.16338 ,-0.35829 ,-0.7002   ,-1.3494  ,-2.4572  ,-3.9417   ,-5.3774  ,-6.3865  ,-6.9482  ,-7.2224  ,-7.3492   ,-7.408    ,-7.4371  ,-7.4538  ,-7.4657  ,-7.4764  ,-7.4876  ,-7.5005  ,-7.5156  ,-7.5332  ,-7.5532  ,-7.5755  ,-7.5996  ,-7.6252  ,-7.6517  ,-7.6787   ,-7.7056  ,-7.732   ,-7.7574  ,-7.7815  ,-7.8039  ,-7.8244 ,-7.8429  ,-7.8592  ,-7.8733 ,-7.8853 ,-7.8952 ,-7.903 ,-7.9087,-7.9124,-7.9141,-7.9137,-7.9111,-7.9062,-7.8987,-7.8883,
                                                                    1.1313 ,0.58334 ,0.25363  ,0.05092   ,-0.091238,-0.23135 ,-0.44319  ,-0.84423 ,-1.6029  ,-2.8341   ,-4.3546  ,-5.6981  ,-6.5761  ,-7.0435  ,-7.2666   ,-7.3692   ,-7.417   ,-7.4412  ,-7.4558  ,-7.467   ,-7.4781  ,-7.4904  ,-7.5049  ,-7.5217  ,-7.541   ,-7.5626  ,-7.5861  ,-7.6112  ,-7.6374  ,-7.6643   ,-7.6913  ,-7.718   ,-7.7439  ,-7.7686  ,-7.7918  ,-7.8131 ,-7.8324  ,-7.8496  ,-7.8645 ,-7.8772 ,-7.8877 ,-7.8959,-7.902 ,-7.9059,-7.9077,-7.9071,-7.9042,-7.8987,-7.8903,-7.8789,
                                                                    1.3473 ,0.71549 ,0.33512  ,0.107     ,-0.038366,-0.15369 ,-0.29285  ,-0.53488 ,-1.0104  ,-1.8902   ,-3.2327  ,-4.7526  ,-5.9825  ,-6.7351  ,-7.1208   ,-7.3018   ,-7.3845  ,-7.4235  ,-7.4439  ,-7.4572  ,-7.4686  ,-7.4807  ,-7.4946  ,-7.5107  ,-7.5293  ,-7.5501  ,-7.573   ,-7.5975  ,-7.6233  ,-7.6499   ,-7.6769  ,-7.7038  ,-7.7301  ,-7.7553  ,-7.7792  ,-7.8013 ,-7.8214  ,-7.8394  ,-7.8551 ,-7.8685 ,-7.8796 ,-7.8883,-7.8948,-7.8989,-7.9006,-7.8999,-7.8965,-7.8904,-7.8811,-7.8684,
                                                                    1.5947 ,0.86749 ,0.42711  ,0.16525   ,0.0063442,-0.10218 ,-0.20435  ,-0.35383 ,-0.6387  ,-1.2037   ,-2.2116  ,-3.6446  ,-5.1268  ,-6.2299  ,-6.8669   ,-7.183    ,-7.3292  ,-7.396   ,-7.428   ,-7.4458  ,-7.4588  ,-7.4711  ,-7.4847  ,-7.5002  ,-7.5181  ,-7.5381  ,-7.5602  ,-7.5841  ,-7.6093  ,-7.6356   ,-7.6625  ,-7.6894  ,-7.716   ,-7.7417  ,-7.7661  ,-7.7889 ,-7.8098  ,-7.8286  ,-7.845  ,-7.8592 ,-7.8709 ,-7.8801,-7.8869,-7.8912,-7.8929,-7.8919,-7.8881,-7.8813,-7.871 ,-7.8569,
                                                                    1.8772 ,1.0427  ,0.53266  ,0.22942   ,0.049411 ,-0.063426,-0.14977  ,-0.24923 ,-0.41907 ,-0.75948  ,-1.428   ,-2.5649  ,-4.0595  ,-5.4704  ,-6.4417   ,-6.9747   ,-7.2323  ,-7.3502  ,-7.4044  ,-7.4315  ,-7.448   ,-7.4615  ,-7.4751  ,-7.4902  ,-7.5073  ,-7.5266  ,-7.5479  ,-7.571   ,-7.5956  ,-7.6215   ,-7.648   ,-7.6749  ,-7.7016  ,-7.7277  ,-7.7526  ,-7.776  ,-7.7976  ,-7.8172  ,-7.8344 ,-7.8492 ,-7.8615 ,-7.8713,-7.8784,-7.8828,-7.8845,-7.8832,-7.8789,-7.8712,-7.8598,-7.8442,
                                                                    2.1982 ,1.2445  ,0.65462  ,0.30233   ,0.094847 ,-0.029713,-0.11283  ,-0.18763 ,-0.29289 ,-0.4927   ,-0.90171 ,-1.6861  ,-2.9454  ,-4.4668  ,-5.779    ,-6.6201   ,-7.062   ,-7.271   ,-7.3664  ,-7.4112  ,-7.435   ,-7.4512  ,-7.4655  ,-7.4805  ,-7.497   ,-7.5155  ,-7.536   ,-7.5583  ,-7.5822  ,-7.6075   ,-7.6337  ,-7.6604  ,-7.6871  ,-7.7134  ,-7.7387  ,-7.7627 ,-7.7849  ,-7.8052  ,-7.8231 ,-7.8386 ,-7.8515 ,-7.8617,-7.8691,-7.8737,-7.8753,-7.8737,-7.8688,-7.8602,-7.8475,-7.8303,
                                                                    2.5612 ,1.4764  ,0.79581  ,0.38637   ,0.14537  ,0.003478 ,-0.084229 ,-0.14916 ,-0.22041 ,-0.33892  ,-0.57855 ,-1.0696  ,-1.9792  ,-3.3457  ,-4.8564   ,-6.0505   ,-6.7684  ,-7.1321  ,-7.3016  ,-7.3794  ,-7.4173  ,-7.4392  ,-7.4555  ,-7.4708  ,-7.487   ,-7.5048  ,-7.5245  ,-7.546   ,-7.5691  ,-7.5937   ,-7.6194  ,-7.6458  ,-7.6724  ,-7.6988  ,-7.7244  ,-7.7489 ,-7.7717  ,-7.7926  ,-7.8112 ,-7.8273 ,-7.8407 ,-7.8514,-7.8591,-7.8638,-7.8653,-7.8633,-7.8577,-7.8481,-7.834 ,-7.8149,
                                                                    2.9689 ,1.7418  ,0.95917  ,0.48381   ,0.20307  ,0.038976 ,-0.058658 ,-0.12248 ,-0.17747 ,-0.25153  ,-0.39043 ,-0.68031 ,-1.2668  ,-2.3063  ,-3.7564   ,-5.2198   ,-6.2849  ,-6.8904  ,-7.1883  ,-7.3263  ,-7.3906  ,-7.4236  ,-7.4443  ,-7.461   ,-7.4773  ,-7.4945  ,-7.5134  ,-7.5341  ,-7.5564  ,-7.5803   ,-7.6054  ,-7.6313  ,-7.6577  ,-7.6841  ,-7.7099  ,-7.7346 ,-7.7579  ,-7.7794  ,-7.7986 ,-7.8153 ,-7.8293 ,-7.8403,-7.8483,-7.8531,-7.8543,-7.852 ,-7.8456,-7.8348,-7.8192,-7.798 ,
                                                                    3.4237 ,2.0442  ,1.1478   ,0.597     ,0.2698   ,0.07875  ,-0.032977 ,-0.10121 ,-0.15019 ,-0.20132  ,-0.28367 ,-0.45025 ,-0.80162 ,-1.4963  ,-2.6643   ,-4.167    ,-5.551   ,-6.4845  ,-6.9904  ,-7.2341  ,-7.3469  ,-7.4009  ,-7.4303  ,-7.4503  ,-7.4674  ,-7.4845  ,-7.5028  ,-7.5226  ,-7.5441  ,-7.5671   ,-7.5915  ,-7.617   ,-7.643   ,-7.6692  ,-7.695   ,-7.72   ,-7.7437  ,-7.7656  ,-7.7853 ,-7.8026 ,-7.817  ,-7.8285,-7.8367,-7.8414,-7.8424,-7.8395,-7.8323,-7.8202,-7.8028,-7.7794,
                                                                    3.9269 ,2.3872  ,1.365    ,0.72844   ,0.3474   ,0.12437  ,-0.0051831,-0.081675,-0.13074 ,-0.17125  ,-0.22326 ,-0.31905 ,-0.52116 ,-0.94614 ,-1.7602   ,-3.0475   ,-4.5671  ,-5.8469  ,-6.6527  ,-7.0725  ,-7.2718  ,-7.3648  ,-7.4108  ,-7.4375  ,-7.457   ,-7.4745  ,-7.4924  ,-7.5115  ,-7.5321  ,-7.5543   ,-7.5779  ,-7.6027  ,-7.6283  ,-7.6542  ,-7.68    ,-7.7051 ,-7.729   ,-7.7512  ,-7.7714 ,-7.7891 ,-7.804  ,-7.8157,-7.8241,-7.8288,-7.8295,-7.8259,-7.8177,-7.8042,-7.7849,-7.7589,
                                                                    4.4788 ,2.7737  ,1.6143   ,0.88087   ,0.43778  ,0.17728  ,0.026155  ,-0.061704,-0.11475 ,-0.1517   ,-0.18834 ,-0.2452  ,-0.35967 ,-0.60606 ,-1.1175   ,-2.0589   ,-3.4483  ,-4.9478  ,-6.1068  ,-6.7937  ,-7.1403  ,-7.3034  ,-7.3806  ,-7.4203  ,-7.4449  ,-7.4641  ,-7.4821  ,-7.5006  ,-7.5205  ,-7.5418   ,-7.5646  ,-7.5887  ,-7.6137  ,-7.6392  ,-7.6647  ,-7.6898 ,-7.7138  ,-7.7363  ,-7.7569 ,-7.7749 ,-7.7901 ,-7.8021,-7.8105,-7.8151,-7.8154,-7.8111,-7.8017,-7.7866,-7.7652,-7.7364,
                                                                    5.0783 ,3.2063  ,1.899    ,1.0572    ,0.54305  ,0.23893  ,0.062221  ,-0.039872,-0.099692,-0.13733  ,-0.16706 ,-0.2033  ,-0.26866 ,-0.40755 ,-0.70817  ,-1.3192   ,-2.3915  ,-3.8576  ,-5.302   ,-6.3322  ,-6.9114  ,-7.1965  ,-7.3303  ,-7.3947  ,-7.4295  ,-7.4526  ,-7.4716  ,-7.49    ,-7.5092  ,-7.5297   ,-7.5517  ,-7.5749  ,-7.5992  ,-7.6242  ,-7.6493  ,-7.6742 ,-7.6983  ,-7.7209  ,-7.7416 ,-7.7599 ,-7.7754 ,-7.7875,-7.7959,-7.8003,-7.8001,-7.7949,-7.7842,-7.7673,-7.7435,-7.7116,
                                                                    5.7225 ,3.6868  ,2.2227   ,1.2606    ,0.66551  ,0.31089  ,0.10413   ,-0.01513 ,-0.084048,-0.12513  ,-0.15281 ,-0.17876 ,-0.21747 ,-0.29508 ,-0.46502  ,-0.8312   ,-1.5546  ,-2.7549  ,-4.2658  ,-5.6248  ,-6.5252  ,-7.0094  ,-7.2432  ,-7.3531  ,-7.4075  ,-7.4384  ,-7.4603  ,-7.4793  ,-7.4982  ,-7.5179   ,-7.539   ,-7.5614  ,-7.5849  ,-7.6092  ,-7.6339  ,-7.6584 ,-7.6823  ,-7.705   ,-7.7258 ,-7.7442 ,-7.7598 ,-7.772 ,-7.7803,-7.7843,-7.7834,-7.7773,-7.7651,-7.7462,-7.7196,-7.6843,
                                                                    6.407  ,4.2161  ,2.5885   ,1.4945    ,0.80773  ,0.39488  ,0.15305   ,0.013417 ,-0.066802,-0.11335  ,-0.14195 ,-0.1634  ,-0.18813 ,-0.23196 ,-0.32609  ,-0.53499  ,-0.97939 ,-1.8262  ,-3.1435  ,-4.6629  ,-5.9135  ,-6.6889  ,-7.0905  ,-7.282   ,-7.3727  ,-7.4191  ,-7.447   ,-7.468   ,-7.4871  ,-7.5064   ,-7.5266  ,-7.5481  ,-7.5708  ,-7.5943  ,-7.6184  ,-7.6425 ,-7.6661  ,-7.6885  ,-7.7093 ,-7.7277 ,-7.7433 ,-7.7554,-7.7635,-7.767 ,-7.7654,-7.758 ,-7.7442,-7.723 ,-7.6935,-7.6543,
                                                                    7.1256 ,4.7937  ,2.9994   ,1.7622    ,0.97252  ,0.49285  ,0.21023   ,0.04664  ,-0.047181,-0.10095  ,-0.13241 ,-0.15273 ,-0.17055 ,-0.19621 ,-0.24802  ,-0.36382  ,-0.62105 ,-1.1571  ,-2.1345  ,-3.5494  ,-5.0397  ,-6.167   ,-6.8262  ,-7.1576  ,-7.3143  ,-7.3897  ,-7.4297  ,-7.4553  ,-7.4757  ,-7.4949   ,-7.5145  ,-7.5351  ,-7.5569  ,-7.5796  ,-7.6029  ,-7.6264 ,-7.6495  ,-7.6716  ,-7.6922 ,-7.7105 ,-7.7259 ,-7.7377,-7.7454,-7.7483,-7.7458,-7.737 ,-7.7213,-7.6976,-7.6647,-7.6214,
                                                                    7.8704 ,5.4175  ,3.4574   ,2.0674    ,1.1629   ,0.60697  ,0.27709   ,0.085453 ,-0.024505,-0.087187 ,-0.123   ,-0.14423 ,-0.15912 ,-0.17533 ,-0.20409  ,-0.26728  ,-0.41098 ,-0.72724 ,-1.3683  ,-2.478   ,-3.9623  ,-5.3889  ,-6.3862  ,-6.9407  ,-7.2128  ,-7.3414  ,-7.4044  ,-7.4394  ,-7.4632  ,-7.4832   ,-7.5025  ,-7.5224  ,-7.5432  ,-7.565   ,-7.5874  ,-7.6102 ,-7.6327  ,-7.6543  ,-7.6745 ,-7.6925 ,-7.7075 ,-7.719 ,-7.7261,-7.7282,-7.7245,-7.7142,-7.6963,-7.6697,-7.6332,-7.5852,
                                                                    8.6323 ,6.084   ,3.964    ,2.4132    ,1.3822   ,0.73965  ,0.35523   ,0.13087  ,0.0019043,-0.07149  ,-0.11297 ,-0.13652 ,-0.1508  ,-0.1624  ,-0.17888  ,-0.21317  ,-0.29181 ,-0.47061 ,-0.85759 ,-1.6156  ,-2.8518  ,-4.3715  ,-5.7054  ,-6.5732  ,-7.0354  ,-7.2584  ,-7.364   ,-7.4172  ,-7.4482  ,-7.4707   ,-7.4904  ,-7.5097  ,-7.5297  ,-7.5505  ,-7.572   ,-7.5939 ,-7.6156  ,-7.6366  ,-7.6562 ,-7.6737 ,-7.6882 ,-7.6991,-7.7055,-7.7066,-7.7015,-7.6893,-7.669 ,-7.6392,-7.5986,-7.5455,
                                                                    9.4014 ,6.7875  ,4.5192   ,2.8028    ,1.6337   ,0.89358  ,0.44645   ,0.18404  ,0.032771 ,-0.053338 ,-0.1018  ,-0.12875 ,-0.14388 ,-0.15362 ,-0.1638   ,-0.18253  ,-0.22511 ,-0.32377 ,-0.54575 ,-1.0158  ,-1.9002  ,-3.249   ,-4.7664  ,-5.9864  ,-6.7307  ,-7.1134  ,-7.2959  ,-7.3831  ,-7.4284  ,-7.4563   ,-7.4777  ,-7.4971  ,-7.5164  ,-7.5362  ,-7.5567  ,-7.5776 ,-7.5984  ,-7.6185  ,-7.6374 ,-7.6541 ,-7.6679 ,-7.678 ,-7.6835,-7.6833,-7.6766,-7.6623,-7.6392,-7.6058,-7.5607,-7.5019,
                                                                    10.167 ,7.5212  ,5.1219   ,3.2386    ,1.921    ,1.0717   ,0.55281   ,0.24626  ,0.068905 ,-0.032204 ,-0.089042,-0.12041 ,-0.1374  ,-0.1469  ,-0.15411  ,-0.16468  ,-0.1877  ,-0.24146 ,-0.36518 ,-0.63938 ,-1.2052  ,-2.2211  ,-3.6601  ,-5.138   ,-6.2317  ,-6.8623  ,-7.1774  ,-7.3267  ,-7.3991  ,-7.4382   ,-7.4636  ,-7.484   ,-7.503   ,-7.522   ,-7.5415  ,-7.5612 ,-7.5809  ,-7.6001  ,-7.6179 ,-7.6338 ,-7.6466 ,-7.6557,-7.6599,-7.6583,-7.6497,-7.6329,-7.6066,-7.5693,-7.5191,-7.454 ,
                                                                    10.919 ,8.2764  ,5.7691   ,3.7225    ,2.2475   ,1.2772   ,0.6766    ,0.31903  ,0.11123  ,-0.0075159,-0.074317,-0.11112 ,-0.13079 ,-0.1411  ,-0.14724  ,-0.15373  ,-0.16635 ,-0.19556 ,-0.26351 ,-0.41795 ,-0.75455 ,-1.4286  ,-2.5755  ,-4.0746  ,-5.4797  ,-6.4424  ,-6.9712  ,-7.2297  ,-7.352   ,-7.4125   ,-7.4465  ,-7.4698  ,-7.4893  ,-7.5078  ,-7.5262  ,-7.5448 ,-7.5633  ,-7.5813  ,-7.598  ,-7.6126 ,-7.6243 ,-7.6321,-7.6348,-7.6314,-7.6206,-7.601 ,-7.5712,-7.5293,-7.4735,-7.4014,
                                                                    11.648 ,9.0434  ,6.4561   ,4.2552    ,2.6164   ,1.5133   ,0.82038   ,0.40404  ,0.1608   ,0.021372  ,-0.05721 ,-0.10055 ,-0.12369 ,-0.13556 ,-0.14179  ,-0.14648  ,-0.15373 ,-0.16973 ,-0.20691 ,-0.29234 ,-0.48414 ,-0.89446 ,-1.6878  ,-2.9578  ,-4.4816  ,-5.787   ,-6.621   ,-7.0608  ,-7.2722  ,-7.3726   ,-7.4235  ,-7.4533  ,-7.4748  ,-7.4933  ,-7.5109  ,-7.5284 ,-7.5456  ,-7.5621  ,-7.5775 ,-7.5907 ,-7.601  ,-7.6071,-7.6081,-7.6026,-7.5891,-7.5663,-7.5325,-7.4856,-7.4236,-7.3438,
                                                                    12.346 ,9.812   ,7.1768   ,4.8361    ,3.0305   ,1.7836   ,0.98696   ,0.50324  ,0.21883  ,0.055197  ,-0.03727 ,-0.088413,-0.11583 ,-0.1299  ,-0.13697  ,-0.14117  ,-0.14582 ,-0.15486 ,-0.17533 ,-0.2223  ,-0.32912 ,-0.56611 ,-1.0624  ,-1.9835  ,-3.3603  ,-4.8711  ,-6.0581  ,-6.7706  ,-7.1338  ,-7.3065   ,-7.3892  ,-7.4321  ,-7.4584  ,-7.4781  ,-7.4954  ,-7.5117 ,-7.5276  ,-7.5426  ,-7.5564 ,-7.568  ,-7.5765 ,-7.5808,-7.5796,-7.5716,-7.5552,-7.5287,-7.4903,-7.4379,-7.3689,-7.2806,
                                                                    13.006 ,10.572  ,7.9231   ,5.463     ,3.492    ,2.0915   ,1.1794    ,0.61878  ,0.28673  ,0.09481   ,-0.013989,-0.074386,-0.10696 ,-0.12386 ,-0.13234  ,-0.13684  ,-0.14044 ,-0.14592 ,-0.15744 ,-0.18337 ,-0.2423  ,-0.37525 ,-0.66661 ,-1.2616  ,-2.3143  ,-3.7733  ,-5.2348  ,-6.293   ,-6.8945  ,-7.1929   ,-7.3338  ,-7.4019  ,-7.4383  ,-7.4615  ,-7.4792  ,-7.4948 ,-7.5093  ,-7.5228  ,-7.5347 ,-7.5445 ,-7.5509 ,-7.553 ,-7.5493,-7.5384,-7.5186,-7.4879,-7.4444,-7.3857,-7.3091,-7.2114,
                                                                    13.622 ,11.313  ,8.6859   ,6.1323    ,4.002    ,2.4404   ,1.401     ,0.75312  ,0.36608  ,0.14119   ,0.013219 ,-0.058123,-0.096882,-0.11725 ,-0.12762  ,-0.13297  ,-0.1364  ,-0.14019 ,-0.14699 ,-0.16153 ,-0.19406 ,-0.26764 ,-0.43253 ,-0.78876 ,-1.4948  ,-2.6766  ,-4.1862  ,-5.5664  ,-6.4933  ,-6.9961   ,-7.2401  ,-7.3549  ,-7.4111  ,-7.4418  ,-7.4618  ,-7.4773 ,-7.4906  ,-7.5024  ,-7.5124 ,-7.52   ,-7.5242 ,-7.5237,-7.5171,-7.5028,-7.479 ,-7.4437,-7.3945,-7.3288,-7.2436,-7.1357,
                                                                    14.192 ,12.026  ,9.4551   ,6.8382    ,4.5607   ,2.8332   ,1.6551    ,0.90894  ,0.45871  ,0.19548   ,0.04503  ,-0.039227,-0.08534 ,-0.1099  ,-0.12263  ,-0.12923  ,-0.13303 ,-0.13618 ,-0.14056 ,-0.14899 ,-0.16714 ,-0.20774 ,-0.29932 ,-0.5031  ,-0.9359  ,-1.7635  ,-3.0645  ,-4.5882  ,-5.8623  ,-6.6616   ,-7.0783  ,-7.277   ,-7.3704  ,-7.4164  ,-7.442   ,-7.4587 ,-7.4713  ,-7.4815  ,-7.4895 ,-7.4947 ,-7.4962 ,-7.4927,-7.4828,-7.4647,-7.4364,-7.3957,-7.3401,-7.2666,-7.1721,-7.0527,
                                                                    14.713 ,12.704  ,10.22    ,7.5738    ,5.1666   ,3.2725   ,1.9453    ,1.0892   ,0.56668  ,0.25898   ,0.082229 ,-0.01724 ,-0.07207 ,-0.10163 ,-0.11724  ,-0.12544  ,-0.12998 ,-0.13307 ,-0.1363  ,-0.14148 ,-0.15184 ,-0.17436 ,-0.22487 ,-0.33859 ,-0.58953 ,-1.1114  ,-2.0679  ,-3.4694  ,-4.9696  ,-6.1211   ,-6.8009  ,-7.1436  ,-7.3048  ,-7.3802  ,-7.4175  ,-7.438  ,-7.4509  ,-7.4598  ,-7.4656 ,-7.4683 ,-7.4668 ,-7.46  ,-7.4463,-7.4239,-7.3905,-7.3438,-7.281 ,-7.1988,-7.0939,-6.9621,
                                                                    15.185 ,13.341  ,10.971   ,8.3302    ,5.8169   ,3.76     ,2.2749    ,1.2971   ,0.69229  ,0.33318   ,0.12573  ,0.0083716,-0.056763,-0.092254,-0.11129  ,-0.12145  ,-0.12701 ,-0.13042 ,-0.1332  ,-0.1367  ,-0.14285 ,-0.15551 ,-0.18338 ,-0.24611 ,-0.38705 ,-0.69471 ,-1.3184  ,-2.4061  ,-3.8813  ,-5.3223   ,-6.3429  ,-6.9142  ,-7.1941  ,-7.324   ,-7.3842  ,-7.4134 ,-7.4287  ,-7.437   ,-7.4408 ,-7.4407 ,-7.436  ,-7.4255,-7.4075,-7.3801,-7.341 ,-7.2876,-7.2167,-7.125 ,-7.0085,-6.863 ,
                                                                    15.61  ,13.933  ,11.699   ,9.0977    ,6.5066   ,4.2963   ,2.6471    ,1.5359   ,0.83813  ,0.4198    ,0.17657  ,0.038227 ,-0.039059,-0.081568,-0.10466  ,-0.11715  ,-0.12399 ,-0.12797 ,-0.13071 ,-0.13338 ,-0.13725 ,-0.14455 ,-0.16003 ,-0.19449 ,-0.27232 ,-0.44658 ,-0.82179 ,-1.5591  ,-2.7736  ,-4.2892   ,-5.6406  ,-6.5295  ,-7.0044  ,-7.2313  ,-7.3348  ,-7.3817 ,-7.4031  ,-7.4123  ,-7.4147 ,-7.4118 ,-7.4036 ,-7.3889,-7.3661,-7.3332,-7.2877,-7.2267,-7.1468,-7.0445,-6.9154,-6.7547,
                                                                    15.989 ,14.477  ,12.394   ,9.8661    ,7.2295   ,4.8807   ,3.0648    ,1.8092   ,1.007    ,0.52077   ,0.23598  ,0.073051 ,-0.018539,-0.069323,-0.097198 ,-0.11243  ,-0.1208  ,-0.12555 ,-0.12851 ,-0.13083 ,-0.13349 ,-0.13784 ,-0.14653 ,-0.16546 ,-0.20806 ,-0.30455 ,-0.5194  ,-0.97405 ,-1.8348  ,-3.1635   ,-4.6821  ,-5.9211  ,-6.6831  ,-7.0739  ,-7.2559  ,-7.3369 ,-7.3715  ,-7.3846  ,-7.3866 ,-7.3812 ,-7.3693 ,-7.3501,-7.322 ,-7.2829,-7.2302,-7.1607,-7.071 ,-6.9568,-6.8138,-6.6366,
                                                                    16.325 ,14.972  ,13.051   ,10.625    ,7.9775   ,5.511    ,3.53      ,2.1203   ,1.202    ,0.63828   ,0.30534  ,0.11368  ,0.0052919,-0.055229,-0.088726 ,-0.10718  ,-0.11737 ,-0.12307 ,-0.12643 ,-0.12868 ,-0.13071 ,-0.13343 ,-0.13836 ,-0.14875 ,-0.17192 ,-0.2246  ,-0.34406 ,-0.60802 ,-1.1546  ,-2.1447   ,-3.5665  ,-5.0504  ,-6.1625  ,-6.8062  ,-7.1242  ,-7.2682 ,-7.3292  ,-7.3518  ,-7.3556 ,-7.3485 ,-7.333  ,-7.3089,-7.275 ,-7.229 ,-7.1682,-7.0894,-6.9886,-6.8614,-6.703 ,-6.5078,
                                                                    16.621 ,15.419  ,13.664   ,11.365    ,8.7413   ,6.1833   ,4.0438    ,2.4727   ,1.4264   ,0.77477   ,0.38627  ,0.1611   ,0.033011 ,-0.03895 ,-0.079045 ,-0.10128  ,-0.11359 ,-0.12045 ,-0.12435 ,-0.12672 ,-0.12844 ,-0.13023 ,-0.13305 ,-0.13869 ,-0.15115 ,-0.17956 ,-0.24471 ,-0.39236 ,-0.7152  ,-1.3661   ,-2.4859  ,-3.9718  ,-5.3861  ,-6.3647  ,-6.901   ,-7.1564 ,-7.2674  ,-7.31    ,-7.32   ,-7.3128 ,-7.294  ,-7.2649,-7.2246,-7.1711,-7.1014,-7.0121,-6.8991,-6.7577,-6.5825,-6.3675,
                                                                    16.88  ,15.819  ,14.231   ,12.076    ,9.5108   ,6.8919   ,4.6063    ,2.8694   ,1.6837   ,0.93295   ,0.48057  ,0.21644  ,0.065294 ,-0.020089,-0.067923 ,-0.09458  ,-0.10939 ,-0.11761 ,-0.12219 ,-0.1248  ,-0.12641 ,-0.12768 ,-0.12926 ,-0.13222 ,-0.13871 ,-0.1537  ,-0.18856 ,-0.26908 ,-0.45118 ,-0.84382  ,-1.6102  ,-2.8529  ,-4.3676  ,-5.6834  ,-6.5287  ,-6.9692 ,-7.1705  ,-7.2521  ,-7.2766 ,-7.2726 ,-7.2517 ,-7.2176,-7.1706,-7.1088,-7.0293,-6.9285,-6.8021,-6.6449,-6.4512,-6.2147,
                                                                    17.106 ,16.174  ,14.748   ,12.752    ,10.276   ,7.6296   ,5.216     ,3.3126   ,1.9772   ,1.1158    ,0.59031  ,0.281    ,0.10293  ,0.0018215,-0.055085 ,-0.086926 ,-0.10466 ,-0.11449 ,-0.1199  ,-0.12285 ,-0.12447 ,-0.12544 ,-0.12628 ,-0.12764 ,-0.13073 ,-0.13824 ,-0.15631 ,-0.19912 ,-0.29856 ,-0.52245  ,-0.99669 ,-1.8871  ,-3.2374  ,-4.7425  ,-5.9386  ,-6.6558 ,-7.0115  ,-7.1654  ,-7.2197 ,-7.2254 ,-7.2048 ,-7.1664,-7.1124,-7.0416,-6.9515,-6.838 ,-6.6968,-6.5223,-6.3085,-6.0487,
                                                                    17.303 ,16.488  ,15.218   ,13.386    ,11.026   ,8.3875   ,5.8697    ,3.8041   ,2.3104   ,1.3264    ,0.71782  ,0.3563   ,0.14684  ,0.027325 ,-0.040213 ,-0.078127 ,-0.099287,-0.11101 ,-0.11741 ,-0.12082 ,-0.12254 ,-0.12333 ,-0.12368 ,-0.12405 ,-0.12512 ,-0.12834 ,-0.13703 ,-0.15884 ,-0.21142 ,-0.33407  ,-0.60827 ,-1.1763  ,-2.1948  ,-3.6286  ,-5.0862  ,-6.1493 ,-6.7469  ,-7.0277  ,-7.1389 ,-7.1664 ,-7.1511 ,-7.1102,-7.0494,-6.9692,-6.8673,-6.7401,-6.5826,-6.3892,-6.1535,-5.8683,
                                                                    17.473 ,16.764  ,15.639   ,13.975    ,11.751   ,9.1559   ,6.5626    ,4.3445   ,2.6865   ,1.5682    ,0.86566  ,0.44404  ,0.19808  ,0.057059 ,-0.022931 ,-0.067964 ,-0.093144,-0.10709 ,-0.11468 ,-0.11865 ,-0.12056 ,-0.12125 ,-0.12126 ,-0.12094 ,-0.12072 ,-0.12139 ,-0.1247  ,-0.13475 ,-0.1611  ,-0.22566  ,-0.37663 ,-0.71079 ,-1.3842  ,-2.5288  ,-4.0143  ,-5.3899 ,-6.3143  ,-6.8021  ,-7.016  ,-7.0872 ,-7.0868 ,-7.047 ,-6.9805,-6.8906,-6.7763,-6.6339,-6.4587,-6.2448,-5.9852,-5.6726,
                                                                    17.62  ,17.006  ,16.015   ,14.515    ,12.445   ,9.9246   ,7.2881    ,4.9329   ,3.1082   ,1.8447    ,1.0367   ,0.54617  ,0.25789  ,0.091766 ,-0.0028013,-0.056181 ,-0.086078,-0.10264 ,-0.11164 ,-0.1163  ,-0.11846 ,-0.11913 ,-0.11887 ,-0.11804 ,-0.11692 ,-0.11594 ,-0.11603 ,-0.11933 ,-0.13094 ,-0.16278  ,-0.24198 ,-0.4273  ,-0.83201 ,-1.6206  ,-2.8813  ,-4.3812 ,-5.6466  ,-6.4326  ,-6.8201 ,-6.9731 ,-7.005  ,-6.9736,-6.904 ,-6.8049,-6.6775,-6.5189,-6.3245,-6.0881,-5.8027,-5.4606,
                                                                    17.747 ,17.216  ,16.349   ,15.007    ,13.099   ,10.683   ,8.0382    ,5.5671   ,3.5775   ,2.1592    ,1.234    ,0.66492  ,0.32769  ,0.13231  ,0.020687  ,-0.042479 ,-0.077916,-0.097564,-0.10822 ,-0.11372 ,-0.11622 ,-0.11691 ,-0.11644 ,-0.11518 ,-0.11336 ,-0.11122 ,-0.10922 ,-0.10845 ,-0.11158 ,-0.12496  ,-0.1634  ,-0.26043 ,-0.4871  ,-0.97353 ,-1.884   ,-3.2418 ,-4.716   ,-5.8504  ,-6.5023 ,-6.7982 ,-6.8935 ,-6.8845,-6.8173,-6.7106,-6.5699,-6.3941,-6.1789,-5.9182,-5.605 ,-5.2312,
                                                                    17.856 ,17.398  ,16.642   ,15.451    ,13.71    ,11.422   ,8.8035    ,6.2429   ,4.0954   ,2.5151    ,1.4608   ,0.80273  ,0.40909  ,0.17969  ,0.04813   ,-0.026507 ,-0.068451,-0.091731,-0.10436 ,-0.11086 ,-0.11379 ,-0.11455 ,-0.11391 ,-0.11226 ,-0.10984 ,-0.10677 ,-0.10327 ,-0.09987 ,-0.097881,-0.1006   ,-0.11592 ,-0.16228 ,-0.2809  ,-0.55687 ,-1.1361  ,-2.1696 ,-3.5969  ,-5.0061  ,-5.9965 ,-6.5206 ,-6.7312 ,-6.7696,-6.7155,-6.6053,-6.4521,-6.2583,-6.0209,-5.7342,-5.391 ,-4.9833,
                                                                    17.949 ,17.556  ,16.899   ,15.848    ,14.273   ,12.131   ,9.5738    ,6.9546   ,4.662    ,2.9154    ,1.7206   ,0.96234  ,0.50396  ,0.23507  ,0.080224  ,-0.0078565,-0.057445,-0.085002,-0.09996 ,-0.10767 ,-0.11113 ,-0.11202 ,-0.11122 ,-0.10922 ,-0.10623 ,-0.10235 ,-0.097652,-0.092341,-0.087006,-0.083265 ,-0.085195,-0.10263 ,-0.15842 ,-0.30301 ,-0.63705 ,-1.3189 ,-2.4698  ,-3.9308  ,-5.2399 ,-6.0794 ,-6.4824 ,-6.6114,-6.5906,-6.4849,-6.3219,-6.1103,-5.8495,-5.535 ,-5.1597,-4.7159,
                                                                    18.029 ,17.692  ,17.124   ,16.201    ,14.788   ,12.805   ,10.339    ,7.6947   ,5.2755   ,3.3623    ,2.0168   ,1.1467   ,0.61438  ,0.29978  ,0.11778   ,0.013948  ,-0.044619,-0.077214,-0.09493 ,-0.10408 ,-0.1082  ,-0.10927 ,-0.10835 ,-0.106   ,-0.10245 ,-0.097795,-0.09203 ,-0.08517 ,-0.077402,-0.069417 ,-0.063195,-0.063794,-0.083441,-0.15039 ,-0.32597 ,-0.72735,-1.5191  ,-2.7729  ,-4.2267 ,-5.4063 ,-6.093  ,-6.3802,-6.428 ,-6.3427,-6.1759,-5.9478,-5.6632,-5.3192,-4.9099,-4.4278,
                                                                    18.097 ,17.809  ,17.319   ,16.513    ,15.255   ,13.437   ,11.088    ,8.4545   ,5.9327   ,3.8576    ,2.3529   ,1.3591   ,0.74271  ,0.37534  ,0.16172   ,0.039465  ,-0.029647,-0.068174,-0.089152,-0.10001 ,-0.10494 ,-0.10627 ,-0.10525 ,-0.10255 ,-0.098439,-0.093007,-0.08621 ,-0.077956,-0.068193,-0.057099 ,-0.045489,-0.035789,-0.034274,-0.056125,-0.13616 ,-0.34831,-0.82623 ,-1.7309  ,-3.0633 ,-4.4665 ,-5.4942 ,-6.0287,-6.203 ,-6.1667,-6.0082,-5.768 ,-5.46  ,-5.0855,-4.6404,-4.1179,
                                                                    18.156 ,17.909  ,17.488   ,16.787    ,15.674   ,14.023   ,11.813    ,9.224    ,6.6286   ,4.4015    ,2.7319   ,1.6028   ,0.89155  ,0.46352  ,0.21315   ,0.069342  ,-0.012146,-0.057658,-0.082491,-0.095393,-0.1013  ,-0.10297 ,-0.10189 ,-0.09883 ,-0.094138,-0.087904,-0.080061,-0.070449,-0.058875,-0.045211 ,-0.029622,-0.013072,0.0014623,0.0061992,-0.017678,-0.11289,-0.36762 ,-0.93034 ,-1.9444 ,-3.3219 ,-4.6319 ,-5.4916,-5.8751,-5.9363,-5.809 ,-5.5658,-5.2371,-4.8321,-4.3498,-3.7853,
                                                                    18.206 ,17.995  ,17.634   ,17.027    ,16.048   ,14.561   ,12.504    ,9.993    ,7.3566   ,4.9933    ,3.1566   ,1.8814   ,1.0638   ,0.56628  ,0.27331   ,0.10433   ,0.0083273,-0.045404,-0.074792,-0.090121,-0.097216,-0.099324,-0.098214,-0.0948  ,-0.089498,-0.082419,-0.073483,-0.062484,-0.049132,-0.033114 ,-0.01421 ,0.0074407,0.030694 ,0.051893 ,0.061396 ,0.035931,-0.076627,-0.38006 ,-1.0336 ,-2.1453 ,-3.526  ,-4.7039,-5.3847,-5.6173,-5.5613,-5.3331,-4.9905,-4.5565,-4.0365,-3.4287,
                                                                    18.248 ,18.069  ,17.759   ,17.236    ,16.379   ,15.051   ,13.157    ,10.751   ,8.1084   ,5.6305    ,3.629    ,2.1982   ,1.2625   ,0.68587  ,0.34364   ,0.14532   ,0.032294 ,-0.031106,-0.065871,-0.084085,-0.092608,-0.095276,-0.094181,-0.090411,-0.084466,-0.076487,-0.06639 ,-0.053931,-0.038751,-0.02041  ,0.0015549,0.027488 ,0.057285 ,0.08959  ,0.11991  ,0.1363  ,0.11007  ,-0.021953,-0.37982,-1.1264 ,-2.3137 ,-3.6502,-4.6626,-5.1569,-5.2363,-5.0557,-4.7133,-4.2553,-3.6986,-3.0469,
                                                                    18.285 ,18.132  ,17.868   ,17.417    ,16.67    ,15.492   ,13.765    ,11.488   ,8.8746   ,6.3089    ,4.15     ,2.5566   ,1.491    ,0.82477  ,0.42578   ,0.19331   ,0.060362 ,-0.014406,-0.055516,-0.077153,-0.087394,-0.090763,-0.089738,-0.085612,-0.078987,-0.070045,-0.058701,-0.044682,-0.02757 ,-0.0068265,0.018176 ,0.048089 ,0.083424 ,0.12413  ,0.1686   ,0.2113  ,0.23747  ,0.21181  ,0.058473,-0.35842,-1.1943 ,-2.4248,-3.6673,-4.4866,-4.7882,-4.7097,-4.3937,-3.9226,-3.333 ,-2.6384,
                                                                    18.316 ,18.186  ,17.961   ,17.574    ,16.926   ,15.887   ,14.326    ,12.196   ,9.645    ,7.0227    ,4.7195   ,2.9596   ,1.7528   ,0.98573  ,0.52161   ,0.24948   ,0.093235 ,0.0051144,-0.043477,-0.069173,-0.081474,-0.085715,-0.084826,-0.080348,-0.073003,-0.063026,-0.050336,-0.034633,-0.015447,0.0078449 ,0.036    ,0.069884 ,0.1104   ,0.15827  ,0.21353  ,0.27425 ,0.33363  ,0.3734   ,0.35037 ,0.17444 ,-0.30385,-1.2174,-2.4483,-3.5481,-4.1528,-4.256 ,-4.0122,-3.5491,-2.9353,-2.2011};

            return arrayCurvatureMin[y * arrayWidth + x];
        }
        case 7:
        {
            static const std::vector<double> arrayCurvatureMin   = {1.9666,1.4305,0.98418,0.60429,0.26735,-0.051368,-0.37854 ,-0.74252,-1.1705  ,-1.6818  ,-2.2772  ,-2.9306   ,-3.5916  ,-4.2031   ,-4.7226  ,-5.134   ,-5.4437  ,-5.67    ,-5.834    ,-5.9543  ,-6.0451    ,-6.1168  ,-6.1764  ,-6.2286  ,-6.2766  ,-6.3223  ,-6.3672  ,-6.4119  ,-6.4572  ,-6.5033  ,-6.5504  ,-6.5987   ,-6.6481  ,-6.6988  ,-6.7507  ,-6.8038  ,-6.8579  ,-6.9131 ,-6.9692  ,-7.0261  ,-7.0838  ,-7.1421  ,-7.2008  ,-7.2599  ,-7.3192  ,-7.3786  ,-7.438  ,-7.4972 ,-7.556  ,-7.6144 ,
                                                                    2.2085,1.6368,1.165  ,0.77089,0.43237,0.12692  ,-0.16932 ,-0.48234,-0.83975 ,-1.2675  ,-1.7822  ,-2.3809   ,-3.0331  ,-3.6859   ,-4.2827  ,-4.7843  ,-5.1781  ,-5.4726  ,-5.6871   ,-5.8424  ,-5.9564    ,-6.0427  ,-6.1112  ,-6.1684  ,-6.2189  ,-6.2655  ,-6.3101  ,-6.354   ,-6.398   ,-6.4425  ,-6.4879  ,-6.5344   ,-6.5821  ,-6.631   ,-6.6812  ,-6.7326  ,-6.7853  ,-6.839  ,-6.8938  ,-6.9496  ,-7.0062  ,-7.0636  ,-7.1216  ,-7.1802  ,-7.2392  ,-7.2984  ,-7.3578 ,-7.4172 ,-7.4764 ,-7.5353 ,
                                                                    2.4638,1.8517,1.3491 ,0.93419,0.5858 ,0.28298  ,0.0041881,-0.27368,-0.57603 ,-0.92973 ,-1.3594  ,-1.879    ,-2.4816  ,-3.1324   ,-3.7766  ,-4.3585  ,-4.8424  ,-5.2191  ,-5.4992   ,-5.7025  ,-5.8495    ,-5.9576  ,-6.0397  ,-6.1052  ,-6.1602  ,-6.209   ,-6.2543  ,-6.2979  ,-6.3409  ,-6.3841  ,-6.4279  ,-6.4727   ,-6.5186  ,-6.5657  ,-6.6141  ,-6.6638  ,-6.7147  ,-6.7668 ,-6.8202  ,-6.8746  ,-6.93    ,-6.9863  ,-7.0434  ,-7.1012  ,-7.1596  ,-7.2185  ,-7.2776 ,-7.3369 ,-7.3963 ,-7.4555 ,
                                                                    2.7346,2.078 ,1.5399 ,1.0988 ,0.73397,0.42534  ,0.15311  ,-0.10335,-0.36651 ,-0.66126 ,-1.0138  ,-1.4473   ,-1.9729  ,-2.5797   ,-3.2289  ,-3.8639  ,-4.4307  ,-4.8971  ,-5.2572   ,-5.5235  ,-5.7162    ,-5.8555  ,-5.958   ,-6.0363  ,-6.0989  ,-6.1519  ,-6.1991  ,-6.2432  ,-6.2858  ,-6.3279  ,-6.3703  ,-6.4135   ,-6.4576  ,-6.5029  ,-6.5494  ,-6.5973  ,-6.6464  ,-6.6969 ,-6.7486  ,-6.8014  ,-6.8554  ,-6.9105  ,-6.9665  ,-7.0233  ,-7.0809  ,-7.139   ,-7.1977 ,-7.2567 ,-7.3159 ,-7.3752 ,
                                                                    3.0227,2.3179,1.7402 ,1.2683 ,0.88164,0.56055  ,0.28629  ,0.040079,-0.1979  ,-0.4496  ,-0.7395  ,-1.0931   ,-1.5321  ,-2.0647   ,-2.6757  ,-3.3228  ,-3.948   ,-4.4995  ,-4.9487   ,-5.2928  ,-5.5459    ,-5.7285  ,-5.8605  ,-5.9578  ,-6.0324  ,-6.0924  ,-6.1434  ,-6.1892  ,-6.2321  ,-6.2737  ,-6.315   ,-6.3567   ,-6.3991  ,-6.4426  ,-6.4873  ,-6.5333  ,-6.5806  ,-6.6292 ,-6.6791  ,-6.7304  ,-6.7828  ,-6.8364  ,-6.891   ,-6.9467  ,-7.0032  ,-7.0605  ,-7.1184 ,-7.1768 ,-7.2357 ,-7.2948 ,
                                                                    3.3298,2.5732,1.9522 ,1.4455 ,1.0325 ,0.69353  ,0.41051  ,0.16578 ,-0.058469,-0.2814  ,-0.52454 ,-0.81201  ,-1.1687  ,-1.6145   ,-2.1547  ,-2.7699  ,-3.4142  ,-4.0291  ,-4.5651   ,-4.9973  ,-5.3258    ,-5.5663  ,-5.7395  ,-5.8646  ,-5.957   ,-6.0281  ,-6.0856  ,-6.1349  ,-6.1792  ,-6.221   ,-6.2617  ,-6.3022   ,-6.3431  ,-6.3849  ,-6.4278  ,-6.4719  ,-6.5173  ,-6.564  ,-6.6121  ,-6.6615  ,-6.7122  ,-6.7642  ,-6.8174  ,-6.8716  ,-6.9269  ,-6.9831  ,-7.04   ,-7.0977 ,-7.1559 ,-7.2146 ,
                                                                    3.6573,2.8457,2.1779 ,1.6326 ,1.1892 ,0.82801  ,0.53086  ,0.28078 ,0.061278 ,-0.14462 ,-0.35554 ,-0.59269  ,-0.87991 ,-1.2414   ,-1.6952  ,-2.2435  ,-2.8625  ,-3.5034  ,-4.1074   ,-4.6276  ,-5.043     ,-5.3566  ,-5.5851  ,-5.7492  ,-5.8678  ,-5.9557  ,-6.0235  ,-6.0787  ,-6.1262  ,-6.1693  ,-6.21    ,-6.2497   ,-6.2895  ,-6.3297  ,-6.3708  ,-6.413   ,-6.4565  ,-6.5013 ,-6.5475  ,-6.5951  ,-6.644   ,-6.6942  ,-6.7457  ,-6.7984  ,-6.8522  ,-6.9071  ,-6.9629 ,-7.0195 ,-7.0768 ,-7.1348 ,
                                                                    4.0064,3.1369,2.4188 ,1.8316 ,1.3542 ,0.96685  ,0.65117  ,0.39036 ,0.16864  ,-0.029458,-0.22021 ,-0.4218   ,-0.65525 ,-0.94415  ,-1.3119  ,-1.7747  ,-2.3312  ,-2.9537  ,-3.5904   ,-4.1828  ,-4.6872    ,-5.0862  ,-5.3852  ,-5.6022  ,-5.7578  ,-5.8703  ,-5.9539  ,-6.0187  ,-6.0717  ,-6.1175  ,-6.1593  ,-6.199    ,-6.2379  ,-6.2768  ,-6.3164  ,-6.3568  ,-6.3984  ,-6.4413 ,-6.4855  ,-6.5311  ,-6.5781  ,-6.6265  ,-6.6762  ,-6.7272  ,-6.7794  ,-6.8327  ,-6.8872 ,-6.9425 ,-6.9988 ,-7.0558 ,
                                                                    4.3783,3.4479,2.6765 ,2.044  ,1.5293 ,1.1123   ,0.77435  ,0.49845 ,0.2691   ,0.071689 ,-0.1084  ,-0.28683  ,-0.48149 ,-0.71328  ,-1.0056  ,-1.3809  ,-1.8535  ,-2.4183  ,-3.0437   ,-3.6753  ,-4.2556    ,-4.744   ,-5.1268  ,-5.4119  ,-5.6179  ,-5.7654  ,-5.8722  ,-5.9517  ,-6.0136  ,-6.0645  ,-6.1088  ,-6.1493   ,-6.188   ,-6.2261  ,-6.2643  ,-6.3031  ,-6.3429  ,-6.3839 ,-6.4261  ,-6.4697  ,-6.5148  ,-6.5612  ,-6.609   ,-6.6581  ,-6.7086  ,-6.7603  ,-6.8131 ,-6.8671 ,-6.922  ,-6.9778 ,
                                                                    4.7739,3.7801,2.9524 ,2.2713 ,1.7161 ,1.2663   ,0.90268  ,0.60804 ,0.36673  ,0.16451  ,-0.012205,-0.17726  ,-0.34587 ,-0.53573  ,-0.76767 ,-1.0649  ,-1.449   ,-1.9319  ,-2.5049   ,-3.1324  ,-3.7582    ,-4.3257  ,-4.798   ,-5.165   ,-5.4366  ,-5.6322  ,-5.772   ,-5.8733  ,-5.949   ,-6.0082  ,-6.0571  ,-6.1      ,-6.1394  ,-6.1771  ,-6.2143  ,-6.2518  ,-6.2899  ,-6.329  ,-6.3694  ,-6.411   ,-6.454   ,-6.4984  ,-6.5442  ,-6.5914  ,-6.64    ,-6.6899  ,-6.741  ,-6.7933 ,-6.8467 ,-6.9011 ,
                                                                    5.1939,4.1346,3.2477 ,2.515  ,1.9161 ,1.4303   ,1.038    ,0.72141 ,0.46457  ,0.25322  ,0.074291 ,-0.084899 ,-0.23756 ,-0.39855  ,-0.58549 ,-0.81923 ,-1.1227  ,-1.5165  ,-2.0102   ,-2.5911  ,-3.22      ,-3.8389  ,-4.3932  ,-4.8495  ,-5.201   ,-5.4596  ,-5.6452  ,-5.7778  ,-5.874   ,-5.946   ,-6.0026  ,-6.0497   ,-6.0911  ,-6.1294  ,-6.1662  ,-6.2026  ,-6.2393  ,-6.2767 ,-6.3152  ,-6.3548  ,-6.3958  ,-6.4382  ,-6.4819  ,-6.5271  ,-6.5737  ,-6.6216  ,-6.6709 ,-6.7214 ,-6.7731 ,-6.8259 ,
                                                                    5.6391,4.5124,3.5636 ,2.7763 ,2.1306 ,1.6058   ,1.1818   ,0.84037 ,0.56492  ,0.34094  ,0.15546  ,-0.0035708,-0.14803 ,-0.2906   ,-0.44595 ,-0.63166 ,-0.86866 ,-1.1796  ,-1.584    ,-2.0886  ,-2.6771    ,-3.3064  ,-3.9177  ,-4.4582  ,-4.8984  ,-5.2348  ,-5.4809  ,-5.657   ,-5.7827  ,-5.874   ,-5.9427  ,-5.9968   ,-6.0421  ,-6.0821  ,-6.1193  ,-6.1552  ,-6.1908  ,-6.2268 ,-6.2635  ,-6.3013  ,-6.3402  ,-6.3805  ,-6.4222  ,-6.4653  ,-6.5098  ,-6.5557  ,-6.6029 ,-6.6515 ,-6.7013 ,-6.7523 ,
                                                                    6.1097,4.9141,3.9011 ,3.0563 ,2.3608 ,1.794    ,1.3355   ,0.96642 ,0.66957  ,0.42999  ,0.2345   ,0.071245  ,-0.070861,-0.20303  ,-0.33757 ,-0.48899 ,-0.67498 ,-0.91656 ,-1.2361   ,-1.6517  ,-2.1674    ,-2.7628  ,-3.3916  ,-3.9943  ,-4.5207  ,-4.9449  ,-5.2665  ,-5.5006  ,-5.6677  ,-5.7869  ,-5.8736  ,-5.939    ,-5.9908  ,-6.0344  ,-6.0731  ,-6.1092  ,-6.1442  ,-6.179  ,-6.2142  ,-6.2501  ,-6.2872  ,-6.3255  ,-6.365   ,-6.406   ,-6.4483  ,-6.4921  ,-6.5372 ,-6.5837 ,-6.6314 ,-6.6805 ,
                                                                    6.6059,5.3405,4.2612 ,3.3562 ,2.6079 ,1.9961   ,1.5002   ,1.1007  ,0.77991  ,0.5221   ,0.3137   ,0.1428    ,-0.001403,-0.12919  ,-0.25122 ,-0.37951 ,-0.52855 ,-0.71619 ,-0.96352  ,-1.2926  ,-1.72      ,-2.2467  ,-2.8484  ,-3.4757  ,-4.069   ,-4.5808  ,-4.9891  ,-5.2963  ,-5.5189  ,-5.6773  ,-5.7903  ,-5.8726   ,-5.9349  ,-5.9846  ,-6.0265  ,-6.0639  ,-6.099   ,-6.133  ,-6.167   ,-6.2014  ,-6.2366  ,-6.2729  ,-6.3103  ,-6.3491  ,-6.3893  ,-6.4308  ,-6.4737 ,-6.518  ,-6.5635 ,-6.6104 ,
                                                                    7.1273,5.7915,4.6442 ,3.6765 ,2.8726 ,2.2129   ,1.6768   ,1.2442  ,0.89689  ,0.61843  ,0.39465  ,0.21329   ,0.063569 ,-0.06436  ,-0.1801  ,-0.29385 ,-0.41748 ,-0.56551 ,-0.756    ,-1.0102  ,-1.3496    ,-1.7891  ,-2.3267  ,-2.9339  ,-3.5586  ,-4.1416  ,-4.6385  ,-5.0311  ,-5.3244  ,-5.5358  ,-5.686   ,-5.7931   ,-5.8713  ,-5.9306  ,-5.9781  ,-6.0184  ,-6.0546  ,-6.0886 ,-6.1217  ,-6.1547  ,-6.1883  ,-6.2227  ,-6.2581  ,-6.2947  ,-6.3327  ,-6.3719  ,-6.4125 ,-6.4544 ,-6.4977 ,-6.5422 ,
                                                                    7.6723,6.2665,5.05   ,4.0175 ,3.1551 ,2.4446   ,1.8654   ,1.3972  ,1.0209   ,0.71946  ,0.47806  ,0.2839    ,0.12594  ,-0.0055521,-0.11961 ,-0.22525 ,-0.33235 ,-0.45268 ,-0.60091  ,-0.79533 ,-1.0573    ,-1.4077  ,-1.8597  ,-2.4078  ,-3.0196  ,-3.6407  ,-4.2126  ,-4.6943  ,-5.0713  ,-5.351   ,-5.5517  ,-5.694    ,-5.7955  ,-5.8696  ,-5.9262  ,-5.9716  ,-6.0103  ,-6.0451 ,-6.078   ,-6.1101  ,-6.1422  ,-6.1748  ,-6.2083  ,-6.2428  ,-6.2785  ,-6.3154  ,-6.3536 ,-6.393  ,-6.4338 ,-6.4758 ,
                                                                    8.2378,6.763 ,5.4767 ,4.3774 ,3.4543 ,2.6902   ,2.0651   ,1.5586  ,1.151    ,0.8244   ,0.56334  ,0.35436   ,0.18598  ,0.048297  ,-0.067491,-0.16971 ,-0.2669  ,-0.36871 ,-0.48689  ,-0.63634 ,-0.83563   ,-1.1061  ,-1.4682  ,-1.9328  ,-2.4908  ,-3.1062  ,-3.7225  ,-4.2825  ,-4.7488  ,-5.1103  ,-5.3767  ,-5.5672   ,-5.7019  ,-5.7979  ,-5.8682  ,-5.9219  ,-5.9652  ,-6.0023 ,-6.0358  ,-6.0675  ,-6.0985  ,-6.1295  ,-6.1611  ,-6.1935  ,-6.2269  ,-6.2614  ,-6.2971 ,-6.334  ,-6.3721 ,-6.4115 ,
                                                                    8.8171,7.2753,5.9193 ,4.7521 ,3.766  ,2.9459   ,2.2724   ,1.7251  ,1.2838   ,0.93004  ,0.64746  ,0.42192   ,0.24135  ,0.095506  ,-0.024483,-0.12662 ,-0.21863 ,-0.30868 ,-0.40627  ,-0.52324 ,-0.67472   ,-0.87962 ,-1.1593  ,-1.5333  ,-2.0104  ,-2.5777  ,-3.1955  ,-3.8058  ,-4.3531  ,-4.8036  ,-5.1497  ,-5.4032   ,-5.5836  ,-5.7109  ,-5.8016  ,-5.8681  ,-5.919   ,-5.9601 ,-5.9954  ,-6.0274  ,-6.0577  ,-6.0874  ,-6.1172  ,-6.1475  ,-6.1786  ,-6.2106  ,-6.2437 ,-6.2779 ,-6.3132 ,-6.3496 ,
                                                                    9.4016,7.7955,6.3704 ,5.1346 ,4.0838 ,3.2053   ,2.4808   ,1.8902  ,1.413    ,1.0301   ,0.72426  ,0.48058   ,0.28636  ,0.13082   ,0.0048412,-0.099559,-0.18968 ,-0.27273 ,-0.35662  ,-0.45078 ,-0.56718   ,-0.72129 ,-0.93232 ,-1.2216  ,-1.6077  ,-2.0969  ,-2.6724  ,-3.2912  ,-3.8942  ,-4.4279  ,-4.8622  ,-5.193    ,-5.4336  ,-5.6041  ,-5.7241  ,-5.8094  ,-5.8719  ,-5.9198 ,-5.9585  ,-5.9918  ,-6.022   ,-6.0506  ,-6.0787  ,-6.1068  ,-6.1355  ,-6.1648  ,-6.1951 ,-6.2263 ,-6.2585 ,-6.2918 ,
                                                                    9.9906,8.3228,6.8294 ,5.5239 ,4.4063 ,3.4665   ,2.688    ,2.0511  ,1.5352   ,1.1205   ,0.78917  ,0.52543   ,0.31582  ,0.14895   ,0.015269 ,-0.09338 ,-0.18417 ,-0.26379 ,-0.33909  ,-0.41774 ,-0.50919   ,-0.62564 ,-0.78291 ,-1.0005  ,-1.2998  ,-1.6979  ,-2.1986  ,-2.7814  ,-3.3997  ,-3.9938  ,-4.5129  ,-4.9306   ,-5.246   ,-5.4737  ,-5.6343  ,-5.7467  ,-5.8265  ,-5.8847 ,-5.9291  ,-5.965   ,-5.9958  ,-6.0237  ,-6.0501  ,-6.0761  ,-6.102   ,-6.1284  ,-6.1555 ,-6.1833 ,-6.212  ,-6.2417 ,
                                                                    10.601,8.8738,7.3128 ,5.9364 ,4.7493 ,3.7448   ,2.9083   ,2.2211  ,1.6626   ,1.2126   ,0.85254  ,0.56577   ,0.33805  ,0.15732   ,0.013479 ,-0.10198 ,-0.1963  ,-0.27599 ,-0.34722  ,-0.41652 ,-0.49144   ,-0.58152 ,-0.69923 ,-0.86077 ,-1.086   ,-1.396   ,-1.8067  ,-2.3188  ,-2.9079  ,-3.5246  ,-4.1089  ,-4.6128   ,-5.0139  ,-5.3138  ,-5.5289  ,-5.6796  ,-5.7846  ,-5.8586 ,-5.9123  ,-5.953   ,-5.9856  ,-6.0134  ,-6.0385  ,-6.0621  ,-6.0853  ,-6.1084  ,-6.1318 ,-6.1558 ,-6.1804 ,-6.2058 ,
                                                                    11.247,9.4653,7.8387 ,6.3915 ,5.1333 ,4.0613   ,3.1636   ,2.4224  ,1.8176   ,1.3288   ,0.93667  ,0.62382   ,0.37518  ,0.17793   ,0.021327 ,-0.10358 ,-0.20431 ,-0.28734 ,-0.35857  ,-0.42376 ,-0.48917   ,-0.56234 ,-0.6529  ,-0.77366 ,-0.94124 ,-1.1758  ,-1.498   ,-1.9223  ,-2.4464  ,-3.0419  ,-3.6566  ,-4.231    ,-4.7201  ,-5.1052  ,-5.3908  ,-5.5943  ,-5.736   ,-5.8343 ,-5.9031  ,-5.9527  ,-5.9899  ,-6.0195  ,-6.0444  ,-6.0666  ,-6.0874  ,-6.1076  ,-6.1275 ,-6.1477 ,-6.1682 ,-6.1892 ,
                                                                    11.926,10.096,8.4076 ,6.8911 ,5.5617 ,4.4208   ,3.4594   ,2.6616  ,2.0078   ,1.4776   ,1.051    ,0.70997   ,0.43852  ,0.22303   ,0.052083 ,-0.083841,-0.19265 ,-0.281   ,-0.3547   ,-0.41913 ,-0.4797    ,-0.54253 ,-0.61513 ,-0.70736 ,-0.83235 ,-1.0072  ,-1.2522  ,-1.5878  ,-2.0265  ,-2.5628  ,-3.1644  ,-3.7767   ,-4.3411  ,-4.8158  ,-5.1858  ,-5.4582  ,-5.6513  ,-5.7852 ,-5.8778  ,-5.9424  ,-5.9888  ,-6.0234  ,-6.0507  ,-6.0734  ,-6.0935  ,-6.112   ,-6.1297 ,-6.1471 ,-6.1644 ,-6.1819 ,
                                                                    12.629,10.757,9.011  ,7.4277 ,6.0276 ,4.8169   ,3.79     ,2.9332  ,2.228    ,1.6539   ,1.1908   ,0.81966   ,0.52383  ,0.28883   ,0.10245  ,-0.045444,-0.16326 ,-0.25796 ,-0.33545  ,-0.40097 ,-0.45944   ,-0.51605 ,-0.57677 ,-0.6492  ,-0.74343 ,-0.87301 ,-1.0555  ,-1.3113  ,-1.6604  ,-2.1133  ,-2.6609  ,-3.2673   ,-3.876   ,-4.4293  ,-4.8892  ,-5.2443  ,-5.5038  ,-5.6868 ,-5.8132  ,-5.9003  ,-5.9609  ,-6.0042  ,-6.0364  ,-6.0615  ,-6.0821  ,-6.1001  ,-6.1165 ,-6.1319 ,-6.1469 ,-6.1617 ,
                                                                    13.346,11.44 ,9.6418 ,7.9947 ,6.525  ,5.2439   ,4.1498   ,3.2317  ,2.4724   ,1.8519   ,1.3498   ,0.9466    ,0.62464  ,0.36867   ,0.16567  ,0.0047921,-0.12292 ,-0.22485 ,-0.30715  ,-0.37507 ,-0.43334   ,-0.48653 ,-0.53962 ,-0.59855 ,-0.67102 ,-0.76742 ,-0.90174 ,-1.0919  ,-1.3587  ,-1.721   ,-2.1875  ,-2.7453   ,-3.3549  ,-3.9581  ,-4.499   ,-4.9433  ,-5.2829  ,-5.5293 ,-5.7018  ,-5.8204  ,-5.9015  ,-5.9576  ,-5.9972  ,-6.0261  ,-6.0483  ,-6.0662  ,-6.0814 ,-6.0949 ,-6.1075 ,-6.1196 ,
                                                                    14.073,12.139,10.296 ,8.5885 ,7.0507 ,5.699    ,4.5363   ,3.5546  ,2.7385   ,2.0689   ,1.5253   ,1.0876    ,0.73752  ,0.45887   ,0.23783  ,0.062789 ,-0.075843,-0.18595 ,-0.27403  ,-0.34549 ,-0.40502   ,-0.4569  ,-0.50546 ,-0.5555  ,-0.61297 ,-0.68574 ,-0.78454 ,-0.9238  ,-1.1219  ,-1.3997  ,-1.7752  ,-2.2545   ,-2.8213  ,-3.4323  ,-4.0283  ,-4.5555  ,-4.9833  ,-5.3071 ,-5.5399  ,-5.7017  ,-5.8119  ,-5.8865  ,-5.9373  ,-5.9724  ,-5.9973  ,-6.0158  ,-6.03   ,-6.0417 ,-6.0517 ,-6.0607 ,
                                                                    14.806,12.852,10.969 ,9.2065 ,7.6031 ,6.1813   ,4.9488   ,3.9013  ,3.0259   ,2.3044   ,1.7165   ,1.2418    ,0.86139  ,0.5582    ,0.31753  ,0.127    ,-0.023681,-0.14296 ,-0.23775  ,-0.31376 ,-0.37575   ,-0.42793 ,-0.47423 ,-0.51873 ,-0.56615 ,-0.62249 ,-0.69583 ,-0.79726 ,-0.94171 ,-1.148   ,-1.4368  ,-1.8253   ,-2.3168  ,-2.8912  ,-3.5019  ,-4.0891  ,-4.6015  ,-5.0122 ,-5.3198  ,-5.5388  ,-5.6895  ,-5.791   ,-5.8586  ,-5.9035  ,-5.9335  ,-5.9539  ,-5.968  ,-5.9781 ,-5.9858 ,-5.9919 ,
                                                                    15.54 ,13.576,11.66  ,9.847  ,8.1811 ,6.6901   ,5.3873   ,4.2723  ,3.335    ,2.5587   ,1.9237   ,1.4095    ,0.99638  ,0.66662   ,0.40462  ,0.19717  ,0.033244 ,-0.096233,-0.19867  ,-0.28012 ,-0.34558   ,-0.39929 ,-0.445   ,-0.48636 ,-0.52727 ,-0.57242 ,-0.62789 ,-0.702   ,-0.80626 ,-0.9561  ,-1.1707  ,-1.4707   ,-1.8718  ,-2.3748  ,-2.9554  ,-3.5641  ,-4.141   ,-4.6374 ,-5.0304  ,-5.3214  ,-5.5265  ,-5.6659  ,-5.7584  ,-5.8186  ,-5.8574  ,-5.8821  ,-5.8976 ,-5.9073 ,-5.9132 ,-5.917  ,
                                                                    16.272,14.305,12.365 ,10.508 ,8.7836 ,7.2253   ,5.8521   ,4.6681  ,3.6665   ,2.8327   ,2.1478   ,1.5913    ,1.1431   ,0.78465   ,0.49951  ,0.27362  ,0.095184 ,-0.045551,-0.15655  ,-0.24429 ,-0.31407   ,-0.37025 ,-0.41661 ,-0.45655 ,-0.49347 ,-0.53114 ,-0.57425 ,-0.629   ,-0.70401 ,-0.81122 ,-0.96656 ,-1.1896   ,-1.5007  ,-1.9141  ,-2.4275  ,-3.0128  ,-3.6178  ,-4.1828 ,-4.6623  ,-5.0371  ,-5.3115  ,-5.5025  ,-5.6306  ,-5.714   ,-5.7669  ,-5.7996  ,-5.819  ,-5.8299 ,-5.8355 ,-5.8378 ,
                                                                    16.999,15.039,13.082 ,11.188 ,9.4093 ,7.7864   ,6.3432   ,5.0893  ,4.0215   ,3.1276   ,2.3899   ,1.7884    ,1.3025   ,0.91308   ,0.60287  ,0.35693  ,0.16265  ,0.0095513,-0.11093  ,-0.20576 ,-0.2806    ,-0.34006 ,-0.38799 ,-0.42774 ,-0.46241 ,-0.49521 ,-0.52984 ,-0.571   ,-0.62508 ,-0.70101 ,-0.81119 ,-0.97208  ,-1.2035  ,-1.5256  ,-1.9508  ,-2.4736  ,-3.062   ,-3.6616 ,-4.2133  ,-4.6751  ,-5.0314  ,-5.2891  ,-5.4663  ,-5.5834  ,-5.6581  ,-5.7041  ,-5.7311 ,-5.7459 ,-5.7531 ,-5.7556 ,
                                                                    17.718,15.773,13.808 ,11.884 ,10.057 ,8.3726   ,6.8609   ,5.5366  ,4.4008   ,3.4443   ,2.6511   ,2.0017    ,1.4755   ,1.0528    ,0.7155   ,0.44781  ,0.23627  ,0.069667 ,-0.061239 ,-0.16395 ,-0.24454   ,-0.30794 ,-0.35816 ,-0.39861 ,-0.43226 ,-0.46198 ,-0.49083 ,-0.52246 ,-0.56164 ,-0.61497 ,-0.69174 ,-0.80484  ,-0.97124 ,-1.2111  ,-1.5439  ,-1.9803  ,-2.5115  ,-3.1015 ,-3.6941  ,-4.2314  ,-4.6749  ,-5.0129  ,-5.2545  ,-5.4185  ,-5.5253  ,-5.5922  ,-5.6322 ,-5.6547 ,-5.6662 ,-5.6709 ,
                                                                    18.427,16.504,14.539 ,12.593 ,10.724 ,8.983    ,7.4049   ,6.0103  ,4.8053   ,3.784    ,2.9326   ,2.2325    ,1.6633   ,1.2048    ,0.83831  ,0.54707  ,0.3168   ,0.13547  ,-0.0068338,-0.11822 ,-0.20522   ,-0.27314 ,-0.32623 ,-0.36803 ,-0.40151 ,-0.42938 ,-0.4543  ,-0.47921 ,-0.50775 ,-0.54482 ,-0.59725 ,-0.67473  ,-0.79068 ,-0.96259 ,-1.2108  ,-1.5542  ,-2.0014  ,-2.5399 ,-3.1303  ,-3.7146  ,-4.2368  ,-4.6622  ,-4.9825  ,-5.2089  ,-5.361   ,-5.4589  ,-5.5193 ,-5.5549 ,-5.5745 ,-5.5843 ,
                                                                    19.122,17.23 ,15.274 ,13.314 ,11.41  ,9.6165   ,7.9747   ,6.5107  ,5.2357   ,4.1477   ,3.2355   ,2.482     ,1.867    ,1.3703    ,0.9723   ,0.65562  ,0.40504  ,0.20773  ,0.052999  ,-0.067863,-0.16193   ,-0.2349  ,-0.29135 ,-0.335   ,-0.36893 ,-0.39578 ,-0.418   ,-0.4381  ,-0.45897 ,-0.48427 ,-0.51904 ,-0.57042  ,-0.64854 ,-0.76736 ,-0.94491 ,-1.2017  ,-1.5558  ,-2.0135 ,-2.5586  ,-3.1485  ,-3.7239  ,-4.2311  ,-4.6389  ,-4.9428  ,-5.1556  ,-5.2975  ,-5.3883 ,-5.4441 ,-5.477  ,-5.4954 ,
                                                                    19.8  ,17.947,16.008 ,14.042 ,12.111 ,10.271   ,8.5697   ,7.0379  ,5.6926   ,4.5364   ,3.5611   ,2.7514    ,2.088    ,1.5503    ,1.1186   ,0.77445  ,0.50192  ,0.28727  ,0.11905   ,-0.012133,-0.11392   ,-0.19246 ,-0.2527  ,-0.29861 ,-0.33341 ,-0.35982 ,-0.38019 ,-0.39678 ,-0.41199 ,-0.42869 ,-0.45061 ,-0.483    ,-0.53333 ,-0.61221 ,-0.7342  ,-0.91784 ,-1.1838  ,-1.549  ,-2.0174  ,-2.5691  ,-3.158   ,-3.7245  ,-4.2171  ,-4.6088  ,-4.8979  ,-5.0991  ,-5.2329 ,-5.3185 ,-5.3715 ,-5.4033 ,
                                                                    20.461,18.652,16.739 ,14.776 ,12.826 ,10.946   ,9.189    ,7.5916  ,6.1765   ,4.951    ,3.9105   ,3.042     ,2.3274   ,1.7462    ,1.2783   ,0.90467  ,0.60843  ,0.37502  ,0.19219   ,0.049802 ,-0.060378  ,-0.14501 ,-0.20945 ,-0.25796 ,-0.29399 ,-0.32037 ,-0.33948 ,-0.35347 ,-0.36445 ,-0.37475 ,-0.38727 ,-0.40587  ,-0.43604 ,-0.48564 ,-0.56581 ,-0.69172 ,-0.88243 ,-1.1585 ,-1.5359  ,-2.0157  ,-2.5742  ,-3.1624  ,-3.7203  ,-4.1995  ,-4.5766  ,-4.8529  ,-5.0445 ,-5.1719 ,-5.254  ,-5.3056 ,
                                                                    21.101,19.343,17.463 ,15.512 ,13.551 ,11.639   ,9.8313   ,8.1716  ,6.6877   ,5.3923   ,4.2849   ,3.3553    ,2.5867   ,1.9593    ,1.4528   ,1.0475   ,0.72568  ,0.47201  ,0.27336   ,0.11884  ,-0.00045781,-0.091732,-0.16078 ,-0.21223 ,-0.24981 ,-0.27651 ,-0.29482 ,-0.30688 ,-0.31466 ,-0.3202  ,-0.32581 ,-0.33445  ,-0.35018 ,-0.37875 ,-0.42845 ,-0.51098 ,-0.64214 ,-0.84143,-1.1293  ,-1.5203  ,-2.0125  ,-2.5784  ,-3.1662  ,-3.7163  ,-4.1831  ,-4.5471  ,-4.8123 ,-4.9959 ,-5.1184 ,-5.1981 ,
                                                                    21.719,20.017,18.178 ,16.248 ,14.283 ,12.347   ,10.495   ,8.7771  ,7.2263   ,5.861    ,4.6855   ,3.6925    ,2.8673   ,2.1911    ,1.6434   ,1.2041   ,0.85487  ,0.57935  ,0.36359   ,0.19592  ,0.066726   ,-0.031787,-0.10591 ,-0.16069 ,-0.20014 ,-0.22752 ,-0.24546 ,-0.25617 ,-0.26162 ,-0.26366 ,-0.26429 ,-0.26586  ,-0.27142 ,-0.28524 ,-0.31339 ,-0.36459 ,-0.45111 ,-0.58934,-0.79918 ,-1.1008  ,-1.507   ,-2.0127  ,-2.5864  ,-3.1741  ,-3.7168  ,-4.172   ,-4.524  ,-4.7794 ,-4.9561 ,-5.0745 ,
                                                                    22.315,20.674,18.881 ,16.979 ,15.021 ,13.068   ,11.179   ,9.4073  ,7.7922   ,6.3578   ,5.1131   ,4.0549    ,3.1707   ,2.443     ,1.8515   ,1.376    ,0.99728  ,0.6982   ,0.46394   ,0.28201  ,0.14205    ,0.03561  ,-0.044147,-0.10271 ,-0.14446 ,-0.17294 ,-0.19099 ,-0.20096 ,-0.20488 ,-0.20455 ,-0.20178 ,-0.19855  ,-0.19725 ,-0.20108 ,-0.2145  ,-0.24388 ,-0.2984  ,-0.39087,-0.53823 ,-0.76065 ,-1.0778  ,-1.5006  ,-2.0205  ,-2.6022  ,-3.1896  ,-3.7248  ,-4.1688 ,-4.5095 ,-4.7555 ,-4.9259 ,
                                                                    22.886,21.31 ,19.57  ,17.704 ,15.76  ,13.799   ,11.88    ,10.061  ,8.3852   ,6.883    ,5.5689   ,4.4439    ,3.4985   ,2.7167    ,2.0789   ,1.5646   ,1.1542   ,0.82978  ,0.57551   ,0.37808  ,0.22634    ,0.11115  ,0.02507  ,-0.037872,-0.082483,-0.11262 ,-0.13139 ,-0.14133 ,-0.14456 ,-0.14295 ,-0.1382  ,-0.13206  ,-0.12652 ,-0.12401 ,-0.12784 ,-0.14266 ,-0.17516 ,-0.23491,-0.33525 ,-0.49361 ,-0.73037 ,-1.0644  ,-1.5047  ,-2.0389  ,-2.6282  ,-3.2145  ,-3.7417 ,-4.1743 ,-4.5038 ,-4.7408 ,
                                                                    23.433,21.924,20.242 ,18.419 ,16.499 ,14.538   ,12.598   ,10.737  ,9.0045   ,7.4368   ,6.0537   ,4.8608    ,3.8521   ,3.0136    ,2.3269   ,1.7714   ,1.3271   ,0.97529  ,0.69932   ,0.48498  ,0.32028    ,0.19533  ,0.10206  ,0.03395  ,-0.014262,-0.046798,-0.067059,-0.077805,-0.081311,-0.079505,-0.074093,-0.066697 ,-0.059007,-0.052975,-0.051081,-0.056713,-0.074691,-0.11199,-0.17859 ,-0.28833 ,-0.45912 ,-0.71141 ,-1.0633  ,-1.5211  ,-2.0691  ,-2.6648  ,-3.2489 ,-3.7671 ,-4.188  ,-4.5061 ,
                                                                    23.955,22.516,20.896 ,19.122 ,17.233 ,15.282   ,13.328   ,11.433  ,9.6496   ,8.0193   ,6.5681   ,5.3066    ,4.2329   ,3.3355    ,2.5972   ,1.9979   ,1.5172   ,1.1359   ,0.8363    ,0.60341  ,0.42433    ,0.2884   ,0.18684  ,0.11255  ,0.059778 ,0.023899 ,0.0011946,-0.01134 ,-0.016154,-0.015277,-0.010445,-0.0032129,0.0049111,0.012319 ,0.017099 ,0.016764 ,0.0078514,-0.01463,-0.057937,-0.13248 ,-0.2526  ,-0.43666 ,-0.70509 ,-1.0749  ,-1.5501  ,-2.1108  ,-2.7113 ,-3.2917 ,-3.8    ,-4.2085 ,
                                                                    24.453,23.084,21.531 ,19.812 ,17.961 ,16.028   ,14.07    ,12.148  ,10.319   ,8.6302   ,7.1126   ,5.7825    ,4.6424   ,3.6838    ,2.8912   ,2.2454   ,1.7258   ,1.3124   ,0.98715   ,0.73381  ,0.53869    ,0.39028  ,0.27908  ,0.19736  ,0.13884  ,0.098489 ,0.072237 ,0.056825 ,0.049616 ,0.048464 ,0.051582 ,0.05744   ,0.064646 ,0.071821 ,0.077451 ,0.079675 ,0.075996 ,0.062867,0.035111 ,-0.014844,-0.097856,-0.2288  ,-0.42647 ,-0.71122 ,-1.0989  ,-1.5905  ,-2.1627 ,-2.7663 ,-3.3414 ,-3.8386 ,
                                                                    24.926,23.629,22.144 ,20.486 ,18.681 ,16.774   ,14.82    ,12.879  ,11.012   ,9.2689   ,7.6873   ,6.2891    ,5.0815   ,4.0597    ,3.2103   ,2.5151   ,1.9537   ,1.5057   ,1.1523    ,0.87632  ,0.66319    ,0.50053  ,0.37809  ,0.28746  ,0.22184  ,0.17573  ,0.14473  ,0.12531  ,0.11464  ,0.11044  ,0.11086  ,0.11437   ,0.11963  ,0.12543  ,0.1305   ,0.13343  ,0.13238  ,0.12483 ,0.10712  ,0.073846 ,0.017087 ,-0.074505,-0.21636 ,-0.42768 ,-0.72864 ,-1.1336  ,-1.6408 ,-2.2228 ,-2.8279 ,-3.3961 ,
                                                                    25.375,24.15 ,22.737 ,21.142 ,19.389 ,17.516   ,15.575   ,13.625  ,11.726   ,9.9345   ,8.2922   ,6.827     ,5.5513   ,4.4643    ,3.5554   ,2.8078   ,2.2016   ,1.7162   ,1.3318    ,1.0307   ,0.79734    ,0.61841  ,0.48289  ,0.38171  ,0.30749  ,0.25429  ,0.21732  ,0.1928   ,0.17769  ,0.16956  ,0.16648  ,0.1669    ,0.1695   ,0.17314  ,0.17675  ,0.17918  ,0.17906  ,0.17456 ,0.1631   ,0.14087  ,0.1022   ,0.038778 ,-0.061271,-0.21395 ,-0.43883 ,-0.75577 ,-1.1776 ,-1.6993 ,-2.2895 ,-2.8943 ,
                                                                    25.802,24.649,23.307 ,21.78  ,20.085 ,18.254   ,16.334   ,14.382  ,12.46    ,10.626   ,8.9265   ,7.3959    ,6.0519   ,4.8982    ,3.9273   ,3.1242   ,2.47     ,1.9438   ,1.5255    ,1.1965   ,0.94038    ,0.74294  ,0.59237  ,0.47889  ,0.39454  ,0.33289  ,0.28881  ,0.25819  ,0.23779  ,0.22502  ,0.21786  ,0.21467   ,0.21414  ,0.21517  ,0.21676  ,0.21797  ,0.2177   ,0.21462 ,0.20687  ,0.19177  ,0.16531  ,0.12154  ,0.05171  ,-0.056605,-0.21998 ,-0.45831 ,-0.79102,-1.2291 ,-1.7642 ,-2.3612 ,
                                                                    26.207,25.124,23.855 ,22.398 ,20.765 ,18.983   ,17.092   ,15.148  ,13.211   ,11.34    ,9.5889   ,7.9954    ,6.5833   ,5.3615    ,4.326    ,3.4644   ,2.7586   ,2.1883   ,1.7329    ,1.373    ,1.0914     ,0.87311  ,0.70543  ,0.57787  ,0.48189  ,0.41054  ,0.35829  ,0.32072  ,0.29435  ,0.27641  ,0.26474  ,0.25763   ,0.25371  ,0.25188  ,0.25119  ,0.25076  ,0.24972  ,0.24702 ,0.24131  ,0.23068  ,0.21232  ,0.18201  ,0.13347  ,0.057483 ,-0.058939,-0.23294 ,-0.48467,-0.83295,-1.2867 ,-1.8341 ,
                                                                    26.592,25.578,24.382 ,22.997 ,21.43  ,19.701   ,17.848   ,15.92   ,13.976   ,12.076   ,10.278   ,8.6242    ,7.1448   ,5.8538    ,4.7515   ,3.8282   ,3.0674   ,2.4493   ,1.9533    ,1.5594   ,1.2496     ,1.008    ,0.82117  ,0.67782  ,0.56877  ,0.48656  ,0.42523  ,0.38001  ,0.34714  ,0.32365  ,0.30722  ,0.29602   ,0.2886   ,0.28382  ,0.28072  ,0.27847  ,0.2763   ,0.27336 ,0.26862  ,0.26069  ,0.2476   ,0.22637  ,0.19255  ,0.13952  ,0.057558 ,-0.066895,-0.25154,-0.51668,-0.88039,-1.3493 ,
                                                                    26.957,26.011,24.886 ,23.574 ,22.076 ,20.407   ,18.598   ,16.694  ,14.751   ,12.83    ,10.99    ,9.2803    ,7.7351   ,6.3743    ,5.203    ,4.2149   ,3.3956   ,2.7263   ,2.1862    ,1.7552   ,1.4144     ,1.1471   ,0.93904  ,0.77824  ,0.65481  ,0.56071  ,0.48951  ,0.43606  ,0.39628  ,0.36698  ,0.34563  ,0.33026   ,0.31932  ,0.31159  ,0.30607  ,0.30194  ,0.29847  ,0.29491 ,0.29046  ,0.28408  ,0.27434  ,0.25915  ,0.23538  ,0.19833  ,0.14099  ,0.053145 ,-0.07936,-0.27474,-0.55339,-0.93241,
                                                                    27.304,26.424,25.37  ,24.131 ,22.704 ,21.099   ,19.339   ,17.467  ,15.533   ,13.599   ,11.724   ,9.9613    ,8.3522   ,6.9216    ,5.6797   ,4.6242   ,3.743    ,3.0188   ,2.4313    ,1.9599   ,1.5854     ,1.29     ,1.0589   ,0.87907  ,0.74003  ,0.63311  ,0.55134  ,0.48917  ,0.44217  ,0.40685  ,0.38048  ,0.36092   ,0.34648  ,0.33584  ,0.32794  ,0.32193  ,0.31707  ,0.31271 ,0.30813  ,0.30254  ,0.29487  ,0.28364  ,0.26665  ,0.24059  ,0.20048  ,0.13892  ,0.045186,-0.09548,-0.30178,-0.59408,
                                                                    27.635,26.817,25.833 ,24.666 ,23.311 ,21.772   ,20.068   ,18.234  ,16.317   ,14.377   ,12.474   ,10.664    ,8.9938   ,7.4941    ,6.1806   ,5.0552   ,4.1092   ,3.3268   ,2.6885    ,2.1738   ,1.7628     ,1.4371   ,1.181    ,0.98064  ,0.82482  ,0.70419  ,0.61122  ,0.53987  ,0.48535  ,0.44385  ,0.41239  ,0.38862   ,0.3707   ,0.3572   ,0.34697  ,0.33911  ,0.33286  ,0.32756 ,0.3226   ,0.3173   ,0.31085  ,0.30217  ,0.2897   ,0.27112  ,0.24293  ,0.19986  ,0.13409 ,0.034383,-0.11463,-0.3321};

            return arrayCurvatureMin[y * arrayWidth + x];
        }
        case 8:
        {
            static const std::vector<double> arrayCurvatureMin   = {5.033 ,4.0224,3.1254,2.3318,1.6312,1.0136,0.46905,-0.01168,-0.43729 ,-0.81582,-1.1545  ,-1.46   ,-1.7376  ,-1.9921 ,-2.2271 ,-2.445   ,-2.6472   ,-2.8346 ,-3.0069 ,-3.164    ,-3.3055 ,-3.4314  ,-3.5418 ,-3.6376  ,-3.7197 ,-3.7896  ,-3.849    ,-3.8999  ,-3.9445  ,-3.9855  ,-4.0264   ,-4.0715   ,-4.1256  ,-4.1937  ,-4.2786   ,-4.3779  ,-4.4829  ,-4.5813   ,-4.6628  ,-4.7231  ,-4.7633  ,-4.7876  ,-4.8002  ,-4.8047   ,-4.8038    ,-4.7991  ,-4.7917  ,-4.7822  ,-4.7711   ,-4.7585  ,
                                                                    5.3617,4.3142,3.3825,2.5565,1.8255,1.1792,0.60766,0.10125 ,-0.34893 ,-0.75101,-1.1123  ,-1.439  ,-1.7363  ,-2.0084 ,-2.2581 ,-2.4874  ,-2.6971   ,-2.8879 ,-3.0598 ,-3.2131   ,-3.3481 ,-3.4657  ,-3.5669 ,-3.6531  ,-3.7261 ,-3.7876  ,-3.8396   ,-3.8841  ,-3.9238  ,-3.9613  ,-4.0004   ,-4.0454   ,-4.1014  ,-4.1728  ,-4.2612   ,-4.362   ,-4.465   ,-4.5579   ,-4.6318  ,-4.6844  ,-4.7181  ,-4.7371  ,-4.7458  ,-4.7473   ,-4.7439    ,-4.7371  ,-4.7278  ,-4.7166  ,-4.7038   ,-4.6896  ,
                                                                    5.7019,4.6161,3.648 ,2.7876,2.0242,1.3472,0.74646,0.21226 ,-0.26436 ,-0.69147,-1.0761  ,-1.424  ,-1.74    ,-2.0274 ,-2.2886 ,-2.5251  ,-2.7379   ,-2.9278 ,-3.0956 ,-3.2423   ,-3.3692 ,-3.4778  ,-3.57   ,-3.6476  ,-3.7127 ,-3.7672  ,-3.8133   ,-3.8531  ,-3.8892  ,-3.9246  ,-3.963    ,-4.0092   ,-4.068   ,-4.1433  ,-4.235    ,-4.3366  ,-4.4366  ,-4.523    ,-4.5892  ,-4.6344  ,-4.6619  ,-4.6762  ,-4.6811  ,-4.6797   ,-4.6739    ,-4.6649  ,-4.6537  ,-4.6406  ,-4.626    ,-4.6101  ,
                                                                    6.0535,4.9275,3.9212,3.0244,2.2263,1.5164,0.88438,0.32056 ,-0.18386 ,-0.63661,-1.0442  ,-1.412  ,-1.7439  ,-2.043  ,-2.3115 ,-2.5509  ,-2.7627   ,-2.9485 ,-3.1099 ,-3.2488   ,-3.3672 ,-3.4674  ,-3.5516 ,-3.6219  ,-3.6806 ,-3.7297  ,-3.7713   ,-3.8076  ,-3.8413  ,-3.8756  ,-3.9145   ,-3.9628   ,-4.0252  ,-4.1047  ,-4.1993   ,-4.3006  ,-4.3964  ,-4.4758   ,-4.5341  ,-4.5722  ,-4.594   ,-4.6037  ,-4.6051  ,-4.6009   ,-4.5926    ,-4.5815  ,-4.5683  ,-4.5532  ,-4.5367   ,-4.5189  ,
                                                                    6.4159,5.2479,4.2013,3.2658,2.4308,1.6859,1.0209 ,0.42625 ,-0.10631 ,-0.58402,-1.0128  ,-1.3974 ,-1.7416  ,-2.0481 ,-2.3195 ,-2.5581  ,-2.766    ,-2.9458 ,-3.1    ,-3.2311   ,-3.3419 ,-3.4349  ,-3.5125 ,-3.5771  ,-3.6308 ,-3.6758  ,-3.7141   ,-3.748   ,-3.7804  ,-3.8146  ,-3.8548   ,-3.9059   ,-3.9725  ,-4.0561  ,-4.1529   ,-4.2528  ,-4.3433  ,-4.4152   ,-4.4657  ,-4.497   ,-4.5133  ,-4.5188  ,-4.5168  ,-4.5097   ,-4.499     ,-4.4857  ,-4.4703  ,-4.4532  ,-4.4347   ,-4.4148  ,
                                                                    6.7884,5.5765,4.4873,3.5109,2.6371,1.8555,1.1565 ,0.53101 ,-0.028658,-0.52914,-0.9759  ,-1.3734 ,-1.7255  ,-2.0354 ,-2.3064 ,-2.5416  ,-2.7443   ,-2.9177 ,-3.065  ,-3.1894   ,-3.2938 ,-3.3811  ,-3.4537 ,-3.5139  ,-3.564  ,-3.606   ,-3.642    ,-3.6745  ,-3.7063  ,-3.7411  ,-3.7834   ,-3.8381   ,-3.9091  ,-3.9966  ,-4.0946   ,-4.1917  ,-4.276   ,-4.3401   ,-4.383   ,-4.4079  ,-4.419   ,-4.4204  ,-4.4152  ,-4.4053   ,-4.3921    ,-4.3765  ,-4.3589  ,-4.3396  ,-4.3189   ,-4.297   ,
                                                                    7.1702,5.9122,4.7783,3.7593,2.8452,2.0263,1.2935 ,0.63857 ,0.054311 ,-0.46548,-0.92615 ,-1.3324 ,-1.6886  ,-1.999  ,-2.2676 ,-2.4985  ,-2.6959   ,-2.8636 ,-3.0053 ,-3.1244   ,-3.224  ,-3.307   ,-3.376  ,-3.4331  ,-3.4806 ,-3.5205  ,-3.555    ,-3.5868  ,-3.6187  ,-3.6549  ,-3.6999   ,-3.7586   ,-3.8341  ,-3.9249  ,-4.023    ,-4.1162  ,-4.1935  ,-4.2496   ,-4.2851  ,-4.3038  ,-4.3101  ,-4.3076  ,-4.2992  ,-4.2865   ,-4.2708    ,-4.2528  ,-4.2329  ,-4.2114  ,-4.1885   ,-4.1643  ,
                                                                    7.5604,6.2544,5.0743,4.0114,3.0566,2.2011,1.4363 ,0.75476 ,0.14956  ,-0.3855 ,-0.85611 ,-1.2677 ,-1.6255  ,-1.9348 ,-2.2005 ,-2.4275  ,-2.6206   ,-2.7839 ,-2.9215 ,-3.0369   ,-3.1332 ,-3.2134  ,-3.2798 ,-3.3349  ,-3.3806 ,-3.4192  ,-3.4529   ,-3.4845  ,-3.5172  ,-3.5553  ,-3.6036   ,-3.6667   ,-3.7466  ,-3.8399  ,-3.9368   ,-4.0249  ,-4.0947  ,-4.1428   ,-4.1713  ,-4.1842  ,-4.1858  ,-4.1797  ,-4.1681  ,-4.1527   ,-4.1345    ,-4.1141  ,-4.0918  ,-4.0681  ,-4.0429   ,-4.0165  ,
                                                                    7.9585,6.6032,5.3761,4.2694,3.275 ,2.3849,1.5913 ,0.88689 ,0.26472  ,-0.28188,-0.75941 ,-1.1743 ,-1.5327  ,-1.8408 ,-2.1043 ,-2.3287  ,-2.5189   ,-2.6795 ,-2.8146 ,-2.9278   ,-3.0222 ,-3.1006  ,-3.1656 ,-3.2193  ,-3.264  ,-3.3019  ,-3.3353   ,-3.3672  ,-3.4012  ,-3.4418  ,-3.4939   ,-3.5617   ,-3.6457  ,-3.7405  ,-3.835    ,-3.917   ,-3.9789  ,-4.0192   ,-4.0409  ,-4.0484  ,-4.0457  ,-4.036   ,-4.0215  ,-4.0034   ,-3.9826    ,-3.9598  ,-3.9353  ,-3.9092  ,-3.8817   ,-3.8531  ,
                                                                    8.3647,6.9599,5.6862,4.5372,3.5058,2.5846,1.766  ,1.0426  ,0.40683  ,-0.14875,-0.63167 ,-1.0493 ,-1.4087  ,-1.7166 ,-1.9794 ,-2.2027  ,-2.3918   ,-2.5514 ,-2.6855 ,-2.7978   ,-2.8913 ,-2.969   ,-3.0332 ,-3.0864  ,-3.1306 ,-3.1682  ,-3.2018   ,-3.2346  ,-3.2704  ,-3.3141  ,-3.3705   ,-3.443    ,-3.5306  ,-3.6258  ,-3.7166   ,-3.7917  ,-3.8456  ,-3.8784   ,-3.8938  ,-3.8962  ,-3.8896  ,-3.8765  ,-3.8591  ,-3.8384   ,-3.8152    ,-3.7901  ,-3.7632  ,-3.7349  ,-3.7053   ,-3.6745  ,
                                                                    8.7807,7.3275,6.0094,4.8211,3.7562,2.8077,1.9678 ,1.2284  ,0.58116  ,0.017636,-0.4706  ,-0.89169,-1.2534  ,-1.5628 ,-1.8266 ,-2.0507  ,-2.2404   ,-2.4004 ,-2.5349 ,-2.6473   ,-2.741  ,-2.8187  ,-2.8829 ,-2.936   ,-2.9802 ,-3.0181  ,-3.0524   ,-3.0865  ,-3.1247  ,-3.1719  ,-3.233    ,-3.3102   ,-3.4007  ,-3.4951  ,-3.5812   ,-3.649   ,-3.6949  ,-3.7206   ,-3.7301  ,-3.728   ,-3.7176  ,-3.7015  ,-3.6813  ,-3.6582   ,-3.6327    ,-3.6054  ,-3.5764  ,-3.546   ,-3.5143   ,-3.4815  ,
                                                                    9.2101,7.7111,6.352 ,5.1282,4.0338,3.0614,2.2027 ,1.449   ,0.7908   ,0.21893 ,-0.27574 ,-0.70192,-1.0677  ,-1.3806 ,-1.6473 ,-1.8739  ,-2.0657   ,-2.2275 ,-2.3634 ,-2.477    ,-2.5716 ,-2.65    ,-2.7147 ,-2.7683  ,-2.813  ,-2.8515  ,-2.8869   ,-2.9229  ,-2.9641  ,-3.0155  ,-3.0816   ,-3.1633   ,-3.2558  ,-3.3484  ,-3.4288   ,-3.489   ,-3.5273  ,-3.5463   ,-3.5506  ,-3.5443  ,-3.5305  ,-3.5117  ,-3.4891  ,-3.4637   ,-3.4361    ,-3.4068  ,-3.3758  ,-3.3436  ,-3.3101   ,-3.2756  ,
                                                                    9.6583,8.1174,6.7212,5.4659,4.345 ,3.351 ,2.4747 ,1.7066  ,1.0367   ,0.45504 ,-0.047922,-0.48121,-0.85318 ,-1.1714 ,-1.4428 ,-1.6734  ,-1.8687   ,-2.0335 ,-2.1718 ,-2.2875   ,-2.3837 ,-2.4634  ,-2.5292 ,-2.5837  ,-2.6294 ,-2.6691  ,-2.7061   ,-2.7446  ,-2.7894  ,-2.8456  ,-2.9169   ,-3.0028   ,-3.0965  ,-3.1861  ,-3.2602   ,-3.3128  ,-3.3439  ,-3.3569   ,-3.3565  ,-3.3465  ,-3.3299  ,-3.3085  ,-3.2838  ,-3.2565   ,-3.2271    ,-3.196   ,-3.1635  ,-3.1296  ,-3.0946   ,-3.0587  ,
                                                                    10.132,8.5534,7.1241,5.84  ,4.6945,3.6795,2.7853 ,2.0017  ,1.3182   ,0.72468 ,0.21124  ,-0.23132,-0.61147 ,-0.93692,-1.2146 ,-1.4507  ,-1.6507   ,-1.8194 ,-1.9611 ,-2.0796   ,-2.1781 ,-2.2597  ,-2.3271 ,-2.383   ,-2.4302 ,-2.4716  ,-2.5109   ,-2.5525  ,-2.6016  ,-2.6632  ,-2.74     ,-2.8296   ,-2.9234  ,-3.0092  ,-3.0767   ,-3.1219  ,-3.1464  ,-3.1542   ,-3.1497  ,-3.1366  ,-3.1175  ,-3.0941  ,-3.0676  ,-3.0388   ,-3.0079    ,-2.9755  ,-2.9416  ,-2.9065  ,-2.8704   ,-2.8335  ,
                                                                    10.638,9.0253,7.5655,6.2541,5.0844,4.0476,3.1339 ,2.3328  ,1.6335   ,1.0258  ,0.49962  ,0.04569 ,-0.34456 ,-0.67891,-0.96437,-1.2072  ,-1.413    ,-1.5867 ,-1.7326 ,-1.8546   ,-1.9562 ,-2.0403  ,-2.1099 ,-2.1678  ,-2.2169 ,-2.2606  ,-2.3029   ,-2.3485  ,-2.4027  ,-2.4701  ,-2.5524   ,-2.6452   ,-2.7383  ,-2.8194  ,-2.8802   ,-2.9186  ,-2.9372  ,-2.9405   ,-2.9326  ,-2.917   ,-2.8959  ,-2.871   ,-2.8433  ,-2.8132   ,-2.7814    ,-2.7479  ,-2.7132  ,-2.6773  ,-2.6405   ,-2.603   ,
                                                                    11.181,9.5371,8.048 ,6.7094,5.5144,4.4541,3.5187 ,2.6976  ,1.9802   ,1.3559  ,0.81477  ,0.34748 ,-0.054637,-0.39945,-0.69407,-0.94489 ,-1.1576   ,-1.3372 ,-1.4881 ,-1.6144   ,-1.7196 ,-1.8068  ,-1.8792 ,-1.9397  ,-1.9916 ,-2.0383  ,-2.0843   ,-2.1347  ,-2.1948  ,-2.2686  ,-2.3564   ,-2.4517   ,-2.5431  ,-2.6192  ,-2.6735   ,-2.7057  ,-2.7193  ,-2.7188   ,-2.7082  ,-2.6907  ,-2.6683  ,-2.6424  ,-2.6137  ,-2.583    ,-2.5506    ,-2.5167  ,-2.4815  ,-2.4454  ,-2.4083   ,-2.3707  ,
                                                                    11.765,10.09 ,8.5718,7.2047,5.9824,4.8964,3.9368 ,3.0932  ,2.3551   ,1.712   ,1.1539   ,0.67141 ,0.25578  ,-0.10096,-0.40605,-0.66598 ,-0.88657  ,-1.073  ,-1.2298 ,-1.3612   ,-1.4707 ,-1.5618  ,-1.6377 ,-1.7015  ,-1.7567 ,-1.8073  ,-1.858    ,-1.914   ,-1.9807  ,-2.0615  ,-2.1544   ,-2.2515   ,-2.3406  ,-2.4114  ,-2.4596   ,-2.4864  ,-2.4958  ,-2.4924   ,-2.4799  ,-2.4612  ,-2.4379  ,-2.4115  ,-2.3825  ,-2.3516   ,-2.319     ,-2.2851  ,-2.2501  ,-2.214   ,-2.1773   ,-2.14    ,
                                                                    12.389,10.684,9.1348,7.7372,6.4853,5.371 ,4.3846 ,3.5162  ,2.7551   ,2.091   ,1.514    ,1.0145  ,0.58384  ,0.21378 ,-0.10302,-0.3732  ,-0.60273  ,-0.7969 ,-0.96048,-1.0977   ,-1.2124 ,-1.3081  ,-1.3881 ,-1.456   ,-1.5155 ,-1.5708  ,-1.6269   ,-1.6896  ,-1.7638  ,-1.8517  ,-1.9495   ,-2.0476   ,-2.1338  ,-2.1994  ,-2.2421   ,-2.2642  ,-2.2703  ,-2.2648   ,-2.2512  ,-2.2318  ,-2.2083  ,-2.1819  ,-2.1531  ,-2.1224   ,-2.0902    ,-2.0568  ,-2.0222  ,-1.9868  ,-1.9508   ,-1.9143  ,
                                                                    13.052,11.315,9.7332,8.303 ,7.019 ,5.8739,4.8583 ,3.9626  ,3.1764   ,2.4894  ,1.8917   ,1.3736  ,0.92637  ,0.54162 ,0.21188 ,-0.069677,-0.30916  ,-0.51206,-0.68326,-0.82718  ,-0.9478 ,-1.0488  ,-1.1339 ,-1.2067  ,-1.2712 ,-1.3321  ,-1.3948   ,-1.4651  ,-1.5475  ,-1.6425  ,-1.7447   ,-1.843    ,-1.9259  ,-1.9865  ,-2.0244   ,-2.0428  ,-2.0466  ,-2.0397   ,-2.0255  ,-2.0061  ,-1.9829  ,-1.9569  ,-1.9287  ,-1.8988   ,-1.8674    ,-1.8348  ,-1.8012  ,-1.7669  ,-1.732    ,-1.6968  ,
                                                                    13.749,11.979,10.363,8.8974,7.5791,6.4007,5.3538 ,4.4287  ,3.6153   ,2.9036  ,2.2834   ,1.7452  ,1.2799   ,0.87917 ,0.53526 ,0.24122  ,-0.0092556,-0.22182,-0.40156,-0.55303  ,-0.68041,-0.78761 ,-0.87846,-0.95693 ,-1.0274 ,-1.0949  ,-1.1652   ,-1.244   ,-1.335   ,-1.4372  ,-1.5431   ,-1.6409   ,-1.7202  ,-1.7762  ,-1.81     ,-1.8256  ,-1.8278  ,-1.8203   ,-1.8062  ,-1.7872  ,-1.7648  ,-1.7397  ,-1.7126  ,-1.6838   ,-1.6536    ,-1.6222  ,-1.59    ,-1.557   ,-1.5236   ,-1.49    ,
                                                                    14.477,12.671,11.018,9.5157,8.1608,6.9472,5.8667 ,4.9103  ,4.068    ,3.3297  ,2.6855   ,2.1256  ,1.6409   ,1.2229  ,0.86361 ,0.55599  ,0.2935    ,0.0703  ,-0.11885,-0.27873  ,-0.41368,-0.52785 ,-0.6253 ,-0.71032 ,-0.78766,-0.86275 ,-0.94159  ,-1.0296  ,-1.1296  ,-1.2386  ,-1.3474   ,-1.4441   ,-1.5197  ,-1.5717  ,-1.6022   ,-1.6158  ,-1.6171  ,-1.6096   ,-1.596   ,-1.578   ,-1.5568  ,-1.533   ,-1.5072  ,-1.4799   ,-1.4511    ,-1.4213  ,-1.3907  ,-1.3594  ,-1.3278   ,-1.296   ,
                                                                    15.228,13.386,11.694,10.153,8.7597,7.5088,6.393  ,5.4034  ,4.5304   ,3.764   ,3.0942   ,2.5112  ,2.0058   ,1.5692  ,1.1935  ,0.87118  ,0.59568   ,0.36094 ,0.16149 ,-0.0076056,-0.15093,-0.27284 ,-0.3777 ,-0.4701  ,-0.5552 ,-0.63882 ,-0.7271   ,-0.82508 ,-0.93411 ,-1.0493  ,-1.1602   ,-1.2552   ,-1.3273  ,-1.3756  ,-1.4036   ,-1.416   ,-1.417   ,-1.4101   ,-1.3975  ,-1.3808  ,-1.361   ,-1.3388  ,-1.3146  ,-1.2889   ,-1.2619    ,-1.2339  ,-1.2051  ,-1.1757  ,-1.146    ,-1.1162  ,
                                                                    15.998,14.118,12.386,10.805,9.3712,8.0815,6.9287 ,5.9043  ,4.999    ,4.2029  ,3.506    ,2.8986  ,2.3712   ,1.9149  ,1.5216  ,1.1836   ,0.89417   ,0.64699 ,0.43643 ,0.25733   ,0.10488 ,-0.025521,-0.13854,-0.23913 ,-0.33287,-0.42592 ,-0.52443  ,-0.63277 ,-0.75066 ,-0.87132 ,-0.98334  ,-1.0762   ,-1.1449  ,-1.1903  ,-1.2165   ,-1.2282  ,-1.2295  ,-1.2235   ,-1.2122  ,-1.197   ,-1.1789  ,-1.1585  ,-1.1362  ,-1.1123   ,-1.0872    ,-1.0611  ,-1.0342  ,-1.0068  ,-0.9792   ,-0.95157 ,
                                                                    16.781,14.862,13.089,11.466,9.9914,8.6614,7.4701 ,6.4095  ,5.4704   ,4.6433  ,3.918    ,3.2848  ,2.7342   ,2.2571  ,1.8452  ,1.4906   ,1.1863    ,0.92584 ,0.7034  ,0.51356   ,0.35129 ,0.2117   ,0.089803,-0.019762,-0.12296,-0.22626 ,-0.33566  ,-0.45456 ,-0.58085 ,-0.706   ,-0.81832  ,-0.90876  ,-0.97435 ,-1.0173  ,-1.0422   ,-1.0538  ,-1.0558  ,-1.051    ,-1.0413  ,-1.0279  ,-1.0116  ,-0.99306 ,-0.97266 ,-0.95073  ,-0.92756   ,-0.90341 ,-0.87856 ,-0.85326 ,-0.82775  ,-0.80229 ,
                                                                    17.573,15.614,13.8  ,12.134,10.617,9.2453,8.0142 ,6.916   ,5.9419   ,5.0824  ,4.3275   ,3.6674  ,3.0925   ,2.5935  ,2.1619  ,1.7898   ,1.4699    ,1.1954  ,0.96038 ,0.75915   ,0.58642 ,0.43699  ,0.30554 ,0.18628  ,0.072849,-0.041426,-0.1622   ,-0.29158 ,-0.42554 ,-0.55412 ,-0.66594  ,-0.75377  ,-0.81655 ,-0.85762 ,-0.88182  ,-0.89368 ,-0.8967  ,-0.89336  ,-0.88533 ,-0.87375 ,-0.85939 ,-0.84277 ,-0.8243  ,-0.80433  ,-0.78312   ,-0.76097 ,-0.73815 ,-0.71491 ,-0.69151  ,-0.6682  ,
                                                                    18.369,16.371,14.514,12.805,11.244,9.8304,8.5585 ,7.4215  ,6.4112   ,5.5181  ,4.7325   ,4.0444  ,3.4441   ,2.9223  ,2.4703  ,2.0798   ,1.7434    ,1.4543  ,1.206   ,0.99278   ,0.809   ,0.64915  ,0.50751 ,0.37788  ,0.25349 ,0.12762  ,-0.0048019,-0.14433 ,-0.28501 ,-0.41592 ,-0.52655  ,-0.61169  ,-0.67201 ,-0.71167 ,-0.73559  ,-0.74806 ,-0.75229 ,-0.75053  ,-0.74428 ,-0.73458 ,-0.72212 ,-0.70743 ,-0.69088 ,-0.67283  ,-0.65356   ,-0.63336 ,-0.61252 ,-0.5913  ,-0.56995  ,-0.54872 ,
                                                                    19.166,17.128,15.23 ,13.477,11.872,10.415,9.1011 ,7.9244  ,6.8768   ,5.9491  ,5.1317   ,4.4145  ,3.788    ,3.2425  ,2.7691  ,2.3595   ,2.0061    ,1.7016  ,1.4395  ,1.2138    ,1.0184  ,0.84755  ,0.69513 ,0.55447  ,0.41849 ,0.28052  ,0.1364    ,-0.012681,-0.15896 ,-0.29111 ,-0.40002  ,-0.48251  ,-0.54074 ,-0.57941 ,-0.6034   ,-0.61672 ,-0.62231 ,-0.62219  ,-0.61774 ,-0.60989 ,-0.5993  ,-0.58648 ,-0.57179 ,-0.5556   ,-0.53819   ,-0.51987 ,-0.50093 ,-0.48164 ,-0.46226  ,-0.44302 ,
                                                                    19.962,17.885,15.944,14.147,12.498,10.997,9.6409 ,8.4236  ,7.3378   ,6.3745  ,5.5243   ,4.7773  ,4.1235   ,3.5535  ,3.0581  ,2.6287   ,2.2575    ,1.9372  ,1.6607  ,1.4219    ,1.2145  ,1.0321   ,0.86834 ,0.71603  ,0.56786 ,0.41744  ,0.26179   ,0.10402  ,-0.046643,-0.17905 ,-0.28584  ,-0.36583  ,-0.42234 ,-0.4604  ,-0.48474  ,-0.49908 ,-0.50611 ,-0.50765  ,-0.50496 ,-0.49891 ,-0.49013 ,-0.47908 ,-0.46617 ,-0.45174  ,-0.4361    ,-0.41956 ,-0.40243 ,-0.38497 ,-0.36745  ,-0.35008 ,
                                                                    20.754,18.639,16.656,14.816,13.122,11.577,10.177 ,8.9187  ,7.7938   ,6.7942  ,5.9104   ,5.1325  ,4.4508   ,3.8555  ,3.3373  ,2.8875   ,2.498     ,2.1612  ,1.8699  ,1.6176    ,1.3975  ,1.2032   ,1.0275  ,0.86293  ,0.70206 ,0.539    ,0.37224   ,0.2068   ,0.052993 ,-0.078824,-0.18328  ,-0.26097  ,-0.31613 ,-0.3539  ,-0.37881  ,-0.3943  ,-0.40279 ,-0.40598  ,-0.405   ,-0.40067 ,-0.39359 ,-0.38423 ,-0.37298 ,-0.3602   ,-0.34622   ,-0.33135 ,-0.31591 ,-0.30017 ,-0.28438  ,-0.26877 ,
                                                                    21.542,19.389,17.365,15.481,13.743,12.153,10.71  ,9.4099  ,8.2453   ,7.2084  ,6.2902   ,5.4808  ,4.7703   ,4.1489  ,3.6073  ,3.1365   ,2.7281    ,2.3743  ,2.0678  ,1.8014    ,1.5683  ,1.3614   ,1.1732  ,0.99586  ,0.82185 ,0.64614  ,0.46891   ,0.29699  ,0.14122  ,0.010662 ,-0.091383 ,-0.16706  ,-0.22121 ,-0.25898 ,-0.28464  ,-0.30134 ,-0.31132 ,-0.3161   ,-0.31676 ,-0.31405 ,-0.30857 ,-0.30078 ,-0.29108 ,-0.27984  ,-0.26741   ,-0.25411 ,-0.24024 ,-0.2261  ,-0.21192  ,-0.19792 ,
                                                                    22.324,20.136,18.072,16.145,14.362,12.728,11.241 ,9.8976  ,8.6926   ,7.6179  ,6.6645   ,5.8228  ,5.0828   ,4.4347  ,3.869   ,3.3766   ,2.9488    ,2.5775  ,2.2552  ,1.9744    ,1.7277  ,1.5078   ,1.3065  ,1.1157   ,0.92827 ,0.74008  ,0.55322   ,0.37607  ,0.21941  ,0.090585 ,-0.0091119,-0.083089 ,-0.13658 ,-0.1746  ,-0.20114  ,-0.21911 ,-0.23055 ,-0.23688  ,-0.23908 ,-0.2379  ,-0.23391 ,-0.22758 ,-0.21932 ,-0.20952  ,-0.19852   ,-0.18666 ,-0.17426 ,-0.1616  ,-0.14891  ,-0.13641 ,
                                                                    23.102,20.878,18.775,16.805,14.979,13.3  ,11.769 ,10.383  ,9.1369   ,8.0237  ,7.0344   ,6.1596  ,5.3895   ,4.714   ,4.1236  ,3.6089   ,3.1612    ,2.7719  ,2.4333  ,2.1375    ,1.8768  ,1.6432   ,1.4283  ,1.2237   ,1.0225  ,0.8222   ,0.62673   ,0.44561  ,0.28895  ,0.16214  ,0.064613  ,-0.0080142,-0.061157,-0.099641,-0.12719  ,-0.14646 ,-0.15934 ,-0.16715  ,-0.17082 ,-0.17107 ,-0.16847 ,-0.16349 ,-0.15657 ,-0.14809  ,-0.13842   ,-0.12791 ,-0.11686 ,-0.10557 ,-0.094263 ,-0.08314 ,
                                                                    23.874,21.617,19.475,17.465,15.595,13.871,12.296 ,10.867  ,9.5793   ,8.4269  ,7.4011   ,6.4926  ,5.6916   ,4.9881  ,4.3723  ,3.8348   ,3.3665    ,2.9588  ,2.6034  ,2.2921    ,2.0168  ,1.769    ,1.5398  ,1.3208   ,1.1058  ,0.89399  ,0.69104   ,0.50714  ,0.35119  ,0.22649  ,0.13086   ,0.05922   ,0.0061209,-0.033025,-0.061666 ,-0.082267,-0.096561,-0.10578  ,-0.11084 ,-0.11243 ,-0.11113 ,-0.10741 ,-0.10173 ,-0.09448  ,-0.086048  ,-0.076782,-0.066998,-0.056973,-0.04694  ,-0.037084,
                                                                    24.641,22.352,20.174,18.123,16.21 ,14.442,12.823 ,11.35   ,10.021   ,8.8289  ,7.766    ,6.823   ,5.9905   ,5.2583  ,4.6165  ,4.0556   ,3.5662    ,3.1394  ,2.7667  ,2.4394    ,2.1489  ,1.8861   ,1.642   ,1.4082   ,1.1796  ,0.957    ,0.74777   ,0.56215  ,0.40739  ,0.28474  ,0.19066   ,0.11963   ,0.066293 ,0.026311 ,-0.0035112,-0.025454,-0.041134,-0.051722 ,-0.058097,-0.060951,-0.060858,-0.058321,-0.053794,-0.047697 ,-0.040419  ,-0.032316,-0.023704,-0.014857,-0.0060002,0.0026897,
                                                                    25.404,23.086,20.872,18.782,16.827,15.015,13.351 ,11.835  ,10.463   ,9.2311  ,8.1304   ,7.1524  ,6.2876   ,5.526   ,4.8576  ,4.2726   ,3.7615    ,3.3151  ,2.9244  ,2.5805    ,2.2741  ,1.9958   ,1.7362  ,1.4873   ,1.2452  ,1.0128   ,0.79845   ,0.61199  ,0.45869  ,0.33789  ,0.24497   ,0.17419   ,0.12034  ,0.079359 ,0.048279  ,0.024978 ,0.0079334,-0.0039758,-0.011609,-0.01566 ,-0.016712,-0.015283,-0.011844,-0.0068289,-0.00063476,0.0063782,0.013894 ,0.021643 ,0.029408  ,0.037022 ,
                                                                    26.163,23.817,21.57 ,19.441,17.445,15.59 ,13.882 ,12.322  ,10.908   ,9.6348  ,8.4957   ,7.482   ,6.5842   ,5.7925  ,5.0968  ,4.4871   ,3.9537    ,3.487   ,3.0778  ,2.7165    ,2.3936  ,2.0991   ,1.8233  ,1.5591   ,1.3041  ,1.0627   ,0.84451   ,0.65789  ,0.50612  ,0.38684  ,0.29467   ,0.22377   ,0.16916  ,0.12703  ,0.094622  ,0.069942 ,0.051548 ,0.038353  ,0.029507 ,0.024313 ,0.02217  ,0.022545 ,0.02495  ,0.02894   ,0.034108   ,0.040092 ,0.046576 ,0.053298 ,0.060047  ,0.066664 ,
                                                                    26.919,24.547,22.269,20.103,18.066,16.168,14.416 ,12.812  ,11.356   ,10.041  ,8.8633   ,7.8132  ,6.8818   ,6.0593  ,5.3355  ,4.7004   ,4.144     ,3.6564  ,3.2279  ,2.8487    ,2.5084  ,2.1969   ,1.9045  ,1.6249   ,1.3575  ,1.1084   ,0.88724   ,0.70092  ,0.55058  ,0.43243  ,0.34056   ,0.2692    ,0.21358  ,0.17016  ,0.13635   ,0.11027  ,0.090527 ,0.076069  ,0.066044 ,0.059743 ,0.056551 ,0.055913 ,0.057327 ,0.060335  ,0.064523   ,0.069527 ,0.075035 ,0.08079  ,0.086588  ,0.092278 ,
                                                                    27.672,25.277,22.968,20.767,18.69 ,16.75 ,14.955 ,13.308  ,11.808   ,10.452  ,9.2343   ,8.1472  ,7.1815   ,6.3274  ,5.5749  ,4.9137   ,4.3335    ,3.8243  ,3.3759  ,2.9778    ,2.6195  ,2.2901   ,1.9806  ,1.6857   ,1.4067  ,1.1511   ,0.92779   ,0.742    ,0.59284  ,0.47535  ,0.38336   ,0.31118   ,0.25434  ,0.20949  ,0.1742    ,0.14669  ,0.1256   ,0.10989   ,0.098703 ,0.091321 ,0.087104 ,0.085482 ,0.085934 ,0.087992  ,0.091236   ,0.095301 ,0.099878 ,0.10472  ,0.10962   ,0.11444  ,
                                                                    28.423,26.007,23.67 ,21.434,19.319,17.338,15.5   ,13.808  ,12.265   ,10.868  ,9.61     ,8.4852  ,7.4845   ,6.5981  ,5.8161  ,5.128    ,4.5233    ,3.9917  ,3.5226  ,3.1049    ,2.7275  ,2.3797   ,2.0528  ,1.7428   ,1.4531  ,1.1919   ,0.96715   ,0.78192  ,0.63355  ,0.51625  ,0.42368   ,0.35038   ,0.29209  ,0.24568  ,0.20884   ,0.17985  ,0.1574   ,0.14044   ,0.1281   ,0.11965  ,0.11442  ,0.11183  ,0.11134  ,0.11247   ,0.11479    ,0.11795  ,0.12163  ,0.1256   ,0.12965   ,0.13367  ,
                                                                    29.172,26.736,24.373,22.106,19.954,17.932,16.051 ,14.316  ,12.73    ,11.29   ,9.9914   ,8.8282  ,7.7917   ,6.8724  ,6.0601  ,5.3443   ,4.7144    ,4.1596  ,3.6688  ,3.2306    ,2.8334  ,2.4664   ,2.1217  ,1.7973   ,1.4978  ,1.2319   ,1.0062    ,0.82133  ,0.6733   ,0.55566  ,0.4621    ,0.38736   ,0.32742  ,0.27932  ,0.24084   ,0.21033  ,0.18649  ,0.16826   ,0.15477  ,0.14525  ,0.13901  ,0.13546  ,0.13403  ,0.13424   ,0.13567    ,0.13794  ,0.14076  ,0.14388  ,0.14714   ,0.15038  ,
                                                                    29.919,27.467,25.079,22.781,20.594,18.532,16.609 ,14.831  ,13.201   ,11.719  ,10.379   ,9.1773  ,8.1043   ,7.1512  ,6.3078  ,5.5636   ,4.9076    ,4.3287  ,3.8153  ,3.3557    ,2.9377  ,2.5509   ,2.1885  ,1.8501   ,1.5419  ,1.2722   ,1.0455    ,0.86078  ,0.71255  ,0.59407  ,0.4991    ,0.42262   ,0.36083  ,0.3109   ,0.27071   ,0.23863  ,0.21336  ,0.19385   ,0.17918  ,0.16857  ,0.16132  ,0.1568   ,0.15443  ,0.15374   ,0.15427    ,0.15567  ,0.15765  ,0.15997  ,0.16245   ,0.16496  ,
                                                                    30.665,28.198,25.788,23.462,21.24 ,19.14 ,17.175 ,15.354  ,13.681   ,12.155  ,10.775   ,9.5332  ,8.4231   ,7.4354  ,6.5601  ,5.7865   ,5.1036    ,4.4997  ,3.9628  ,3.4807    ,3.0411  ,2.634    ,2.2539  ,1.9023   ,1.5864  ,1.3133   ,1.0858    ,0.90071  ,0.75172  ,0.63188  ,0.5351    ,0.45661   ,0.39277  ,0.34088  ,0.29889   ,0.26517  ,0.23843  ,0.21759   ,0.20173  ,0.19002  ,0.18173  ,0.17622  ,0.17292  ,0.1713    ,0.17095    ,0.17149  ,0.17265  ,0.17418  ,0.17592   ,0.17773  ,
                                                                    31.408,28.93 ,26.5  ,24.147,21.892,19.755,17.75  ,15.886  ,14.169   ,12.601  ,11.178   ,9.8967  ,8.7487   ,7.7257  ,6.8177  ,6.0139   ,5.3031    ,4.6732  ,4.1118  ,3.6062    ,3.1442  ,2.7163   ,2.3189  ,1.955    ,1.6321  ,1.3561   ,1.1275    ,0.9415   ,0.79116  ,0.66945  ,0.57049   ,0.4897    ,0.42363  ,0.36965  ,0.32575   ,0.29033  ,0.26207  ,0.23986   ,0.22276  ,0.20991  ,0.20056  ,0.19405  ,0.18978  ,0.18725   ,0.186      ,0.1857   ,0.18604  ,0.18681  ,0.18782   ,0.18897  ,
                                                                    32.15 ,29.662,27.215,24.837,22.551,20.378,18.333 ,16.427  ,14.667   ,13.056  ,11.591   ,10.269  ,9.082    ,8.0229  ,7.0813  ,6.2464   ,5.5067    ,4.8498  ,4.2629  ,3.7328    ,3.2475  ,2.7987   ,2.3843  ,2.0089   ,1.6798  ,1.401    ,1.1709    ,0.98346  ,0.83116  ,0.70709  ,0.60558   ,0.52224   ,0.45374  ,0.39754  ,0.35163   ,0.31443  ,0.28458  ,0.26096   ,0.24257  ,0.22854  ,0.2181   ,0.21056  ,0.20531  ,0.20183   ,0.19969    ,0.19854  ,0.19808  ,0.19809  ,0.19841   ,0.19892  ,
                                                                    32.89 ,30.395,27.933,25.531,23.216,21.008,18.925 ,16.978  ,15.175   ,13.52   ,12.013   ,10.649  ,9.4236   ,8.3276  ,7.3515  ,6.4846   ,5.715     ,5.03    ,4.4164  ,3.8608    ,3.3516  ,2.8817   ,2.4509  ,2.065    ,1.7302  ,1.4485   ,1.2164    ,1.0268   ,0.87198  ,0.74509  ,0.64067   ,0.55451   ,0.4834   ,0.42483  ,0.37682   ,0.33775  ,0.30624  ,0.28114   ,0.26141  ,0.24616  ,0.23458  ,0.22597  ,0.21971  ,0.21528   ,0.21224    ,0.21023  ,0.20898  ,0.20825  ,0.20789   ,0.20777  ,
                                                                    33.628,31.128,28.653,26.231,23.888,21.647,19.526 ,17.538  ,15.693   ,13.995  ,12.445   ,11.039  ,9.774    ,8.6403  ,7.6289  ,6.729    ,5.9284    ,5.2143  ,4.5728  ,3.9907    ,3.4569  ,2.9661   ,2.5198  ,2.1239   ,1.7837  ,1.499    ,1.2642    ,1.0718   ,0.91385  ,0.78367  ,0.676     ,0.58679   ,0.51287  ,0.4518   ,0.40156   ,0.36052  ,0.32728  ,0.30062   ,0.2795   ,0.26297  ,0.25021  ,0.24049  ,0.23319  ,0.22779   ,0.22383    ,0.22097  ,0.21892  ,0.21746  ,0.21642   ,0.2157   ,
                                                                    34.363,31.861,29.375,26.935,24.567,22.294,20.136 ,18.108  ,16.221   ,14.48   ,12.887   ,11.44   ,10.134   ,8.9615  ,7.9139  ,6.98     ,6.1474    ,5.4029  ,4.7324  ,4.1228    ,3.5642  ,3.0528   ,2.5915  ,2.1865   ,1.8409  ,1.5525   ,1.3146    ,1.1186   ,0.95696  ,0.82305  ,0.71181   ,0.61928   ,0.54238  ,0.47864  ,0.42606   ,0.38296  ,0.34789  ,0.31961   ,0.29702  ,0.27915  ,0.26515  ,0.25428  ,0.24591  ,0.23951   ,0.23462    ,0.2309   ,0.22805  ,0.22586  ,0.22417   ,0.22285  ,
                                                                    35.095,32.593,30.1  ,27.644,25.252,22.949,20.756 ,18.689  ,16.76    ,14.976  ,13.34    ,11.85   ,10.503   ,9.2917  ,8.207   ,7.238    ,6.3723    ,5.5962  ,4.8956  ,4.2578    ,3.674   ,3.1423   ,2.6669  ,2.2531   ,1.9021  ,1.6095   ,1.3675    ,1.1674   ,1.0015   ,0.86342  ,0.74828   ,0.6522    ,0.57212  ,0.50557  ,0.45052   ,0.40524  ,0.36825  ,0.33825   ,0.31412  ,0.29484  ,0.27956  ,0.2675   ,0.25802  ,0.25059   ,0.24475    ,0.24016  ,0.23651  ,0.2336   ,0.23126   ,0.22935  ,
                                                                    35.823,33.325,30.826,28.356,25.943,23.612,21.385 ,19.28   ,17.31    ,15.483  ,13.803   ,12.271  ,10.882   ,9.6312  ,8.5085  ,7.5035   ,6.6034    ,5.7945  ,5.0627  ,4.3959    ,3.7869  ,3.2356   ,2.7467  ,2.3243   ,1.9675  ,1.6701   ,1.4233    ,1.2183   ,1.0476   ,0.90496  ,0.7856    ,0.68572   ,0.60226  ,0.53274  ,0.47508   ,0.42751  ,0.3885   ,0.3567    ,0.33094  ,0.31019  ,0.29355  ,0.28025  ,0.26963  ,0.26115   ,0.25434    ,0.24887  ,0.24442  ,0.24079  ,0.2378    ,0.23531  ,
                                                                    36.548,34.055,31.554,29.072,26.639,24.282,22.023 ,19.881  ,17.87    ,16.001  ,14.279   ,12.703  ,11.272   ,9.9805  ,8.8188  ,7.7766   ,6.8411    ,5.9982  ,5.2341  ,4.5377    ,3.9036  ,3.3333   ,2.8314  ,2.4005   ,2.0374  ,1.7343   ,1.4819    ,1.2715   ,1.0954   ,0.94781  ,0.82393   ,0.71999   ,0.63296  ,0.56031  ,0.4999    ,0.44991  ,0.40877  ,0.37507   ,0.3476   ,0.3253   ,0.30725  ,0.29266  ,0.28085  ,0.27128   ,0.26349    ,0.25712  ,0.25187  ,0.24752  ,0.24388   ,0.24082};
        }
        case 9:
        {
            static const std::vector<double> arrayCurvatureMin   = {9.3553,6.9589,4.9189,3.227 ,1.867 ,0.81473,0.039801,-0.49344,-0.82387 ,-0.99168,-1.0364  ,-0.99498,-0.90061,-0.78149,-0.66049 ,-0.55503  ,-0.47743 ,-0.43549 ,-0.43323 ,-0.47168,-0.54961,-0.66429 ,-0.81199 ,-0.98846 ,-1.1892  ,-1.4098 ,-1.6459 ,-1.893    ,-2.1472  ,-2.4041  ,-2.6594  ,-2.9089  ,-3.1481  ,-3.3727  ,-3.5783  ,-3.7609  ,-3.9169  ,-4.0433  ,-4.1376  ,-4.198   ,-4.2231  ,-4.2124   ,-4.1653  ,-4.0821  ,-3.9628  ,-3.808    ,-3.6181  ,-3.3938  ,-3.1358  ,-2.8447    ,
                                                                    10.224,7.7644,5.6607,3.9056,2.4839,1.3726 ,0.54214 ,-0.04254,-0.41992 ,-0.63007,-0.7126  ,-0.70486,-0.64047,-0.54829,-0.45187 ,-0.36933  ,-0.31369 ,-0.29338 ,-0.31296 ,-0.3739 ,-0.47532,-0.61468 ,-0.78827 ,-0.99168 ,-1.2201  ,-1.4684 ,-1.7313 ,-2.0036   ,-2.2799  ,-2.5547  ,-2.8225  ,-3.0783  ,-3.317   ,-3.534   ,-3.7255  ,-3.8881  ,-4.0192  ,-4.1169  ,-4.18    ,-4.2076  ,-4.1994  ,-4.1554   ,-4.0759  ,-3.9611  ,-3.8115  ,-3.6277   ,-3.4101  ,-3.1595  ,-2.8765  ,-2.5618    ,
                                                                    11.103,8.5806,6.4128,4.594 ,3.1099,1.9386 ,1.0513  ,0.41358 ,-0.012589,-0.26717,-0.3899  ,-0.41848,-0.38704,-0.32508,-0.25686 ,-0.20124  ,-0.17189 ,-0.17783 ,-0.2241  ,-0.31247,-0.44218,-0.61055 ,-0.81352 ,-1.0461  ,-1.3024  ,-1.5765 ,-1.8619 ,-2.152    ,-2.4404  ,-2.7209  ,-2.9877  ,-3.2355  ,-3.4599  ,-3.6571  ,-3.8242  ,-3.9592  ,-4.0605  ,-4.1273  ,-4.159   ,-4.1556  ,-4.1172  ,-4.044    ,-3.9363  ,-3.7946  ,-3.6193  ,-3.4109   ,-3.1699  ,-2.8969  ,-2.5927  ,-2.258     ,
                                                                    11.992,9.406 ,7.1739,5.2906,3.7431,2.5104 ,1.5645  ,0.87177 ,0.39439  ,0.092581,-0.073502,-0.14191,-0.14731,-0.11984,-0.084432,-0.060598 ,-0.062588,-0.099837,-0.17759 ,-0.29758,-0.45874,-0.6578  ,-0.88983 ,-1.1488  ,-1.4277  ,-1.7194 ,-2.0164 ,-2.3117   ,-2.5983  ,-2.8705  ,-3.1228  ,-3.3511  ,-3.5519  ,-3.7227  ,-3.8616  ,-3.9676  ,-4.0398  ,-4.0779  ,-4.0818  ,-4.0516  ,-3.9874  ,-3.8896   ,-3.7583  ,-3.5939  ,-3.3968  ,-3.1675   ,-2.9065  ,-2.6144  ,-2.292   ,-1.9402    ,
                                                                    12.888,10.239,7.9417,5.993 ,4.3806,3.0846 ,2.0779  ,1.3273  ,0.79546  ,0.44278 ,0.22925  ,0.1165  ,0.069408,0.057279,0.054644 ,0.041545  ,0.0034537,-0.069109,-0.18115 ,-0.33389,-0.52547,-0.75163 ,-1.0064  ,-1.2826  ,-1.5727  ,-1.8688 ,-2.1635 ,-2.45     ,-2.7223  ,-2.9753  ,-3.205   ,-3.4082  ,-3.5826  ,-3.7264  ,-3.8386  ,-3.9183  ,-3.9652  ,-3.979   ,-3.9596  ,-3.907   ,-3.8214  ,-3.7029   ,-3.5517  ,-3.368   ,-3.1524  ,-2.9051   ,-2.6268  ,-2.3183  ,-1.9803  ,-1.6142    ,
                                                                    13.79 ,11.077,8.7132,6.6974,5.018 ,3.6562 ,2.5854  ,1.7733  ,1.1829   ,0.7747  ,0.5087   ,0.34634 ,0.25219 ,0.19527 ,0.14989  ,0.096069  ,0.019467 ,-0.088947,-0.23349 ,-0.41453,-0.62928,-0.87265 ,-1.1381  ,-1.4181  ,-1.7053  ,-1.9924 ,-2.2728 ,-2.5407   ,-2.7914  ,-3.0208  ,-3.226   ,-3.4044  ,-3.5543  ,-3.6745  ,-3.7639  ,-3.8218  ,-3.8479  ,-3.8418  ,-3.8032  ,-3.7322  ,-3.6287  ,-3.4927   ,-3.3245  ,-3.1243  ,-2.8925  ,-2.6297   ,-2.3365  ,-2.0138  ,-1.6628  ,-1.2848    ,
                                                                    14.694,11.915,9.4836,7.3983,5.6491,4.2178 ,3.079   ,2.2007  ,1.5467   ,1.0777  ,0.75382  ,0.53669 ,0.39093 ,0.28571 ,0.19564  ,0.1012    ,-0.011348,-0.15029 ,-0.31924 ,-0.51794,-0.74328,-0.99017 ,-1.2524  ,-1.5234  ,-1.7965  ,-2.0658 ,-2.3257 ,-2.5717   ,-2.7997  ,-3.0066  ,-3.1897  ,-3.347   ,-3.4768  ,-3.5778  ,-3.6489  ,-3.6895  ,-3.6988  ,-3.6765  ,-3.6221  ,-3.5356  ,-3.4169  ,-3.2659   ,-3.083   ,-2.8683  ,-2.6224  ,-2.3461   ,-2.04    ,-1.7054  ,-1.3436  ,-0.95625   ,
                                                                    15.594,12.748,10.246,8.0882,6.2652,4.7601 ,3.5482  ,2.5987  ,1.8758   ,1.3409  ,0.95498  ,0.67998 ,0.48121 ,0.3285  ,0.19712  ,0.068224  ,-0.071354,-0.22956 ,-0.40995 ,-0.61258,-0.83502,-1.0731  ,-1.3219  ,-1.576   ,-1.83    ,-2.079  ,-2.3183 ,-2.544    ,-2.7525  ,-2.9408  ,-3.1064  ,-3.2471  ,-3.3613  ,-3.4474  ,-3.5042  ,-3.5309  ,-3.5266  ,-3.4909  ,-3.4233  ,-3.3236  ,-3.1916  ,-3.0276   ,-2.8317  ,-2.6044  ,-2.3463  ,-2.0582   ,-1.7414  ,-1.3971  ,-1.0269  ,-0.63275   ,
                                                                    16.484,13.568,10.992,8.757 ,6.8557,5.2718 ,3.982   ,2.9566  ,2.161    ,1.5579  ,1.1092   ,0.7779  ,0.53037 ,0.33719 ,0.17418  ,0.022638  ,-0.13077 ,-0.29462 ,-0.47341 ,-0.66845,-0.87866,-1.1013  ,-1.3326  ,-1.5683  ,-1.8039  ,-2.0349 ,-2.2572 ,-2.4668   ,-2.6603  ,-2.8347  ,-2.9873  ,-3.1157  ,-3.218   ,-3.2926  ,-3.3382  ,-3.3537  ,-3.3382  ,-3.2911  ,-3.2121  ,-3.1008  ,-2.9573  ,-2.7817   ,-2.5745  ,-2.3362  ,-2.0676  ,-1.7698   ,-1.4443  ,-1.0925  ,-0.71642 ,-0.31825   ,
                                                                    17.355,14.364,11.71 ,9.3936,7.4094,5.7429 ,4.3721  ,3.2691  ,2.4011   ,1.7322  ,1.2259   ,0.84628 ,0.56039 ,0.33902 ,0.15766  ,-0.0031405,-0.15795 ,-0.31678 ,-0.48571 ,-0.66755,-0.86258,-1.0692  ,-1.2845  ,-1.5047  ,-1.7257  ,-1.9433 ,-2.1532 ,-2.3515   ,-2.5345  ,-2.6991  ,-2.8423  ,-2.9616  ,-3.055   ,-3.1207  ,-3.1572  ,-3.1634  ,-3.1384  ,-3.0816  ,-2.9925  ,-2.8711  ,-2.7174  ,-2.5318   ,-2.3147  ,-2.067   ,-1.7897  ,-1.4842   ,-1.1521  ,-0.79528 ,-0.41598 ,-0.016673  ,
                                                                    18.194,15.125,12.389,9.9882,7.919 ,6.169  ,4.7185  ,3.5417  ,2.6075   ,1.8819  ,1.3291   ,0.91383 ,0.60292 ,0.36636 ,0.1784   ,0.018015  ,-0.13102 ,-0.28034 ,-0.43748 ,-0.60644,-0.78841,-0.98241 ,-1.1859  ,-1.3954  ,-1.6067  ,-1.8154 ,-2.0174 ,-2.2084   ,-2.3845  ,-2.5425  ,-2.6791  ,-2.7918  ,-2.8784  ,-2.9369  ,-2.9659  ,-2.9643  ,-2.9311  ,-2.8657  ,-2.7679  ,-2.6376  ,-2.475   ,-2.2806   ,-2.0552  ,-1.7998  ,-1.5158  ,-1.2046   ,-0.8683  ,-0.50904 ,-0.1293  ,0.26812    ,
                                                                    18.992,15.843,13.023,10.538,8.3862,6.5578 ,5.0353  ,3.7946  ,2.8062   ,2.0368  ,1.4509   ,1.0126  ,0.68765 ,0.44458 ,0.25608  ,0.099501  ,-0.042839,-0.18382 ,-0.33195 ,-0.49194,-0.66543,-0.85169 ,-1.0483  ,-1.2517  ,-1.4576  ,-1.6616 ,-1.8591 ,-2.0458   ,-2.2178  ,-2.3714  ,-2.5035  ,-2.6112  ,-2.6925  ,-2.7452  ,-2.7679  ,-2.7595  ,-2.7192  ,-2.6464  ,-2.5409  ,-2.4029  ,-2.2327  ,-2.0311   ,-1.799   ,-1.5377  ,-1.2489  ,-0.93435  ,-0.59633 ,-0.23731 ,0.13994  ,0.53236    ,
                                                                    19.743,16.514,13.615,11.052,8.8272,6.9316 ,5.3497  ,4.0587  ,3.0295   ,2.2286  ,1.6197   ,1.1661  ,0.83202 ,0.58478 ,0.39569  ,0.24096   ,0.10191  ,-0.035078,-0.17906 ,-0.3351 ,-0.50506,-0.68827 ,-0.88228 ,-1.0834  ,-1.2874  ,-1.4896 ,-1.6853 ,-1.87     ,-2.0397  ,-2.1906  ,-2.3195  ,-2.4236  ,-2.5006  ,-2.5485  ,-2.566   ,-2.5518  ,-2.5052  ,-2.4259  ,-2.3139  ,-2.1694  ,-1.9931  ,-1.7857   ,-1.5487  ,-1.2836  ,-0.99216 ,-0.67665  ,-0.33953 ,0.016425 ,0.38818  ,0.77248    ,
                                                                    20.452,17.152,14.183,11.556,9.2708,7.322  ,5.6941  ,4.3646  ,3.304    ,2.4783  ,1.8506   ,1.383   ,1.039   ,0.78488 ,0.59127  ,0.43365   ,0.29279  ,0.15463  ,0.0097939,-0.14699,-0.31767,-0.50164 ,-0.69643 ,-0.89831 ,-1.1029  ,-1.3054 ,-1.5011 ,-1.6854   ,-1.8542  ,-2.0036  ,-2.1304  ,-2.2318  ,-2.3054  ,-2.3493  ,-2.3623  ,-2.3432  ,-2.2914  ,-2.2067  ,-2.0893  ,-1.9397  ,-1.7586  ,-1.5473   ,-1.3073  ,-1.0403  ,-0.74865 ,-0.43468  ,-0.10118 ,0.24883  ,0.6121   ,0.98515    ,
                                                                    21.14 ,17.782,14.759,12.08 ,9.7488,7.7581 ,6.0932  ,4.7312  ,3.6425   ,2.7927  ,2.1445   ,1.6597  ,1.3014  ,1.0354  ,0.83197  ,0.66618   ,0.51846  ,0.37443  ,0.2245   ,0.063211,-0.11151,-0.29912 ,-0.49715 ,-0.70188 ,-0.90884 ,-1.1132 ,-1.3103 ,-1.4953   ,-1.6641  ,-1.813   ,-1.9385  ,-2.0379  ,-2.1089  ,-2.1497  ,-2.1589  ,-2.1358  ,-2.0798  ,-1.9909  ,-1.8694  ,-1.716   ,-1.5319  ,-1.3184   ,-1.0774  ,-0.81086 ,-0.52133 ,-0.2115   ,0.11563  ,0.45682  ,0.80864  ,1.1675     ,
                                                                    21.838,18.437,15.372,12.653,10.284,8.2568 ,6.5574  ,5.163   ,4.0441   ,3.1667  ,2.4934   ,1.9862  ,1.6081  ,1.3249  ,1.1065   ,0.92791   ,0.7691   ,0.61546  ,0.45716  ,0.28851 ,0.10729 ,-0.086044,-0.28911 ,-0.49821 ,-0.70886 ,-0.91627,-1.1156 ,-1.3022   ,-1.4719  ,-1.6209  ,-1.7458  ,-1.8439  ,-1.913   ,-1.9514  ,-1.9578  ,-1.9316  ,-1.8725  ,-1.7805  ,-1.6563  ,-1.5008  ,-1.3153  ,-1.1016   ,-0.86169 ,-0.59793 ,-0.31303 ,-0.0099652,0.30807  ,0.63764  ,0.9752   ,1.317      ,
                                                                    22.575,19.142,16.043,13.29 ,10.884,8.8203 ,7.0841  ,5.6535  ,4.4999   ,3.5895  ,2.8857   ,2.3508  ,1.948   ,1.6429  ,1.4055   ,1.2103    ,1.0371   ,0.87105  ,0.70191  ,0.52381 ,0.33432 ,0.13372  ,-0.075686,-0.29027 ,-0.50559 ,-0.71686,-0.91926,-1.1082   ,-1.2793  ,-1.429   ,-1.554   ,-1.6515  ,-1.7194  ,-1.7562  ,-1.7608  ,-1.7325  ,-1.6714  ,-1.5777  ,-1.4523  ,-1.2963  ,-1.1113  ,-0.89941  ,-0.6628  ,-0.40414 ,-0.12636 ,0.16737   ,0.47367  ,0.78899  ,1.1097   ,1.4321     ,
                                                                    23.371,19.909,16.778,13.989,11.545,9.441  ,7.6635  ,6.1918  ,4.9981   ,4.0498  ,3.3107   ,2.7436  ,2.3119  ,1.9814  ,1.7216   ,1.507     ,1.317    ,1.1363   ,0.95455  ,0.76543 ,0.56633 ,0.35733  ,0.1406   ,-0.080327,-0.30107 ,-0.51689,-0.72298,-0.91477  ,-1.088   ,-1.2391  ,-1.3647  ,-1.4623  ,-1.5299  ,-1.5659  ,-1.5695  ,-1.5404  ,-1.4785  ,-1.3845  ,-1.2594  ,-1.1047  ,-0.92219 ,-0.71412  ,-0.48304 ,-0.23182 ,0.036423 ,0.31835   ,0.61046  ,0.90913  ,1.2107   ,1.5115     ,
                                                                    24.227,20.737,17.572,14.743,12.257,10.107 ,8.2839  ,6.7663  ,5.5282   ,4.5376  ,3.7594   ,3.1566  ,2.6929  ,2.334   ,2.0495   ,1.8133    ,1.6046   ,1.4078   ,1.212    ,1.0106  ,0.8009  ,0.5826   ,0.35776  ,0.1298   ,-0.097004,-0.31796,-0.52831,-0.72354  ,-0.89948 ,-1.0525  ,-1.1794  ,-1.2779  ,-1.3459  ,-1.3821  ,-1.3859  ,-1.3569  ,-1.2957  ,-1.2028  ,-1.0797  ,-0.92808 ,-0.74996 ,-0.54783  ,-0.32447 ,-0.082932,0.17349  ,0.44135   ,0.71706  ,0.99698  ,1.2775   ,1.555      ,
                                                                    25.139,21.616,18.412,15.541,13.007,10.808 ,8.9349  ,7.3676  ,6.0815   ,5.0455  ,4.2251   ,3.584   ,3.0858  ,2.6965  ,2.3852   ,2.1258    ,1.897    ,1.6828   ,1.4718   ,1.2574  ,1.0361  ,0.8078   ,0.57419  ,0.33858  ,0.10516  ,-0.12147,-0.33663,-0.53585  ,-0.71504 ,-0.87063 ,-0.99962 ,-1.0996  ,-1.1689  ,-1.2063  ,-1.2113  ,-1.1839  ,-1.1246  ,-1.0345  ,-0.91506 ,-0.76827 ,-0.59645 ,-0.40228  ,-0.18871 ,0.041032 ,0.28355  ,0.53531   ,0.79272  ,1.0522   ,1.3101   ,1.5632     ,
                                                                    26.093,22.534,19.288,16.37 ,13.786,11.535 ,9.6079  ,7.9883  ,6.6515   ,5.5676  ,4.7029   ,4.0213  ,3.4869  ,3.0655  ,2.7261   ,2.4421    ,2.192    ,1.9593   ,1.7323   ,1.5039  ,1.2704  ,1.0314   ,0.78848  ,0.54467  ,0.3041   ,0.071276,-0.14921,-0.35297  ,-0.536   ,-0.69482 ,-0.82657 ,-0.92899 ,-1.0005  ,-1.0401  ,-1.0475  ,-1.0228  ,-0.96688 ,-0.88101 ,-0.76691 ,-0.62677 ,-0.4631  ,-0.27879  ,-0.076934,0.13912  ,0.36591  ,0.59988   ,0.83746  ,1.0752   ,1.3095   ,1.5374     ,
                                                                    27.079,23.481,20.19 ,17.222,14.584,12.279 ,10.297  ,8.6227  ,7.2334   ,6.1     ,5.1892   ,4.4657  ,3.8935  ,3.4386  ,3.0697   ,2.7601    ,2.4877   ,2.2356   ,1.9918   ,1.7487  ,1.5025  ,1.2522   ,0.99936  ,0.74685  ,0.49863  ,0.25911 ,0.032777,-0.17607  ,-0.36354 ,-0.52627 ,-0.6615  ,-0.76716 ,-0.84179 ,-0.88463 ,-0.89553 ,-0.87492 ,-0.8238  ,-0.74368 ,-0.63651 ,-0.50468 ,-0.35089 ,-0.17816  ,0.01027  ,0.21101  ,0.42056  ,0.6354    ,0.85205  ,1.0671   ,1.2774   ,1.4799     ,
                                                                    28.087,24.447,21.109,18.089,15.398,13.036 ,10.997  ,9.267   ,7.8238   ,6.6394  ,5.6813   ,4.9145  ,4.3034  ,3.8138  ,3.4145   ,3.0782    ,2.7826   ,2.5104   ,2.2491   ,1.9907  ,1.731   ,1.469    ,1.2057   ,0.944    ,0.68766  ,0.44096 ,0.20827 ,-0.0062193,-0.19872 ,-0.36601 ,-0.50548 ,-0.6152  ,-0.6939  ,-0.74102 ,-0.75659 ,-0.74128 ,-0.6963  ,-0.62337 ,-0.52463 ,-0.40264 ,-0.26027 ,-0.10065  ,0.072934 ,0.25704  ,0.4482   ,0.64297   ,0.838    ,1.0301   ,1.2162   ,1.3936     ,
                                                                    29.108,25.426,22.04 ,18.968,16.22 ,13.801 ,11.705  ,9.9181  ,8.42     ,7.1835  ,6.177    ,5.3658  ,4.7147  ,4.1896  ,3.7588   ,3.3951    ,3.0755   ,2.7824   ,2.5028   ,2.2285  ,1.9549  ,1.6806   ,1.4065   ,1.1351   ,0.87022  ,0.61587 ,0.37634 ,0.15569   ,-0.042433,-0.21493 ,-0.35935 ,-0.47394 ,-0.55762 ,-0.61001 ,-0.63135 ,-0.62253 ,-0.58494 ,-0.5205  ,-0.43155 ,-0.32076 ,-0.19113 ,-0.045844 ,0.11175  ,0.27825  ,0.45025  ,0.62442   ,0.79757  ,0.9667   ,1.129    ,1.2821     ,
                                                                    30.138,26.412,22.978,19.853,17.049,14.573 ,12.418  ,10.574  ,9.0198   ,7.7304  ,6.6745   ,5.8181  ,5.1261  ,4.5644  ,4.1014   ,3.7093    ,3.365    ,3.0503   ,2.752    ,2.4613  ,2.1732  ,1.8862   ,1.6007   ,1.3194   ,1.0455   ,0.78308 ,0.53626 ,0.30896   ,0.10468  ,-0.073634,-0.22369 ,-0.34391 ,-0.43341 ,-0.492   ,-0.52013 ,-0.51886 ,-0.4898  ,-0.43502 ,-0.357   ,-0.25854 ,-0.14271 ,-0.012728 ,0.12807  ,0.27634  ,0.42879  ,0.58223   ,0.73365  ,0.88029  ,1.0196   ,1.1495     ,
                                                                    31.171,27.403,23.919,20.741,17.882,15.347 ,13.134  ,11.231  ,9.6213   ,8.2783  ,7.1723   ,6.2698  ,5.536   ,4.937   ,4.4409   ,4.0198    ,3.6502   ,3.3133   ,2.9957   ,2.688   ,2.3851  ,2.0849   ,1.7878   ,1.4961   ,1.2129   ,0.94203 ,0.68754 ,0.45315   ,0.24224  ,0.057569 ,-0.098738,-0.22528 ,-0.32136 ,-0.38698 ,-0.42278 ,-0.43002 ,-0.41045 ,-0.3663  ,-0.30016 ,-0.21491 ,-0.11367 ,0.00036037,0.12389  ,0.25366  ,0.38652  ,0.51945   ,0.64968  ,0.77464  ,0.89209  ,1.0001     ,
                                                                    32.205,28.394,24.862,21.631,18.716,16.123 ,13.851  ,11.89   ,10.223   ,8.8257  ,7.6689   ,6.7196  ,5.9434  ,5.3062  ,4.7764   ,4.3256    ,3.9299   ,3.5704   ,3.233    ,2.908   ,2.5899  ,2.2762   ,1.9672   ,1.6648   ,1.372    ,1.0925  ,0.82998 ,0.58817   ,0.37023  ,0.17875  ,0.015674 ,-0.11777 ,-0.22108 ,-0.29444 ,-0.33868 ,-0.3552  ,-0.34591 ,-0.31316 ,-0.25961 ,-0.18821 ,-0.10206 ,-0.0043353,0.10175  ,0.21308  ,0.32665  ,0.43964   ,0.54947  ,0.65385  ,0.75077  ,0.83853    ,
                                                                    33.237,29.383,25.803,22.52 ,19.549,16.898 ,14.567  ,12.547  ,10.823   ,9.3712  ,8.163    ,7.1663  ,6.3469  ,5.671   ,5.1068   ,4.6257    ,4.2036   ,3.8209   ,3.4633   ,3.1207  ,2.7871  ,2.4598   ,2.1387   ,1.8254   ,1.5229   ,1.2344  ,0.96377 ,0.71429   ,0.48902  ,0.29039  ,0.12014  ,-0.020673,-0.13173 ,-0.2134  ,-0.26667 ,-0.29308 ,-0.29467 ,-0.27387 ,-0.2334  ,-0.17622 ,-0.1054  ,-0.024064 ,0.064706 ,0.15792  ,0.25276  ,0.34661   ,0.43712  ,0.52221  ,0.60011  ,0.66935    ,
                                                                    34.265,30.369,26.742,23.406,20.38 ,17.671 ,15.281  ,13.202  ,11.421   ,9.9135  ,8.6535   ,7.6087  ,6.7456  ,6.0304  ,5.4313   ,4.9195    ,4.4704   ,4.0643   ,3.6862   ,3.3258  ,2.9766  ,2.6355   ,2.3023   ,1.9781   ,1.6658   ,1.3684  ,1.0895  ,0.83222   ,0.59946  ,0.39344  ,0.21575  ,0.067227 ,-0.051967,-0.14236 ,-0.20509 ,-0.24183 ,-0.2547  ,-0.2462  ,-0.21908 ,-0.17627 ,-0.12079 ,-0.055674 ,0.016132 ,0.091799 ,0.16869  ,0.2444    ,0.31678  ,0.38399  ,0.44446  ,0.49694    ,
                                                                    35.285,31.348,27.674,24.288,21.206,18.439 ,15.99   ,13.853  ,12.014   ,10.451  ,9.1391   ,8.0459  ,7.1386  ,6.3837  ,5.7492   ,5.2063    ,4.7301   ,4.3003   ,3.9017   ,3.5233  ,3.1585  ,2.8038   ,2.4585   ,2.1236   ,1.8015   ,1.4953  ,1.2082  ,0.94313   ,0.70284  ,0.48935  ,0.30406  ,0.14763  ,0.020063 ,-0.079327,-0.1518  ,-0.19914 ,-0.22353 ,-0.2275  ,-0.2138  ,-0.18532 ,-0.14499 ,-0.095735 ,-0.040359,0.018482 ,0.078341 ,0.13701   ,0.19256  ,0.24334  ,0.28798  ,0.3254     ,
                                                                    36.297,32.32 ,28.6  ,25.163,22.026,19.202 ,16.695  ,14.499  ,12.602   ,10.984  ,9.6189   ,8.477   ,7.5252  ,6.7302  ,6.0602   ,5.486     ,4.9825   ,4.529    ,4.1099   ,3.7138  ,3.3336  ,2.9654   ,2.6082   ,2.2628   ,1.9313   ,1.6165  ,1.3214  ,1.0487    ,0.80097  ,0.58004  ,0.38712  ,0.22272  ,0.086654 ,-0.021876,-0.10425 ,-0.16231 ,-0.19831 ,-0.21477 ,-0.21441 ,-0.20006 ,-0.17456 ,-0.14068  ,-0.10109 ,-0.058251,-0.014438,0.028343  ,0.068347 ,0.10411  ,0.13443  ,0.1584     ,
                                                                    37.299,33.282,29.517,26.03 ,22.839,19.959 ,17.393  ,15.138  ,13.183   ,11.51   ,10.092   ,8.9014  ,7.9049  ,7.0698  ,6.3641   ,5.7586    ,5.2279   ,4.7509   ,4.3115   ,3.898   ,3.5027  ,3.1215   ,2.7529   ,2.3975   ,2.0569   ,1.7339  ,1.4311  ,1.1511    ,0.89612  ,0.6679   ,0.46744  ,0.29508  ,0.1505   ,0.032786 ,-0.059527,-0.12835 ,-0.17594 ,-0.2048  ,-0.21762 ,-0.21711 ,-0.20601 ,-0.18696  ,-0.16244 ,-0.13476 ,-0.106   ,-0.077994 ,-0.052303,-0.030233,-0.012826,-0.00087356,
                                                                    38.288,34.233,30.424,26.887,23.643,20.706 ,18.082  ,15.769  ,13.757   ,12.028  ,10.558   ,9.3187  ,8.2777  ,7.4025  ,6.6613   ,6.0248    ,5.4671   ,4.967    ,4.5078   ,4.0773  ,3.6675  ,3.2739   ,2.8945   ,2.5296   ,2.1806   ,1.8499  ,1.5398  ,1.2529    ,0.99099  ,0.75571  ,0.54786  ,0.36764  ,0.2146   ,0.087717 ,-0.014527,-0.094075,-0.15318 ,-0.19431 ,-0.22007 ,-0.23308 ,-0.23594 ,-0.23114  ,-0.221   ,-0.20766 ,-0.19302 ,-0.17875  ,-0.16625 ,-0.15669 ,-0.15097 ,-0.14978   ,
                                                                    39.263,35.17 ,31.318,27.733,24.437,21.444 ,18.763  ,16.393  ,14.323   ,12.539  ,11.017   ,9.7291  ,8.6437  ,7.7288  ,6.9525   ,6.2853    ,5.7013   ,5.1786   ,4.7003   ,4.2536  ,3.8301  ,3.4248   ,3.0354   ,2.6618   ,2.3051   ,1.9673  ,1.6506  ,1.3571    ,1.0887   ,0.84659  ,0.63154  ,0.44357  ,0.28213  ,0.14612  ,0.033973 ,-0.056261,-0.1268  ,-0.18006 ,-0.21855 ,-0.24476 ,-0.26116 ,-0.27008  ,-0.27368 ,-0.27394 ,-0.27259 ,-0.27114  ,-0.27085 ,-0.27277 ,-0.27769 ,-0.28623   ,
                                                                    40.221,36.093,32.199,28.566,25.219,22.172 ,19.434  ,17.007  ,14.881   ,13.043  ,11.469   ,10.133  ,9.0038  ,8.0497  ,7.2389   ,6.5418    ,5.9322   ,5.3878   ,4.8913   ,4.4294  ,3.9932  ,3.5772   ,3.1787   ,2.7972   ,2.4336   ,2.0894  ,1.7666  ,1.4671    ,1.1925   ,0.94387  ,0.72179  ,0.52617  ,0.35637  ,0.21123  ,0.089159 ,-0.011759,-0.093702,-0.15899 ,-0.21003 ,-0.24919 ,-0.27878 ,-0.30097  ,-0.31778 ,-0.33101 ,-0.34224 ,-0.35285  ,-0.36397 ,-0.37654 ,-0.39126 ,-0.40868   ,
                                                                    41.163,37    ,33.065,29.386,25.988,22.888 ,20.095  ,17.612  ,15.431   ,13.539  ,11.914   ,10.531  ,9.359   ,8.3666  ,7.5223   ,6.7963    ,6.1621   ,5.5972   ,5.0836   ,4.6076  ,4.1598  ,3.7342   ,3.3277   ,2.9393   ,2.5695   ,2.2197  ,1.8914  ,1.5863    ,1.3059   ,1.051    ,0.82197  ,0.61872  ,0.4405   ,0.28614  ,0.15406  ,0.042371 ,-0.051018,-0.12833 ,-0.19185 ,-0.24379 ,-0.28633 ,-0.32147  ,-0.35106 ,-0.37677 ,-0.40004 ,-0.42211  ,-0.44401 ,-0.46657 ,-0.49045 ,-0.51612   ,
                                                                    42.085,37.889,33.915,30.192,26.745,23.593 ,20.746  ,18.208  ,15.973   ,14.028  ,12.354   ,10.925  ,9.7111  ,8.6817  ,7.8051   ,7.0515    ,6.3941   ,5.8099   ,5.2805   ,4.7918  ,4.3336  ,3.8996   ,3.4862   ,3.0919   ,2.7168   ,2.362   ,2.0287  ,1.7185    ,1.4324   ,1.1713   ,0.93533  ,0.72432  ,0.53752  ,0.37373  ,0.23142  ,0.10876  ,0.003757 ,-0.085694,-0.16172 ,-0.22642 ,-0.28178 ,-0.32967  ,-0.37178 ,-0.40964 ,-0.44455 ,-0.47764  ,-0.50986 ,-0.54196 ,-0.57452 ,-0.60798   ,
                                                                    42.988,38.76 ,34.748,30.982,27.488,24.286 ,21.387  ,18.795  ,16.508   ,14.512  ,12.79    ,11.317  ,10.063  ,8.9977  ,8.0903   ,7.3106    ,6.6315   ,6.0296   ,5.4859   ,4.9858  ,4.5187  ,4.0775   ,3.6582   ,3.2589   ,2.8794   ,2.5201  ,2.1823  ,1.8671    ,1.5754   ,1.308    ,1.0649   ,0.84586  ,0.65012  ,0.47653  ,0.32361  ,0.18963  ,0.072705 ,-0.029128,-0.11785 ,-0.1954  ,-0.2636  ,-0.32418  ,-0.37869 ,-0.42849 ,-0.47481 ,-0.51866  ,-0.5609  ,-0.60221 ,-0.64316 ,-0.68414   ,
                                                                    43.869,39.611,35.564,31.757,28.218,24.968 ,22.018  ,19.376  ,17.038   ,14.994  ,13.225   ,11.709  ,10.417  ,9.318   ,8.3816   ,7.5775    ,6.8784   ,6.2605   ,5.7041   ,5.1942  ,4.7194  ,4.2722   ,3.848    ,3.4446   ,3.0611   ,2.6978  ,2.3556  ,2.0354    ,1.7381   ,1.464    ,1.2134   ,0.9858   ,0.78059  ,0.59664  ,0.43256  ,0.28674  ,0.15744  ,0.042827 ,-0.058907,-0.14954 ,-0.23075 ,-0.30409  ,-0.37099 ,-0.4327  ,-0.49033 ,-0.54479  ,-0.5969  ,-0.64727 ,-0.69644 ,-0.74478   ,
                                                                    44.729,40.444,36.363,32.518,28.936,25.64  ,22.643  ,19.952  ,17.565   ,15.474  ,13.662   ,12.105  ,10.777  ,9.6463  ,8.683    ,7.8565    ,7.1393   ,6.507    ,5.9396   ,5.4212  ,4.9401  ,4.488    ,4.0598   ,3.6528   ,3.2657   ,2.8986  ,2.552   ,2.2266    ,1.923    ,1.6418   ,1.383    ,1.1462   ,0.93073  ,0.73567  ,0.5597   ,0.40137  ,0.25906  ,0.13114  ,0.01595  ,-0.088127,-0.18262 ,-0.26893  ,-0.34836 ,-0.42204 ,-0.49099 ,-0.55607  ,-0.618   ,-0.67739 ,-0.73473 ,-0.79039   ,
                                                                    45.568,41.257,37.145,33.265,29.643,26.304 ,23.262  ,20.525  ,18.093   ,15.958  ,14.104   ,12.509  ,11.147  ,9.9874  ,8.9992   ,8.1524    ,7.4189   ,6.7741   ,6.1972   ,5.6717  ,5.1853  ,4.7291   ,4.2975   ,3.8872   ,3.4967   ,3.1256  ,2.7742  ,2.443     ,2.1326   ,1.8434   ,1.5754   ,1.3284   ,1.1018   ,0.8947   ,0.70595  ,0.53427  ,0.37822  ,0.23633  ,0.1071   ,-0.010896,-0.11905 ,-0.21863  ,-0.31082 ,-0.39665 ,-0.47705 ,-0.55282  ,-0.62464 ,-0.69308 ,-0.75863 ,-0.82167   ,
                                                                    46.385,42.052,37.913,33.999,30.341,26.962 ,23.878  ,21.098  ,18.624   ,16.448  ,14.556   ,12.926  ,11.533  ,10.346  ,9.3352   ,8.4702    ,7.7224   ,7.0666   ,6.4815   ,5.9501  ,5.4592  ,4.9995   ,4.5648   ,4.1513   ,3.7571   ,3.3814  ,3.0245  ,2.6867    ,2.3684   ,2.0701   ,1.7919   ,1.5335   ,1.2946   ,1.0743   ,0.87173  ,0.68573  ,0.51507  ,0.35843  ,0.21452  ,0.082038 ,-0.040247,-0.15348  ,-0.25874 ,-0.35696 ,-0.449   ,-0.53562  ,-0.61746 ,-0.69508 ,-0.76895 ,-0.83947   ,
                                                                    47.183,42.831,38.667,34.725,31.033,27.618 ,24.496  ,21.677  ,19.163   ,16.95   ,15.022   ,13.36   ,11.938  ,10.727  ,9.6963   ,8.815     ,8.0546   ,7.3893   ,6.7973   ,6.2607  ,5.7659  ,5.3029   ,4.865    ,4.4478   ,4.0492   ,3.6682  ,3.3047  ,2.9591    ,2.6316   ,2.3229   ,2.033    ,1.7619   ,1.5093   ,1.2746   ,1.057    ,0.8556   ,0.66936  ,0.49715  ,0.33781  ,0.19022  ,0.053277 ,-0.074054 ,-0.19273 ,-0.30365 ,-0.40759 ,-0.50526  ,-0.5973  ,-0.68426 ,-0.7666  ,-0.84474   ,
                                                                    47.963,43.595,39.411,35.443,31.723,28.275 ,25.118  ,22.264  ,19.716   ,17.468  ,15.508   ,13.817  ,12.37   ,11.137  ,10.088   ,9.1919    ,8.4204   ,7.7468   ,7.1486   ,6.6074  ,6.1088  ,5.6423   ,5.2005   ,4.779    ,4.375    ,3.9874  ,3.6159  ,3.2609    ,2.9227   ,2.6019   ,2.2988   ,2.0134   ,1.7456   ,1.4951   ,1.2612   ,1.0433   ,0.8405   ,0.65182  ,0.47629  ,0.31293  ,0.16077  ,0.018872  ,-0.11363 ,-0.23756 ,-0.35368 ,-0.46265  ,-0.5651  ,-0.66157 ,-0.75256 ,-0.83851   ,
                                                                    48.727,44.348,40.147,36.159,32.414,28.938 ,25.752  ,22.866  ,20.286   ,18.008  ,16.02    ,14.303  ,12.832  ,11.58   ,10.515   ,9.6059    ,8.8244   ,8.1432   ,7.5392   ,6.9934  ,6.4906  ,6.0199   ,5.5735   ,5.1463   ,4.7355   ,4.3396  ,3.9585  ,3.5923    ,3.2416   ,2.9069   ,2.5887   ,2.2874   ,2.0029   ,1.7351   ,1.4837   ,1.2481   ,1.0276   ,0.82154  ,0.62903  ,0.44923  ,0.28128  ,0.12435   ,-0.022385,-0.15968 ,-0.28825 ,-0.40877  ,-0.52185 ,-0.62805 ,-0.72787 ,-0.8218    ,
                                                                    49.48 ,45.094,40.881,36.876,33.111,29.612 ,26.4    ,23.488  ,20.881   ,18.576  ,16.563   ,14.822  ,13.331  ,12.061  ,10.982   ,10.061    ,9.2707   ,8.5822   ,7.9723   ,7.4213  ,6.9136  ,6.4375   ,5.985    ,5.5505   ,5.131    ,4.7251  ,4.3322  ,3.9528    ,3.5875   ,3.237    ,2.902    ,2.5829   ,2.2801   ,1.9935   ,1.7232   ,1.4687   ,1.2296   ,1.0052   ,0.79496  ,0.59806  ,0.41377  ,0.24132   ,0.079964 ,-0.071025,-0.21235 ,-0.34466  ,-0.46859 ,-0.58471 ,-0.69356 ,-0.79564   ,
                                                                    50.224,45.836,41.617,37.601,33.821,30.303 ,27.071  ,24.136  ,21.506   ,19.178  ,17.142   ,15.382  ,13.872  ,12.586  ,11.493   ,10.562    ,9.7625   ,9.0666   ,8.4502   ,7.893   ,7.379   ,6.896    ,6.4354   ,5.9917   ,5.5614   ,5.1431  ,4.7363  ,4.3415    ,3.9595   ,3.5911   ,3.2372   ,2.8986   ,2.5758   ,2.2691   ,1.9785   ,1.704    ,1.4452   ,1.2016   ,0.97288  ,0.75826  ,0.55709  ,0.36868   ,0.19233  ,0.027328 ,-0.12701 ,-0.27135  ,-0.40634 ,-0.53257 ,-0.65063 ,-0.76104   ,
                                                                    50.966,46.581,42.359,38.338,34.548,31.017 ,27.768  ,24.815  ,22.166   ,19.819  ,17.764   ,15.985  ,14.459  ,13.159  ,12.054   ,11.112    ,10.303   ,9.5984   ,8.9743   ,8.4095  ,7.8873  ,7.3953   ,6.9246   ,6.4691   ,6.0257   ,5.5925  ,5.1695  ,4.7569    ,4.3559   ,3.9675   ,3.5929   ,3.233    ,2.8886   ,2.5602   ,2.2481   ,1.9524   ,1.673    ,1.4095   ,1.1616   ,0.92863  ,0.71009  ,0.50532   ,0.31363  ,0.13433  ,-0.033279,-0.18987  ,-0.3361  ,-0.47262 ,-0.60003 ,-0.71893   ,
                                                                    51.709,47.333,43.116,39.093,35.299,31.76  ,28.499  ,25.533  ,22.868   ,20.504  ,18.433   ,16.638  ,15.097  ,13.783  ,12.665   ,11.712    ,10.893   ,10.179   ,9.5454   ,8.9709  ,8.4383  ,7.9349   ,7.4514   ,6.9817   ,6.5224   ,6.0718  ,5.6299  ,5.1972    ,4.775    ,4.3645   ,3.9672   ,3.5843   ,3.2167   ,2.8653   ,2.5305   ,2.2127   ,1.9117   ,1.6275   ,1.3598   ,1.108    ,0.87165  ,0.65014   ,0.44281  ,0.24895  ,0.06786  ,-0.10116  ,-0.2588  ,-0.40576 ,-0.54268 ,-0.67019   ,
                                                                    52.462,48.099,43.891,39.874,36.08 ,32.537 ,29.27   ,26.293  ,23.616   ,21.238  ,19.153   ,17.343  ,15.788  ,14.46   ,13.33    ,12.364    ,11.533   ,10.808   ,10.163   ,9.5765  ,9.031   ,8.5134   ,8.0142   ,7.5275   ,7.0495   ,6.5789  ,6.1156  ,5.6603    ,5.2146   ,4.78     ,4.3582   ,3.9506   ,3.5585   ,3.1828   ,2.8242   ,2.4832   ,2.16     ,1.8544   ,1.5662   ,1.2952   ,1.0407   ,0.8021    ,0.57886  ,0.37023  ,0.17547  ,-0.0061352,-0.17533 ,-0.33285 ,-0.47939 ,-0.61565};
        }
        case 10:
        {
            static const std::vector<double> arrayCurvatureMin   = {15.285,13.216,11.458,9.9932,8.7933,7.8225,7.0439,6.4226,5.9274,5.532 ,5.2143,4.9563,4.7439,4.5657,4.4127,4.2781,4.1563,4.0431,3.9351,3.8298,3.725 ,3.6191,3.5108,3.3991,3.283 ,3.1618,3.0351,2.9024,2.7636,2.6184,2.4669,2.3097,2.1472 ,1.9808 ,1.8123 ,1.6443 ,1.4809 ,1.3276 ,1.1918 ,1.0824 ,1.0099 ,0.98353,1.0086 ,1.082  ,1.189  ,1.3048 ,1.4009    ,1.4536    ,1.4499   ,1.3883   ,
                                                                    15.817,13.674,11.842,10.309,9.0484,8.0254,7.203 ,6.5456,6.0212,5.602 ,5.2653,4.992 ,4.7672,4.5789,4.4177,4.2761,4.1485,4.0302,3.9177,3.8083,3.6998,3.5903,3.4786,3.3636,3.2442,3.1199,2.9901,2.8543,2.7124,2.5642,2.4099,2.25  ,2.0853 ,1.9171 ,1.7475 ,1.5796 ,1.4177 ,1.2678 ,1.138  ,1.0376 ,0.97708,0.96496,1.0045 ,1.0892 ,1.2013 ,1.3139 ,1.3991    ,1.436     ,1.415    ,1.3373   ,
                                                                    16.37 ,14.153,12.246,10.643,9.3197,8.2421,7.3737,6.6783,6.1227,5.6785,5.3214,5.0318,4.7939,4.5948,4.4248,4.2758,4.1419,4.0182,3.9009,3.7872,3.6747,3.5616,3.4463,3.3279,3.2052,3.0776,2.9445,2.8056,2.6605,2.5093,2.3523,2.1898,2.0229 ,1.853  ,1.6826 ,1.515  ,1.3551 ,1.2094 ,1.0864 ,0.99595,0.94821,0.95065,1.0041 ,1.0987 ,1.2136 ,1.3206 ,1.3931    ,1.4133    ,1.3752   ,1.2823   ,
                                                                    16.945,14.654,12.672,10.997,9.608 ,8.4736,7.5569,6.8213,6.2329,5.7619,5.3833,5.0763,4.8243,4.6138,4.4342,4.2774,4.1368,4.0073,3.885 ,3.7667,3.6501,3.533 ,3.4141,3.292 ,3.1658,3.0349,2.8985,2.7563,2.6081,2.4539,2.294 ,2.129 ,1.9599 ,1.7887 ,1.6177 ,1.4509 ,1.2935 ,1.1527 ,1.0374 ,0.95772,0.92347,0.94058,1.0071 ,1.1098 ,1.225  ,1.3242 ,1.3825    ,1.3857    ,1.3308   ,1.2235   ,
                                                                    17.542,15.178,13.119,11.371,9.9146,8.7209,7.7536,6.9757,6.3524,5.8531,5.4516,5.126 ,4.8589,4.6361,4.4465,4.2813,4.1335,3.9979,3.8702,3.747 ,3.6259,3.5047,3.3819,3.2561,3.1264,2.9919,2.8521,2.7065,2.5551,2.3979,2.2352,2.0677,1.8967 ,1.7241 ,1.5529 ,1.3873 ,1.2331 ,1.0979 ,0.99117,0.9232 ,0.90299,0.93466,1.0132 ,1.122  ,1.2351 ,1.3242 ,1.3671    ,1.3531    ,1.282    ,1.1612   ,
                                                                    18.162,15.725,13.59 ,11.766,10.241,8.9853,7.9649,7.1424,6.4824,5.953 ,5.5269,5.1816,4.8984,4.6624,4.462 ,4.2877,4.1323,3.9901,3.8566,3.7284,3.6025,3.477 ,3.3501,3.2204,3.0869,2.9487,2.8054,2.6565,2.5018,2.3415,2.1759,2.006 ,1.8332 ,1.6596 ,1.4885 ,1.3245 ,1.1741 ,1.0454 ,0.94813,0.89264,0.88688,0.93271,1.0218 ,1.1344 ,1.243  ,1.3201 ,1.3469    ,1.3156    ,1.2291   ,1.0959   ,
                                                                    18.805,16.297,14.085,12.185,10.588,9.2681,8.1921,7.3227,6.6236,6.0623,5.6102,5.2437,4.9433,4.6932,4.4812,4.2972,4.1336,3.9844,3.8447,3.7109,3.5801,3.45  ,3.3187,3.1849,3.0475,2.9056,2.7586,2.6061,2.4481,2.2847,2.1164,1.9441,1.7695 ,1.5951 ,1.4245 ,1.2628 ,1.1168 ,0.99539,0.90857,0.86627,0.87516,0.93448,1.0325 ,1.1465 ,1.2481 ,1.3118 ,1.3217    ,1.2736    ,1.1724   ,1.0278   ,
                                                                    19.472,16.894,14.606,12.628,10.957,9.5706,8.4365,7.5177,6.7774,6.1821,5.7023,5.3132,4.9943,4.7292,4.5046,4.3103,4.1379,3.9811,3.8348,3.6951,3.559 ,3.4239,3.2881,3.15  ,3.0084,2.8626,2.7118,2.5557,2.3943,2.2278,2.0566,1.8821,1.706  ,1.5311 ,1.3613 ,1.2024 ,1.0616 ,0.94833,0.87277,0.84427,0.86779,0.93961,1.0446 ,1.1575 ,1.25   ,1.2989 ,1.2916    ,1.2272    ,1.1123   ,0.95718  ,
                                                                    20.164,17.518,15.153,13.096,11.35 ,9.8944,8.6996,7.7287,6.9449,6.3135,5.8042,5.3909,5.0522,4.7709,4.533 ,4.3274,4.1456,3.9806,3.8272,3.6813,3.5394,3.3992,3.2585,3.1158,2.9699,2.8199,2.6652,2.5054,2.3405,2.1708,1.9969,1.8202,1.6427 ,1.4675 ,1.2991 ,1.1435 ,1.0088 ,0.9045 ,0.84103,0.82679,0.86465,0.94768,1.0574 ,1.1669 ,1.2483 ,1.2812 ,1.2569    ,1.1768    ,1.0491   ,0.8845   ,
                                                                    20.881,18.169,15.728,13.592,11.768,10.241,8.9828,7.9574,7.1275,6.4578,5.917 ,5.4778,5.118 ,4.8193,4.5669,4.3493,4.1573,3.9835,3.8225,3.6698,3.5218,3.376 ,3.2301,3.0826,2.9322,2.7779,2.619 ,2.4554,2.2869,2.114 ,1.9373,1.7585,1.5799 ,1.4048 ,1.2381 ,1.0865 ,0.9586 ,0.86427,0.81362,0.81392,0.86554,0.9582 ,1.0704 ,1.1739 ,1.2424 ,1.2588 ,1.2175    ,1.1225    ,0.98307  ,0.81002  ,
                                                                    21.623,18.848,16.332,14.116,12.213,10.612,9.288 ,8.2051,7.3266,6.6164,6.0419,5.5752,5.1927,4.8752,4.6073,4.3767,4.1736,3.9905,3.8212,3.6612,3.5067,3.3548,3.2035,3.0508,2.8955,2.7366,2.5736,2.4059,2.2337,2.0575,1.8782,1.6974,1.5177 ,1.3431 ,1.1787 ,1.0318 ,0.91153,0.82798,0.79079,0.80567,0.87019,0.97061,1.0827 ,1.1781 ,1.2323 ,1.2317 ,1.1739    ,1.0649    ,0.91465  ,0.73406  ,
                                                                    22.391,19.556,16.966,14.67 ,12.686,11.01 ,9.6166,8.4737,7.5439,6.7906,6.1804,5.6841,5.2773,4.9397,4.6551,4.4104,4.1955,4.0021,3.8239,3.656 ,3.4944,3.3362,3.1789,3.0207,2.8603,2.6966,2.529 ,2.3572,2.1812,2.0017,1.8196,1.637 ,1.4566 ,1.2827 ,1.1212 ,0.97977,0.86794,0.79597,0.77274,0.80202,0.87825,0.98431,1.0939 ,1.1789 ,1.2177 ,1.2    ,1.1263    ,1.0041    ,0.84416  ,0.65694  ,
                                                                    23.185,20.293,17.631,15.255,13.19 ,11.435,9.9707,8.7649,7.781 ,6.9821,6.3338,5.806 ,5.3731,5.0139,4.7113,4.4514,4.2236,4.0191,3.8313,3.6548,3.4856,3.3205,3.1568,2.9928,2.8269,2.6581,2.4858,2.3096,2.1297,1.9468,1.762 ,1.5776,1.3968 ,1.2241 ,1.066  ,0.93075,0.82823,0.76859,0.75964,0.80283,0.88928,0.99868,1.1033 ,1.1761 ,1.1986 ,1.1639 ,1.075     ,0.94072   ,0.77195  ,0.57897  ,
                                                                    24.005,21.06 ,18.328,15.873,13.725,11.89 ,10.352,9.0804,8.0396,7.1925,6.5038,5.9424,5.4815,5.099 ,4.7769,4.5007,4.2589,4.0425,3.8442,3.6585,3.481 ,3.3084,3.1379,2.9676,2.7959,2.6217,2.4443,2.2636,2.0796,1.8932,1.7057,1.5197,1.3386 ,1.1675 ,1.0135 ,0.8852 ,0.79282,0.74616,0.7516 ,0.80792,0.9028 ,1.0131 ,1.1104 ,1.1694 ,1.1751 ,1.1236 ,1.0203    ,0.87497   ,0.69836  ,0.50045  ,
                                                                    24.849,21.856,19.058,16.524,14.293,12.377,10.762,9.4221,8.3217,7.4236,6.692 ,6.0947,5.6039,5.1964,4.8534,4.5595,4.3026,4.0731,3.8635,3.6678,3.4814,3.3007,3.1228,2.9457,2.7677,2.5877,2.405 ,2.2194,2.0312,1.8412,1.6511,1.4635,1.2825 ,1.1135 ,0.96417,0.84362,0.76215,0.72897,0.74868,0.81703,0.9183 ,1.027  ,1.1147 ,1.1587 ,1.1472 ,1.0796 ,0.9628    ,0.80729   ,0.62376  ,0.42169  ,
                                                                    25.718,22.681,19.819,17.209,14.895,12.896,11.203,9.7917,8.6289,7.6772,6.9001,6.2646,5.7418,5.3076,4.942 ,4.629 ,4.3558,4.1122,3.8903,3.6837,3.4876,3.2981,3.1122,2.9277,2.743 ,2.5568,2.3684,2.1777,1.985 ,1.7914,1.5986,1.4097,1.229  ,1.0626 ,0.91855,0.80651,0.73667,0.71731,0.75086,0.82984,0.93522,1.0398 ,1.1159 ,1.1438 ,1.1152 ,1.032  ,0.90272   ,0.73806   ,0.54848  ,0.34297  ,
                                                                    26.609,23.534,20.613,17.928,15.532,13.449,11.675,10.191,8.963 ,7.9549,7.1298,6.4538,5.8969,5.434 ,5.0442,4.7106,4.4198,4.161 ,3.9256,3.7072,3.5005,3.3015,3.1069,2.9145,2.7225,2.5296,2.3352,2.1391,1.9417,1.7442,1.5487,1.3586,1.1785 ,1.0152 ,0.87721,0.77441,0.71682,0.71142,0.7581 ,0.84598,0.95303,1.051  ,1.1138 ,1.1249 ,1.0793 ,0.98139,0.84055   ,0.66767   ,0.47287  ,0.26462  ,
                                                                    27.52 ,24.413,21.437,18.681,16.204,14.036,12.181,10.621,9.3256,8.2585,7.3829,6.6639,6.0708,5.5773,5.1615,4.8058,4.4959,4.2206,3.9708,3.7395,3.5213,3.3119,3.108 ,2.907 ,2.7071,2.507 ,2.3061,2.1041,1.9018,1.7003,1.5022,1.3109,1.1318 ,0.97211,0.84076,0.74792,0.70304,0.71149,0.77027,0.86505,0.97121,1.0602 ,1.1083 ,1.1022 ,1.04   ,0.92815,0.77671   ,0.59653   ,0.39729  ,0.18693  ,
                                                                    28.447,25.315,22.291,19.467,16.91 ,14.659,12.721,11.083,9.718 ,8.5895,7.6608,6.8966,6.2651,5.739 ,5.2955,4.9161,4.5857,4.2925,4.027 ,3.7819,3.5512,3.3305,3.1164,2.906 ,2.6976,2.4897,2.2818,2.0737,1.8661,1.6605,1.4595,1.2673,1.0895 ,0.93391,0.80988,0.72763,0.69579,0.7177 ,0.78721,0.88661,0.9893 ,1.0672 ,1.0993 ,1.0758 ,0.99766,0.87274,0.71167   ,0.52506   ,0.3221   ,0.11022  ,
                                                                    29.388,26.238,23.17 ,20.283,17.65 ,15.316,13.295,11.579,10.141,8.9489,7.9649,7.1532,6.4812,5.9206,5.4476,5.0429,4.6906,4.3782,4.0958,3.8356,3.5914,3.3585,3.1332,2.9128,2.6951,2.4788,2.2633,2.0485,1.8353,1.6253,1.4216,1.2285,1.0524 ,0.90139,0.78532,0.71418,0.69548,0.73016,0.80872,0.91028,1.0069 ,1.0718 ,1.087  ,1.0462 ,0.95272,0.81568,0.64589   ,0.45365   ,0.24766  ,0.034834 ,
                                                                    30.337,27.177,24.073,21.128,18.423,16.007,13.903,12.106,10.595,9.3376,8.2963,7.4351,6.7206,6.1236,5.6194,5.1877,4.812 ,4.4791,4.1786,3.9021,3.6434,3.3972,3.1599,2.9284,2.7007,2.4753,2.2516,2.0297,1.8105,1.5959,1.3893,1.1955,1.0213 ,0.87544,0.76788,0.70825,0.70255,0.74898,0.83457,0.93568,1.0238 ,1.074  ,1.0717 ,1.0138 ,0.90574,0.75748,0.57987   ,0.38277   ,0.17437  ,-0.03889 ,
                                                                    31.29 ,28.127,24.995,21.999,19.224,16.73 ,14.543,12.667,11.081,9.756 ,8.6554,7.7428,6.9841,6.3491,5.8121,5.352 ,4.9515,4.5967,4.2767,3.9829,3.7084,3.448 ,3.1976,2.9543,2.7157,2.4804,2.2479,2.0183,1.7927,1.5732,1.3636,1.1692,0.99726,0.85699,0.75844,0.71053,0.71742,0.77423,0.86455,0.9625 ,1.0397 ,1.0738 ,1.0537 ,0.97911,0.85729,0.69872,0.5141    ,0.31285   ,0.10263  ,-0.11058 ,
                                                                    32.241,29.084,25.932,22.89 ,20.052,17.482,15.215,13.258,11.597,10.204,9.0425,8.077 ,7.2725,6.5979,6.0267,5.5368,5.1101,4.7323,4.3917,4.0793,3.788 ,3.5123,3.2479,2.9917,2.7413,2.4955,2.2535,2.0155,1.783 ,1.5584,1.3456,1.1507,0.98133,0.84708,0.75791,0.72175,0.74048,0.80595,0.89848,0.99054,1.0547 ,1.0716 ,1.0335 ,0.94276,0.80801,0.63999,0.44914   ,0.24438   ,0.032868 ,-0.17983 ,
                                                                    33.184,30.042,26.877,23.798,20.901,18.259,15.914,13.878,12.142,10.68 ,9.4568,8.4372,7.5857,6.8704,6.2638,5.7429,5.289 ,4.8869,4.5246,4.1926,3.8835,3.5914,3.3121,3.0421,2.7791,2.5218,2.2695,2.0227,1.7828,1.5526,1.3366,1.1412,0.97472,0.84685,0.7673 ,0.74264,0.77215,0.84421,0.93622,1.0197 ,1.0689 ,1.0677 ,1.0117 ,0.90544,0.75856,0.58192,0.38555   ,0.17788   ,-0.034441,-0.24622 ,
                                                                    34.116,30.996,27.826,24.717,21.769,19.06 ,16.639,14.525,12.714,11.182,9.8972,8.8228,7.9234,7.1664,6.5235,5.9708,5.4888,5.0616,4.6767,4.3241,3.9962,3.6869,3.3915,3.1069,2.8304,2.5608,2.2975,2.0411,1.7933,1.5573,1.338 ,1.1422,0.97873,0.85752,0.78765,0.77398,0.81283,0.88911,0.97773,1.0501 ,1.0826 ,1.0627 ,0.98898,0.86789,0.7097 ,0.52518,0.32395   ,0.11389   ,-0.098794,-0.30927 ,
                                                                    35.03 ,31.94 ,28.774,25.642,22.649,19.877,17.384,15.195,13.31 ,11.709,10.362,9.232 ,8.2842,7.4851,6.8054,6.2203,5.7096,5.2566,4.8484,4.4746,4.1271,3.7997,3.4876,3.1874,2.8966,2.6139,2.3389,2.0723,1.8162,1.5739,1.3511,1.155 ,0.99478,0.8804 ,0.82008,0.81655,0.86293,0.94078,1.0231 ,1.0818 ,1.0963 ,1.0574 ,0.96622,0.83096,0.66221,0.4705 ,0.26498   ,0.053004  ,-0.15964 ,-0.36845 ,
                                                                    35.922,32.869,29.715,26.569,23.537,20.708,18.147,15.885,13.927,12.257,10.848,9.6627,8.6665,7.8251,7.1084,6.4907,5.951 ,5.4721,5.0402,4.6446,4.277 ,3.9309,3.6014,3.2849,2.9791,2.6827,2.3952,2.1179,1.8528,1.6041,1.3776,1.1812,1.0244 ,0.91688,0.86572,0.87115,0.92291,0.99943,1.0725 ,1.1154 ,1.1107 ,1.0524 ,0.94433,0.79557,0.61695,0.41868,0.20936   ,-0.0041141,-0.21637 ,-0.42319 ,
                                                                    36.788,33.779,30.645,27.491,24.427,21.548,18.922,16.59 ,14.561,12.823,11.352,10.112,9.0676,8.1841,7.4306,6.7806,6.2121,5.7072,5.2516,4.8342,4.4462,4.0811,3.7337,3.4005,3.079 ,2.7682,2.4678,2.1791,1.9048,1.6493,1.419 ,1.2224,1.0691 ,0.9684 ,0.92575,0.93862,0.99324,1.0654 ,1.1263 ,1.1514 ,1.1265 ,1.0489 ,0.92434,0.76272,0.57487,0.37055,0.15788   ,-0.05675  ,-0.26831 ,-0.47285 ,
                                                                    37.626,34.667,31.559,28.405,25.317,22.392,19.706,17.306,15.207,13.404,11.871,10.577,9.4845,8.5594,7.7697,7.0879,6.4911,5.9606,5.4818,5.0427,4.6345,4.2503,3.8849,3.5347,3.1973,2.8717,2.5579,2.2575,1.9736,1.7112,1.4771,1.2803,1.1306 ,1.0364 ,1.0014 ,1.0198 ,1.0745 ,1.1391 ,1.1851 ,1.1907 ,1.1449 ,1.048  ,0.90741,0.73352,0.53697,0.32706,0.11137   ,-0.10411  ,-0.31471 ,-0.51669 ,
                                                                    38.433,35.529,32.453,29.307,26.201,23.236,20.494,18.029,15.863,13.994,12.401,11.053,9.9137,8.9477,8.1225,7.4098,6.7855,6.2305,5.729 ,5.269 ,4.8411,4.4382,4.0549,3.6878,3.3344,2.9939,2.6665,2.3543,2.0607,1.7913,1.5535,1.3566,1.2105 ,1.1224 ,1.0937 ,1.1154 ,1.1672 ,1.2211 ,1.2497 ,1.2341 ,1.167  ,1.051  ,0.8948 ,0.70916,0.50436,0.28921,0.070798  ,-0.14531  ,-0.35474 ,-0.55391 ,
                                                                    39.208,36.363,33.326,30.194,27.075,24.076,21.283,18.757,16.525,14.592,12.939,11.537,10.351,9.3453,8.4856,7.7429,7.0925,6.5139,5.9911,5.5112,5.0645,4.6436,4.2431,3.8595,3.4904,3.1352,2.7944,2.4704,2.1672,1.8909,1.6496,1.453 ,1.3103 ,1.2277 ,1.2039 ,1.2264 ,1.272  ,1.3121 ,1.3209 ,1.283  ,1.1941 ,1.0592 ,0.88789,0.69092,0.47823,0.25812,0.037182  ,-0.17938  ,-0.38745 ,-0.58361 ,
                                                                    39.95 ,37.167,34.174,31.062,27.938,24.91 ,22.07 ,19.485,17.189,15.193,13.481,12.026,10.794,9.7484,8.8552,8.0838,7.4084,6.8078,6.265 ,5.7666,5.3024,4.8648,4.4482,4.0489,3.6648,3.2954,2.9417,2.6063,2.294 ,2.0113,1.7669,1.5707,1.4315 ,1.3537 ,1.3329 ,1.3534 ,1.3897 ,1.4129 ,1.4    ,1.3386 ,1.2278 ,1.0742 ,0.88816,0.68022,0.4599 ,0.235  ,0.011669  ,-0.20524  ,-0.41181 ,-0.60476 ,
                                                                    40.658,37.942,34.997,31.911,28.787,25.736,22.853,20.211,17.854,15.794,14.023,12.516,11.238,10.153,9.2274,8.4284,7.7295,7.1084,6.5473,6.0321,5.5521,5.0993,4.668 ,4.2544,3.8565,3.474 ,3.1082,2.7622,2.4414,2.1531,1.9063,1.711 ,1.5752 ,1.5013 ,1.4816 ,1.4973 ,1.5209 ,1.5246 ,1.488  ,1.4025 ,1.2697 ,1.0977 ,0.89718,0.67853,0.45077,0.22118,-0.0044872,-0.22167  ,-0.42665 ,-0.61624 ,
                                                                    41.333,38.686,35.794,32.738,29.621,26.551,23.629,20.934,18.516,16.395,14.565,13.004,11.68 ,10.557,9.5987,8.773 ,8.0518,7.4118,6.8341,6.304 ,5.8101,5.344 ,4.8998,4.4736,4.0635,3.6693,3.2928,2.9376,2.6096,2.3167,2.0685,1.8748,1.7424 ,1.6713 ,1.6506 ,1.6584 ,1.6663 ,1.648  ,1.5865 ,1.4762 ,1.3215 ,1.1313 ,0.91663,0.68746,0.45234,0.21809,-0.009916 ,-0.22737  ,-0.4307  ,-0.61676 ,
                                                                    41.976,39.399,36.564,33.544,30.438,27.355,24.398,21.653,19.175,16.992,15.103,13.488,12.118,10.956,9.9657,9.114 ,8.3716,7.714 ,7.1215,6.5784,6.0727,5.5954,5.1404,4.7037,4.2834,3.8794,3.494 ,3.1312,2.7976,2.5018,2.2535,2.0622,1.9333 ,1.864  ,1.8402 ,1.8373 ,1.8267 ,1.7843 ,1.6967 ,1.5615 ,1.3849 ,1.1769 ,0.94822,0.70866,0.46621,0.22726,-0.0031286,-0.22088  ,-0.42252 ,-0.60493 ,
                                                                    42.587,40.083,37.306,34.327,31.238,28.146,25.159,22.366,19.831,17.586,15.636,13.968,12.55 ,11.348,10.326,9.4483,8.6854,8.0114,7.4056,6.8513,6.3358,5.8496,5.3861,4.9412,4.5129,4.1015,3.7094,3.3413,3.0043,2.7075,2.4608,2.2731,2.1479 ,2.0794 ,2.0502 ,2.0339 ,2.0024 ,1.9344 ,1.82   ,1.6598 ,1.4617 ,1.2361 ,0.99373,0.74386,0.49407,0.25035,0.017474  ,-0.20063  ,-0.40056 ,-0.57921 ,
                                                                    43.168,40.736,38.022,35.087,32.02 ,28.925,25.911,23.074,20.482,18.175,16.166,14.441,12.975,11.732,10.677,9.7734,8.9902,8.3007,7.6827,7.1188,6.5954,6.1024,5.6328,5.1822,4.7486,4.3323,3.9362,3.5654,3.2276,2.9321,2.6892,2.5064,2.3852 ,2.3166 ,2.28   ,2.248  ,2.1937 ,2.099  ,1.9576 ,1.7727 ,1.5535 ,1.3107 ,1.0549 ,0.79479,0.53762,0.28904,0.053587  ,-0.16493  ,-0.36315 ,-0.53789 ,
                                                                    43.719,41.361,38.712,35.825,32.784,29.691,26.655,23.776,21.129,18.761,16.69 ,14.908,13.392,12.107,11.018,10.087,9.2839,8.5791,7.9498,7.3775,6.8478,6.35  ,5.8765,5.4227,4.9864,4.5681,4.1709,3.8004,3.4646,3.1735,2.9366,2.7604,2.6437 ,2.5741 ,2.5285 ,2.4788 ,2.4004 ,2.2787 ,2.1105 ,1.9015 ,1.6618 ,1.4023 ,1.1333 ,0.86313,0.5986 ,0.34509,0.10698   ,-0.112    ,-0.30846 ,-0.47914 ,
                                                                    44.242,41.959,39.376,36.541,33.53 ,30.444,27.391,24.473,21.773,19.344,17.211,15.37 ,13.802,12.473,11.348,10.389,9.5648,8.8445,8.2042,7.6245,7.0898,6.5888,6.1135,5.6588,5.2224,4.8049,4.4097,4.0426,3.7123,3.4285,3.2004,3.0326,2.9209 ,2.8497 ,2.7937 ,2.7252 ,2.622  ,2.4734 ,2.2792 ,2.0471 ,1.7877 ,1.5123 ,1.2305 ,0.95048,0.67866,0.42023,0.17944   ,-0.039992 ,-0.23461 ,-0.40103 ,
                                                                    44.737,42.529,40.014,37.234,34.259,31.184,28.118,25.165,22.414,19.925,17.728,15.827,14.204,12.83 ,11.667,10.679,9.832 ,9.0957,8.4443,7.8574,7.3186,6.8158,6.3403,5.8867,5.4528,5.0389,4.6487,4.2885,3.967 ,3.6939,3.4773,3.3197,3.2137 ,3.1404 ,3.0732 ,2.9852 ,2.8575 ,2.6827 ,2.4639 ,2.2099 ,1.9321 ,1.6416 ,1.3477 ,1.0582 ,0.77935,0.51612,0.27274   ,0.052966  ,-0.13964 ,-0.3015  ,
                                                                    45.207,43.073,40.628,37.906,34.97 ,31.912,28.837,25.853,23.053,20.504,18.244,16.281,14.602,13.178,11.975,10.956,10.085,9.3319,8.6689,8.0747,7.5322,7.0284,6.5541,6.1035,5.6741,5.2664,4.8842,4.5341,4.2248,3.9656,3.7634,3.6178,3.5181 ,3.4426 ,3.3637 ,3.2563 ,3.1049 ,2.9057 ,2.6639 ,2.3899 ,2.0953 ,1.7908 ,1.4857 ,1.1874 ,0.90198,0.63427,0.38857   ,0.16871   ,-0.021569,-0.17845 ,
                                                                    45.652,43.592,41.218,38.556,35.664,32.626,29.548,26.537,23.691,21.084,18.759,16.733,14.995,13.52 ,12.275,11.222,10.325,9.5533,8.8776,8.2756,7.7293,7.2249,6.7526,6.3064,5.8833,5.4842,5.1127,4.7757,4.4819,4.2399,4.0545,3.9226,3.8299 ,3.7522 ,3.6617 ,3.5357 ,3.3622 ,3.1406 ,2.8782 ,2.5864 ,2.2768 ,1.96   ,1.645  ,1.3389 ,1.0475 ,0.77595,0.52841   ,0.30895   ,0.12152  ,-0.029794,
                                                                    46.074,44.087,41.784,39.185,36.34 ,33.328,30.251,27.217,24.328,21.664,19.275,17.184,15.385,13.857,12.566,11.477,10.553,9.7603,9.0706,8.4599,7.9093,7.4043,6.9346,6.4935,6.0784,5.6896,5.3313,5.0102,4.7349,4.5127,4.3465,4.2296,4.1443 ,4.0645 ,3.9629 ,3.8197 ,3.6262 ,3.385  ,3.105  ,2.798  ,2.4758 ,2.1487 ,1.8253 ,1.5127 ,1.2166 ,0.94204,0.6935    ,0.47521   ,0.29139  ,0.14649  ,
                                                                    46.472,44.558,42.327,39.793,36.998,34.017,30.946,27.894,24.965,22.246,19.793,17.636,15.774,14.19 ,12.852,11.724,10.769,9.954 ,9.2485,8.6279,8.0722,7.5663,7.0991,6.6639,6.2575,5.8807,5.5376,5.2348,4.9805,4.7806,4.6352,4.5342,4.4566 ,4.3749 ,4.2631 ,4.1044 ,3.8938 ,3.6361 ,3.3417 ,3.0227 ,2.6906 ,2.3555 ,2.0259 ,1.7085 ,1.4092 ,1.133  ,0.88468   ,0.66868   ,0.4896   ,0.35226  ,
                                                                    46.85 ,45.007,42.848,40.38 ,37.639,34.693,31.633,28.568,25.603,22.831,20.315,18.09 ,16.164,14.521,13.132,11.963,10.975,10.136,9.4126,8.7805,8.2187,7.711 ,7.2461,6.8169,6.4201,6.0564,5.7299,5.4475,5.2161,5.0402,4.9168,4.832 ,4.762  ,4.6786 ,4.5577 ,4.3858 ,4.1612 ,3.8906 ,3.5853 ,3.2576 ,2.9188 ,2.5786 ,2.2452 ,1.9253 ,1.6247 ,1.3487 ,1.1023    ,0.89016   ,0.71735  ,0.58911  ,
                                                                    47.206,45.434,43.347,40.946,38.263,35.355,32.312,29.238,26.242,23.419,20.841,18.549,16.556,14.852,13.41 ,12.196,11.174,10.307,9.5641,8.9189,8.3495,7.8392,7.3761,6.9527,6.5657,6.2159,5.9075,5.6467,5.4399,5.2891,5.188 ,5.1192,5.0563 ,4.9713 ,4.8426 ,4.66   ,4.4245 ,4.1447 ,3.8323 ,3.4996 ,3.1575 ,2.8153 ,2.481  ,2.1613 ,1.862  ,1.5885 ,1.3461    ,1.1399    ,0.9754   ,0.85825  ,
                                                                    47.543,45.84 ,43.824,41.492,38.868,36.004,32.983,29.905,26.881,24.011,21.372,19.013,16.952,15.185,13.687,12.426,11.366,10.469,9.7047,9.0446,8.4661,7.9521,7.49  ,7.0721,6.695 ,6.3596,6.0699,5.832 ,5.6506,5.5253,5.4461,5.3923,5.3359 ,5.2494 ,5.1141 ,4.9232 ,4.6802 ,4.3948 ,4.0792 ,3.7451 ,3.4033 ,3.0626 ,2.7307 ,2.4142 ,2.1191 ,1.851  ,1.6153    ,1.4177    ,1.264    ,1.1604   ,
                                                                    47.861,46.225,44.28 ,42.018,39.456,36.639,33.644,30.568,27.522,24.607,21.909,19.483,17.354,15.521,13.965,12.655,11.553,10.625,9.8361,9.1591,8.57  ,8.051 ,7.589 ,7.176 ,6.8088,6.488 ,6.2176,6.003 ,5.8474,5.7472,5.6891,5.6489,5.5977 ,5.5095 ,5.3689 ,5.1722 ,4.9248 ,4.6374 ,4.3221 ,3.9905 ,3.6526 ,3.3169 ,2.991  ,2.6812 ,2.3937 ,2.1342 ,1.9084    ,1.7225    ,1.5827   ,1.4958   ,
                                                                    48.161,46.591,44.716,42.524,40.026,37.26 ,34.295,31.227,28.162,25.208,22.454,19.961,17.762,15.863,14.246,12.883,11.738,10.775,9.9601,9.2642,8.6629,8.1375,7.6746,7.266 ,6.9083,6.6021,6.3514,6.1603,6.0301,5.9542,5.9153,5.8868,5.8395 ,5.7495 ,5.6048 ,5.4043 ,5.1554 ,4.8693 ,4.5579 ,4.2322 ,3.9017 ,3.5748 ,3.2584 ,2.959  ,2.6827 ,2.4355 ,2.2234    ,2.0529    ,1.9307   ,1.864    ,
                                                                    48.443,46.938,45.132,43.01 ,40.578,37.865,34.937,31.881,28.804,25.814,23.005,20.447,18.179,16.21 ,14.531,13.113,11.922,10.922,10.078,9.3616,8.7463,8.2132,7.7485,7.3436,6.995 ,6.7035,6.4723,6.3046,6.1992,6.1459,6.124 ,6.105 ,6.0599 ,5.9678 ,5.8198 ,5.6176 ,5.3698 ,5.0878 ,4.7834 ,4.4669 ,4.1473 ,3.8325 ,3.5293 ,3.2441 ,2.9829 ,2.752  ,2.5576    ,2.4066    ,2.3061   ,2.2639};
        }
        case 11:
        {
            static const std::vector<double> arrayCurvatureMin   = {18.186,15.925,13.891,12.07 ,10.447,9.0065,7.7317,6.6075,5.6189,4.7519,3.9933,3.3313,2.7549,2.2543,1.8208,1.4466,1.1247,0.84902,0.61413,0.41511,0.24754,0.10736,-0.0091509,-0.10541 ,-0.18449  ,-0.24911,-0.30154 ,-0.34359 ,-0.37663 ,-0.40154 ,-0.41879 ,-0.42846 ,-0.43029 ,-0.42368 ,-0.40774 ,-0.3813  ,-0.34295 ,-0.29109 ,-0.22404,-0.1402  ,-0.03827,0.082378,0.22138,0.37694,0.54559,0.72202,0.89935,1.0697,1.2253,1.3593,
                                                                    18.79 ,16.47 ,14.38 ,12.507,10.836,9.3508,8.0355,6.8745,5.8527,4.9558,4.1703,3.4841,2.8859,2.3657,1.9143,1.5238,1.1871,0.89774,0.65037,0.44004,0.26244,0.11366,-0.0098515,-0.11134 ,-0.1938   ,-0.25995,-0.3122  ,-0.35257 ,-0.38266 ,-0.40358 ,-0.41598 ,-0.42003 ,-0.41548 ,-0.40166 ,-0.37757 ,-0.34189 ,-0.29308 ,-0.22949 ,-0.14951,-0.051799,0.064391,0.19884 ,0.34999,0.51459,0.68761,0.86241,1.0313 ,1.1865,1.3211,1.4306,
                                                                    19.41 ,17.03 ,14.884,12.958,11.237,9.7064,8.3497,7.1511,6.0954,5.1679,4.3551,3.6443,3.0242,2.4842,2.0151,1.6084,1.257 ,0.95422,0.69447,0.47277,0.28477,0.12667,-0.0049341,-0.11309 ,-0.20059  ,-0.26997,-0.32357 ,-0.36342 ,-0.39123 ,-0.40831 ,-0.41549 ,-0.41312 ,-0.40101 ,-0.37856 ,-0.3447  ,-0.2981  ,-0.23717 ,-0.16034 ,-0.06625,0.045949 ,0.1762  ,0.32314 ,0.4838 ,0.65343,0.82564,0.99291,1.1476 ,1.2828,1.3938,1.4787,
                                                                    20.045,17.605,15.401,13.421,11.65 ,10.073,8.6742,7.4372,6.3468,5.3881,4.5474,3.8117,3.1692,2.6093,2.1223,1.6996,1.3337,1.0177 ,0.74587,0.51304,0.31478,0.14726,0.0071513 ,-0.10846 ,-0.20213  ,-0.27618,-0.33267 ,-0.37348 ,-0.40021 ,-0.4142  ,-0.41637 ,-0.40721 ,-0.3867  ,-0.35435 ,-0.3092  ,-0.25    ,-0.17535 ,-0.083943,0.025093,0.15182  ,0.2951  ,0.45221 ,0.61868,0.78841,0.95413,1.1082 ,1.244  ,1.3566,1.4438,1.5064,
                                                                    20.696,18.195,15.933,13.898,12.076,10.452,9.0092,7.7329,6.607 ,5.6164,4.7471,3.9858,3.3207,2.7405,2.2355,1.7967,1.4164,1.0874 ,0.80377,0.56017,0.35202,0.17539,0.026927  ,-0.096234,-0.19652  ,-0.27602,-0.33652 ,-0.37959 ,-0.40659 ,-0.41862 ,-0.41653 ,-0.40082 ,-0.37159 ,-0.32851 ,-0.27088 ,-0.19772 ,-0.10805 ,-0.001127,0.1231  ,0.26362  ,0.41792 ,0.58183 ,0.74951,0.91395,1.0677 ,1.2042 ,1.3184 ,1.4081,1.4738,1.5185,
                                                                    21.363,18.8  ,16.479,14.388,12.514,10.842,9.355 ,8.0384,6.8761,5.8528,4.9541,4.1668,3.4784,2.8776,2.3543,1.8993,1.5044,1.1626 ,0.86739,0.61337,0.39577,0.21051,0.054157  ,-0.076179,-0.18288  ,-0.26789,-0.3328  ,-0.37886 ,-0.40708 ,-0.41823 ,-0.41287 ,-0.39133 ,-0.35367 ,-0.29971 ,-0.22906 ,-0.14122 ,-0.035901,0.086715 ,0.22551 ,0.37805  ,0.54034 ,0.70679 ,0.87062,1.0246 ,1.1621 ,1.2783 ,1.3707 ,1.4398,1.4883,1.5204,
                                                                    22.045,19.421,17.04 ,14.892,12.965,11.243,9.7116,8.3538,7.1542,6.0973,5.1686,4.3545,3.6422,3.0204,2.4784,2.0069,1.5975,1.2427 ,0.93612,0.67194,0.44526,0.25188,0.088225  ,-0.048656,-0.16117  ,-0.2512 ,-0.3202  ,-0.3692  ,-0.39891 ,-0.40974 ,-0.40186 ,-0.3753  ,-0.32995 ,-0.26565 ,-0.18232 ,-0.080085,0.040413 ,0.17772  ,0.32922 ,0.49085  ,0.65711 ,0.82131 ,0.97638,1.1158 ,1.2346 ,1.3304 ,1.4034 ,1.4561,1.4925,1.5169,
                                                                    22.744,20.057,17.615,15.41 ,13.429,11.657,10.079,8.6793,7.4414,6.3501,5.3906,4.549 ,3.8123,3.1688,2.6077,2.1194,1.6952,1.3275 ,1.0095 ,0.73532,0.49989,0.29882,0.12845   ,-0.014272,-0.13179  ,-0.22601,-0.29829 ,-0.34956 ,-0.38028 ,-0.39059 ,-0.38034 ,-0.34924 ,-0.29695 ,-0.22331 ,-0.12852 ,-0.013384,0.12036  ,0.2698   ,0.43062 ,0.59713  ,0.76256 ,0.91975 ,1.0621 ,1.1847 ,1.2848 ,1.3626 ,1.4203 ,1.4617,1.491 ,1.5118,
                                                                    23.459,20.708,18.206,15.943,13.906,12.083,10.458,9.0152,7.7381,6.6115,5.6203,4.7504,3.9886,3.3229,2.7423,2.2368,1.7975,1.4167 ,1.0872 ,0.80317,0.5592 ,0.35083,0.17427   ,0.026392 ,-0.095299 ,-0.19272,-0.26722 ,-0.31963 ,-0.35033 ,-0.35926 ,-0.3461  ,-0.31035 ,-0.25159 ,-0.16973 ,-0.065313,0.060098 ,0.20364  ,0.36092  ,0.52607 ,0.6921   ,0.85157 ,0.99759 ,1.1249 ,1.2306 ,1.3144 ,1.3783 ,1.4258 ,1.4608,1.4869,1.5071,
                                                                    24.191,21.376,18.812,16.489,14.397,12.522,10.849,9.3616,8.0444,6.8816,5.8578,4.9589,4.1713,3.4828,2.8821,2.359 ,1.9044,1.5102 ,1.1693 ,0.87531,0.62296,0.40759,0.22531   ,0.072921 ,-0.052118 ,-0.15174,-0.22729 ,-0.27957 ,-0.30889 ,-0.31518 ,-0.29805 ,-0.25706 ,-0.19192 ,-0.10285 ,0.0090083,0.14125  ,0.28989  ,0.4493   ,0.61257 ,0.77209  ,0.92066 ,1.0525  ,1.1643 ,1.2551 ,1.3265 ,1.3815 ,1.4236 ,1.4562,1.4821,1.5032,
                                                                    24.944,22.062,19.434,17.051,14.902,12.974,11.251,9.7189,8.3606,7.1605,6.1034,5.1747,4.3606,3.6487,3.0274,2.4862,2.0159,1.6081 ,1.2555 ,0.95169,0.69109,0.46899,0.28141   ,0.12511  ,-0.0024902,-0.10334,-0.17877 ,-0.2296  ,-0.25614 ,-0.25838 ,-0.23609 ,-0.18907 ,-0.11751 ,-0.022262,0.094662 ,0.22983  ,0.37822  ,0.53343  ,0.6882  ,0.83541  ,0.96909 ,1.0853  ,1.1826 ,1.2618 ,1.3251 ,1.3756 ,1.4161 ,1.4491,1.4765,1.4994,
                                                                    25.727,22.77 ,20.074,17.629,15.421,13.439,11.666,10.087,8.6868,7.4487,6.3573,5.3979,4.5567,3.8207,3.1782,2.6185,2.1322,1.7105 ,1.3461 ,1.0324 ,0.76367,0.53508,0.34262   ,0.18298  ,0.05356   ,-0.04758,-0.1218  ,-0.16992 ,-0.19237 ,-0.1893  ,-0.1608  ,-0.10721 ,-0.029431,0.070673 ,0.18999  ,0.32398  ,0.46681  ,0.61184  ,0.75243 ,0.88288  ,0.99926 ,1.0998  ,1.1846 ,1.2553 ,1.3141 ,1.3632 ,1.4046 ,1.4396,1.4693,1.4945,
                                                                    26.565,23.511,20.738,18.225,15.956,13.917,12.093,10.467,9.0235,7.7463,6.6197,5.6289,4.7597,3.999 ,3.3347,2.7561,2.2534,1.8176 ,1.4413 ,1.1177 ,0.84088,0.60607,0.40912   ,0.2467   ,0.1162    ,0.015652,-0.056328,-0.10066 ,-0.11792 ,-0.1086  ,-0.073304,-0.013134,0.069992 ,0.17303  ,0.29165  ,0.42036  ,0.55296  ,0.68324  ,0.80596 ,0.9175   ,1.0162  ,1.1021  ,1.1764 ,1.2409 ,1.297  ,1.3462 ,1.3892 ,1.4267,1.4591,1.4868,
                                                                    27.529,24.32 ,21.44 ,18.846,16.51 ,14.412,12.534,10.859,9.371 ,8.0536,6.891 ,5.8678,4.9699,4.1838,3.4973,2.8993,2.3797,1.9297 ,1.5413 ,1.2077 ,0.92302,0.68225,0.48124   ,0.31661  ,0.18574   ,0.08661 ,0.017763 ,-0.021882,-0.033217,-0.017204,0.024808 ,0.090745 ,0.17751  ,0.28082  ,0.39531  ,0.51497  ,0.63383  ,0.74684  ,0.85058 ,0.94356  ,1.0259  ,1.099   ,1.1644 ,1.2235 ,1.2773 ,1.3262 ,1.3704 ,1.4099,1.4446,1.4748,
                                                                    28.762,25.286,22.227,19.514,17.092,14.926,12.99 ,11.264,9.7299,8.3712,7.1714,6.1151,5.1876,4.3755,3.6661,3.0482,2.5115,2.0469 ,1.6464 ,1.3029 ,1.0105 ,0.76404,0.5594    ,0.39314  ,0.26256   ,0.16559 ,0.1006   ,0.066218 ,0.061102 ,0.083604 ,0.13145  ,0.20143  ,0.28926  ,0.38961  ,0.49654  ,0.60414  ,0.70739  ,0.80281  ,0.88883 ,0.96556  ,1.0342  ,1.0965  ,1.1539 ,1.2075 ,1.2578 ,1.305  ,1.3488 ,1.389 ,1.4253,1.4573,
                                                                    30.361,26.56 ,23.204,20.285,17.73 ,15.473,13.467,11.684,10.102,8.6996,7.4614,6.3709,5.4132,4.5743,3.8415,3.2033,2.6491,2.1697 ,1.7569 ,1.4036 ,1.1037 ,0.8519 ,0.64409   ,0.47674  ,0.34708   ,0.25285 ,0.19217  ,0.16323  ,0.16404  ,0.19208  ,0.24401  ,0.31549  ,0.4012   ,0.49525  ,0.59173  ,0.68555  ,0.77309  ,0.8525   ,0.92358 ,0.98726  ,1.045   ,1.0983  ,1.1483 ,1.1959 ,1.2414 ,1.2849 ,1.3263 ,1.3653,1.4015,1.4345,
                                                                    32.06 ,28.176,24.525,21.286,18.496,16.087,13.982,12.128,10.489,9.0407,7.762 ,6.636 ,5.647 ,4.7806,4.0238,3.3648,2.7929,2.2985 ,1.8734 ,1.5104 ,1.2031 ,0.94634,0.7358    ,0.56787  ,0.43961   ,0.34849 ,0.29221  ,0.26836  ,0.27413  ,0.30597  ,0.35948  ,0.42933  ,0.50964  ,0.59453  ,0.67882  ,0.7587   ,0.83199  ,0.89808  ,0.95743 ,1.0111   ,1.0603  ,1.106   ,1.149  ,1.1901 ,1.2296 ,1.2678 ,1.3049 ,1.3407,1.375 ,1.4075,
                                                                    33.471,29.775,26.131,22.652,19.533,16.858,14.58 ,12.616,10.903,9.3989,8.0754,6.9114,5.8897,4.9949,4.2135,3.5332,2.9432,2.4337 ,1.9963 ,1.6236 ,1.3092 ,1.0479 ,0.83501   ,0.66689  ,0.54031   ,0.45233 ,0.40007  ,0.38035  ,0.38937  ,0.4226   ,0.47467  ,0.53972  ,0.61185  ,0.68582  ,0.75762  ,0.82476  ,0.88617  ,0.94182  ,0.99224 ,1.0382   ,1.0803  ,1.1193  ,1.1557 ,1.1902 ,1.2232 ,1.2552 ,1.2866 ,1.3176,1.3483,1.3784,
                                                                    34.587,31.047,27.607,24.217,20.934,17.94 ,15.37 ,13.208,11.371,9.787 ,8.4073,7.2   ,6.1428,5.218 ,4.411 ,3.709 ,3.1005,2.5758 ,2.1261 ,1.7439 ,1.4227 ,1.157  ,0.94208   ,0.77398  ,0.64905   ,0.56381 ,0.51462  ,0.49739  ,0.50737  ,0.53912  ,0.58674  ,0.64432  ,0.70655  ,0.76922  ,0.82953  ,0.88595  ,0.93794  ,0.9855   ,1.0289  ,1.0684   ,1.1043  ,1.1371  ,1.1673 ,1.1953 ,1.2218 ,1.2474 ,1.2727 ,1.298 ,1.3238,1.35  ,
                                                                    35.565,32.064,28.746,25.556,22.426,19.36 ,16.504,14.03 ,11.966,10.241,8.7746,7.5094,6.4099,5.4519,4.6176,3.8928,3.2655,2.7254 ,2.2633 ,1.8717 ,1.5439 ,1.274  ,1.0572    ,0.88904  ,0.76533   ,0.68189 ,0.63423  ,0.61731  ,0.62556  ,0.653    ,0.69367  ,0.74216  ,0.79409  ,0.8463   ,0.89683  ,0.9446   ,0.98908  ,1.03     ,1.0673  ,1.1009   ,1.131   ,1.1579  ,1.1821 ,1.2041 ,1.2245 ,1.2441 ,1.2634 ,1.2829,1.3033,1.3247,
                                                                    36.498,32.974,29.673,26.574,23.624,20.751,17.916,15.216,12.834,10.853,9.2242,7.8618,6.7011,5.701 ,4.8355,4.0861,3.439 ,2.883  ,2.4086 ,2.0076 ,1.6731 ,1.3992 ,1.1803    ,1.0116   ,0.88822   ,0.80511 ,0.75691  ,0.73781  ,0.7417   ,0.76253  ,0.79477  ,0.8338   ,0.87613  ,0.91934  ,0.96181  ,1.0024   ,1.0405   ,1.0753   ,1.1067  ,1.1346   ,1.159   ,1.1802  ,1.1988 ,1.2153 ,1.2303 ,1.2444 ,1.2583 ,1.2725,1.2875,1.3039,
                                                                    37.427,33.855,30.52 ,27.419,24.532,21.812,19.187,16.59 ,14.063,11.777,9.8664,8.3166,7.0451,5.9789,5.0708,4.2918,3.6226,3.0497 ,2.5624 ,2.1519 ,1.8107 ,1.5324 ,1.311     ,1.1409   ,1.0164    ,0.93165 ,0.8806   ,0.85691  ,0.85436  ,0.86726  ,0.89073  ,0.92092  ,0.95491  ,0.99053  ,1.0261   ,1.0603   ,1.0921   ,1.1209   ,1.1464  ,1.1684   ,1.1871  ,1.2029  ,1.2163 ,1.2278 ,1.2379 ,1.2473 ,1.2565 ,1.2661,1.2764,1.288 ,
                                                                    38.363,34.737,31.349,28.207,25.304,22.623,20.121,17.729,15.369,13.032,10.853,9.0038,7.5171,6.3227,5.3411,4.5178,3.8201,3.2274 ,2.7258 ,2.3051 ,1.9567 ,1.6734 ,1.4485    ,1.2755   ,1.1481    ,1.0597  ,1.0036   ,0.97345  ,0.96331  ,0.96801  ,0.98329  ,1.0057   ,1.0325   ,1.0613   ,1.0903   ,1.1181   ,1.1435   ,1.166    ,1.1852  ,1.2013   ,1.2145  ,1.2251  ,1.2337 ,1.2407 ,1.2467 ,1.2521 ,1.2573 ,1.2629,1.2693,1.2768,
                                                                    39.311,35.628,32.183,28.985,26.035,23.328,20.845,18.548,16.376,14.245,12.109,10.054,8.2632,6.8253,5.694 ,4.7869,4.0418,3.4207 ,2.901  ,2.468  ,2.1111 ,1.8216 ,1.5917    ,1.4142   ,1.2819    ,1.1877  ,1.125    ,1.0874   ,1.0695   ,1.0665   ,1.0745   ,1.09     ,1.11     ,1.1321   ,1.1543   ,1.1752   ,1.1938   ,1.2096   ,1.2225  ,1.2327   ,1.2405  ,1.2463  ,1.2505 ,1.2535 ,1.2559 ,1.2578 ,1.2598 ,1.2622,1.2653,1.2694,
                                                                    40.273,36.534,33.028,29.77 ,26.763,24.006,21.491,19.197,17.094,15.126,13.212,11.281,9.3677,7.642 ,6.242 ,5.16  ,4.317 ,3.643  ,3.0936 ,2.6429 ,2.2743 ,1.9763 ,1.7395    ,1.5555   ,1.4164    ,1.3152  ,1.245    ,1.1999   ,1.1747   ,1.1648   ,1.166    ,1.1747   ,1.1878   ,1.2027   ,1.2175   ,1.2309   ,1.2422   ,1.2512   ,1.2577  ,1.2622   ,1.2649  ,1.2662  ,1.2664 ,1.2659 ,1.2651 ,1.2641 ,1.2634 ,1.2632,1.2636,1.265 ,
                                                                    41.249,37.453,33.889,30.569,27.501,24.686,22.119,19.789,17.678,15.757,13.977,12.265,10.538,8.7837,7.1359,5.7681,4.722 ,3.9318 ,3.3208 ,2.8368 ,2.4488 ,2.1381 ,1.8914    ,1.6988   ,1.5516    ,1.4425  ,1.3648   ,1.3127   ,1.2808   ,1.2643   ,1.2586   ,1.26     ,1.2653   ,1.2721   ,1.2788   ,1.2843   ,1.2881   ,1.2902   ,1.2906  ,1.2896   ,1.2875  ,1.2846  ,1.2812 ,1.2776 ,1.274  ,1.2706 ,1.2677 ,1.2653,1.2636,1.2628,
                                                                    42.237,38.388,34.764,31.383,28.253,25.377,22.752,20.371,18.221,16.284,14.534,12.928,11.403,9.8729,8.2901,6.7386,5.4037,4.381  ,3.6303 ,3.0716 ,2.6439 ,2.3103 ,2.0486    ,1.8448   ,1.6884    ,1.5713  ,1.4865   ,1.4277   ,1.3891   ,1.3654   ,1.3521   ,1.3451   ,1.3417   ,1.3396   ,1.3376   ,1.3349   ,1.3312   ,1.3265   ,1.321   ,1.3148   ,1.3082  ,1.3015  ,1.2949 ,1.2885 ,1.2826 ,1.2771 ,1.2722 ,1.2681,1.2647,1.2622,
                                                                    43.239,39.337,35.656,32.213,29.021,26.082,23.398,20.96 ,18.76 ,16.784,15.013,13.424,11.98 ,10.626,9.2812,7.8755,6.4401,5.1451 ,4.135  ,3.4082 ,2.8883 ,2.506  ,2.2174    ,1.9971   ,1.8296    ,1.704   ,1.6118   ,1.5459   ,1.4998   ,1.4679   ,1.4455   ,1.429    ,1.4158   ,1.4042   ,1.3932   ,1.3823   ,1.3712   ,1.3601   ,1.3489  ,1.3379   ,1.3273  ,1.3172  ,1.3076 ,1.2988 ,1.2907 ,1.2834 ,1.2769 ,1.2713,1.2665,1.2626,
                                                                    44.253,40.3  ,36.562,33.06 ,29.805,26.805,24.059,21.563,19.308,17.284,15.474,13.862,12.425,11.131,9.9331,8.7597,7.529 ,6.2262 ,4.9835 ,3.9785 ,3.2596 ,2.7641 ,2.4166    ,2.1655   ,1.9805    ,1.8436  ,1.7424   ,1.6678   ,1.6124   ,1.5705   ,1.5376   ,1.5105   ,1.4869   ,1.4654   ,1.4454   ,1.4264   ,1.4083   ,1.391    ,1.3746  ,1.3592   ,1.3449  ,1.3316  ,1.3194 ,1.3083 ,1.2984 ,1.2895 ,1.2816 ,1.2747,1.2688,1.2638,
                                                                    45.28 ,41.278,37.485,33.923,30.606,27.544,24.736,22.181,19.871,17.794,15.939,14.289,12.827,11.533,10.378,9.323 ,8.3044,7.2399 ,6.08   ,4.9065 ,3.9063 ,3.182  ,2.6977    ,2.375    ,2.1532    ,1.9954  ,1.88     ,1.7932   ,1.726    ,1.6722   ,1.6274   ,1.5887   ,1.5543   ,1.523    ,1.4941   ,1.4674   ,1.4425   ,1.4195   ,1.3983  ,1.3788   ,1.361   ,1.3449  ,1.3304 ,1.3173 ,1.3057 ,1.2953 ,1.2862 ,1.2782,1.2713,1.2655,
                                                                    46.319,42.269,38.422,34.802,31.424,28.3  ,25.431,22.816,20.449,18.319,16.415,14.722,13.224,11.905,10.744,9.7177,8.7912,7.9109 ,6.9997 ,5.9873 ,4.9031 ,3.9167 ,3.1779    ,2.6902   ,2.3785    ,2.173   ,2.0296   ,1.9233   ,1.8402   ,1.7721   ,1.7141   ,1.6633   ,1.6179   ,1.5768   ,1.5394   ,1.5053   ,1.4742   ,1.4458   ,1.4202  ,1.3969   ,1.3761  ,1.3573  ,1.3406 ,1.3258 ,1.3126 ,1.301  ,1.2908 ,1.2818,1.2741,1.2674,
                                                                    47.369,43.273,39.374,35.697,32.259,29.073,26.142,23.467,21.043,18.859,16.904,15.166,13.627,12.274,11.088,10.052,9.1433,8.3336 ,7.5777 ,6.8071 ,5.9419 ,4.9655 ,4.0068    ,3.2457   ,2.735     ,2.4142  ,2.2076   ,2.0646   ,1.9572   ,1.8708   ,1.7979   ,1.7343   ,1.6778   ,1.6272   ,1.5816   ,1.5406   ,1.5036   ,1.4703   ,1.4405  ,1.4139   ,1.3901  ,1.369   ,1.3503 ,1.3338 ,1.3192 ,1.3064 ,1.2952 ,1.2854,1.2769,1.2696,
                                                                    48.429,44.291,40.341,36.608,33.111,29.863,26.871,24.136,21.653,19.414,17.408,15.622,14.041,12.649,11.431,10.37 ,9.4508,8.6529 ,7.9517 ,7.3096 ,6.6662 ,5.9402 ,5.0795    ,4.1608   ,3.37      ,2.8145  ,2.4629   ,2.2394   ,2.0867   ,1.9725   ,1.8804   ,1.8025   ,1.7346   ,1.6746   ,1.6212   ,1.5735   ,1.531    ,1.4932   ,1.4596  ,1.4298   ,1.4034  ,1.3801  ,1.3595 ,1.3415 ,1.3256 ,1.3117 ,1.2996 ,1.289 ,1.2799,1.2719,
                                                                    49.501,45.321,41.322,37.534,33.979,30.671,27.617,24.821,22.279,19.985,17.926,16.091,14.465,13.034,11.781,10.691,9.7506,8.9428 ,8.2504 ,7.6506 ,7.1099 ,6.575  ,5.9686    ,5.2182   ,4.3484    ,3.5272  ,2.9108   ,2.511    ,2.2592   ,2.0911   ,1.9682   ,1.871    ,1.7898   ,1.7198   ,1.6586   ,1.6046   ,1.557    ,1.5148   ,1.4777  ,1.4449   ,1.416   ,1.3906  ,1.3684 ,1.3489 ,1.3318 ,1.3169 ,1.3039 ,1.2926,1.2828,1.2743,
                                                                    50.582,46.363,42.318,38.477,34.864,31.495,28.38 ,25.523,22.922,20.57 ,18.458,16.573,14.901,13.429,12.141,11.022,10.058,9.2333 ,8.5329 ,7.9379 ,7.4262 ,6.9674 ,6.5157    ,6.0027   ,5.3492    ,4.5374  ,3.6967   ,3.0144   ,2.5551   ,2.2669   ,2.0801   ,1.9485   ,1.8475   ,1.765    ,1.695    ,1.6346   ,1.5818   ,1.5355   ,1.4949  ,1.4594   ,1.4282  ,1.4009  ,1.377  ,1.3561 ,1.3379 ,1.3221 ,1.3083 ,1.2963,1.2858,1.2768,
                                                                    51.672,47.416,43.326,39.434,35.764,32.336,29.16 ,26.242,23.58 ,21.17 ,19.003,17.068,15.35 ,13.837,12.514,11.367,10.379,9.5363 ,8.8211 ,8.2158 ,7.7014 ,7.2568 ,6.8561    ,6.4626   ,6.0191    ,5.4493  ,4.7052   ,3.866    ,3.1245   ,2.5994   ,2.2683   ,2.0599   ,1.9192   ,1.8154   ,1.733    ,1.6645   ,1.6061   ,1.5556   ,1.5117  ,1.4734   ,1.44    ,1.4108  ,1.3854 ,1.3632 ,1.3439 ,1.3272 ,1.3126 ,1.2999,1.2889,1.2794,
                                                                    52.77 ,48.481,44.348,40.405,36.68 ,33.193,29.956,26.975,24.253,21.785,19.563,17.576,15.813,14.261,12.905,11.73 ,10.719,9.856  ,9.1219 ,8.4983 ,7.9671 ,7.5108 ,7.1117    ,6.7497   ,6.3962    ,6.0056  ,5.509    ,4.8407   ,4.0287   ,3.2443   ,2.652    ,2.2714   ,2.0371   ,1.8859   ,1.7793   ,1.6977   ,1.6316   ,1.576    ,1.5284  ,1.4873   ,1.4517  ,1.4207  ,1.3937 ,1.3703 ,1.3499 ,1.3322 ,1.3169 ,1.3035,1.292 ,1.282 ,
                                                                    53.877,49.556,45.382,41.391,37.611,34.065,30.766,27.724,24.941,22.414,20.137,18.1  ,16.294,14.704,13.316,12.114,11.079,10.193 ,9.4357 ,8.788  ,8.2325 ,7.7532 ,7.336     ,6.9676   ,6.6326    ,6.3094  ,5.9616   ,5.5298   ,4.9415   ,4.1806   ,3.3765   ,2.7206   ,2.2841   ,2.0181   ,1.8533   ,1.743    ,1.6621   ,1.5984   ,1.5459  ,1.5015   ,1.4634  ,1.4305  ,1.4021 ,1.3774 ,1.3559 ,1.3374 ,1.3212 ,1.3073,1.2951,1.2847,
                                                                    54.99 ,50.64 ,46.427,42.389,38.556,34.951,31.592,28.488,25.644,23.059,20.728,18.643,16.795,15.169,13.751,12.52 ,11.459,10.546 ,9.7605 ,9.0836 ,8.4985 ,7.9907 ,7.5477    ,7.1588   ,6.8139    ,6.5011  ,6.2038   ,5.8931   ,5.5189   ,5.0096   ,4.318    ,3.5209   ,2.8113   ,2.3134   ,2.0083   ,1.8255   ,1.7091   ,1.628    ,1.5665  ,1.517    ,1.4757  ,1.4407  ,1.4105 ,1.3845 ,1.362  ,1.3425 ,1.3257 ,1.311 ,1.2984,1.2874,
                                                                    56.109,51.732,47.483,43.399,39.513,35.851,32.432,29.267,26.363,23.721,21.339,19.209,17.321,15.661,14.21 ,12.949,11.857,10.912 ,10.093 ,9.3826 ,8.7644 ,8.2249 ,7.7528    ,7.3386   ,6.9738    ,6.6502  ,6.3587   ,6.0856   ,5.808    ,5.4845   ,5.0495   ,4.4377   ,3.6738   ,2.9268   ,2.365    ,2.0126   ,1.8054   ,1.6798   ,1.5969  ,1.5368   ,1.49    ,1.4517  ,1.4194 ,1.3919 ,1.3682 ,1.3478 ,1.3302 ,1.3149,1.3017,1.2902,
                                                                    57.233,52.831,48.547,44.42 ,40.483,36.765,33.286,30.062,27.1  ,24.405,21.975,19.802,17.875,16.179,14.694,13.398,12.269,11.287 ,10.431 ,9.6828 ,9.0288 ,8.4559 ,7.9533    ,7.5121   ,7.1244    ,6.7831  ,6.4813   ,6.211    ,5.9613   ,5.7137   ,5.4345   ,5.0663   ,4.5373   ,3.8289   ,3.0664   ,2.4431   ,2.0351   ,1.7962   ,1.6568  ,1.5699   ,1.5101  ,1.4653  ,1.4295 ,1.3999 ,1.3748 ,1.3533 ,1.3348 ,1.3188,1.3051,1.2931,
                                                                    58.36 ,53.936,49.619,45.45 ,41.464,37.691,34.156,30.875,27.86 ,25.115,22.639,20.424,18.458,16.724,15.199,13.863,12.693,11.668 ,10.77  ,9.9825 ,9.291  ,8.6837 ,8.1502    ,7.6817   ,7.2704    ,6.9093  ,6.5921   ,6.3126   ,6.0638   ,5.8366   ,5.616    ,5.375    ,5.0652   ,4.616    ,3.9788   ,3.226    ,2.5499   ,2.0798   ,1.8004  ,1.6417   ,1.5477  ,1.4867  ,1.4431 ,1.4094 ,1.3821 ,1.3592 ,1.3397 ,1.323 ,1.3086,1.2961,
                                                                    59.488,55.045,50.698,46.489,42.456,38.632,35.044,31.711,28.646,25.854,23.334,21.078,19.07 ,17.292,15.723,14.341,13.124,12.053 ,11.11  ,10.281 ,9.5508 ,8.9087 ,8.3443    ,7.8486   ,7.4137    ,7.0324  ,6.6984   ,6.4057   ,6.1487   ,5.9211   ,5.7152   ,5.519    ,5.3108   ,5.0508   ,4.6746   ,4.1167   ,3.3982   ,2.6853   ,2.1499  ,1.821    ,1.6361  ,1.5312  ,1.467  ,1.4236 ,1.3915 ,1.366  ,1.345  ,1.3273,1.3122,1.2992,
                                                                    60.617,56.157,51.782,47.537,43.461,39.591,35.954,32.574,29.463,26.627,24.063,21.762,19.707,17.881,16.261,14.827,13.559,12.439 ,11.45  ,10.577 ,9.8084 ,9.1315 ,8.5363    ,8.0137   ,7.5555    ,7.1542  ,6.8031   ,6.4961   ,6.2278   ,5.9928   ,5.7857   ,5.5998   ,5.4254   ,5.2451   ,5.0269   ,4.7149   ,4.2377   ,3.5744   ,2.8462  ,2.2479   ,1.861   ,1.6419  ,1.5213 ,1.4514 ,1.4069 ,1.3756 ,1.3517 ,1.3323,1.3162,1.3026,
                                                                    61.745,57.271,52.872,48.595,44.481,40.569,36.891,33.467,30.313,27.433,24.824,22.473,20.367,18.485,16.81 ,15.32 ,13.997,12.825 ,11.789 ,10.872 ,10.064 ,9.3528 ,8.7272    ,8.1781   ,7.6968    ,7.2756  ,6.9074   ,6.586    ,6.3055   ,6.0608   ,5.847    ,5.6593   ,5.4918   ,5.3367   ,5.1803   ,4.9968   ,4.7396   ,4.3391   ,3.7454  ,3.0267   ,2.3747  ,1.9233  ,1.6612 ,1.519  ,1.4403 ,1.3932 ,1.3619 ,1.339 ,1.321 ,1.3062,
                                                                    62.871,58.387,53.969,49.665,45.52 ,41.574,37.858,34.395,31.199,28.274,25.614,23.209,21.044,19.102,17.365,15.816,14.437,13.212 ,12.127 ,11.167 ,10.32  ,9.5737 ,8.9179    ,8.3425   ,7.8384    ,7.3975  ,7.0123   ,6.6762   ,6.3833   ,6.1282   ,5.9059   ,5.7122   ,5.5425   ,5.3919   ,5.254    ,5.1178   ,4.9631   ,4.7516   ,4.4204  ,3.9037   ,3.2189  ,2.5294  ,2.0106 ,1.6965 ,1.5258 ,1.4344 ,1.3826 ,1.3503,1.3279,1.3109,
                                                                    63.995,59.506,55.075,50.752,46.582,42.608,38.859,35.359,32.12 ,29.144,26.428,23.963,21.734,19.727,17.926,16.315,14.877,13.599 ,12.465 ,11.461 ,10.575 ,9.795  ,9.1093    ,8.5078   ,7.981     ,7.5205  ,7.1183   ,6.7676   ,6.4621   ,6.1961   ,5.9648   ,5.7636   ,5.5885   ,5.4354   ,5.3003   ,5.1775   ,5.0587   ,4.9276   ,4.7538  ,4.4827   ,4.0437  ,3.4131  ,2.7083 ,2.1252 ,1.7505 ,1.5434 ,1.4345 ,1.3755,1.3409,1.3184,
                                                                    65.118,60.631,56.195,51.86 ,47.673,43.675,39.896,36.357,33.071,30.04 ,27.263,24.731,22.434,20.358,18.49 ,16.816,15.319,13.987 ,12.804 ,11.756 ,10.832 ,10.018 ,9.3021    ,8.6746   ,8.1253    ,7.645   ,7.2259   ,6.8605   ,6.5422   ,6.2653   ,6.0246   ,5.8154   ,5.6337   ,5.4757   ,5.3379   ,5.2166   ,5.1074   ,5.0033   ,4.8918  ,4.7486   ,4.5284  ,4.1627  ,3.6003 ,2.9051 ,2.2675 ,1.8262 ,1.5738 ,1.4415,1.3723,1.3339,
                                                                    66.243,61.766,57.333,52.994,48.795,44.776,40.966,37.386,34.048,30.957,28.113,25.51 ,23.141,20.994,19.058,17.319,15.763,14.376 ,13.144 ,12.054 ,11.091 ,10.243 ,9.4973    ,8.8438   ,8.2717    ,7.7718  ,7.3355   ,6.9552   ,6.6241   ,6.3361   ,6.0857   ,5.8683   ,5.6795   ,5.5157   ,5.3734   ,5.2494   ,5.1407   ,5.0433   ,4.9519  ,4.8566   ,4.7383  ,4.5601  ,4.2602 ,3.7728 ,3.1114 ,2.4364 ,1.9263 ,1.6196,1.457 ,1.3737,
                                                                    67.374,62.916,58.493,54.156,49.949,45.91 ,42.066,38.441,35.046,31.89 ,28.974,26.298,23.853,21.634,19.628,17.824,16.209,14.768 ,13.488 ,12.354 ,11.353 ,10.471 ,9.6955    ,9.0158   ,8.4209    ,7.9011  ,7.4474   ,7.0521   ,6.7079   ,6.4086   ,6.1484   ,5.9225   ,5.7264   ,5.5563   ,5.4088   ,5.2807   ,5.1694   ,5.0719   ,4.9851  ,4.9046   ,4.8227  ,4.7245  ,4.5805 ,4.3375 ,3.9254 ,3.3179 ,2.628  ,2.0527,1.6836,1.4827};
        }
        }

        return 0;
    }
