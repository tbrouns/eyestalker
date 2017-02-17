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
    else                       { index = *itr; }
    return index;
}

// Scores

inline double calculateCertainty(double x, double midpoint)
{
    double a = certaintyAsymptoteX * midpoint;
    double k = - log((1 / certaintyAsymptoteY) - 1) / a;
    return (1 - 2 / (1 + exp(-k * (x - midpoint)))); // Logistic
}

inline double calculateScoreIntensity(double x)
{
    double a1 = 1.0, b1 = 0, c1 = 12;
    return (a1 * exp(-pow((x - b1) / c1, 2)));  // Gaussian
}

inline double calculateScoreRadius(double x)
{
    double a1 = 0.54, b1 = 1.006, c1 = 0.044, a2 = 0.52, b2 = 0.97, c2 = 0.085;
    return (a1 * exp(-pow((x - b1) / c1, 2)) + a2 * exp(-pow((x - b2) / c2, 2))); // Two-term Gaussian
}

inline double calculateScoreRadiusVar(double x)
{
    double a = 2.73, b = -0.13;
    return (exp(-a * (x - b))); // Exponential decay
}

inline double calculateScoreCurvature(double x)
{
    double a1 = 0.90, b1 = 0.0, c1 = 8.0;
    return (a1 * exp(-pow((x - b1) / c1, 2))); // Gaussian
}

inline double calculateScoreCircumference(double x)
{
    double a = 10, b = 0.3;
    return (1 / (1 + exp(-a * (x - b)))); // Logistic
}

inline double calculateScoreGradient(double x)
{
    double a1 = 0.98, b1 = 0, c1 = 6.4;
    return (a1 * exp(-pow((x - b1) / c1, 2))); // Gaussian
}

double calculateScoreTotal(const detectionVariables& mDetectionVariables, std::vector<double>& inputVector, bool USE_CERTAINTY)
{
    for (int i = 0, vSize = inputVector.size(); i < vSize; i++) // check for NaNs or Infs
    {
        double val = inputVector[i];
        if (!std::isfinite(val)) { inputVector[i] = 0; }
        else if (val < 0)        { inputVector[i] = std::abs(val); }
    }

    // Do score calculation

    double changeRadius        = inputVector[0];
    double changeCircumference = inputVector[1];
    double changeCurvature     = inputVector[2];
    double changeIntensity     = inputVector[3];
    double changeGradient      = inputVector[4];
    double varianceRadius      = inputVector[5];
    double edgeLength          = inputVector[6];

    changeIntensity = std::abs(changeIntensity);
    changeGradient  = std::abs(changeGradient);
    changeCurvature = std::abs(changeCurvature);

    double certaintyFactorPosition = 0.5 * (mDetectionVariables.certaintyPosition + 1);
    double certaintyFactorFeatures = 0.5 * (mDetectionVariables.certaintyFeatures + 1);
    double certaintyFactorAverages = 0.5 * (mDetectionVariables.certaintyAverages + 1);

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

    double factorLength = (1 / (mDetectionVariables.predictedCircumference - 2 * curvatureWindowLength)) * (edgeLength - 2 * curvatureWindowLength);
    if (factorLength < 0) { factorLength = 0; }

    double factorRadiusVar = certaintyFactorPosition  * scoreFactorRadiusVar * factorLength;
    double factorCurvature = certaintyFactorCurvature * scoreFactorCurvature * factorLength;

    // Calculate scores

    double scoreCircumference = factorCircumference * calculateScoreCircumference(changeCircumference);
    double scoreRadiusVar     = factorRadiusVar     * calculateScoreRadiusVar(varianceRadius);
    double scoreCurvature     = factorCurvature     * calculateScoreCurvature(changeCurvature);
    double scoreRadius        = factorRadius        * calculateScoreRadius   (changeRadius);
    double scoreIntensity     = factorIntensity     * calculateScoreIntensity(changeIntensity);
    double scoreGradient      = factorGradient      * calculateScoreGradient (changeGradient);

    const double norm =  factorRadius + factorRadiusVar + factorCurvature + factorCircumference + factorIntensity + factorGradient;

    double scoreTotal;
    if (norm > 0) { scoreTotal = (scoreRadius + scoreRadiusVar + scoreGradient + scoreCurvature + scoreCircumference + scoreIntensity) / norm; }
    else          { scoreTotal = 0; }
    return scoreTotal;
}

// Detection

std::vector<unsigned int> calculateIntImg(const cv::Mat& img, int imgWidth, AOIProperties searchAOI)
{
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

AOIProperties detectGlint(const cv::Mat& img, int imgWidth, AOIProperties searchAOI, AOIProperties glintAOI)
{
    int stepSize        = 2;
    int numDirections   = 4;
    int glintThreshold  = 200;

    int gradientWindowLength = glintAOI.wdth;
    int glintRadius = round(0.5 * glintAOI.wdth);

    uchar *ptr = img.data;

    std::vector<double> imageGradient(searchAOI.wdth * searchAOI.hght, 0.0);

    std::vector<int> dZ(numDirections);
    dZ[0] = -searchAOI.wdth - 1;
    dZ[1] = -searchAOI.wdth + 1;
    dZ[2] =  searchAOI.wdth + 1;
    dZ[3] =  searchAOI.wdth - 1;

    for (int y = gradientWindowLength; y < searchAOI.hght - gradientWindowLength; y = y + stepSize)
    {
        for (int x = gradientWindowLength; x < searchAOI.wdth - gradientWindowLength; x = x + stepSize)
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

AOIProperties detectPupilApprox(const std::vector<unsigned int>& I, AOIProperties searchAOI, AOIProperties haarAOI, AOIProperties glintAOI)
{
    int stepSize = 2;

    double intensityInnerMin = std::numeric_limits<double>::max(); // set to maximum double value;

    haarAOI.xPos = 0;
    haarAOI.yPos = 0;

    int innerArea = (haarAOI.wdth - 1) * (haarAOI.hght - 1);

    int wdth = searchAOI.wdth - haarAOI.wdth;
    int hght = searchAOI.hght - haarAOI.hght;

    for (int iRow = 0; iRow < hght; iRow = iRow + stepSize)
    {
        for (int iCol = 0; iCol < wdth; iCol = iCol + stepSize)
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
                    if
                            (X[i] > xTopLeft &&
                             X[i] < xBtmRght &&
                             Y[j] > yTopLeft &&
                             Y[j] < yBtmRght) { GLINT_OVERLAP = true; }
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

            if (intensityInner < intensityInnerMin)
            {
                haarAOI.xPos = xTopLeft;
                haarAOI.yPos = yTopLeft;
                intensityInnerMin = intensityInner;
            }
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

    double pupilRadius = (mDetectionVariables.predictedCircumference * mDetectionVariables.thresholdChangeCircumference) / (2 * M_PI);

    // Find a starting edge point using Starburst-like algorithm

    std::vector<int> startIndices;

    for (int m = 0; m < 8; m++)
    {
        bool STOP_SEARCH = false;

        int x = round(mDetectionVariables.predictedXPosRelative + dX[m]);
        int y = round(mDetectionVariables.predictedYPosRelative + dY[m]);

        double R = 0;

        while (!STOP_SEARCH)
        {
            int dx = dX[m];
            int dy = dY[m];

            for (int n = 0; n < 2 && !STOP_SEARCH; n++)
            {
                x = x + dx * (1 - n);
                y = y + dy * (n);

                if (R > pupilRadius || x < 0 || x >= mAOI.wdth || y < 0 || y >= mAOI.hght)
                {
                    STOP_SEARCH = true;
                    break;
                }

                int centreIndex = y * mAOI.wdth + x;

                int tag = cannyEdgeVector[centreIndex];

                if (tag > 0)
                {
                    if (tag == 1) { startIndices.push_back(centreIndex); }
                }
            }

            R += sqrt((double) (dx * dx + dy * dy));
        }
    }

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

    double ratio;
    if (x != 0) { ratio = std::abs(y / x); }
    else        { ratio = 0; } // arbitrary choice

    int dir = 0;

    if (x >= 0) // Right half
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

    return dir;
}

int connectEdges(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1};
    std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1};

    std::vector<int> edgePoints = {startIndex};
    int edgePointNew = {startIndex};

    for (int iEdgePoint = 0; iEdgePoint < curvatureWindowLength; iEdgePoint++) // move back through edge
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

    if (edgeLength >= curvatureWindowLength)
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

int findEdgePoints(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
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
                int edgePointNew = connectEdges(cannyEdgeVector, mAOI, centreIndex); // connect possible edge terminals

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

void constructGraphTree(std::vector<vertexProperties>& vVertexProperties, std::vector<branchProperties>& vBranchProperties, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int pointIndexStart)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    vertexProperties vertexStart; // first vertex
    vertexStart.pointIndex = pointIndexStart;
    vertexStart.index      = 0;

    { // Find all branches connected to start vertex
        int centreXPos =  vertexStart.pointIndex % mAOI.wdth;
        int centreYPos = (vertexStart.pointIndex - centreXPos) / mAOI.wdth;

        for (int m = 0; m < 8; m++) // loop through 8-connected environment of the vertex
        {
            int neighbourXPos = centreXPos + dX[m];
            int neighbourYPos = centreYPos + dY[m];

            if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

            int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

            if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point is filled ...
            {
                cannyEdgeVector[neighbourIndex] = 3; // ... then tag it
                vertexStart.connectedPoints.push_back(neighbourIndex);
            }
        }
    }

    std::vector<vertexProperties> verticesOld = {vertexStart}; // store first vertex

    int branchNumber = 0; // counters
    int vertexNumber = 0;

    int numVertices;

    do
    {
        std::vector<vertexProperties> verticesNew; // store new found vertices

        numVertices = verticesOld.size();

        for (int iVertex = 0; iVertex < numVertices; iVertex++)
        {
            vertexProperties vertexCurrent = verticesOld[iVertex]; // run through all found vertices
            cannyEdgeVector[vertexCurrent.pointIndex] = 3; // tag vertex

            // Loop through all found branches (if any) to find all edge points belonging to each branch

            int numEdges = vertexCurrent.connectedPoints.size();

            for (int iBranch = 0; iBranch < numEdges || branchNumber == 0; iBranch++)
            {
                branchProperties branchNew;
                if (branchNumber == 0) { branchNew.pointIndices.push_back(vertexCurrent.pointIndex); } // first vertex should be included in first branch

                if (numEdges > 0)
                {
                    int edgePointIndexOld = vertexCurrent.connectedPoints[iBranch];
                    branchNew.pointIndices.push_back(edgePointIndexOld);

                    bool VERTEX_FOUND = false;

                    while(!VERTEX_FOUND) // keep going until vertex is encountered
                    {
                        std::vector<int> edgePointsNew;

                        int centreIndex = edgePointIndexOld;
                        int centreXPos  = centreIndex % mAOI.wdth;
                        int centreYPos  = (centreIndex - centreXPos) / mAOI.wdth;

                        for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                        {
                            int neighbourXPos = centreXPos + dX[m];
                            int neighbourYPos = centreYPos + dY[m];

                            if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                            int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                            if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point was tagged previously ...
                            {
                                cannyEdgeVector[neighbourIndex] = 3; // ... give it a new tag
                                edgePointsNew.push_back(neighbourIndex);
                            }
                        }

                        int numEdgePoints = edgePointsNew.size();

                        if (numEdgePoints == 1) // edge continues
                        {
                            edgePointIndexOld = edgePointsNew[0]; // update index
                            branchNew.pointIndices.push_back(edgePointIndexOld); // check new index
                            edgePointsNew.clear();
                        }
                        else // new vertex found
                        {
                            vertexProperties vertexNew;
                            vertexNew.pointIndex = centreIndex;
                            vertexNew.connectedBranches.push_back(branchNumber); // add current branch to new vertex connections
                            if (edgePointsNew.size() > 0) { vertexNew.connectedPoints = edgePointsNew; }
                            verticesNew.push_back(vertexNew); // new vertices to be checked
                            VERTEX_FOUND = true;
                        }
                    }
                }

                vertexCurrent.connectedBranches.push_back(branchNumber); // add new found branches

                int branchLength = branchNew.pointIndices.size();

                if (branchLength > 0)
                {
                    branchNew.length = branchLength;
                    branchNew.index  = branchNumber;
                    vBranchProperties.push_back(branchNew); // record all branches
                }

                branchNumber++;
            }

            vertexCurrent.index = vertexNumber;
            vVertexProperties.push_back(vertexCurrent); // record all vertices
            vertexNumber++;
        }

        numVertices = verticesNew.size();
        verticesOld = verticesNew;

    } while (numVertices > 0);

    // Find vertices that branches are connected to

    numVertices = vVertexProperties.size();

    for (int iVertex = 0; iVertex < numVertices; iVertex++)
    {
        std::vector<int> vertexConnections = vVertexProperties[iVertex].connectedBranches;

        for (int iEdge = 0, numConnections = vertexConnections.size(); iEdge < numConnections; iEdge++)
        {
            vBranchProperties[vertexConnections[iEdge]].connectedVertices.push_back(iVertex);
        }
    }
}

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

void findLongestPath(std::vector<int>& allPoints, std::vector<int>& pathPoints, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int edgePointStart)
{
    std::vector<branchProperties> vBranchPropertiesAll;
    std::vector<vertexProperties> vVertexPropertiesAll;

    // Find all vertices and all connected branches (i.e. obtain graph tree)

    constructGraphTree(vVertexPropertiesAll, vBranchPropertiesAll, cannyEdgeVector, mAOI, edgePointStart);

    int numBranchesAll = vBranchPropertiesAll.size();
    int numVerticesAll = vVertexPropertiesAll.size();

    // Find longest path in edge collection

    std::vector<std::vector<int>> pathAll(numVerticesAll);

    for (int iVertex = 0; iVertex < numVerticesAll; iVertex++) // loop through all vertices
    {
        vertexProperties vertexStart = vVertexPropertiesAll[iVertex]; // starting vertex

        if (vertexStart.connectedBranches.size() > 1) { continue; } // only start with terminal vertices

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
                    std::vector<int> connectedVertices = vBranchPropertiesAll[branchNumber].connectedVertices; // get the vertices the branch is connected to (almost always 2, but can be 1)

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

std::vector<edgeProperties> edgeSelection(const detectionVariables& mDetectionVariables, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
{
    std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all edges

    std::vector<int> startIndicesRaw = findEdges(mDetectionVariables, cannyEdgeVector, mAOI);

    int numEdges = startIndicesRaw.size();
    std::vector<int> startIndices(numEdges);

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        startIndices[iEdge] = findEdgePoints(cannyEdgeVector, mAOI, startIndicesRaw[iEdge]); // tag all edge points of found edges
    }

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
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
            {  cannyEdgeVector[edgePointIndex] = 1; }
        }

        edgeProperties mEdgeProperties;
        mEdgeProperties.pointIndices = pathIndices;

        vEdgePropertiesAll.push_back(mEdgeProperties);
    }

    return vEdgePropertiesAll;
}

double calculateCurvatureAverage(double circumference, double aspectRatio)
{
    int arrayWdth = arrayCircumferences.size();

    int arrayXPos = *std::upper_bound(arrayCircumferences.begin(), arrayCircumferences.end(), circumference);
    int arrayYPos = *std::upper_bound(arrayAspectRatios.begin(),   arrayAspectRatios.end(),   aspectRatio);

    double curvatureAverage = 0.5 * (arrayCurvatureMax[arrayYPos * arrayWdth + arrayXPos] + arrayCurvatureMin[arrayYPos * arrayWdth + arrayXPos]);

    return curvatureAverage;
}

void calculateCurvatureLimits(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, double& curvatureUpperLimit, double& curvatureLowerLimit)
{
    // Calculate curvature limits

    double circumferenceLow  = mDetectionVariables.predictedCircumference * (1 / mDetectionVariables.thresholdChangeCircumference);
    double circumferenceHigh = mDetectionVariables.predictedCircumference *      mDetectionVariables.thresholdChangeCircumference;

    double aspectRatioLow  = mDetectionVariables.predictedAspectRatio * (1 / mDetectionVariables.thresholdChangeAspectRatio);
    double aspectRatioHigh = mDetectionVariables.predictedAspectRatio *      mDetectionVariables.thresholdChangeAspectRatio;

    std::vector<int> arrayXPos(2);
    arrayXPos[0] = findUpperBound(arrayCircumferences, circumferenceLow);
    arrayXPos[1] = findUpperBound(arrayCircumferences, circumferenceHigh);

    std::vector<int> arrayYPos(2);
    arrayYPos[0] = findUpperBound(arrayAspectRatios, aspectRatioLow);
    arrayYPos[1] = findUpperBound(arrayAspectRatios, aspectRatioHigh);

    // Calculate limits

    std::vector<double> curvaturesMax(4);
    std::vector<double> curvaturesMin(4);

    int arrayWidth = arrayCircumferences.size();

    for (int i = 0; i < 2; i++)
    {
        int x = arrayXPos[i];

        for (int j = 0; j < 2; j++)
        {
            int y = arrayYPos[j];
            curvaturesMax[2 * j + i] = arrayCurvatureMax[y * arrayWidth + x];
            curvaturesMin[2 * j + i] = arrayCurvatureMin[y * arrayWidth + x];
        }
    }

    curvatureUpperLimit = *std::max_element(std::begin(curvaturesMax), std::end(curvaturesMax)) + mDetectionParameters.curvatureOffset;
    curvatureLowerLimit = *std::min_element(std::begin(curvaturesMin), std::end(curvaturesMin)) - mDetectionParameters.curvatureOffset;
}

std::vector<double> calculateCurvatures(std::vector<double>& xNormals, std::vector<double>& yNormals, const std::vector<double>& xTangentsAll, const std::vector<double>& yTangentsAll)
{
    int edgeSize = xTangentsAll.size();

    std::vector<double> curvatures(edgeSize, 0.0);

    int numPos = 0;
    int numNeg = 0;

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
    {
        // calculate window tangents

        // first window

        std::vector<double> tangentsX_1(xTangentsAll.begin() + iEdgePoint - curvatureWindowLength, xTangentsAll.begin() + iEdgePoint - 1);
        std::vector<double> tangentsY_1(yTangentsAll.begin() + iEdgePoint - curvatureWindowLength, yTangentsAll.begin() + iEdgePoint - 1);

        double tangentXMean_1 = calculateMean(tangentsX_1);
        double tangentYMean_1 = calculateMean(tangentsY_1);

        // second window

        std::vector<double> tangentsX_2(xTangentsAll.begin() + iEdgePoint + 1, xTangentsAll.begin() + iEdgePoint + curvatureWindowLength);
        std::vector<double> tangentsY_2(yTangentsAll.begin() + iEdgePoint + 1, yTangentsAll.begin() + iEdgePoint + curvatureWindowLength);

        double tangentXMean_2 = calculateMean(tangentsX_2);
        double tangentYMean_2 = calculateMean(tangentsY_2);

        // calculate vector difference

        double vectorAngle = atan2(tangentYMean_2, tangentXMean_2) - atan2(tangentYMean_1, tangentXMean_1);

        if      (vectorAngle >  M_PI) { vectorAngle = vectorAngle - 2 * M_PI; }
        else if (vectorAngle < -M_PI) { vectorAngle = vectorAngle + 2 * M_PI; }

        if      (vectorAngle > 0) { numPos++; }
        else if (vectorAngle < 0) { numNeg++; }

        curvatures[iEdgePoint] = 180 * vectorAngle / M_PI; // in degrees

        xNormals[iEdgePoint] = tangentXMean_2 - tangentXMean_1;
        yNormals[iEdgePoint] = tangentYMean_2 - tangentYMean_1;
    }

    if (numNeg > numPos) // if majority sign is negative, then swap all signs
    {
        for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
        {
            curvatures[iEdgePoint] = -curvatures[iEdgePoint];
        }
    }

    return curvatures;
}

std::vector<edgeProperties> edgeSegmentationCurvature(const edgeProperties& mEdgeProperties, const double curvatureUpperLimit, const double curvatureLowerLimit)
{
    int edgeSize = mEdgeProperties.curvatures.size();

    // find breakpoints based on curvature thresholding

    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(-1); // add first point (+ 1 is added later)

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
    {
        double curvature = mEdgeProperties.curvatures[iEdgePoint];

        if (curvature >= curvatureUpperLimit || curvature <= curvatureLowerLimit)
        {
            breakPoints.push_back(iEdgePoint);
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

        if (diffLeft > minimumEdgeLength || diffRght > minimumEdgeLength)
        {
            breakPointsNew.push_back(iEdgePointCntr);
        }
    }

    breakPointsNew.push_back(edgeSize - 1); // add last point again

    // cut edge at breakpoints

    std::vector<edgeProperties> vEdgePropertiesNew;

    int numBreakPointsNew = breakPointsNew.size();

    for (int iBreakPoint = 0; iBreakPoint < numBreakPointsNew - 1; iBreakPoint++)
    {
        int iStartBreakPoint = breakPointsNew[iBreakPoint] + 1;
        int iEndBreakPoint   = breakPointsNew[iBreakPoint  + 1];

        edgeProperties mEdgePropertiesNew;
        mEdgePropertiesNew.pointIndices = std::vector<int>    (mEdgeProperties.pointIndices.begin() + iStartBreakPoint, mEdgeProperties.pointIndices.begin() + iEndBreakPoint);
        mEdgePropertiesNew.curvatures   = std::vector<double> (mEdgeProperties.curvatures.begin()   + iStartBreakPoint, mEdgeProperties.curvatures.begin()   + iEndBreakPoint);
        mEdgePropertiesNew.xnormals     = std::vector<double> (mEdgeProperties.xnormals.begin()     + iStartBreakPoint, mEdgeProperties.xnormals.begin()     + iEndBreakPoint);
        mEdgePropertiesNew.ynormals     = std::vector<double> (mEdgeProperties.ynormals.begin()     + iStartBreakPoint, mEdgeProperties.ynormals.begin()     + iEndBreakPoint);

        vEdgePropertiesNew.push_back(mEdgePropertiesNew);
    }

    return vEdgePropertiesNew;
}

double calculateEdgeLength(const std::vector<int>& edgePoints, const AOIProperties& mAOI)
{
    // calculates edge length without need of edge continuity

    double lengthTotal = 0;
    int edgeSize = edgePoints.size();

    int iEdgePoint  = 0;
    bool BREAK_LOOP = false;

    do // run through all edge points
    {
        int jEdgePoint = iEdgePoint + lengthWindowLength;

        if (jEdgePoint >= edgeSize)
        {
            BREAK_LOOP = true;
            jEdgePoint = edgeSize - 1;
        }

        int edgePointIndex_1 = edgePoints[iEdgePoint];
        int edgePointIndex_2 = edgePoints[jEdgePoint];

        int edgePointXPos_1  = edgePointIndex_1 % mAOI.wdth;
        int edgePointXPos_2  = edgePointIndex_2 % mAOI.wdth;

        int edgePointYPos_1  = (edgePointIndex_1 - edgePointXPos_1) / mAOI.wdth;
        int edgePointYPos_2  = (edgePointIndex_2 - edgePointXPos_2) / mAOI.wdth;

        int dX = edgePointXPos_1 - edgePointXPos_2;
        int dY = edgePointYPos_1 - edgePointYPos_2;

        lengthTotal += sqrt(dX * dX + dY * dY);

        iEdgePoint = iEdgePoint + lengthWindowLength;
    } while (!BREAK_LOOP);

    return lengthTotal;
}

std::vector<edgeProperties> edgeSegmentationLength(const detectionVariables& mDetectionVariables, const edgeProperties& mEdgeProperties, const AOIProperties& mAOI)
{
    // This functions cuts edge terminals to make the edge shorter, if the edge is significantly longer than predicted

    int edgeSize = mEdgeProperties.curvatures.size();

    // find breakpoints based on length thresholding

    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(-1); // add first point (+ 1 is added later)

    double lengthDifference = mEdgeProperties.length - mDetectionVariables.predictedCircumference;

    if (lengthDifference > lengthWindowLength) // difference should be large enough
    {
        if (edgeSize > 2 * lengthDifference + lengthWindowLength) // should be enough space between breakpoints
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
            mEdgePropertiesNew.curvature = calculateMean(mEdgePropertiesNew.curvatures);
            mEdgePropertiesNew.length    = calculateEdgeLength(mEdgePropertiesNew.pointIndices, mAOI);

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
        inputVector_1[6] = vEdgeProperties[2].length;

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

std::vector<edgeProperties> edgeTerminalFilter(const detectionVariables& mDetectionVariables, const edgeProperties& mEdgeProperties, const AOIProperties& mAOI, double scoreThresholdDiff)
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

            if (edgeSize <= 3 * curvatureWindowLength) { continue; } // don't segment short edges

            int iEdgePoint  = 0;
            bool BREAK_LOOP = false;

            //            int numBreakPoints = ceil(edgeSize / curvatureWindowLength) + 1;
            std::vector<int> breakPointsAll;
            std::vector<double> scoreDifference;
            int iBreakPoint = 0;

            do // run through all edge points
            {
                iEdgePoint = iEdgePoint + curvatureWindowLength;

                if (iEdgePoint >= edgeSize - curvatureWindowLength - 1)
                {
                    BREAK_LOOP = true;
                    iEdgePoint = edgeSize - curvatureWindowLength - 1;
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

                    mEdgePropertiesNew.radius    = calculateMean(mEdgePropertiesNew.radii);
                    mEdgePropertiesNew.curvature = calculateMean(mEdgePropertiesNew.curvatures);
                    mEdgePropertiesNew.intensity = calculateMeanInt(mEdgePropertiesNew.intensities);
                    mEdgePropertiesNew.gradient  = calculateMeanInt(mEdgePropertiesNew.gradients);
                    mEdgePropertiesNew.length    = calculateEdgeLength(mEdgePropertiesNew.pointIndices, mAOI);

                    vEdgePropertiesTemp[i] = mEdgePropertiesNew;
                }

                // Shortest edge length

                double edgeLength;
                double edgeLength_1 = vEdgePropertiesTemp[0].length;
                double edgeLength_2 = vEdgePropertiesTemp[1].length;
                if (edgeLength_1 < edgeLength_2) { edgeLength = edgeLength_1; }
                else                             { edgeLength = edgeLength_2; }

                // calculate score difference

                std::vector<double> inputVector(7);
                inputVector[0] = vEdgePropertiesTemp[0].radius    / vEdgePropertiesTemp[1].radius;
                inputVector[1] = 0; // don't use circumference difference
                inputVector[2] = vEdgePropertiesTemp[0].curvature - vEdgePropertiesTemp[1].curvature;
                inputVector[3] = vEdgePropertiesTemp[0].intensity - vEdgePropertiesTemp[1].intensity;
                inputVector[4] = vEdgePropertiesTemp[0].gradient  - vEdgePropertiesTemp[1].gradient;
                inputVector[5] = 0; // don't use radial variance
                inputVector[6] = edgeLength;

                scoreDifference.push_back(calculateScoreTotal(mDetectionVariables,inputVector, false));
                breakPointsAll.push_back(iEdgePoint);
                iBreakPoint++;
            } while (!BREAK_LOOP);

            // Find point of maximum score difference

            auto itr = std::min_element(scoreDifference.begin(), scoreDifference.end());
            double scoreDifferenceMax = *itr;
            int indexMax = std::distance(scoreDifference.begin(), itr);

            std::vector<edgeProperties> vEdgePropertiesTemp;

            if (scoreDifferenceMax < scoreThresholdDiff)
            {
                std::vector<int> breakPoints = {0, indexMax, indexMax + 1, edgeSize - 1};

                for (int i = 0; i < 2; i++) // split edge up
                {
                    edgeProperties mEdgePropertiesNew;

                    mEdgePropertiesNew.pointIndices = std::vector<int>   (mEdgePropertiesOld.pointIndices.begin() + breakPoints[2 * i], mEdgePropertiesOld.pointIndices.begin() + breakPoints[2 * i + 1]);
                    mEdgePropertiesNew.intensities  = std::vector<int>   (mEdgePropertiesOld.intensities.begin()  + breakPoints[2 * i], mEdgePropertiesOld.intensities.begin()  + breakPoints[2 * i + 1]);
                    mEdgePropertiesNew.gradients    = std::vector<int>   (mEdgePropertiesOld.gradients.begin()    + breakPoints[2 * i], mEdgePropertiesOld.gradients.begin()    + breakPoints[2 * i + 1]);
                    mEdgePropertiesNew.radii        = std::vector<double>(mEdgePropertiesOld.radii.begin()        + breakPoints[2 * i], mEdgePropertiesOld.radii.begin()        + breakPoints[2 * i + 1]);
                    mEdgePropertiesNew.curvatures   = std::vector<double>(mEdgePropertiesOld.curvatures.begin()   + breakPoints[2 * i], mEdgePropertiesOld.curvatures.begin()   + breakPoints[2 * i + 1]);

                    vEdgePropertiesTemp.push_back(mEdgePropertiesOld);
                }
            }
            else { vEdgePropertiesTemp.push_back(mEdgePropertiesOld); }

            vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
        }

        vEdgePropertiesOld = vEdgePropertiesNew;
        vEdgePropertiesAll.insert(vEdgePropertiesAll.end(), vEdgePropertiesOld.begin(), vEdgePropertiesOld.end());

    } while (vEdgePropertiesOld.size() > 1);

    return vEdgePropertiesAll;
}

void calculateCurvatureStats(edgeProperties& mEdgeProperties)
{
    // Calculate min, max and mean curvature

    double curvature    = 0;
    double curvatureMax = 0;
    double curvatureMin = 0;

    std::vector<double> edgeCurvaturesNew;

    int edgeSize = mEdgeProperties.curvatures.size();

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
    {
        double curvature = mEdgeProperties.curvatures[iEdgePoint];
        if (curvature < 180.0) { edgeCurvaturesNew.push_back(curvature); }
    }

    int edgeLengthNew = edgeCurvaturesNew.size();

    if (edgeLengthNew > 0)
    {
        curvature    = calculateMean(edgeCurvaturesNew);
        curvatureMax = *std::max_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
        curvatureMin = *std::min_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
    }
    else { curvature = 360; curvatureMax = 360; curvatureMin = 360; }

    mEdgeProperties.curvature    = curvature;
    mEdgeProperties.curvatureMax = curvatureMax;
    mEdgeProperties.curvatureMin = curvatureMin;
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

std::vector<edgeProperties> removeShortEdges(const std::vector<edgeProperties>& vEdgeProperties)
{
    std::vector<edgeProperties> vEdgePropertiesNew;
    for (int iEdge = 0, numEdges = vEdgeProperties.size(); iEdge < numEdges; iEdge++) // ignore short edges
    {
        edgeProperties mEdgeProperties = vEdgeProperties[iEdge];
        int edgeSize = mEdgeProperties.pointIndices.size();
        if (edgeSize > minimumEdgeLength) { vEdgePropertiesNew.push_back(mEdgeProperties); }
    }

    return vEdgePropertiesNew;
}

std::vector<int> edgeClassification(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgePropertiesAll)
{
    int numEdgesMax = mDetectionParameters.ellipseFitNumberMaximum;
    int numEdges    = vEdgePropertiesAll.size();
    if (numEdgesMax > numEdges) { numEdgesMax = numEdges; }

    // Classify edges based on score

    std::vector<double> totalScores(numEdges);
    std::vector<int> pupilEdges;

    if (mDetectionVariables.certaintyAverages > 0 || mDetectionVariables.certaintyFeatures > 0 || mDetectionVariables.certaintyPosition > 0) // need to have enough certainty
    {
        for (int iEdge = 0; iEdge < numEdges; iEdge++)
        {
            std::vector<double> inputVector(7);

            inputVector[0] = vEdgePropertiesAll[iEdge].radius    / (mDetectionVariables.predictedCircumference / (2 * M_PI));
            inputVector[1] = vEdgePropertiesAll[iEdge].length    / mDetectionVariables.predictedCircumference;
            inputVector[2] = vEdgePropertiesAll[iEdge].curvature - mDetectionVariables.predictedCurvature;
            inputVector[3] = vEdgePropertiesAll[iEdge].intensity - mDetectionVariables.averageIntensity;
            inputVector[4] = vEdgePropertiesAll[iEdge].gradient  - mDetectionVariables.averageGradient;
            inputVector[5] = vEdgePropertiesAll[iEdge].radiusVar;
            inputVector[6] = vEdgePropertiesAll[iEdge].length;

            totalScores[iEdge] = calculateScoreTotal(mDetectionVariables,inputVector, true);
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

    std::vector<double> v(6);

    v[0] = semiMajor;
    v[1] = semiMinor;
    v[2] = x;
    v[3] = y;
    v[4] = w;
    v[5] = h;
    v[6] = alpha;

    return v;
}

ellipseProperties fitEllipse(std::vector<int> edgeIndices, int edgeSetSize, int haarWidth)
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

        double edgePointX = edgePointIndex % haarWidth;
        double edgePointY = (edgePointIndex - edgePointX) / haarWidth;

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
    mEllipseProperties.radius        = 0.5 * (semiMinor + semiMajor);
    mEllipseProperties.xPos          = ellipseParameters[2];
    mEllipseProperties.yPos          = ellipseParameters[3];
    mEllipseProperties.width         = ellipseParameters[4];
    mEllipseProperties.height        = ellipseParameters[5];
    mEllipseProperties.coefficients  = ellipseFitCoefficients;

    for (int iParameter = 0; iParameter < 6; iParameter++)
    {
        if (std::isnan(ellipseParameters[iParameter])) { mEllipseProperties.DETECTED = false; } // ERROR
    }

    return mEllipseProperties;
}

std::vector<ellipseProperties> getEllipseFits(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI)
{
    ellipseProperties mEllipseProperties;
    mEllipseProperties.DETECTED = false;

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

            double circumferenceUpperLimit = mDetectionVariables.averageCircumference * mDetectionVariables.offsetCircumference;
            double circumferenceLowerLimit = mDetectionVariables.averageCircumference / mDetectionVariables.offsetCircumference;

            if (circumferenceUpperLimit > mDetectionParameters.circumferenceMax) { circumferenceUpperLimit = mDetectionParameters.circumferenceMax; }
            if (circumferenceLowerLimit < mDetectionParameters.circumferenceMin) { circumferenceLowerLimit = mDetectionParameters.circumferenceMin; }

            if (circumferenceApprox > circumferenceUpperLimit) { continue; } // no large ellipse
            if (circumferenceApprox < circumferenceLowerLimit) { continue; } // no small ellipse

            // Fit ellipse

            ellipseProperties mEllipsePropertiesNew = fitEllipse(edgeIndices, edgeSetSize, mAOI.wdth);

            if (!mEllipsePropertiesNew.DETECTED) { continue; } // error

            // Absolute displacement filter

            double dX = mEllipsePropertiesNew.xPos - (mDetectionVariables.predictedXPos - mAOI.xPos);
            double dY = mEllipsePropertiesNew.yPos - (mDetectionVariables.predictedYPos - mAOI.yPos);
            double dR = sqrt(dX * dX + dY * dY);
            if (dR > mDetectionVariables.thresholdChangePosition) { continue; } // no large ellipse displacements

            // Absolute size and shape filter

            if (mEllipsePropertiesNew.circumference > circumferenceUpperLimit) { continue; } // no large ellipse
            if (mEllipsePropertiesNew.circumference < circumferenceLowerLimit) { continue; } // no small ellipse
            if (mEllipsePropertiesNew.aspectRatio   < mDetectionParameters.aspectRatioMin)   { continue; } // no extreme deviations from circular shape

            // Relative change in size and shape filter

            double relativeChangeAspectRatio   = mEllipsePropertiesNew.aspectRatio   / mDetectionVariables.predictedAspectRatio;
            double relativeChangeCircumference = mEllipsePropertiesNew.circumference / mDetectionVariables.predictedCircumference;

            if (relativeChangeAspectRatio   < 1.0) { relativeChangeAspectRatio   = 1 / relativeChangeAspectRatio; }
            if (relativeChangeCircumference < 1.0) { relativeChangeCircumference = 1 / relativeChangeCircumference; }

            if (relativeChangeAspectRatio   > mDetectionVariables.thresholdChangeAspectRatio  ) { continue; } // no large ellipse shape changes
            if (relativeChangeCircumference > mDetectionVariables.thresholdChangeCircumference) { continue; } // no large ellipse size changes

            // Calculate error between fit and every edge point

            double A = mEllipsePropertiesNew.coefficients[0];
            double B = mEllipsePropertiesNew.coefficients[1];
            double C = mEllipsePropertiesNew.coefficients[2];
            double D = mEllipsePropertiesNew.coefficients[3];
            double E = mEllipsePropertiesNew.coefficients[4];
            double F = mEllipsePropertiesNew.coefficients[5];

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
            std::vector<double> fitErrorsMax(fitErrorsSorted.begin(), fitErrorsSorted.begin() + round(fitErrorFraction * edgeSetLength));
            double fitErrorMax = calculateMean(fitErrorsMax);

            // Fit error filter

            if (fitErrorMax > mDetectionParameters.ellipseFitErrorMaximum) { continue; } // no large fit errors

            mEllipsePropertiesNew.fitError = fitErrorMax;

            // Save parameters of accepted fit

            mEllipsePropertiesNew.edgeIndices = combiEdgeIndices;
            mEllipsePropertiesNew.edgeLength  = edgeSetLength;
            vEllipsePropertiesAll.push_back(mEllipsePropertiesNew);
        }
        while (std::next_permutation(edgeCombination.begin(), edgeCombination.end()));
    }

    return vEllipsePropertiesAll;
}

int ellipseFitFilter(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, std::vector<ellipseProperties> vEllipseProperties)
{
    int numFits = vEllipseProperties.size();

    std::vector<double> featureChange(numFits); // new type of fit error

    double maxScoreAspectRatio   = 1;
    double maxScoreCircumference = 1;
    double maxScoreDisplacement  = 1;
    double maxScoreFitError      = 2;
    double maxScoreLength        = 2;

    double scoreAspectRatio   = 0;
    double scoreCircumference = 0;
    double scoreDisplacement  = 0;

    double certaintyFactorPosition = 0.5 * (mDetectionVariables.certaintyPosition + 1);
    double certaintyFactorFeatures = 0.5 * (mDetectionVariables.certaintyFeatures + 1);

    for (int iFit = 0; iFit < numFits; iFit++)
    {
        ellipseProperties mEllipseProperties = vEllipseProperties[iFit];

        double dx = mEllipseProperties.xPos - mDetectionVariables.predictedXPos;
        double dy = mEllipseProperties.yPos - mDetectionVariables.predictedYPos;
        double displacementChange  = sqrt(dx * dx + dy * dy);
        double circumferenceChange = (std::abs(mEllipseProperties.circumference - mDetectionVariables.predictedCircumference));
        double aspectRatioChange   = (std::abs(mEllipseProperties.aspectRatio   - mDetectionVariables.predictedAspectRatio));
        double lengthDifference    = (std::abs(mEllipseProperties.edgeLength    - mDetectionVariables.predictedCircumference));
        double fitError            = mEllipseProperties.fitError;

        scoreAspectRatio   = (-maxScoreAspectRatio   / mDetectionParameters.thresholdChangeAspectRatio)   * aspectRatioChange   + maxScoreAspectRatio;
        scoreCircumference = (-maxScoreCircumference / mDetectionParameters.thresholdChangeCircumference) * circumferenceChange + maxScoreCircumference;
        scoreDisplacement  = (-maxScoreDisplacement  / mDetectionParameters.thresholdChangePosition)      * displacementChange  + maxScoreDisplacement;

        double scoreFitError = (-maxScoreFitError   / mDetectionParameters.ellipseFitErrorMaximum)  * fitError         + maxScoreFitError;
        double scoreLength   = (-maxScoreLength * 2 / mDetectionVariables.predictedCircumference) * lengthDifference + maxScoreLength;

        if (scoreCircumference  < 0) { scoreCircumference   = 0; }
        if (scoreAspectRatio    < 0) { scoreAspectRatio     = 0; }
        if (scoreFitError       < 0) { scoreFitError        = 0; }
        if (scoreLength         < 0) { scoreLength          = 0; }

        featureChange[iFit] = certaintyFactorFeatures * (scoreAspectRatio + scoreCircumference) + certaintyFactorPosition * scoreDisplacement + scoreFitError + scoreLength;
    }

    int acceptedFitIndex = std::distance(featureChange.begin(), std::max_element(featureChange.begin(), featureChange.end()));

    return acceptedFitIndex;
}

void checkVariableLimits(detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters)
{
    if (mDetectionVariables.thresholdChangeAspectRatio      < mDetectionParameters.thresholdChangeAspectRatio)
    {        mDetectionVariables.thresholdChangeAspectRatio = mDetectionParameters.thresholdChangeAspectRatio; }

    if (mDetectionVariables.thresholdChangeCircumference      < mDetectionParameters.thresholdChangeCircumference)
    {        mDetectionVariables.thresholdChangeCircumference = mDetectionParameters.thresholdChangeCircumference; }

    if (mDetectionVariables.thresholdChangePosition      < mDetectionParameters.thresholdChangePosition)
    {        mDetectionVariables.thresholdChangePosition = mDetectionParameters.thresholdChangePosition; }

    if      (mDetectionVariables.certaintyPosition < -1.0) { mDetectionVariables.certaintyPosition = -1.0; }
    else if (mDetectionVariables.certaintyPosition >  1.0) { mDetectionVariables.certaintyPosition =  1.0; }

    if      (mDetectionVariables.certaintyFeatures < -1.0) { mDetectionVariables.certaintyFeatures = -1.0; }
    else if (mDetectionVariables.certaintyFeatures >  1.0) { mDetectionVariables.certaintyFeatures =  1.0; }

    if      (mDetectionVariables.certaintyAverages  < -1.0) { mDetectionVariables.certaintyAverages  = -1.0; }
    else if (mDetectionVariables.certaintyAverages  >  1.0) { mDetectionVariables.certaintyAverages  =  1.0; }
}

detectionVariables eyeStalker(const cv::Mat& imageOriginalBGR, const AOIProperties& mAOI, detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, dataVariables& mDataVariables, drawVariables& mDrawVariables)
{
    // Keep variables within limits

    checkVariableLimits(mDetectionVariables, mDetectionParameters);

    // Define search area

    int imageWdth = imageOriginalBGR.cols;

    AOIProperties searchAOI;

    searchAOI.xPos = round(mDetectionVariables.predictedXPos - mDetectionVariables.thresholdChangePosition);
    searchAOI.yPos = round(mDetectionVariables.predictedYPos - mDetectionVariables.thresholdChangePosition);
    int searchEndX = round(mDetectionVariables.predictedXPos + mDetectionVariables.thresholdChangePosition);
    int searchEndY = round(mDetectionVariables.predictedYPos + mDetectionVariables.thresholdChangePosition);
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

    if (searchAOI.wdth > haarAOI.wdth || searchAOI.hght > haarAOI.hght)
    {
        std::vector<unsigned int> integralImage = calculateIntImg(imageOriginalGray, imageWdth, searchAOI);

        glintAOI.wdth = mDetectionParameters.glintWdth;
        glintAOI.hght = glintAOI.wdth;
        glintAOI      = detectGlint(imageOriginalGray, imageWdth, searchAOI, glintAOI);
        glintAOI.xPos = searchAOI.xPos + glintAOI.xPos;
        glintAOI.yPos = searchAOI.yPos + glintAOI.yPos;

        haarAOI      = detectPupilApprox(integralImage, searchAOI, haarAOI, glintAOI);
        haarAOI.xPos = searchAOI.xPos + haarAOI.xPos;
        haarAOI.yPos = searchAOI.yPos + haarAOI.yPos;
    }
    else // search not required
    {
        haarAOI.xPos = searchAOI.xPos;
        haarAOI.yPos = searchAOI.yPos;
        haarAOI.flag = false; glintAOI.flag = false;
    }

    // Create new AOI for Canny edge deteciton

    double AOIXPos   = haarAOI.xPos + 0.5 * haarAOI.wdth; // centre of haar-like detector
    double AOIYPos   = haarAOI.yPos + 0.5 * haarAOI.hght;
    double AOIWdth   = mDetectionVariables.predictedWidth;
    double AOIHght   = mDetectionVariables.predictedHeight;
    double AOIOffset = mDetectionVariables.thresholdChangePosition + mDetectionVariables.thresholdChangeCircumference / (2 * M_PI);

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
    std::vector<int> cannyEdges = cannyConversion(imageCannyEdges, cannyAOI); // convert to binary vector
    std::vector<int> cannyEdgesSharpened = sharpenEdges(cannyEdges, cannyAOI); // morphological operation
    std::vector<int> edgeIndices = getEdgeIndices(cannyEdgesSharpened, 1); // used for drawing function

    /////////////////////////////////////////////////////////////////////////////
    //////////////////////////// EDGE SELECTION   ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    { // calculate predicted positions relative to Canny edge AOI
        double certaintyFactorPosition = 0.5 * (mDetectionVariables.certaintyPosition + 1);
        mDetectionVariables.predictedXPosRelative = (1 - certaintyFactorPosition) * AOIXPos + certaintyFactorPosition * mDetectionVariables.predictedXPos - cannyAOI.xPos;
        mDetectionVariables.predictedYPosRelative = (1 - certaintyFactorPosition) * AOIYPos + certaintyFactorPosition * mDetectionVariables.predictedYPos - cannyAOI.yPos;
    }

    std::vector<edgeProperties> vEdgePropertiesAll = edgeSelection(mDetectionVariables, cannyEdgesSharpened, cannyAOI);
    vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

    /////////////////////////////////////////////////////////////////////////////
    //////////////////////////// EDGE SEGMENTATION   ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    // Calculate curvature limits

    double curvatureUpperLimit;
    double curvatureLowerLimit;

    calculateCurvatureLimits(mDetectionVariables, mDetectionParameters, curvatureUpperLimit, curvatureLowerLimit);

    // Curvature calculation and edge segmentation

    { std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            int edgeSize = mEdgeProperties.pointIndices.size();

            std::vector<double> edgeXTangents(edgeSize, 0.0);
            std::vector<double> edgeYTangents(edgeSize, 0.0);

            calculateEdgeDirections(mEdgeProperties.pointIndices, edgeXTangents, edgeYTangents, cannyAOI);

            std::vector<double> edgeXNormals(edgeSize, 0.0);
            std::vector<double> edgeYNormals(edgeSize, 0.0);

            mEdgeProperties.curvatures = calculateCurvatures(edgeXNormals, edgeYNormals, edgeXTangents, edgeYTangents);
            mEdgeProperties.xnormals   = edgeXNormals;
            mEdgeProperties.ynormals   = edgeYNormals;

            // Do segmentation
            std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationCurvature(mEdgeProperties, curvatureUpperLimit, curvatureLowerLimit);
            vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
    }

    vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

    // Calculate additional edge properties and do length segmentation

    { std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            mEdgeProperties.length      = calculateEdgeLength(mEdgeProperties.pointIndices, cannyAOI);
            mEdgeProperties.radii       = calculateEdgeRadii(mEdgeProperties, cannyAOI, mDetectionVariables.predictedXPosRelative, mDetectionVariables.predictedYPosRelative);
            mEdgeProperties.gradients   = calculateRadialGradients(mDetectionVariables, mDetectionParameters, imageAOIGray, mEdgeProperties.pointIndices);
            mEdgeProperties.intensities = findEdgeIntensities(imageAOIGray, mDetectionParameters, mEdgeProperties, cannyAOI);

            // Do segmentation
            std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationLength(mDetectionVariables, mEdgeProperties, cannyAOI);
            vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
    }

    vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// EDGE TERMINAL FILTER   ////////////////////////
    ////////////////////////////////////////////////////////////////////////

    //    { std::vector<edgeProperties> vEdgePropertiesNew;

    //        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
    //        {
    //            edgeProperties mEdgeProperties    = vEdgePropertiesAll[iEdge];
    //            std::vector<edgeProperties> vEdgePropertiesTemp = edgeTerminalFilter(mDetectionVariables, mEdgeProperties, cannyAOI, mDetectionParameters.scoreThresholdDiff);
    //            vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
    //        }

    //        vEdgePropertiesAll = vEdgePropertiesNew;
    //    }

    //    vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

    ///////////////////////////////////////////////////////////////////////////
    ////////////////////////// EDGE CLASSIFICATION  ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Calculate some edge properties

    { std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, jEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            // re-evaluate edge properties

            mEdgeProperties.length    = calculateEdgeLength(mEdgeProperties.pointIndices, cannyAOI);

            mEdgeProperties.intensity = calculateMeanInt(mEdgeProperties.intensities);
            mEdgeProperties.gradient  = calculateMeanInt(mEdgeProperties.gradients);

            mEdgeProperties.radius    = calculateMean(mEdgeProperties.radii);
            mEdgeProperties.radiusVar = calculateVariance(mEdgeProperties.radii) / mEdgeProperties.radius; // relative variance

            calculateCurvatureStats(mEdgeProperties);

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

    int acceptedFitIndex = 0;

    if (numFits > 0)
    {
        mEllipseProperties.DETECTED = true;

        if (numFits > 1) { acceptedFitIndex = ellipseFitFilter(mDetectionVariables, mDetectionParameters, vEllipsePropertiesAll); } // grab best fit

        mEllipseProperties = vEllipsePropertiesAll[acceptedFitIndex];

        // Tag fits

        for (int iFit = 0; iFit < numFits; iFit++)
        {
            if (iFit == acceptedFitIndex) { vEllipsePropertiesAll[iFit].tag = 1; }
            else                          { vEllipsePropertiesAll[iFit].tag = 0; }
        }

        // Tag edges

        for (int iEdge = 0, numEdges = mEllipseProperties.edgeIndices.size(); iEdge < numEdges; iEdge++)
        {
            int jEdge = mEllipseProperties.edgeIndices[iEdge];
            vEdgePropertiesAll[jEdge].tag = 2;
        }

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
            for (iEdge = 0, numEdges = intensities.size(); iEdge < numEdges; iEdge++)
            {
                double weight  = numEdgePoints[iEdge] / numEdgePointsTotal_1;
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
            for (iEdge = 0, numEdges = curvatures.size(); iEdge < numEdges; iEdge++)
            {
                double weight = numEdgePoints[iEdge] / numEdgePointsTotal_2;
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

    // Logistic functions are used to add response latency to changes in certainty.
    // Multiple detections are required until a reduction in certainty leads to changes in parameter limits or variable predicteds

    double certaintyFactorPosition = 1 - 1 / (1 + exp(-certaintyLatency * mDetectionVariablesNew.certaintyPosition));
    double certaintyFactorFeatures = 1 - 1 / (1 + exp(-certaintyLatency * mDetectionVariablesNew.certaintyFeatures));
    double certaintyFactorAverages = 1 - 1 / (1 + exp(-certaintyLatency * mDetectionVariablesNew.certaintyAverages));

    // Calculate new threshold limits. Thresholds are harsher with higher certainties (= lower certainty factors)

    int AOISize;
    if (mAOI.wdth > mAOI.hght) { AOISize = mAOI.wdth; }
    else                       { AOISize = mAOI.hght; }

    double maxChangeThresholdAspectRatio   = 1.0 / mDetectionParameters.aspectRatioMin;
    double maxChangeThresholdCircumference = mDetectionParameters.circumferenceMax / mDetectionParameters.circumferenceMin;
    double maxChangeThresholdPosition      = AOISize;
    double maxOffsetCircumference          = mDetectionParameters.circumferenceMax / mDetectionParameters.circumferenceMin;

    double rangeChangeThresholdAspectRatio   = maxChangeThresholdAspectRatio   - mDetectionParameters.thresholdChangeAspectRatio;
    double rangeChangeThresholdCircumference = maxChangeThresholdCircumference - mDetectionParameters.thresholdChangeCircumference;
    double rangeChangeThresholdPosition      = maxChangeThresholdPosition      - mDetectionParameters.thresholdChangePosition;
    double rangeOffsetCircumference          = maxOffsetCircumference          - mDetectionParameters.circumferenceOffset;

    mDetectionVariablesNew.thresholdChangeAspectRatio   = rangeChangeThresholdAspectRatio   * certaintyFactorFeatures + mDetectionParameters.thresholdChangeAspectRatio;
    mDetectionVariablesNew.thresholdChangeCircumference = rangeChangeThresholdCircumference * certaintyFactorFeatures + mDetectionParameters.thresholdChangeCircumference;
    mDetectionVariablesNew.thresholdChangePosition      = rangeChangeThresholdPosition      * certaintyFactorPosition + mDetectionParameters.thresholdChangePosition;
    mDetectionVariablesNew.offsetCircumference          = rangeOffsetCircumference          * certaintyFactorAverages + mDetectionParameters.circumferenceOffset;

    std::vector<double> certaintyFactors = {certaintyFactorAverages, certaintyFactorFeatures, certaintyFactorAverages};
    double certaintyFactorMin = *std::min_element(certaintyFactors.begin(), certaintyFactors.end());
    mDetectionVariablesNew.thresholdScore = (1 - certaintyFactorMin) * mDetectionParameters.scoreThreshold;

    if (!mEllipseProperties.DETECTED) // pupil not detected
    {
        // Running averages

        double meanAspectRatio   = initialAspectRatio;
        double meanCircumference = 0.5 * (mDetectionParameters.circumferenceMax + mDetectionParameters.circumferenceMin);
        double meanWidth         = meanCircumference / M_PI;
        double meanHeight        = meanCircumference / M_PI;
        double meanIntensity     = initialIntensity;
        double meanGradient      = 0;
        double meanCurvature     = calculateCurvatureAverage(meanCircumference,meanAspectRatio);

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

        mDetectionVariablesNew.averageAspectRatio   = mDetectionVariables.averageAspectRatio   + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanAspectRatio   - mDetectionVariables.averageAspectRatio);
        mDetectionVariablesNew.averageCircumference = mDetectionVariables.averageCircumference + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanCircumference - mDetectionVariables.averageCircumference);
        mDetectionVariablesNew.averageWidth         = mDetectionVariables.averageWidth         + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanWidth         - mDetectionVariables.averageWidth);
        mDetectionVariablesNew.averageHeight        = mDetectionVariables.averageHeight        + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanHeight        - mDetectionVariables.averageHeight);
        mDetectionVariablesNew.averageCurvature     = mDetectionVariables.averageCurvature     + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanCurvature     - mDetectionVariables.averageCurvature);
        mDetectionVariablesNew.averageIntensity     = mDetectionVariables.averageIntensity     + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanIntensity     - mDetectionVariables.averageIntensity);
        mDetectionVariablesNew.averageGradient      = mDetectionVariables.averageGradient      + mDetectionParameters.alphaAverages * certaintyFactorAverages * (meanGradient      - mDetectionVariables.averageGradient);

        // Feature predictions should decay to average terms. Certainty term gives some latency.

        mDetectionVariablesNew.predictedAspectRatio   = mDetectionVariables.predictedAspectRatio   + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageAspectRatio   - mDetectionVariables.predictedAspectRatio)   + certaintyFactorFeatures * mDetectionVariables.momentumAspectRatio;
        mDetectionVariablesNew.predictedCircumference = mDetectionVariables.predictedCircumference + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageCircumference - mDetectionVariables.predictedCircumference) + certaintyFactorFeatures * mDetectionVariables.momentumCircumference;
        mDetectionVariablesNew.predictedWidth         = mDetectionVariables.predictedWidth         + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageWidth         - mDetectionVariables.predictedWidth)         + certaintyFactorFeatures * mDetectionVariables.momentumWidth;
        mDetectionVariablesNew.predictedHeight        = mDetectionVariables.predictedHeight        + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageHeight        - mDetectionVariables.predictedHeight)        + certaintyFactorFeatures * mDetectionVariables.momentumHeight;
        mDetectionVariablesNew.predictedCurvature     = mDetectionVariables.predictedCurvature     + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageCurvature     - mDetectionVariables.predictedCurvature)     + certaintyFactorFeatures * mDetectionVariables.momentumCurvature;
        mDetectionVariablesNew.predictedIntensity     = mDetectionVariables.predictedIntensity     + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageIntensity     - mDetectionVariables.predictedIntensity)     + certaintyFactorFeatures * mDetectionVariables.momentumIntensity;
        mDetectionVariablesNew.predictedGradient      = mDetectionVariables.predictedGradient      + mDetectionParameters.alphaFeatures * certaintyFactorFeatures * (mDetectionVariablesNew.averageGradient      - mDetectionVariables.predictedGradient)      + certaintyFactorFeatures * mDetectionVariables.momentumGradient;

        // Position predictions should decay to approximate detection position. Certainty term gives some latency.

        mDetectionVariablesNew.predictedXPos = mDetectionVariables.predictedXPos + mDetectionParameters.alphaPosition * certaintyFactorPosition * (haarAOI.xPos + 0.5 * haarAOI.wdth - mDetectionVariables.predictedXPos) + certaintyFactorPosition * mDetectionVariables.momentumXPos;
        mDetectionVariablesNew.predictedYPos = mDetectionVariables.predictedYPos + mDetectionParameters.alphaPosition * certaintyFactorPosition * (haarAOI.yPos + 0.5 * haarAOI.hght - mDetectionVariables.predictedYPos) + certaintyFactorPosition * mDetectionVariables.momentumYPos;

        // Certainty decays to minimum value

        mDetectionVariablesNew.certaintyPosition = mDetectionVariables.certaintyPosition - mDetectionParameters.alphaCertainty * mDetectionParameters.alphaPosition;
        mDetectionVariablesNew.certaintyFeatures = mDetectionVariables.certaintyFeatures - mDetectionParameters.alphaCertainty * mDetectionParameters.alphaFeatures;
        mDetectionVariablesNew.certaintyAverages = mDetectionVariables.certaintyAverages - mDetectionParameters.alphaCertainty * mDetectionParameters.alphaAverages;
    }
    else // pupil detected
    {
        // Exact values

        mDataVariables.exactAspectRatio   = mEllipseProperties.aspectRatio;
        mDataVariables.exactCircumference = mEllipseProperties.circumference;

        mDataVariables.exactXPos = mEllipseProperties.xPos + cannyAOI.xPos;
        mDataVariables.exactYPos = mEllipseProperties.yPos + cannyAOI.yPos;

        // Delta error

        double errorAspectRatio   = mEllipseProperties.aspectRatio   - mDetectionVariables.predictedAspectRatio;
        double errorCircumference = mEllipseProperties.circumference - mDetectionVariables.predictedCircumference;
        double errorWidth         = mEllipseProperties.width         - mDetectionVariables.predictedWidth;
        double errorHeight        = mEllipseProperties.height        - mDetectionVariables.predictedHeight;
        double errorCurvature     = mEllipseProperties.curvature     - mDetectionVariables.predictedCurvature;
        double errorIntensity     = mEllipseProperties.intensity     - mDetectionVariables.predictedIntensity;
        double errorGradient      = mEllipseProperties.gradient      - mDetectionVariables.predictedGradient;
        double errorXPosition     = mDataVariables.exactXPos         - mDetectionVariables.predictedXPos;
        double errorYPosition     = mDataVariables.exactYPos         - mDetectionVariables.predictedYPos;

        // Determine certainty of current measurement


        double displacement                = sqrt(errorXPosition * errorXPosition + errorYPosition * errorYPosition);
        double relativeChangeAspectRatio   = mEllipseProperties.aspectRatio   / mDetectionVariables.predictedAspectRatio;
        double relativeChangeCircumference = mEllipseProperties.circumference / mDetectionVariables.predictedCircumference;

        if (relativeChangeAspectRatio   < 1) { relativeChangeAspectRatio   = 1 / relativeChangeAspectRatio;   }
        if (relativeChangeCircumference < 1) { relativeChangeCircumference = 1 / relativeChangeCircumference; }

        double certaintyPosition      = calculateCertainty(displacement,                    mDetectionParameters.thresholdChangePosition);
        double certaintyAspectRatio   = calculateCertainty(relativeChangeAspectRatio   - 1, mDetectionParameters.thresholdChangeAspectRatio   - 1); // offset by -1 to give maximum at zero
        double certaintyCircumference = calculateCertainty(relativeChangeCircumference - 1, mDetectionParameters.thresholdChangeCircumference - 1);
        double certaintyFeatures      = 0.5 * (certaintyAspectRatio + certaintyCircumference);

        // Increase or reduce certainty with a fraction of a constant step size (product of the two learning rate parameters)

        mDetectionVariablesNew.certaintyPosition = mDetectionVariables.certaintyPosition + certaintyPosition * mDetectionParameters.alphaCertainty * mDetectionParameters.alphaPosition;
        mDetectionVariablesNew.certaintyFeatures = mDetectionVariables.certaintyFeatures + certaintyFeatures * mDetectionParameters.alphaCertainty * mDetectionParameters.alphaFeatures;
        mDetectionVariablesNew.certaintyAverages = mDetectionVariables.certaintyAverages + certaintyFeatures * mDetectionParameters.alphaCertainty * mDetectionParameters.alphaAverages;

        // Momentum. Rate of change approximation.

        mDetectionVariablesNew.momentumAspectRatio   =  mDetectionVariablesNew.momentumAspectRatio   + mDetectionParameters.alphaFeatures * (errorAspectRatio   - mDetectionVariablesNew.momentumAspectRatio);
        mDetectionVariablesNew.momentumCircumference =  mDetectionVariablesNew.momentumCircumference + mDetectionParameters.alphaFeatures * (errorCircumference - mDetectionVariablesNew.momentumCircumference);
        mDetectionVariablesNew.momentumWidth         =  mDetectionVariablesNew.momentumWidth         + mDetectionParameters.alphaFeatures * (errorWidth         - mDetectionVariablesNew.momentumWidth);
        mDetectionVariablesNew.momentumHeight        =  mDetectionVariablesNew.momentumHeight        + mDetectionParameters.alphaFeatures * (errorHeight        - mDetectionVariablesNew.momentumHeight);
        mDetectionVariablesNew.momentumGradient      =  mDetectionVariablesNew.momentumGradient      + mDetectionParameters.alphaFeatures * (errorGradient      - mDetectionVariablesNew.momentumGradient);
        mDetectionVariablesNew.momentumIntensity     =  mDetectionVariablesNew.momentumIntensity     + mDetectionParameters.alphaFeatures * (errorIntensity     - mDetectionVariablesNew.momentumIntensity);
        mDetectionVariablesNew.momentumCurvature     =  mDetectionVariablesNew.momentumCurvature     + mDetectionParameters.alphaFeatures * (errorCurvature     - mDetectionVariablesNew.momentumCurvature);
        mDetectionVariablesNew.momentumXPos          =  mDetectionVariablesNew.momentumXPos          + mDetectionParameters.alphaPosition * (errorXPosition     - mDetectionVariablesNew.momentumXPos);
        mDetectionVariablesNew.momentumYPos          =  mDetectionVariablesNew.momentumYPos          + mDetectionParameters.alphaPosition * (errorYPosition     - mDetectionVariablesNew.momentumYPos);

        // Update predictions. Only add momentum when certainty is high

        mDetectionVariablesNew.predictedAspectRatio   = mDetectionVariables.predictedAspectRatio   + mDetectionParameters.alphaFeatures * errorAspectRatio   + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumAspectRatio;
        mDetectionVariablesNew.predictedCircumference = mDetectionVariables.predictedCircumference + mDetectionParameters.alphaFeatures * errorCircumference + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumCircumference;
        mDetectionVariablesNew.predictedWidth         = mDetectionVariables.predictedWidth         + mDetectionParameters.alphaFeatures * errorWidth         + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumWidth;
        mDetectionVariablesNew.predictedHeight        = mDetectionVariables.predictedHeight        + mDetectionParameters.alphaFeatures * errorHeight        + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumHeight;
        mDetectionVariablesNew.predictedCurvature     = mDetectionVariables.predictedCurvature     + mDetectionParameters.alphaFeatures * errorCurvature     + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumCurvature;
        mDetectionVariablesNew.predictedIntensity     = mDetectionVariables.predictedIntensity     + mDetectionParameters.alphaFeatures * errorIntensity     + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumIntensity;
        mDetectionVariablesNew.predictedGradient      = mDetectionVariables.predictedGradient      + mDetectionParameters.alphaFeatures * errorGradient      + (1 - certaintyFactorFeatures) * mDetectionVariables.momentumGradient;
        mDetectionVariablesNew.predictedXPos          = mDetectionVariables.predictedXPos          + mDetectionParameters.alphaPosition * errorXPosition     + (1 - certaintyFactorPosition) * mDetectionVariables.momentumXPos;
        mDetectionVariablesNew.predictedYPos          = mDetectionVariables.predictedYPos          + mDetectionParameters.alphaPosition * errorYPosition     + (1 - certaintyFactorPosition) * mDetectionVariables.momentumYPos;

        // Averages

        mDetectionVariablesNew.averageAspectRatio   = mDetectionVariables.averageAspectRatio   + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedAspectRatio   - mDetectionVariables.averageAspectRatio);
        mDetectionVariablesNew.averageCircumference = mDetectionVariables.averageCircumference + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedCircumference - mDetectionVariables.averageCircumference);
        mDetectionVariablesNew.averageWidth         = mDetectionVariables.averageWidth         + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedWidth         - mDetectionVariables.averageWidth);
        mDetectionVariablesNew.averageHeight        = mDetectionVariables.averageHeight        + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedHeight        - mDetectionVariables.averageHeight);
        mDetectionVariablesNew.averageCurvature     = mDetectionVariables.averageCurvature     + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedCurvature     - mDetectionVariables.averageCurvature);
        mDetectionVariablesNew.averageIntensity     = mDetectionVariables.averageIntensity     + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedIntensity     - mDetectionVariables.averageIntensity);
        mDetectionVariablesNew.averageGradient      = mDetectionVariables.averageGradient      + mDetectionParameters.alphaAverages * (mDetectionVariables.predictedGradient      - mDetectionVariables.averageGradient);
    }

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
