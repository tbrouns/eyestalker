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

    double factorLength = scoreFactorLength * (edgeLength / mDetectionVariables.predictedCircumference) + (1 - scoreFactorLength);

    double factorRadiusVar = certaintyFactorPosition  * scoreFactorRadiusVar * factorLength;
    double factorCurvature = certaintyFactorCurvature * scoreFactorCurvature * factorLength;

    // Calculate scores

    double scoreRadius        = factorRadius        * calculateGaussian2(changeRadius,        parametersEdgeRadius       );
    double scoreCircumference = factorCircumference * calculateGaussian2(changeCircumference, parametersEdgeCircumference);
    double scoreCurvature     = factorCurvature     * calculateGaussian2(changeCurvature,     parametersEdgeCurvature    );
    double scoreIntensity     = factorIntensity     * calculateGaussian2(changeIntensity,     parametersEdgeIntensity    );
    double scoreGradient      = factorGradient      * calculateGaussian2(changeGradient,      parametersEdgeGradient     );
    double scoreRadiusVar     = factorRadiusVar     * calculateGaussian2(varianceRadius,      parametersEdgeRadiusVar    );

    const double norm =  factorRadius + factorRadiusVar + factorCurvature + factorCircumference + factorIntensity + factorGradient;

    double scoreTotal = 0;
    if (norm > 0) { scoreTotal = (scoreRadius + scoreRadiusVar + scoreGradient + scoreCurvature + scoreCircumference + scoreIntensity) / norm; }
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

                if (x < 0 || x >= mAOI.wdth || y < 0 || y >= mAOI.hght)
                {
                    STOP_SEARCH = true;
                    break;
                }

                int centreIndex = y * mAOI.wdth + x;

                int tag = cannyEdgeVector[centreIndex];

                if (tag > 0)
                {
                    if (tag == 1)
                    {
                        startIndices.push_back(centreIndex);
                        if (R > pupilRadius) { STOP_SEARCH = true; }
                    }
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

int connectEdges(const detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1};
    std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1};

    std::vector<int> edgePoints = {startIndex};
    int edgePointNew = {startIndex};

    for (int iEdgePoint = 0; iEdgePoint < mDetectionParameters.windowLengthCurvature; iEdgePoint++) // move back through edge
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

    if (edgeLength >= mDetectionParameters.windowLengthCurvature)
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
        circumferences[0] = mDetectionVariables.averageCircumference * mDetectionVariables.offsetCircumference;
        circumferences[1] = mDetectionVariables.averageCircumference / mDetectionVariables.offsetCircumference;
    }
    else
    {
        circumferences[0] = mDetectionVariables.predictedCircumference * (1 / mDetectionVariables.thresholdChangeCircumference);
        circumferences[1] = mDetectionVariables.predictedCircumference *      mDetectionVariables.thresholdChangeCircumference;
    }

    aspectRatios[0] = mDetectionVariables.predictedAspectRatio * (1 / mDetectionVariables.thresholdChangeAspectRatio);
    aspectRatios[1] = mDetectionVariables.predictedAspectRatio *      mDetectionVariables.thresholdChangeAspectRatio;

    // Calculate limits

    std::vector<double> curvaturesMax(4);
    std::vector<double> curvaturesMin(4);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            curvaturesMax[2 * j + i] = getCurvatureUpperLimit(circumferences[i], aspectRatios[j], mDetectionParameters.windowLengthCurvature);
            curvaturesMin[2 * j + i] = getCurvatureLowerLimit(circumferences[i], aspectRatios[j], mDetectionParameters.windowLengthCurvature);
        }
    }

    curvatureUpperLimit = *std::max_element(std::begin(curvaturesMax), std::end(curvaturesMax)) + mDetectionParameters.curvatureOffset;
    curvatureLowerLimit = *std::min_element(std::begin(curvaturesMin), std::end(curvaturesMin)) - mDetectionParameters.curvatureOffset;
}

std::vector<double> calculateCurvatures(const detectionParameters& mDetectionParameters, std::vector<double>& xNormals, std::vector<double>& yNormals, const std::vector<double>& xTangentsAll, const std::vector<double>& yTangentsAll)
{
    int edgeSize = xTangentsAll.size();

    std::vector<double> curvatures(edgeSize, 0.0);

    int numPos = 0;
    int numNeg = 0;

    for (int iEdgePoint = mDetectionParameters.windowLengthCurvature; iEdgePoint < edgeSize - mDetectionParameters.windowLengthCurvature; iEdgePoint++)
    {
        // calculate window tangents

        // first window

        std::vector<double> tangentsX_1(xTangentsAll.begin() + iEdgePoint - mDetectionParameters.windowLengthCurvature, xTangentsAll.begin() + iEdgePoint);
        std::vector<double> tangentsY_1(yTangentsAll.begin() + iEdgePoint - mDetectionParameters.windowLengthCurvature, yTangentsAll.begin() + iEdgePoint);

        double tangentXMean_1 = calculateMean(tangentsX_1);
        double tangentYMean_1 = calculateMean(tangentsY_1);

        // second window

        std::vector<double> tangentsX_2(xTangentsAll.begin() + iEdgePoint + 1, xTangentsAll.begin() + iEdgePoint + 1 + mDetectionParameters.windowLengthCurvature);
        std::vector<double> tangentsY_2(yTangentsAll.begin() + iEdgePoint + 1, yTangentsAll.begin() + iEdgePoint + 1 + mDetectionParameters.windowLengthCurvature);

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
        for (int iEdgePoint = mDetectionParameters.windowLengthCurvature; iEdgePoint < edgeSize - mDetectionParameters.windowLengthCurvature; iEdgePoint++)
        {
            curvatures[iEdgePoint] = -curvatures[iEdgePoint];
        }
    }

    return curvatures;
}

std::vector<edgeProperties> edgeSegmentationCurvature(const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties, const double curvatureUpperLimit, const double curvatureLowerLimit)
{
    int edgeSize = mEdgeProperties.curvatures.size();

    // find breakpoints based on curvature thresholding

    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(-1); // add first point (+ 1 is added later)

    for (int iEdgePoint = mDetectionParameters.windowLengthCurvature; iEdgePoint < edgeSize - mDetectionParameters.windowLengthCurvature; iEdgePoint++)
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

        if (diffLeft > mDetectionParameters.windowLengthCurvature || diffRght > mDetectionParameters.windowLengthCurvature)
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

double calculateEdgeLength(const detectionParameters& mDetectionParameters, const std::vector<int>& edgePoints, const AOIProperties& mAOI)
{
    // calculates edge length without need of edge continuity

    double lengthTotal = 0;
    int edgeSize = edgePoints.size();

    int iEdgePoint  = 0;
    bool BREAK_LOOP = false;

    do // run through all edge points
    {
        int jEdgePoint = iEdgePoint + mDetectionParameters.windowLengthEdge;

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

        iEdgePoint = iEdgePoint + mDetectionParameters.windowLengthEdge;
    } while (!BREAK_LOOP);

    return lengthTotal;
}

void calculateCurvatureStats(const detectionParameters& mDetectionParameters, edgeProperties& mEdgeProperties)
{
    // Calculate min, max and mean curvature

    std::vector<double> edgeCurvaturesNew;

    int edgeSize = mEdgeProperties.curvatures.size();

    for (int iEdgePoint = mDetectionParameters.windowLengthCurvature; iEdgePoint < edgeSize - mDetectionParameters.windowLengthCurvature; iEdgePoint++)
    {
        edgeCurvaturesNew.push_back(mEdgeProperties.curvatures[iEdgePoint]);
    }

    int edgeLengthNew = edgeCurvaturesNew.size();

    double curvature    = 360;
    double curvatureMax = 360;
    double curvatureMin = 360;

    if (edgeLengthNew > 0)
    {
        curvature    = calculateMean(edgeCurvaturesNew);
        curvatureMax = *std::max_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
        curvatureMin = *std::min_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
    }

    mEdgeProperties.curvature    = curvature;
    mEdgeProperties.curvatureMax = curvatureMax;
    mEdgeProperties.curvatureMin = curvatureMin;
}

std::vector<edgeProperties> edgeSegmentationLength(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties, const AOIProperties& mAOI)
{
    // This functions cuts edge terminals to make the edge shorter, if the edge is significantly longer than predicted

    int edgeSize = mEdgeProperties.curvatures.size();

    // find breakpoints based on length thresholding

    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(-1); // add first point (+ 1 is added later)

    double lengthDifference = mEdgeProperties.length - mDetectionVariables.predictedCircumference;

    if (lengthDifference > mDetectionParameters.fitEdgeFraction * mDetectionVariables.predictedCircumference) // difference should be large enough
    {
        if (edgeSize > 2 * lengthDifference + mDetectionParameters.windowLengthCurvature) // should be enough space between breakpoints
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
            mEdgePropertiesNew.length    = calculateEdgeLength(mDetectionParameters, mEdgePropertiesNew.pointIndices, mAOI);

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

std::vector<edgeProperties> edgeTerminalFilter(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties, const AOIProperties& mAOI, double scoreThresholdDiff)
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

            if (edgeSize <= 3 * mDetectionParameters.windowLengthCurvature) // don't segment short edges
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
                iEdgePoint = iEdgePoint + mDetectionParameters.windowLengthCurvature;

                if (iEdgePoint >= edgeSize - mDetectionParameters.windowLengthCurvature - 1)
                {
                    BREAK_LOOP = true;
                    iEdgePoint = edgeSize - mDetectionParameters.windowLengthCurvature - 1;
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
                    mEdgePropertiesNew.length    = calculateEdgeLength(mDetectionParameters, mEdgePropertiesNew.pointIndices, mAOI);

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

            if (scoreDifferenceMax > scoreThresholdDiff)
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
        if (edgeSize >= mDetectionParameters.windowLengthCurvature) { vEdgePropertiesNew.push_back(mEdgeProperties); }
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

    if (mDetectionVariables.certaintyAverages > 0 || mDetectionVariables.certaintyFeatures > 0 || mDetectionVariables.certaintyPosition > 0) // need to have enough certainty
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

            ellipseProperties mEllipseProperties = fitEllipse(edgeIndices, edgeSetSize, mAOI.wdth);

            if (!mEllipseProperties.DETECTED) { continue; } // error

            // Absolute displacement filter

            double dX = mEllipseProperties.xPos - (mDetectionVariables.predictedXPos - mAOI.xPos);
            double dY = mEllipseProperties.yPos - (mDetectionVariables.predictedYPos - mAOI.yPos);
            double dR = sqrt(dX * dX + dY * dY);
            if (dR > mDetectionVariables.thresholdChangePosition) { continue; } // no large ellipse displacements

            // Absolute size and shape filter

            if (mEllipseProperties.circumference > circumferenceUpperLimit) { continue; } // no large ellipse
            if (mEllipseProperties.circumference < circumferenceLowerLimit) { continue; } // no small ellipse
            if (mEllipseProperties.aspectRatio   < mDetectionParameters.aspectRatioMin)   { continue; } // no extreme deviations from circular shape

            // Relative change in size and shape filter

            double relativeChangeAspectRatio   = mEllipseProperties.aspectRatio   / mDetectionVariables.predictedAspectRatio;
            double relativeChangeCircumference = mEllipseProperties.circumference / mDetectionVariables.predictedCircumference;

            if (relativeChangeAspectRatio   < 1.0) { relativeChangeAspectRatio   = 1 / relativeChangeAspectRatio; }
            if (relativeChangeCircumference < 1.0) { relativeChangeCircumference = 1 / relativeChangeCircumference; }

            if (relativeChangeAspectRatio   > mDetectionVariables.thresholdChangeAspectRatio  ) { continue; } // no large ellipse shape changes
            if (relativeChangeCircumference > mDetectionVariables.thresholdChangeCircumference) { continue; } // no large ellipse size changes

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

int ellipseFitFilter(const detectionVariables& mDetectionVariables, std::vector<ellipseProperties> vEllipseProperties)
{
    static const double scoreFactorFitAspectRatio   = 0.41;
    static const double scoreFactorFitCircumference = 0.98;
    static const double scoreFactorFitDisplacement  = 0.00;
    static const double scoreFactorFitLength        = 3.71;
    static const double scoreFactorFitError         = 0.18;

    static const std::vector<double> parametersFitAspectRatio   = {0.73173, 0.0005102, 0.0093987, 0.23184, 0.0005102, 0.033331};
    static const std::vector<double> parametersFitCircumference = {0.81525,   0.99994, 0.0051646, 0.17624,    1.0089, 0.034909};
    static const std::vector<double> parametersFitDisplacement  = {0.37568,  0.040816,   0.66938, 0.59611,  0.040816,   3.1728};
    static const std::vector<double> parametersFitLength        = {0.12018,   0.76041,  0.033629, 0.97603,   0.99961,  0.26897};
    static const std::vector<double> parametersFitError         = {0.83881,   0.20388,  0.075325, 0.28425,   0.41703,  0.24919};

    int numFits = vEllipseProperties.size();

    std::vector<double> scoreFits(numFits);

    double certaintyFactorPosition = 0.5 * (mDetectionVariables.certaintyPosition + 1);
    double certaintyFactorFeatures = 0.5 * (mDetectionVariables.certaintyFeatures + 1);
    double certaintyFactorAverages = 0.5 * (mDetectionVariables.certaintyAverages + 1);

    if (certaintyFactorAverages > certaintyFactorFeatures + certaintyOffset) { certaintyFactorFeatures = certaintyFactorAverages; }

    for (int iFit = 0; iFit < numFits; iFit++)
    {
        ellipseProperties mEllipseProperties = vEllipseProperties[iFit];

        double dx = mEllipseProperties.xPos - mDetectionVariables.predictedXPos;
        double dy = mEllipseProperties.yPos - mDetectionVariables.predictedYPos;
        double changePosition      = sqrt(dx * dx + dy * dy);
        double changeAspectRatio   = std::abs(mEllipseProperties.aspectRatio - mDetectionVariables.predictedAspectRatio);
        double changeCircumference = mEllipseProperties.circumference / mDetectionVariables.predictedCircumference;
        double changeLength        = mEllipseProperties.edgeLength    / mDetectionVariables.predictedCircumference;
        double fitError            = mEllipseProperties.fitError;

        double factorAspectRatio   = certaintyFactorFeatures * scoreFactorFitAspectRatio;
        double factorCircumference = certaintyFactorFeatures * scoreFactorFitCircumference;
        double factorDisplacement  = certaintyFactorPosition * scoreFactorFitDisplacement;
        double factorLength        = certaintyFactorFeatures * scoreFactorFitLength;
        double factorFitError      =                           scoreFactorFitError;

        // Calculate scores
        double scoreAspectRatio   = factorAspectRatio   * calculateGaussian2(changeAspectRatio,   parametersFitAspectRatio   );
        double scoreCircumference = factorCircumference * calculateGaussian2(changeCircumference, parametersFitCircumference );
        double scoreDisplacement  = factorDisplacement  * calculateGaussian2(changePosition,      parametersFitDisplacement  );
        double scoreLength        = factorLength        * calculateGaussian2(changeLength,        parametersFitLength        );
        double scoreFitError      = factorFitError      * calculateGaussian2(fitError,            parametersFitError         );

        double norm =  factorAspectRatio + factorCircumference + factorDisplacement + factorLength + factorFitError;

        double scoreTotal = 0;
        if (norm > 0) { scoreTotal = (scoreAspectRatio + scoreCircumference + scoreDisplacement + scoreLength + scoreFitError) / norm; }

        scoreFits[iFit] = scoreTotal;
    }

    int acceptedFitIndex = std::distance(scoreFits.begin(), std::max_element(scoreFits.begin(), scoreFits.end()));

    return acceptedFitIndex;
}

void getWindowLengthEdge(detectionParameters& mDetectionParameters)
{
    std::vector<int> windowLengthCurvatures = {4, 7, 11};

    // Functions finds closest discrete approximation to the desired exact window length

    int numWindowLengths = windowLengthCurvatures.size();

    std::vector<double> diffs(numWindowLengths);

    for (int iLength = 0; iLength < numWindowLengths; iLength++)
    {
        diffs[iLength] = std::abs(windowLengthCurvatures[iLength] - mDetectionParameters.windowLengthEdge);
    }

    int windowLengthIndex = std::distance(diffs.begin(), std::min_element(diffs.begin(), diffs.end()));

    mDetectionParameters.windowLengthCurvature = windowLengthCurvatures[windowLengthIndex]; // closest approximation
}

void setCurvatureMeasurement(detectionParameters& mDetectionParameters, int imgWdth)
{
    // For development

    mDetectionParameters.windowLengthCurvature = mDetectionParameters.windowLengthEdge;

    mDetectionParameters.aspectRatioMin   = 0.0;
    mDetectionParameters.circumferenceMax = M_PI * imgWdth;
    mDetectionParameters.circumferenceMin = 1.0; // 1.0 avoids inf

    mDetectionParameters.curvatureOffset = 360;

    mDetectionParameters.scoreThreshold     = 0.0;
    mDetectionParameters.scoreThresholdDiff = 1.0;

    mDetectionParameters.glintWdth = 0.0;
}

void checkVariableLimits(detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters)
{
    if (mDetectionVariables.thresholdScore > mDetectionParameters.scoreThreshold)
    {   mDetectionVariables.thresholdScore = mDetectionParameters.scoreThreshold; }

    if (mDetectionVariables.thresholdChangeAspectRatio < mDetectionParameters.thresholdChangeAspectRatio)
    {   mDetectionVariables.thresholdChangeAspectRatio = mDetectionParameters.thresholdChangeAspectRatio; }

    if (mDetectionVariables.thresholdChangeCircumference < mDetectionParameters.thresholdChangeCircumference)
    {   mDetectionVariables.thresholdChangeCircumference = mDetectionParameters.thresholdChangeCircumference; }

    if (mDetectionVariables.thresholdChangePosition < mDetectionParameters.thresholdChangePosition)
    {   mDetectionVariables.thresholdChangePosition = mDetectionParameters.thresholdChangePosition; }

    if      (mDetectionVariables.certaintyPosition < -1.0) { mDetectionVariables.certaintyPosition = -1.0; }
    else if (mDetectionVariables.certaintyPosition >  1.0) { mDetectionVariables.certaintyPosition =  1.0; }

    if      (mDetectionVariables.certaintyFeatures < -1.0) { mDetectionVariables.certaintyFeatures = -1.0; }
    else if (mDetectionVariables.certaintyFeatures >  1.0) { mDetectionVariables.certaintyFeatures =  1.0; }

    if      (mDetectionVariables.certaintyAverages < -1.0) { mDetectionVariables.certaintyAverages = -1.0; }
    else if (mDetectionVariables.certaintyAverages >  1.0) { mDetectionVariables.certaintyAverages =  1.0; }
}

detectionVariables eyeStalker(const cv::Mat& imageOriginalBGR, const AOIProperties& mAOI, detectionVariables& mDetectionVariables, detectionParameters& mDetectionParameters, dataVariables& mDataVariables, drawVariables& mDrawVariables, const developmentOptions& mAdvancedOptions)
{
    int imageWdth = imageOriginalBGR.cols;

    checkVariableLimits(mDetectionVariables, mDetectionParameters); // keep variables within limits

    if (mAdvancedOptions.CURVATURE_MEASUREMENT) { setCurvatureMeasurement(mDetectionParameters, imageWdth); }
    else                                           { getWindowLengthEdge(mDetectionParameters); }

    // Define search area

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

        mEdgeProperties.length      = calculateEdgeLength(mDetectionParameters, mEdgeProperties.pointIndices, cannyAOI);
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
            std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationLength(mDetectionVariables, mDetectionParameters, vEdgePropertiesAll[iEdge], cannyAOI);
            vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end());
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
    }

    vEdgePropertiesAll = removeShortEdges(mDetectionParameters, vEdgePropertiesAll);

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// EDGE TERMINAL FILTER   ////////////////////////
    ////////////////////////////////////////////////////////////////////////

    { std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties    = vEdgePropertiesAll[iEdge];
            std::vector<edgeProperties> vEdgePropertiesTemp = edgeTerminalFilter(mDetectionVariables, mDetectionParameters, mEdgeProperties, cannyAOI, mDetectionParameters.scoreThresholdDiff);
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

            mEdgeProperties.length    = calculateEdgeLength(mDetectionParameters, mEdgeProperties.pointIndices, cannyAOI);

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

    int acceptedFitIndex = 0;

    if (numFits > 0)
    {
        mEllipseProperties.DETECTED = true;

        if (numFits > 1) { acceptedFitIndex = ellipseFitFilter(mDetectionVariables, vEllipsePropertiesAll); } // grab best fit

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

        double curvatureMax = getCurvatureUpperLimit(meanCircumference, meanAspectRatio, mDetectionParameters.windowLengthCurvature);
        double curvatureMin = getCurvatureLowerLimit(meanCircumference, meanAspectRatio, mDetectionParameters.windowLengthCurvature);

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

    case 4:
    {
        static const std::vector<double> arrayCurvatureMax =    {111.8,108.2,104.5,100.5,96.5,92.5,88.4,84.4,80.6,76.9,73.5,70.3,67.4,64.7,62.3,60.1,58.2,56.5,55.1,53.8,52.7,51.7,50.9,50.2,49.6,49  ,48.6,48.2,47.8,47.5,47.2,47  ,46.8,46.6,46.4,46.3,46.1,46  ,45.9,45.7,45.6,45.5,45.3,45.2,45  ,44.9,44.7,44.5,44.4,44.2,
                                                                 110.2,106.6,102.8,98.8 ,94.8,90.7,86.7,82.8,79  ,75.5,72.1,69  ,66.2,63.6,61.3,59.2,57.4,55.8,54.4,53.2,52.2,51.3,50.5,49.8,49.2,48.7,48.3,47.9,47.6,47.3,47.1,46.9,46.7,46.5,46.3,46.2,46  ,45.9,45.8,45.6,45.5,45.3,45.2,45.1,44.9,44.7,44.6,44.4,44.2,44  ,
                                                                 108.5,104.8,100.9,97   ,92.9,88.9,84.9,81.1,77.4,73.9,70.6,67.6,64.9,62.4,60.2,58.3,56.6,55.1,53.7,52.6,51.6,50.8,50  ,49.4,48.9,48.4,48  ,47.6,47.3,47.1,46.8,46.6,46.5,46.3,46.1,46  ,45.8,45.7,45.6,45.5,45.3,45.2,45  ,44.9,44.7,44.6,44.4,44.2,44  ,43.8,
                                                                 106.6,102.8,99   ,95   ,91  ,87  ,83  ,79.3,75.6,72.2,69.1,66.2,63.6,61.2,59.1,57.2,55.6,54.2,52.9,51.9,50.9,50.1,49.5,48.9,48.4,47.9,47.6,47.2,47  ,46.7,46.5,46.3,46.2,46  ,45.9,45.7,45.6,45.5,45.3,45.2,45.1,44.9,44.8,44.6,44.5,44.3,44.1,43.9,43.7,43.5,
                                                                 104.6,100.8,96.9 ,92.9 ,88.9,84.9,81  ,77.3,73.8,70.5,67.4,64.6,62.1,59.9,57.8,56.1,54.5,53.2,52  ,51  ,50.1,49.4,48.8,48.2,47.8,47.4,47  ,46.7,46.5,46.3,46.1,45.9,45.7,45.6,45.5,45.3,45.2,45.1,45  ,44.9,44.7,44.6,44.4,44.3,44.1,43.9,43.7,43.5,43.3,43.1,
                                                                 102.5,98.7 ,94.7 ,90.7 ,86.7,82.8,79  ,75.3,71.9,68.7,65.7,63  ,60.6,58.4,56.5,54.8,53.4,52.1,51  ,50  ,49.2,48.6,48  ,47.5,47  ,46.7,46.4,46.1,45.9,45.7,45.5,45.3,45.2,45.1,45  ,44.8,44.7,44.6,44.5,44.4,44.3,44.1,44  ,43.8,43.6,43.5,43.3,43  ,42.8,42.6,
                                                                 100.4,96.6 ,92.6 ,88.6 ,84.6,80.7,76.9,73.3,69.9,66.8,64  ,61.4,59  ,57  ,55.1,53.5,52.1,50.9,49.9,49  ,48.2,47.6,47.1,46.6,46.2,45.9,45.6,45.3,45.1,45  ,44.8,44.7,44.5,44.4,44.3,44.2,44.1,44  ,43.9,43.8,43.6,43.5,43.4,43.2,43  ,42.8,42.6,42.4,42.2,42  ,
                                                                 98.4 ,94.5 ,90.5 ,86.5 ,82.5,78.6,74.9,71.4,68.1,65  ,62.2,59.7,57.5,55.5,53.7,52.2,50.9,49.7,48.8,47.9,47.2,46.6,46.1,45.7,45.3,45  ,44.7,44.5,44.3,44.1,44  ,43.9,43.8,43.7,43.6,43.5,43.4,43.2,43.1,43  ,42.9,42.8,42.6,42.5,42.3,42.1,41.9,41.7,41.4,41.2,
                                                                 96.5 ,92.5 ,88.5 ,84.5 ,80.5,76.7,73  ,69.5,66.3,63.3,60.6,58.2,56  ,54.1,52.4,51  ,49.7,48.6,47.7,46.9,46.2,45.6,45.1,44.7,44.3,44  ,43.8,43.6,43.4,43.2,43.1,43  ,42.9,42.8,42.7,42.6,42.5,42.4,42.3,42.2,42  ,41.9,41.7,41.6,41.4,41.2,41  ,40.8,40.6,40.3,
                                                                 94.7 ,90.7 ,86.6 ,82.6 ,78.7,74.9,71.2,67.8,64.7,61.8,59.1,56.8,54.7,52.8,51.2,49.8,48.6,47.5,46.6,45.9,45.2,44.6,44.2,43.8,43.4,43.1,42.9,42.7,42.5,42.3,42.2,42.1,42  ,41.9,41.8,41.7,41.5,41.4,41.3,41.2,41.1,40.9,40.7,40.6,40.4,40.2,40  ,39.8,39.5,39.3,
                                                                 93   ,89   ,84.9 ,80.9 ,77  ,73.2,69.6,66.3,63.2,60.4,57.8,55.5,53.5,51.7,50.2,48.8,47.6,46.6,45.7,45  ,44.3,43.8,43.3,42.9,42.6,42.3,42.1,41.9,41.7,41.5,41.3,41.2,41.1,41  ,40.8,40.7,40.6,40.5,40.3,40.2,40  ,39.9,39.7,39.5,39.3,39.1,38.9,38.7,38.4,38.2,
                                                                 91.4 ,87.4 ,83.3 ,79.3 ,75.4,71.7,68.2,64.9,61.9,59.1,56.7,54.4,52.5,50.8,49.2,47.9,46.8,45.8,45  ,44.2,43.6,43.1,42.6,42.2,41.9,41.6,41.3,41.1,40.9,40.7,40.6,40.4,40.3,40.1,40  ,39.8,39.7,39.5,39.4,39.2,39  ,38.9,38.7,38.5,38.2,38  ,37.8,37.5,37.3,37.1,
                                                                 89.9 ,85.8 ,81.8 ,77.8 ,73.9,70.3,66.8,63.6,60.7,58  ,55.6,53.5,51.6,49.9,48.5,47.2,46.1,45.2,44.3,43.6,43  ,42.5,42  ,41.6,41.3,41  ,40.7,40.5,40.3,40.1,39.9,39.7,39.6,39.4,39.2,39.1,38.9,38.7,38.5,38.3,38.1,37.9,37.7,37.4,37.2,37  ,36.7,36.5,36.2,36  ,
                                                                 88.5 ,84.4 ,80.3 ,76.4 ,72.6,69  ,65.6,62.5,59.6,57  ,54.7,52.7,50.8,49.2,47.8,46.6,45.5,44.6,43.8,43.1,42.5,42  ,41.6,41.2,40.8,40.5,40.3,40  ,39.8,39.6,39.4,39.2,39  ,38.8,38.6,38.4,38.2,38  ,37.8,37.5,37.3,37.1,36.8,36.5,36.3,36  ,35.7,35.5,35.2,35  ,
                                                                 87.1 ,83   ,79   ,75.1 ,71.3,67.8,64.5,61.4,58.7,56.2,53.9,51.9,50.2,48.6,47.3,46.1,45.1,44.2,43.4,42.8,42.2,41.7,41.2,40.8,40.5,40.2,39.9,39.6,39.4,39.2,38.9,38.7,38.5,38.3,38.1,37.8,37.6,37.4,37.1,36.8,36.6,36.3,36  ,35.8,35.5,35.2,34.9,34.6,34.3,34.1,
                                                                 85.8 ,81.7 ,77.7 ,73.8 ,70.1,66.7,63.4,60.5,57.8,55.4,53.2,51.3,49.6,48.1,46.8,45.7,44.7,43.9,43.1,42.5,41.9,41.4,41  ,40.6,40.2,39.9,39.6,39.3,39.1,38.8,38.6,38.3,38.1,37.9,37.6,37.4,37.1,36.8,36.6,36.3,36  ,35.7,35.4,35.1,34.8,34.5,34.2,33.9,33.6,33.3,
                                                                 84.4 ,80.4 ,76.4 ,72.6 ,69  ,65.6,62.5,59.6,57  ,54.6,52.6,50.7,49.1,47.7,46.4,45.3,44.4,43.6,42.8,42.2,41.7,41.2,40.7,40.4,40  ,39.7,39.4,39.1,38.8,38.6,38.3,38  ,37.8,37.5,37.2,37  ,36.7,36.4,36.1,35.8,35.5,35.2,34.8,34.5,34.2,33.9,33.6,33.3,33  ,32.7,
                                                                 83.1 ,79.1 ,75.2 ,71.4 ,67.9,64.6,61.5,58.7,56.2,54  ,52  ,50.2,48.6,47.3,46.1,45  ,44.1,43.3,42.6,42  ,41.5,41  ,40.6,40.2,39.8,39.5,39.2,38.9,38.6,38.3,38.1,37.8,37.5,37.2,36.9,36.6,36.3,36  ,35.7,35.4,35  ,34.7,34.4,34.1,33.7,33.4,33.1,32.8,32.5,32.2,
                                                                 81.9 ,77.9 ,74   ,70.3 ,66.8,63.6,60.6,58  ,55.5,53.4,51.4,49.7,48.2,46.9,45.8,44.8,43.9,43.1,42.4,41.8,41.3,40.9,40.4,40  ,39.7,39.3,39  ,38.7,38.4,38.1,37.8,37.5,37.2,36.9,36.6,36.3,36  ,35.7,35.3,35  ,34.7,34.3,34  ,33.7,33.3,33  ,32.7,32.4,32.1,31.8,
                                                                 80.6 ,76.7 ,72.8 ,69.2 ,65.8,62.7,59.8,57.2,54.9,52.8,50.9,49.3,47.8,46.6,45.5,44.5,43.7,42.9,42.3,41.7,41.2,40.7,40.3,39.9,39.6,39.2,38.9,38.6,38.3,38  ,37.6,37.3,37  ,36.7,36.4,36  ,35.7,35.4,35  ,34.7,34.3,34  ,33.6,33.3,33  ,32.7,32.4,32.1,31.7,31.4,
                                                                 79.4 ,75.5 ,71.7 ,68.2 ,64.9,61.8,59  ,56.5,54.2,52.2,50.4,48.9,47.5,46.3,45.2,44.3,43.5,42.8,42.1,41.6,41.1,40.6,40.2,39.8,39.4,39.1,38.8,38.4,38.1,37.8,37.5,37.1,36.8,36.5,36.1,35.8,35.4,35.1,34.7,34.4,34  ,33.7,33.4,33  ,32.7,32.4,32.1,31.8,31.4,31.1,
                                                                 78.2 ,74.3 ,70.6 ,67.1 ,63.9,60.9,58.2,55.8,53.6,51.7,50  ,48.5,47.2,46  ,45  ,44.1,43.3,42.6,42  ,41.4,41  ,40.5,40.1,39.7,39.3,39  ,38.6,38.3,37.9,37.6,37.3,36.9,36.6,36.2,35.9,35.5,35.2,34.8,34.5,34.1,33.8,33.4,33.1,32.8,32.5,32.2,31.8,31.5,31.2,30.8,
                                                                 77   ,73.2 ,69.5 ,66.1 ,63  ,60.1,57.5,55.2,53.1,51.2,49.6,48.1,46.9,45.7,44.8,43.9,43.1,42.5,41.9,41.3,40.8,40.4,40  ,39.6,39.2,38.8,38.5,38.1,37.8,37.4,37.1,36.7,36.4,36  ,35.6,35.3,34.9,34.6,34.2,33.9,33.5,33.2,32.9,32.6,32.3,31.9,31.6,31.3,30.9,30.5,
                                                                 75.8 ,72.1 ,68.5 ,65.2 ,62.1,59.3,56.8,54.6,52.5,50.7,49.2,47.8,46.6,45.5,44.5,43.7,43  ,42.3,41.7,41.2,40.7,40.3,39.9,39.5,39.1,38.7,38.3,38  ,37.6,37.2,36.9,36.5,36.1,35.8,35.4,35  ,34.7,34.3,34  ,33.6,33.3,33  ,32.7,32.4,32.1,31.8,31.4,31.1,30.7,30.1,
                                                                 74.7 ,71   ,67.5 ,64.3 ,61.3,58.6,56.2,54  ,52  ,50.3,48.8,47.5,46.3,45.2,44.3,43.5,42.8,42.2,41.6,41.1,40.6,40.1,39.7,39.3,38.9,38.5,38.2,37.8,37.4,37  ,36.7,36.3,35.9,35.5,35.2,34.8,34.5,34.1,33.8,33.4,33.1,32.8,32.5,32.2,31.9,31.6,31.2,30.8,30.4,29.8,
                                                                 73.5 ,69.9 ,66.5 ,63.4 ,60.5,57.9,55.5,53.4,51.5,49.9,48.4,47.1,46  ,45  ,44.1,43.4,42.7,42  ,41.5,41  ,40.5,40  ,39.6,39.2,38.8,38.4,38  ,37.6,37.2,36.8,36.5,36.1,35.7,35.3,35  ,34.6,34.2,33.9,33.6,33.3,32.9,32.6,32.3,32.1,31.7,31.4,31  ,30.6,30  ,29.4,
                                                                 72.4 ,68.9 ,65.5 ,62.5 ,59.7,57.2,54.9,52.9,51.1,49.5,48.1,46.8,45.8,44.8,43.9,43.2,42.5,41.9,41.3,40.8,40.3,39.9,39.4,39  ,38.6,38.2,37.8,37.4,37  ,36.6,36.2,35.9,35.5,35.1,34.7,34.4,34  ,33.7,33.4,33.1,32.8,32.5,32.2,31.9,31.6,31.2,30.8,30.3,29.7,28.9,
                                                                 71.3 ,67.8 ,64.6 ,61.6 ,58.9,56.5,54.3,52.3,50.6,49.1,47.7,46.6,45.5,44.6,43.7,43  ,42.3,41.7,41.2,40.7,40.2,39.7,39.3,38.8,38.4,38  ,37.6,37.2,36.8,36.4,36  ,35.6,35.3,34.9,34.5,34.2,33.9,33.5,33.2,32.9,32.7,32.4,32.1,31.8,31.4,31  ,30.6,30  ,29.2,28.3,
                                                                 70.3 ,66.9 ,63.7 ,60.8 ,58.2,55.8,53.7,51.9,50.2,48.7,47.4,46.3,45.3,44.3,43.5,42.8,42.2,41.6,41  ,40.5,40  ,39.5,39.1,38.7,38.2,37.8,37.4,37  ,36.6,36.2,35.8,35.4,35.1,34.7,34.4,34  ,33.7,33.4,33.1,32.8,32.5,32.3,32  ,31.6,31.3,30.8,30.3,29.6,28.7,27.6,
                                                                 69.2 ,65.9 ,62.8 ,60   ,57.5,55.2,53.2,51.4,49.8,48.3,47.1,46  ,45  ,44.1,43.3,42.6,42  ,41.4,40.8,40.3,39.8,39.4,38.9,38.5,38  ,37.6,37.2,36.8,36.4,36  ,35.6,35.2,34.9,34.5,34.2,33.9,33.6,33.3,33  ,32.7,32.4,32.1,31.8,31.5,31.1,30.6,29.9,29.1,28.1,26.9,
                                                                 68.2 ,65   ,62   ,59.3 ,56.8,54.6,52.7,50.9,49.4,48  ,46.8,45.7,44.8,43.9,43.1,42.4,41.8,41.2,40.7,40.1,39.6,39.2,38.7,38.2,37.8,37.4,37  ,36.6,36.2,35.8,35.4,35  ,34.7,34.3,34  ,33.7,33.4,33.2,32.9,32.6,32.3,32  ,31.7,31.3,30.8,30.2,29.5,28.5,27.4,26.1,
                                                                 67.2 ,64   ,61.2 ,58.5 ,56.2,54  ,52.1,50.5,49  ,47.6,46.5,45.4,44.5,43.7,42.9,42.2,41.6,41  ,40.5,39.9,39.4,38.9,38.5,38  ,37.6,37.2,36.7,36.3,35.9,35.6,35.2,34.8,34.5,34.2,33.9,33.6,33.3,33.1,32.8,32.5,32.2,31.9,31.6,31.1,30.5,29.8,28.9,27.8,26.6,25.4,
                                                                 66.2 ,63.2 ,60.4 ,57.8 ,55.5,53.5,51.6,50  ,48.6,47.3,46.2,45.2,44.3,43.4,42.7,42  ,41.4,40.8,40.2,39.7,39.2,38.7,38.3,37.8,37.4,36.9,36.5,36.1,35.7,35.4,35  ,34.7,34.4,34.1,33.8,33.5,33.2,33  ,32.7,32.4,32.1,31.8,31.4,30.8,30.2,29.3,28.3,27.1,25.9,24.6,
                                                                 65.3 ,62.3 ,59.6 ,57.1 ,54.9,52.9,51.2,49.6,48.2,47  ,45.9,44.9,44  ,43.2,42.5,41.8,41.2,40.6,40  ,39.5,39  ,38.5,38  ,37.6,37.1,36.7,36.3,35.9,35.6,35.2,34.9,34.5,34.2,33.9,33.7,33.4,33.1,32.9,32.6,32.3,32  ,31.6,31.1,30.5,29.7,28.7,27.6,26.4,25.1,24  ,
                                                                 64.4 ,61.5 ,58.8 ,56.4 ,54.3,52.4,50.7,49.2,47.8,46.6,45.6,44.6,43.7,42.9,42.2,41.6,40.9,40.3,39.8,39.3,38.8,38.3,37.8,37.4,36.9,36.5,36.1,35.7,35.4,35  ,34.7,34.4,34.1,33.8,33.6,33.3,33.1,32.8,32.5,32.2,31.9,31.4,30.8,30.1,29.2,28.1,26.9,25.6,24.4,23.4,
                                                                 63.5 ,60.7 ,58.1 ,55.8 ,53.7,51.9,50.2,48.8,47.5,46.3,45.3,44.3,43.5,42.7,42  ,41.3,40.7,40.1,39.5,39  ,38.5,38  ,37.6,37.1,36.7,36.3,35.9,35.6,35.2,34.9,34.6,34.3,34  ,33.8,33.5,33.3,33  ,32.8,32.5,32.1,31.7,31.2,30.5,29.6,28.6,27.4,26.1,24.9,23.8,23  ,
                                                                 62.6 ,59.9 ,57.4 ,55.1 ,53.2,51.4,49.8,48.4,47.1,46  ,44.9,44  ,43.2,42.4,41.7,41  ,40.4,39.8,39.3,38.8,38.3,37.8,37.4,36.9,36.5,36.1,35.8,35.4,35.1,34.8,34.5,34.2,33.9,33.7,33.5,33.2,33  ,32.7,32.4,32  ,31.5,30.8,30  ,29  ,27.9,26.6,25.4,24.3,23.3,22.6,
                                                                 61.7 ,59.1 ,56.7 ,54.5 ,52.6,50.9,49.3,47.9,46.7,45.6,44.6,43.7,42.9,42.1,41.4,40.8,40.2,39.6,39  ,38.5,38  ,37.6,37.1,36.7,36.3,36  ,35.6,35.3,35  ,34.7,34.4,34.1,33.9,33.6,33.4,33.2,32.9,32.6,32.2,31.8,31.2,30.4,29.5,28.4,27.1,25.9,24.7,23.7,23  ,22.4,
                                                                 60.9 ,58.3 ,56   ,53.9 ,52  ,50.4,48.9,47.5,46.3,45.3,44.3,43.4,42.6,41.8,41.2,40.5,39.9,39.3,38.8,38.3,37.8,37.4,36.9,36.5,36.2,35.8,35.5,35.2,34.9,34.6,34.3,34.1,33.8,33.6,33.4,33.1,32.8,32.5,32.1,31.5,30.8,29.9,28.8,27.6,26.4,25.2,24.2,23.3,22.7,22.2,
                                                                 60.1 ,57.6 ,55.3 ,53.3 ,51.5,49.9,48.4,47.1,46  ,44.9,44  ,43.1,42.3,41.6,40.9,40.2,39.6,39.1,38.6,38.1,37.6,37.2,36.7,36.4,36  ,35.7,35.3,35.1,34.8,34.5,34.3,34  ,33.8,33.6,33.3,33.1,32.7,32.3,31.8,31.2,30.3,29.3,28.1,26.9,25.7,24.6,23.7,23  ,22.5,22.1,
                                                                 59.3 ,56.9 ,54.7 ,52.7 ,51  ,49.4,48  ,46.7,45.6,44.6,43.6,42.8,42  ,41.3,40.6,40  ,39.4,38.8,38.3,37.8,37.4,37  ,36.6,36.2,35.9,35.6,35.3,35  ,34.7,34.5,34.2,34  ,33.8,33.5,33.3,33  ,32.6,32.1,31.5,30.7,29.8,28.6,27.4,26.2,25  ,24.1,23.3,22.8,22.3,22  ,
                                                                 58.5 ,56.2 ,54.1 ,52.2 ,50.4,48.9,47.6,46.3,45.2,44.2,43.3,42.4,41.7,41  ,40.3,39.7,39.1,38.6,38.1,37.6,37.2,36.8,36.4,36.1,35.8,35.5,35.2,34.9,34.7,34.4,34.2,34  ,33.8,33.5,33.2,32.9,32.4,31.9,31.1,30.2,29.1,27.9,26.7,25.5,24.5,23.7,23.1,22.6,22.3,22  ,
                                                                 57.8 ,55.5 ,53.4 ,51.6 ,49.9,48.4,47.1,45.9,44.8,43.8,42.9,42.1,41.4,40.7,40  ,39.4,38.9,38.4,37.9,37.4,37  ,36.7,36.3,36  ,35.7,35.4,35.1,34.9,34.6,34.4,34.2,34  ,33.7,33.5,33.1,32.7,32.2,31.5,30.7,29.6,28.4,27.2,26  ,24.9,24.1,23.4,22.9,22.5,22.3,22.1,
                                                                 57   ,54.8 ,52.8 ,51   ,49.4,48  ,46.7,45.5,44.4,43.5,42.6,41.8,41.1,40.4,39.8,39.2,38.7,38.2,37.7,37.3,36.9,36.5,36.2,35.9,35.6,35.3,35.1,34.8,34.6,34.4,34.2,34  ,33.7,33.4,33  ,32.5,31.9,31.1,30.1,29  ,27.7,26.5,25.4,24.5,23.7,23.2,22.8,22.5,22.3,22.1,
                                                                 56.3 ,54.1 ,52.2 ,50.4 ,48.9,47.5,46.2,45.1,44.1,43.1,42.3,41.5,40.8,40.1,39.5,39  ,38.4,38  ,37.5,37.1,36.8,36.4,36.1,35.8,35.5,35.3,35.1,34.8,34.6,34.4,34.2,33.9,33.6,33.3,32.8,32.3,31.5,30.6,29.5,28.2,27  ,25.9,24.9,24.1,23.5,23.1,22.7,22.5,22.3,22.2,
                                                                 55.5 ,53.5 ,51.6 ,49.9 ,48.4,47  ,45.8,44.7,43.7,42.8,41.9,41.2,40.5,39.8,39.3,38.7,38.2,37.8,37.4,37  ,36.6,36.3,36  ,35.8,35.5,35.3,35  ,34.8,34.6,34.4,34.2,33.9,33.6,33.1,32.6,31.9,31  ,30  ,28.8,27.5,26.3,25.3,24.5,23.8,23.4,23  ,22.7,22.5,22.4,22.3,
                                                                 54.8 ,52.8 ,51   ,49.3 ,47.9,46.5,45.3,44.3,43.3,42.4,41.6,40.9,40.2,39.6,39  ,38.5,38.1,37.6,37.2,36.9,36.6,36.3,36  ,35.7,35.5,35.3,35  ,34.8,34.6,34.4,34.1,33.8,33.4,32.9,32.3,31.5,30.4,29.3,28  ,26.8,25.8,24.9,24.2,23.6,23.3,23  ,22.8,22.6,22.5,22.4,
                                                                 54.1 ,52.1 ,50.4 ,48.8 ,47.4,46.1,44.9,43.9,42.9,42.1,41.3,40.6,40  ,39.4,38.8,38.4,37.9,37.5,37.1,36.8,36.5,36.2,36  ,35.7,35.5,35.3,35.1,34.8,34.6,34.4,34.1,33.7,33.3,32.7,31.9,30.9,29.8,28.6,27.3,26.2,25.3,24.5,24  ,23.5,23.2,23  ,22.8,22.7,22.6,22.5,
                                                                 53.4 ,51.5 ,49.8 ,48.2 ,46.9,45.6,44.5,43.5,42.6,41.7,41  ,40.3,39.7,39.2,38.7,38.2,37.8,37.4,37.1,36.7,36.4,36.2,35.9,35.7,35.5,35.3,35.1,34.9,34.6,34.3,34  ,33.6,33  ,32.3,31.4,30.3,29.1,27.9,26.7,25.7,24.9,24.3,23.8,23.5,23.2,23  ,22.9,22.8,22.7,22.6,
                                                                 52.7 ,50.9 ,49.2 ,47.7 ,46.4,45.2,44.1,43.1,42.2,41.4,40.7,40.1,39.5,39  ,38.5,38.1,37.7,37.3,37  ,36.7,36.4,36.2,35.9,35.7,35.5,35.3,35.1,34.9,34.6,34.3,33.9,33.4,32.7,31.8,30.8,29.6,28.4,27.2,26.1,25.3,24.6,24.1,23.7,23.5,23.3,23.1,23  ,22.9,22.8,22.7};

        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 7:
    {
        static const std::vector<double> arrayCurvatureMax =    {140.93,135.33,129.81,124.67,119.97,115.04,110.11,105.4 ,100.98,96.877,93.126,89.752,86.771,84.179,81.94 ,79.981,78.205,76.517,74.848,73.166,71.467,69.765,68.079,66.428,64.824,63.278,61.794,60.374,59.019,57.728,56.5  ,55.331,54.22 ,53.164,52.161,51.208,50.302,49.441,48.623,47.846,47.106,46.402,45.73 ,45.088,44.469,43.87 ,43.283,42.698,42.102,41.48 ,
                                                                 137.63,132.07,126.56,121.19,116.25,111.79,107.15,102.54,98.199,94.179,90.515,87.236,84.354,81.86 ,79.712,77.831,76.118,74.482,72.858,71.22 ,69.569,67.919,66.291,64.699,63.158,61.675,60.254,58.898,57.605,56.375,55.206,54.095,53.039,52.038,51.086,50.183,49.326,48.511,47.738,47.003,46.303,45.637,45    ,44.389,43.799,43.224,42.654,42.078,41.482,40.85 ,
                                                                 134.39,128.9 ,123.45,118.09,112.92,108.22,104.05,99.742,95.511,91.579,88.006,84.822,82.039,79.643,77.584,75.778,74.125,72.536,70.956,69.362,67.757,66.16 ,64.587,63.054,61.574,60.152,58.793,57.497,56.264,55.092,53.98 ,52.924,51.922,50.972,50.071,49.215,48.404,47.634,46.902,46.207,45.545,44.914,44.309,43.727,43.162,42.605,42.047,41.475,40.872,40.223,
                                                                 131.23,125.81,120.42,115.13,109.97,105.05,100.64,96.818,92.886,89.074,85.598,82.511,79.828,77.527,75.554,73.819,72.222,70.679,69.14 ,67.588,66.031,64.483,62.965,61.49 ,60.068,58.706,57.406,56.168,54.993,53.878,52.82 ,51.817,50.866,49.966,49.112,48.302,47.534,46.805,46.114,45.456,44.829,44.231,43.656,43.099,42.555,42.012,41.46 ,40.884,40.267,39.596,
                                                                 128.14,122.79,117.48,112.26,107.17,102.27,97.647,93.581,90.14 ,86.636,83.288,80.303,77.719,75.512,73.621,71.953,70.409,68.909,67.409,65.899,64.386,62.889,61.423,60.004,58.639,57.334,56.091,54.91 ,53.79 ,52.728,51.723,50.77 ,49.869,49.015,48.206,47.44 ,46.713,46.024,45.37 ,44.747,44.154,43.585,43.036,42.502,41.973,41.439,40.887,40.3  ,39.663,38.969,
                                                                 125.13,119.86,114.62,109.48,104.46,99.635,95.038,90.763,87.086,84.077,81.045,78.194,75.711,73.596,71.784,70.178,68.682,67.223,65.761,64.29 ,62.822,61.373,59.959,58.593,57.284,56.034,54.846,53.719,52.652,51.642,50.686,49.782,48.927,48.118,47.352,46.627,45.939,45.287,44.668,44.079,43.515,42.973,42.448,41.931,41.413,40.882,40.323,39.719,39.058,38.341,
                                                                 122.2 ,117.01,111.86,106.79,101.86,97.105,92.581,88.335,84.455,81.213,78.68 ,76.15 ,73.799,71.777,70.039,68.492,67.041,65.619,64.193,62.761,61.336,59.933,58.57 ,57.256,55.999,54.803,53.668,52.592,51.576,50.615,49.707,48.849,48.038,47.271,46.546,45.86 ,45.21 ,44.593,44.006,43.447,42.911,42.393,41.887,41.383,40.871,40.336,39.763,39.136,38.451,37.716,
                                                                 119.36,114.26,109.19,104.2 ,99.346,94.673,90.229,86.06 ,82.214,78.776,76.011,73.981,71.95 ,70.048,68.385,66.892,65.482,64.095,62.703,61.309,59.925,58.568,57.253,55.989,54.783,53.638,52.553,51.528,50.559,49.645,48.782,47.968,47.199,46.473,45.786,45.137,44.522,43.938,43.382,42.85 ,42.339,41.841,41.35 ,40.854,40.342,39.796,39.203,38.551,37.843,37.101,
                                                                 116.61,111.59,106.61,101.7 ,96.933,92.341,87.978,83.892,80.126,76.727,73.772,71.508,69.969,68.376,66.815,65.376,64.004,62.648,61.29 ,59.931,58.587,57.274,56.005,54.79 ,53.633,52.537,51.5  ,50.522,49.599,48.73 ,47.91 ,47.138,46.409,45.721,45.071,44.456,43.873,43.32 ,42.792,42.285,41.795,41.314,40.833,40.34 ,39.821,39.258,38.64 ,37.962,37.238,36.502,
                                                                 113.96,109.03,104.13,99.307,94.617,90.108,85.827,81.824,78.144,74.829,71.918,69.465,67.681,66.564,65.292,63.935,62.602,61.277,59.949,58.625,57.319,56.049,54.825,53.656,52.546,51.496,50.505,49.572,48.694,47.867,47.088,46.355,45.664,45.012,44.396,43.814,43.262,42.736,42.234,41.75 ,41.278,40.809,40.333,39.836,39.303,38.717,38.072,37.372,36.644,35.928,
                                                                 111.4 ,106.55,101.74,97.006,92.4  ,87.973,83.775,79.856,76.262,73.037,70.211,67.8  ,65.821,64.439,63.621,62.535,61.27 ,59.977,58.679,57.388,56.12 ,54.89 ,53.709,52.584,51.519,50.513,49.567,48.676,47.839,47.053,46.313,45.617,44.962,44.344,43.761,43.209,42.685,42.185,41.706,41.241,40.783,40.321,39.845,39.337,38.783,38.171,37.5  ,36.787,36.069,35.387,
                                                                 108.93,104.18,99.452,94.803,90.281,85.936,81.821,77.985,74.479,71.343,68.607,66.278,64.331,62.738,61.631,60.976,59.972,58.743,57.477,56.218,54.985,53.794,52.655,51.572,50.549,49.586,48.681,47.832,47.034,46.286,45.583,44.922,44.301,43.715,43.162,42.638,42.14 ,41.664,41.205,40.755,40.306,39.847,39.363,38.839,38.26 ,37.619,36.928,36.214,35.52 ,34.886,
                                                                 106.56,101.9 ,97.26 ,92.698,88.259,83.997,79.964,76.212,72.791,69.744,67.098,64.854,62.977,61.397,60.058,59.094,58.507,57.536,56.335,55.111,53.914,52.76 ,51.659,50.617,49.635,48.712,47.847,47.036,46.276,45.564,44.896,44.268,43.678,43.123,42.598,42.1  ,41.626,41.171,40.728,40.289,39.845,39.382,38.885,38.337,37.729,37.064,36.362,35.662,35.006,34.426,
                                                                 104.29,99.714,95.165,90.688,86.333,82.153,78.202,74.533,71.198,68.239,65.68 ,63.519,61.713,60.183,58.836,57.624,56.715,56.157,55.217,54.061,52.901,51.783,50.72 ,49.716,48.773,47.888,47.061,46.286,45.562,44.884,44.248,43.652,43.092,42.565,42.066,41.592,41.139,40.701,40.271,39.84 ,39.395,38.922,38.404,37.829,37.192,36.509,35.81 ,35.138,34.531,34.011,
                                                                 102.11,97.627,93.167,88.775,84.502,80.404,76.534,72.948,69.698,66.824,64.351,62.27 ,60.531,59.051,57.732,56.496,55.333,54.448,53.921,53.03 ,51.941,50.862,49.835,48.868,47.961,47.113,46.32 ,45.581,44.89 ,44.244,43.639,43.073,42.541,42.039,41.564,41.112,40.678,40.254,39.833,39.403,38.951,38.461,37.917,37.312,36.652,35.961,35.28 ,34.65 ,34.099,33.638,
                                                                 100.03,95.635,91.262,86.955,82.765,78.747,74.958,71.454,68.287,65.499,63.109,61.103,59.427,57.992,56.702,55.48 ,54.286,53.147,52.292,51.815,50.996,49.989,49.001,48.068,47.196,46.382,45.624,44.917,44.257,43.642,43.066,42.527,42.021,41.544,41.091,40.658,40.238,39.825,39.408,38.974,38.509,37.996,37.423,36.789,36.113,35.429,34.78 ,34.201,33.71 ,33.305,
                                                                 98.038,93.737,89.451,85.228,81.119,77.181,73.472,70.049,66.965,64.259,61.95 ,60.016,58.396,57.002,55.737,54.53 ,53.344,52.174,51.07 ,50.265,49.86 ,49.126,48.21 ,47.316,46.477,45.695,44.968,44.292,43.662,43.076,42.527,42.014,41.532,41.076,40.643,40.226,39.818,39.411,38.993,38.55 ,38.064,37.523,36.919,36.263,35.583,34.921,34.317,33.795,33.361,33.008,
                                                                 96.144,91.932,87.732,83.591,79.563,75.704,72.074,68.731,65.727,63.103,60.871,59.004,57.437,56.077,54.833,53.639,52.464,51.303,50.171,49.122,48.388,48.066,47.425,46.602,45.799,45.049,44.352,43.705,43.103,42.543,42.02 ,41.531,41.07 ,40.635,40.218,39.814,39.415,39.009,38.584,38.124,37.613,37.041,36.409,35.74 ,35.071,34.445,33.894,33.43 ,33.051,32.743,
                                                                 94.342,90.217,86.102,82.043,78.094,74.313,70.761,67.496,64.573,62.028,59.869,58.066,56.545,55.215,53.988,52.804,51.638,50.489,49.368,48.295,47.321,46.669,46.437,45.89 ,45.157,44.441,43.773,43.153,42.578,42.042,41.543,41.076,40.635,40.217,39.814,39.42 ,39.024,38.615,38.177,37.694,37.152,36.549,35.897,35.228,34.585,34.006,33.512,33.105,32.775,32.507,
                                                                 92.629,88.592,84.56 ,80.581,76.71 ,73.006,69.531,66.344,63.499,61.031,58.942,57.197,55.718,54.412,53.197,52.021,50.863,49.725,48.619,47.561,46.564,45.675,45.11 ,44.969,44.512,43.865,43.228,42.634,42.084,41.572,41.094,40.646,40.224,39.82 ,39.428,39.039,38.643,38.223,37.766,37.255,36.681,36.051,35.389,34.736,34.132,33.608,33.171,32.816,32.529,32.296,
                                                                 91.005,87.053,83.103,79.203,75.409,71.781,68.381,65.27 ,62.502,60.108,58.087,56.394,54.951,53.666,52.46 ,51.288,50.136,49.008,47.917,46.877,45.896,44.986,44.186,43.707,43.653,43.281,42.711,42.147,41.619,41.129,40.671,40.241,39.834,39.442,39.057,38.67 ,38.266,37.831,37.348,36.804,36.2  ,35.552,34.895,34.271,33.717,33.25 ,32.868,32.56 ,32.31 ,32.106,
                                                                 89.467,85.598,81.728,77.907,74.188,70.634,67.308,64.273,61.58 ,59.258,57.299,55.655,54.243,52.972,51.772,50.603,49.455,48.336,47.259,46.236,45.275,44.382,43.56 ,42.847,42.451,42.477,42.183,41.683,41.182,40.712,40.272,39.858,39.463,39.08 ,38.699,38.307,37.891,37.433,36.919,36.343,35.713,35.059,34.422,33.841,33.342,32.931,32.599,32.332,32.114,31.934,
                                                                 88.012,84.225,80.435,76.689,73.044,69.564,66.311,63.35 ,60.73 ,58.476,56.576,54.975,53.589,52.329,51.131,49.963,48.818,47.707,46.642,45.635,44.693,43.82 ,43.015,42.28 ,41.65 ,41.33 ,41.429,41.204,40.765,40.318,39.895,39.495,39.109,38.731,38.348,37.947,37.512,37.026,36.479,35.872,35.228,34.582,33.979,33.45 ,33.007,32.649,32.361,32.129,31.938,31.778,
                                                                 86.639,82.932,79.218,75.547,71.976,68.567,65.387,62.498,59.949,57.761,55.915,54.352,52.986,51.733,50.535,49.365,48.222,47.118,46.064,45.072,44.148,43.294,42.509,41.79 ,41.137,40.583,40.331,40.494,40.329,39.942,39.538,39.148,38.77 ,38.392,38.003,37.586,37.125,36.607,36.026,35.397,34.75 ,34.129,33.571,33.097,32.71 ,32.399,32.15 ,31.947,31.779,31.636,
                                                                 85.343,81.715,78.077,74.479,70.979,67.641,64.532,61.714,59.234,57.109,55.313,53.783,52.432,51.182,49.981,48.808,47.666,46.567,45.524,44.546,43.638,42.802,42.036,41.336,40.697,40.118,39.632,39.439,39.657,39.541,39.193,38.817,38.441,38.06 ,37.658,37.219,36.728,36.175,35.565,34.924,34.291,33.707,33.202,32.784,32.447,32.178,31.962,31.784,31.635,31.506,
                                                                 84.124,80.572,77.008,73.482,70.052,66.784,63.745,60.997,58.583,56.517,54.766,53.264,51.924,50.673,49.467,48.29 ,47.148,46.053,45.019,44.053,43.161,42.342,41.593,40.911,40.289,39.723,39.21 ,38.782,38.639,38.899,38.819,38.491,38.12 ,37.73 ,37.309,36.842,36.317,35.731,35.101,34.462,33.857,33.321,32.871,32.506,32.215,31.982,31.793,31.637,31.503,31.387,
                                                                 82.978,79.501,76.009,72.553,69.192,65.992,63.022,60.342,57.993,55.982,54.272,52.792,51.459,50.204,48.99 ,47.808,46.664,45.574,44.547,43.593,42.715,41.911,41.179,40.513,39.908,39.357,38.854,38.397,38.016,37.909,38.196,38.13 ,37.797,37.397,36.952,36.453,35.892,35.279,34.641,34.02 ,33.456,32.973,32.578,32.262,32.01 ,31.808,31.642,31.503,31.383,31.277,
                                                                 81.902,78.498,75.077,71.69 ,68.396,65.264,62.362,59.748,57.461,55.501,53.827,52.365,51.034,49.773,48.55 ,47.36 ,46.215,45.127,44.108,43.164,42.298,41.509,40.792,40.141,39.551,39.014,38.524,38.074,37.661,37.312,37.221,37.511,37.428,37.053,36.583,36.049,35.456,34.825,34.195,33.606,33.091,32.664,32.32 ,32.047,31.829,31.652,31.506,31.381,31.273,31.176,
                                                                 80.894,77.561,74.209,70.889,67.663,64.597,61.761,59.212,56.984,55.072,53.429,51.98 ,50.647,49.376,48.142,46.945,45.797,44.711,43.698,42.764,41.91 ,41.133,40.429,39.793,39.216,38.692,38.213,37.773,37.361,36.975,36.636,36.534,36.791,36.652,36.194,35.63 ,35.013,34.379,33.77 ,33.225,32.764,32.39 ,32.093,31.857,31.668,31.512,31.382,31.27 ,31.171,31.083,
                                                                 79.95 ,76.687,73.403,70.149,66.988,63.988,61.217,58.731,56.559,54.691,53.074,51.633,50.296,49.013,47.767,46.561,45.409,44.325,43.317,42.391,41.547,40.782,40.091,39.467,38.902,38.389,37.921,37.487,37.08 ,36.687,36.301,35.94 ,35.788,35.971,35.743,35.195,34.57 ,33.948,33.375,32.881,32.474,32.15 ,31.893,31.689,31.523,31.386,31.269,31.168,31.077,30.996,
                                                                 79.069,75.874,72.656,69.467,66.371,63.435,60.728,58.303,56.185,54.355,52.76 ,51.323,49.977,48.681,47.421,46.206,45.05 ,43.966,42.963,42.044,41.21 ,40.455,39.775,39.161,38.607,38.104,37.643,37.215,36.808,36.41 ,36.006,35.584,35.16 ,34.92 ,35.009,34.707,34.13 ,33.54 ,33.014,32.574,32.22 ,31.939,31.717,31.54 ,31.394,31.272,31.166,31.074,30.991,30.915,
                                                                 78.248,75.119,71.965,68.84 ,65.807,62.936,60.292,57.926,55.857,54.063,52.485,51.047,49.69 ,48.378,47.104,45.879,44.718,43.633,42.634,41.722,40.895,40.15 ,39.479,38.876,38.331,37.835,37.379,36.953,36.543,36.135,35.713,35.258,34.759,34.242,33.914,33.951,33.66 ,33.157,32.689,32.303,31.996,31.754,31.562,31.406,31.277,31.167,31.072,30.987,30.91 ,30.84 ,
                                                                 77.483,74.419,71.328,68.266,65.297,62.488,59.906,57.597,55.574,53.811,52.245,50.803,49.432,48.102,46.813,45.578,44.411,43.326,42.329,41.422,40.603,39.865,39.203,38.608,38.07 ,37.58 ,37.127,36.699,36.282,35.859,35.413,34.925,34.383,33.792,33.206,32.858,32.944,32.761,32.395,32.065,31.801,31.592,31.424,31.287,31.171,31.072,30.985,30.907,30.836,30.771,
                                                                 76.773,73.772,70.744,67.743,64.836,62.089,59.567,57.313,55.334,53.596,52.038,50.589,49.2  ,47.852,46.548,45.301,44.128,43.041,42.046,41.144,40.331,39.6  ,38.945,38.356,37.823,37.336,36.883,36.45 ,36.02 ,35.577,35.1  ,34.575,33.995,33.372,32.741,32.177,31.913,32.126,32.086,31.851,31.629,31.448,31.301,31.179,31.075,30.984,30.904,30.832,30.766,30.706,
                                                                 76.115,73.175,70.208,67.269,64.422,61.737,59.274,57.072,55.132,53.417,51.862,50.402,48.994,47.627,46.306,45.047,43.868,42.779,41.785,40.886,40.078,39.353,38.703,38.119,37.59 ,37.103,36.646,36.202,35.755,35.285,34.773,34.209,33.596,32.957,32.33 ,31.762,31.312,31.184,31.533,31.614,31.471,31.319,31.19 ,31.081,30.987,30.904,30.83 ,30.763,30.702,30.645,
                                                                 75.507,72.627,69.72 ,66.841,64.055,61.43 ,59.024,56.872,54.968,53.27 ,51.714,50.24 ,48.812,47.424,46.086,44.816,43.63 ,42.538,41.544,40.647,39.843,39.122,38.477,37.895,37.366,36.877,36.411,35.953,35.482,34.979,34.428,33.827,33.19 ,32.549,31.947,31.418,30.982,30.672,30.668,31.12 ,31.282,31.198,31.089,30.992,30.906,30.83 ,30.761,30.699,30.641,30.589,
                                                                 74.945,72.124,69.277,66.457,63.731,61.166,58.816,56.71 ,54.839,53.154,51.593,50.102,48.652,47.242,45.887,44.605,43.411,42.316,41.321,40.426,39.624,38.906,38.263,37.683,37.152,36.655,36.177,35.698,35.197,34.656,34.067,33.434,32.785,32.158,31.593,31.114,30.726,30.424,30.226,30.309,30.826,31.039,30.992,30.909,30.831,30.761,30.697,30.639,30.585,30.536,
                                                                 74.428,71.666,68.876,66.116,63.449,60.942,58.646,56.585,54.741,53.066,51.497,49.986,48.513,47.081,45.708,44.413,43.212,42.112,41.116,40.222,39.421,38.704,38.061,37.479,36.943,36.436,35.939,35.433,34.897,34.316,33.69 ,33.036,32.39 ,31.792,31.273,30.846,30.506,30.241,30.037,29.914,30.049,30.607,30.851,30.827,30.761,30.697,30.638,30.583,30.533,30.486,
                                                                 73.954,71.249,68.517,65.815,63.207,60.758,58.514,56.494,54.674,53.004,51.423,49.89 ,48.393,46.939,45.548,44.24 ,43.03 ,41.926,40.928,40.032,39.231,38.514,37.869,37.282,36.737,36.214,35.694,35.155,34.579,33.959,33.303,32.642,32.014,31.456,30.988,30.612,30.316,30.085,29.904,29.763,29.685,29.852,30.434,30.7  ,30.69 ,30.637,30.582,30.531,30.484,30.44 ,
                                                                 73.521,70.871,68.196,65.552,63.003,60.61 ,58.417,56.435,54.635,52.965,51.37 ,49.814,48.291,46.814,45.405,44.084,42.865,41.755,40.754,39.857,39.054,38.335,37.685,37.09 ,36.53 ,35.987,35.436,34.86 ,34.243,33.587,32.914,32.26 ,31.664,31.154,30.738,30.409,30.152,29.951,29.79 ,29.661,29.558,29.508,29.693,30.291,30.572,30.574,30.53 ,30.483,30.439,30.397,
                                                                 73.126,70.532,67.913,65.326,62.835,60.498,58.353,56.405,54.622,52.949,51.336,49.755,48.205,46.706,45.278,43.944,42.716,41.6  ,40.595,39.694,38.888,38.164,37.507,36.899,36.321,35.75 ,35.164,34.546,33.889,33.207,32.53 ,31.899,31.346,30.888,30.521,30.234,30.01 ,29.833,29.69 ,29.574,29.476,29.397,29.364,29.562,30.17 ,30.463,30.474,30.437,30.396,30.357,
                                                                 72.768,70.228,67.665,65.136,62.702,60.419,58.32 ,56.403,54.632,52.952,51.319,49.712,48.136,46.613,45.168,43.82 ,42.582,41.459,40.448,39.542,38.731,37.999,37.331,36.706,36.104,35.5  ,34.874,34.214,33.523,32.827,32.162,31.567,31.063,30.655,30.333,30.082,29.885,29.729,29.601,29.495,29.406,29.328,29.265,29.244,29.449,30.066,30.369,30.387,30.355,30.32 ,
                                                                 72.445,69.959,67.452,64.979,62.602,60.372,58.317,56.428,54.664,52.975,51.319,49.684,48.081,46.535,45.071,43.709,42.461,41.33 ,40.312,39.4  ,38.581,37.839,37.155,36.509,35.875,35.233,34.564,33.865,33.151,32.456,31.818,31.266,30.813,30.453,30.172,29.951,29.777,29.637,29.522,29.425,29.342,29.269,29.205,29.153,29.141,29.352,29.974,30.286,30.31 ,30.284,
                                                                 72.155,69.722,67.27 ,64.854,62.534,60.355,58.341,56.475,54.716,53.014,51.334,49.67 ,48.04 ,46.471,44.988,43.612,42.353,41.212,40.187,39.266,38.436,37.68 ,36.976,36.302,35.632,34.947,34.236,33.505,32.781,32.102,31.502,31.001,30.597,30.28 ,30.032,29.837,29.682,29.555,29.45 ,29.36 ,29.283,29.215,29.155,29.101,29.057,29.051,29.267,29.894,30.214,30.243,
                                                                 71.896,69.517,67.12 ,64.76 ,62.495,60.367,58.39 ,56.545,54.787,53.068,51.362,49.669,48.012,46.419,44.918,43.527,42.256,41.105,40.069,39.138,38.294,37.519,36.79 ,36.082,35.371,34.642,33.893,33.141,32.422,31.774,31.22 ,30.769,30.411,30.131,29.912,29.738,29.598,29.482,29.385,29.302,29.229,29.165,29.108,29.057,29.011,28.973,28.973,29.192,29.823,30.15 ,
                                                                 71.667,69.341,66.998,64.695,62.485,60.405,58.463,56.635,54.874,53.136,51.402,49.68 ,47.996,46.38 ,44.86 ,43.453,42.169,41.006,39.959,39.013,38.152,37.353,36.593,35.845,35.09 ,34.318,33.54 ,32.782,32.084,31.476,30.972,30.568,30.251,30.003,29.808,29.651,29.523,29.416,29.326,29.248,29.179,29.119,29.065,29.016,28.972,28.932,28.9  ,28.905,29.125,29.76 ,
                                                                 71.467,69.193,66.906,64.658,62.502,60.468,58.559,56.743,54.976,53.217,51.455,49.703,47.991,46.351,44.812,43.389,42.091,40.915,39.853,38.89 ,38.007,37.179,36.381,35.59 ,34.789,33.981,33.185,32.438,31.773,31.212,30.757,30.397,30.115,29.893,29.717,29.574,29.456,29.357,29.272,29.198,29.134,29.076,29.025,28.978,28.936,28.898,28.864,28.836,28.845,29.067,
                                                                 71.294,69.072,66.84 ,64.648,62.544,60.555,58.674,56.868,55.091,53.309,51.517,49.735,47.996,46.333,44.774,43.335,42.021,40.829,39.75 ,38.766,37.855,36.992,36.152,35.314,34.471,33.636,32.839,32.115,31.493,30.981,30.572,30.251,29.998,29.798,29.637,29.505,29.396,29.303,29.223,29.153,29.091,29.037,28.988,28.943,28.903,28.867,28.833,28.803,28.78 ,28.792,
                                                                 71.147,68.977,66.799,64.663,62.611,60.663,58.809,57.008,55.218,53.411,51.589,49.777,48.01 ,46.324,44.745,43.287,41.956,40.747,39.647,38.637,37.693,36.789,35.902,35.018,34.141,33.292,32.509,31.821,31.247,30.782,30.415,30.126,29.898,29.716,29.567,29.445,29.341,29.253,29.177,29.111,29.052,29    ,28.953,28.911,28.872,28.837,28.806,28.776,28.75 ,28.73 ,
                                                                 71.025,68.907,66.783,64.702,62.701,60.792,58.96 ,57.162,55.356,53.523,51.67 ,49.827,48.033,46.323,44.723,43.247,41.896,40.666,39.541,38.499,37.517,36.567,35.632,34.707,33.805,32.959,32.203,31.56 ,31.033,30.613,30.282,30.02 ,29.812,29.644,29.505,29.39 ,29.292,29.209,29.136,29.072,29.016,28.966,28.921,28.88 ,28.843,28.81 ,28.78 ,28.752,28.726,28.703};

        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 11:
    {
        static const std::vector<double> arrayCurvatureMax =    {151.19,150.15,148.46,146.19,143.47,140.47,137.34,134.2 ,131.15,128.24,125.5 ,122.94,120.58,118.39,116.36,114.48,112.71,111.03,109.42,107.84,106.28,104.71,103.1 ,101.44,99.72 ,97.933,96.086,94.194,92.282,90.376,88.507,86.7  ,84.975,83.342,81.806,80.365,79.008,77.727,76.509,75.342,74.216,73.122,72.051,70.997,69.956,68.925,67.902,66.886,65.876,64.873,
                                                                 151.61,150.24,148.24,145.68,142.74,139.55,136.28,133.03,129.89,126.91,124.11,121.51,119.1 ,116.87,114.8 ,112.87,111.05,109.31,107.64,105.99,104.35,102.69,100.99,99.241,97.43 ,95.562,93.648,91.71 ,89.775,87.87 ,86.024,84.258,82.584,81.01 ,79.532,78.145,76.837,75.598,74.415,73.276,72.172,71.094,70.036,68.993,67.96 ,66.936,65.919,64.909,63.906,62.91 ,
                                                                 151.86,150.16,147.85,145.03,141.87,138.51,135.1 ,131.75,128.52,125.48,122.63,119.98,117.53,115.25,113.14,111.16,109.28,107.49,105.75,104.03,102.3 ,100.55,98.765,96.924,95.03 ,93.093,91.129,89.164,87.225,85.341,83.534,81.82 ,80.206,78.692,77.273,75.94 ,74.679,73.48 ,72.329,71.217,70.134,69.074,68.03 ,66.998,65.976,64.961,63.953,62.951,61.957,60.972,
                                                                 151.95,149.92,147.31,144.24,140.86,137.34,133.81,130.35,127.06,123.95,121.04,118.35,115.86,113.54,111.38,109.35,107.42,105.56,103.75,101.95,100.14,98.305,96.428,94.505,92.54 ,90.549,88.554,86.582,84.66 ,82.812,81.057,79.404,77.854,76.404,75.044,73.762,72.547,71.385,70.266,69.179,68.118,67.075,66.047,65.028,64.018,63.015,62.018,61.029,60.047,59.075,
                                                                 151.88,149.53,146.63,143.32,139.74,136.06,132.41,128.86,125.49,122.32,119.37,116.63,114.1 ,111.74,109.53,107.45,105.46,103.53,101.65,99.77 ,97.879,95.958,93.999,92.004,89.984,87.958,85.951,83.992,82.104,80.308,78.616,77.03 ,75.547,74.16 ,72.857,71.626,70.453,69.328,68.239,67.178,66.138,65.114,64.102,63.098,62.102,61.113,60.131,59.156,58.19 ,57.233,
                                                                 151.66,148.99,145.81,142.26,138.5 ,134.68,130.91,127.28,123.84,120.62,117.62,114.84,112.26,109.85,107.6 ,105.46,103.41,101.42,99.461,97.503,95.53 ,93.529,91.499,89.446,87.387,85.346,83.349,81.421,79.584,77.851,76.228,74.713,73.298,71.974,70.726,69.543,68.411,67.321,66.262,65.226,64.208,63.204,62.209,61.222,60.243,59.27 ,58.304,57.346,56.398,55.46 ,
                                                                 151.29,148.32,144.87,141.1 ,137.16,133.2 ,129.33,125.62,122.12,118.85,115.8 ,112.98,110.36,107.91,105.6 ,103.41,101.3 ,99.238,97.2  ,95.164,93.113,91.041,88.951,86.857,84.779,82.743,80.775,78.896,77.123,75.461,73.912,72.469,71.122,69.858,68.663,67.526,66.434,65.377,64.346,63.336,62.341,61.357,60.381,59.413,58.452,57.497,56.55 ,55.611,54.682,53.765,
                                                                 150.78,147.52,143.81,139.83,135.73,131.65,127.68,123.9 ,120.34,117.01,113.93,111.07,108.4 ,105.9 ,103.55,101.3 ,99.123,96.992,94.88 ,92.768,90.647,88.515,86.381,84.265,82.188,80.178,78.257,76.441,74.74 ,73.155,71.682,70.311,69.029,67.823,66.679,65.586,64.531,63.507,62.505,61.52 ,60.547,59.584,58.629,57.68 ,56.738,55.803,54.876,53.958,53.051,52.156,
                                                                 150.16,146.61,142.66,138.48,134.23,130.03,125.97,122.12,118.51,115.14,112.02,109.11,106.4 ,103.86,101.45,99.143,96.902,94.7  ,92.517,90.336,88.154,85.976,83.817,81.698,79.644,77.678,75.818,74.075,72.453,70.948,69.551,68.25 ,67.03 ,65.879,64.784,63.731,62.713,61.719,60.745,59.785,58.836,57.895,56.961,56.033,55.111,54.197,53.291,52.395,51.51 ,50.639,
                                                                 149.44,145.6 ,141.43,137.06,132.67,128.36,124.22,120.31,116.65,113.24,110.08,107.13,104.38,101.79,99.327,96.957,94.649,92.379,90.128,87.887,85.658,83.451,81.285,79.184,77.171,75.264,73.478,71.815,70.276,68.85 ,67.528,66.293,65.134,64.034,62.983,61.97 ,60.986,60.022,59.075,58.14 ,57.214,56.295,55.382,54.475,53.575,52.683,51.799,50.925,50.064,49.216,
                                                                 148.62,144.53,140.13,135.6 ,131.07,126.66,122.45,118.48,114.78,111.33,108.12,105.14,102.34,99.704,97.184,94.753,92.38 ,90.045,87.734,85.444,83.182,80.965,78.812,76.748,74.792,72.958,71.253,69.675,68.218,66.871,65.619,64.448,63.344,62.293,61.284,60.308,59.355,58.421,57.5  ,56.589,55.686,54.789,53.897,53.012,52.133,51.263,50.401,49.551,48.713,47.889,
                                                                 147.74,143.39,138.79,134.1 ,129.45,124.95,120.68,116.65,112.9 ,109.42,106.17,103.15,100.31,97.614,95.037,92.544,90.108,87.714,85.353,83.028,80.751,78.542,76.423,74.412,72.527,70.774,69.155,67.663,66.287,65.015,63.83 ,62.718,61.665,60.659,59.689,58.747,57.825,56.918,56.022,55.134,54.253,53.378,52.508,51.644,50.787,49.938,49.099,48.271,47.456,46.657,
                                                                 146.81,142.22,137.43,132.59,127.83,123.25,118.91,114.84,111.04,107.52,104.23,101.17,98.277,95.531,92.896,90.342,87.848,85.403,83.004,80.66 ,78.387,76.206,74.136,72.194,70.389,68.723,67.191,65.783,64.486,63.284,62.161,61.104,60.098,59.133,58.199,57.288,56.395,55.513,54.641,53.776,52.917,52.062,51.213,50.371,49.535,48.707,47.89 ,47.085,46.293,45.517,
                                                                 145.84,141.03,136.06,131.08,126.21,121.56,117.16,113.04,109.21,105.64,102.32,99.206,96.266,93.465,90.772,88.16 ,85.615,83.129,80.707,78.362,76.111,73.975,71.97 ,70.106,68.388,66.811,65.365,64.039,62.815,61.678,60.613,59.605,58.643,57.715,56.814,55.932,55.065,54.207,53.358,52.514,51.675,50.841,50.012,49.19 ,48.374,47.568,46.772,45.989,45.219,44.466,
                                                                 144.86,139.83,134.69,129.58,124.63,119.9 ,115.44,111.28,107.41,103.8 ,100.43,97.272,94.28 ,91.423,88.673,86.009,83.421,80.908,78.479,76.151,73.94 ,71.865,69.936,68.158,66.529,65.04 ,63.678,62.428,61.272,60.195,59.182,58.219,57.296,56.402,55.531,54.675,53.832,52.997,52.168,51.344,50.525,49.71 ,48.901,48.098,47.302,46.516,45.741,44.979,44.231,43.5  ,
                                                                 143.88,138.65,133.35,128.12,123.07,118.28,113.77,109.56,105.65,102   ,98.584,95.371,92.325,89.412,86.609,83.9  ,81.281,78.756,76.338,74.043,71.888,69.886,68.041,66.353,64.813,63.41 ,62.127,60.947,59.854,58.831,57.864,56.942,56.053,55.19 ,54.346,53.514,52.693,51.879,51.069,50.264,49.463,48.666,47.875,47.091,46.314,45.547,44.792,44.05 ,43.323,42.612,
                                                                 142.91,137.49,132.03,126.69,121.56,116.7 ,112.14,107.89,103.93,100.23,96.771,93.506,90.404,87.438,84.588,81.844,79.207,76.686,74.296,72.052,69.966,68.045,66.29 ,64.692,63.239,61.917,60.707,59.591,58.554,57.579,56.655,55.768,54.91 ,54.074,53.254,52.445,51.643,50.847,50.056,49.268,48.484,47.704,46.93 ,46.163,45.404,44.656,43.919,43.196,42.489,41.799,
                                                                 141.96,136.35,130.76,125.31,120.1 ,115.18,110.57,106.27,102.26,98.513,94.998,91.678,88.522,85.506,82.617,79.851,77.212,74.713,72.367,70.187,68.18 ,66.347,64.681,63.172,61.802,60.555,59.411,58.354,57.366,56.435,55.547,54.692,53.861,53.049,52.25 ,51.46 ,50.677,49.898,49.122,48.35 ,47.582,46.818,46.061,45.31 ,44.568,43.837,43.118,42.413,41.725,41.054,
                                                                 141.04,135.26,129.53,123.98,118.69,113.71,109.04,104.69,100.63,96.834,93.263,89.888,86.681,83.623,80.706,77.932,75.308,72.846,70.559,68.454,66.532,64.789,63.212,61.787,60.494,59.316,58.232,57.227,56.283,55.389,54.533,53.706,52.9  ,52.109,51.329,50.556,49.788,49.024,48.263,47.506,46.752,46.003,45.26 ,44.525,43.799,43.084,42.382,41.695,41.024,40.371,
                                                                 140.16,134.21,128.35,122.7 ,117.34,112.29,107.57,103.16,99.046,95.193,91.566,88.138,84.885,81.794,78.863,76.096,73.503,71.094,68.877,66.855,65.022,63.368,61.877,60.531,59.309,58.192,57.162,56.202,55.297,54.436,53.608,52.804,52.019,51.246,50.482,49.725,48.971,48.221,47.473,46.729,45.988,45.253,44.524,43.803,43.092,42.392,41.706,41.035,40.381,39.745,
                                                                 139.31,133.2 ,127.22,121.48,116.03,110.92,106.13,101.67,97.499,93.588,89.906,86.428,83.137,80.026,77.095,74.352,71.805,69.462,67.325,65.39 ,63.646,62.078,60.668,59.394,58.236,57.174,56.191,55.27 ,54.399,53.566,52.762,51.979,51.212,50.455,49.705,48.961,48.219,47.481,46.746,46.013,45.285,44.562,43.846,43.138,42.441,41.756,41.085,40.429,39.791,39.171,
                                                                 138.49,132.23,126.13,120.29,114.77,109.59,104.74,100.22,95.986,92.016,88.281,84.761,81.444,78.326,75.411,72.708,70.221,67.953,65.901,64.054,62.398,60.912,59.576,58.367,57.266,56.253,55.31 ,54.424,53.582,52.773,51.99 ,51.224,50.472,49.728,48.991,48.258,47.528,46.8  ,46.075,45.353,44.636,43.925,43.221,42.526,41.841,41.17 ,40.513,39.872,39.249,38.644,
                                                                 137.71,131.3 ,125.08,119.15,113.55,108.3 ,103.38,98.795,94.502,90.475,86.693,83.14 ,79.809,76.7  ,73.818,71.168,68.752,66.567,64.602,62.843,61.27 ,59.86 ,58.591,57.441,56.39 ,55.419,54.512,53.655,52.837,52.049,51.283,50.532,49.793,49.061,48.334,47.61 ,46.89 ,46.171,45.456,44.744,44.037,43.336,42.643,41.96 ,41.288,40.629,39.985,39.358,38.749,38.159,
                                                                 136.96,130.4 ,124.07,118.04,112.36,107.03,102.05,97.397,93.043,88.965,85.144,81.57 ,78.24 ,75.156,72.322,69.738,67.401,65.302,63.424,61.749,60.253,58.913,57.705,56.606,55.599,54.664,53.787,52.955,52.158,51.387,50.635,49.897,49.168,48.445,47.727,47.012,46.3  ,45.59 ,44.883,44.18 ,43.482,42.791,42.108,41.436,40.775,40.129,39.498,38.884,38.288,37.712,
                                                                 136.24,129.54,123.09,116.96,111.19,105.79,100.74,96.019,91.608,87.484,83.636,80.056,76.744,73.7  ,70.926,68.417,66.164,64.152,62.359,60.763,59.339,58.061,56.906,55.853,54.883,53.979,53.127,52.316,51.537,50.78 ,50.04 ,49.312,48.592,47.877,47.167,46.459,45.753,45.051,44.351,43.656,42.967,42.285,41.612,40.95 ,40.3  ,39.665,39.046,38.444,37.862,37.299,
                                                                 135.53,128.69,122.12,115.9 ,110.05,104.56,99.441,94.658,90.195,86.037,82.175,78.605,75.325,72.336,69.633,67.206,65.039,63.112,61.4  ,59.877,58.517,57.294,56.186,55.172,54.234,53.356,52.526,51.733,50.967,50.222,49.492,48.772,48.059,47.351,46.647,45.945,45.245,44.549,43.856,43.169,42.487,41.814,41.15 ,40.498,39.858,39.234,38.627,38.037,37.467,36.917,
                                                                 134.84,127.85,121.17,114.85,108.91,103.34,98.152,93.312,88.807,84.628,80.767,77.222,73.99 ,71.067,68.443,66.1  ,64.019,62.174,60.537,59.081,57.778,56.604,55.537,54.556,53.645,52.789,51.977,51.197,50.443,49.707,48.985,48.271,47.564,46.862,46.163,45.466,44.772,44.081,43.395,42.714,42.039,41.374,40.719,40.076,39.446,38.833,38.236,37.658,37.099,36.561,
                                                                 134.15,127.02,120.22,113.8 ,107.77,102.13,96.872,91.982,87.448,83.261,79.418,75.914,72.743,69.895,67.354,65.097,63.098,61.33 ,59.762,58.365,57.113,55.981,54.949,53.996,53.108,52.27 ,51.472,50.705,49.959,49.231,48.514,47.806,47.104,46.405,45.71 ,45.018,44.329,43.643,42.962,42.287,41.62 ,40.962,40.315,39.681,39.061,38.457,37.871,37.304,36.756,36.23 ,
                                                                 133.45,126.19,119.27,112.75,106.63,100.91,95.599,90.671,86.121,81.944,78.134,74.685,71.585,68.818,66.362,64.189,62.269,60.571,59.065,57.721,56.513,55.418,54.415,53.486,52.617,51.794,51.008,50.249,49.511,48.788,48.076,47.372,46.673,45.978,45.286,44.597,43.912,43.231,42.555,41.886,41.225,40.574,39.935,39.31 ,38.699,38.105,37.529,36.972,36.436,35.921,
                                                                 132.74,125.34,118.31,111.68,105.48,99.699,94.337,89.384,84.835,80.683,76.922,73.539,70.518,67.835,65.463,63.37 ,61.523,59.889,58.438,57.141,55.971,54.907,53.929,53.019,52.165,51.355,50.577,49.826,49.093,48.374,47.666,46.964,46.268,45.576,44.887,44.201,43.519,42.842,42.171,41.508,40.853,40.209,39.578,38.96 ,38.358,37.774,37.208,36.661,36.136,35.631,
                                                                 132.01,124.47,117.32,110.6 ,104.32,98.486,93.089,88.128,83.596,79.486,75.785,72.478,69.539,66.941,64.65 ,62.632,60.852,59.276,57.873,56.616,55.479,54.441,53.483,52.591,51.749,50.948,50.177,49.431,48.702,47.986,47.28 ,46.581,45.886,45.196,44.509,43.826,43.147,42.474,41.808,41.149,40.501,39.863,39.239,38.63 ,38.037,37.461,36.905,36.368,35.853,35.359,
                                                                 131.25,123.58,116.32,109.51,103.16,97.278,91.863,86.91 ,82.411,78.356,74.727,71.5  ,68.646,66.131,63.917,61.968,60.248,58.723,57.363,56.14 ,55.031,54.015,53.074,52.194,51.362,50.568,49.803,49.061,48.334,47.621,46.916,46.217,45.524,44.835,44.15 ,43.47 ,42.794,42.124,41.462,40.809,40.166,39.535,38.918,38.316,37.732,37.165,36.618,36.092,35.586,35.103,
                                                                 130.46,122.66,115.29,108.4 ,101.99,96.083,90.665,85.737,81.287,77.299,73.748,70.606,67.835,65.398,63.256,61.37 ,59.704,58.224,56.9  ,55.706,54.62 ,53.622,52.696,51.826,51.002,50.213,49.452,48.711,47.987,47.274,46.57 ,45.872,45.18 ,44.492,43.808,43.13 ,42.457,41.79 ,41.132,40.484,39.847,39.222,38.612,38.018,37.442,36.884,36.346,35.829,35.334,34.86 ,
                                                                 129.64,121.71,114.24,107.28,100.83,94.907,89.505,84.617,80.229,76.316,72.848,69.79 ,67.1  ,64.737,62.661,60.831,59.211,57.77 ,56.478,55.309,54.243,53.259,52.344,51.482,50.663,49.878,49.119,48.38 ,47.657,46.944,46.24 ,45.542,44.85 ,44.163,43.481,42.804,42.134,41.471,40.817,40.173,39.541,38.923,38.32 ,37.734,37.166,36.617,36.088,35.581,35.095,34.631,
                                                                 128.77,120.72,113.17,106.15,99.678,93.759,88.389,83.557,79.24 ,75.408,72.024,69.048,66.435,64.14 ,62.123,60.343,58.765,57.357,56.091,54.943,53.893,52.921,52.014,51.159,50.344,49.561,48.804,48.065,47.341,46.628,45.924,45.226,44.534,43.847,43.166,42.491,41.823,41.163,40.513,39.875,39.248,38.637,38.041,37.462,36.902,36.362,35.842,35.343,34.867,34.413,
                                                                 127.86,119.71,112.08,105.02,98.545,92.648,87.326,82.56 ,78.321,74.574,71.273,68.375,65.833,63.601,61.636,59.9  ,58.358,56.978,55.735,54.605,53.567,52.605,51.704,50.853,50.041,49.26 ,48.502,47.763,47.038,46.325,45.619,44.921,44.229,43.543,42.862,42.189,41.524,40.867,40.221,39.587,38.967,38.361,37.772,37.201,36.649,36.117,35.606,35.117,34.65 ,34.205,
                                                                 126.92,118.67,110.99,103.91,97.439,91.581,86.32 ,81.63 ,77.474,73.81 ,70.59 ,67.766,65.289,63.112,61.194,59.496,57.985,56.63 ,55.405,54.289,53.261,52.307,51.411,50.562,49.752,48.97 ,48.213,47.473,46.747,46.032,45.325,44.626,43.934,43.248,42.569,41.897,41.234,40.581,39.939,39.31 ,38.695,38.096,37.514,36.95 ,36.406,35.883,35.381,34.901,34.443,34.007,
                                                                 125.94,117.61,109.9 ,102.81,96.372,90.567,85.376,80.767,76.695,73.113,69.97 ,67.214,64.795,62.668,60.791,59.126,57.641,56.306,55.097,53.992,52.973,52.024,51.131,50.284,49.474,48.692,47.933,47.192,46.464,45.748,45.04 ,44.34 ,43.647,42.962,42.284,41.614,40.953,40.304,39.666,39.042,38.432,37.84 ,37.265,36.709,36.173,35.658,35.164,34.693,34.244,33.818,
                                                                 124.93,116.55,108.81,101.75,95.349,89.609,84.497,79.971,75.983,72.479,69.407,66.712,64.346,62.263,60.422,58.785,57.322,56.005,54.808,53.712,52.699,51.753,50.863,50.017,49.206,48.423,47.662,46.919,46.19 ,45.471,44.762,44.061,43.368,42.683,42.006,41.338,40.68 ,40.034,39.401,38.781,38.178,37.592,37.024,36.475,35.947,35.441,34.956,34.494,34.054,33.636,
                                                                 123.9 ,115.49,107.75,100.72,94.378,88.712,83.682,79.239,75.331,71.902,68.894,66.256,63.936,61.892,60.081,58.469,57.025,55.721,54.535,53.446,52.437,51.494,50.605,49.758,48.946,48.161,47.399,46.653,45.922,45.202,44.491,43.789,43.096,42.411,41.735,41.069,40.414,39.772,39.143,38.529,37.931,37.351,36.791,36.25 ,35.73 ,35.232,34.755,34.302,33.871,33.462,
                                                                 122.87,114.43,106.72,99.733,93.462,87.876,82.929,78.569,74.737,71.376,68.427,65.838,63.56 ,61.549,59.765,58.174,56.745,55.453,54.275,53.191,52.186,51.244,50.354,49.507,48.693,47.906,47.141,46.393,45.659,44.937,44.225,43.523,42.829,42.145,41.47 ,40.806,40.154,39.516,38.891,38.282,37.691,37.118,36.564,36.031,35.519,35.029,34.562,34.117,33.694,33.294,
                                                                 121.84,113.41,105.73,98.801,92.606,87.102,82.238,77.956,74.195,70.895,68    ,65.456,63.214,61.231,59.47 ,57.896,56.48 ,55.198,54.026,52.946,51.942,51.001,50.111,49.261,48.445,47.656,46.888,46.138,45.402,44.678,43.965,43.261,42.567,41.884,41.211,40.549,39.9  ,39.265,38.646,38.043,37.457,36.891,36.344,35.819,35.315,34.834,34.374,33.938,33.524,33.133,
                                                                 120.82,112.41,104.78,97.923,91.808,86.386,81.602,77.394,73.699,70.456,67.608,65.102,62.892,60.935,59.193,57.633,56.228,54.953,53.786,52.708,51.706,50.764,49.872,49.021,48.202,47.41 ,46.64 ,45.887,45.148,44.422,43.708,43.004,42.31 ,41.627,40.956,40.297,39.651,39.021,38.406,37.808,37.229,36.67 ,36.131,35.613,35.117,34.644,34.193,33.765,33.36 ,32.978,
                                                                 119.82,111.45,103.88,97.103,91.068,85.728,81.02 ,76.879,73.244,70.051,67.246,64.775,62.592,60.656,58.93 ,57.383,55.986,54.717,53.553,52.477,51.475,50.532,49.638,48.784,47.963,47.168,46.394,45.639,44.898,44.17 ,43.455,42.75 ,42.057,41.375,40.705,40.049,39.407,38.781,38.171,37.579,37.007,36.454,35.922,35.413,34.925,34.46 ,34.018,33.598,33.202,32.828,
                                                                 118.86,110.54,103.04,96.338,90.385,85.121,80.484,76.407,72.825,69.678,66.909,64.468,62.309,60.391,58.68 ,57.142,55.753,54.488,53.326,52.251,51.248,50.303,49.407,48.551,47.726,46.928,46.152,45.394,44.651,43.921,43.204,42.499,41.807,41.126,40.459,39.805,39.167,38.545,37.941,37.355,36.789,36.244,35.72 ,35.218,34.738,34.281,33.848,33.437,33.048,32.682,
                                                                 117.94,109.68,102.25,95.629,89.753,84.563,79.991,75.971,72.438,69.331,66.595,64.18 ,62.041,60.14 ,58.439,56.911,55.527,54.264,53.104,52.029,51.025,50.078,49.179,48.32 ,47.492,46.691,45.911,45.151,44.405,43.675,42.957,42.252,41.56 ,40.881,40.216,39.566,38.931,38.314,37.715,37.136,36.576,36.038,35.522,35.028,34.556,34.108,33.682,33.28 ,32.9  ,32.542,
                                                                 117.06,108.86,101.52,94.971,89.17 ,84.048,79.537,75.568,72.078,69.007,66.299,63.907,61.786,59.898,58.207,56.685,55.306,54.046,52.886,51.81 ,50.804,49.855,48.953,48.09 ,47.259,46.455,45.673,44.909,44.162,43.43 ,42.712,42.007,41.316,40.638,39.976,39.329,38.699,38.087,37.494,36.921,36.368,35.837,35.328,34.842,34.379,33.939,33.522,33.128,32.756,32.406,
                                                                 116.23,108.11,100.83,94.362,88.631,83.572,79.115,75.193,71.741,68.701,66.02 ,63.647,61.542,59.665,57.982,56.466,55.089,53.83 ,52.671,51.593,50.585,49.633,48.728,47.862,47.028,46.22 ,45.435,44.67 ,43.921,43.187,42.469,41.764,41.074,40.399,39.739,39.096,38.471,37.864,37.277,36.71 ,36.164,35.641,35.14 ,34.662,34.207,33.775,33.366,32.98 ,32.617,32.275,
                                                                 115.45,107.4 ,100.2 ,93.798,88.131,83.129,78.722,74.841,71.424,68.412,65.753,63.398,61.306,59.438,57.762,56.25 ,54.876,53.617,52.457,51.378,50.367,49.412,48.504,47.635,46.797,45.987,45.199,44.431,43.681,42.946,42.227,41.524,40.835,40.162,39.506,38.867,38.246,37.645,37.063,36.503,35.965,35.449,34.955,34.485,34.039,33.615,33.215,32.837,32.481,32.147,
                                                                 114.72,106.74,99.607,93.273,87.667,82.717,78.354,74.51 ,71.124,68.137,65.497,63.157,61.076,59.217,57.547,56.037,54.664,53.406,52.245,51.164,50.15 ,49.192,48.281,47.408,46.567,45.754,44.963,44.193,43.442,42.707,41.988,41.285,40.598,39.928,39.275,38.64 ,38.024,37.429,36.854,36.3  ,35.769,35.26 ,34.775,34.313,33.875,33.459,33.067,32.697,32.35 ,32.024};

        return arrayCurvatureMax[y * arrayWidth + x];
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

    case 4:
    {
        static const std::vector<double> arrayCurvatureMin   = {-0.3,-1.4,-3  ,-4.9,-6.9,-8.7,-10 ,-10.8,-11.2,-11.4,-11.4,-11.3,-11.2,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.1 ,-0.8,-2.1,-3.8,-5.8,-7.8,-9.3,-10.4,-11  ,-11.3,-11.4,-11.3,-11.3,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.4,-11.6,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.3 ,-0.3,-1.3,-2.8,-4.7,-6.7,-8.5,-9.9 ,-10.7,-11.2,-11.3,-11.4,-11.3,-11.2,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.5 ,0   ,-0.7,-1.9,-3.6,-5.7,-7.6,-9.2 ,-10.3,-11  ,-11.3,-11.4,-11.3,-11.3,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.3 ,-0.3,-1.2,-2.7,-4.5,-6.6,-8.4 ,-9.8 ,-10.7,-11.1,-11.3,-11.4,-11.3,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.4 ,0.1 ,-0.6,-1.8,-3.5,-5.5,-7.4 ,-9.1 ,-10.2,-10.9,-11.2,-11.3,-11.3,-11.3,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.5 ,0.3 ,-0.2,-1.1,-2.5,-4.3,-6.4 ,-8.2 ,-9.7 ,-10.6,-11.1,-11.3,-11.3,-11.3,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.3,-11.4,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.5 ,0.4 ,0.1 ,-0.5,-1.7,-3.3,-5.3 ,-7.3 ,-8.9 ,-10.1,-10.8,-11.2,-11.3,-11.3,-11.3,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.5 ,0.4 ,0.3 ,-0.1,-1  ,-2.3,-4.2 ,-6.2 ,-8.1 ,-9.5 ,-10.5,-11  ,-11.3,-11.3,-11.3,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.4 ,0.4 ,0.3 ,0.1 ,-0.5,-1.5,-3.1 ,-5.1 ,-7.1 ,-8.8 ,-10  ,-10.8,-11.2,-11.3,-11.3,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.3,-11.4,-11.6,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.4 ,0.4 ,0.4 ,0.3 ,-0.1,-0.9,-2.2 ,-4   ,-6   ,-7.9 ,-9.4 ,-10.4,-11  ,-11.2,-11.3,-11.3,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.6 ,0.3 ,0.3 ,0.4 ,0.3 ,0.2 ,-0.4,-1.4 ,-2.9 ,-4.9 ,-6.9 ,-8.6 ,-9.9 ,-10.7,-11.1,-11.3,-11.3,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.7 ,0.3 ,0.3 ,0.3 ,0.4 ,0.3 ,-0  ,-0.8 ,-2   ,-3.8 ,-5.8 ,-7.7 ,-9.3 ,-10.3,-10.9,-11.2,-11.3,-11.2,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.3,-11.4,-11.6,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.7 ,0.3 ,0.2 ,0.2 ,0.3 ,0.4 ,0.2 ,-0.3 ,-1.3 ,-2.8 ,-4.7 ,-6.7 ,-8.5 ,-9.8 ,-10.6,-11  ,-11.2,-11.2,-11.2,-11.1,-11  ,-11  ,-11  ,-11  ,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.8 ,0.3 ,0.2 ,0.2 ,0.3 ,0.4 ,0.3 ,0    ,-0.7 ,-1.9 ,-3.6 ,-5.6 ,-7.5 ,-9.1 ,-10.2,-10.8,-11.1,-11.2,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                0.9 ,0.3 ,0.1 ,0.1 ,0.2 ,0.3 ,0.4 ,0.2  ,-0.2 ,-1.1 ,-2.6 ,-4.5 ,-6.5 ,-8.3 ,-9.6 ,-10.5,-11  ,-11.2,-11.2,-11.2,-11.1,-11  ,-11  ,-11  ,-11.1,-11.1,-11.3,-11.4,-11.6,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                1.1 ,0.3 ,0.1 ,0   ,0.1 ,0.3 ,0.4 ,0.3  ,0.1  ,-0.5 ,-1.7 ,-3.4 ,-5.4 ,-7.3 ,-8.9 ,-10.1,-10.7,-11.1,-11.2,-11.2,-11.1,-11.1,-11  ,-11  ,-11  ,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                1.2 ,0.4 ,0.1 ,-0  ,0.1 ,0.2 ,0.3 ,0.4  ,0.3  ,-0.1 ,-1   ,-2.4 ,-4.3 ,-6.3 ,-8.1 ,-9.5 ,-10.4,-10.9,-11.1,-11.2,-11.2,-11.1,-11.1,-11  ,-11  ,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                1.4 ,0.5 ,0.1 ,-0.1,-0  ,0.1 ,0.2 ,0.4  ,0.4  ,0.2  ,-0.4 ,-1.5 ,-3.2 ,-5.1 ,-7.1 ,-8.8 ,-10  ,-10.7,-11  ,-11.2,-11.2,-11.1,-11.1,-11  ,-11  ,-11.1,-11.2,-11.3,-11.4,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                1.7 ,0.5 ,0.1 ,-0.1,-0.1,0   ,0.2 ,0.3  ,0.4  ,0.3  ,-0   ,-0.9 ,-2.2 ,-4   ,-6   ,-7.9 ,-9.4 ,-10.3,-10.9,-11.1,-11.2,-11.2,-11.1,-11.1,-11  ,-11.1,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                1.9 ,0.7 ,0.1 ,-0.1,-0.1,-0.1,0.1 ,0.2  ,0.4  ,0.4  ,0.2  ,-0.3 ,-1.4 ,-3   ,-4.9 ,-6.9 ,-8.6 ,-9.8 ,-10.6,-11  ,-11.1,-11.2,-11.1,-11.1,-11.1,-11.1,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                2.3 ,0.8 ,0.1 ,-0.1,-0.2,-0.1,0   ,0.2  ,0.3  ,0.4  ,0.4  ,0.1  ,-0.7 ,-2   ,-3.8 ,-5.8 ,-7.7 ,-9.2 ,-10.2,-10.8,-11.1,-11.2,-11.2,-11.1,-11.1,-11.1,-11.1,-11.2,-11.3,-11.4,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                2.6 ,0.9 ,0.2 ,-0.1,-0.2,-0.2,-0.1,0.1  ,0.2  ,0.4  ,0.5  ,0.3  ,-0.2 ,-1.2 ,-2.8 ,-4.7 ,-6.7 ,-8.5 ,-9.7 ,-10.5,-10.9,-11.1,-11.2,-11.1,-11.1,-11.1,-11.1,-11.1,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                3   ,1.1 ,0.3 ,-0.1,-0.2,-0.2,-0.1,-0   ,0.2  ,0.3  ,0.5  ,0.5  ,0.2  ,-0.6 ,-1.8 ,-3.6 ,-5.6 ,-7.6 ,-9.1 ,-10.1,-10.7,-11  ,-11.1,-11.2,-11.1,-11.1,-11.1,-11.1,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                3.5 ,1.4 ,0.3 ,-0.1,-0.3,-0.3,-0.2,-0.1 ,0.1  ,0.3  ,0.4  ,0.5  ,0.4  ,-0.1 ,-1.1 ,-2.6 ,-4.5 ,-6.5 ,-8.3 ,-9.6 ,-10.5,-10.9,-11.1,-11.2,-11.1,-11.1,-11.1,-11.1,-11.2,-11.3,-11.4,-11.6,-11.7,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                4   ,1.6 ,0.4 ,-0.1,-0.3,-0.3,-0.3,-0.1 ,0    ,0.2  ,0.4  ,0.5  ,0.5  ,0.2  ,-0.4 ,-1.7 ,-3.4 ,-5.4 ,-7.4 ,-9   ,-10.1,-10.7,-11  ,-11.1,-11.2,-11.1,-11.1,-11.1,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                4.6 ,1.9 ,0.6 ,-0  ,-0.3,-0.3,-0.3,-0.2 ,-0.1 ,0.1  ,0.3  ,0.5  ,0.6  ,0.5  ,0    ,-0.9 ,-2.4 ,-4.3 ,-6.3 ,-8.1 ,-9.5 ,-10.4,-10.9,-11.1,-11.1,-11.1,-11.1,-11.1,-11.2,-11.2,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                5.2 ,2.2 ,0.7 ,0   ,-0.3,-0.3,-0.3,-0.3 ,-0.1 ,0    ,0.2  ,0.4  ,0.5  ,0.6  ,0.3  ,-0.3 ,-1.5 ,-3.2 ,-5.2 ,-7.2 ,-8.8 ,-10  ,-10.6,-11  ,-11.1,-11.2,-11.1,-11.1,-11.1,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                6   ,2.6 ,0.9 ,0.1 ,-0.2,-0.4,-0.4,-0.3 ,-0.2 ,-0.1 ,0.1  ,0.3  ,0.5  ,0.6  ,0.5  ,0.1  ,-0.8 ,-2.2 ,-4.1 ,-6.1 ,-8   ,-9.4 ,-10.3,-10.8,-11.1,-11.1,-11.1,-11.1,-11.1,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                6.8 ,3   ,1.1 ,0.2 ,-0.2,-0.4,-0.4,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.6  ,0.4  ,-0.2 ,-1.4 ,-3   ,-5   ,-7   ,-8.7 ,-9.9 ,-10.6,-11  ,-11.1,-11.2,-11.1,-11.1,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                7.7 ,3.5 ,1.3 ,0.3 ,-0.2,-0.4,-0.4,-0.4 ,-0.3 ,-0.2 ,-0   ,0.1  ,0.3  ,0.5  ,0.6  ,0.6  ,0.2  ,-0.7 ,-2.1 ,-3.9 ,-6   ,-7.8 ,-9.3 ,-10.2,-10.8,-11  ,-11.1,-11.2,-11.2,-11.2,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                8.6 ,4   ,1.6 ,0.4 ,-0.1,-0.3,-0.4,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.6  ,0.4  ,-0.1 ,-1.2 ,-2.9 ,-4.8 ,-6.8 ,-8.5 ,-9.8 ,-10.5,-10.9,-11.1,-11.2,-11.2,-11.2,-11.2,-11.3,-11.4,-11.6,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                9.7 ,4.6 ,1.9 ,0.6 ,-0.1,-0.3,-0.4,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0   ,0.1  ,0.3  ,0.5  ,0.6  ,0.6  ,0.2  ,-0.5 ,-1.9 ,-3.7 ,-5.8 ,-7.7 ,-9.2 ,-10.2,-10.7,-11  ,-11.1,-11.2,-11.2,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                10.8,5.3 ,2.2 ,0.7 ,0   ,-0.3,-0.4,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0.1  ,0.2  ,0.4  ,0.6  ,0.6  ,0.5  ,-0   ,-1.1 ,-2.7 ,-4.6 ,-6.7 ,-8.4 ,-9.7 ,-10.5,-10.9,-11.1,-11.2,-11.2,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                12  ,6   ,2.6 ,0.9 ,0.1 ,-0.2,-0.4,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0   ,0.1  ,0.3  ,0.5  ,0.6  ,0.6  ,0.3  ,-0.4 ,-1.7 ,-3.5 ,-5.6 ,-7.5 ,-9   ,-10.1,-10.7,-11  ,-11.1,-11.2,-11.2,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                13.3,6.8 ,3   ,1.1 ,0.2 ,-0.2,-0.4,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.7  ,0.5  ,0.1  ,-0.9 ,-2.5 ,-4.4 ,-6.5 ,-8.2 ,-9.6 ,-10.4,-10.9,-11.1,-11.2,-11.2,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                14.7,7.7 ,3.5 ,1.4 ,0.3 ,-0.2,-0.4,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.3 ,-0.2 ,-0   ,0.1  ,0.3  ,0.5  ,0.7  ,0.7  ,0.4  ,-0.3 ,-1.6 ,-3.3 ,-5.4 ,-7.3 ,-8.9 ,-10  ,-10.6,-11  ,-11.1,-11.2,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                16.2,8.7 ,4.1 ,1.6 ,0.4 ,-0.1,-0.3,-0.4 ,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.7  ,0.6  ,0.1  ,-0.8 ,-2.3 ,-4.2 ,-6.3 ,-8.1 ,-9.5 ,-10.3,-10.8,-11.1,-11.2,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                17.7,9.8 ,4.7 ,1.9 ,0.6 ,-0  ,-0.3,-0.4 ,-0.4 ,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0   ,0.1  ,0.3  ,0.5  ,0.7  ,0.7  ,0.4  ,-0.2 ,-1.4 ,-3.2 ,-5.2 ,-7.2 ,-8.8 ,-9.9 ,-10.6,-11  ,-11.1,-11.2,-11.3,-11.3,-11.4,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                19.3,10.9,5.4 ,2.3 ,0.8 ,0.1 ,-0.3,-0.4 ,-0.4 ,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.7  ,0.6  ,0.2  ,-0.7 ,-2.1 ,-4.1 ,-6.1 ,-7.9 ,-9.3 ,-10.3,-10.8,-11.1,-11.2,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                20.9,12.1,6.1 ,2.7 ,1   ,0.1 ,-0.2,-0.4 ,-0.4 ,-0.5 ,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0   ,0.1  ,0.3  ,0.5  ,0.7  ,0.7  ,0.5  ,-0.1 ,-1.3 ,-3   ,-5   ,-7   ,-8.6 ,-9.8 ,-10.5,-10.9,-11.1,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                22.6,13.5,6.9 ,3.1 ,1.2 ,0.2 ,-0.2,-0.4 ,-0.4 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.7  ,0.7  ,0.3  ,-0.6 ,-2   ,-3.9 ,-5.9 ,-7.8 ,-9.2 ,-10.2,-10.8,-11  ,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                24.3,14.9,7.8 ,3.6 ,1.4 ,0.4 ,-0.1,-0.3 ,-0.4 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0.1  ,0.3  ,0.5  ,0.7  ,0.7  ,0.6  ,-0   ,-1.1 ,-2.8 ,-4.8 ,-6.8 ,-8.5 ,-9.7 ,-10.5,-10.9,-11.1,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                26  ,16.3,8.8 ,4.2 ,1.7 ,0.5 ,-0.1,-0.3 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.7  ,0.7  ,0.4  ,-0.4 ,-1.8 ,-3.7 ,-5.7 ,-7.6 ,-9.1 ,-10.1,-10.7,-11  ,-11.2,-11.3,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                27.7,17.9,9.9 ,4.8 ,2   ,0.6 ,-0  ,-0.3 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0.1  ,0.3  ,0.5  ,0.7  ,0.8  ,0.6  ,0.1  ,-1   ,-2.6 ,-4.6 ,-6.6 ,-8.4 ,-9.6 ,-10.4,-10.9,-11.1,-11.3,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                29.3,19.5,11  ,5.4 ,2.4 ,0.8 ,0.1 ,-0.3 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.8  ,0.7  ,0.4  ,-0.3 ,-1.7 ,-3.5 ,-5.5 ,-7.5 ,-9   ,-10  ,-10.7,-11  ,-11.2,-11.3,-11.5,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                30.9,21.1,12.3,6.2 ,2.7 ,1   ,0.2 ,-0.2 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0.1  ,0.3  ,0.5  ,0.7  ,0.8  ,0.7  ,0.2  ,-0.8 ,-2.4 ,-4.4 ,-6.4 ,-8.2 ,-9.5 ,-10.4,-10.9,-11.1,-11.3,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                32.5,22.8,13.6,7   ,3.2 ,1.2 ,0.3 ,-0.2 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0    ,0.2  ,0.4  ,0.6  ,0.8  ,0.8  ,0.5  ,-0.2 ,-1.5 ,-3.3 ,-5.3 ,-7.3 ,-8.9 ,-10  ,-10.6,-11  ,-11.2,-11.4,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                34  ,24.5,15  ,7.9 ,3.7 ,1.5 ,0.4 ,-0.1 ,-0.3 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.2 ,-0.1 ,0.1  ,0.3  ,0.5  ,0.7  ,0.8  ,0.7  ,0.3  ,-0.7 ,-2.2 ,-4.2 ,-6.3 ,-8.1 ,-9.4 ,-10.3,-10.8,-11.1,-11.3,-11.5,-11.6,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,
                                                                35.4,26.1,16.5,8.9 ,4.2 ,1.7 ,0.5 ,-0.1 ,-0.3 ,-0.4 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.5 ,-0.4 ,-0.4 ,-0.3 ,-0.3 ,-0.1 ,-0   ,0.2  ,0.4  ,0.6  ,0.8  ,0.8  ,0.6  ,-0.1 ,-1.3 ,-3.1 ,-5.1 ,-7.1 ,-8.7 ,-9.9 ,-10.6,-11  ,-11.2,-11.4,-11.6,-11.7,-11.7,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8,-11.8};

        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 7:
    {
        static const std::vector<double> arrayCurvatureMin   = {6.6983,5.3967,4.2588,3.2846,2.4655,1.7872,1.2317,0.78085,0.4169 ,0.12404,-0.1115  ,-0.30161   ,-0.45669  ,-0.58645  ,-0.70118 ,-0.8138   ,-0.94326  ,-1.1196  ,-1.3896  ,-1.8169  ,-2.4607  ,-3.3181   ,-4.2713  ,-5.1336   ,-5.7779  ,-6.1929 ,-6.4338  ,-6.5632  ,-6.6263  ,-6.6503  ,-6.6492  ,-6.6289  ,-6.5894  ,-6.5265  ,-6.4313  ,-6.2897  ,-6.0812  ,-5.778   ,-5.3461  ,-4.75    ,-3.9654  ,-2.9962  ,-1.8902  ,-0.73727,0.35682  ,1.3075   ,2.0743   ,2.6586   ,3.0874  ,3.3961  ,
                                                                7.0648,5.7174,4.5352,3.5197,2.6634,1.9521,1.3685,0.89375,0.50997,0.20088,-0.047706,-0.2479    ,-0.41007  ,-0.54346  ,-0.65701 ,-0.76085  ,-0.86861  ,-1.0012  ,-1.1924  ,-1.4936  ,-1.9708  ,-2.673    ,-3.5696  ,-4.5142   ,-5.3235  ,-5.9013 ,-6.2608  ,-6.4636  ,-6.5683  ,-6.6149  ,-6.6261  ,-6.6129  ,-6.5784  ,-6.5205  ,-6.432   ,-6.3003  ,-6.1067  ,-5.8251  ,-5.4227  ,-4.8637  ,-4.1198  ,-3.1871  ,-2.1026  ,-0.94827,0.17033  ,1.161    ,1.9724   ,2.5976   ,3.0596  ,3.3927  ,
                                                                7.4474,6.0531,4.8254,3.7672,2.872 ,2.1266,1.5133,1.0135 ,0.60878,0.28242,0.019824 ,-0.19148   ,-0.36202  ,-0.50088  ,-0.61631 ,-0.71676  ,-0.81251  ,-0.91826 ,-1.0574  ,-1.2678  ,-1.6057  ,-2.1368   ,-2.8958  ,-3.8207   ,-4.7427  ,-5.4919 ,-6.0042  ,-6.3124  ,-6.4805  ,-6.5624  ,-6.5928  ,-6.5901  ,-6.5622  ,-6.5098  ,-6.4276  ,-6.305   ,-6.1249  ,-5.8632  ,-5.4884  ,-4.9649  ,-4.2616  ,-3.3673  ,-2.3088  ,-1.1589 ,-0.021008,1.0066   ,1.8619   ,2.5289   ,3.0253  ,3.3842  ,
                                                                7.8466,6.4043,5.1298,4.0276,3.0922,2.3109,1.6667,1.1405 ,0.71371,0.36907,0.091525 ,-0.13181   ,-0.31174  ,-0.45741  ,-0.57677 ,-0.67735  ,-0.76735  ,-0.85747 ,-0.96392 ,-1.1131  ,-1.3473  ,-1.727    ,-2.3151  ,-3.1266   ,-4.067   ,-4.9534 ,-5.6376  ,-6.0866  ,-6.3475  ,-6.4839  ,-6.5442  ,-6.5578  ,-6.5393  ,-6.4935  ,-6.4178  ,-6.3037  ,-6.1361  ,-5.8926  ,-5.5436  ,-5.0541  ,-4.3908  ,-3.5364  ,-2.508   ,-1.3684 ,-0.21669 ,0.84417  ,1.7422   ,2.4515   ,2.9835  ,3.3696  ,
                                                                8.2629,6.7716,5.4492,4.3014,3.3243,2.5058,1.8291,1.2753 ,0.82519,0.46122,0.16779  ,-0.068435  ,-0.25863  ,-0.41214  ,-0.53687 ,-0.63988  ,-0.72819  ,-0.81    ,-0.89674 ,-1.0065  ,-1.1691  ,-1.4318   ,-1.8582  ,-2.5048   ,-3.3619  ,-4.3037 ,-5.1431  ,-5.7598  ,-6.1485  ,-6.3659  ,-6.4727  ,-6.5118  ,-6.5075  ,-6.4705  ,-6.4021  ,-6.2962  ,-6.1402  ,-5.9137  ,-5.5888  ,-5.1316  ,-4.5076  ,-3.6942  ,-2.6995  ,-1.5758 ,-0.41617 ,0.67375  ,1.6128   ,2.3645   ,2.9334  ,3.348   ,
                                                                8.6967,7.1556,5.784 ,4.5894,3.5689,2.7117,2.0011,1.4182 ,0.94363,0.55926,0.249    ,-0.00096952,-0.20223  ,-0.36443  ,-0.49558 ,-0.60261  ,-0.69186  ,-0.77004 ,-0.84559 ,-0.93102 ,-1.0467  ,-1.2265   ,-1.5225  ,-1.9999   ,-2.7045  ,-3.5976 ,-4.5263  ,-5.3093  ,-5.8577  ,-6.1897  ,-6.3668  ,-6.4456  ,-6.4632  ,-6.4389  ,-6.3794  ,-6.2823  ,-6.1373  ,-5.9267  ,-5.6245  ,-5.1981  ,-4.6123  ,-3.8407  ,-2.8828  ,-1.7802 ,-0.61884 ,0.49539  ,1.473    ,2.267    ,2.874   ,3.3186  ,
                                                                9.1486,7.5568,6.1349,4.892 ,3.8267,2.9291,2.1832,1.5698 ,1.0695 ,0.66358,0.33552  ,0.070962   ,-0.14213  ,-0.31376  ,-0.45218 ,-0.56433  ,-0.65626  ,-0.73385 ,-0.80369 ,-0.87469 ,-0.96078 ,-1.0851   ,-1.2859  ,-1.6199   ,-2.1518  ,-2.9116 ,-3.829   ,-4.7303  ,-5.4498  ,-5.9309  ,-6.2097  ,-6.3492  ,-6.4009  ,-6.3958  ,-6.3484  ,-6.2612  ,-6.1275  ,-5.932   ,-5.6513  ,-5.2542  ,-4.7057  ,-3.9761  ,-3.0575  ,-1.9812 ,-0.82408 ,0.30927  ,1.3226   ,2.1584   ,2.8043  ,3.2805  ,
                                                                9.6188,7.9757,6.5024,5.2098,4.0983,3.1588,2.3759,1.7306 ,1.2032 ,0.77462,0.42774  ,0.14773    ,-0.077965 ,-0.25972  ,-0.40612 ,-0.52422  ,-0.62004  ,-0.699   ,-0.76662 ,-0.82959 ,-0.89758 ,-0.98626  ,-1.1221  ,-1.3478   ,-1.7244  ,-2.3131 ,-3.1229  ,-4.0508  ,-4.9116  ,-5.5632  ,-5.9789  ,-6.208   ,-6.3121  ,-6.3366  ,-6.3068  ,-6.2321  ,-6.1104  ,-5.9298  ,-5.6696  ,-5.3005  ,-4.7882  ,-4.1007  ,-3.2236  ,-2.178  ,-1.0312  ,0.11575  ,1.1614   ,2.038    ,2.7237  ,3.2331  ,
                                                                10.108,8.4127,6.8871,5.5436,4.3842,3.4012,2.5798,1.9011 ,1.3453 ,0.8928 ,0.52607  ,0.22969    ,-0.0093609,-0.20193  ,-0.35693 ,-0.48166  ,-0.58223  ,-0.66386 ,-0.73151 ,-0.79052 ,-0.84788 ,-0.91427  ,-1.0075  ,-1.1577   ,-1.4127  ,-1.8363 ,-2.4824  ,-3.334   ,-4.258   ,-5.0671  ,-5.6483  ,-6.0016  ,-6.1841  ,-6.2542  ,-6.2508  ,-6.1932  ,-6.0855  ,-5.9203  ,-5.68    ,-5.3379  ,-4.8607  ,-4.215   ,-3.3809  ,-2.37   ,-1.2394  ,-0.084656,0.98943  ,1.9055   ,2.6316  ,3.1756  ,
                                                                10.616,8.8684,7.2895,5.8938,4.6852,3.6571,2.7955,2.0819 ,1.4962 ,1.0186 ,0.63092  ,0.31726    ,0.064045  ,-0.14001  ,-0.30422 ,-0.43617  ,-0.54213  ,-0.62735 ,-0.69648 ,-0.75411 ,-0.80559 ,-0.85841  ,-0.92454 ,-1.0244   ,-1.1922  ,-1.4807 ,-1.9551  ,-2.6573  ,-3.5406  ,-4.4459  ,-5.1944  ,-5.7048  ,-5.999   ,-6.1379  ,-6.1747  ,-6.1415  ,-6.0516  ,-5.9032  ,-5.6829  ,-5.3671  ,-4.924   ,-4.3196  ,-3.5296  ,-2.5569 ,-1.4478  ,-0.29126 ,0.80688  ,1.7606   ,2.5273  ,3.1075  ,
                                                                11.144,9.3432,7.7101,6.2611,5.0018,3.9269,3.0236,2.2735 ,1.6565 ,1.1525 ,0.74274  ,0.41082    ,0.14263   ,-0.073604 ,-0.24761 ,-0.38731  ,-0.49922  ,-0.58867 ,-0.66025 ,-0.71815 ,-0.76674 ,-0.81155  ,-0.86078 ,-0.92804  ,-1.0367  ,-1.2256 ,-1.5522  ,-2.0804  ,-2.8352  ,-3.7378  ,-4.6104  ,-5.2921  ,-5.733   ,-5.9717  ,-6.0692  ,-6.0725  ,-6.0065  ,-5.8779  ,-5.6786  ,-5.3886  ,-4.9788  ,-4.4151  ,-3.6698  ,-2.7381 ,-1.6557  ,-0.50321 ,0.61424  ,1.6032   ,2.4107  ,3.0285  ,
                                                                11.691,9.8373,8.1496,6.6462,5.3346,4.2115,3.2647,2.4766 ,1.8268 ,1.295  ,0.862    ,0.51082    ,0.22681   ,-0.0023223,-0.18672 ,-0.33469  ,-0.45301  ,-0.54721 ,-0.62191 ,-0.68114 ,-0.72875 ,-0.76907  ,-0.80793 ,-0.8545   ,-0.92438 ,-1.0444 ,-1.2583  ,-1.6275  ,-2.2115  ,-3.013   ,-3.9213  ,-4.7491  ,-5.3603  ,-5.7343  ,-5.9208  ,-5.9784  ,-5.9465  ,-5.8427  ,-5.6665  ,-5.4029  ,-5.0257  ,-4.502   ,-3.8017  ,-2.9133 ,-1.8622  ,-0.7195  ,0.41217  ,1.4336   ,2.2814  ,2.938   ,
                                                                12.258,10.351,8.6083,7.0494,5.6844,4.5115,3.5196,2.6917 ,2.0076 ,1.4467 ,0.9892   ,0.61771    ,0.31699   ,0.074232  ,-0.12117 ,-0.27792  ,-0.40309  ,-0.50247 ,-0.58079 ,-0.64206 ,-0.68989 ,-0.72792  ,-0.76055 ,-0.79415  ,-0.83909 ,-0.91335,-1.0477  ,-1.2909  ,-1.7071  ,-2.3475  ,-3.1879  ,-4.0878  ,-4.8606  ,-5.4     ,-5.7103  ,-5.8477  ,-5.8654  ,-5.7949  ,-5.6458  ,-5.4098  ,-5.0651  ,-4.5808  ,-3.9256  ,-3.0821 ,-2.0662  ,-0.93899 ,0.20156  ,1.2521   ,2.1394  ,2.8359  ,
                                                                12.845,10.885,9.0867,7.4715,6.0518,4.8275,3.7889,2.9196 ,2.1996 ,1.6081 ,1.1249   ,0.73198    ,0.41363   ,0.15648   ,-0.050557,-0.21659  ,-0.34906  ,-0.45399 ,-0.53633 ,-0.60015 ,-0.64898 ,-0.68613  ,-0.71515 ,-0.74066  ,-0.76975 ,-0.81432,-0.89514 ,-1.0473  ,-1.3248  ,-1.7923  ,-2.488   ,-3.3572  ,-4.2349  ,-4.9453  ,-5.4134  ,-5.6634  ,-5.7532  ,-5.7292  ,-5.614   ,-5.4086  ,-5.0971  ,-4.6519  ,-4.0414  ,-3.2441 ,-2.2669  ,-1.1604  ,-0.016516,1.0594   ,1.9849  ,2.722   ,
                                                                13.452,11.44 ,9.5852,7.913 ,6.4373,5.1602,4.0732,3.1609 ,2.4034 ,1.7798 ,1.2696   ,0.85417    ,0.51722   ,0.24487   ,0.025542 ,-0.15032  ,-0.29051  ,-0.40137 ,-0.48806 ,-0.55482 ,-0.60521 ,-0.64238  ,-0.66945 ,-0.69001  ,-0.70907 ,-0.73465,-0.78051 ,-0.87066 ,-1.045   ,-1.362   ,-1.8846  ,-2.6326  ,-3.519   ,-4.3619  ,-5.0047  ,-5.4029  ,-5.5952  ,-5.6373  ,-5.5669  ,-5.3974  ,-5.1211  ,-4.715   ,-4.1492  ,-3.3989 ,-2.4632  ,-1.3825  ,-0.24082 ,0.85613  ,1.818   ,2.5961  ,
                                                                14.078,12.015,10.104,8.3744,6.8417,5.5103,4.3733,3.4163 ,2.6197 ,1.9625 ,1.4239   ,0.98481    ,0.62826   ,0.33989   ,0.10757  ,-0.078671 ,-0.22704  ,-0.3442  ,-0.43557 ,-0.5056  ,-0.55796 ,-0.59579  ,-0.62199 ,-0.63962  ,-0.65242 ,-0.66597,-0.68945 ,-0.73884 ,-0.84184 ,-1.0434  ,-1.4055  ,-1.9861  ,-2.7813  ,-3.6718  ,-4.4686  ,-5.0406  ,-5.3708  ,-5.5062  ,-5.4975  ,-5.3727  ,-5.1354  ,-4.7697  ,-4.2486  ,-3.5459 ,-2.6542  ,-1.6038  ,-0.46999 ,0.64332  ,1.6391  ,2.4583  ,
                                                                14.724,12.61 ,10.644,8.8561,7.2655,5.8784,4.69  ,3.6866 ,2.8492 ,2.1569 ,1.5886   ,1.1245     ,0.74731   ,0.44202   ,0.19601  ,-0.0012127,-0.15823  ,-0.28208 ,-0.37848 ,-0.45209 ,-0.50677 ,-0.54574  ,-0.57186 ,-0.58789  ,-0.5969  ,-0.60296,-0.61235 ,-0.63569 ,-0.69155 ,-0.81171 ,-1.0462  ,-1.4593  ,-2.0991  ,-2.9337  ,-3.8147  ,-4.5555  ,-5.0546  ,-5.3174  ,-5.3944  ,-5.3284  ,-5.1371  ,-4.8146  ,-4.339   ,-3.6846 ,-2.8389  ,-1.8229  ,-0.70252 ,0.42206  ,1.4489  ,2.3086  ,
                                                                15.39 ,13.226,11.205,9.3587,7.7093,6.2654,5.0239,3.9724 ,3.0927 ,2.3637 ,1.7641   ,1.2738     ,0.87494   ,0.55182   ,0.29133  ,0.082507  ,-0.083685 ,-0.21465 ,-0.31642 ,-0.39395 ,-0.45129 ,-0.49183  ,-0.51847 ,-0.53388  ,-0.54077 ,-0.54235,-0.54311 ,-0.55024 ,-0.57603 ,-0.64198 ,-0.78433 ,-1.0579  ,-1.5271  ,-2.2253  ,-3.0889  ,-3.9462  ,-4.6226  ,-5.0468  ,-5.2411  ,-5.2548  ,-5.1211  ,-4.8471  ,-4.4191  ,-3.814  ,-3.0162  ,-2.0385  ,-0.93686 ,0.19362  ,1.2479  ,2.1471  ,
                                                                16.074,13.862,11.787,9.8824,8.1737,6.6718,5.3758,4.2746 ,3.3508 ,2.5836 ,1.9513   ,1.4335     ,1.0117    ,0.66982   ,0.39405  ,0.17295   ,-0.0029728,-0.14152 ,-0.24909 ,-0.3309  ,-0.39128 ,-0.4338   ,-0.4615  ,-0.47707  ,-0.48308 ,-0.48225,-0.47803 ,-0.47544 ,-0.4827  ,-0.51407 ,-0.59432 ,-0.76438 ,-1.0833  ,-1.6125  ,-2.3653  ,-3.245   ,-4.0643  ,-4.6687  ,-5.0153  ,-5.1373  ,-5.0789  ,-4.8628  ,-4.4866  ,-3.9327 ,-3.1848  ,-2.2489  ,-1.1713  ,-0.040579,1.0371  ,1.9742  ,
                                                                16.777,14.518,12.39 ,10.428,8.6592,7.0984,5.7464,4.594  ,3.6244 ,2.8173 ,2.1508   ,1.6042     ,1.1584    ,0.79662   ,0.5047   ,0.27058   ,0.084308  ,-0.062357,-0.17618 ,-0.26271 ,-0.32656 ,-0.37152  ,-0.40082 ,-0.41726  ,-0.42337 ,-0.42166,-0.41497 ,-0.40701 ,-0.40341 ,-0.4136  ,-0.45404 ,-0.55311 ,-0.75655 ,-1.1265  ,-1.7177  ,-2.518   ,-3.3984  ,-4.1656  ,-4.6907  ,-4.9553  ,-4.9977  ,-4.8544  ,-4.5376  ,-4.0386 ,-3.3431  ,-2.4525  ,-1.4042  ,-0.279   ,0.81745 ,1.7902  ,
                                                                17.498,15.195,13.014,10.995,9.1662,7.5456,6.1365,4.9313 ,3.9144 ,3.0656 ,2.3634   ,1.7866     ,1.3154    ,0.93279   ,0.62379  ,0.37586   ,0.17855   ,0.023161 ,-0.097461,-0.18923 ,-0.25705 ,-0.30498  ,-0.33647 ,-0.35445  ,-0.36154 ,-0.36016,-0.35282 ,-0.34245 ,-0.33307 ,-0.33101 ,-0.34707 ,-0.40016 ,-0.52261 ,-0.7649  ,-1.1905  ,-1.8428  ,-2.68    ,-3.5439  ,-4.2454  ,-4.6835  ,-4.8592  ,-4.8109  ,-4.566   ,-4.1282 ,-3.489   ,-2.6475  ,-1.6336  ,-0.51994 ,0.59012 ,1.5957  ,
                                                                18.236,15.891,13.66 ,11.584,9.6953,8.0142,6.5468,5.2873 ,4.2214 ,3.3295 ,2.59     ,1.9814     ,1.4836    ,1.0789    ,0.75185  ,0.48924   ,0.28011   ,0.11531  ,-0.012755,-0.11037 ,-0.18275 ,-0.23427  ,-0.26858 ,-0.28883  ,-0.29772 ,-0.2977 ,-0.2911  ,-0.28041 ,-0.26865 ,-0.26021 ,-0.26219 ,-0.28692 ,-0.3561  ,-0.50629 ,-0.79233 ,-1.2767  ,-1.9861  ,-2.8457  ,-3.6749  ,-4.297   ,-4.6392  ,-4.7159  ,-4.5621  ,-4.1963 ,-3.6193  ,-2.8315  ,-1.8574  ,-0.76152 ,0.35647 ,1.3913  ,
                                                                18.99 ,16.606,14.326,12.195,10.247,8.5048,6.978 ,5.6629 ,4.5463 ,3.6096 ,2.8312   ,2.1893     ,1.6636    ,1.2357    ,0.88939  ,0.61112   ,0.38931   ,0.2143   ,0.078054 ,-0.026107,-0.10376 ,-0.15956  ,-0.19743 ,-0.22069  ,-0.2322  ,-0.23446,-0.22974 ,-0.22025 ,-0.20845 ,-0.19748 ,-0.19208 ,-0.20026 ,-0.23614 ,-0.32452 ,-0.50642 ,-0.8403  ,-1.3843  ,-2.1431  ,-3.0077  ,-3.7829  ,-4.3118  ,-4.5473  ,-4.5117  ,-4.2345 ,-3.7291  ,-3.0015  ,-2.0732  ,-1.0016  ,0.1181  ,1.1779  ,
                                                                19.759,17.339,15.013,12.828,10.821,9.0178,7.4308,6.0587 ,4.89   ,3.9068 ,3.0878   ,2.4112     ,1.8561    ,1.4035    ,1.0369   ,0.74189   ,0.50641   ,0.32029  ,0.17501  ,0.063475 ,-0.020245,-0.081122 ,-0.12333 ,-0.15041  ,-0.16536 ,-0.17075,-0.16885 ,-0.16175 ,-0.15152 ,-0.14054 ,-0.13211 ,-0.13142 ,-0.14753 ,-0.19661 ,-0.30694 ,-0.52405 ,-0.90873 ,-1.5107  ,-2.3073  ,-3.1565  ,-3.8579  ,-4.2791  ,-4.3947  ,-4.2302 ,-3.8112  ,-3.1529  ,-2.2777  ,-1.2377  ,-0.12312,0.95675 ,
                                                                20.542,18.09 ,15.721,13.483,11.418,9.5537,7.9059,6.4756 ,5.2533 ,4.2219 ,3.3607   ,2.6477     ,2.0616    ,1.5831    ,1.1948   ,0.88186   ,0.63163   ,0.43336  ,0.27807  ,0.15823  ,0.067538 ,0.00070469,-0.046706,-0.078429 ,-0.097617,-0.10694,-0.1087  ,-0.10489 ,-0.097372,-0.088085,-0.079381,-0.074717,-0.079867,-0.10514 ,-0.16907 ,-0.30374 ,-0.559   ,-0.99602 ,-1.6512  ,-2.4698  ,-3.281   ,-3.8883  ,-4.186   ,-4.1656 ,-3.8547  ,-3.2792  ,-2.4667  ,-1.4668  ,-0.36492,0.72909 ,
                                                                21.338,18.858,16.448,14.16 ,12.038,10.113,8.4038,6.9144 ,5.637  ,4.5558 ,3.6507   ,2.8995     ,2.2809    ,1.7748    ,1.3635   ,1.0313    ,0.76509   ,0.55352  ,0.38713  ,0.25794  ,0.15928  ,0.085543  ,0.032034 ,-0.0051845,-0.029413,-0.04342,-0.049535,-0.04974 ,-0.045776,-0.039292,-0.032084,-0.026519,-0.026307,-0.037931,-0.07318 ,-0.15329 ,-0.31435 ,-0.61003 ,-1.099   ,-1.7988  ,-2.6197  ,-3.3682  ,-3.8603  ,-4.0171 ,-3.8434  ,-3.3707  ,-2.6342  ,-1.6846  ,-0.60443,0.49678 ,
                                                                22.146,19.641,17.195,14.859,12.681,10.696,8.9251,7.3756 ,6.0418 ,4.9092 ,3.9584   ,3.1674     ,2.5145    ,1.9792    ,1.5433   ,1.1904    ,0.90686   ,0.68072  ,0.50201  ,0.36234  ,0.25464  ,0.17299   ,0.11246  ,0.068896  ,0.038862 ,0.019488,0.0083933,0.0035916,0.0034026,0.0063378,0.010937 ,0.015498 ,0.017585 ,0.013134 ,-0.005197,-0.050923,-0.14826 ,-0.33741 ,-0.6748  ,-1.213   ,-1.9445  ,-2.744   ,-3.4034  ,-3.7586 ,-3.7559  ,-3.4132  ,-2.7713  ,-1.8856  ,-0.8379 ,0.26218 ,
                                                                22.965,20.439,17.959,15.579,13.347,11.303,9.4703,7.8599 ,6.4684 ,5.2828 ,4.2845   ,3.4518     ,2.7628    ,2.1967    ,1.7344   ,1.3593    ,1.0569    ,0.81481  ,0.62249  ,0.47111  ,0.35324  ,0.26266   ,0.19417  ,0.14343   ,0.10687  ,0.081508,0.064902 ,0.05504  ,0.05026  ,0.049155 ,0.050454 ,0.052848 ,0.054674 ,0.053364 ,0.044402 ,0.01943  ,-0.037022,-0.15239 ,-0.37083 ,-0.74982 ,-1.331   ,-2.0768  ,-2.8275  ,-3.3702 ,-3.5662  ,-3.3867  ,-2.8653  ,-2.0618  ,-1.0602 ,0.028508,
                                                                23.792,21.251,18.741,16.319,14.036,11.933,10.04 ,8.3678 ,6.9174 ,5.6773 ,4.6297   ,3.7535     ,3.0265    ,2.4275    ,1.9372   ,1.5381    ,1.2152    ,0.95565  ,0.74832  ,0.58395  ,0.45475  ,0.35419   ,0.27682  ,0.21812   ,0.17436  ,0.14245 ,0.11988  ,0.10459  ,0.094904 ,0.089437 ,0.087018 ,0.086552 ,0.086834 ,0.086198 ,0.081881 ,0.068825 ,0.037513 ,-0.029677,-0.16354 ,-0.41169 ,-0.83015 ,-1.4443  ,-2.1816  ,-2.8528 ,-3.2507  ,-3.2665  ,-2.8982  ,-2.2016  ,-1.264  ,-0.19977,
                                                                24.627,22.074,19.54 ,17.08 ,14.748,12.588,10.633,8.8998 ,7.3894 ,6.0933 ,4.9946   ,4.0727     ,3.3057    ,2.672     ,2.1516   ,1.7267    ,1.3817    ,1.103    ,0.87923  ,0.70056  ,0.55886  ,0.44728   ,0.36014  ,0.29273   ,0.24116  ,0.20222 ,0.1733   ,0.15229  ,0.13747  ,0.12744  ,0.12105  ,0.11733  ,0.11532  ,0.11387  ,0.1113   ,0.10459  ,0.088075 ,0.050963 ,-0.026673,-0.17896 ,-0.45599 ,-0.90898 ,-1.5409  ,-2.242  ,-2.8008  ,-3.027   ,-2.8465  ,-2.2886  ,-1.439  ,-0.41629,
                                                                25.467,22.909,20.354,17.86 ,15.481,13.267,11.252,9.4561 ,7.8848 ,6.5311 ,5.3794   ,4.41       ,3.6008    ,2.9302    ,2.3778   ,1.9252    ,1.5562    ,1.2568   ,1.015    ,0.82069  ,0.66531  ,0.54172   ,0.44397  ,0.36715   ,0.3072   ,0.2608  ,0.22522  ,0.19826  ,0.17813  ,0.16341  ,0.15292  ,0.14573  ,0.14099  ,0.13787  ,0.13528  ,0.13147  ,0.12318  ,0.10408  ,0.062019 ,-0.025296,-0.19507 ,-0.49833 ,-0.97715 ,-1.6063 ,-2.2387  ,-2.6519  ,-2.6836  ,-2.3009  ,-1.5704 ,-0.61188,
                                                                26.312,23.753,21.182,18.657,16.236,13.969,11.894,10.037 ,8.4039 ,6.9911 ,5.7847   ,4.7655     ,3.912     ,3.2023    ,2.6157   ,2.1334    ,1.7385    ,1.4167   ,1.1555   ,0.94418  ,0.77397  ,0.63738   ,0.52822  ,0.44134   ,0.37252  ,0.31827 ,0.27575  ,0.24265  ,0.2171   ,0.1976   ,0.18294  ,0.17216  ,0.16447  ,0.15917  ,0.15554  ,0.15254  ,0.14837  ,0.13953  ,0.11903  ,0.07335  ,-0.022135,-0.20712 ,-0.53135 ,-1.0227 ,-1.6226  ,-2.1509  ,-2.3882  ,-2.2121  ,-1.6378 ,-0.77332,
                                                                27.159,24.604,22.023,19.472,17.011,14.693,12.561,10.642 ,8.9467 ,7.4735 ,6.2105   ,5.1396     ,4.2395    ,3.4885    ,2.8655   ,2.3513    ,1.9287    ,1.5828   ,1.3005   ,1.0709   ,0.88476  ,0.73425   ,0.61292  ,0.51538   ,0.43722  ,0.37478 ,0.32508  ,0.28568  ,0.25461  ,0.23028  ,0.21142  ,0.197    ,0.18624  ,0.17849  ,0.1732   ,0.16975  ,0.16721  ,0.16377  ,0.15574  ,0.13551  ,0.088245 ,-0.012772,-0.20878 ,-0.54529,-1.0304  ,-1.5697  ,-1.9584  ,-1.9972  ,-1.6158 ,-0.8822 ,
                                                                28.007,25.462,22.875,20.303,17.806,15.441,13.252,11.271 ,9.5132 ,7.9784 ,6.6571   ,5.5322     ,4.5833    ,3.7886    ,3.127    ,2.579     ,2.1267    ,1.7549   ,1.4501   ,1.2009   ,0.99775  ,0.83242   ,0.69819  ,0.58943   ,0.50149  ,0.43055 ,0.37344  ,0.3276   ,0.29093  ,0.26173  ,0.23863  ,0.22056  ,0.20667  ,0.19633  ,0.18904  ,0.1844   ,0.18197  ,0.18095  ,0.17964  ,0.17424  ,0.15666  ,0.11091  ,0.0086169,-0.19163,-0.52749 ,-0.98237 ,-1.4265  ,-1.6452  ,-1.4774 ,-0.91444,
                                                                28.854,26.324,23.736,21.147,18.619,16.209,13.966,11.925 ,10.103 ,8.5058 ,7.1244   ,5.9436     ,4.9435    ,4.1029    ,3.4005   ,2.8164    ,2.3326    ,1.9332   ,1.6044   ,1.3344   ,1.1131   ,0.93207   ,0.78424  ,0.66372   ,0.5656   ,0.48584 ,0.42112  ,0.36869  ,0.32632  ,0.2922   ,0.26485  ,0.24312  ,0.22608  ,0.21306  ,0.20358  ,0.19733  ,0.19411  ,0.19369  ,0.19556  ,0.19823  ,0.19798  ,0.18647  ,0.14685  ,0.049714,-0.1447  ,-0.46241 ,-0.85864 ,-1.1745  ,-1.2043 ,-0.84254,
                                                                29.699,27.188,24.605,22.005,19.449,16.997,14.702,12.601 ,10.716 ,9.0553 ,7.6123   ,6.3735     ,5.3202    ,4.4314    ,3.6859   ,3.0637    ,2.5465    ,2.1178   ,1.7636   ,1.4714   ,1.231    ,1.0335    ,0.87138  ,0.73855   ,0.62984  ,0.54097 ,0.46841  ,0.40924  ,0.36108  ,0.32198  ,0.29035  ,0.26493  ,0.24473  ,0.229    ,0.21722  ,0.20911  ,0.20458  ,0.20369  ,0.20657  ,0.21312  ,0.22229  ,0.23074  ,0.23015  ,0.20327 ,0.12049  ,-0.054371,-0.33232 ,-0.64021 ,-0.80223,-0.64312,
                                                                30.539,28.053,25.48 ,22.872,20.293,17.804,15.459,13.3   ,11.352 ,9.6266 ,8.1206   ,6.8221     ,5.7134    ,4.7742    ,3.9836   ,3.3212    ,2.7686    ,2.3091   ,1.9279   ,1.6125   ,1.3519   ,1.137     ,0.95995  ,0.81429   ,0.69458  ,0.59627 ,0.51564  ,0.44956  ,0.39549  ,0.35133  ,0.31538  ,0.28624  ,0.26285  ,0.24438  ,0.23025  ,0.22014  ,0.21399  ,0.21199  ,0.21461  ,0.22249  ,0.23615  ,0.25532  ,0.27739  ,0.29453 ,0.28954  ,0.23343  ,0.095262 ,-0.11926 ,-0.31303,-0.31087,
                                                                31.373,28.917,26.358,23.749,21.151,18.628,16.236,14.02  ,12.01  ,10.219 ,8.6492   ,7.2893     ,6.1233    ,5.1316    ,4.2936   ,3.5892    ,2.9995    ,2.5074   ,2.0979   ,1.7579   ,1.4762   ,1.243     ,1.0504   ,0.89134   ,0.76018  ,0.65211 ,0.56314  ,0.48996  ,0.42984  ,0.38053  ,0.34017  ,0.30728  ,0.28067  ,0.25942  ,0.2429   ,0.23071  ,0.22277  ,0.21929  ,0.22089  ,0.22858  ,0.24374  ,0.26788  ,0.30189  ,0.34438 ,0.38847  ,0.41745  ,0.40329  ,0.32084  ,0.1914  ,0.1264  ,
                                                                32.199,29.777,27.237,24.632,22.021,19.467,17.031,14.76  ,12.688 ,10.833 ,9.1979   ,7.775      ,6.5498    ,5.5038    ,4.6165   ,3.8681    ,3.2394    ,2.7132   ,2.274    ,1.9083   ,1.6043   ,1.3521    ,1.1431   ,0.97012   ,0.82705  ,0.70885 ,0.61126  ,0.53076  ,0.46441  ,0.40981  ,0.36497  ,0.32825  ,0.29837  ,0.27432  ,0.25536  ,0.24104  ,0.2312   ,0.22604  ,0.22619  ,0.23281  ,0.24768  ,0.27328  ,0.31251  ,0.36799 ,0.44009  ,0.52325  ,0.6011   ,0.64587  ,0.63673 ,0.60583 ,
                                                                33.015,30.631,28.116,25.519,22.899,20.32 ,17.842,15.519 ,13.387 ,11.467 ,9.7662   ,8.2793     ,6.9934    ,5.891     ,4.9527   ,4.1584    ,3.489     ,2.9271   ,2.4568   ,2.0641   ,1.7369   ,1.4647    ,1.2387   ,1.0511    ,0.89562  ,0.76687 ,0.66034  ,0.57226  ,0.49949  ,0.43945  ,0.38999  ,0.34935  ,0.31613  ,0.28922  ,0.26779  ,0.25129  ,0.23949  ,0.23254  ,0.23099  ,0.23604  ,0.24961  ,0.2746   ,0.31497  ,0.37563 ,0.46151  ,0.57525  ,0.7125   ,0.85599  ,0.97554 ,1.0504  ,
                                                                33.82 ,31.478,28.992,26.408,23.784,21.183,18.669,16.296 ,14.105 ,12.121 ,10.354   ,8.8021     ,7.4541    ,6.2939    ,5.3026   ,4.4607    ,3.749     ,3.1498   ,2.6469   ,2.226    ,1.8745   ,1.5814    ,1.3375   ,1.1347    ,0.9663   ,0.82657 ,0.71073  ,0.61477  ,0.53534  ,0.46967  ,0.41545  ,0.37078  ,0.33413  ,0.30429  ,0.28033  ,0.26161  ,0.24781  ,0.23898  ,0.23562  ,0.23883  ,0.25053  ,0.27373  ,0.31286  ,0.37396 ,0.46458  ,0.59262  ,0.76305  ,0.97161  ,1.1969  ,1.401   ,
                                                                34.612,32.315,29.862,27.296,24.673,22.055,19.508,17.088 ,14.841 ,12.794 ,10.961   ,9.3434     ,7.9322    ,6.7127    ,5.6669   ,4.7757    ,4.02      ,3.382    ,2.8452   ,2.3948   ,2.0178   ,1.7028    ,1.4402   ,1.2215    ,1.0395   ,0.88832 ,0.76277  ,0.6586   ,0.57224  ,0.50072  ,0.44155  ,0.39271  ,0.35252  ,0.31966  ,0.2931   ,0.27211  ,0.25627  ,0.24552  ,0.24028  ,0.24151  ,0.25106  ,0.27185  ,0.30839  ,0.36719 ,0.45706  ,0.58893  ,0.77396  ,1.0187   ,1.3159  ,1.6341  ,
                                                                35.389,33.14 ,30.725,28.181,25.564,22.934,20.358,17.895 ,15.593 ,13.485 ,11.587   ,9.9034     ,8.4281    ,7.148     ,6.0462   ,5.1041    ,4.3028    ,3.6244   ,3.0522   ,2.5709   ,2.1673   ,1.8295    ,1.5473   ,1.3119    ,1.1158   ,0.95253 ,0.81681  ,0.70405  ,0.61045  ,0.53282  ,0.46851  ,0.41532  ,0.37146  ,0.33548  ,0.30623  ,0.2829   ,0.26497  ,0.25228  ,0.24511  ,0.24432  ,0.25157  ,0.26965  ,0.30292  ,0.3579  ,0.44387  ,0.57321  ,0.76071  ,1.0205   ,1.3581  ,1.7575  ,
                                                                36.149,33.952,31.577,29.061,26.454,23.817,21.217,18.715 ,16.362 ,14.193 ,12.231   ,10.482     ,8.9421    ,7.6004    ,6.4412   ,5.4466    ,4.5982    ,3.8779   ,3.2687   ,2.7554   ,2.3239   ,1.9621    ,1.6594   ,1.4064    ,1.1954   ,1.0196  ,0.87319  ,0.75143  ,0.65023  ,0.56621  ,0.49651  ,0.43878  ,0.39109  ,0.35185  ,0.31983  ,0.29409  ,0.27401  ,0.25934  ,0.25022  ,0.24738  ,0.2523   ,0.26755  ,0.29725  ,0.34767 ,0.42805  ,0.55121  ,0.73377  ,0.99441  ,1.3483  ,1.7948  ,
                                                                36.891,34.748,32.417,29.932,27.341,24.703,22.083,19.546 ,17.144 ,14.919 ,12.894   ,11.079     ,9.4743    ,8.0703    ,6.8525   ,5.8041    ,4.907     ,4.1432   ,3.4957   ,2.9487   ,2.4881   ,2.1012    ,1.7769   ,1.5056    ,1.279    ,1.0899  ,0.93225  ,0.80102  ,0.69184  ,0.6011   ,0.52574  ,0.46326  ,0.41155  ,0.36891  ,0.334    ,0.30575  ,0.28347  ,0.26676  ,0.25568  ,0.2508   ,0.2534   ,0.26582  ,0.29185  ,0.33743 ,0.41141  ,0.5265   ,0.69988  ,0.95265  ,1.3061  ,1.7722  ,
                                                                37.613,35.526,33.243,30.793,28.222,25.587,22.954,20.386 ,17.94  ,15.66  ,13.574   ,11.695     ,10.025    ,8.5581    ,7.2808   ,6.1773    ,5.2299    ,4.4211   ,3.7337   ,3.1518   ,2.6607   ,2.2475    ,1.9006   ,1.6099    ,1.3668   ,1.1638  ,0.99433  ,0.85313  ,0.73554  ,0.63772  ,0.55641  ,0.48891  ,0.43298  ,0.38678  ,0.34883  ,0.31799  ,0.29341  ,0.27463  ,0.26156  ,0.25463  ,0.25494  ,0.26458  ,0.28701  ,0.32772 ,0.39506  ,0.50122  ,0.66323  ,0.90311  ,1.2458  ,1.7117  ,
                                                                38.315,36.285,34.052,31.642,29.096,26.47 ,23.827,21.234 ,18.748 ,16.416 ,14.271   ,12.329     ,10.595    ,9.0645    ,7.7268   ,6.5668    ,5.5679    ,4.7125   ,3.9837   ,3.3653   ,2.8424   ,2.4016    ,2.0309   ,1.7199    ,1.4595   ,1.2417  ,1.0598   ,0.90805  ,0.78159  ,0.67629  ,0.58869  ,0.5159   ,0.45553  ,0.40558  ,0.36444  ,0.33087  ,0.30392  ,0.283    ,0.26791  ,0.25894  ,0.257    ,0.26393  ,0.28286  ,0.31883 ,0.37958  ,0.47661  ,0.62634  ,0.85076  ,1.1764  ,1.6294  ,
                                                                38.995,37.024,34.844,32.476,29.96 ,27.348,24.702,22.088 ,19.566 ,17.187 ,14.985   ,12.981     ,11.183    ,9.5896    ,8.1908   ,6.9735    ,5.9216    ,5.0182   ,4.2464   ,3.59     ,3.0338   ,2.5641    ,2.1685   ,1.8361    ,1.5574   ,1.324   ,1.129    ,0.96609  ,0.83024  ,0.71703  ,0.62278  ,0.54441  ,0.47933  ,0.42541  ,0.38093  ,0.34449  ,0.31506  ,0.29192  ,0.27478  ,0.26375  ,0.25959  ,0.2639   ,0.27949  ,0.31091 ,0.3653   ,0.45338  ,0.59069  ,0.79858  ,1.1041  ,1.5364  ,
                                                                39.652,37.741,35.616,33.294,30.812,28.219,25.575,22.945 ,20.393 ,17.97  ,15.715   ,13.651     ,11.791    ,10.134    ,8.6736   ,7.3979    ,6.2918    ,5.3389   ,4.5227   ,3.8267   ,3.2357   ,2.7357    ,2.3139   ,1.959     ,1.661    ,1.4112  ,1.2022   ,1.0276   ,0.88176  ,0.76018  ,0.65888  ,0.57458  ,0.50452  ,0.44641  ,0.39838  ,0.35893  ,0.32689  ,0.30146  ,0.2822   ,0.26911  ,0.26276  ,0.26454  ,0.27691  ,0.30403 ,0.35238  ,0.43188  ,0.55707  ,0.74832  ,1.0322  ,1.4396  ,
                                                                40.286,38.436,36.367,34.095,31.651,29.081,26.445,23.806 ,21.228 ,18.765 ,16.46    ,14.339     ,12.417    ,10.698    ,9.1755   ,7.8407    ,6.6792    ,5.6755   ,4.8131   ,4.0761   ,3.4489   ,2.9171    ,2.4678   ,2.0891    ,1.7708   ,1.5037  ,1.2799   ,1.0928   ,0.93643  ,0.80595  ,0.69718  ,0.60659  ,0.53124  ,0.46869  ,0.41691  ,0.37427  ,0.3395   ,0.31167  ,0.29023  ,0.27504  ,0.26653  ,0.26583  ,0.27515  ,0.29819 ,0.34087  ,0.41227  ,0.52588  ,0.70091  ,0.96311 ,1.3436};

        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 11:
    {
        static const std::vector<double> arrayCurvatureMin   = {23.102,21.94 ,20.879,19.894,18.943,17.95 ,16.786,15.261,13.178,10.531,7.6944,5.2435,3.506 ,2.4377,1.835 ,1.5114,1.344 ,1.2612,1.223 ,1.2067,1.199 ,1.1916,1.179 ,1.1577,1.1252,1.0801,1.0212,0.94804,0.8603,0.75783,0.64056,0.50843,0.36133,0.19905,0.021263,-0.17245,-0.38265,-0.60996,-0.85503,-1.1185 ,-1.401  ,-1.7028 ,-2.024  ,-2.3644 ,-2.7231 ,-3.0987  ,-3.4894 ,-3.8926 ,-4.3052 ,-4.7237 ,
                                                                24.385,23.141,21.996,20.926,19.886,18.795,17.519,15.862,13.637,10.873,7.9831,5.5375,3.8237,2.77  ,2.1696,1.8425,1.672 ,1.5904,1.5593,1.5563,1.5674,1.5834,1.5977,1.6054,1.6028,1.5874,1.5572,1.5108 ,1.4472,1.366  ,1.2667 ,1.1489 ,1.0126 ,0.85729,0.68285 ,0.48885 ,0.27489 ,0.040538,-0.2146 ,-0.49082,-0.78822,-1.1066 ,-1.4455 ,-1.8038 ,-2.1799 ,-2.5717  ,-2.9764 ,-3.3907 ,-3.8109 ,-4.2332 ,
                                                                25.74 ,24.417,23.193,22.041,20.914,19.726,18.337,16.547,14.179,11.3  ,8.358 ,5.9138,4.2133,3.1605,2.5475,2.2015,2.0125,1.9169,1.8786,1.8758,1.8946,1.9248,1.9588,1.9906,2.0152,2.0286,2.0276,2.0098 ,1.9732,1.9166 ,1.8388 ,1.7393 ,1.6174 ,1.4727 ,1.3049  ,1.1136  ,0.89868 ,0.65984 ,0.39708 ,0.11059 ,-0.19916,-0.53131,-0.88451,-1.2569 ,-1.646  ,-2.0488  ,-2.4618 ,-2.8811 ,-3.3028 ,-3.7225 ,
                                                                27.157,25.758,24.458,23.227,22.016,20.734,19.233,17.309,14.801,11.812,8.8231,6.3797,4.686 ,3.6243,2.9872,2.6096,2.3882,2.264 ,2.2034,2.186 ,2.1982,2.2299,2.2727,2.3196,2.3644,2.4018,2.4274,2.4371 ,2.4281,2.3977 ,2.3441 ,2.2659 ,2.1619 ,2.0314 ,1.8738  ,1.6888  ,1.4763  ,1.2364  ,0.96954 ,0.67641 ,0.35824 ,0.016674,-0.34605,-0.72716,-1.1233 ,-1.5308  ,-1.9455 ,-2.3631 ,-2.7793 ,-3.1903 ,
                                                                28.624,27.152,25.779,24.473,23.18 ,21.806,20.194,18.139,15.495,12.403,9.376 ,6.9355,5.2454,4.169 ,3.5004,3.0822,2.8179,2.6529,2.5564,2.5098,2.5004,2.5186,2.556 ,2.6053,2.6594,2.712 ,2.7573,2.79   ,2.8056,2.8006 ,2.7716 ,2.7162 ,2.6325 ,2.5193 ,2.3754  ,2.2005  ,1.9946  ,1.758   ,1.4916  ,1.1969  ,0.87565 ,0.53056 ,0.1647  ,-0.21828,-0.61429,-1.0189  ,-1.4277 ,-1.8362 ,-2.2403 ,-2.6361 ,
                                                                30.133,28.59 ,27.146,25.767,24.395,22.93 ,21.21 ,19.026,16.25 ,13.066,10.009,7.5757,5.8888,4.7952,4.0916,3.6279,3.3141,3.0997,2.9567,2.8684,2.8235,2.8134,2.8303,2.867 ,2.9163,2.9712,3.0252,3.0719 ,3.1056,3.1211 ,3.114  ,3.0805 ,3.0178 ,2.9235 ,2.7963  ,2.6352  ,2.4403  ,2.2121  ,1.9519  ,1.6618  ,1.3445  ,1.0034  ,0.64235 ,0.2658  ,-0.12158,-0.51502 ,-0.90986,-1.3017 ,-1.6868 ,-2.0617 ,
                                                                31.676,30.063,28.549,27.099,25.65 ,24.096,22.268,19.957,17.055,13.789,10.712,8.2908,6.6083,5.4974,4.7582,4.2477,3.8816,3.6133,3.4171,3.2782,3.1869,3.1357,3.1177,3.1265,3.1555,3.1976,3.2461,3.2941 ,3.3348,3.3621 ,3.3702 ,3.3542 ,3.3101 ,3.2345 ,3.1251  ,2.9808  ,2.8011  ,2.5868  ,2.3394  ,2.0615  ,1.7563  ,1.428   ,1.0811  ,0.72048 ,0.3511  ,-0.022119,-0.39463,-0.76236,-1.1219 ,-1.4705 ,
                                                                33.246,31.564,29.982,28.461,26.936,25.294,23.358,20.921,17.9  ,14.562,11.475,9.0695,7.3933,6.2664,5.4929,4.9367,4.5186,4.1951,3.943 ,3.7487,3.6038,3.502 ,3.4376,3.4052,3.399 ,3.4127,3.4399,3.4736 ,3.5071,3.5334 ,3.5461 ,3.5393 ,3.5078 ,3.4474 ,3.355   ,3.2286  ,3.0675  ,2.8722  ,2.6444  ,2.3867  ,2.1029  ,1.7973  ,1.4748  ,1.1403  ,0.79879 ,0.45505  ,0.11325 ,-0.22301,-0.55084,-0.86807,
                                                                34.837,33.087,31.436,29.845,28.244,26.513,24.47 ,21.909,18.774,15.374,12.284,9.8999,8.2318,7.0908,6.2853,5.6859,5.2179,4.8408,4.5328,4.2817,4.0799,3.922 ,3.8035,3.7197,3.6661,3.6374,3.6279,3.6316 ,3.6417,3.6519 ,3.6552 ,3.6456 ,3.6171 ,3.565  ,3.4854  ,3.3757  ,3.2346  ,3.0623  ,2.86    ,2.6305  ,2.3772  ,2.1044  ,1.8166  ,1.5186  ,1.215   ,0.90985  ,0.6069  ,0.30912 ,0.018818,-0.26234,
                                                                36.445,34.626,32.907,31.246,29.568,27.748,25.596,22.912,19.669,16.216,13.131,10.771,9.1119,7.9584,7.1235,6.484 ,5.9694,5.5415,5.1798,4.873 ,4.614 ,4.3981,4.2213,4.0801,3.9706,3.8887,3.8298,3.7889 ,3.7604,3.7384 ,3.717  ,3.6901 ,3.6521 ,3.5981 ,3.5238  ,3.4264  ,3.304   ,3.1563  ,2.9839  ,2.789   ,2.5743  ,2.3434  ,2.1001  ,1.8483  ,1.5917  ,1.3338   ,1.0773  ,0.82457 ,0.57733 ,0.33675 ,
                                                                38.066,36.178,34.391,32.658,30.903,28.991,26.729,23.923,20.577,17.079,14.006,11.671,10.023,8.8578,7.9957,7.319 ,6.761 ,6.2857,5.8736,5.5136,5.1993,4.926 ,4.6902,4.4889,4.3189,4.177 ,4.0597,3.9632 ,3.8831,3.8149 ,3.7538 ,3.695  ,3.6336 ,3.5653 ,3.4864  ,3.3939  ,3.2857  ,3.161   ,3.0199  ,2.8632  ,2.6926  ,2.5106  ,2.3194  ,2.1218  ,1.9203  ,1.7169   ,1.5137  ,1.312   ,1.113   ,0.91738 ,
                                                                39.699,37.74 ,35.883,34.078,32.243,30.238,27.864,24.935,21.492,17.957,14.9  ,12.593,10.954,9.7784,8.8907,8.1792,7.5809,7.0614,6.6021,6.1923,5.8254,5.497 ,5.2035,4.9422,4.7104,4.5055,4.3249,4.1658 ,4.0253,3.9003 ,3.7877 ,3.6841 ,3.5863 ,3.4913 ,3.3964  ,3.2993  ,3.1981  ,3.0915  ,2.979   ,2.8604  ,2.7358  ,2.6059  ,2.4714  ,2.3332  ,2.1921  ,2.0489   ,1.9043  ,1.7587  ,1.6126  ,1.4662  ,
                                                                41.339,39.31 ,37.381,35.503,33.586,31.485,28.996,25.946,22.409,18.845,15.808,13.527,11.897,10.711,9.7986,9.054 ,8.4176,7.8567,7.3532,6.8966,6.4802,6.0995,5.7508,5.4315,5.139 ,4.8711,4.6258,4.4014 ,4.1961,4.0083 ,3.8365 ,3.6792 ,3.5352 ,3.4031 ,3.2818  ,3.17    ,3.0666  ,2.9704  ,2.8803  ,2.7951  ,2.7136  ,2.6345  ,2.5567  ,2.4791  ,2.4007  ,2.3205   ,2.2378  ,2.152   ,2.0626  ,1.9693  ,
                                                                42.988,40.886,38.884,36.93 ,34.929,32.729,30.123,26.95 ,23.325,19.738,16.724,14.468,12.846,11.648,10.711,9.9339,9.2608,8.6604,8.115 ,7.6142,7.1511,6.7209,6.32  ,5.9454,5.5947,5.2659,4.9576,4.6685 ,4.3982,4.1462 ,3.9125 ,3.6976 ,3.502  ,3.3261 ,3.1706  ,3.0355  ,2.9208  ,2.8257  ,2.7489  ,2.6887  ,2.6427  ,2.6081  ,2.582   ,2.5615  ,2.5435  ,2.5256   ,2.5052  ,2.4806  ,2.4502  ,2.4129  ,
                                                                44.642,42.466,40.39 ,38.358,36.27 ,33.967,31.24 ,27.946,24.238,20.633,17.642,15.411,13.794,12.582,11.62 ,10.811,10.102,9.4629,8.877 ,8.3337,7.8261,7.349 ,6.8985,6.4717,6.066 ,5.6798,5.3118,4.9616 ,4.6293,4.3156 ,4.0221 ,3.7504 ,3.5026 ,3.2809 ,3.0874  ,2.9234  ,2.7899  ,2.6869  ,2.6133  ,2.5671  ,2.5454  ,2.5446  ,2.5604  ,2.5884  ,2.6241  ,2.6632   ,2.7019  ,2.7368  ,2.7652  ,2.7851  ,
                                                                46.302,44.05 ,41.898,39.785,37.607,35.198,32.348,28.931,25.145,21.527,18.561,16.351,14.735,13.508,12.519,11.678,10.933,10.256,9.6299,9.0453,8.4945,7.9723,7.4747,6.9986,6.5415,6.1019,5.6789,5.2727 ,4.884 ,4.5145 ,4.1665 ,3.8429 ,3.5469 ,3.2818 ,3.0506  ,2.8557  ,2.699   ,2.5809  ,2.501   ,2.4575  ,2.4474  ,2.4669  ,2.5111  ,2.5749  ,2.6529  ,2.7394   ,2.8295  ,2.9184  ,3.0022  ,3.0777  ,
                                                                47.967,45.638,43.407,41.211,38.94 ,36.421,33.443,29.906,26.046,22.419,19.477,17.285,15.667,14.422,13.404,12.529,11.747,11.031,10.366,9.7401,9.1469,8.581 ,8.0383,7.5156,7.0107,6.5222,6.0497,5.5938 ,5.1559,4.7382 ,4.3437 ,3.9761 ,3.6389 ,3.3361 ,3.0712  ,2.8471  ,2.6659  ,2.5286  ,2.4352  ,2.3844  ,2.3737  ,2.3994  ,2.4572  ,2.5416  ,2.647   ,2.7675   ,2.8972  ,3.0304  ,3.1622  ,3.2881  ,
                                                                49.637,47.23 ,44.918,42.636,40.268,37.634,34.526,30.868,26.939,23.307,20.387,18.209,16.586,15.319,14.271,13.36 ,12.539,11.784,11.078,10.411,9.7758,9.1672,8.581 ,8.0144,7.4652,6.9326,6.4166,5.9181 ,5.439 ,4.982  ,4.5504 ,4.1479 ,3.7785 ,3.4459 ,3.1536  ,2.9044  ,2.7004  ,2.5428  ,2.4318  ,2.3666  ,2.3453  ,2.365   ,2.4219  ,2.5115  ,2.6284  ,2.767    ,2.9214  ,3.0857  ,3.2541  ,3.4213  ,
                                                                51.313,48.825,46.431,44.06 ,41.591,38.838,35.596,31.819,27.825,24.189,21.29 ,19.122,17.489,16.197,15.115,14.166,13.306,12.509,11.761,11.052,10.375,9.7245,9.0965,8.4886,7.8991,7.3274,6.7741,6.2404 ,5.7285,5.2414 ,4.7824 ,4.355  ,3.9627 ,3.6089 ,3.2964  ,3.0278  ,2.8047  ,2.6281  ,2.4984  ,2.4149  ,2.3765  ,2.381   ,2.4256  ,2.5067  ,2.6201  ,2.7609   ,2.924   ,3.1036  ,3.2941  ,3.4898  ,
                                                                52.995,50.425,47.945,45.482,42.908,40.031,36.654,32.757,28.703,25.066,22.184,20.021,18.374,17.053,15.935,14.945,14.042,13.202,12.412,11.66 ,10.941,10.249,9.5809,8.9346,8.3087,7.7033,7.1191,6.5577 ,6.0216,5.5134 ,5.0363 ,4.5931 ,4.1868 ,3.82   ,3.4948  ,3.2129  ,2.9753  ,2.7827  ,2.6353  ,2.5327  ,2.474   ,2.4577  ,2.4821  ,2.5447  ,2.6424  ,2.7716   ,2.9281  ,3.1074  ,3.3042  ,3.5134  ,
                                                                54.683,52.029,49.462,46.903,44.219,41.215,37.698,33.684,29.573,25.936,23.069,20.905,19.239,17.885,16.727,15.693,14.746,13.861,13.026,12.231,11.47 ,10.738,10.032,9.3508,8.6931,8.0595,7.451 ,6.8694 ,6.3169,5.7959 ,5.3089 ,4.8582 ,4.4457 ,4.0731 ,3.7416  ,3.452   ,3.2047  ,2.9998  ,2.8373  ,2.7167  ,2.6373  ,2.5984  ,2.5986  ,2.6365  ,2.7103  ,2.8175   ,2.9553  ,3.1204  ,3.3087  ,3.5158  ,
                                                                56.378,53.638,50.98 ,48.323,45.526,42.388,38.729,34.6  ,30.436,26.799,23.942,21.773,20.082,18.692,17.489,16.409,15.415,14.484,13.603,12.765,11.962,11.192,10.451,9.738 ,9.0533,8.3972,7.7709,7.176  ,6.6144,6.0879 ,5.5981 ,5.1465 ,4.734  ,4.3613 ,4.0287  ,3.7361  ,3.4833  ,3.2699  ,3.0955  ,2.9594  ,2.861   ,2.7999  ,2.7752  ,2.7861  ,2.8318  ,2.9109   ,3.0217  ,3.1622  ,3.3297  ,3.5212  ,
                                                                58.081,55.253,52.502,49.742,46.827,43.552,39.748,35.506,31.293,27.655,24.804,22.623,20.903,19.471,18.222,17.093,16.048,15.069,14.143,13.261,12.418,11.611,10.839,10.099,9.3923,8.7195,8.0815,7.4796 ,6.915 ,6.3888 ,5.9018 ,5.4543 ,5.0464 ,4.6778 ,4.3479  ,4.056   ,3.8012  ,3.5826  ,3.3993  ,3.2507  ,3.1361  ,3.0549  ,3.0068  ,2.9914  ,3.0082  ,3.0568   ,3.1365  ,3.2463  ,3.3849  ,3.5502  ,
                                                                59.791,56.873,54.027,51.162,48.123,44.706,40.756,36.402,32.143,28.503,25.652,23.455,21.7  ,20.223,18.923,17.742,16.647,15.619,14.646,13.722,12.841,12.001,11.2  ,10.438,9.7144,9.0303,8.3861,7.7825 ,7.2199,6.6985 ,6.2182 ,5.7784 ,5.3781 ,5.0163 ,4.6917  ,4.4029  ,4.1485  ,3.9273  ,3.7381  ,3.5799  ,3.452   ,3.3539  ,3.2852  ,3.2458  ,3.2355  ,3.2544   ,3.3024  ,3.3792  ,3.4846  ,3.6177  ,
                                                                61.508,58.5  ,55.556,52.581,49.414,45.85 ,41.752,37.289,32.986,29.343,26.488,24.268,22.473,20.947,19.594,18.359,17.211,16.134,15.116,14.151,13.234,12.364,11.54 ,10.76 ,10.025,9.3342,8.6884,8.0872 ,7.5303,7.0169 ,6.5457 ,6.1155 ,5.7246 ,5.371  ,5.0531  ,4.7688  ,4.5165  ,4.2946  ,4.1017  ,3.9366  ,3.7984  ,3.6866  ,3.6007  ,3.5406  ,3.5063  ,3.4981   ,3.5162  ,3.561   ,3.6327  ,3.7312  ,
                                                                63.234,60.133,57.088,53.999,50.7  ,46.984,42.736,38.168,33.823,30.174,27.31 ,25.062,23.223,21.644,20.235,18.945,17.744,16.618,15.556,14.552,13.604,12.708,11.864,11.071,10.329,9.6356,8.9919,8.3962 ,7.8472,7.3434 ,6.8826 ,6.4628 ,6.0817 ,5.7369 ,5.4261  ,5.147   ,4.8977  ,4.6764  ,4.4814  ,4.3114  ,4.1655  ,4.0429  ,3.9431  ,3.8659  ,3.8114  ,3.7798   ,3.7716  ,3.7871  ,3.827   ,3.8919  ,
                                                                64.967,61.772,58.624,55.417,51.981,48.109,43.71 ,39.038,34.655,30.997,28.118,25.836,23.95 ,22.315,20.849,19.502,18.25 ,17.076,15.972,14.933,13.957,13.039,12.18 ,11.378,10.631,9.9389,9.2995,8.711  ,8.1712,7.6774 ,7.2272 ,6.8176 ,6.446  ,6.1095 ,5.8056  ,5.5317  ,5.2857  ,5.0656  ,4.8696  ,4.6964  ,4.5448  ,4.414   ,4.3033  ,4.2125  ,4.1415  ,4.0905   ,4.0598  ,4.0501  ,4.0619  ,4.0961  ,
                                                                66.708,63.417,60.162,56.834,53.256,49.224,44.674,39.902,35.481,31.812,28.913,26.593,24.655,22.963,21.438,20.036,18.732,17.513,16.371,15.301,14.299,13.364,12.493,11.685,10.937,10.247,9.6137,9.033  ,8.5022,8.0182 ,7.5777 ,7.1774 ,6.8143 ,6.4852 ,6.1874  ,5.9183  ,5.6754  ,5.4567  ,5.2605  ,5.0851  ,4.9294  ,4.7924  ,4.6733  ,4.5719  ,4.4878  ,4.4211   ,4.3721  ,4.3413  ,4.3293  ,4.3368  ,
                                                                68.456,65.066,61.703,58.25 ,54.524,50.329,45.628,40.76 ,36.302,32.62 ,29.696,27.334,25.341,23.59 ,22.007,20.55 ,19.198,17.936,16.759,15.661,14.638,13.689,12.809,11.997,11.25 ,10.564,9.9361,9.3627 ,8.8401,8.3646 ,7.9324 ,7.54   ,7.184  ,6.8611 ,6.5683  ,6.303   ,6.0628  ,5.8455  ,5.6492  ,5.4723  ,5.3136  ,5.1719  ,5.0465  ,4.9369  ,4.8426  ,4.7636   ,4.7     ,4.6521  ,4.6205  ,4.6058  ,
                                                                70.209,66.72 ,63.246,59.663,55.786,51.423,46.572,41.612,37.12 ,33.421,30.469,28.061,26.013,24.203,22.562,21.053,19.654,18.354,17.145,16.022,14.981,14.02 ,13.134,12.319,11.573,10.89 ,10.267,9.7002 ,9.1843,8.7156 ,8.2899 ,7.9036 ,7.5529 ,7.2346 ,6.9456  ,6.6832  ,6.4449  ,6.2285  ,6.0321  ,5.8541  ,5.6931  ,5.5479  ,5.4177  ,5.3017  ,5.1996  ,5.1111   ,5.0361  ,4.9747  ,4.9275  ,4.8948  ,
                                                                71.967,68.378,64.789,61.073,57.04 ,52.508,47.509,42.461,37.936,34.218,31.234,28.778,26.674,24.806,23.11 ,21.551,20.109,18.772,17.535,16.39 ,15.334,14.362,13.47 ,12.654,11.908,11.227,10.608,10.045 ,9.5341,9.07   ,8.6488 ,8.2665 ,7.9194 ,7.604  ,7.3174  ,7.0568  ,6.8195  ,6.6035  ,6.4068  ,6.2277  ,6.0647  ,5.9167  ,5.7828  ,5.662   ,5.5539  ,5.458    ,5.3742  ,5.3024  ,5.2429  ,5.196   ,
                                                                73.727,70.036,66.331,62.479,58.286,53.582,48.438,43.308,38.751,35.013,31.995,29.491,27.33 ,25.406,23.657,22.051,20.569,19.2  ,17.936,16.772,15.702,14.72 ,13.823,13.003,12.256,11.576,10.958,10.398 ,9.8887,9.427  ,9.0079 ,8.6275 ,8.2821 ,7.968  ,7.6823  ,7.4222  ,7.1851  ,6.9687  ,6.7713  ,6.5909  ,6.4262  ,6.2759  ,6.1389  ,6.0143  ,5.9015  ,5.8      ,5.7095  ,5.6297  ,5.5607  ,5.5027  ,
                                                                75.489,71.695,67.872,63.88 ,59.523,54.648,49.363,44.155,39.571,35.811,32.758,30.205,27.99 ,26.011,24.212,22.562,21.042,19.643,18.355,17.172,16.088,15.097,14.192,13.368,12.618,11.936,11.318,10.757 ,10.247,9.7855 ,9.3664 ,8.9858 ,8.6401 ,8.3257 ,8.0395  ,7.7786  ,7.5406  ,7.3232  ,7.1245  ,6.9426  ,6.7761  ,6.6236  ,6.484   ,6.3563  ,6.2399  ,6.134    ,6.0383  ,5.9524  ,5.8762  ,5.8097  ,
                                                                77.25 ,73.353,69.409,65.274,60.752,55.706,50.285,45.008,40.397,36.616,33.527,30.927,28.659,26.629,24.782,23.091,21.536,20.108,18.797,17.595,16.497,15.494,14.58 ,13.749,12.994,12.308,11.686,11.121 ,10.609,10.145 ,9.7235 ,9.3408 ,8.9929 ,8.6764 ,8.3883  ,8.1255  ,7.8856  ,7.6664  ,7.4658  ,7.282   ,7.1134  ,6.9587  ,6.8168  ,6.6865  ,6.567   ,6.4577   ,6.3579  ,6.2673  ,6.1856  ,6.1127  ,
                                                                79.009,75.007,70.942,66.663,61.974,56.758,51.209,45.869,41.236,37.435,34.311,31.665,29.346,27.266,25.375,23.644,22.056,20.599,19.265,18.044,16.929,15.913,14.988,14.147,13.384,12.691,12.062,11.492 ,10.974,10.505 ,10.079 ,9.692  ,9.3402 ,9.0201 ,8.7286  ,8.4627  ,8.22    ,7.9981  ,7.7949  ,7.6087  ,7.4378  ,7.2808  ,7.1365  ,7.0038  ,6.8817  ,6.7695   ,6.6665  ,6.5722  ,6.4863  ,6.4085  ,
                                                                80.762,76.658,72.47 ,68.046,63.19 ,57.809,52.139,46.745,42.094,38.273,35.116,32.425,30.059,27.932,25.997,24.227,22.606,21.121,19.761,18.52 ,17.387,16.355,15.416,14.562,13.787,13.084,12.446,11.867 ,11.342,10.865 ,10.432 ,10.039 ,9.682  ,9.3568 ,9.0606  ,8.7905  ,8.5438  ,8.3184  ,8.112   ,7.9229  ,7.7493  ,7.5897  ,7.4429  ,7.3078  ,7.1833  ,7.0685   ,6.9629  ,6.8657  ,6.7765  ,6.6949  ,
                                                                82.51 ,78.303,73.994,69.424,64.403,58.861,53.081,47.641,42.977,39.139,35.95 ,33.216,30.803,28.631,26.654,24.846,23.191,21.675,20.289,19.024,17.87 ,16.819,15.863,14.994,14.204,13.488,12.837,12.247 ,11.712,11.226 ,10.784 ,10.383 ,10.018 ,9.6867 ,9.3845  ,9.109   ,8.8575  ,8.6277  ,8.4175  ,8.2248  ,8.048   ,7.8855  ,7.7361  ,7.5985  ,7.4716  ,7.3545   ,7.2465  ,7.1468  ,7.0549  ,6.9704  ,
                                                                84.25 ,79.943,75.514,70.8  ,65.615,59.921,54.042,48.565,43.892,40.039,36.819,34.044,31.586,29.369,27.349,25.503,23.812,22.265,20.85 ,19.558,18.379,17.306,16.329,15.441,14.634,13.901,13.236,12.632 ,12.084,11.586 ,11.134 ,10.723 ,10.35  ,10.01  ,9.7008  ,9.4189  ,9.1616  ,8.9266  ,8.7117  ,8.5149  ,8.3344  ,8.1687  ,8.0163  ,7.876   ,7.7467  ,7.6273   ,7.5171  ,7.4152  ,7.3211  ,7.2342  ,
                                                                85.982,81.578,77.031,72.176,66.833,60.994,55.027,49.524,44.847,40.98 ,37.73 ,34.915,32.412,30.15 ,28.087,26.2  ,24.472,22.89 ,21.442,20.12 ,18.914,17.815,16.814,15.904,15.077,14.325,13.642,13.022 ,12.459,11.947 ,11.483 ,11.061 ,10.677 ,10.328 ,10.01   ,9.7207  ,9.4567  ,9.2157  ,8.9955  ,8.7939  ,8.6093  ,8.4398  ,8.2842  ,8.1409  ,8.009   ,7.8872   ,7.7748  ,7.6708  ,7.5746  ,7.4857  ,
                                                                87.706,83.209,78.549,73.556,68.06 ,62.087,56.045,50.524,45.847,41.967,38.688,35.833,33.285,30.977,28.87 ,26.94 ,25.171,23.551,22.068,20.713,19.475,18.347,17.319,16.383,15.532,14.759,14.055,13.417 ,12.837,12.31  ,11.831 ,11.396 ,11     ,10.641 ,10.313  ,10.015  ,9.7437  ,9.4959  ,9.2695  ,9.0626  ,8.8732  ,8.6996  ,8.5402  ,8.3938  ,8.2589  ,8.1346   ,8.0199  ,7.9138  ,7.8157  ,7.7249  ,
                                                                89.424,84.839,80.071,74.945,69.303,63.206,57.102,51.572,46.897,43.006,39.698,36.802,34.208,31.852,29.698,27.723,25.91 ,24.248,22.726,21.334,20.062,18.9  ,17.842,16.878,16    ,15.202,14.477,13.817 ,13.218,12.673 ,12.179 ,11.729 ,11.321 ,10.949 ,10.611  ,10.304  ,10.023  ,9.7678  ,9.5347  ,9.3218  ,9.1271  ,8.9488  ,8.7853  ,8.6352  ,8.4972  ,8.3701   ,8.253   ,8.1448  ,8.0447  ,7.9521  ,
                                                                91.137,86.472,81.601,76.348,70.567,64.358,58.204,52.673,48.004,44.101,40.762,37.824,35.182,32.777,30.572,28.548,26.689,24.982,23.417,21.984,20.673,19.476,18.384,17.388,16.482,15.656,14.906,14.223 ,13.603,13.039 ,12.527 ,12.062 ,11.639 ,11.254 ,10.904  ,10.586  ,10.296  ,10.033  ,9.7919  ,9.5723  ,9.3717  ,9.1882  ,9.0202  ,8.8661  ,8.7246  ,8.5945   ,8.4746  ,8.3641  ,8.262   ,8.1675  ,
                                                                92.849,88.11 ,83.143,77.769,71.858,65.549,59.357,53.831,49.168,45.253,41.883,38.901,36.209,33.75 ,31.493,29.417,27.507,25.751,24.14 ,22.662,21.31 ,20.074,18.945,17.915,16.976,16.122,15.344,14.636 ,13.993,13.408 ,12.877 ,12.394 ,11.955 ,11.556 ,11.194  ,10.864  ,10.564  ,10.291  ,10.042  ,9.8151  ,9.608   ,9.4188  ,9.2457  ,9.0872  ,8.9419  ,8.8084   ,8.6856  ,8.5725  ,8.4682  ,8.3718  ,
                                                                94.563,89.758,84.702,79.214,73.181,66.784,60.565,55.049,50.394,46.464,43.06 ,40.032,37.287,34.773,32.46 ,30.328,28.364,26.556,24.894,23.369,21.972,20.693,19.525,18.458,17.485,16.598,15.791,15.056 ,14.388,13.781 ,13.229 ,12.727 ,12.271 ,11.857 ,11.481  ,11.138  ,10.827  ,10.544  ,10.286  ,10.051  ,9.8369  ,9.6414  ,9.4628  ,9.2994  ,9.1498  ,9.0126   ,8.8866  ,8.7707  ,8.664   ,8.5655  ,
                                                                96.284,91.421,86.281,80.686,74.539,68.067,61.83 ,56.329,51.681,47.735,44.295,41.217,38.417,35.844,33.472,31.281,29.259,27.396,25.681,24.105,22.659,21.335,20.124,19.018,18.008,17.087,16.248,15.485 ,14.79 ,14.158 ,13.584 ,13.062 ,12.588 ,12.157 ,11.766  ,11.41   ,11.086  ,10.792  ,10.525  ,10.281  ,10.059  ,9.8569  ,9.6723  ,9.5036  ,9.3493  ,9.208    ,9.0785  ,8.9595  ,8.8501  ,8.7493  ,
                                                                98.014,93.101,87.885,82.188,75.937,69.399,63.155,57.671,53.029,49.064,45.585,42.456,39.598,36.963,34.528,32.275,30.193,28.27 ,26.499,24.869,23.372,22    ,20.744,19.595,18.546,17.589,16.716,15.922 ,15.199,14.541 ,13.943 ,13.4   ,12.906 ,12.457 ,12.05   ,11.679  ,11.343  ,11.037  ,10.759  ,10.506  ,10.276  ,10.066  ,9.8751  ,9.7006  ,9.5412  ,9.3955   ,9.262   ,9.1397  ,9.0273  ,8.924   ,
                                                                99.757,94.802,89.516,83.722,77.374,70.782,64.538,59.074,54.436,50.451,46.93 ,43.748,40.828,38.129,35.629,33.311,31.164,29.18 ,27.349,25.662,24.111,22.688,21.384,20.191,19.1  ,18.104,17.196,16.369 ,15.616,14.93  ,14.307 ,13.741 ,13.226 ,12.759 ,12.334  ,11.948  ,11.598  ,11.28   ,10.99   ,10.727  ,10.488  ,10.27   ,10.072  ,9.8913  ,9.7264  ,9.5758   ,9.438   ,9.3119  ,9.1963  ,9.0902  ,
                                                                101.52,96.525,91.174,85.289,78.853,72.217,65.979,60.538,55.902,51.893,48.329,45.09 ,42.107,39.341,36.773,34.387,32.174,30.124,28.231,26.484,24.877,23.4  ,22.046,20.805,19.671,18.634,17.689,16.827 ,16.042,15.327 ,14.677 ,14.087 ,13.55  ,13.062 ,12.619  ,12.217  ,11.852  ,11.52   ,11.219  ,10.945  ,10.696  ,10.47   ,10.264  ,10.076  ,9.9056  ,9.7497   ,9.6073  ,9.4771  ,9.358   ,9.2487  ,
                                                                103.29,98.271,92.859,86.888,80.371,73.701,67.478,62.058,57.424,53.389,49.779,46.481,43.433,40.598,37.959,35.503,33.22 ,31.104,29.145,27.337,25.67 ,24.137,22.73 ,21.44 ,20.259,19.18 ,18.195,17.296 ,16.478,15.733 ,15.054 ,14.438 ,13.878 ,13.369 ,12.907  ,12.487  ,12.106  ,11.76   ,11.446  ,11.161  ,10.902  ,10.666  ,10.452  ,10.257  ,10.08   ,9.918    ,9.7706  ,9.636   ,9.5129  ,9.4003  ,
                                                                105.09,100.04,94.57 ,88.517,81.928,75.233,69.03 ,63.634,58.999,54.936,51.278,47.92 ,44.805,41.899,39.188,36.659,34.305,32.118,30.092,28.219,26.491,24.899,23.437,22.095,20.867,19.743,18.716,17.779 ,16.925,16.147 ,15.44  ,14.796 ,14.211 ,13.68  ,13.197  ,12.759  ,12.361  ,12      ,11.672  ,11.375  ,11.105  ,10.859  ,10.636  ,10.434  ,10.249  ,10.082   ,9.9286  ,9.7892  ,9.6619  ,9.5456};

        return arrayCurvatureMin[y * arrayWidth + x];
    }
    }

    return 0;
}
