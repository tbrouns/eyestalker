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

#include "headers/eyestalker.h"

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

double ceil2( double value )
{
    if (value < 0.0) { return floor(value); }
    else             { return  ceil(value); }
}

inline double calculateScoreIntensity(double d_intensity)
{
    return (1 - 1 / (1 + exp(-0.10 * (d_intensity - 25))));
}

inline double calculateScoreRadius(double d_radius)
{
    return (exp(-0.05 * pow(d_radius, 2)));
}

inline double calculateScoreCurvature(double d_curvature)
{
    return (exp(-0.05 * pow(d_curvature, 2)));
}

inline double calculateScoreCircumference(double d_circumference)
{
    return (1 - 1 / (1 + exp(-0.10 * (d_circumference - 75))));
}

inline double calculateScoreGradient(double d_gradient)
{
    return d_gradient;
}

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

AOIProperties detectPupilApprox(const std::vector<unsigned int>& I, AOIProperties searchAOI, AOIProperties innerAOI, AOIProperties glintAOI)
{
    int stepSize = 2;

    innerAOI.xPos = 0;
    innerAOI.yPos = 0;
    
    double minPupilIntensity = std::numeric_limits<double>::max(); // set to maximum
    
    int innerAOIArea = (innerAOI.wdth - 1) * (innerAOI.hght - 1);
    
    for (int y = 0; y < searchAOI.hght - innerAOI.hght; y = y + stepSize)
    {
        for (int x = 0; x < searchAOI.wdth - innerAOI.wdth; x = x + stepSize)
        {
            // vertices of inner square
            
            int topLeftX = x;
            int topLeftY = y;
            
            int backRightX = topLeftX + (innerAOI.wdth  - 1);
            int backRightY = topLeftY + (innerAOI.hght - 1);
            
            int topLeftIndex  = searchAOI.wdth * topLeftY + topLeftX;
            int topRghtIndex  = topLeftIndex + (innerAOI.wdth  - 1);
            int backLeftIndex = topLeftIndex + (innerAOI.hght - 1) * searchAOI.wdth;
            int backRghtIndex = topRghtIndex + backLeftIndex - topLeftIndex;
            
            // calculate glint intensity
            
            double glintIntensity = 0.0;
            double glintArea = 0.0;
            
            bool glintWithinHaarDetector = false; // flag for glint overlap
            
            std::vector<int> z(4);
            z[0] = 0;
            z[1] = 0;
            z[2] = glintAOI.wdth;
            z[3] = glintAOI.wdth;
            
            // check if glint overlaps with Haar detector
            
            for (int m = 0; m < 4; m++)
            {
                int n = (m + 1) % 4;
                
                if (glintAOI.xPos + z[m] >= topLeftX && glintAOI.xPos + z[m] <= backRightX)
                {
                    if (glintAOI.yPos + z[n] >= topLeftY && glintAOI.yPos + z[n] <= backRightY)
                    {
                        glintWithinHaarDetector = true;
                        break;
                    }
                }
            }
            
            if (glintWithinHaarDetector) // if yes, check how much it overlaps
            {
                bool glintOverlapsLeftEdge   = false;
                bool glintOverlapsRightEdge  = false;
                bool glintOverlapsTopEdge    = false;
                bool glintOverlapsBottomEdge = false;
                
                if (glintAOI.xPos < topLeftX)
                {
                    glintOverlapsLeftEdge = true;
                }
                else if (glintAOI.xPos + glintAOI.wdth > backRightX)
                {
                    glintOverlapsRightEdge = true;
                }
                
                if (glintAOI.yPos < topLeftY)
                {
                    glintOverlapsTopEdge = true;
                }
                else if (glintAOI.yPos + glintAOI.wdth > backRightY)
                {
                    glintOverlapsBottomEdge = true;
                }
                
                // coordinates of corners of glint square
                int glintTopLeftIndex  = searchAOI.wdth *  glintAOI.yPos + glintAOI.xPos;
                int glintTopRghtIndex  = searchAOI.wdth *  glintAOI.yPos + glintAOI.xPos  + glintAOI.wdth;
                int glintBackLeftIndex = searchAOI.wdth * (glintAOI.yPos + glintAOI.wdth) + glintAOI.xPos;
                int glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                
                // check if glint square overlaps with edge or corner of pupil square
                
                // check edge overlap
                
                if (!glintOverlapsLeftEdge && !glintOverlapsRightEdge && glintOverlapsTopEdge) // top edge
                {
                    glintTopLeftIndex  = searchAOI.wdth * topLeftY + glintAOI.xPos;
                    glintTopRghtIndex  = searchAOI.wdth * topLeftY + glintAOI.xPos + glintAOI.wdth;
                    glintBackLeftIndex = searchAOI.wdth * (glintAOI.yPos + glintAOI.wdth) + glintAOI.xPos;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                
                if (!glintOverlapsLeftEdge && !glintOverlapsRightEdge && glintOverlapsBottomEdge) // bottom edge
                {
                    glintTopLeftIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos;
                    glintTopRghtIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos + glintAOI.wdth;
                    glintBackLeftIndex = searchAOI.wdth * backRightY + glintAOI.xPos;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                
                if (glintOverlapsLeftEdge && !glintOverlapsTopEdge && !glintOverlapsBottomEdge) // left edge
                {
                    glintTopLeftIndex  = searchAOI.wdth * glintAOI.yPos + topLeftX;
                    glintTopRghtIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos + glintAOI.wdth;
                    glintBackLeftIndex = searchAOI.wdth * (glintAOI.yPos + glintAOI.wdth) + topLeftX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && !glintOverlapsTopEdge && !glintOverlapsBottomEdge) // right edge
                {
                    glintTopLeftIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos;
                    glintTopRghtIndex  = searchAOI.wdth * glintAOI.yPos + backRightX;
                    glintBackLeftIndex = searchAOI.wdth * (glintAOI.yPos + glintAOI.wdth) + glintAOI.xPos;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                // check corner overlap
                
                if (glintOverlapsLeftEdge && glintOverlapsTopEdge) // top left corner
                {
                    glintTopLeftIndex  = topLeftIndex;
                    glintTopRghtIndex  = searchAOI.wdth * topLeftY + glintAOI.xPos + glintAOI.wdth;
                    glintBackLeftIndex = searchAOI.wdth * (glintAOI.yPos + glintAOI.wdth) + topLeftX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && glintOverlapsTopEdge) // top right corner
                {
                    glintTopLeftIndex  = searchAOI.wdth * topLeftY + glintAOI.xPos;
                    glintTopRghtIndex  = topRghtIndex ;
                    glintBackLeftIndex = searchAOI.wdth * (glintAOI.yPos + glintAOI.wdth) + glintAOI.xPos;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsLeftEdge && glintOverlapsBottomEdge) // bottom left corner
                {
                    glintTopLeftIndex  = searchAOI.wdth * glintAOI.yPos + topLeftX;
                    glintTopRghtIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos + glintAOI.wdth;
                    glintBackLeftIndex = backLeftIndex;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && glintOverlapsBottomEdge) // bottom right corner
                {
                    glintTopLeftIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos;
                    glintTopRghtIndex  = searchAOI.wdth * glintAOI.yPos + backRightX;
                    glintBackLeftIndex = searchAOI.wdth * backRightY + glintAOI.xPos;
                    glintBackRghtIndex = backRghtIndex;
                }
                
                // calculate area and intensity of glint
                glintIntensity = I[glintBackRghtIndex] - I[glintBackLeftIndex] - I[glintTopRghtIndex] + I[glintTopLeftIndex];
                glintArea      = (glintTopRghtIndex - glintTopLeftIndex) * ((glintBackLeftIndex - glintTopLeftIndex) / searchAOI.wdth);
            }
            
            // calculate average pupil intensity, adjusting for glint
            
            double pupilIntensity = ((I[backRghtIndex] - I[backLeftIndex] - I[topRghtIndex] + I[topLeftIndex]) - glintIntensity) / (innerAOIArea - glintArea);
            
            if (pupilIntensity < minPupilIntensity)
            {
                innerAOI.xPos = topLeftX;
                innerAOI.yPos = topLeftY;
                minPupilIntensity = pupilIntensity;
            }
        }
    }

    return innerAOI;
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

bool findEdge(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int& startEdgeIndex, double pupilXCentre, double pupilYCentre)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    // Find a starting edge point using Starburst type algorithm

    bool NEW_EDGE_FOUND = false;

    for (int m = 0; m < 8 && !NEW_EDGE_FOUND; m++)
    {
        bool STOP_SEARCH = false;

        int x = round(pupilXCentre);
        int y = round(pupilYCentre);

        while (!STOP_SEARCH)
        {
            for (int n = 0; n < 2 && !STOP_SEARCH; n++)
            {
                x = x + dX[m] * (1 - n);
                y = y + dY[m] * (n);

                if (x < 0 || x >= mAOI.wdth || y < 0 || y >= mAOI.hght) { STOP_SEARCH = true; break; }

                int centreIndex = y * mAOI.wdth + x;

                if (cannyEdgeVector[centreIndex] == 1)
                {
                    startEdgeIndex = centreIndex;
                    NEW_EDGE_FOUND = true;
                    STOP_SEARCH    = true;
                }
                else if (cannyEdgeVector[centreIndex] > 1) { STOP_SEARCH = true; }
            }
        }
    }

    return NEW_EDGE_FOUND;
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

            for (int iEdge = 0; iEdge < numEdges || branchNumber == 0; iEdge++)
            {
                branchProperties branchNew;
                if (branchNumber == 0) { branchNew.pointIndices.push_back(vertexCurrent.pointIndex); } // first vertex should be included in first branch

                if (numEdges > 0)
                {
                    int edgePointIndexOld = vertexCurrent.connectedPoints[iEdge];
                    branchNew.pointIndices.push_back(edgePointIndexOld);

                    while(true) // keep going until vertex is encountered
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
                            vertexNew.connectedPoints = edgePointsNew;
                            verticesNew.push_back(vertexNew); // new vertices to be checked
                            break;
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

std::vector<edgeProperties> edgeSelection(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, double pupilXCentre, double pupilYCentre)
{
    std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all edges

    std::vector<int> startIndices;
    int startIndex = 0;

    while(true)
    {
        if (!findEdge (cannyEdgeVector, mAOI, startIndex, pupilXCentre, pupilYCentre)) { break; } // no (more) edges found
        startIndices.push_back(findEdgePoints(cannyEdgeVector, mAOI, startIndex)); // tag all desired edges
    }

    int numEdges = startIndices.size();

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

std::vector<double> calculateCurvatures(edgeProperties& mEdgeProperties, const std::vector<double>& edgeXTangents, const std::vector<double>& edgeYTangents)
{
    int edgeSize = edgeXTangents.size();

    mEdgeProperties.xnormals.resize(edgeSize);
    mEdgeProperties.ynormals.resize(edgeSize);
    mEdgeProperties.curvatures.resize(edgeSize);

    int numPos = 0;
    int numNeg = 0;

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
    {
        // calculate window tangents

        // first window

        double meanXTangent_1 = calculateMean(std::vector<double>(&edgeXTangents[iEdgePoint - curvatureWindowLength],&edgeXTangents[iEdgePoint - 1]));
        double meanYTangent_1 = calculateMean(std::vector<double>(&edgeYTangents[iEdgePoint - curvatureWindowLength],&edgeYTangents[iEdgePoint - 1]));

        // second window

        double meanXTangent_2 = calculateMean(std::vector<double>(&edgeXTangents[iEdgePoint + 1],&edgeXTangents[iEdgePoint + curvatureWindowLength]));
        double meanYTangent_2 = calculateMean(std::vector<double>(&edgeYTangents[iEdgePoint + 1],&edgeYTangents[iEdgePoint + curvatureWindowLength]));

        // calculate vector difference

        double vectorAngle = atan2(meanYTangent_2, meanXTangent_2) - atan2(meanYTangent_1, meanXTangent_1);

        if      (vectorAngle >  M_PI) { vectorAngle = vectorAngle - 2 * M_PI; }
        else if (vectorAngle < -M_PI) { vectorAngle = vectorAngle + 2 * M_PI; }

        if      (vectorAngle > 0) { numPos++; }
        else if (vectorAngle < 0) { numNeg++; }

        mEdgeProperties.curvatures[iEdgePoint] = 180 * vectorAngle / M_PI; // in degrees

        mEdgeProperties.xnormals[iEdgePoint] = meanXTangent_2 - meanXTangent_1;
        mEdgeProperties.ynormals[iEdgePoint] = meanYTangent_2 - meanYTangent_1;
    }

    if (numNeg > numPos) // if majority sign is negative, then swap all signs
    {
        for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
        {
            mEdgeProperties.curvatures[iEdgePoint] = -mEdgeProperties.curvatures[iEdgePoint];
        }
    }

    return mEdgeProperties.curvatures;
}

std::vector<edgeProperties> edgeSegmentationCurvature(edgeProperties mEdgeProperties, const double curvatureLowerLimit, const double curvatureUpperLimit)
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

    // cut edge at breakpoints

    std::vector<edgeProperties> vEdgePropertiesNew;

    int numBreakPoints = breakPoints.size();

    for (int iBreakPoint = 0; iBreakPoint < numBreakPoints - 1; iBreakPoint++)
    {
        int iStartBreakPoint = breakPoints[iBreakPoint] + 1;
        int iEndBreakPoint   = breakPoints[iBreakPoint  + 1];

        edgeProperties mEdgePropertiesNew;
        mEdgePropertiesNew.pointIndices = std::vector<int>    (mEdgeProperties.pointIndices.begin() + iStartBreakPoint, mEdgeProperties.pointIndices.begin() + iEndBreakPoint);
        mEdgePropertiesNew.curvatures   = std::vector<double> (mEdgeProperties.curvatures.begin()   + iStartBreakPoint, mEdgeProperties.curvatures.begin()   + iEndBreakPoint);
        mEdgePropertiesNew.xnormals     = std::vector<double> (mEdgeProperties.xnormals.begin()     + iStartBreakPoint, mEdgeProperties.xnormals.begin()     + iEndBreakPoint);
        mEdgePropertiesNew.ynormals     = std::vector<double> (mEdgeProperties.ynormals.begin()     + iStartBreakPoint, mEdgeProperties.ynormals.begin()     + iEndBreakPoint);


        vEdgePropertiesNew.push_back(mEdgePropertiesNew);
    }

    return vEdgePropertiesNew;
}

std::vector<edgeProperties> edgeSegmentationLength(const detectionProperties& mDetectionProperties, const edgeProperties& mEdgeProperties, bool USE_PRIOR_POSITION)
{
    int edgeSize = mEdgeProperties.curvatures.size();

    // find breakpoints based on length thresholding

    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(-1); // add first point (+ 1 is added later)

    if (edgeSize > lengthWindowLength)
    {
        double lengthDifference = mEdgeProperties.length - mDetectionProperties.v.predictionCircumference;

        if (lengthDifference > lengthWindowLength) // if edge is significantly longer than prediction, cut terminals to make it shorter
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
            //            mEdgePropertiesNew.xnormals     = std::vector<double>(mEdgeProperties.xnormals.begin()     + iStartBreakPoint, mEdgeProperties.xnormals.begin()     + iEndBreakPoint);
            //            mEdgePropertiesNew.ynormals     = std::vector<double>(mEdgeProperties.ynormals.begin()     + iStartBreakPoint, mEdgeProperties.ynormals.begin()     + iEndBreakPoint);

            mEdgePropertiesNew.intensity = calculateMeanInt(mEdgePropertiesNew.intensities);
            mEdgePropertiesNew.gradient  = calculateMeanInt(mEdgePropertiesNew.gradients);
            mEdgePropertiesNew.radius    = calculateMean(mEdgePropertiesNew.radii);
            mEdgePropertiesNew.curvature = calculateMean(mEdgePropertiesNew.curvatures);

            vEdgeProperties[iBreakPoint] = mEdgePropertiesNew;
        }

        // only remove one of the two edge terminals - re-attach other one

        // compare edge characteristics of terminals with middle section

        // 0 = start terminal
        // 1 = main section
        // 2 = end terminal

        double dIntensity_0 = std::abs(vEdgeProperties[0].intensity - vEdgeProperties[1].intensity);
        double dIntensity_2 = std::abs(vEdgeProperties[2].intensity - vEdgeProperties[1].intensity);

        double dCurvature_0 = std::abs(vEdgeProperties[0].curvature - vEdgeProperties[1].curvature);
        double dCurvature_2 = std::abs(vEdgeProperties[2].curvature - vEdgeProperties[1].curvature);

        double dRadius_0 = std::abs(vEdgeProperties[0].radius - vEdgeProperties[1].radius);
        double dRadius_2 = std::abs(vEdgeProperties[2].radius - vEdgeProperties[1].radius);

        double dGradient_0 = std::abs(vEdgeProperties[0].gradient - vEdgeProperties[1].gradient);
        double dGradient_2 = std::abs(vEdgeProperties[2].gradient - vEdgeProperties[1].gradient);

        double scoreIntensity_0 = calculateScoreIntensity(dIntensity_0);
        double scoreIntensity_2 = calculateScoreIntensity(dIntensity_2);

        double scoreCurvature_0 = calculateScoreCurvature(dCurvature_0);
        double scoreCurvature_2 = calculateScoreCurvature(dCurvature_2);

        double scoreRadius_0 = 0;
        double scoreRadius_2 = 0;

        double scoreGradient_0 = 0;
        double scoreGradient_2 = 0;

        if (USE_PRIOR_POSITION)
        {
            scoreRadius_0 = calculateScoreRadius(dRadius_0);
            scoreRadius_2 = calculateScoreRadius(dRadius_2);

            scoreGradient_0 = calculateScoreGradient(dGradient_0);
            scoreGradient_2 = calculateScoreGradient(dGradient_2);
        }

        double scoreTotal_0 = scoreFactorIntensity * scoreIntensity_0 + scoreFactorCurvature * scoreCurvature_0 + scoreFactorRadius * scoreRadius_0 + scoreFactorGradient * scoreGradient_0;
        double scoreTotal_2 = scoreFactorIntensity * scoreIntensity_2 + scoreFactorCurvature * scoreCurvature_2 + scoreFactorRadius * scoreRadius_2 + scoreFactorGradient * scoreGradient_2;

        int indexStart = 0;
        int indexEnd   = 0;

        if      (scoreTotal_0 > scoreTotal_2) { indexStart = 0; indexEnd = 1; }
        else if (scoreTotal_0 < scoreTotal_2) { indexStart = 1; indexEnd = 2; }

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

            //            mEdgePropertiesNew.xnormals.reserve(vEdgeProperties[indexStart].xnormals.size() + vEdgeProperties[indexEnd].xnormals.size());
            //            mEdgePropertiesNew.xnormals.insert(mEdgePropertiesNew.xnormals.end(), vEdgeProperties[indexStart].xnormals.begin(), vEdgeProperties[indexStart].xnormals.end() );
            //            mEdgePropertiesNew.xnormals.insert(mEdgePropertiesNew.xnormals.end(), vEdgeProperties[indexEnd].xnormals.begin(),   vEdgeProperties[indexEnd].xnormals.end() );

            //            mEdgePropertiesNew.ynormals.reserve(vEdgeProperties[indexStart].ynormals.size() + vEdgeProperties[indexEnd].ynormals.size());
            //            mEdgePropertiesNew.ynormals.insert(mEdgePropertiesNew.ynormals.end(), vEdgeProperties[indexStart].ynormals.begin(), vEdgeProperties[indexStart].ynormals.end() );
            //            mEdgePropertiesNew.ynormals.insert(mEdgePropertiesNew.ynormals.end(), vEdgeProperties[indexEnd].ynormals.begin(),   vEdgeProperties[indexEnd].ynormals.end() );

            std::vector<edgeProperties> vEdgePropertiesNew(2);
            vEdgePropertiesNew[0] = mEdgePropertiesNew;
            vEdgePropertiesNew[1] = vEdgeProperties[(indexStart + 2) % 3];

            vEdgeProperties = vEdgePropertiesNew; // update return vector
        }
    }
    else
    {
        vEdgeProperties[0] = mEdgeProperties;
    }

    return vEdgeProperties;
}

std::vector<int> calculateRadialGradients(const cv::Mat& img, const std::vector<int> edgeIndices, double pupilXCentre, double pupilYCentre, int kernelSize)
{
    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1 };
    std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1 };

    uchar *ptr = img.data;
    int width  = img.cols;
    int height = img.rows;

    int kernelRadius = (kernelSize - 1) / 2;

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

        if (xPos < 0 || xPos >= width || yPos < 0 || yPos >= height) { continue; }
        if (xNeg < 0 || xNeg >= width || yNeg < 0 || yNeg >= height) { continue; }

        int neighbourIndexPositive = width * yPos + xPos;
        int neighbourIndexNegative = width * yNeg + xNeg;

        double positiveIntensity = ptr[neighbourIndexPositive];
        double negativeIntensity = ptr[neighbourIndexNegative];

        gradientVector.push_back(positiveIntensity - negativeIntensity);
    }

    return gradientVector;
}

std::vector<int> findEdgeIntensities(const cv::Mat& img, const edgeProperties& mEdgeProperties, AOIProperties mAOI)
{
    int edgeSize = mEdgeProperties.pointIndices.size();

    // calculate pixel intensities within inner curve of edge points

    std::vector<int> edgeIntensities(edgeSize);

    uchar *ptr_img = img.data;

    for (int iEdgePoint = 0; iEdgePoint < edgeSize; iEdgePoint++)
    {
        int edgePointIndex = mEdgeProperties.pointIndices[iEdgePoint];
        int edgePointXPos  = edgePointIndex % mAOI.wdth;
        int edgePointYPos  = (edgePointIndex - edgePointXPos) / mAOI.wdth;

        int offsetXPos = edgePointXPos + edgeIntensitiesPositionOffset * ceil2(mEdgeProperties.xnormals[iEdgePoint]);
        int offsetYPos = edgePointYPos + edgeIntensitiesPositionOffset * ceil2(mEdgeProperties.ynormals[iEdgePoint]);

        if (edgePointXPos < 0 || edgePointXPos > mAOI.wdth || edgePointYPos < 0 || edgePointYPos > mAOI.hght)
        {       edgeIntensities[iEdgePoint] = (int) ptr_img[edgePointXPos + edgePointYPos * mAOI.wdth]; }
        else {  edgeIntensities[iEdgePoint] = (int) ptr_img[   offsetXPos +    offsetYPos * mAOI.wdth]; }
    }

    return edgeIntensities;
}

edgeProperties edgePointFilter(const edgeProperties& mEdgeProperties, double pointScoreThreshold)
{
    edgeProperties mEdgePropertiesNew = mEdgeProperties;

    for (int iEdgePoint = 0, edgeSize =  mEdgeProperties.pointIndices.size(); iEdgePoint < edgeSize; iEdgePoint++)
    {
        double score;

        // calculate some score here

        if (score > pointScoreThreshold)
        {

        }

    }

    return mEdgePropertiesNew;
}

double calculateEdgeLength(const std::vector<int>& edgePoints, AOIProperties mAOI)
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

void calculateCurvatureRange(edgeProperties& mEdgeProperties)
{
    // Calculate min, max and mean curvature

    double curvature    = 0;
    double curvatureMax = 0;
    double curvatureMin = 0;

    std::vector<double> edgeCurvaturesNew;

    int edgeSize = mEdgeProperties.curvatures.size();

    for (int iEdgePoint = 0; iEdgePoint < edgeSize; iEdgePoint++)
    {
        double curvature = mEdgeProperties.curvatures[iEdgePoint];
        if (curvature < 180.0) { edgeCurvaturesNew.push_back(curvature); }
    }

    int edgeLengthNew = edgeCurvaturesNew.size();

    if (edgeLengthNew > 0)
    {
        for (int iEdgePoint = 0; iEdgePoint < edgeLengthNew; iEdgePoint++)
        {
            curvature += std::abs(edgeCurvaturesNew[iEdgePoint] / edgeLengthNew);
        }

        curvatureMax = *std::max_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
        curvatureMin = *std::min_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
    }
    else { curvature = 360; curvatureMax = 360; curvatureMin = 360; }

    mEdgeProperties.curvature = curvature;
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

std::vector<int> edgeClassification(detectionProperties mDetectionProperties, const std::vector<edgeProperties>& vEdgePropertiesAll, double curvaturePrediction, bool USE_PRIOR_POSITION, bool USE_PRIOR_FEATURES)
{
    int numEdgesMax = mDetectionProperties.p.ellipseFitNumberMaximum;
    int numEdges    = vEdgePropertiesAll.size();
    if (numEdgesMax > numEdges) { numEdgesMax = numEdges; }

    // Classify edges based on score

    std::vector<double> totalScores(numEdges);

    const double norm = scoreFactorRadius + scoreFactorCurvature + scoreFactorCircumference + scoreFactorIntensity;

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        double dRadius          = mDetectionProperties.v.predictionCircumference / (2 * M_PI) - vEdgePropertiesAll[iEdge].radius;
        double dCircumference   = mDetectionProperties.v.predictionCircumference - vEdgePropertiesAll[iEdge].length;
        double dCurvature       = curvaturePrediction - vEdgePropertiesAll[iEdge].curvature;
        double dIntensity       = vEdgePropertiesAll[iEdge].intensity - mDetectionProperties.v.averageIntensity;

        double scoreRadius          = 0;
        double scoreCurvature       = 0;
        double scoreCircumference   = 0;
        double scoreIntensity       = 0;

        if (USE_PRIOR_POSITION)
        {
            scoreRadius = scoreFactorRadius * calculateScoreRadius(dRadius);
        }

        if (USE_PRIOR_FEATURES)
        {
            scoreCurvature     = scoreFactorCurvature     * calculateScoreCurvature(dCurvature);
            scoreCircumference = scoreFactorCircumference * calculateScoreCircumference(dCircumference);
        }

        scoreIntensity = scoreFactorIntensity * calculateScoreIntensity(dIntensity);

        totalScores[iEdge] = (scoreRadius + scoreCurvature + scoreCircumference + scoreIntensity) / norm;
    }

    // Only pick edges above threshold

    std::vector<int> pupilEdges;

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        double factor = scoreFactorIntensity + USE_PRIOR_POSITION * scoreFactorRadius + USE_PRIOR_FEATURES * (scoreFactorCurvature + scoreFactorCircumference);
        factor = factor / norm; // equal 1 if all priors are active

        if (totalScores[iEdge] >= factor * mDetectionProperties.p.scoreThreshold) { pupilEdges.push_back(iEdge); }
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
    mEllipseProperties.PUPIL_DETECTED = true;
    double h = pow((semiMajor - semiMinor), 2) / pow((semiMajor + semiMinor), 2);
    mEllipseProperties.circumference = M_PI * (semiMajor + semiMinor) * (1 + (3 * h) / (10 + sqrt(4 - 3 * h))); // ramanujans 2nd approximation
    mEllipseProperties.aspectRatio   = semiMinor / semiMajor;
    mEllipseProperties.radius        = 0.5 * (semiMinor + semiMajor);
    mEllipseProperties.xPos          = ellipseParameters[2];
    mEllipseProperties.yPos          = ellipseParameters[3];
    mEllipseProperties.width         = ellipseParameters[4];
    mEllipseProperties.height        = ellipseParameters[5];
    mEllipseProperties.coefficients  = ellipseFitCoefficients;

    for (int iParameter = 0; iParameter < 6; iParameter++)
    {
        if (std::isnan(ellipseParameters[iParameter])) { mEllipseProperties.PUPIL_DETECTED = false; } // ERROR
    }

    return mEllipseProperties;
}

std::vector<ellipseProperties> getEllipseFits(const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI, detectionProperties mDetectionProperties, bool USE_PRIOR_POSITION, bool USE_PRIOR_FEATURES)
{
    ellipseProperties mEllipseProperties;
    mEllipseProperties.PUPIL_DETECTED = false;

    int totalNumberOfEdges = vEdgePropertiesAll.size(); // total number of edges

    std::vector<ellipseProperties> vEllipsePropertiesAll; // vector to record information for each accepted ellipse fit

    for (int combiNumEdges = totalNumberOfEdges; combiNumEdges >= 1; combiNumEdges--) // loop through all possible edge set sizes
    {
        std::vector<bool> edgeCombination(totalNumberOfEdges);
        std::fill(edgeCombination.begin() + totalNumberOfEdges - combiNumEdges, edgeCombination.end(), true);

        do // loop through all possible edge combinations for the current set size
        {
            std::vector<int> combiEdgeIndices (combiNumEdges);
            std::vector<int> combiEdgeLengths (combiNumEdges);
            std::vector<int> combiEdgeSizes   (combiNumEdges);

            std::vector<std::vector<int>> combiEdgePoints(combiNumEdges);

            for (int iEdge = 0, jEdge = 0; iEdge < totalNumberOfEdges; ++iEdge)
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

            if (USE_PRIOR_FEATURES)
            {
                if (edgeSetLength < mDetectionProperties.v.predictionCircumference * mDetectionProperties.p.edgeLengthFraction)
                { continue; }
            }
            else
            {
                if (edgeSetLength <= mDetectionProperties.p.circumferenceMin * mDetectionProperties.p.edgeLengthFraction)
                { continue; }
            }

            int edgeSetSize = std::accumulate(combiEdgeSizes.begin(), combiEdgeSizes.end(), 0);

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

            int combiWdth = XPosMax - XPosMin;
            int combiHght = YPosMax - YPosMin;

            if (combiWdth * M_PI > mDetectionProperties.p.circumferenceMax || combiHght * M_PI > mDetectionProperties.p.circumferenceMax) { continue; } // no large ellipse
            if (combiWdth * M_PI < mDetectionProperties.p.circumferenceMin && combiHght * M_PI < mDetectionProperties.p.circumferenceMin) { continue; } // no small ellipse

            // Fit ellipse

            ellipseProperties mEllipsePropertiesNew = fitEllipse(edgeIndices, edgeSetSize, mAOI.wdth);

            if (!mEllipsePropertiesNew.PUPIL_DETECTED) { continue; } // error

            // Size and shape filters

            if (mEllipsePropertiesNew.circumference > mDetectionProperties.p.circumferenceMax) { continue; } // no large ellipse
            if (mEllipsePropertiesNew.circumference < mDetectionProperties.p.circumferenceMin) { continue; } // no small ellipse
            if (mEllipsePropertiesNew.aspectRatio   < mDetectionProperties.p.aspectRatioMin)   { continue; } // no extreme deviations from circular shape

            if (USE_PRIOR_POSITION)
            {
                double dX = mEllipsePropertiesNew.xPos - (mDetectionProperties.v.predictionXPos - mAOI.xPos);
                double dY = mEllipsePropertiesNew.yPos - (mDetectionProperties.v.predictionYPos - mAOI.yPos);
                double dR = sqrt(dX * dX + dY * dY);
                if (dR > mDetectionProperties.v.changeThresholdPosition) { continue; } // no large ellipse displacements
            }

            if (USE_PRIOR_FEATURES)
            {
                if (std::abs(mEllipsePropertiesNew.circumference - mDetectionProperties.v.predictionCircumference) > mDetectionProperties.v.changeThresholdCircumference) { continue; } // no large ellipse size changes
                if (std::abs(mEllipsePropertiesNew.aspectRatio   - mDetectionProperties.v.predictionAspectRatio  ) > mDetectionProperties.v.changeThresholdAspectRatio  ) { continue; } // no large ellipse shape changes
            }

            // calculate error between fit and every edge point

            double A = mEllipsePropertiesNew.coefficients[0];
            double B = mEllipsePropertiesNew.coefficients[1];
            double C = mEllipsePropertiesNew.coefficients[2];
            double D = mEllipsePropertiesNew.coefficients[3];
            double E = mEllipsePropertiesNew.coefficients[4];
            double F = mEllipsePropertiesNew.coefficients[5];

            // calculate errors

            std::vector<double> fitErrors(edgeSetSize);

            for (int iEdgePoint = 0; iEdgePoint < edgeSetSize; iEdgePoint++)
            {
                int edgePointIndex = edgeIndices[iEdgePoint];

                double x = edgePointIndex % mAOI.wdth;
                double y = (edgePointIndex - x) / mAOI.wdth;

                fitErrors[iEdgePoint] = std::abs(A * x * x + B * x * y + C * y * y + D * x + E * y + F);
            }

            std::vector<double> fitErrorsSorted = fitErrors;
            std::sort(fitErrorsSorted.begin(), fitErrorsSorted.end());
            std::reverse(fitErrorsSorted.begin(), fitErrorsSorted.end());
            std::vector<double> fitErrorsMax(fitErrorsSorted.begin(), fitErrorsSorted.begin() + round(fitErrorFraction * edgeSetLength));
            double fitErrorMax = calculateMean(fitErrorsMax);

            if (fitErrorMax > mDetectionProperties.p.ellipseFitErrorMaximum) { continue; } // no large fit errors

            mEllipsePropertiesNew.fitError = fitErrorMax;

            // save parameters of accepted fit

            mEllipsePropertiesNew.edgeIndices = combiEdgeIndices;
            mEllipsePropertiesNew.edgeLength  = edgeSetLength;
            vEllipsePropertiesAll.push_back(mEllipsePropertiesNew);
        }
        while (std::next_permutation(edgeCombination.begin(), edgeCombination.end()));
    }

    return vEllipsePropertiesAll;
}

int ellipseFitFilter(detectionProperties mDetectionProperties, std::vector<ellipseProperties> vEllipseProperties, bool USE_PRIOR_POSITION, bool USE_PRIOR_FEATURES)
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

    for (int iFit = 0; iFit < numFits; iFit++)
    {
        ellipseProperties mEllipseProperties = vEllipseProperties[iFit];

        double dx = mEllipseProperties.xPos - mDetectionProperties.v.predictionXPos;
        double dy = mEllipseProperties.yPos - mDetectionProperties.v.predictionYPos;
        double displacementChange  = sqrt(dx * dx + dy * dy);
        double circumferenceChange = (std::abs(mEllipseProperties.circumference - mDetectionProperties.v.predictionCircumference));
        double aspectRatioChange   = (std::abs(mEllipseProperties.aspectRatio   - mDetectionProperties.v.predictionAspectRatio));
        double lengthDifference    = (std::abs(mEllipseProperties.edgeLength    - mDetectionProperties.v.predictionCircumference));
        double fitError            = mEllipseProperties.fitError;

        if (USE_PRIOR_POSITION)
        {
            scoreDisplacement  = (-maxScoreDisplacement  / mDetectionProperties.p.changeThresholdPosition)  * displacementChange  + maxScoreDisplacement;
        }

        if (USE_PRIOR_FEATURES)
        {
            scoreAspectRatio   = (-maxScoreAspectRatio   / mDetectionProperties.p.changeThresholdAspectRatio)   * aspectRatioChange   + maxScoreAspectRatio;
            scoreCircumference = (-maxScoreCircumference / mDetectionProperties.p.changeThresholdCircumference) * circumferenceChange + maxScoreCircumference;
        }

        double scoreFitError = (-maxScoreFitError   / mDetectionProperties.p.ellipseFitErrorMaximum)  * fitError         + maxScoreFitError;
        double scoreLength   = (-maxScoreLength * 2 / mDetectionProperties.v.predictionCircumference) * lengthDifference + maxScoreLength;

        if (scoreCircumference  < 0) { scoreCircumference   = 0; }
        if (scoreAspectRatio    < 0) { scoreAspectRatio     = 0; }
        if (scoreFitError       < 0) { scoreFitError        = 0; }
        if (scoreLength         < 0) { scoreLength          = 0; }

        featureChange[iFit] = scoreAspectRatio + scoreCircumference + scoreDisplacement + scoreFitError + scoreLength;
    }

    int acceptedFitIndex = std::distance(featureChange.begin(), std::max_element(featureChange.begin(), featureChange.end()));

    return acceptedFitIndex;
}

void checkVariableLimits(detectionProperties& mDetectionProperties)
{
    if (mDetectionProperties.v.changeThresholdAspectRatio      < mDetectionProperties.p.changeThresholdAspectRatio)
    {        mDetectionProperties.v.changeThresholdAspectRatio = mDetectionProperties.p.changeThresholdAspectRatio; }

    if (mDetectionProperties.v.changeThresholdCircumference      < mDetectionProperties.p.changeThresholdCircumference)
    {        mDetectionProperties.v.changeThresholdCircumference = mDetectionProperties.p.changeThresholdCircumference; }

    if (mDetectionProperties.v.changeThresholdPosition      < mDetectionProperties.p.changeThresholdPosition)
    {        mDetectionProperties.v.changeThresholdPosition = mDetectionProperties.p.changeThresholdPosition; }

    if (mDetectionProperties.v.curvatureOffset      < mDetectionProperties.p.curvatureOffset)
    {        mDetectionProperties.v.curvatureOffset = mDetectionProperties.p.curvatureOffset; }
}

detectionProperties pupilDetection(const cv::Mat& imageOriginalBGR, detectionProperties mDetectionProperties)
{
    // Define some variables

    detectionProperties mDetectionPropertiesNew = mDetectionProperties; // new properties for new frame
    mDetectionPropertiesNew.m.ERROR_DETECTED = false;

    ellipseProperties mEllipseProperties;
    mEllipseProperties.PUPIL_DETECTED = false;

    checkVariableLimits(mDetectionProperties);

    // Define search area

    int imageWdth = imageOriginalBGR.cols;
    int imageHght = imageOriginalBGR.rows;

    AOIProperties searchAOI;

    searchAOI.xPos = round(mDetectionProperties.v.predictionXPos - mDetectionProperties.v.searchRadius);
    searchAOI.yPos = round(mDetectionProperties.v.predictionYPos - mDetectionProperties.v.searchRadius);

    int searchEndX = round(mDetectionProperties.v.predictionXPos + mDetectionProperties.v.searchRadius);
    int searchEndY = round(mDetectionProperties.v.predictionYPos + mDetectionProperties.v.searchRadius);

    if (searchAOI.xPos < mDetectionProperties.p.AOIXPos)
    {   searchAOI.xPos = mDetectionProperties.p.AOIXPos; }

    if (searchAOI.yPos < mDetectionProperties.p.AOIYPos)
    {   searchAOI.yPos = mDetectionProperties.p.AOIYPos; }

    if (searchEndX > mDetectionProperties.p.AOIXPos + mDetectionProperties.p.AOIWdth - 1)
    {   searchEndX = mDetectionProperties.p.AOIXPos + mDetectionProperties.p.AOIWdth - 1; }

    if (searchEndY > mDetectionProperties.p.AOIYPos + mDetectionProperties.p.AOIHght - 1)
    {   searchEndY = mDetectionProperties.p.AOIYPos + mDetectionProperties.p.AOIHght - 1; }

    searchAOI.wdth = searchEndX - searchAOI.xPos + 1;
    searchAOI.hght = searchEndY - searchAOI.yPos + 1;

    AOIProperties innerAOI;
    innerAOI.wdth = round(pupilHaarReductionFactor * mDetectionProperties.v.predictionWidth);
    innerAOI.hght = round(pupilHaarReductionFactor * mDetectionProperties.v.predictionHeight);
    if (innerAOI.wdth > searchAOI.wdth) { innerAOI.wdth = searchAOI.wdth; }
    if (innerAOI.hght > searchAOI.hght) { innerAOI.hght = searchAOI.hght; }

    AOIProperties outerAOI = mDetectionProperties.m.outerAOI;

    if (innerAOI.wdth > 0 && innerAOI.hght > 0)
    {
        // Check if prior pupil information should be used

        bool USE_PRIOR_POSITION = false;
        bool USE_PRIOR_FEATURES = false;

        if (mDetectionProperties.v.certaintyPosition >= 0.5) { USE_PRIOR_POSITION = true; }
        if (mDetectionProperties.v.certaintyFeatures >= 0.5) { USE_PRIOR_FEATURES = true; }

        // Convert to grayscale

        cv::Mat imageOriginalGray;
        cv::cvtColor(imageOriginalBGR, imageOriginalGray, cv::COLOR_BGR2GRAY);

        ////////////////////////////////////////////////////////////////////
        /////////////////////// INITIAL DETECTION  /////////////////////////
        ////////////////////////////////////////////////////////////////////

        std::vector<unsigned int> integralImage = calculateIntImg(imageOriginalGray, imageWdth, searchAOI);

        AOIProperties glintAOI;
        glintAOI.wdth = mDetectionProperties.p.glintWdth;
        glintAOI.hght = glintAOI.wdth;

        glintAOI      = detectGlint(imageOriginalGray, imageWdth, searchAOI, glintAOI);
        glintAOI.xPos = searchAOI.xPos + glintAOI.xPos;
        glintAOI.yPos = searchAOI.yPos + glintAOI.yPos;

        innerAOI      = detectPupilApprox(integralImage, searchAOI, innerAOI, glintAOI);
        innerAOI.xPos = searchAOI.xPos + innerAOI.xPos;
        innerAOI.yPos = searchAOI.yPos + innerAOI.yPos;

        outerAOI.xPos = innerAOI.xPos -     mDetectionProperties.p.haarOffset;
        outerAOI.yPos = innerAOI.yPos -     mDetectionProperties.p.haarOffset;
        outerAOI.wdth = innerAOI.wdth + 2 * mDetectionProperties.p.haarOffset;
        outerAOI.hght = innerAOI.hght + 2 * mDetectionProperties.p.haarOffset;

        // Check limits

        if (outerAOI.wdth >= imageWdth) { outerAOI.wdth = imageWdth - 1; }
        if (outerAOI.hght >= imageHght) { outerAOI.hght = imageHght - 1; }
        if (outerAOI.xPos < 0) { outerAOI.xPos = 0; }
        if (outerAOI.yPos < 0) { outerAOI.yPos = 0; }
        if (outerAOI.xPos + outerAOI.wdth >= imageWdth) { outerAOI.wdth = imageWdth - outerAOI.xPos - 1; }
        if (outerAOI.yPos + outerAOI.hght >= imageHght) { outerAOI.hght = imageHght - outerAOI.yPos - 1; }

        // Crop image to outer pupil Haar region

        cv::Rect outerRect(outerAOI.xPos, outerAOI.yPos, outerAOI.wdth, outerAOI.hght);
        cv::Mat imageAOIBGR = imageOriginalBGR(outerRect);

        // Convert back to grayscale

        cv::Mat imageAOIGray;
        cv::cvtColor(imageAOIBGR, imageAOIGray, cv::COLOR_BGR2GRAY);

        ///////////////////////////////////////////////////////////////////////
        /////////////////////// CANNY EDGE DETECTION  /////////////////////////
        ///////////////////////////////////////////////////////////////////////

        cv::Mat imageAOIGrayBlurred;
        int cannyBlurLevel = 2 * mDetectionProperties.p.cannyBlurLevel - 1; // should be odd
        if (cannyBlurLevel > 0) { cv::GaussianBlur(imageAOIGray, imageAOIGrayBlurred, cv::Size(cannyBlurLevel, cannyBlurLevel), 0, 0);
        } else                  { imageAOIGrayBlurred = imageAOIGray; }

        // binary vector

        cv::Mat imageCannyEdges;
        cv::Canny(imageAOIGrayBlurred, imageCannyEdges, mDetectionProperties.p.cannyThresholdHigh, mDetectionProperties.p.cannyThresholdLow, 5);
        std::vector<int> cannyEdges = cannyConversion(imageCannyEdges, outerAOI);
        std::vector<int> cannyEdgesSharpened = sharpenEdges(cannyEdges, outerAOI); // Morphological operation
        std::vector<int> edgeIndices = getEdgeIndices(cannyEdgesSharpened, 1); // used for drawing function

        /////////////////////////////////////////////////////////////////////////////
        //////////////////////////// EDGE SELECTION   ///////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        double xPosPredictionRelative;
        double yPosPredictionRelative;

        if (USE_PRIOR_POSITION)
        {
            xPosPredictionRelative = mDetectionProperties.v.predictionXPos - outerAOI.xPos;
            yPosPredictionRelative = mDetectionProperties.v.predictionYPos - outerAOI.yPos;
        }
        else // use centre of haar-like detector
        {
            xPosPredictionRelative = 0.5 * outerAOI.wdth;
            yPosPredictionRelative = 0.5 * outerAOI.hght;
        }

        std::vector<edgeProperties> vEdgePropertiesAll = edgeSelection(cannyEdgesSharpened, outerAOI, xPosPredictionRelative, yPosPredictionRelative);
        vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

        /////////////////////////////////////////////////////////////////////////////
        //////////////////////////// EDGE SEGMENTATION   ////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        // Calculate curvature limits

        int arrayWidth = arrayCircumferences.size();

        int arrayXPos = 0;
        for (int x = 0; x < arrayWidth; x++)
        {
            if (arrayCircumferences[x] < mDetectionProperties.v.predictionCircumference)
            {
                arrayXPos = x;
                break;
            }
        }

        int arrayHeight = arrayAspectRatios.size();

        int arrayYPos = 0;
        for (int y = 0; y < arrayHeight; y++)
        {
            if (arrayAspectRatios[y] < mDetectionProperties.v.predictionAspectRatio)
            {
                arrayYPos = y;
                break;
            }
        }

        double curvatureUpperLimit = arrayCurvatureMax[arrayXPos * arrayWidth + arrayYPos] + mDetectionProperties.v.curvatureOffset;
        double curvatureLowerLimit = arrayCurvatureMin[arrayXPos * arrayWidth + arrayYPos] - mDetectionProperties.v.curvatureOffset;

        // Curvature calculation and segmentation

        std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            int edgeSize = mEdgeProperties.pointIndices.size();

            std::vector<double> edgeXTangents(edgeSize);
            std::vector<double> edgeYTangents(edgeSize);

            calculateEdgeDirections(mEdgeProperties.pointIndices, edgeXTangents, edgeYTangents, outerAOI);
            calculateCurvatures(mEdgeProperties, edgeXTangents, edgeYTangents);

            if (USE_PRIOR_FEATURES) // only do curvature segmentation if circumference and aspect ratio predictions are accurate
            {
                std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationCurvature(mEdgeProperties, curvatureLowerLimit, curvatureUpperLimit);
                vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end()); // record edges
            }
            else
            {
                vEdgePropertiesNew.push_back(mEdgeProperties);
            }
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
        vEdgePropertiesNew.clear();

        vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

        // Calculate additional edge properties and do length segmentation

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            mEdgeProperties.length      = calculateEdgeLength(mEdgeProperties.pointIndices, outerAOI);
            mEdgeProperties.radii       = calculateEdgeRadii(mEdgeProperties, outerAOI, xPosPredictionRelative, yPosPredictionRelative);
            mEdgeProperties.gradients   = calculateRadialGradients(imageAOIGray, mEdgeProperties.pointIndices, xPosPredictionRelative, yPosPredictionRelative, mDetectionProperties.p.cannyKernelSize);
            mEdgeProperties.intensities = findEdgeIntensities(imageAOIGray, mEdgeProperties, outerAOI);

            if (USE_PRIOR_FEATURES) // Prior circumference is required
            {
                std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationLength(mDetectionProperties, mEdgeProperties, USE_PRIOR_POSITION);
                vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end()); // record edges
            }
            else
            {
                vEdgePropertiesNew.push_back(mEdgeProperties);
            }
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
        vEdgePropertiesNew.clear();

        vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

        ////////////////////////////////////////////////////////////////////////
        //////////////////////// EDGE POINT FILTER   ///////////////////////////
        ////////////////////////////////////////////////////////////////////////

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties    = vEdgePropertiesAll[iEdge];
            edgeProperties mEdgePropertiesNew = edgePointFilter(mEdgeProperties, mDetectionProperties.p.scoreThresholdPoints);
            vEdgePropertiesNew.push_back(mEdgePropertiesNew);
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
        vEdgePropertiesNew.clear();

        vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

        ///////////////////////////////////////////////////////////////////////////
        ////////////////////////// EDGE CLASSIFICATION  ///////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        // Calculate some edge properties

        for (int iEdge = 0, jEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            // re-evaluate edge properties

            mEdgeProperties.length    = calculateEdgeLength(mEdgeProperties.pointIndices, outerAOI);

            mEdgeProperties.intensity = calculateMeanInt(mEdgeProperties.intensities);
            mEdgeProperties.radius    = calculateMean(mEdgeProperties.radii);
            mEdgeProperties.gradient  = calculateMeanInt(mEdgeProperties.gradients);

            calculateCurvatureRange(mEdgeProperties);

            mEdgeProperties.index     = iEdge;
            mEdgeProperties.tag       = 0;

            restoreEdgePoints(mEdgeProperties, cannyEdgesSharpened, outerAOI); // Restore some points
            mEdgeProperties.size      = mEdgeProperties.pointIndices.size();

            vEdgePropertiesNew.push_back(mEdgeProperties);

            jEdge++;
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
        vEdgePropertiesNew.clear();

        // Do edge classification

        std::vector<ellipseProperties> vEllipsePropertiesAll;

        if (vEdgePropertiesAll.size() > 0) // ignore empty edge collections
        {
            double curvaturePrediction = 0.5 * (curvatureUpperLimit + curvatureLowerLimit);

            std::vector<int> acceptedEdges = edgeClassification(mDetectionProperties, vEdgePropertiesAll, curvaturePrediction, USE_PRIOR_POSITION, USE_PRIOR_FEATURES);

            for (int iEdge = 0, numEdges = acceptedEdges.size(); iEdge < numEdges; iEdge++) // grab accepted edges
            {
                int jEdge = acceptedEdges[iEdge];
                edgeProperties mEdgeProperties = vEdgePropertiesAll[jEdge];
                mEdgeProperties.tag = 1; // new tag
                vEdgePropertiesNew.push_back(mEdgeProperties);
            }

            //////////////////////////////////////////////////////////////////
            /////////////////////// ELLIPSE FITTING  /////////////////////////
            //////////////////////////////////////////////////////////////////

            vEllipsePropertiesAll = getEllipseFits(vEdgePropertiesNew, outerAOI, mDetectionProperties, USE_PRIOR_POSITION, USE_PRIOR_FEATURES); // ellipse fitting

            int numFits = vEllipsePropertiesAll.size();

            int acceptedFitIndex = 0;

            if (numFits > 0)
            {
                mEllipseProperties.PUPIL_DETECTED = true;

                if (numFits > 1) { acceptedFitIndex = ellipseFitFilter(mDetectionProperties, vEllipsePropertiesAll, USE_PRIOR_POSITION, USE_PRIOR_FEATURES); } // grab best fit

                mEllipseProperties = vEllipsePropertiesAll[acceptedFitIndex];

                // Classify fits

                for (int iFit = 0; iFit < numFits; iFit++)
                {
                    if (iFit == acceptedFitIndex) { vEllipsePropertiesAll[iFit].tag = 1; }
                    else                          { vEllipsePropertiesAll[iFit].tag = 0; }
                }

                // Classify edges

                for (int iEdge = 0, numEdges = mEllipseProperties.edgeIndices.size(); iEdge < numEdges; iEdge++)
                {
                    int jEdge = mEllipseProperties.edgeIndices[iEdge];
                    vEdgePropertiesAll[jEdge].tag = 2;
                }
            }
        }

        // Calculate average pixel brightness of accepted edge points

        double intensitySum = 0;
        int numEdges        = 0;

        for (int iEdge = 0, numEdgesAll = vEdgePropertiesAll.size(); iEdge < numEdgesAll; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            if (mEdgeProperties.tag == 2)
            {
                intensitySum += (double) mEdgeProperties.intensity;
                numEdges++;
            }
        }

        if (numEdges > 0) { mDetectionPropertiesNew.v.averageIntensity = intensitySum / numEdges; }

        /////////////////////////////////////////////////////////////////
        /////////////////////// SAVING DATA  ////////////////////////////
        /////////////////////////////////////////////////////////////////

        mDetectionPropertiesNew.m.edgePropertiesAll    = vEdgePropertiesAll; // edge data
        mDetectionPropertiesNew.m.ellipsePropertiesAll = vEllipsePropertiesAll; // ellipse data

        // Save parameters

        mDetectionPropertiesNew.v.PUPIL_DETECTED = mEllipseProperties.PUPIL_DETECTED;

        // For draw functions

        mDetectionPropertiesNew.m.outerAOI = outerAOI;
        mDetectionPropertiesNew.m.innerAOI = innerAOI;
        mDetectionPropertiesNew.m.glintAOI = glintAOI;

        mDetectionPropertiesNew.m.cannyEdgeIndices    = edgeIndices;
        mDetectionPropertiesNew.m.ellipseCoefficients = mEllipseProperties.coefficients;
    }
    else
    {
        mDetectionPropertiesNew.m.ERROR_DETECTED = true;
    }

    mDetectionProperties.p.alphaCertainty = 0.2; // turn this into adjustable parameter

    if (!mEllipseProperties.PUPIL_DETECTED) // pupil not detected
    {
        // Running averages

        mDetectionPropertiesNew.v.momentumAspectRatio   = mDetectionProperties.v.momentumAspectRatio   * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumCircumference = mDetectionProperties.v.momentumCircumference * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumWidth         = mDetectionProperties.v.momentumWidth         * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumHeight        = mDetectionProperties.v.momentumHeight        * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumXPos          = mDetectionProperties.v.momentumXPos          * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumYPos          = mDetectionProperties.v.momentumYPos          * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.predictionAspectRatio   = mDetectionProperties.v.predictionAspectRatio   + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.averageAspectRatio   - mDetectionProperties.v.predictionAspectRatio);
        mDetectionPropertiesNew.v.predictionCircumference = mDetectionProperties.v.predictionCircumference + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.averageCircumference - mDetectionProperties.v.predictionCircumference);
        mDetectionPropertiesNew.v.predictionWidth         = mDetectionProperties.v.predictionWidth         + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.averageWidth         - mDetectionProperties.v.predictionWidth);
        mDetectionPropertiesNew.v.predictionHeight        = mDetectionProperties.v.predictionHeight        + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.averageHeight        - mDetectionProperties.v.predictionHeight);

        mDetectionPropertiesNew.v.predictionXPos          = mDetectionProperties.v.predictionXPos          + mDetectionProperties.p.alphaPosition * (outerAOI.xPos + 0.5 * outerAOI.wdth - mDetectionProperties.v.predictionXPos) + mDetectionProperties.v.momentumXPos;
        mDetectionPropertiesNew.v.predictionYPos          = mDetectionProperties.v.predictionYPos          + mDetectionProperties.p.alphaPosition * (outerAOI.yPos + 0.5 * outerAOI.hght - mDetectionProperties.v.predictionYPos) + mDetectionProperties.v.momentumYPos;

        mDetectionPropertiesNew.v.certaintyPosition = mDetectionPropertiesNew.v.certaintyPosition - mDetectionProperties.v.certaintyPosition * mDetectionProperties.p.alphaCertainty;
        mDetectionPropertiesNew.v.certaintyFeatures = mDetectionPropertiesNew.v.certaintyFeatures - mDetectionProperties.v.certaintyFeatures * mDetectionProperties.p.alphaCertainty;
    }
    else // pupil detected
    {
        // Exact values

        mDetectionPropertiesNew.v.exactAspectRatio   = mEllipseProperties.aspectRatio;
        mDetectionPropertiesNew.v.exactCircumference = mEllipseProperties.circumference;

        mDetectionPropertiesNew.v.exactXPos = mEllipseProperties.xPos + outerAOI.xPos;
        mDetectionPropertiesNew.v.exactYPos = mEllipseProperties.yPos + outerAOI.yPos;

        // Running averages

        double changeAspectRatio   = mDetectionPropertiesNew.v.predictionAspectRatio    - mDetectionProperties.v.predictionAspectRatio;
        double changeCircumference = mDetectionPropertiesNew.v.predictionCircumference  - mDetectionProperties.v.predictionCircumference;
        double changeWidth         = mDetectionPropertiesNew.v.predictionWidth          - mDetectionProperties.v.predictionWidth;
        double changeHeight        = mDetectionPropertiesNew.v.predictionHeight         - mDetectionProperties.v.predictionHeight;
        double changeXPosition     = mDetectionPropertiesNew.v.predictionXPos           - mDetectionProperties.v.predictionXPos;
        double changeYPosition     = mDetectionPropertiesNew.v.predictionYPos           - mDetectionProperties.v.predictionYPos;

        double errorAspectRatio   = mEllipseProperties.aspectRatio      - mDetectionProperties.v.predictionAspectRatio;
        double errorCircumference = mEllipseProperties.circumference    - mDetectionProperties.v.predictionCircumference;
        double errorWidth         = mEllipseProperties.width            - mDetectionProperties.v.predictionWidth;
        double errorHeight        = mEllipseProperties.height           - mDetectionProperties.v.predictionHeight;
        double errorXPosition     = mDetectionPropertiesNew.v.exactXPos - mDetectionProperties.v.predictionXPos;
        double errorYPosition     = mDetectionPropertiesNew.v.exactYPos - mDetectionProperties.v.predictionYPos;

        mDetectionPropertiesNew.v.predictionAspectRatio   = mDetectionProperties.v.predictionAspectRatio   + mDetectionProperties.p.alphaFeatures * errorAspectRatio   + mDetectionProperties.v.momentumAspectRatio;
        mDetectionPropertiesNew.v.predictionCircumference = mDetectionProperties.v.predictionCircumference + mDetectionProperties.p.alphaFeatures * errorCircumference + mDetectionProperties.v.momentumCircumference;
        mDetectionPropertiesNew.v.predictionWidth         = mDetectionProperties.v.predictionWidth         + mDetectionProperties.p.alphaFeatures * errorWidth         + mDetectionProperties.v.momentumWidth;
        mDetectionPropertiesNew.v.predictionHeight        = mDetectionProperties.v.predictionHeight        + mDetectionProperties.p.alphaFeatures * errorHeight        + mDetectionProperties.v.momentumHeight;

        mDetectionPropertiesNew.v.predictionXPos          = mDetectionProperties.v.predictionXPos + mDetectionProperties.p.alphaPosition * errorXPosition     + mDetectionProperties.v.momentumXPos;
        mDetectionPropertiesNew.v.predictionYPos          = mDetectionProperties.v.predictionYPos + mDetectionProperties.p.alphaPosition * errorYPosition     + mDetectionProperties.v.momentumYPos;

        mDetectionPropertiesNew.v.averageAspectRatio   = mDetectionProperties.v.averageAspectRatio   + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.predictionAspectRatio   - mDetectionProperties.v.averageAspectRatio);
        mDetectionPropertiesNew.v.averageCircumference = mDetectionProperties.v.averageCircumference + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.predictionCircumference - mDetectionProperties.v.averageCircumference);
        mDetectionPropertiesNew.v.averageWidth         = mDetectionProperties.v.averageWidth         + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.predictionWidth         - mDetectionProperties.v.averageWidth);
        mDetectionPropertiesNew.v.averageHeight        = mDetectionProperties.v.averageHeight        + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.predictionHeight        - mDetectionProperties.v.averageHeight);
        mDetectionPropertiesNew.v.averageIntensity     = mDetectionProperties.v.averageIntensity     + mDetectionProperties.p.alphaAverage * (mDetectionPropertiesNew.v.averageIntensity     - mDetectionProperties.v.averageIntensity);

        mDetectionPropertiesNew.v.momentumAspectRatio   =  changeAspectRatio   * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumCircumference =  changeCircumference * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumWidth         =  changeWidth         * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumHeight        =  changeHeight        * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumXPos          =  changeXPosition     * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.momentumYPos          =  changeYPosition     * mDetectionProperties.p.alphaMomentum;

        // Determine whether prior pupil characteristics should be used in next frame

        double displacement       = sqrt(changeXPosition * changeXPosition + changeYPosition * changeYPosition);
        double certaintyPosition  = 1 - 1 / (1 + exp(-displacement + mDetectionProperties.p.changeThresholdPosition));

        double factor = mDetectionProperties.p.changeThresholdCircumference / mDetectionProperties.p.changeThresholdAspectRatio;
        double certaintyAspectRatio   = 1 - 1 / (1 + exp(factor * (mDetectionProperties.p.changeThresholdAspectRatio   - std::abs(changeAspectRatio))));
        double certaintyCircumference = 1 - 1 / (1 + exp(         (mDetectionProperties.p.changeThresholdCircumference - std::abs(changeCircumference))));
        double certaintyFeatures      = 0.5 * (certaintyAspectRatio + certaintyCircumference);

        mDetectionPropertiesNew.v.certaintyPosition = mDetectionProperties.v.certaintyPosition + (certaintyPosition - mDetectionProperties.v.certaintyPosition) * mDetectionProperties.p.alphaCertainty;
        mDetectionPropertiesNew.v.certaintyFeatures = mDetectionProperties.v.certaintyFeatures + (certaintyFeatures - mDetectionProperties.v.certaintyFeatures) * mDetectionProperties.p.alphaCertainty;
    }

    // Calculate new limits

    int imgSize;
    if (imageWdth > imageHght) { imgSize = imageWdth; }
    else                       { imgSize = imageHght; }

    int AOISize;
    if (outerAOI.wdth > outerAOI.hght) { AOISize = outerAOI.wdth; }
    else                               { AOISize = outerAOI.hght; }

    double rangeSearchRadius                 = imgSize - AOISize;
    double rangeChangeThresholdAspectRatio   = 1.0 - mDetectionProperties.p.changeThresholdAspectRatio;
    double rangeChangeThresholdCircumference = mDetectionProperties.p.circumferenceMax - mDetectionProperties.p.changeThresholdCircumference;
    double rangeChangeThresholdPosition      = imgSize - mDetectionProperties.p.changeThresholdPosition;
    double rangeCurvatureOffset              = 180 - mDetectionProperties.p.curvatureOffset;

    mDetectionPropertiesNew.v.changeThresholdAspectRatio   = rangeChangeThresholdAspectRatio   * (1 - mDetectionPropertiesNew.v.certaintyFeatures) + mDetectionProperties.p.changeThresholdAspectRatio;
    mDetectionPropertiesNew.v.changeThresholdCircumference = rangeChangeThresholdCircumference * (1 - mDetectionPropertiesNew.v.certaintyFeatures) + mDetectionProperties.p.changeThresholdCircumference;
    mDetectionPropertiesNew.v.curvatureOffset              = rangeCurvatureOffset              * (1 - mDetectionPropertiesNew.v.certaintyFeatures) + mDetectionProperties.p.curvatureOffset;

    mDetectionPropertiesNew.v.changeThresholdPosition = rangeChangeThresholdPosition * (1 - mDetectionPropertiesNew.v.certaintyPosition)  + mDetectionProperties.p.changeThresholdPosition;
    mDetectionPropertiesNew.v.searchRadius            = round(rangeSearchRadius      * (1 - mDetectionPropertiesNew.v.certaintyPosition)) + AOISize;

    return mDetectionPropertiesNew;
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
