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
    
    double minPupilIntensity = pow(10, 10); // arbitrarily large value
    
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

bool findEdge(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int& startEdgeIndex, double pupilXCentre, double pupilYCentre, bool USE_PRIOR_INFORMATION)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    // Find a starting edge point

    bool NEW_EDGE_FOUND = false;

    if (USE_PRIOR_INFORMATION)
    {
        for (int m = 0; m < 8 && !NEW_EDGE_FOUND; m++) // starburst
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
    }
    else
    {
        int mAOIArea = mAOI.wdth * mAOI.hght;

        for (int iEdgePoint = startEdgeIndex; iEdgePoint < mAOIArea && !NEW_EDGE_FOUND; iEdgePoint++)
        {
            if (cannyEdgeVector[iEdgePoint] == 1)
            {
                startEdgeIndex = iEdgePoint;
                NEW_EDGE_FOUND = true;
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

std::vector<edgeProperties> edgeSelection(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, double pupilXCentre, double pupilYCentre, bool USE_PRIOR_INFORMATION)
{
    std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all edges

    std::vector<int> startIndices;
    int startIndex = 0;

    while(true)
    {
        if (!findEdge (cannyEdgeVector, mAOI, startIndex, pupilXCentre, pupilYCentre, USE_PRIOR_INFORMATION)) { break; } // no (more) edges found
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

std::vector<int> findBreakPoints(std::vector<double> curvatures, double curvatureLowerLimit, double curvatureUpperLimit)
{
    int edgeSize = curvatures.size();

    // find breakpoints

    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(0); // add first point

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeSize - curvatureWindowLength; iEdgePoint++)
    {
        double curvature = curvatures[iEdgePoint];

        if (curvature >= curvatureUpperLimit || curvature <= curvatureLowerLimit)
        {
            breakPoints.push_back(iEdgePoint);
        }
    }

    breakPoints.push_back(edgeSize - 1); // add last point

    // evaluate each partial edge

    std::sort(breakPoints.begin(), breakPoints.end());

    return breakPoints;
}

std::vector<edgeProperties> edgeSegmentation(edgeProperties mEdgeProperties, std::vector<int> breakPoints)
{
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

std::vector<int> findEdgeIntensities(const cv::Mat& img, const std::vector<int>& edgePoints, const std::vector<double>& edgePointXNormals, const std::vector<double>& edgePointYNormals,  AOIProperties mAOI)
{
    int edgeLength = edgePoints.size();

    // calculate pixel intensities within inner curve of edge points

    std::vector<int> edgeIntensities(edgeLength);

    uchar *ptr_img = img.data;

    for (int iEdgePoint = 0; iEdgePoint < edgeLength; iEdgePoint++)
    {
        int edgePointIndex = edgePoints[iEdgePoint];
        int edgePointXPos  = edgePointIndex % mAOI.wdth;
        int edgePointYPos  = (edgePointIndex - edgePointXPos) / mAOI.wdth;

        int offsetXPos = edgePointXPos + edgeIntensitiesPositionOffset * ceil2(edgePointXNormals[iEdgePoint]);
        int offsetYPos = edgePointYPos + edgeIntensitiesPositionOffset * ceil2(edgePointYNormals[iEdgePoint]);

        if (offsetXPos < 0 || offsetXPos > mAOI.wdth || offsetYPos < 0 || offsetYPos > mAOI.hght)
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

        lengthTotal += sqrt(pow(dX,2) + pow(dY,2));

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

        edgePointRadii[iEdgePoint] = sqrt(pow(dX, 2) + pow(dY, 2));
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

std::vector<int> edgeClassification(detectionProperties mDetectionProperties, const std::vector<edgeProperties>& vEdgePropertiesAll, bool USE_PRIOR_INFORMATION)
{
    int numEdgesMax = mDetectionProperties.p.ellipseFitNumberMaximum;
    int numEdges    = vEdgePropertiesAll.size();
    if (numEdgesMax > numEdges) { numEdgesMax = numEdges; }

    // Classify edges based on score

    std::vector<double> totalScores(numEdges);

    const double norm = scoreFactorRadius + scoreFactorCurvature + scoreFactorCircumference + scoreFactorIntensity;

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        double dRadius          = mDetectionProperties.v.circumferencePrediction / (2 * M_PI) - vEdgePropertiesAll[iEdge].radius;
        double dCircumference   = mDetectionProperties.v.circumferencePrediction - vEdgePropertiesAll[iEdge].length;
        double dCurvature       = mDetectionProperties.v.edgeCurvaturePrediction - vEdgePropertiesAll[iEdge].curvature;
        double intensity        = vEdgePropertiesAll[iEdge].intensity;

        double scoreRadius          = 0;
        double scoreCurvature       = 0;
        double scoreCircumference   = 0;
        double scoreIntensity       = 0;

        if (USE_PRIOR_INFORMATION)
        {
            scoreRadius        = scoreFactorRadius        * (exp(-(pow(   dRadius, 2) / 20)));
            scoreCurvature     = scoreFactorCurvature     * (exp(-(pow(dCurvature, 2) / 20)));
            scoreCircumference = scoreFactorCircumference * (1 - 1 / (1 + exp(-0.10 * (dCircumference - 75))));
        }

        scoreIntensity = scoreFactorIntensity * (1 - 1 / (1 + exp(-0.25 * (intensity - 50))));

        totalScores[iEdge] = (scoreRadius + scoreCurvature + scoreCircumference + scoreIntensity) / norm;
    }

    // Only pick edges above threshold

    std::vector<int> pupilEdges;

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        if (totalScores[iEdge] >= mDetectionProperties.p.scoreThreshold || !USE_PRIOR_INFORMATION)
        {
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

    double minEigenValue   = pow(10, 50); // arbitrarily large number
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

std::vector<ellipseProperties> getEllipseFits(const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI, detectionProperties mDetectionProperties, bool USE_PRIOR_INFORMATION)
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

            if (USE_PRIOR_INFORMATION)
            {
                if (edgeSetLength < mDetectionProperties.v.circumferencePrediction * mDetectionProperties.p.edgeLengthFraction)
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

            if (USE_PRIOR_INFORMATION)
            {
                double dX = mEllipsePropertiesNew.xPos - (mDetectionProperties.v.xPosPrediction - mAOI.xPos);
                double dY = mEllipsePropertiesNew.yPos - (mDetectionProperties.v.yPosPrediction - mAOI.yPos);
                double dR = sqrt(pow(dX,2) + pow(dY,2));

                if (dR > mDetectionProperties.v.thresholdDisplacementChange) { continue; } // no large ellipse displacements
                if (std::abs(mEllipsePropertiesNew.circumference - mDetectionProperties.v.circumferencePrediction) > mDetectionProperties.v.thresholdCircumferenceChange) { continue; } // no large ellipse size changes
                if (std::abs(mEllipsePropertiesNew.aspectRatio   - mDetectionProperties.v.aspectRatioPrediction  ) > mDetectionProperties.v.thresholdAspectRatioChange  ) { continue; } // no large ellipse shape changes
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

int ellipseFitFilter(detectionProperties mDetectionProperties, std::vector<ellipseProperties> vEllipseProperties, bool USE_PRIOR_INFORMATION)
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

        double dx = mEllipseProperties.xPos - mDetectionProperties.v.xPosPrediction;
        double dy = mEllipseProperties.yPos - mDetectionProperties.v.yPosPrediction;
        double displacementChange  = sqrt(pow(dx,2) + pow(dy,2));
        double circumferenceChange = (std::abs(mEllipseProperties.circumference - mDetectionProperties.v.circumferencePrediction));
        double aspectRatioChange   = (std::abs(mEllipseProperties.aspectRatio   - mDetectionProperties.v.aspectRatioPrediction));
        double lengthDifference    = (std::abs(mEllipseProperties.edgeLength    - mDetectionProperties.v.circumferencePrediction));
        double fitError            = mEllipseProperties.fitError;

        if (USE_PRIOR_INFORMATION)
        {
            scoreAspectRatio   = (-maxScoreAspectRatio   / mDetectionProperties.p.aspectRatioChangeThreshold)   * aspectRatioChange   + maxScoreAspectRatio;
            scoreCircumference = (-maxScoreCircumference / mDetectionProperties.p.circumferenceChangeThreshold) * circumferenceChange + maxScoreCircumference;
            scoreDisplacement  = (-maxScoreDisplacement  / mDetectionProperties.p.displacementChangeThreshold)  * displacementChange  + maxScoreDisplacement;
        }

        double scoreFitError = (-maxScoreFitError   / mDetectionProperties.p.ellipseFitErrorMaximum)  * fitError         + maxScoreFitError;
        double scoreLength   = (-maxScoreLength * 2 / mDetectionProperties.v.circumferencePrediction) * lengthDifference + maxScoreLength;

        if (scoreCircumference  < 0) { scoreCircumference   = 0; }
        if (scoreAspectRatio    < 0) { scoreAspectRatio     = 0; }
        if (scoreFitError       < 0) { scoreFitError        = 0; }
        if (scoreLength         < 0) { scoreLength          = 0; }

        featureChange[iFit] = scoreAspectRatio + scoreCircumference + scoreDisplacement + scoreFitError + scoreLength;
    }

    int acceptedFitIndex = std::distance(featureChange.begin(), std::max_element(featureChange.begin(), featureChange.end()));

    return acceptedFitIndex;
}

void checkVariableLimits(detectionProperties& mDetectionProperties, int imgWdth, int imgHght, AOIProperties mAOI)
{
    if (mDetectionProperties.v.searchRadius > 0.5 * imgWdth)
    {   mDetectionProperties.v.searchRadius = 0.5 * imgWdth; }
    else if (mDetectionProperties.v.searchRadius < (0.5 * mAOI.wdth))
    {        mDetectionProperties.v.searchRadius = ceil(0.5 * mAOI.wdth); }

    if (mDetectionProperties.v.searchRadius > 0.5 * imgHght)
    {   mDetectionProperties.v.searchRadius = 0.5 * imgHght; }
    else if (mDetectionProperties.v.searchRadius < (0.5 * mAOI.hght))
    {        mDetectionProperties.v.searchRadius = ceil(0.5 * mAOI.hght); }

    if (mDetectionProperties.v.thresholdAspectRatioChange > 1.0)
    {   mDetectionProperties.v.thresholdAspectRatioChange = 1.0; }
    else if (mDetectionProperties.v.thresholdAspectRatioChange < mDetectionProperties.p.aspectRatioChangeThreshold)
    {        mDetectionProperties.v.thresholdAspectRatioChange = mDetectionProperties.p.aspectRatioChangeThreshold; }

    if (mDetectionProperties.v.thresholdCircumferenceChange > mDetectionProperties.p.circumferenceMax)
    {   mDetectionProperties.v.thresholdCircumferenceChange = mDetectionProperties.p.circumferenceMax; }
    else if (mDetectionProperties.v.thresholdCircumferenceChange < mDetectionProperties.p.circumferenceChangeThreshold)
    {        mDetectionProperties.v.thresholdCircumferenceChange = mDetectionProperties.p.circumferenceChangeThreshold; }

    if (mDetectionProperties.v.thresholdDisplacementChange > imgWdth)
    {   mDetectionProperties.v.thresholdDisplacementChange = imgWdth; }
    else if (mDetectionProperties.v.thresholdDisplacementChange < mDetectionProperties.p.displacementChangeThreshold)
    {        mDetectionProperties.v.thresholdDisplacementChange = mDetectionProperties.p.displacementChangeThreshold; }

    if (mDetectionProperties.v.curvatureOffset > 180)
    {   mDetectionProperties.v.curvatureOffset = 180; }
    else if (mDetectionProperties.v.curvatureOffset < mDetectionProperties.p.curvatureOffset)
    {        mDetectionProperties.v.curvatureOffset = mDetectionProperties.p.curvatureOffset; }

    if (mDetectionProperties.v.priorCertainty > certaintyUpperLimit)
    {   mDetectionProperties.v.priorCertainty = certaintyUpperLimit; }
    else if (mDetectionProperties.v.priorCertainty < certaintyLowerLimit)
    {        mDetectionProperties.v.priorCertainty = certaintyLowerLimit; }
}

detectionProperties pupilDetection(const cv::Mat& imageOriginalBGR, detectionProperties mDetectionProperties)
{
    // Define some variables

    detectionProperties mDetectionPropertiesNew = mDetectionProperties; // new properties for new frame
    mDetectionPropertiesNew.m.ERROR_DETECTED = false;

    ellipseProperties mEllipseProperties;
    mEllipseProperties.PUPIL_DETECTED = false;

    // Define search area

    int imageWdth = imageOriginalBGR.cols;
    int imageHght = imageOriginalBGR.rows;

    AOIProperties searchAOI;

    searchAOI.xPos = round(mDetectionProperties.v.xPosPrediction - mDetectionProperties.v.searchRadius);
    searchAOI.yPos = round(mDetectionProperties.v.yPosPrediction - mDetectionProperties.v.searchRadius);

    int searchEndX = round(mDetectionProperties.v.xPosPrediction + mDetectionProperties.v.searchRadius);
    int searchEndY = round(mDetectionProperties.v.yPosPrediction + mDetectionProperties.v.searchRadius);

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
    innerAOI.wdth = round(pupilHaarReductionFactor * mDetectionProperties.v.widthPrediction);
    innerAOI.hght = round(pupilHaarReductionFactor * mDetectionProperties.v.heightPrediction);
    if (innerAOI.wdth > searchAOI.wdth) { innerAOI.wdth = searchAOI.wdth; }
    if (innerAOI.hght > searchAOI.hght) { innerAOI.hght = searchAOI.hght; }

    AOIProperties outerAOI = mDetectionProperties.m.outerAOI;

    if (innerAOI.wdth > 0 && innerAOI.hght > 0)
    {
        // Needed for offline mode when threshold is too low

        checkVariableLimits(mDetectionProperties, imageWdth, imageHght, outerAOI);

        // Check if prior pupil information should be used

        bool USE_PRIOR_INFORMATION = false;
        if (mDetectionProperties.v.priorCertainty >= certaintyThreshold) { USE_PRIOR_INFORMATION = true; }

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
        std::vector<int> edgeIndices = getEdgeIndices(cannyEdgesSharpened, 1);

        /////////////////////////////////////////////////////////////////////////////
        //////////////////////////// EDGE SELECTION   ///////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        double xPosPredictionRelative = mDetectionProperties.v.xPosPrediction - outerAOI.xPos;
        double yPosPredictionRelative = mDetectionProperties.v.yPosPrediction - outerAOI.yPos;

        std::vector<edgeProperties> vEdgePropertiesAll = edgeSelection(cannyEdgesSharpened, outerAOI, xPosPredictionRelative, yPosPredictionRelative, USE_PRIOR_INFORMATION);
        vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

        /////////////////////////////////////////////////////////////////////////////
        //////////////////////////// EDGE SEGMENTATION   ////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        // Calculate curvature limits

        int arrayWidth = arrayCircumferences.size();

        int arrayXPos = 0;
        for (int x = 0; x < arrayWidth; x++)
        {
            if (arrayCircumferences[x] < mDetectionProperties.v.circumferencePrediction)
            {
                arrayXPos = x;
                break;
            }
        }

        int arrayHeight = arrayAspectRatios.size();

        int arrayYPos = 0;
        for (int y = 0; y < arrayHeight; y++)
        {
            if (arrayAspectRatios[y] < mDetectionProperties.v.aspectRatioPrediction)
            {
                arrayYPos = y;
                break;
            }
        }

        double curvatureUpperLimit = arrayCurvatureMax[arrayXPos * arrayWidth + arrayYPos] + mDetectionProperties.v.curvatureOffset;
        double curvatureLowerLimit = arrayCurvatureMin[arrayXPos * arrayWidth + arrayYPos] - mDetectionProperties.v.curvatureOffset;

        // Do segmentation

        std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];

            int edgeSize = mEdgeProperties.pointIndices.size();

            std::vector<double> edgeXTangents(edgeSize);
            std::vector<double> edgeYTangents(edgeSize);

            calculateEdgeDirections(mEdgeProperties.pointIndices, edgeXTangents, edgeYTangents, outerAOI);
            calculateCurvatures(mEdgeProperties, edgeXTangents, edgeYTangents);

            // Do segmentation

            std::vector<int> breakPoints = findBreakPoints(mEdgeProperties.curvatures, curvatureLowerLimit, curvatureUpperLimit);
            std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentation(mEdgeProperties, breakPoints);

            vEdgePropertiesNew.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end()); // record edges
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
        vEdgePropertiesNew.clear();

        vEdgePropertiesAll = removeShortEdges(vEdgePropertiesAll);

        ////////////////////////////////////////////////////////////////////////
        //////////////////////// EDGE POINT FILTER   ///////////////////////////
        ////////////////////////////////////////////////////////////////////////

        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties  = vEdgePropertiesAll[iEdge];
            mEdgeProperties.intensities     = findEdgeIntensities(imageAOIGray, mEdgeProperties.pointIndices, mEdgeProperties.xnormals, mEdgeProperties.ynormals, outerAOI);
            mEdgeProperties.radii           = calculateEdgeRadii(mEdgeProperties, outerAOI, xPosPredictionRelative, yPosPredictionRelative);
            mEdgeProperties.gradients       = calculateRadialGradients(imageAOIGray, mEdgeProperties.pointIndices, xPosPredictionRelative, yPosPredictionRelative, mDetectionProperties.p.cannyKernelSize);

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

            mEdgeProperties.length    = calculateEdgeLength(mEdgeProperties.pointIndices, outerAOI);

            restoreEdgePoints(mEdgeProperties, cannyEdgesSharpened, outerAOI); // Restore some points

            mEdgeProperties.intensity = calculateMeanInt(mEdgeProperties.intensities);
            mEdgeProperties.radius    = calculateMean(mEdgeProperties.radii);
            mEdgeProperties.gradient  = calculateMeanInt(mEdgeProperties.gradients);

            calculateCurvatureRange(mEdgeProperties);

            mEdgeProperties.size      = mEdgeProperties.pointIndices.size();

            mEdgeProperties.index     = jEdge;
            mEdgeProperties.tag       = 0;

            vEdgePropertiesNew.push_back(mEdgeProperties);

            jEdge++;
        }

        vEdgePropertiesAll = vEdgePropertiesNew;
        vEdgePropertiesNew.clear();

        // Do edge classification

        std::vector<ellipseProperties> vEllipsePropertiesAll;

        if (vEdgePropertiesAll.size() > 0) // ignore empty edge collections
        {
            mDetectionProperties.v.edgeCurvaturePrediction = 0.5 * (curvatureUpperLimit + curvatureLowerLimit);

            std::vector<int> acceptedEdges = edgeClassification(mDetectionProperties, vEdgePropertiesAll, USE_PRIOR_INFORMATION);

            for (int iEdge = 0, numEdges = acceptedEdges.size(); iEdge < numEdges; iEdge++) // grab accepted edges
            {
                int jEdge = acceptedEdges[iEdge];
                vEdgePropertiesAll[jEdge].tag = 1; // new tag
                vEdgePropertiesNew.push_back(vEdgePropertiesAll[jEdge]);
            }

            vEdgePropertiesAll = vEdgePropertiesNew;
            vEdgePropertiesNew.clear();


            //////////////////////////////////////////////////////////////////
            /////////////////////// ELLIPSE FITTING  /////////////////////////
            //////////////////////////////////////////////////////////////////

            vEllipsePropertiesAll = getEllipseFits(vEdgePropertiesAll, outerAOI, mDetectionProperties, USE_PRIOR_INFORMATION); // ellipse fitting

            int numFits = vEllipsePropertiesAll.size();

            int acceptedFitIndex = 0;

            if (numFits > 0)
            {
                mEllipseProperties.PUPIL_DETECTED = true;

                if (numFits > 1) { acceptedFitIndex = ellipseFitFilter(mDetectionProperties, vEllipsePropertiesAll, USE_PRIOR_INFORMATION); } // grab best fit

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

        /////////////////////////////////////////////////////////////////
        /////////////////////// SAVING DATA  ////////////////////////////
        /////////////////////////////////////////////////////////////////

        mDetectionPropertiesNew.m.edgePropertiesAll    = vEdgePropertiesAll; // edge data
        mDetectionPropertiesNew.m.ellipsePropertiesAll = vEllipsePropertiesAll; // ellipse data

        // Save parameters

        mDetectionPropertiesNew.v.edgeCurvaturePrediction = mDetectionProperties.v.edgeCurvaturePrediction;

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

    // For running averages

    if (!mEllipseProperties.PUPIL_DETECTED) // pupil not detected
    {
        mDetectionPropertiesNew.v.aspectRatioAverage    = mDetectionProperties.v.aspectRatioAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.aspectRatioPrediction - mDetectionProperties.v.aspectRatioAverage);
        mDetectionPropertiesNew.v.aspectRatioMomentum   = mDetectionProperties.v.aspectRatioMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.aspectRatioPrediction = mDetectionProperties.v.aspectRatioPrediction + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.aspectRatioAverage - mDetectionProperties.v.aspectRatioPrediction);

        mDetectionPropertiesNew.v.circumferenceAverage    = mDetectionProperties.v.circumferenceAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.circumferencePrediction - mDetectionProperties.v.circumferenceAverage);
        mDetectionPropertiesNew.v.circumferenceMomentum   = mDetectionProperties.v.circumferenceMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.circumferencePrediction = mDetectionProperties.v.circumferencePrediction + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.circumferenceAverage - mDetectionProperties.v.circumferencePrediction);

        mDetectionPropertiesNew.v.curvatureOffset = mDetectionProperties.v.curvatureOffset * (1 / mDetectionProperties.p.alphaFeatures);

        mDetectionPropertiesNew.v.heightAverage    = mDetectionProperties.v.heightAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.heightPrediction - mDetectionProperties.v.heightAverage);
        mDetectionPropertiesNew.v.heightMomentum   = mDetectionProperties.v.heightMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.heightPrediction = mDetectionProperties.v.heightPrediction + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.heightAverage - mDetectionProperties.v.heightPrediction);

        mDetectionPropertiesNew.v.radiusMomentum   = mDetectionProperties.v.radiusMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.radiusPrediction = mDetectionProperties.v.circumferencePrediction / (2 * M_PI);

        mDetectionPropertiesNew.v.searchRadius = mDetectionProperties.v.searchRadius * (1 / mDetectionProperties.p.alphaPosition);

        mDetectionPropertiesNew.v.thresholdAspectRatioChange   = mDetectionProperties.v.thresholdAspectRatioChange   * (1 / mDetectionProperties.p.alphaFeatures);
        mDetectionPropertiesNew.v.thresholdCircumferenceChange = mDetectionProperties.v.thresholdCircumferenceChange * (1 / mDetectionProperties.p.alphaFeatures);
        mDetectionPropertiesNew.v.thresholdDisplacementChange  = mDetectionProperties.v.thresholdDisplacementChange  * (1 / mDetectionProperties.p.alphaPosition);

        mDetectionPropertiesNew.v.widthAverage    = mDetectionProperties.v.widthAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.widthPrediction - mDetectionProperties.v.widthAverage);
        mDetectionPropertiesNew.v.widthMomentum   = mDetectionProperties.v.widthMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.widthPrediction = mDetectionProperties.v.widthPrediction + mDetectionProperties.p.alphaFeatures * (mDetectionPropertiesNew.v.widthAverage - mDetectionProperties.v.widthPrediction);

        mDetectionPropertiesNew.v.xPosPrediction = mDetectionProperties.v.xPosPrediction + mDetectionProperties.p.alphaPosition * (outerAOI.xPos + 0.5 * outerAOI.wdth - mDetectionProperties.v.xPosPrediction) + mDetectionProperties.v.xVelocity;
        mDetectionPropertiesNew.v.yPosPrediction = mDetectionProperties.v.yPosPrediction + mDetectionProperties.p.alphaPosition * (outerAOI.yPos + 0.5 * outerAOI.hght - mDetectionProperties.v.yPosPrediction) + mDetectionProperties.v.yVelocity;

        mDetectionPropertiesNew.v.xVelocity     = mDetectionProperties.v.xVelocity * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.yVelocity     = mDetectionProperties.v.yVelocity * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.priorCertainty = mDetectionProperties.v.priorCertainty * mDetectionProperties.p.alphaFeatures;
    }
    else // pupil detected
    {
        mDetectionPropertiesNew.v.aspectRatioExact   = mEllipseProperties.aspectRatio;
        mDetectionPropertiesNew.v.circumferenceExact = mEllipseProperties.circumference;

        mDetectionPropertiesNew.v.aspectRatioPrediction =  mDetectionProperties.v.aspectRatioPrediction + mDetectionProperties.p.alphaFeatures * (mEllipseProperties.aspectRatio - mDetectionProperties.v.aspectRatioPrediction) + mDetectionProperties.v.aspectRatioMomentum;
        mDetectionPropertiesNew.v.aspectRatioAverage    =  mDetectionProperties.v.aspectRatioAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.aspectRatioPrediction - mDetectionProperties.v.aspectRatioAverage);
        mDetectionPropertiesNew.v.aspectRatioMomentum   = (mDetectionProperties.v.aspectRatioMomentum   + mDetectionPropertiesNew.v.aspectRatioPrediction - mDetectionProperties.v.aspectRatioPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.circumferencePrediction =  mDetectionProperties.v.circumferencePrediction + mDetectionProperties.p.alphaFeatures * (mEllipseProperties.circumference - mDetectionProperties.v.circumferencePrediction) + mDetectionProperties.v.circumferenceMomentum;
        mDetectionPropertiesNew.v.circumferenceAverage    =  mDetectionProperties.v.circumferenceAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.circumferencePrediction - mDetectionProperties.v.circumferenceAverage);
        mDetectionPropertiesNew.v.circumferenceMomentum   = (mDetectionProperties.v.circumferenceMomentum   + mDetectionPropertiesNew.v.circumferencePrediction - mDetectionProperties.v.circumferencePrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.curvatureOffset = mDetectionProperties.v.curvatureOffset * mDetectionProperties.p.alphaFeatures;

        mDetectionPropertiesNew.v.heightPrediction =  mDetectionProperties.v.heightPrediction + mDetectionProperties.p.alphaFeatures * (mEllipseProperties.height - mDetectionProperties.v.heightPrediction) + mDetectionProperties.v.heightMomentum;
        mDetectionPropertiesNew.v.heightAverage    =  mDetectionProperties.v.heightAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.heightPrediction - mDetectionProperties.v.heightAverage);
        mDetectionPropertiesNew.v.heightMomentum   = (mDetectionProperties.v.heightMomentum   + mDetectionPropertiesNew.v.heightPrediction - mDetectionProperties.v.heightPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.radiusPrediction =  mDetectionProperties.v.radiusPrediction + mDetectionProperties.p.alphaFeatures * (mEllipseProperties.radius - mDetectionProperties.v.radiusPrediction) + mDetectionProperties.v.radiusMomentum;
        mDetectionPropertiesNew.v.radiusMomentum   = (mDetectionProperties.v.radiusMomentum + (mDetectionPropertiesNew.v.radiusPrediction - mDetectionProperties.v.radiusPrediction)) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.searchRadius = mDetectionProperties.v.searchRadius * mDetectionProperties.p.alphaPosition;

        mDetectionPropertiesNew.v.thresholdAspectRatioChange   = mDetectionProperties.v.thresholdAspectRatioChange   * mDetectionProperties.p.alphaFeatures;
        mDetectionPropertiesNew.v.thresholdCircumferenceChange = mDetectionProperties.v.thresholdCircumferenceChange * mDetectionProperties.p.alphaFeatures;
        mDetectionPropertiesNew.v.thresholdDisplacementChange  = mDetectionProperties.v.thresholdDisplacementChange  * mDetectionProperties.p.alphaPosition;

        mDetectionPropertiesNew.v.widthPrediction =  mDetectionProperties.v.widthPrediction + mDetectionProperties.p.alphaFeatures * (mEllipseProperties.width - mDetectionProperties.v.widthPrediction) + mDetectionProperties.v.widthMomentum;
        mDetectionPropertiesNew.v.widthAverage    =  mDetectionProperties.v.widthAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.widthPrediction - mDetectionProperties.v.widthAverage);
        mDetectionPropertiesNew.v.widthMomentum   = (mDetectionProperties.v.widthMomentum   + mDetectionPropertiesNew.v.widthPrediction - mDetectionProperties.v.widthPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.xPosExact = mEllipseProperties.xPos + outerAOI.xPos;
        mDetectionPropertiesNew.v.yPosExact = mEllipseProperties.yPos + outerAOI.yPos;

        mDetectionPropertiesNew.v.xPosPrediction = mDetectionProperties.v.xPosPrediction + mDetectionProperties.p.alphaPosition * (mDetectionPropertiesNew.v.xPosExact - mDetectionProperties.v.xPosPrediction) + mDetectionProperties.v.xVelocity;
        mDetectionPropertiesNew.v.yPosPrediction = mDetectionProperties.v.yPosPrediction + mDetectionProperties.p.alphaPosition * (mDetectionPropertiesNew.v.yPosExact - mDetectionProperties.v.yPosPrediction) + mDetectionProperties.v.yVelocity;

        mDetectionPropertiesNew.v.xVelocity = (mDetectionProperties.v.xVelocity + mDetectionPropertiesNew.v.xPosPrediction - mDetectionProperties.v.xPosPrediction) * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.yVelocity = (mDetectionProperties.v.yVelocity + mDetectionPropertiesNew.v.yPosPrediction - mDetectionProperties.v.yPosPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.priorCertainty = mDetectionProperties.v.priorCertainty * (1 / mDetectionProperties.p.alphaFeatures);

        // Grab pupil image

        if (SAVE_PUPIL_IMAGE)
        {
            int wdth = mEllipseProperties.width  * pupilImageFactor;
            int hght = mEllipseProperties.height * pupilImageFactor;

            if (wdth <= 0 || hght <= 0)
            {
                wdth = mEllipseProperties.radius * pupilImageFactor;
                hght = mEllipseProperties.radius * pupilImageFactor;
            }

            int xPos = round(mDetectionPropertiesNew.v.xPosExact - 0.5 * wdth);
            int yPos = round(mDetectionPropertiesNew.v.yPosExact - 0.5 * hght);

            if (xPos < 0)
            {
                imageWdth = imageWdth + xPos;
                xPos      = 0;
            }

            if (yPos < 0)
            {
                imageHght = imageHght + yPos;
                yPos      = 0;
            }

            if (xPos + wdth >= imageWdth)
            {
                int dx = xPos + wdth - (imageWdth - 1);
                xPos   = xPos + dx;
                wdth   = wdth - 2 * dx;
            }

            if (yPos + hght >= imageHght)
            {
                int dy = yPos + hght - (imageHght - 1);
                yPos   = yPos + dy;
                hght   = hght - 2 * dy;
            }

            cv::Rect pupilROI(xPos, yPos, wdth, hght);
            cv::Mat imageAOIBGR = imageOriginalBGR(pupilROI);
            cv::Mat imageAOIGray;
            cv::cvtColor(imageAOIBGR, imageAOIGray, cv::COLOR_BGR2GRAY);

            mDetectionPropertiesNew.m.imagePupil = imageAOIGray;
        }
    }

    checkVariableLimits(mDetectionPropertiesNew, imageWdth, imageHght, outerAOI);

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
